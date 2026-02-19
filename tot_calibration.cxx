// tot_calibration_v2.cxx
//
// Refactor of Tristan's tot_calibration.C to match the event/mapping style used in
// led_analysis.C / draw_waveform_reversed.cxx / gain_match.cxx / adc_calibration.cxx:
//  - Uses common_led.{h,cxx} mapping + signal helpers
//  - Reads adc/tot as flat buffers using TLeaf::GetLen() (dynamic sizes)
//  - Uses mapping CSV produced by generate_mapping.py
//  - Keeps the original purpose: calibrate ToT response vs "(ref - missing)" for ToT events
//
// IMPORTANT about the x-axis:
//  - In the original Tristan code this was (beam_energy[GeV] - noncentral_energy[GeV]).
//  - Here we generalize it to (label_value - noncentral_equivalent_value).
//  - If you pass label_value as LED voltage (V), we convert the non-central ADC sum into a
//    voltage-equivalent using adc_calibration output (ADC per 1.27V reference unit).
//
// Run example:
//   root -l -q 'common_led.cxx+ tot_calibration.cxx+ tot_calib_scan("data", "eeemcal_desy_dec2025_mapping_v2.csv", "outputs")'
//
// Or edit the run list inside tot_calib_scan() using:
//   std::vector<std::pair<int,float>> runs = { {48,1.29f}, {51,1.30f}, ... };
// Interactive example for one run:
// root -l -b
// .L common_led.cxx+
// .L tot_calibration.cxx+
// tot_calib_scan("data", "eeemcal_desy_dec2025_mapping_v2.csv", "outputs");

#include "common_led.h"

#include <map>
#include <array>
#include <vector>
#include <unordered_map>
#include <utility>
#include <string>
#include <iostream>
#include <cctype>
#include <cmath>
#include <algorithm>

#include <sys/stat.h>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLatex.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TParameter.h>
#include <TKey.h>
#include <TError.h>

static constexpr int SAMPLES_PER_CHANNEL = 20;
static constexpr int SIPMS_PER_CRYSTAL = 16;
static constexpr int MAX_NUM_CRYSTALS = 25;

// ----- knobs (keep same defaults as Tristan's original) -----
static float g_energy_fraction_cut = 0.30f; // used in COG selection
static long g_max_events = 1000000;
static float g_tot_min_cut = 0.0f; // if ToT sum < this -> set to 0

// crystals to ignore (e.g., Tristan ignored crystal 9)
// Edit this list to ignore additional crystals during the analysis.
static std::vector<int> g_skip_crystals = {9, 15};

static inline bool is_skipped_crystal(int cr)
{
    return std::find(g_skip_crystals.begin(), g_skip_crystals.end(), cr) != g_skip_crystals.end();
}

// COG ellipse cut (adjust if needed)
static bool g_do_cog_ellipse_cut = true;
static float g_cog_cx = 2.0f;
static float g_cog_cy = 2.0f;
static float g_cog_sx = 0.80f;
static float g_cog_sy = 0.80f;

// ADC calibration reference voltage used by adc_calibration.cxx
static float g_adc_ref_voltage_V = 1.27f;

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

static void mkdir_p(const char *dir)
{
    if (!dir)
        return;
    struct stat st;
    if (stat(dir, &st) != 0)
    {
        mkdir(dir, 0755);
    }
}

static void print_progress_25(int progress)
{
    std::cout << " [";
    for (int i = 0; i < MAX_NUM_CRYSTALS; i++)
        std::cout << (i < progress ? "*" : " ");
    std::cout << "]\r" << std::flush;
    if (progress >= MAX_NUM_CRYSTALS)
        std::cout << "\n";
}

// 5x5 geometry mapping used in your other macros.
static void build_xy_for_crystal(int crystal_id, float &x, float &y)
{
    // index -> crystal_id (row-major, top-left to bottom-right)
    int crystal_mapping[MAX_NUM_CRYSTALS] = {
        4, 9, 14, 19, 24,
        3, 8, 13, 18, 23,
        2, 7, 12, 17, 22,
        1, 6, 11, 16, 21,
        0, 5, 10, 15, 20};

    int found = -1;
    for (int i = 0; i < MAX_NUM_CRYSTALS; ++i)
    {
        if (crystal_mapping[i] == crystal_id)
        {
            found = i;
            break;
        }
    }

    if (found < 0)
    {
        x = -999;
        y = -999;
        return;
    }

    int row = found / 5; // 0..4
    int col = found % 5; // 0..4
    x = (float)col;
    y = (float)row;
}

static bool point_in_ellipse(float x, float y, float cx, float cy, float sx, float sy)
{
    if (sx <= 0 || sy <= 0)
        return true;
    const float dx = (x - cx) / sx;
    const float dy = (y - cy) / sy;
    return (dx * dx + dy * dy) <= 1.0f;
}

// Compute and optionally cut on COG, and optionally cut on central9 fraction.
static bool calculate_cog(TH2 *distribution,
                          const std::array<float, MAX_NUM_CRYSTALS> &signals)
{
    // total
    double sumE = 0.0;
    double sumX = 0.0;
    double sumY = 0.0;

    // fraction in central 3x3 around crystal 12
    const int central9[9] = {6, 7, 8, 11, 12, 13, 16, 17, 18};
    double sumE9 = 0.0;

    for (int cr = 0; cr < MAX_NUM_CRYSTALS; ++cr)
    {
        if (is_skipped_crystal(cr))
            continue;
        const float E = signals[cr];
        if (E <= 0)
            continue;

        float x, y;
        build_xy_for_crystal(cr, x, y);
        if (x < -100)
            continue;

        sumE += E;
        sumX += E * x;
        sumY += E * y;
    }

    for (int i = 0; i < 9; ++i)
    {
        const int cr = central9[i];
        if (is_skipped_crystal(cr))
            continue;
        if (cr < 0 || cr >= MAX_NUM_CRYSTALS)
            continue;
        const float E = signals[cr];
        if (E > 0)
            sumE9 += E;
    }

    if (sumE <= 0)
        return false;

    const float x_cog = (float)(sumX / sumE);
    const float y_cog = (float)(sumY / sumE);

    if (distribution)
        distribution->Fill(x_cog, y_cog);

    // central energy fraction cut
    const float frac9 = (sumE > 0) ? (float)(sumE9 / sumE) : 0.0f;
    if (frac9 < g_energy_fraction_cut)
        return false;

    // ellipse cut
    if (g_do_cog_ellipse_cut)
    {
        if (!point_in_ellipse(x_cog, y_cog, g_cog_cx, g_cog_cy, g_cog_sx, g_cog_sy))
            return false;
    }

    return true;
}

// Load ADC->ref calibration produced by adc_calibration.cxx
// Returns ADC counts per "1 ref unit" (where 1 ref unit == g_adc_ref_voltage_V).
static float load_adc_to_ref(const char *adc_calib_root, float ref_voltage)
{
    if (!adc_calib_root)
        return 0.0f;

    TFile *f = TFile::Open(adc_calib_root, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "Error: could not open adc calibration file: " << adc_calib_root << "\n";
        if (f)
        {
            f->Close();
            delete f;
        }
        return 0.0f;
    }

    // Try exact name first
    std::string exact = std::string(Form("mean_adc_to_ref_calibration_%.2fV", ref_voltage));
    TParameter<float> *p = (TParameter<float> *)f->Get(exact.c_str());

    // Fallback: search for any key with prefix
    if (!p)
    {
        TIter nextkey(f->GetListOfKeys());
        TKey *key;
        while ((key = (TKey *)nextkey()))
        {
            std::string name = key->GetName();
            if (name.rfind("mean_adc_to_ref_calibration_", 0) == 0)
            {
                p = (TParameter<float> *)f->Get(name.c_str());
                if (p)
                    break;
            }
        }
    }

    float val = 0.0f;
    if (p)
        val = p->GetVal();

    f->Close();
    delete f;

    return val;
}

// Load gain factors produced by gain_match.cxx
static void load_gain_factors(const char *gain_root,
                              TH1 *&gain_factors,
                              TH1 *&crystal_factor,
                              TFile *&gain_file_handle)
{
    gain_factors = nullptr;
    crystal_factor = nullptr;
    gain_file_handle = nullptr;

    if (gain_root)
    {
        gain_file_handle = TFile::Open(gain_root, "READ");
        if (gain_file_handle && !gain_file_handle->IsZombie())
        {
            gain_factors = (TH1 *)gain_file_handle->Get("gain_factors");
            crystal_factor = (TH1 *)gain_file_handle->Get("crystal_factor");
        }
    }

    // Fallback to unity
    if (!gain_factors)
    {
        gain_factors = new TH1F("gain_factors_unity", "Gain Factors;crystal*16+sipm;Gain", MAX_NUM_CRYSTALS * SIPMS_PER_CRYSTAL, 0, MAX_NUM_CRYSTALS * SIPMS_PER_CRYSTAL);
        for (int i = 1; i <= gain_factors->GetNbinsX(); ++i)
            gain_factors->SetBinContent(i, 1.0);
    }

    if (!crystal_factor)
    {
        crystal_factor = new TH1F("crystal_factor_unity", "Crystal Factors;crystal;Gain", MAX_NUM_CRYSTALS, 0, MAX_NUM_CRYSTALS);
        for (int i = 1; i <= crystal_factor->GetNbinsX(); ++i)
            crystal_factor->SetBinContent(i, 1.0);
    }
}

// -----------------------------------------------------------------------------
// Core: process one run
// -----------------------------------------------------------------------------

static bool tot_calibration_one(const char *filename,
                                float label_value,
                                const std::map<int, std::vector<int>> &mapping,
                                int central_crystal,
                                TH1 *gain_factors,
                                TH1 *crystal_factor,
                                float adc_per_refunit,
                                float ref_voltage,
                                const char *outdir,
                                TH1F *h_missing,
                                TH1F *h_tot,
                                TH2F *h_tot_vs_label_minus_missing,
                                TH2F *h_cog,
                                std::array<TH1F *, SIPMS_PER_CRYSTAL> &h_per_sipm_tot,
                                TH1 *widths_hist,
                                bool use_widths)
{
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie())
    {
        std::cerr << "Error: could not open " << filename << "\n";
        if (f)
        {
            f->Close();
            delete f;
        }
        return false;
    }

    TTree *tree = (TTree *)f->Get("events");
    if (!tree)
    {
        std::cerr << "Error: could not find TTree 'events' in " << filename << "\n";
        f->Close();
        delete f;
        return false;
    }

    TBranch *br_adc = tree->GetBranch("adc");
    TBranch *br_tot = tree->GetBranch("tot");
    if (!br_adc || !br_tot)
    {
        std::cerr << "Missing adc/tot branches in " << filename << "\n";
        f->Close();
        delete f;
        return false;
    }

    TLeaf *leaf_adc = br_adc->GetLeaf("adc");
    TLeaf *leaf_tot = br_tot->GetLeaf("tot");
    if (!leaf_adc || !leaf_tot)
    {
        std::cerr << "Missing adc/tot leaves in " << filename << "\n";
        f->Close();
        delete f;
        return false;
    }

    const int n_adc = leaf_adc->GetLen();
    const int n_tot = leaf_tot->GetLen();
    if (n_adc != n_tot || (n_adc % SAMPLES_PER_CHANNEL) != 0)
    {
        std::cerr << "Unexpected adc/tot shape in " << filename
                  << " n_adc=" << n_adc << " n_tot=" << n_tot
                  << " samples=" << SAMPLES_PER_CHANNEL << "\n";
        f->Close();
        delete f;
        return false;
    }

    const int num_channels = n_adc / SAMPLES_PER_CHANNEL;
    std::vector<uint32_t> adc_buf((size_t)n_adc);
    std::vector<uint32_t> tot_buf((size_t)n_tot);

    tree->SetBranchAddress("adc", adc_buf.data());
    tree->SetBranchAddress("tot", tot_buf.data());

    Long64_t nentries = tree->GetEntries();
    if (g_max_events > 0 && (Long64_t)g_max_events < nentries)
        nentries = (Long64_t)g_max_events;

    if (nentries <= 0)
    {
        std::cerr << "No entries in tree.\n";
        f->Close();
        delete f;
        return false;
    }

    long skipped_noncentral_tot = 0;

    for (Long64_t entry = 0; entry < nentries; ++entry)
    {
        if (entry * 25 / nentries > (Long64_t)(entry == 0 ? -1 : (entry - 1) * 25 / nentries))
        {
            int complete = (int)(entry * 25 / nentries);
            print_progress_25(complete);
        }

        tree->GetEntry(entry);

        // ----- Non-central ADC signals (reject if any non-central ToT) -----
        std::array<float, MAX_NUM_CRYSTALS> signals;
        signals.fill(0.0f);

        bool has_noncentral_tot = false;

        for (int cr = 0; cr < MAX_NUM_CRYSTALS; ++cr)
        {
            if (cr == central_crystal)
                continue;
            if (is_skipped_crystal(cr))
                continue;

            auto it = mapping.find(cr);
            if (it == mapping.end())
                continue;

            const std::vector<int> &chans = it->second;
            float crystal_signal = 0.0f;

            for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; ++sipm)
            {
                if (sipm >= (int)chans.size())
                    break;
                const int ch = chans[sipm];
                if (ch < 0 || ch >= num_channels)
                    continue;

                uint32_t *adc_vals = &adc_buf[(size_t)ch * SAMPLES_PER_CHANNEL];
                uint32_t *tot_vals = &tot_buf[(size_t)ch * SAMPLES_PER_CHANNEL];

                if (!has_noncentral_tot && has_tot(tot_vals))
                    has_noncentral_tot = true;

                const int idx = cr * SIPMS_PER_CRYSTAL + sipm;
                const float gch = (gain_factors ? (float)gain_factors->GetBinContent(idx + 1) : 1.0f);

                const float sig = calculate_signal_adc(adc_vals, gch);
                crystal_signal += sig;
            }

            const float gcr = (crystal_factor ? (float)crystal_factor->GetBinContent(cr + 1) : 1.0f);
            signals[cr] = crystal_signal * gcr;
        }

        if (has_noncentral_tot)
        {
            skipped_noncentral_tot++;
            continue;
        }

        // ----- Central crystal: require ToT on all SiPMs -----
        auto itc = mapping.find(central_crystal);
        if (itc == mapping.end())
            continue;

        const std::vector<int> &center_chans = itc->second;

        float tot_sum = 0.0f;
        int used = 0;

        for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; ++sipm)
        {
            if (sipm >= (int)center_chans.size())
                break;
            const int ch = center_chans[sipm];
            if (ch < 0 || ch >= num_channels)
                continue;

            uint32_t *tot_vals = &tot_buf[(size_t)ch * SAMPLES_PER_CHANNEL];

            if (!has_tot(tot_vals))
            {
                used = -999; // flag reject
                break;
            }

            const uint32_t raw_tot = get_tot_first(tot_vals);

            // optional width cuts per SiPM (based on histogram made by this macro)
            if (use_widths && widths_hist)
            {
                const float mean = (float)widths_hist->GetBinContent(sipm + 1);
                const float sigma = (float)widths_hist->GetBinError(sipm + 1);
                if (sigma > 0)
                {
                    if ((float)raw_tot < (mean - 2.0f * sigma) || (float)raw_tot > (mean + 2.0f * sigma))
                    {
                        used = -999;
                        break;
                    }
                }
            }

            // Fill per-SiPM raw ToT distribution (like Tristan)
            if (h_per_sipm_tot[sipm])
                h_per_sipm_tot[sipm]->Fill((float)raw_tot);

            // Sum with channel gain factor (match the old: center_signal += signal * gain)
            const int idx = central_crystal * SIPMS_PER_CRYSTAL + sipm;
            const float gch = (gain_factors ? (float)gain_factors->GetBinContent(idx + 1) : 1.0f);
            tot_sum += (float)raw_tot * gch;
            used++;
        }

        if (used != SIPMS_PER_CRYSTAL)
            continue;

        // normalize to 16 channels
        if (used > 0)
            tot_sum *= ((float)SIPMS_PER_CRYSTAL / (float)used);

        if (tot_sum < g_tot_min_cut)
            tot_sum = 0.0f;

        // COG cut (use ToT proxy for central crystal so ToT-only / ADC-empty runs still work)
        std::array<float, MAX_NUM_CRYSTALS> cog_w = signals; // noncentral ADC (maybe 0)
        cog_w[central_crystal] = tot_sum;                    // central ToT proxy
        // COG cut (also fills h_cog)
        if (!calculate_cog(h_cog, cog_w))
            continue;

        // total ToT distribution
        if (h_tot)
            h_tot->Fill(tot_sum);

        // Sum non-central energy (ADC-based) and convert to voltage-equivalent using adc_per_refunit
        float noncentral_adc = 0.0f;
        for (int cr = 0; cr < MAX_NUM_CRYSTALS; ++cr)
        {
            if (cr == central_crystal)
                continue;
            if (is_skipped_crystal(cr))
                continue;
            noncentral_adc += signals[cr];
        }

        // Convert ADC -> ref-units -> volts-equivalent
        float noncentral_volt = 0.0f;
        if (adc_per_refunit > 0.0f)
        {
            const float ref_units = noncentral_adc / adc_per_refunit; // 1 unit == ref_voltage
            noncentral_volt = ref_units * ref_voltage;
        }

        if (h_missing)
            h_missing->Fill(noncentral_volt);

        const float x = label_value - noncentral_volt;
        if (h_tot_vs_label_minus_missing)
            h_tot_vs_label_minus_missing->Fill(x, tot_sum);
    }

    std::cout << "\nSkipped events due to non-central ToT: " << skipped_noncentral_tot << "\n";

    f->Close();
    delete f;

    return true;
}

// -----------------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------------

void tot_calib_scan(const char *data_dir = "data",
                    const char *mapping_csv = "eeemcal_desy_dec2025_mapping_v2.csv",
                    const char *outdir = "outputs",
                    const char *gain_root = "outputs/gain_match_1.27V.root",
                    const char *adc_calib_root = "outputs/adc_to_ref_calibration_1.27V.root",
                    const char *widths_root = "outputs/tot_widths.root",
                    int central_crystal = 12)
{
    TH1::AddDirectory(kFALSE); // prevents histograms from attaching to gDirectory
    gROOT->cd();
    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(0);

    mkdir_p(outdir);

    // ------------------ Run list: ------------------
    // Format exactly like your adc_calibration scan list:
    std::vector<std::pair<int, float>> runs = {
        {48, 1.29f},
        {51, 1.30f},
        {54, 1.32f},
        {57, 1.33f},
        {60, 1.34f}};

    // ------------------ Mapping ------------------
    auto mapping = read_mapping_csv(mapping_csv, SIPMS_PER_CRYSTAL);
    if (mapping.empty())
    {
        std::cerr << "Error: mapping is empty. Check mapping_csv='" << mapping_csv << "'\n";
        return;
    }

    // ------------------ Calibration dependencies ------------------
    // adc_calibration.cxx writes mean_adc_to_ref_calibration_1.27V into outputs/adc_to_ref_calibration_1.27V.root
    const float adc_per_refunit = load_adc_to_ref(adc_calib_root, g_adc_ref_voltage_V);
    if (adc_per_refunit <= 0.0f)
    {
        std::cerr << "Error: could not load ADC->ref calibration from '" << adc_calib_root << "'.\n";
        std::cerr << "       Run adc_calibration first (it must write mean_adc_to_ref_calibration_*.\n";
        return;
    }
    std::cout << "Loaded ADC per ref-unit: " << adc_per_refunit << " (1 unit = " << g_adc_ref_voltage_V << " V)\n";

    TH1 *gain_factors = nullptr;
    TH1 *crystal_factor = nullptr;
    TFile *gain_file = nullptr;
    load_gain_factors(gain_root, gain_factors, crystal_factor, gain_file);
    std::cout << "Loaded gain factors from: " << gain_root << "\n";

    // ------------------ Width cuts file (optional) ------------------
    bool use_widths = false;
    TFile *widths_file = TFile::Open(widths_root, "READ");
    if (widths_file && !widths_file->IsZombie())
        use_widths = true;

    // ------------------ Global output histograms ------------------
    TH2F *total_dist = new TH2F("tot_vs_label_minus_missing",
                                "ToT vs (label - missing);label - missing (V);ToT sum (a.u.)",
                                120, 0, 2.5, 1200, 0, 60000);

    TH2F *total_invt = new TH2F("label_minus_missing_vs_tot",
                                "label-missing as f(ToT);ToT sum (a.u.);label - missing (V)",
                                1200, 0, 60000, 120, 0, 2.5);

    // Per-run containers (so we can also build widths if needed)
    struct RunOut
    {
        int run;
        float label;
        TH1F *missing;
        TH1F *tot;
        TH2F *tot2d;
        TH2F *cog;
        std::array<TH1F *, SIPMS_PER_CRYSTAL> per_sipm;

        // --- add these ---
        double tot_mu = 0.0;
        double tot_sigma = 0.0;
        double tot_res = 0.0; // sigma/mu
        bool tot_fit_ok = false;
    };

    std::vector<RunOut> outs;
    outs.reserve(runs.size());

    // ------------------ Process runs ------------------
    for (const auto &rp : runs)
    {
        const int runnum = rp.first;
        const float label = rp.second;

        std::string fn = std::string(Form("%s/Run%03d.root", data_dir, runnum));

        RunOut o;
        o.run = runnum;
        o.label = label;

        o.missing = new TH1F(Form("run%03d_missing_equivV", runnum),
                             Form("Run %03d non-central equiv (V);equiv V;Events", runnum),
                             120, 0, 2.5);

        o.tot = new TH1F(Form("run%03d_tot_sum", runnum),
                         Form("Run %03d ToT sum;ToT sum (a.u.);Events", runnum),
                         1200, 0, 60000);

        o.tot2d = new TH2F(Form("run%03d_tot_vs_label_minus_missing", runnum),
                           Form("Run %03d ToT vs (label-missing);label-missing (V);ToT sum (a.u.)", runnum),
                           120, 0, 2.5, 1200, 0, 60000);

        o.cog = new TH2F(Form("run%03d_cog", runnum),
                         Form("Run %03d COG;X (crystals);Y (crystals)", runnum),
                         100, -0.5, 4.5, 100, -0.5, 4.5);

        for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; ++sipm)
        {
            o.per_sipm[sipm] = new TH1F(Form("run%03d_sipm%02d_rawtot", runnum, sipm),
                                        Form("Run %03d SiPM %02d raw ToT;raw ToT;Events", runnum, sipm),
                                        1024, 0, 4096);
        }

        TH1 *widths_hist = nullptr;
        if (use_widths)
            widths_hist = (TH1 *)widths_file->Get(Form("run%d_tot_widths", runnum));

        std::cout << "\n=== Run " << runnum << " label=" << label << " ===\n";

        bool ok = tot_calibration_one(fn.c_str(), label,
                                      mapping, central_crystal,
                                      gain_factors, crystal_factor,
                                      adc_per_refunit, g_adc_ref_voltage_V,
                                      outdir,
                                      o.missing, o.tot, o.tot2d, o.cog,
                                      o.per_sipm,
                                      widths_hist,
                                      use_widths && (widths_hist != nullptr));

        if (!ok)
        {
            std::cerr << "Warning: failed processing " << fn << "\n";
        }

        // Accumulate into global
        total_dist->Add(o.tot2d);

        // Also fill inverted view (for fitting label as f(ToT) if you want)
        // This is just a re-binning copy:
        for (int ix = 1; ix <= o.tot2d->GetNbinsX(); ++ix)
        {
            for (int iy = 1; iy <= o.tot2d->GetNbinsY(); ++iy)
            {
                const double c = o.tot2d->GetBinContent(ix, iy);
                if (c <= 0)
                    continue;
                const double x = o.tot2d->GetXaxis()->GetBinCenter(ix);
                const double y = o.tot2d->GetYaxis()->GetBinCenter(iy);
                // total_invt has x=ToT, y=label-missing
                total_invt->Fill(y, x, c);
            }
        }

        // ---- Fit ToT-sum peak once per run and store resolution ----
        {
            TH1F *h = o.tot;
            // if (h && h->GetEntries() > 50 && h->GetRMS() > 0)
            if (h && h->GetRMS() > 0)
            {
                const double mean0 = h->GetBinCenter(h->GetMaximumBin());
                const double rms0 = std::max(50.0, h->GetRMS());
                const double fmin = std::max(0.0, mean0 - 1.5 * rms0);
                const double fmax = mean0 + 1.5 * rms0;

                TF1 fit(Form("tot_sum_gaus_run%03d", runnum), "gaus", fmin, fmax);
                fit.SetParameters(h->GetMaximum(), mean0, rms0);

                TFitResultPtr rr = h->Fit(&fit, "RQS"); // stores fit on the hist

                if (rr.Get() && rr->Status() == 0)
                {
                    o.tot_mu = fit.GetParameter(1);
                    o.tot_sigma = std::abs(fit.GetParameter(2));
                }
                else
                {
                    o.tot_mu = h->GetMean();
                    o.tot_sigma = h->GetRMS();
                }
                o.tot_res = (o.tot_mu != 0.0) ? (o.tot_sigma / o.tot_mu) : 0.0;
                o.tot_fit_ok = true;

                std::cout << "Run " << runnum
                          << " ToT sum peak: mean=" << o.tot_mu
                          << " sigma=" << o.tot_sigma
                          << " resolution=" << 100.0 * o.tot_res << " %\n";
            }
        }

        outs.push_back(o);
    }

    if (widths_file)
    {
        std::cout << "\nClosed widths file: " << widths_root << "\n";
        widths_file->Close();
        delete widths_file;
    }

    // ------------------ PDF summary ------------------
    std::cout << "\nInit PDF summary to " << outdir << "\n";
    std::string pdf = std::string(Form("%s/tot_calib.pdf", outdir));
    TCanvas *canvas = new TCanvas("tot_calib", "", 900, 700);
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.04);
    text.SetTextFont(42);

    canvas->SaveAs((pdf + "[").c_str());

    // store peak mean/sigma per run to optionally write widths file
    std::vector<std::vector<double>> peak_mean(outs.size(), std::vector<double>(SIPMS_PER_CRYSTAL, 0.0));
    std::vector<std::vector<double>> peak_sigma(outs.size(), std::vector<double>(SIPMS_PER_CRYSTAL, 0.0));

    for (size_t ir = 0; ir < outs.size(); ++ir)
    {
        std::cout << "\nPreparing summary plots for run " << outs[ir].run << "...\n";
        auto &o = outs[ir];

        canvas->Clear();
        canvas->Divide(3, 2);

        canvas->cd(1);
        o.missing->Draw();
        text.SetTextAlign(31);
        text.DrawLatex(0.88, 0.85, Form("Mean: %.3f V", o.missing->GetMean()));
        std::cout << "Run " << o.run << " missing (non-central) mean (V): " << o.missing->GetMean() << "\n";
        text.DrawLatex(0.88, 0.80, Form("Label: %.2f", o.label));

        canvas->cd(2);
        o.tot->Draw();

        // Draw the fit already done in the run loop
        if (TF1 *f = o.tot->GetFunction(Form("tot_sum_gaus_run%03d", o.run)))
            f->Draw("same");

        canvas->cd(3);
        o.tot2d->Draw("colz");

        canvas->cd(4);
        o.cog->Draw("colz");

        // per-sipm fits
        auto pad = canvas->cd(5);
        pad->Divide(4, 4, 0.002, 0.002);
        for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; ++sipm)
        {
            pad->cd(sipm + 1);
            TH1F *h = o.per_sipm[sipm];
            h->SetTitle("");
            h->Draw("hist");

            // Gaussian fit around max bin (same spirit as old macro)
            const float mean0 = h->GetBinCenter(h->GetMaximumBin());
            const float rms0 = std::max(50.0f, (float)h->GetRMS());
            const float fmin = std::max(0.0f, mean0 - 1.5f * rms0);
            const float fmax = mean0 + 1.5f * rms0;

            TF1 fit(Form("tot_fit_run%03d_sipm%02d", o.run, sipm), "gaus", fmin, fmax);
            fit.SetParameter(0, h->GetMaximum());
            fit.SetParameter(1, mean0);
            fit.SetParameter(2, rms0);

            TFitResultPtr r = h->Fit(&fit, "RQS");
            if (!r.Get() || r->Status() != 0 || fit.GetNDF() < 20)
            {
                // fallback
                TF1 fit2(Form("tot_fit2_run%03d_sipm%02d", o.run, sipm), "gaus", std::max(0.0f, mean0 - 300.0f), mean0 + 300.0f);
                fit2.SetParameter(0, h->GetMaximum());
                fit2.SetParameter(1, mean0);
                fit2.SetParameter(2, 200.0f);
                h->Fit(&fit2, "RQS");
                peak_mean[ir][sipm] = fit2.GetParameter(1);
                peak_sigma[ir][sipm] = fit2.GetParameter(2);
                fit2.Draw("same");
            }
            else
            {
                peak_mean[ir][sipm] = fit.GetParameter(1);
                peak_sigma[ir][sipm] = fit.GetParameter(2);
                fit.Draw("same");
            }
        }

        canvas->cd(6);
        text.SetTextAlign(13);
        text.DrawLatex(0.12, 0.85, Form("Run %03d", o.run));
        text.DrawLatex(0.12, 0.78, Form("Label value: %.2f (treated as V in plots)", o.label));
        text.DrawLatex(0.12, 0.71, Form("ADC->ref: %.1f ADC per %.2fV", adc_per_refunit, g_adc_ref_voltage_V));
        text.DrawLatex(0.12, 0.64, Form("COG frac9 cut: %.2f; ellipse: %s", g_energy_fraction_cut, g_do_cog_ellipse_cut ? "ON" : "OFF"));
        text.DrawLatex(0.12, 0.57,
                       Form("ToT sum: mean=%.1f  sigma=%.1f  res=%.2f%%",
                            o.tot_mu, o.tot_sigma, 100.0 * o.tot_res));

        canvas->SaveAs(pdf.c_str());
    }

    // global page + fits
    canvas->Clear();
    total_dist->Draw("colz");
    canvas->SaveAs(pdf.c_str());

    // Fits of label as a function of ToT (optional, keep same idea as Tristan's pol1/pol2)
    // We fit on total_invt (x=ToT, y=label-missing)
    canvas->Clear();
    total_invt->Draw("colz");

    float scale_factor = 16.0f / (float)SIPMS_PER_CRYSTAL;

    TF1 *fit_hi = new TF1("fit_hi", "pol1", 35000.0f / scale_factor, 50000.0f / scale_factor);
    total_invt->Fit(fit_hi, "RQ");
    fit_hi->SetLineColor(kRed);

    TF1 *fit_mid = new TF1("fit_mid", "pol1", 20000.0f / scale_factor, 26000.0f / scale_factor);
    total_invt->Fit(fit_mid, "RQ");
    fit_mid->SetLineColor(kGreen + 2);

    TF1 *fit_all = new TF1("fit_all", "pol2", 20000.0f / scale_factor, 55000.0f / scale_factor);
    total_invt->Fit(fit_all, "RQ");
    fit_all->SetLineColor(kMagenta);

    fit_hi->Draw("same");
    fit_mid->Draw("same");
    fit_all->Draw("same");

    canvas->SaveAs(pdf.c_str());
    canvas->SaveAs((pdf + "]").c_str());

    // ------------------ ROOT outputs ------------------
    std::string root_out = std::string(Form("%s/tot_calibration_values.root", outdir));
    TFile *out = TFile::Open(root_out.c_str(), "RECREATE");

    // store fit params (same names as old macro)
    TParameter<float> tot_c0("tot_c0", (float)fit_all->GetParameter(0));
    TParameter<float> tot_c1("tot_c1", (float)fit_all->GetParameter(1));
    TParameter<float> tot_c2("tot_c2", (float)fit_all->GetParameter(2));
    TParameter<float> tot_a0("tot_a0", (float)fit_mid->GetParameter(0));
    TParameter<float> tot_a1("tot_a1", (float)fit_mid->GetParameter(1));

    tot_c0.Write();
    tot_c1.Write();
    tot_c2.Write();
    tot_a0.Write();
    tot_a1.Write();

    // store config
    TParameter<float> adc_per_unit("adc_per_refunit", adc_per_refunit);
    TParameter<float> refV("adc_ref_voltage_V", g_adc_ref_voltage_V);
    adc_per_unit.Write();
    refV.Write();

    total_dist->Write();
    total_invt->Write();

    for (auto &o : outs)
    {
        o.missing->Write();

        {
            // ---- Save per-run ToT-sum resolution parameters ----
            TParameter<float> p_mu(Form("run%03d_tot_sum_mu", o.run), (float)o.tot_mu);
            TParameter<float> p_si(Form("run%03d_tot_sum_sigma", o.run), (float)o.tot_sigma);
            TParameter<float> p_re(Form("run%03d_tot_sum_res", o.run), (float)o.tot_res);
            TParameter<float> p_re_pct(Form("run%03d_tot_sum_res_pct", o.run), (float)(100.0 * o.tot_res));

            p_mu.Write();
            p_si.Write();
            p_re.Write();
            p_re_pct.Write();
        }

        o.tot->Write();
        o.tot2d->Write();
        o.cog->Write();
        for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; ++sipm)
            if (o.per_sipm[sipm])
                o.per_sipm[sipm]->Write();
    }

    out->Close();
    delete out;

    std::cout << "Saved: " << pdf << "\n";
    std::cout << "Saved: " << root_out << "\n";

    // ------------------ Create widths file if missing ------------------
    if (!use_widths)
    {
        std::string wout = std::string(Form("%s/tot_widths.root", outdir));
        TFile *wf = TFile::Open(wout.c_str(), "RECREATE");
        for (size_t ir = 0; ir < outs.size(); ++ir)
        {
            const int runnum = outs[ir].run;
            TH1F *hw = new TH1F(Form("run%d_tot_widths", runnum),
                                Form("Run %d ToT widths;SiPM;raw ToT", runnum),
                                SIPMS_PER_CRYSTAL, 0, SIPMS_PER_CRYSTAL);
            for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; ++sipm)
            {
                hw->SetBinContent(sipm + 1, (float)peak_mean[ir][sipm]);
                hw->SetBinError(sipm + 1, (float)peak_sigma[ir][sipm]);
            }
            hw->Write();
        }
        wf->Close();
        delete wf;
        std::cout << "Widths file did not exist; wrote: " << wout << "\n";
    }

    // cleanup gain file handle
    if (gain_file)
    {
        gain_file->Close();
        delete gain_file;
    }

    delete canvas;
}

// float energy_fraction_cut = 0.3;
// long n_events = 1000000;
// // n_events = 100;

// float tot_min_cut = 0000;

// std::vector<int> run_numbers = {330};
// std::vector<float> energies  = {3.8};

// void tot_calibration() {
//     gErrorIgnoreLevel = kWarning;

//     run_numbers.clear();
//     energies.clear();
//     for (int r = 319; r <= 338; r++) {
//         run_numbers.push_back(r);
//         energies.push_back(1.6 + 0.2 * (r - 319));
//     }

//     gStyle->SetOptStat(0);
//     auto mapping = read_mapping("eeemcal_desy_dec2025_mapping.csv");

//     float adc_calib = 0;
//     TFile *adc_calib_file = TFile::Open("output/adc_to_gev_calibration.root", "READ");
//     if (adc_calib_file && !adc_calib_file->IsZombie()) {
//         TParameter<float>* adc_calib_param = (TParameter<float>*)adc_calib_file->Get("mean_adc_to_gev_calibration");
//         if (adc_calib_param) {
//             adc_calib = adc_calib_param->GetVal();
//             std::cout << "Loaded ADC to GeV calibration: " << adc_calib << std::endl;
//         }
//     }
//     if (adc_calib == 0) {
//         std::cerr << "Error: ADC to GeV calibration not found!" << std::endl;
//         return;
//     }
//     adc_calib_file->Close();

//     TH1* gain_factor = nullptr;
//     TH1* crystal_gain_factor;
//     TFile* gain_file = TFile::Open("output/gain_factors.root");
//     if (gain_file && !gain_file->IsZombie()) {
//         gain_factor = (TH1*)gain_file->Get("gain_factors");
//         crystal_gain_factor = (TH1*)gain_file->Get("crystal_factor");
//         std::cout << "Loaded gain factors from file." << std::endl;
//     }
//     if (!gain_factor) {
//         gain_factor = new TH1F("gain_factors", "Gain Factor", 400, 0, 400);
//         for (int i = 1; i <= 400; i++) {
//             gain_factor->SetBinContent(i, 1.0);
//         }
//     }
//     if (!crystal_gain_factor) {
//         crystal_gain_factor = new TH1F("crystal_factor", "Crystal Gain", 25, 0, 25);
//         for (int i = 1; i <= 25; i++) {
//             crystal_gain_factor->SetBinContent(i, 1);
//         }
//     }

//     std::vector<TH1*> missing_energies;
//     std::vector<TH1*> tot_distributions;
//     std::vector<TH1*> tot_energies;
//     std::vector<TH2*> cog_distributions;
//     std::vector<TH1*> tot_bin_distribution;
//     std::vector<TH2*> tot_bin_vs_value;
//     std::vector<std::vector<TH1*>> per_sipm_tot;
//     std::vector<std::vector<std::vector<TH2*>>> tot_covariance;

//     TH2 *total_distribution =new TH2F("run%d_tot_energy", "ToT Energy;Energy (GeV);ToT Signal", 100, 0, 5, 1000, 0, 60000);
//     TH2 *total_distribution_invt =new TH2F("run%d_tot_energy_invt", "Energy as a function of ToT;ToT Signal;Energy (GeV);", 1000, 0, 60000, 100, 0, 5);

//     bool use_widths = false;
//     TFile *input_file = TFile::Open("output/tot_widths.root", "READ");
//     if (!input_file || input_file->IsZombie()) {
//         std::cerr << "Error: Could not open input file 'output/tot_width.root'" << std::endl;
//     } else {
//         use_widths = true;
//     }

//     for (int run = 0; run < run_numbers.size(); run++) {
//         TH1 *peak_widths = nullptr;
//         if (use_widths) {
//             peak_widths = (TH1*)input_file->Get(Form("run%d_tot_widths", run_numbers[run]));
//             if (!peak_widths) {
//                 std::cerr << "Warning: Could not find histogram 'run" << run_numbers[run] << "_tot_widths' in input file." << std::endl;
//                 use_widths = false;
//             }
//             else {
//                 std::cout << "Using ToT widths for run " << run_numbers[run] << std::endl;
//             }
//         }

//         int run_number = run_numbers[run];
//         float energy = energies[run];
//         // Process data file
//         TFile* root_file = TFile::Open(Form("/Users/tristan/dropbox/eeemcal_desy_dec_2025/prod_0/Run%03d.root", run_number));
//         TTree* tree = (TTree*)root_file->Get("events");
//         uint32_t adc[576][20];
//         uint32_t tot[576][20];
//         int tot_events = 0;
//         tree->SetBranchAddress("adc", &adc);
//         tree->SetBranchAddress("tot", &tot);
//         tree->SetBranchStatus("*", 0);
//         tree->SetBranchStatus("adc", 1);
//         tree->SetBranchStatus("tot", 1);

//         missing_energies.push_back(new TH1F(Form("run%d_missing_energy", run_number), Form("Run %d ADC portion of energy", run_number), 100, 0, 6));
//         tot_distributions.push_back(new TH1F(Form("run%d_tot_distribution", run_number), Form("Run %d ToT;ToT;Events", run_number), 1024, 0, 16 * 4096));
//         tot_energies.push_back(new TH2F(Form("run%d_tot_energy", run_number), Form("Run %d ToT Energy;Energy (GeV);(ToT Signal)", run_number), 100, 0, 5, 1000, 0, 60000));
//         cog_distributions.push_back(new TH2F(Form("run%d_cog_distribution", run_number), Form("Run %d Center of Gravity Distribution;X (# Crystals);Y (# Crystals)", run_number), 100, -0.5, 4.5, 100, -0.5, 4.5));
//         tot_bin_distribution.push_back(new TH1F(Form("run%d_tot_bin_distribution", run_number), Form("Run %d ToT Bin Distribution;ToT Bin;Events", run_number), 20, 0, 20));
//         tot_bin_vs_value.push_back(new TH2F(Form("run%d_tot_bin_vs_value", run_number), Form("Run %d ToT Bin vs ToT Value;ToT Value;ToT Bin", run_number), 20, 0, 20, 1024, 0, 4096));

//         per_sipm_tot.push_back(std::vector<TH1*>());
//         tot_covariance.push_back(std::vector<std::vector<TH2*>>());
//         for (int i = 0; i < sipms_to_use; i++) {
//             per_sipm_tot[run].push_back(new TH1F(Form("run%d_sipm%d_tot_distribution", run_number, i), Form("Run %d SiPM %d ToT;ToT;Events", run_number, i), 1024, 0, 4096));
//             tot_covariance[run].push_back(std::vector<TH2*>());
//             for (int j = 0; j < sipms_to_use; j++) {
//                 tot_covariance[run][i].push_back(new TH2F(Form("run%d_sipm%d_sipm%d_tot_covariance", run_number, i, j),
//                                                         Form("Run %d SiPM %d vs SiPM %d ToT;SiPM %d ToT;SiPM %d ToT", run_number, i, j, i, j),
//                                                         1024, 0, 4096, 1024, 0, 4096));
//             }
//         }

//         Long64_t nentries = tree->GetEntries();
//         if (n_events < nentries) {
//             nentries = n_events;
//         }
//         std::cout << "Run " << run_number << std::endl;
//         std::cout << "Processing " << nentries << " events" << std::endl;
//         int complete = 0;
//         for (Long64_t entry = 0; entry < nentries; ++entry) {
//             if (entry * 25 / nentries > complete) {
//                 complete = entry * 25 / nentries;
//                 print_progress(complete);
//             }

//             tree->GetEntry(entry);
//             float central_signal = 0.0f;
//             float central_nine_signal = 0.0f;
//             float total_signal = 0.0f;
//             bool is_tot_event = false;

//             float signals[25];
//             for (int crystal = 0; crystal < 25; crystal++) {
//                 if (crystal == 9) {
//                     continue;
//                 }
//                 if (crystal == 12) {
//                     continue;
//                 }
//                 float crystal_signal = 0.0f;
//                 for (int sipm = 0; sipm < sipms_to_use; sipm++) {
//                     int channel = mapping[crystal][sipm];
//                     float gain = gain_factor->GetBinContent(crystal * 16 + sipm + 1);
//                     float channel_signal = calculate_signal(adc[channel], gain);
//                     crystal_signal += channel_signal;
//                     if (!is_tot_event && is_tot(tot[channel])) {
//                         is_tot_event = true;
//                     }
//                 }
//                 signals[crystal] = crystal_signal * crystal_gain_factor->GetBinContent(crystal + 1);
//             }
//             if (is_tot_event) {
//                 tot_events++;
//                 continue;
//             }

//             // Get the ToT
//             float center_signal = 0;
//             std::vector<float> sipm_signals(16, 0);
//             int channels_used = 0;
//             bool all_tot = true;
//             for (int sipm = 0; sipm < sipms_to_use; sipm++) {
//                 float signal = 0;
//                 int channel = mapping[12][sipm];
//                 float gain = gain_factor->GetBinContent(12 * 16 + sipm + 1);
//                 all_tot &= calculate_signal(adc[channel], tot[channel], gain, signal);
//                 if (!all_tot) {
//                     break;
//                 }
//                 if (use_widths) {
//                     float channel_mean = peak_widths->GetBinContent(sipm + 1);
//                     float channel_sigma = peak_widths->GetBinError(sipm + 1);
//                     // std::cout << "signal: " << signal << ", mean: " << channel_mean << ", sigma: " << channel_sigma << std::endl;
//                     if (signal < channel_mean - (2 * channel_sigma) || signal > channel_mean + (2 * channel_sigma)) {
//                         all_tot = false;
//                         break;
//                         continue;
//                     }
//                 }
//                 sipm_signals[sipm] = signal;
//                     // tot_bin_distribution[run]->Fill(timebin);
//                     // tot_bin_vs_value[run]->Fill(timebin, signal);
//                 per_sipm_tot[run][sipm]->Fill(signal);
//                 center_signal += signal * gain;
//                 channels_used++;
//             }
//             if (!all_tot || channels_used != sipms_to_use) {
//                 continue;
//             }
//             // normalize the center signal to 16 channels used
//             if (channels_used > 0) {
//                 center_signal *= ((float)sipms_to_use / channels_used);
//             }

//             if (center_signal < tot_min_cut) {
//                 center_signal = 0;
//             }
//             // Fill the covariance matrices
//             for (int i = 0; i < sipms_to_use; i++) {
//                 for (int j = 0; j < sipms_to_use; j++) {
//                     tot_covariance[run][i][j]->Fill(sipm_signals[i], sipm_signals[j]);
//                 }
//             }

//             float x_cog, y_cog;
//             bool keep = calculate_cog(cog_distributions[run], signals);
//             keep &= (!is_tot_event);    // Get rid of events with a tot in a non-central channel
//             if (!keep) {
//                 continue;
//             }

//             tot_distributions[run]->Fill(center_signal);

//             // Sum up the total energy found in the non-central crystals
//             float non_central_energy = 0;
//             for (int crystal = 0; crystal < 25; crystal++) {
//                 if (crystal == 12) {continue;}
//                 non_central_energy += signals[crystal];
//             }
//             // Calibrate to GeV
//             non_central_energy /= adc_calib;
//             missing_energies[run]->Fill(non_central_energy);
//             tot_energies[run]->Fill(energy - non_central_energy, center_signal);
//             total_distribution->Fill(energy - non_central_energy, center_signal);
//             total_distribution_invt->Fill(center_signal, energy - non_central_energy);
//         }
//         // std::cout << std::endl;
//         std::cout << "\rTotal ToT events skipped: " << tot_events << std::endl;
//     }

//     if (use_widths) {
//         input_file->Close();
//     }

//     TLatex *text = new TLatex();
//     text->SetNDC();
//     text->SetTextSize(0.04);
//     text->SetTextFont(42);

//     std::vector<std::vector<double>> peak_mean;
//     std::vector<std::vector<double>> peak_sigma;

//     TCanvas* canvas = new TCanvas("tot_calib", "", 800, 600);
//     canvas->SaveAs("output/tot_calib.pdf(");
//     for (int run = 0; run < run_numbers.size(); run++) {
//         peak_mean.push_back(std::vector<double>());
//         peak_sigma.push_back(std::vector<double>());
//         canvas->Clear();
//         canvas->Divide(3, 2);
//         canvas->cd(1);
//         missing_energies[run]->Draw();
//         text->SetTextAlign(31);
//         text->DrawLatex(0.85, 0.85, Form("Mean: %.2f", missing_energies[run]->GetMean()));
//         text->DrawLatex(0.85, 0.80, Form("ADC percentage: %.2f", 100 * missing_energies[run]->GetMean() / energies[run]));

//         canvas->cd(2);
//         tot_distributions[run]->Draw();
//         canvas->cd(3);
//         tot_energies[run]->Draw("colz");
//         canvas->cd(4);
//         cog_distributions[run]->Draw("colz");

//         auto pad = canvas->cd(5);
//         pad->Divide(4, 4, 0, 0);
//         for (int sipm = 0; sipm < sipms_to_use; sipm++) {
//             pad->cd(sipm + 1);
//             // pad->SetMargin(0, 0, 0, 0);
//             per_sipm_tot[run][sipm]->SetTitle("");
//             float mean = per_sipm_tot[run][sipm]->GetBinCenter(per_sipm_tot[run][sipm]->GetMaximumBin());
//             float rms = 500;
//             float fit_min = mean - 1.5 * rms;
//             float fit_max = mean + 1.5 * rms;
//             auto fit = new TF1("tot_fit", "gaus", fit_min, fit_max);
//             fit->SetParameter(0, per_sipm_tot[run][sipm]->GetBinContent(per_sipm_tot[run][sipm]->GetMaximumBin()));
//             fit->SetParameter(1, mean);
//             fit->SetParameter(2, rms);
//             auto result = per_sipm_tot[run][sipm]->Fit(fit, "RQS");
//             // check if the fit failed
//             if (result->Status() != 0 || fit->GetNDF() < 20) {
//                 std::cerr << "Warning: Fit failed for run " << run_numbers[run] << " SiPM " << sipm << std::endl;
//                 std::cout << "       Using default fit parameters." << std::endl;
//                 delete fit;
//                 fit = new TF1("tot_fit", "gaus", 2000, 3000);
//                 fit->SetParameter(0, 200);
//                 fit->SetParameter(1, 2200);
//                 fit->SetParameter(2, 500);
//                 per_sipm_tot[run][sipm]->Fit(fit, "RQS");
//             }
//             per_sipm_tot[run][sipm]->Draw("hist");
//             fit->Draw("same");
//             peak_mean[run].push_back(fit->GetParameter(1));
//             peak_sigma[run].push_back(fit->GetParameter(2));
//         }

//         canvas->cd(6);
//         tot_bin_vs_value[run]->Draw("colz");

//         canvas->SaveAs("output/tot_calib.pdf");

//         bool do_covariance = false;
//         if (do_covariance) {
//             canvas->Clear();
//             canvas->Divide(16, 16, 0, 0);
//             for (int i = 0; i < sipms_to_use; i++) {
//                 for (int j = 0; j < sipms_to_use; j++) {
//                     if (j <= i) {
//                         canvas->cd(i * 16 + j + 1);
//                         tot_covariance[run][i][j]->Draw("colz");
//                     }
//                 }
//             }
//             canvas->SaveAs("output/tot_calib.pdf");
//         }
//     }
//     canvas->Clear();
//     total_distribution->Draw("colz");
//     canvas->SaveAs("output/tot_calib.pdf");

//     float scale_factor = 16.0f / sipms_to_use;

//     TF1 *range_one_fit = new TF1("range_one", "pol1", 35000 / scale_factor, 50000 / scale_factor);
//     total_distribution_invt->Fit(range_one_fit, "RQ");
//     range_one_fit->SetLineColor(kRed);

//     TF1 *range_two_fit = new TF1("range_two", "pol1", 20000 / scale_factor, 26000 / scale_factor);
//     total_distribution_invt->Fit(range_two_fit, "RQ");
//     range_two_fit->SetLineColor(kGreen + 2);

//     TF1 *total_fit = new TF1("total_fit", "pol2", 20000 / scale_factor, 55000 / scale_factor);
//     total_distribution_invt->Fit(total_fit, "RQ");
//     total_fit->SetLineColor(kMagenta);

//     total_distribution_invt->Draw("colz");
//     range_one_fit->Draw("same");
//     range_two_fit->Draw("same");
//     total_fit->Draw("same");
//     canvas->SaveAs("output/tot_calib.pdf)");

//     TFile *tot_calib_file = TFile::Open("output/tot_calibration_values.root", "RECREATE");
//     TParameter<float>* tot_c0 = new TParameter<float>("tot_c0", total_fit->GetParameter(0));
//     TParameter<float>* tot_c1 = new TParameter<float>("tot_c1", total_fit->GetParameter(1));
//     TParameter<float>* tot_c2 = new TParameter<float>("tot_c2", total_fit->GetParameter(2));
//     TParameter<float>* tot_a0 = new TParameter<float>("tot_a0", range_two_fit->GetParameter(0));
//     TParameter<float>* tot_a1 = new TParameter<float>("tot_a1", range_two_fit->GetParameter(1));
//     tot_c0->Write();
//     tot_c1->Write();
//     tot_c2->Write();
//     tot_a0->Write();
//     tot_a1->Write();
//     tot_calib_file->Close();

//     if (!use_widths) {
//         TFile *tot_widths = TFile::Open("output/tot_widths.root", "RECREATE");
//         for (int run = 0; run < run_numbers.size(); run++) {
//             TH1* tot_width = new TH1F(Form("run%d_tot_widths", run_numbers[run]), Form("Run %d ToT Widths;SiPM;Width (ADC counts)", run_numbers[run]), 16, 0, 16);
//             for (int sipm = 0; sipm < sipms_to_use; sipm++) {
//                 float mean = peak_mean[run][sipm];
//                 float sigma = peak_sigma[run][sipm];
//                 tot_width->SetBinContent(sipm + 1, mean);
//                 tot_width->SetBinError(sipm + 1, sigma);
//             }
//             tot_width->Write();
//         }
//     }
// }