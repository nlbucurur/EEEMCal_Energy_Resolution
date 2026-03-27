// tot_calibration_v2.cxx
//
// Refactor of Tristan's tot_calibration.C to match the event/mapping style used in
// led_analysis.C / draw_waveform_reversed.cxx / gain_match.cxx / adc_calibration.cxx:
//  - Uses common_led.{h,cxx} mapping + signal helpers
//  - Reads adc/tot as flat buffers using TLeaf::GetLen() (dynamic sizes)
//  - Uses mapping CSV produced by generate_mapping.py
//  - Keeps the original purpose: calibrate ToT response vs "(ref - missing energy)" for ToT events
//
// IMPORTANT about the x-axis:
//  - In the original Tristan code this was (beam_energy[GeV] - noncentral_energy[GeV]).
//  - Here we generalize it to (label_value - noncentral_equivalent_value).
//  - If you pass label_value as LED voltage (V), we convert the non-central ADC sum into a
//    voltage-equivalent using adc_calibration output (ADC per 1.27V or 1.25V reference unit).
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

#include <TString.h>
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
#include <set>

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
#include <TROOT.h>
#include <TGraphErrors.h>
#include <TProfile.h>
//========================================================//
//==============       Change this        ================//
//========================================================//
// ADC calibration reference voltage used by adc_calibration.cxx
static float g_adc_ref_voltage_V = 1.259f;
static bool g_scale_subset_to_16 = false; // if true, then if all 16 central channels are chosen, we scale their sum to 16 before conversion, then scale back after conversion. This mimics the behavior of adc_calibration.cxx where we always scale to 16 channels before conversion.
//========================================================//
//==============       Change this        ================//
//========================================================//

// ----- knobs (keep same defaults as Tristan's original) -----
static float g_energy_fraction_cut = 0.30f; // used in COG selection
static float g_tot_min_cut = 0.0f; // if ToT sum < this -> set to 0

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

// Compute and optionally cut on COG, and optionally cut on central9 fraction.
static bool calculate_cog_tot(TH2 *distribution,
                          const std::array<float, LED_MAX_NUM_CRYSTALS> &signals)
{
    // total
    double sumE = 0.0;
    double sumX = 0.0;
    double sumY = 0.0;

    // fraction in central 3x3 around crystal 12
    const int central9[9] = {6, 7, 8, 11, 12, 13, 16, 17, 18};
    double sumE9 = 0.0;

    for (int cr = 0; cr < LED_MAX_NUM_CRYSTALS; ++cr)
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
        if (cr < 0 || cr >= LED_MAX_NUM_CRYSTALS)
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
    std::string exact = std::string(Form("mean_adc_to_ref_calibration_%.3fV", ref_voltage));
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

static bool compute_resolution_pct(const TH1 *h,
                                   double &mean, double &mean_err,
                                   double &sigma, double &sigma_err,
                                   double &reso_pct, double &reso_pct_err)
{
    if (!h)
        return false;
    if (h->GetEntries() < 20)
        return false;

    mean = h->GetMean();
    sigma = h->GetRMS();
    mean_err = h->GetMeanError();
    sigma_err = h->GetRMSError();

    if (mean <= 0.0 || sigma <= 0.0)
        return false;

    reso_pct = 100.0 * sigma / mean;

    double rel2 = 0.0;
    if (sigma_err > 0.0)
        rel2 += (sigma_err / sigma) * (sigma_err / sigma);
    if (mean_err > 0.0)
        rel2 += (mean_err / mean) * (mean_err / mean);

    reso_pct_err = reso_pct * std::sqrt(rel2);
    return true;
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
                                std::array<TH1F *, LED_SIPMS_PER_CRYSTAL> &h_per_sipm_tot,
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

    tree->GetEntry(0);
    const int n_adc = leaf_adc->GetLen();
    const int n_tot = leaf_tot->GetLen();
    if (n_adc != n_tot || (n_adc % LED_SAMPLES_PER_CHANNEL) != 0)
    {
        std::cerr << "Unexpected adc/tot shape in " << filename
                  << " n_adc=" << n_adc << " n_tot=" << n_tot
                  << " samples=" << LED_SAMPLES_PER_CHANNEL << "\n";
        f->Close();
        delete f;
        return false;
    }

    const int num_channels = n_adc / LED_SAMPLES_PER_CHANNEL;
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
        std::array<float, LED_MAX_NUM_CRYSTALS> signals;
        signals.fill(0.0f);

        bool has_noncentral_tot = false;

        for (int cr = 0; cr < LED_MAX_NUM_CRYSTALS; ++cr)
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

            for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
            {
                if (sipm >= (int)chans.size())
                    break;
                const int ch = chans[sipm];
                if (ch < 0 || ch >= num_channels)
                    continue;

                uint32_t *adc_vals = &adc_buf[(size_t)ch * LED_SAMPLES_PER_CHANNEL];
                uint32_t *tot_vals = &tot_buf[(size_t)ch * LED_SAMPLES_PER_CHANNEL];

                if (!has_noncentral_tot && has_tot(tot_vals))
                    has_noncentral_tot = true;

                const int idx = cr * LED_SIPMS_PER_CRYSTAL + sipm;
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

        // ----- Central crystal: require ToT on all selected SiPMs -----
        auto itc = mapping.find(central_crystal);
        if (itc == mapping.end())
            continue;

        const std::vector<int> &center_chans = itc->second;

        float tot_sum = 0.0f;
        int used = 0;
        const int n_expected =
            count_selected_sipms_in_mapping(mapping, central_crystal, LED_SIPMS_PER_CRYSTAL);

        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
        {
            if (!use_selected_sipm(central_crystal, sipm))
                continue;

            if (sipm >= (int)center_chans.size())
                break;

            const int ch = center_chans[sipm];

            if (ch < 0 || ch >= num_channels)
                continue;

            uint32_t *tot_vals = &tot_buf[(size_t)ch * LED_SAMPLES_PER_CHANNEL];

            if (!has_tot(tot_vals))
            {
                used = -999; // flag reject
                break;
            }

            // const uint32_t raw_tot = get_tot_first(tot_vals);
            const uint32_t raw_tot = get_tot_max(tot_vals);

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
            const int idx = central_crystal * LED_SIPMS_PER_CRYSTAL + sipm;
            const float gch = (gain_factors ? (float)gain_factors->GetBinContent(idx + 1) : 1.0f);
            tot_sum += (float)raw_tot * gch;
            used++;
        }

        if (used != n_expected)
            continue;

        // if (used != LED_SIPMS_PER_CRYSTAL)
        //     continue;

        // normalize to 16 channels
        if (g_scale_subset_to_16 && used > 0)
            tot_sum *= ((float)LED_SIPMS_PER_CRYSTAL / (float)used);

        if (tot_sum < g_tot_min_cut)
            tot_sum = 0.0f;

        // COG cut (use ToT proxy for central crystal so ToT-only / ADC-empty runs still work)
        std::array<float, LED_MAX_NUM_CRYSTALS> cog_w = signals; // noncentral ADC (maybe 0)
        cog_w[central_crystal] = tot_sum;                        // central ToT proxy
        // COG cut (also fills h_cog)
        if (!calculate_cog_tot(h_cog, cog_w))
            continue;

        // total ToT distribution
        if (h_tot)
            h_tot->Fill(tot_sum);

        // Sum non-central energy (ADC-based) and convert to voltage-equivalent using adc_per_refunit
        float noncentral_adc = 0.0f;
        for (int cr = 0; cr < LED_MAX_NUM_CRYSTALS; ++cr)
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
                    const char *gain_root = "outputs/gain_match_1.259V.root",
                    const char *adc_calib_root = "outputs/adc_to_ref_calibration_1.259V.root",
                    const char *widths_root = "outputs/tot_widths.root",
                    int central_crystal = 12)
{
    TH1::AddDirectory(kFALSE); // prevents histograms from attaching to gDirectory
    gROOT->cd();
    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(0);
    // gStyle->SetPadGridX(true);
    // gStyle->SetPadGridY(true);

    mkdir_p(outdir);

    // ------------------ Run list: ------------------
    //========================================================//
    //==============       Change this        ================//
    //========================================================//

    // Format exactly like your adc_calibration scan list:
    // (run, voltage)
    std::vector<std::pair<int, float>> runs = {
        // {23, 0.0f},
        // {26, 1.2f},
        // {30, 1.22f},
        // {33, 1.24f},
        // {36, 1.25},
        // {39, 1.26f},
        // {42, 1.27f},
        // {45, 1.28f},
        // {48, 1.29f},
        // {51, 1.3f},
        // {54, 1.32f},
        // {57, 1.33f},
        // {60, 1.34f},
        // {170, 1.25f}//,
        // {171, 1.259f},
        // {172, 1.268f},
        // {173, 1.277f},
        // {174, 1.286f},
        {175, 1.295f},
        {176, 1.304f},
        {177, 1.313f},
        {178, 1.322f},
        {179, 1.331f},
        {180, 1.34f} //,
        // {181, 1.349f}//,
        // {182, 1.358f},
        // {183, 1.367f}
    };

    // ------------------ Mapping ------------------
    auto mapping = read_mapping_csv(mapping_csv, LED_SIPMS_PER_CRYSTAL);
    if (mapping.empty())
    {
        std::cerr << "Error: mapping is empty. Check mapping_csv='" << mapping_csv << "'\n";
        return;
    }

    std::vector<double> v_voltage;
    std::vector<double> v_reso_pct;
    std::vector<double> v_reso_pct_err;

    // ------------------ Calibration dependencies ------------------
    // adc_calibration.cxx writes mean_adc_to_ref_calibration_1.27V or 1.25V into outputs/adc_to_ref_calibration_1.27V.root or 1.25
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
                                // 120, 0, 2.5, 1200, 0, 60000);
                                120, 0, 2.5, 1200, 0, 4000);

    TH2F *total_invt = new TH2F("label_minus_missing_vs_tot",
                                "label-missing as f(ToT);ToT sum (a.u.);label - missing (V)",
                                // 1200, 0, 60000, 120, 0, 2.5);
                                1200, 0, 4000, 120, 0, 2.5);

    // Per-run containers (so we can also build widths if needed)
    struct RunOut
    {
        int run;
        float label;
        TH1F *missing;
        TH1F *tot;
        TH2F *tot2d;
        TH2F *cog;
        std::array<TH1F *, LED_SIPMS_PER_CRYSTAL> per_sipm;

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
                        //  1200, 0, 60000);
                         1200, 0, 4000);

        o.tot2d = new TH2F(Form("run%03d_tot_vs_label_minus_missing", runnum),
                           Form("Run %03d ToT vs (label-missing);label-missing (V);ToT sum (a.u.)", runnum),
                           120, 0, 2.5, 1200, 0, 4000);

        o.cog = new TH2F(Form("run%03d_cog", runnum),
                         Form("Run %03d COG;X (crystals);Y (crystals)", runnum),
                         100, -0.5, 4.5, 100, -0.5, 4.5);

        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
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

        // double mean = 0, mean_err = 0, sigma = 0, sigma_err = 0, reso = 0, reso_err = 0;
        // if (compute_resolution_pct(o.tot, mean, mean_err, sigma, sigma_err, reso, reso_err))
        // {
        //     v_voltage.push_back((double)label);
        //     v_reso_pct.push_back(reso);
        //     v_reso_pct_err.push_back(reso_err);
        // }
        // else
        // {
        //     std::cerr << "Warning: could not compute resolution for label=" << label
        //               << " (entries=" << (o.tot ? o.tot->GetEntries() : 0) << ")\n";
        // }

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

            const double reso_pct = (o.tot_mu != 0.0) ? 100.0 * (o.tot_sigma / o.tot_mu) : 0.0;

            // simple error propagation using fit parameter errors (if fit exists)
            double mu_err = 0.0, si_err = 0.0;
            if (TF1 *ffit = o.tot->GetFunction(Form("tot_sum_gaus_run%03d", runnum)))
            {
                mu_err = ffit->GetParError(1);
                si_err = ffit->GetParError(2);
            }
            double reso_err = 0.0;
            if (o.tot_mu > 0.0 && o.tot_sigma > 0.0)
            {
                const double rel2 = (mu_err > 0 ? (mu_err / o.tot_mu) * (mu_err / o.tot_mu) : 0.0) + (si_err > 0 ? (si_err / o.tot_sigma) * (si_err / o.tot_sigma) : 0.0);
                reso_err = reso_pct * std::sqrt(rel2);
            }

            v_voltage.push_back((double)label);
            v_reso_pct.push_back(reso_pct);
            v_reso_pct_err.push_back(reso_err);
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
    TCanvas *canvas = new TCanvas("tot_calib", "", 1680, 720);
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.04);
    text.SetTextFont(42);

    canvas->SaveAs((pdf + "[").c_str());

    // store peak mean/sigma per run to optionally write widths file
    std::vector<std::vector<double>> peak_mean(outs.size(), std::vector<double>(LED_SIPMS_PER_CRYSTAL, 0.0));
    std::vector<std::vector<double>> peak_sigma(outs.size(), std::vector<double>(LED_SIPMS_PER_CRYSTAL, 0.0));

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
        text.DrawLatex(0.88, 0.80, Form("Label: %.3f", o.label));

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
        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
        {
            pad->cd(sipm + 1);
            TH1F *h = o.per_sipm[sipm];
            if (!h)
                continue;

            h->SetTitle("");
            h->Draw("hist");

            TLatex lab;
            lab.SetNDC();
            lab.SetTextSize(0.08);

            // Skip SiPMs that are not part of the selected subset
            if (!use_selected_sipm(central_crystal, sipm))
            {
                lab.SetTextColor(kRed);
                lab.DrawLatex(0.5, 0.5, "Not used");
                continue;
            }

            // Skip empty or pathological histograms
            if (h->GetEntries() < 5 || h->GetRMS() <= 0.0 || h->GetMaximum() <= 0.0)
            {
                lab.DrawLatex(0.25, 0.82, "empty");
                peak_mean[ir][sipm] = 0.0;
                peak_sigma[ir][sipm] = 0.0;
                continue;
            }

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

                TFitResultPtr r2 = h->Fit(&fit2, "RQS");
                if (r2.Get() && r2->Status() == 0)
                {
                    peak_mean[ir][sipm] = fit2.GetParameter(1);
                    peak_sigma[ir][sipm] = fit2.GetParameter(2);
                    fit2.Draw("same");
                }
                else
                {
                    lab.DrawLatex(0.28, 0.82, "fit failed");
                    peak_mean[ir][sipm] = 0.0;
                    peak_sigma[ir][sipm] = 0.0;
                }
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
        text.DrawLatex(0.12, 0.78, Form("Label value: %.3f (treated as V in plots)", o.label));
        text.DrawLatex(0.12, 0.71, Form("ADC->ref: %.1f ADC per %.3fV", adc_per_refunit, g_adc_ref_voltage_V));
        text.DrawLatex(0.12, 0.64, Form("COG frac9 cut: %.3f; ellipse: %s", g_energy_fraction_cut, g_do_cog_ellipse_cut ? "ON" : "OFF"));
        text.DrawLatex(0.12, 0.57,
                       Form("ToT sum: mean=%.1f  sigma=%.1f  res=%.2f%%",
                            o.tot_mu, o.tot_sigma, 100.0 * o.tot_res));

        canvas->SaveAs(pdf.c_str());
    }

    if (v_voltage.size() >= 2)
    {
        TCanvas *cRes = new TCanvas("cRes", "Resolution vs Voltage", 900, 700);
        cRes->SetGrid();

        TGraphErrors *gRes = new TGraphErrors(
            (int)v_voltage.size(),
            v_voltage.data(),
            v_reso_pct.data(),
            nullptr,
            v_reso_pct_err.data());

        gRes->SetTitle("ToT resolution vs LED voltage;LED voltage (V);Resolution (#sigma/mean) [%]");
        gRes->SetMarkerStyle(20);
        gRes->Draw("AP");
        cRes->SaveAs(pdf.c_str());
    }
    else
    {
        std::cerr << "Warning: not enough points to build Resolution vs Voltage graph.\n";
    }

    // global page + fits
    canvas->Clear();
    // grid to this page
    gPad->SetGridx();
    gPad->SetGridy();
    total_dist->Draw("colz");
    canvas->SaveAs(pdf.c_str());

    // Fits of label as a function of ToT (optional, keep same idea as Tristan's pol1/pol2)
    // We fit on total_invt (x=ToT, y=label-missing)
    canvas->Clear();
    total_invt->Draw("colz");

    // float scale_factor = 16.0f / (float)LED_SIPMS_PER_CRYSTAL;

    // TF1 *fit_hi = new TF1("fit_hi", "pol1", 35000.0f / scale_factor, 50000.0f / scale_factor);
    // total_invt->Fit(fit_hi, "RQ");
    // fit_hi->SetLineColor(kRed);

    // TF1 *fit_mid = new TF1("fit_mid", "pol1", 20000.0f / scale_factor, 26000.0f / scale_factor);
    // total_invt->Fit(fit_mid, "RQ");
    // fit_mid->SetLineColor(kGreen + 2);

    // TF1 *fit_all = new TF1("fit_all", "pol2", 20000.0f / scale_factor, 55000.0f / scale_factor);
    // total_invt->Fit(fit_all, "RQ");
    // fit_all->SetLineColor(kMagenta);

    // fit_hi->Draw("same");
    // fit_mid->Draw("same");
    // fit_all->Draw("same");

    // grid this pad
    gPad->SetGridx();
    gPad->SetGridy();

    TProfile *p = total_invt->ProfileX("p_label_minus_missing_vs_tot");
    p->SetTitle("Profile: label-missing vs ToT;ToT;#LT label-missing#GT (V)");

    // TF1 *fit_all = new TF1("fit_all", "pol2", 20000, 55000);
    TF1 *fit_all = new TF1("fit_all", "pol2", 0, 6000);
    // TF1 *fit_mid = new TF1("fit_mid", "pol1", 14000, 23200);
    TF1 *fit_mid = new TF1("fit_mid", "pol1", 0, 5000);
    // TF1 *fit_hi = new TF1("fit_hi", "pol1", 30000, 60000);
    TF1 *fit_hi = new TF1("fit_hi", "pol1", 0, 5000);

    fit_hi->SetLineColor(kRed);
    fit_mid->SetLineColor(kGreen + 2);
    fit_all->SetLineColor(kOrange);

    p->Fit(fit_hi, "RQ");
    p->Fit(fit_mid, "RQ+");
    p->Fit(fit_all, "RQ+");

    p->Draw();
    fit_mid->Draw("same");
    fit_all->Draw("same");
    fit_hi->Draw("same");
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
        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
            if (o.per_sipm[sipm])
                o.per_sipm[sipm]->Write();
    }

    p->Write();

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
                                LED_SIPMS_PER_CRYSTAL, 0, LED_SIPMS_PER_CRYSTAL);
            for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
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