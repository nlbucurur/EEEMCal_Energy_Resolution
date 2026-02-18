// adc_calibration.cxx
//
// Refactor to match led_analysis.C / draw_waveform_reversed.cxx style:
// - Uses common_led.{h,cxx}
// - Uses mapping CSV from generate_mapping.py (FPGA,ASIC,Connector,Crystal,SiPM,Channel)
// - Builds reverse[channel] -> (crystal, sipm)
// - Loops over active channels and routes signals like led_analysis
// - Reads adc/tot/toa as flat buffers using TLeaf::GetLen()
// - Keeps adc_calibration purpose: energy hists, COG, pedestals, energy shares,
//   and produces an "ADC per reference unit" calibration.
//
// Reference unit: voltage = 1.27 V (instead of 1 GeV)
//
// Run examples:
//   root -l -q 'common_led.cxx+ adc_calibration.cxx+ adc_calibration_one("data/Run042.root", 1.27)'
//   root -l -q 'common_led.cxx+ adc_calibration.cxx+ adc_calibration_one("data/Run042.root", 1.27, "eeemcal_desy_dec2025_mapping_v2.csv", "outputs/gain_match_1.27V.root", "outputs")'
// root -l -b
// .L common_led.cxx+
// .L adc_calibration.cxx+
// adc_calibration_one("data/Run042.root", 1.27, "eeemcal_desy_dec2025_mapping_v2.csv", "outputs/gain_match_1.27V.root", "outputs")
//

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

#include <TCanvas.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>
#include <TProfile.h>
#include <TParameter.h>
#include <TError.h>
#include <TSystem.h>
#include <TLeaf.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

static constexpr int SAMPLES_PER_CHANNEL = 20;
static constexpr int SIPMS_PER_CRYSTAL = 16;
static constexpr int MAX_NUM_CRYSTALS = 25;

// ---------------------- small helpers ----------------------

static int extract_run_number(const char *filename)
{
    std::string s(filename);
    int run = -1;
    size_t pos = s.find("Run");
    if (pos != std::string::npos)
    {
        pos += 3;
        std::string digits;
        while (pos < s.size() && std::isdigit((unsigned char)s[pos]))
        {
            digits += s[pos];
            pos++;
        }
        if (!digits.empty())
            run = std::stoi(digits);
    }
    return run;
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

static int32_t get_toa_first_nonzero(uint32_t *toa_values)
{
    for (int i = 0; i < SAMPLES_PER_CHANNEL; ++i)
    {
        if (toa_values[i] > 0)
            return (int32_t)toa_values[i];
    }
    return -1;
}

static bool is_tot_event(uint32_t *tot_values)
{
    for (int i = 0; i < SAMPLES_PER_CHANNEL; ++i)
    {
        if (tot_values[i] > (uint32_t)g_tot_min)
            return true;
    }
    return false;
}

// 5x5 geometry mapping used for COG (same layout you used before).
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

// COG: simple energy-weighted position (only crystals present in event_energy).
// Returns true always, unless do_cut is enabled and ellipse cut fails.
static bool calculate_cog(TH2 *distribution,
                          const std::map<int, float> &event_energy,
                          bool do_cut,
                          float cx, float cy, float sx, float sy)
{
    double sumE = 0.0;
    double sumX = 0.0;
    double sumY = 0.0;

    for (const auto &kv : event_energy)
    {
        int cr = kv.first;
        float E = kv.second;
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

    if (sumE <= 0)
        return false;

    float x_cog = (float)(sumX / sumE);
    float y_cog = (float)(sumY / sumE);

    if (distribution)
        distribution->Fill(x_cog, y_cog);

    if (!do_cut)
        return true;

    float dx = (x_cog - cx) / sx;
    float dy = (y_cog - cy) / sy;
    return (dx * dx + dy * dy) <= 1.0f;
}

// Safe common-mode subtraction into a temporary waveform buffer
static void subtract_common_mode(uint32_t out_adc[SAMPLES_PER_CHANNEL],
                                 uint32_t *in_adc,
                                 float cm)
{
    for (int i = 0; i < SAMPLES_PER_CHANNEL; ++i)
    {
        float v = (float)in_adc[i] - cm;
        if (v < 0)
            v = 0;
        out_adc[i] = (uint32_t)(v + 0.5f);
    }
}

// ---------------------- main function ----------------------

void adc_calibration_one(const char *filename,
                         float voltage,
                         const char *mapping_csv = "eeemcal_desy_dec2025_mapping_v2.csv",
                         const char *gain_root_file = "outputs/gain_match_1.27V.root",
                         const char *outdir = "outputs",
                         Long64_t max_events = 1000000,
                         float energy_fraction_cut = 0.0f,
                         bool use_hybrid_tot = true,
                         bool use_common_mode = true)
{
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;
    gSystem->mkdir(outdir, true);

    // Use same signal config style as your other scripts
    g_signal_method = 3;
    g_tot_min = 50;

    const float REF_VOLTAGE = voltage;

    // Per-crystal common-mode reference channels
    // int common_mode_channels_by_crystal[MAX_NUM_CRYSTALS] = {
    //     224, 296, 188, 386, 152,
    //     512, 548, 242, 134, 404,
    //     206, 368, 332, 458, 8,
    //     440, 116, 350, 62, 98,
    //     422, 170, 260, 314, 278}; // used in DESY

    int common_mode_channels_by_crystal[MAX_NUM_CRYSTALS];
    for (int cr = 0; cr < MAX_NUM_CRYSTALS; ++cr)
        common_mode_channels_by_crystal[cr] = -1;

    // Only crystal 12 has CM on channel 44
    common_mode_channels_by_crystal[12] = 44; // Tedt LED one crystal with CM channel 44. Set to -1 for no CM.

    /*
        a - 8
        b - 26
        c - 62
        d - 44
    */

    // Geometry helpers
    const int central_crystal_id = 12;
    int center_nine_ids[8] = {7, 8, 9, 11, 13, 17, 18, 19};
    int remaining_ids[16] = {0, 1, 2, 3, 4, 5, 6, 10, 11, 15, 16, 20, 21, 22, 23, 24};

    // ---------------- mapping ----------------
    auto mapping = read_mapping_csv(mapping_csv, SIPMS_PER_CRYSTAL);
    if (mapping.empty())
    {
        std::cerr << "ERROR: mapping empty. Check: " << mapping_csv << "\n";
        return;
    }

    // reverse[channel] = (crystal, sipm)
    std::unordered_map<int, std::pair<int, int>> reverse;
    reverse.reserve(mapping.size() * SIPMS_PER_CRYSTAL);

    for (const auto &kv : mapping)
    {
        int cr = kv.first;
        auto chans = get_crystal_channels(mapping, cr, SIPMS_PER_CRYSTAL);
        for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; ++sipm)
        {
            int ch = chans[sipm];
            if (ch >= 0)
                reverse[ch] = {cr, sipm};
        }
    }

    // ---------------- load gain factors ----------------
    TH1 *gain_factor = nullptr;
    TH1 *crystal_gain_factor = nullptr;

    TFile *gf = TFile::Open(gain_root_file);
    if (gf && !gf->IsZombie())
    {
        gain_factor = (TH1 *)gf->Get("gain_factors");
        crystal_gain_factor = (TH1 *)gf->Get("crystal_factor");
        std::cout << "Loaded gain factors from: " << gain_root_file << "\n";
    }

    // Fallbacks
    if (!gain_factor)
    {
        gain_factor = new TH1F("gain_factors_fallback", "gain_factors_fallback", MAX_NUM_CRYSTALS * SIPMS_PER_CRYSTAL, 0, MAX_NUM_CRYSTALS * SIPMS_PER_CRYSTAL);
        for (int i = 1; i <= MAX_NUM_CRYSTALS * SIPMS_PER_CRYSTAL; ++i)
            gain_factor->SetBinContent(i, 1.0);
    }
    if (!crystal_gain_factor)
    {
        crystal_gain_factor = new TH1F("crystal_factor_fallback", "crystal_factor_fallback", MAX_NUM_CRYSTALS, 0, MAX_NUM_CRYSTALS);
        for (int i = 1; i <= MAX_NUM_CRYSTALS; ++i)
            crystal_gain_factor->SetBinContent(i, 1.0);
    }

    // ---------------- open data ----------------
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie())
    {
        std::cerr << "ERROR: cannot open " << filename << "\n";
        if (f)
        {
            f->Close();
            delete f;
        }
        if (gf)
        {
            gf->Close();
            delete gf;
        }
        return;
    }

    TTree *tree = (TTree *)f->Get("events");
    if (!tree)
    {
        std::cerr << "ERROR: cannot find TTree 'events' in " << filename << "\n";
        f->Close();
        delete f;
        if (gf)
        {
            gf->Close();
            delete gf;
        }
        return;
    }

    // Branch leaves for dynamic sizing
    TBranch *br_adc = tree->GetBranch("adc");
    TBranch *br_tot = tree->GetBranch("tot");
    TBranch *br_toa = tree->GetBranch("toa");

    if (!br_adc)
    {
        std::cerr << "ERROR: missing branch adc\n";
        f->Close();
        delete f;
        if (gf)
        {
            gf->Close();
            delete gf;
        }
        return;
    }

    TLeaf *leaf_adc = br_adc->GetLeaf("adc");
    if (!leaf_adc)
    {
        std::cerr << "ERROR: missing leaf adc\n";
        f->Close();
        delete f;
        if (gf)
        {
            gf->Close();
            delete gf;
        }
        return;
    }

    TLeaf *leaf_tot = (br_tot ? br_tot->GetLeaf("tot") : nullptr);
    TLeaf *leaf_toa = (br_toa ? br_toa->GetLeaf("toa") : nullptr);

    const int n_adc = leaf_adc->GetLen();
    if (n_adc <= 0 || (n_adc % SAMPLES_PER_CHANNEL) != 0)
    {
        std::cerr << "ERROR: unexpected adc length " << n_adc
                  << " (not divisible by " << SAMPLES_PER_CHANNEL << ")\n";
        f->Close();
        delete f;
        if (gf)
        {
            gf->Close();
            delete gf;
        }
        return;
    }

    const int n_channels = n_adc / SAMPLES_PER_CHANNEL;

    bool have_tot = (leaf_tot && leaf_tot->GetLen() == n_adc);
    bool have_toa = (leaf_toa && leaf_toa->GetLen() == n_adc);

    if (use_hybrid_tot && !have_tot)
    {
        std::cerr << "Warning: requested hybrid ToT but tot branch/len not compatible. Using ADC-only.\n";
        use_hybrid_tot = false;
    }

    std::vector<uint32_t> adc_buf((size_t)n_adc);
    std::vector<uint32_t> tot_buf;
    std::vector<uint32_t> toa_buf;

    if (use_hybrid_tot)
        tot_buf.resize((size_t)n_adc);
    if (have_toa)
        toa_buf.resize((size_t)n_adc);

    tree->SetBranchAddress("adc", adc_buf.data());
    if (use_hybrid_tot)
        tree->SetBranchAddress("tot", tot_buf.data());
    if (have_toa)
        tree->SetBranchAddress("toa", toa_buf.data());

    auto adc_ptr = [&](int ch) -> uint32_t *
    { return &adc_buf[(size_t)ch * SAMPLES_PER_CHANNEL]; };
    auto tot_ptr = [&](int ch) -> uint32_t *
    { return &tot_buf[(size_t)ch * SAMPLES_PER_CHANNEL]; };
    auto toa_ptr = [&](int ch) -> uint32_t *
    { return &toa_buf[(size_t)ch * SAMPLES_PER_CHANNEL]; };

    // Active channels based on mapping, bounded by file channels
    auto active_channels = get_active_channels_from_mapping(mapping, SIPMS_PER_CRYSTAL, MAX_NUM_CRYSTALS, n_channels);

    // ---------------- Histograms (global) ----------------
    TH1 *central_crystal_energy = new TH1F("central_crystal_energy",
                                           Form("Central Crystal Energy (%.2f V);Energy (ADC);Events", voltage),
                                           500, 0, 75000);

    TH1 *central_nine_energy = new TH1F("central_nine_energy",
                                        Form("Central 3x3 Energy (%.2f V);Energy (ADC);Events", voltage),
                                        500, 0, 75000);

    TH1 *total_energy = new TH1F("total_energy",
                                 Form("Total Energy (%.2f V);Energy (ADC);Events", voltage),
                                 500, 0, 75000);

    TH2 *cog_distribution = new TH2F("cog_distribution",
                                     "Center of Gravity Distribution;X (# Crystals);Y (# Crystals)",
                                     100, -0.5, 4.5, 100, -0.5, 4.5);

    TH2 *pedestals = new TH2F("pedestals", "Pedestals;Channel;Pedestal sample (ADC)",
                              n_channels, 0, n_channels, 1024, 0, 1024);

    TH1 *toa_distribution = new TH1F("toa_distribution", "ToA Distribution;ToA;Events", 1024, 0, 1024);
    TH1 *toa_sample = new TH1F("toa_sample", "ToA Sample Distribution;Sample;Events", 20, 0, 20);

    TH1 *max_sample_index = new TH1F("max_sample_index", "Max Sample Index;Sample;Events", 20, 0, 20);
    TH1 *second_max_sample_index = new TH1F("second_max_sample_index", "Second Max Sample Index;Sample;Events", 20, 0, 20);
    TH1 *third_max_sample_index = new TH1F("third_max_sample_index", "Third Max Sample Index;Sample;Events", 20, 0, 20);

    // E vs ToA for central crystal SiPMs
    std::vector<TH2 *> E_vs_toa;
    E_vs_toa.reserve(SIPMS_PER_CRYSTAL);
    for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; sipm++)
    {
        E_vs_toa.push_back(new TH2F(Form("E_vs_toa_sipm_%d", sipm),
                                    Form("Central Crystal SiPM %d: Energy vs ToA;ToA;Energy (ADC)", sipm),
                                    1024, 0, 1024, 500, 0, 3000));
    }

    // ToA correlations (central crystal only)
    std::vector<std::vector<TH2 *>> toa_correlations(SIPMS_PER_CRYSTAL);
    for (int i = 0; i < SIPMS_PER_CRYSTAL; i++)
    {
        toa_correlations[i].reserve(SIPMS_PER_CRYSTAL);
        for (int j = 0; j < SIPMS_PER_CRYSTAL; j++)
        {
            toa_correlations[i].push_back(new TH2F(Form("toa_correlation_sipm_%02d_sipm_%02d", i, j),
                                                   Form("Central crystal: SiPM %02d vs %02d ToA;SiPM %02d ToA;SiPM %02d ToA", i, j, i, j),
                                                   1024, 0, 1024, 1024, 0, 1024));
        }
    }

    // ---------------- per-crystal/per-sipm histograms ----------------
    // Allocate for 25 crystals, but only fill when that crystal exists in mapping.
    std::vector<std::array<TH1F *, SIPMS_PER_CRYSTAL>> sipm_energy(MAX_NUM_CRYSTALS);
    std::vector<std::array<TH2F *, SIPMS_PER_CRYSTAL>> sipm_waveform(MAX_NUM_CRYSTALS);
    std::vector<TH1F *> crystal_energy(MAX_NUM_CRYSTALS, nullptr);
    std::vector<TH1F *> crystal_energy_shares(MAX_NUM_CRYSTALS, nullptr);
    std::vector<TH2F *> crystal_common_mode(MAX_NUM_CRYSTALS, nullptr);

    for (int cr = 0; cr < MAX_NUM_CRYSTALS; ++cr)
    {
        for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; ++sipm)
        {
            sipm_energy[cr][sipm] = new TH1F(Form("crystal_%02d_sipm_%02d_energy", cr, sipm),
                                             Form("Crystal %02d SiPM %02d Energy;Energy (ADC);Events", cr, sipm),
                                             500, 0, 4000);

            sipm_waveform[cr][sipm] = new TH2F(Form("crystal_%02d_sipm_%02d_waveform", cr, sipm),
                                               Form("Crystal %02d SiPM %02d Waveform;Sample;ADC", cr, sipm),
                                               20, 0, 20, 1024, 0, 1024);
        }

        crystal_energy[cr] = new TH1F(Form("crystal_%02d_energy", cr),
                                      Form("Crystal %02d Energy;Energy (ADC);Events", cr),
                                      500, 0, 40000);

        crystal_energy_shares[cr] = new TH1F(Form("crystal_%02d_energy_share", cr),
                                             Form("Crystal %02d Energy Share;Energy Share;Events", cr),
                                             10000, 0, 1.05);

        crystal_common_mode[cr] = new TH2F(Form("crystal_%02d_common_mode", cr),
                                           Form("Crystal %02d Common Mode;Sample;ADC", cr),
                                           20, 0, 20, 1024, 0, 1024);
    }

    // COG cut parameters (kept but disabled by default for robustness)
    const bool DO_COG_CUT = false; // set true if you want the ellipse selection
    const float COG_CX = 2.0f;
    const float COG_CY = 2.0f;
    const float COG_SX = 0.60f;
    const float COG_SY = 0.60f;

    // ---------------- event loop ----------------
    Long64_t nentries = tree->GetEntries();
    if (max_events > 0 && max_events < nentries)
        nentries = max_events;

    int run = extract_run_number(filename);

    std::cout << "Processing " << nentries << " events from " << filename
              << " (run=" << run << ", voltage=" << voltage << ", channels=" << n_channels << ")\n";

    int tot_events = 0;
    int complete = 0;

    for (Long64_t entry = 0; entry < nentries; ++entry)
    {
        if ((int)(entry * MAX_NUM_CRYSTALS / nentries) > complete)
        {
            complete = (int)(entry * MAX_NUM_CRYSTALS / nentries);
            print_progress_25(complete);
        }

        tree->GetEntry(entry);

        std::unordered_map<int, float> cm_cache;
        std::unordered_map<int, bool> cm_ok;  

        bool is_tot = false;

        // accumulate per-crystal signal for this event
        std::map<int, float> signals;
        for (const auto &kv : mapping)
            signals[kv.first] = 0.0f;

        // central crystal ToA window check (kept; disabled by default)
        bool outside_toa_range = false;
        if (have_toa && mapping.count(central_crystal_id))
        {
            // use central crystal sipm 0 as reference ToA channel
            int ch0 = get_crystal_channels(mapping, central_crystal_id, SIPMS_PER_CRYSTAL)[0];
            if (ch0 >= 0 && ch0 < n_channels)
            {
                uint32_t *toa0 = toa_ptr(ch0);
                int32_t central_toa = get_toa_first_nonzero(toa0);
                if (central_toa > 0)
                {
                    if (central_toa < 200 || central_toa > 800)
                        outside_toa_range = true;
                }
            }
        }
        // if (outside_toa_range) continue; // optional

        // loop active channels
        for (int ch : active_channels)
        {
            if (ch < 0 || ch >= n_channels)
                continue;

            auto it = reverse.find(ch);
            if (it == reverse.end())
                continue;

            int cr = it->second.first;
            int sipm = it->second.second;

            if (cr < 0 || cr >= MAX_NUM_CRYSTALS || sipm < 0 || sipm >= SIPMS_PER_CRYSTAL)
                continue;
            if (!mapping.count(cr))
                continue; // ignore crystals not present

            // gain
            float gain = gain_factor->GetBinContent(cr * SIPMS_PER_CRYSTAL + sipm + 1);

            // common-mode (optional, per crystal)
            float cm = 0.0f;
            bool do_cm = false;
            if (use_common_mode)
            {
                auto itcm = cm_ok.find(cr);
                if (itcm == cm_ok.end())
                {
                    int cm_ch = common_mode_channels_by_crystal[cr];
                    if (cm_ch >= 0 && cm_ch < n_channels)
                    {
                        uint32_t *cm_adc = adc_ptr(cm_ch);
                        cm = (cm_adc[0] + cm_adc[1] /*+ cm_adc[2]*/) / 2.0f;
                        cm_cache[cr] = cm;
                        cm_ok[cr] = true;

                        // record CM waveform
                        for (int s = 0; s < SAMPLES_PER_CHANNEL; ++s)
                        {
                            crystal_common_mode[cr]->Fill(s, cm_adc[s]);
                        }
                    }
                    else
                    {
                        cm_ok[cr] = false;
                    }
                }

                if (cm_ok[cr])
                {
                    cm = cm_cache[cr];
                    do_cm = true;
                }
            }

            // prepare waveform pointer (with optional CM subtraction)
            uint32_t adc_corr[SAMPLES_PER_CHANNEL];
            uint32_t *adc_used = adc_ptr(ch);
            if (do_cm)
            {
                subtract_common_mode(adc_corr, adc_used, cm);
                adc_used = adc_corr;
            }

            float channel_signal = 0.0f;
            bool used_tot = false;

            if (use_hybrid_tot && have_tot)
            {
                // Hybrid decides per-event whether ToT is used
                calculate_signal_hybrid(adc_used, tot_ptr(ch), gain, channel_signal, used_tot);

                // IMPORTANT: only tag ToT when hybrid actually used it
                if (used_tot)
                    is_tot = true;
            }
            else
            {
                // ADC-only: do NOT look at ToT at all (and do NOT call tot_ptr)
                channel_signal = calculate_signal_adc(adc_used, gain);
                used_tot = false;
            }

            // per sipm energy
            sipm_energy[cr][sipm]->Fill(channel_signal);

            // waveform fill
            float pedestal = (adc_used[0] + adc_used[1] /*+ adc_used[2]*/) / 2.0f;
            int i1 = -1, i2 = -1, i3 = -1;
            float max1 = -1e9f, max2 = -1e9f, max3 = -1e9f;

            for (int s = 0; s < SAMPLES_PER_CHANNEL; ++s)
            {
                sipm_waveform[cr][sipm]->Fill(s, adc_used[s]);

                float sample = (float)adc_used[s] - pedestal;
                if (sample > max1)
                {
                    i3 = i2;
                    max3 = max2;
                    i2 = i1;
                    max2 = max1;
                    i1 = s;
                    max1 = sample;
                }
                else if (sample > max2)
                {
                    i3 = i2;
                    max3 = max2;
                    i2 = s;
                    max2 = sample;
                }
                else if (sample > max3)
                {
                    i3 = s;
                    max3 = sample;
                }
            }

            if (i1 >= 0)
                max_sample_index->Fill(i1);
            if (i2 >= 0)
                second_max_sample_index->Fill(i2);
            if (i3 >= 0)
                third_max_sample_index->Fill(i3);

            // pedestal samples
            pedestals->Fill(ch, adc_used[0]);
            pedestals->Fill(ch, adc_used[1]);
            // pedestals->Fill(ch, adc_used[2]);

            // sum into crystal (apply per-crystal gain factor later)
            signals[cr] += channel_signal;

            // TOA logic (central crystal only, like your old code)
            if (have_toa && cr == central_crystal_id)
            {
                int32_t toa_val = get_toa_first_nonzero(toa_ptr(ch));
                if (toa_val >= 0)
                {
                    // record distribution
                    toa_distribution->Fill(toa_val);

                    // which sample had the first non-zero TOA?
                    int timebin = -1;
                    uint32_t *tv = toa_ptr(ch);
                    for (int s = 0; s < SAMPLES_PER_CHANNEL; ++s)
                    {
                        if (tv[s] > 0)
                        {
                            timebin = s;
                            break;
                        }
                    }
                    toa_sample->Fill(timebin);

                    // E vs ToA uses sipm index (not channel index)
                    E_vs_toa[sipm]->Fill(toa_val, channel_signal);

                    // correlations with other central SiPMs
                    auto chans = get_crystal_channels(mapping, central_crystal_id, SIPMS_PER_CRYSTAL);
                    for (int other = 0; other < SIPMS_PER_CRYSTAL; ++other)
                    {
                        if (other == sipm)
                            continue;
                        int ch_other = chans[other];
                        if (ch_other < 0 || ch_other >= n_channels)
                            continue;
                        int32_t toa_other = get_toa_first_nonzero(toa_ptr(ch_other));
                        if (toa_other < 0)
                            continue;
                        toa_correlations[sipm][other]->Fill(toa_val, toa_other);
                    }
                }
            }
        } // end channel loop

        if (is_tot)
        {
            tot_events++;
            // if you want to reject TOT events: continue;
        }

        // apply per-crystal gain factor
        for (auto &kv : signals)
        {
            int cr = kv.first;
            kv.second *= crystal_gain_factor->GetBinContent(cr + 1);
        }

        // COG selection
        bool keep = calculate_cog(cog_distribution, signals, DO_COG_CUT, COG_CX, COG_CY, COG_SX, COG_SY);
        keep &= (!is_tot); // same as your old keep &= (!is_tot_event)
        if (!keep)
            continue;

        // build energies
        float central_signal = 0.0f;
        float central_nine_signal = 0.0f;
        float total_signal = 0.0f;

        if (signals.count(central_crystal_id))
            central_signal = signals[central_crystal_id];

        for (int i = 0; i < 8; i++)
        {
            int cr = center_nine_ids[i];
            if (signals.count(cr))
                central_nine_signal += signals[cr];
        }

        for (int i = 0; i < SIPMS_PER_CRYSTAL; i++)
        {
            int cr = remaining_ids[i];
            if (signals.count(cr))
                total_signal += signals[cr];
        }

        // In your old code: if signals[12]/total_signal < cut skip.
        if (total_signal > 0 && signals.count(central_crystal_id))
        {
            if (signals[central_crystal_id] / total_signal < energy_fraction_cut)
                continue;
        }

        central_nine_signal += central_signal;
        total_signal += central_nine_signal;

        central_crystal_energy->Fill(central_signal);
        central_nine_energy->Fill(central_nine_signal);
        total_energy->Fill(total_signal);

        // per-crystal energy + share
        if (total_signal > 0)
        {
            for (int cr = 0; cr < MAX_NUM_CRYSTALS; cr++)
            {
                if (!signals.count(cr))
                    continue;
                crystal_energy[cr]->Fill(signals[cr]);
                float share = signals[cr] / total_signal;
                if (share >= 1.0f)
                    share = 0.999999f;
                crystal_energy_shares[cr]->Fill(share);
            }
        }
    }

    print_progress_25(MAX_NUM_CRYSTALS);

    std::cout << "Total ToT-tagged events: " << tot_events << " out of " << nentries << "\n";

    // ---------------- calibration: ADC per reference unit (1.27 V) ----------------
    // Same formula you used:
    //   signal_for_1gev = mean(crystal_energy) / mean(crystal_energy_share)
    // Here renamed:
    //   signal_for_ref = mean(crystal_energy) / mean(crystal_energy_share)
    // where "ref" is defined by voltage 1.27 V.
    int crystal_mapping_for_plots[MAX_NUM_CRYSTALS] = {
        4, 9, 14, 19, 24,
        3, 8, 13, 18, 23,
        2, 7, 12, 17, 22,
        1, 6, 11, 16, 21,
        0, 5, 10, 15, 20};

    std::vector<float> ref_calib;
    float mean_calib = 0.0f;
    int n_used = 0;

    for (int i = 0; i < MAX_NUM_CRYSTALS; i++)
    {
        int cr = crystal_mapping_for_plots[i];
        if (!mapping.count(cr))
            continue; // only crystals present
        if (cr == 9)
            continue; // keep exclusion

        float meanE = (float)crystal_energy[cr]->GetMean();
        float meanS = (float)crystal_energy_shares[cr]->GetMean();
        if (meanE <= 0 || meanS <= 0)
            continue;

        float adc_per_ref = meanE / meanS;
        ref_calib.push_back(adc_per_ref);
        mean_calib += adc_per_ref;
        n_used++;
    }

    if (n_used > 0)
        mean_calib /= n_used;
    std::cout << "Reference calibration (" << REF_VOLTAGE << " V defined as 1 unit): " << mean_calib << "\n";

    // ---------------- save outputs ----------------
    int run_out = extract_run_number(filename);

    std::string pdf = std::string(Form("%s/adc_calibration_Run%03d_%.2fV.pdf", outdir, run_out, voltage));
    std::string root_out = std::string(Form("%s/adc_to_ref_calibration_%.2fV.root", outdir, REF_VOLTAGE));

    TCanvas *canvas = new TCanvas("adc_calibration_canvas", "", 900, 700);
    canvas->SetRightMargin(0.05);

    canvas->SaveAs((pdf + "[").c_str());

    central_crystal_energy->Draw("hist e");
    canvas->SaveAs(pdf.c_str());

    central_nine_energy->Draw("hist e");
    canvas->SaveAs(pdf.c_str());

    total_energy->Draw("hist e");
    canvas->SaveAs(pdf.c_str());

    // COG plot + grid + ellipse visualization if you want
    canvas->SetRightMargin(0.10);
    cog_distribution->Draw("colz");
    gPad->SetLogz();

    for (int i = 1; i < 5; i++)
    {
        float loc = i - 0.5f;
        TLine *line_x = new TLine(loc, -0.5, loc, 4.5);
        line_x->SetLineStyle(2);
        line_x->Draw();
        TLine *line_y = new TLine(-0.5, loc, 4.5, loc);
        line_y->SetLineStyle(2);
        line_y->Draw();
    }

    TLine *line_x = new TLine(2, -0.5, 2, 4.5);
    line_x->SetLineStyle(2);
    line_x->SetLineColor(kRed);
    line_x->Draw();
    TLine *line_y = new TLine(-0.5, 2, 4.5, 2);
    line_y->SetLineStyle(2);
    line_y->SetLineColor(kRed);
    line_y->Draw();

    TEllipse *ellipse = new TEllipse(COG_CX, COG_CY, COG_SX, COG_SY);
    ellipse->SetLineColor(kRed);
    ellipse->SetLineWidth(2);
    ellipse->SetFillStyle(0);
    ellipse->Draw();

    canvas->SaveAs(pdf.c_str());

    // pedestals
    canvas->Clear();
    gPad->SetLogz(0);
    pedestals->Draw("colz");
    // pedestals->GetYaxis()->SetRangeUser(0, 200);
    canvas->SaveAs(pdf.c_str());

    // pedestal widths
    canvas->Clear();
    TH1F *pedestals_width = new TH1F("pedestals_width", "Pedestal Width vs Channel;Channel;Width (ADC)",
                                     n_channels, 0, n_channels);

    for (int ch = 0; ch < n_channels; ++ch)
    {
        int xbin = pedestals->GetXaxis()->FindBin(ch + 0.5); // channel center
        TH1D *proj = pedestals->ProjectionY(Form("_py_ch%d", ch), xbin, xbin);
        if (proj && proj->GetEntries() > 1)
            pedestals_width->SetBinContent(xbin, proj->GetStdDev());
        delete proj;
    }

    std::cout << "Pedestals entries: " << pedestals->GetEntries() << "\n";
    std::cout << "Pedestals non-empty x-bins: ";
    for (int ch = 0; ch < n_channels; ++ch)
    {
        int xbin = pedestals->GetXaxis()->FindBin(ch + 0.5);
        auto *p = pedestals->ProjectionY(Form("tmp_ch%03d", ch), xbin, xbin);
        p->SetDirectory(nullptr);
        if (p->GetEntries() > 1)
            std::cout << ch << " ";
        delete p;
    }
    std::cout << "\n";

    // pedestals_width->Draw("P");
    pedestals_width->Draw();
    // pedestals_width->GetYaxis()->SetRangeUser(0, 5);
    canvas->SaveAs(pdf.c_str());

    // crystal energy (5x5 pages)
    canvas->Clear();
    canvas->Divide(5, 5);
    for (int i = 0; i < MAX_NUM_CRYSTALS; i++)
    {
        canvas->cd(i + 1);
        int cr = crystal_mapping_for_plots[i];
        crystal_energy[cr]->Draw("hist e");
    }
    canvas->SaveAs(pdf.c_str());

    // crystal shares + calibration text
    canvas->Clear();
    canvas->Divide(5, 5);
    float mean_calib_check = 0.0f;
    int used_check = 0;

    for (int i = 0; i < MAX_NUM_CRYSTALS; i++)
    {
        canvas->cd(i + 1);
        int cr = crystal_mapping_for_plots[i];

        crystal_energy_shares[cr]->Draw("hist e");
        crystal_energy_shares[cr]->GetXaxis()->SetRangeUser(0.001, 1);
        gPad->SetLogx();

        if (mapping.count(cr) && cr != 9)
        {
            float adc_per_ref = (float)crystal_energy[cr]->GetMean() / (float)crystal_energy_shares[cr]->GetMean();
            mean_calib_check += adc_per_ref;
            used_check++;

            TLatex t;
            t.SetNDC();
            t.SetTextSize(0.04);
            float x_coord = (i == 12 ? 0.15f : 0.40f);
            t.DrawLatex(x_coord, 0.82, Form("ADC @ %.2fV: %.1f", REF_VOLTAGE, adc_per_ref));
        }
    }
    if (used_check > 0)
        mean_calib_check /= used_check;

    canvas->SaveAs(pdf.c_str());
    gPad->SetLogx(0);

    // per-crystal pages: sipm energy + waveform + CM waveform
    for (int cr = 0; cr < MAX_NUM_CRYSTALS; cr++)
    {
        if (!mapping.count(cr))
            continue;

        canvas->Clear();
        canvas->Divide(4, 4);
        for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; sipm++)
        {
            canvas->cd(sipm + 1);
            sipm_energy[cr][sipm]->Draw("hist e");
        }
        canvas->SaveAs(pdf.c_str());

        canvas->Clear();
        canvas->Divide(4, 4);
        for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; sipm++)
        {
            canvas->cd(sipm + 1);
            sipm_waveform[cr][sipm]->Draw("colz");
            gPad->SetLogz();
        }
        canvas->SaveAs(pdf.c_str());
        gPad->SetLogz(0);

        canvas->Clear();
        crystal_common_mode[cr]->Draw("colz");
        gPad->SetLogz();
        canvas->SaveAs(pdf.c_str());
        gPad->SetLogz(0);
    }

    // ToA plots
    canvas->Clear();
    toa_distribution->Draw("hist");
    canvas->SaveAs(pdf.c_str());

    canvas->Clear();
    toa_sample->Draw("hist");
    canvas->SaveAs(pdf.c_str());

    canvas->Clear();
    canvas->Divide(4, 4);
    for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; sipm++)
    {
        canvas->cd(sipm + 1);
        E_vs_toa[sipm]->Draw("colz");
        auto *prof = E_vs_toa[sipm]->ProfileX(Form("prof_%02d", sipm));
        prof->SetDirectory(nullptr);
        prof->Draw("same");
        delete prof;
        gPad->SetLogz();
    }
    canvas->SaveAs(pdf.c_str());
    gPad->SetLogz(0);

    // ToA correlation matrix (central crystal)
    TCanvas *canvas2 = new TCanvas("toa_correlations_canvas", "", 8000, 8000);
    canvas2->Divide(16, 16, 0.0005, 0.0005);
    for (int i = 0; i < SIPMS_PER_CRYSTAL; i++)
    {
        for (int j = 0; j < SIPMS_PER_CRYSTAL; j++)
        {
            if (j <= i)
            {
                canvas2->cd(i * SIPMS_PER_CRYSTAL + j + 1);
                toa_correlations[i][j]->Draw("colz");
            }
        }
    }
    std::string corr_png = std::string(Form("%s/adc_toa_correlation_Run%03d.png", outdir, run_out));
    canvas2->SaveAs(corr_png.c_str());
    delete canvas2;

    // Max sample index plots
    canvas->Clear();
    max_sample_index->Draw("hist");
    canvas->SaveAs(pdf.c_str());

    canvas->Clear();
    second_max_sample_index->Draw("hist");
    canvas->SaveAs(pdf.c_str());

    canvas->Clear();
    third_max_sample_index->Draw("hist");
    canvas->SaveAs(pdf.c_str());

    // Final page with summary text
    canvas->Clear();
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.04);
    t.DrawLatex(0.12, 0.82, Form("Run: %03d   Voltage label: %.2f V", run_out, voltage));
    t.DrawLatex(0.12, 0.76, Form("Reference unit: %.2f V (defined as 1 unit)", REF_VOLTAGE));
    t.DrawLatex(0.12, 0.70, Form("Mean ADC per ref unit: %.2f", mean_calib));
    t.DrawLatex(0.12, 0.64, Form("ToT mode used: %s", use_hybrid_tot ? "hybrid ADC/ToT" : "ADC-only"));
    t.DrawLatex(0.12, 0.58, Form("Common-mode subtraction: %s", use_common_mode ? "enabled" : "disabled"));
    canvas->SaveAs(pdf.c_str());

    canvas->SaveAs((pdf + "]").c_str());
    std::cout << "Saved PDF: " << pdf << "\n";
    std::cout << "Saved ToA correlation PNG: " << corr_png << "\n";

    // ROOT output for the reference calibration
    TFile *out = TFile::Open(root_out.c_str(), "RECREATE");
    TParameter<float> p(Form("mean_adc_to_ref_calibration_%.2fV", REF_VOLTAGE), mean_calib);
    p.Write();

    central_crystal_energy->Write();
    central_nine_energy->Write();
    total_energy->Write();
    cog_distribution->Write();
    pedestals->Write();
    pedestals_width->Write();
    toa_distribution->Write();
    toa_sample->Write();
    max_sample_index->Write();
    second_max_sample_index->Write();
    third_max_sample_index->Write();

    for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; ++sipm)
    {
        E_vs_toa[sipm]->Write();
    }

    for (int i = 0; i < SIPMS_PER_CRYSTAL; i++)
    {
        for (int j = 0; j < SIPMS_PER_CRYSTAL; j++)
        {
            toa_correlations[i][j]->Write();
        }
    }

    for (int cr = 0; cr < MAX_NUM_CRYSTALS; ++cr)
    {
        if (!mapping.count(cr))
            continue;
        crystal_energy[cr]->Write();
        crystal_energy_shares[cr]->Write();
        crystal_common_mode[cr]->Write();
        for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; ++sipm)
        {
            sipm_energy[cr][sipm]->Write();
            sipm_waveform[cr][sipm]->Write();
        }
    }

    out->Close();
    delete out;
    std::cout << "Saved ROOT: " << root_out << "\n";

    // cleanup files
    if (gf)
    {
        gf->Close();
        delete gf;
    }
    f->Close();
    delete f;

    delete canvas;
}