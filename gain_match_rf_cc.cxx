// gain_match_rf_cc.cxx
//
// Purpose: find per-channel gain factors by fitting the peak per SiPM,
// then compute gain_factors[channel] = (mean peak of crystal) / (peak of that SiPM).
// Optionally compute per-crystal equalization factors.
//
// Run interactive:
//   root -l -b
//   .L common_led.cxx+
//   .L gain_match_rf_cc.cxx+
//   gain_scan_rf_cc("data", "eeemcal_desy_dec2025_mapping_v2.csv", "outputs_rf_cc", false);

// Notes:
// - This code reads adc/tot with dynamic sizes (TLeaf::GetLen()) like led_analysis.C.
// - It uses mapping + reverse[channel] to route each channel into the correct (crystal,sipm) histogram.
// - It works for mapping files that contain only 1 crystal OR many crystals.

#include "common_led.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TError.h>
#include <TParameter.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

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

// ------------------------------------------------------------
// Robust peak finder
//   1) gaussian in mean±1.5*rms
//   2) gaussian in peak±sigma
//   3) crystalball in peak±sigma
// Returns true if final fit succeeded.
// ------------------------------------------------------------
static bool fit_peak_crystalball(TH1 *hist, float &peak_out, float &sigma_out)
{
    if (!hist || hist->GetEntries() < 30)
        return false;

    const double xMin = hist->GetXaxis()->GetXmin();
    const double xMax = hist->GetXaxis()->GetXmax();

    const double mean = hist->GetMean();
    const double rms = hist->GetRMS();
    if (!(rms > 0))
        return false;

    // Fit window
    double fmin = mean - 1.5 * rms;
    double fmax = mean + 1.5 * rms;

    // Clamp to histogram limits
    if (fmin < xMin)
        fmin = xMin;
    if (fmax > xMax)
        fmax = xMax;

    // If the window is too small, bail out
    if (fmax - fmin < 5.0)
        return false;

    // Require some entries inside the window
    int bmin = hist->GetXaxis()->FindBin(fmin);
    int bmax = hist->GetXaxis()->FindBin(fmax);
    if (hist->Integral(bmin, bmax) < 20)
        return false;

    // 1) rough gaussian
    TF1 rough("rough_fit", "gaus", fmin, fmax);
    if (hist->Fit(&rough, "RQS") != 0)
        return false;
    rough.SetLineColor(kBlue);
    rough.SetLineStyle(2);

    double peak = rough.GetParameter(1);
    double sigma = rough.GetParameter(2);
    if (!(sigma > 0))
        return false;

    // 2) refined gaussian (clamped)
    double gmin = std::max(xMin, peak - sigma);
    double gmax = std::min(xMax, peak + sigma);
    if (gmax - gmin < 5.0)
    {
        peak_out = (float)peak;
        sigma_out = (float)sigma;
        return true;
    }

    TF1 second("second_fit", "gaus", gmin, gmax);
    if (hist->Fit(&second, "RQ") == 0)
    {
        peak = second.GetParameter(1);
        sigma = second.GetParameter(2);
    }

    if (!(sigma > 0))
    {
        peak_out = (float)peak;
        sigma_out = (float)sigma;
        return true;
    }

    // 3) crystalball (clamped)
    double cmin = std::max(xMin, peak - sigma);
    double cmax = std::min(xMax, peak + sigma);

    if (cmax - cmin < 5.0)
    {
        peak_out = (float)peak;
        sigma_out = (float)sigma;
        return true;
    }

    TF1 cb("final_fit", "crystalball", cmin, cmax);
    cb.SetParameters(second.GetParameter(0), peak, sigma, 1.5, 2.0);
    if (hist->Fit(&cb, "RMQ") != 0)
    {
        peak_out = (float)peak;
        sigma_out = (float)sigma;
        return true;
    }

    peak_out = (float)cb.GetParameter(1);
    sigma_out = (float)cb.GetParameter(2);

    if (hist->GetFunction("final_fit"))
        std::cout << "CB final_fit is stored\n";
    if (hist->GetFunction("second_fit"))
        std::cout << "Gaussian second_fit is stored\n";
    if (hist->GetFunction("rough_fit"))
        std::cout << "Gaussian rough_fit is stored\n";
    return true;
}

// ------------------------------------------------------------
// Dead-channel masks
// Return true if that (crystal,sipm) should be excluded from crystal mean.
// ------------------------------------------------------------
static bool exclude_from_crystal_mean(int crystal, int sipm)
{
    // Crystal 21 dead channels
    if (crystal == 21 && (sipm == 1 || sipm == 4 || sipm == 10 || sipm == 11 || sipm == 13 || sipm == 14 || sipm == 15))
        return true;

    // Crystal 15 weird channels
    if (crystal == 15 && (sipm == 10 || sipm == 0 || sipm == 4 || sipm == 13 || sipm == 11))
        return true;

    return false;
}

// ------------------------------------------------------------
// Core: gain match for ONE file (one run, one voltage label)
// ------------------------------------------------------------
void gain_match_rf_cc_one(const char *filename,
                          float voltage,
                          const char *mapping_csv = "eeemcal_desy_dec2025_mapping_v2.csv",
                          const char *outdir = "outputs",
                          bool use_hybrid_tot = true,
                          int cc = -1,
                          int rf = -1)
{
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    g_signal_method = 3; // same as led_analysis default
    g_tot_min = 50;

    gSystem->mkdir(outdir, true);

    // ---- Read mapping
    auto mapping = read_mapping_csv(mapping_csv, LED_SIPMS_PER_CRYSTAL);
    if (mapping.empty())
    {
        std::cerr << "Error: mapping empty. Check mapping CSV: " << mapping_csv << "\n";
        return;
    }

    // ---- Open ROOT file
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: could not open " << filename << "\n";
        if (file)
        {
            file->Close();
            delete file;
        }
        return;
    }

    TTree *tree = (TTree *)file->Get("events");
    if (!tree)
    {
        std::cerr << "Error: no TTree 'events' in " << filename << "\n";
        file->Close();
        delete file;
        return;
    }

    // ---- Get adc/tot shape dynamically
    TBranch *br_adc = tree->GetBranch("adc");
    if (!br_adc)
    {
        std::cerr << "Error: missing branch 'adc'\n";
        file->Close();
        delete file;
        return;
    }
    TLeaf *leaf_adc = br_adc->GetLeaf("adc");
    if (!leaf_adc)
    {
        std::cerr << "Error: missing leaf 'adc'\n";
        file->Close();
        delete file;
        return;
    }

    TBranch *br_tot = tree->GetBranch("tot");
    TLeaf *leaf_tot = nullptr;
    if (use_hybrid_tot)
    {
        if (!br_tot)
        {
            std::cerr << "Warning: use_hybrid_tot=true but branch 'tot' missing. Falling back to ADC-only.\n";
            use_hybrid_tot = false;
        }
        else
        {
            leaf_tot = br_tot->GetLeaf("tot");
            if (!leaf_tot)
            {
                std::cerr << "Warning: use_hybrid_tot=true but leaf 'tot' missing. Falling back to ADC-only.\n";
                use_hybrid_tot = false;
            }
        }
    }

    const int n_adc = leaf_adc->GetLen();
    if ((n_adc % LED_SAMPLES_PER_CHANNEL) != 0)
    {
        std::cerr << "Error: unexpected adc length n_adc=" << n_adc
                  << " (not divisible by " << LED_SAMPLES_PER_CHANNEL << ")\n";
        file->Close();
        delete file;
        return;
    }

    int n_tot = 0;
    if (use_hybrid_tot)
    {
        n_tot = leaf_tot->GetLen();
        if (n_tot != n_adc)
        {
            std::cerr << "Warning: tot length differs from adc (n_tot=" << n_tot << " n_adc=" << n_adc
                      << "). Falling back to ADC-only.\n";
            use_hybrid_tot = false;
        }
    }

    const int nch_from_file = n_adc / LED_SAMPLES_PER_CHANNEL;

    std::vector<uint32_t> adc_buf((size_t)n_adc);
    std::vector<uint32_t> tot_buf;
    if (use_hybrid_tot)
        tot_buf.resize((size_t)n_adc);

    tree->SetBranchAddress("adc", adc_buf.data());
    if (use_hybrid_tot)
        tree->SetBranchAddress("tot", tot_buf.data());

    auto adc_at = [&](int ch, int t) -> uint32_t &
    {
        return adc_buf[(size_t)ch * LED_SAMPLES_PER_CHANNEL + (size_t)t];
    };
    auto tot_at = [&](int ch, int t) -> uint32_t &
    {
        return tot_buf[(size_t)ch * LED_SAMPLES_PER_CHANNEL + (size_t)t];
    };

    // ---- Build reverse map: channel -> (crystal, sipm)
    auto reverse = build_reverse_channel_map(mapping, LED_SIPMS_PER_CRYSTAL);

    // ---- Create per-(crystal,sipm) histograms only for crystals in mapping
    std::map<int, std::array<TH1F *, LED_SIPMS_PER_CRYSTAL>> h_adc;
    std::map<int, std::array<TH1F *, LED_SIPMS_PER_CRYSTAL>> h_tot;

    for (const auto &kv : mapping)
    {
        const int crystal_id = kv.first;
        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
        {
            h_adc[crystal_id][sipm] = new TH1F(
                Form("h_adc_%.3fV_cr%d_sipm%d", voltage, crystal_id, sipm),
                Form("ADC-only | %.3f V | crystal %d | sipm %d;Signal;Counts", voltage, crystal_id, sipm),
                200, 0, 1024);

            h_tot[crystal_id][sipm] = new TH1F(
                Form("h_tot_%.3fV_cr%d_sipm%d", voltage, crystal_id, sipm),
                Form("ToT-used | %.3f V | crystal %d | sipm %d;Signal;Counts", voltage, crystal_id, sipm),
                200, 0, 1024);
        }
    }

    // ---- Active channels (like led_analysis)
    // Use max_channels = nch_from_file so it adapts to 1 KCU vs many
    auto active_channels = get_active_channels_from_mapping(mapping, LED_SIPMS_PER_CRYSTAL, LED_MAX_NUM_CRYSTALS, nch_from_file);

    // ---- Event loop: loop active channels, route via reverse map
    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t ev = 0; ev < nEntries; ++ev)
    {
        tree->GetEntry(ev);

        for (int ch : active_channels)
        {
            if (ch < 0 || ch >= nch_from_file)
                continue;

            auto it = reverse.find(ch);
            if (it == reverse.end())
                continue;

            const int crystal_id = it->second.first;
            const int sipm = it->second.second;

            if (!use_selected_sipm(crystal_id, sipm))
                continue;

            float sig = 0.0f;
            bool used_tot = false;

            if (use_hybrid_tot)
            {
                // Build pointers to the start of this channel waveform (20 samples)
                uint32_t *adc_ptr = &adc_at(ch, 0);
                uint32_t *tot_ptr = &tot_at(ch, 0);

                calculate_signal_hybrid(adc_ptr, tot_ptr, 1.0f, sig, used_tot);

                if (used_tot)
                    h_tot[crystal_id][sipm]->Fill(sig);
                else
                    h_adc[crystal_id][sipm]->Fill(sig);
            }
            else
            {
                // Pure ADC mode: everything is ADC
                uint32_t *adc_ptr = &adc_at(ch, 0);

                sig = calculate_signal_adc(adc_ptr, 1.0f);

                // In ADC-only mode
                h_adc[crystal_id][sipm]->Fill(sig);
            }
        }
    }

    // ---- Fit peaks, compute gain factors
    // We'll store results in arrays indexed by (crystal_id,sipm), and also in a "global index" = crystal*16+sipm for legacy plots.
    std::map<int, std::array<float, LED_SIPMS_PER_CRYSTAL>> peak;
    std::map<int, std::array<float, LED_SIPMS_PER_CRYSTAL>> sigma;

    for (const auto &kv : mapping)
    {
        const int crystal_id = kv.first;
        for (int sipm_i = 0; sipm_i < LED_SIPMS_PER_CRYSTAL; ++sipm_i)
        {
            float pk = 0, sg = 0;
            TH1F *hist = h_adc[crystal_id][sipm_i];
            bool ok = fit_peak_crystalball(hist, pk, sg);
            peak[crystal_id][sipm_i] = ok ? pk : 0.0f;
            sigma[crystal_id][sipm_i] = ok ? sg : 0.0f;
        }
    }

    // Crystal mean peak
    std::map<int, float> crystal_mean;
    for (const auto &kv : mapping)
    {
        const int crystal_id = kv.first;
        float sum = 0.0f;
        int n = 0;

        for (int sipm_i = 0; sipm_i < LED_SIPMS_PER_CRYSTAL; ++sipm_i)
        {
            if (!use_selected_sipm(crystal_id, sipm_i))
                continue;
            if (exclude_from_crystal_mean(crystal_id, sipm_i))
                continue;

            const float pk = peak[crystal_id][sipm_i];
            if (pk <= 0)
                continue;

            sum += pk;
            n++;
        }

        crystal_mean[crystal_id] = (n > 0) ? (sum / n) : 0.0f;
        std::cout << Form("Crystal %d mean peak: %.1f\n", crystal_id, crystal_mean[crystal_id]);
    }

    // Channel gain factors: crystal_mean / channel_peak
    std::map<int, std::array<float, LED_SIPMS_PER_CRYSTAL>> gain_factor;
    for (const auto &kv : mapping)
    {
        const int crystal_id = kv.first;
        for (int sipm_i = 0; sipm_i < LED_SIPMS_PER_CRYSTAL; ++sipm_i)
        {
            const float pk = peak[crystal_id][sipm_i];
            const float cm = crystal_mean[crystal_id];
            if (pk > 0 && cm > 0)
                gain_factor[crystal_id][sipm_i] = cm / pk;
            else
                gain_factor[crystal_id][sipm_i] = 1.0f;
        }
    }

    // Optional: per-crystal equalization
    float mean_across_crystals = 0.0f;
    int ncr_good = 0;
    for (const auto &kv : crystal_mean)
    {
        const int crystal_id = kv.first;
        const float cm = kv.second;
        if (cm <= 0)
            continue;
        if (crystal_id == 9)
            continue; // keep your "bad crystal" exclusion
        mean_across_crystals += cm;
        ncr_good++;
    }
    if (ncr_good > 0)
        mean_across_crystals /= ncr_good;

    std::map<int, float> crystal_gain;
    for (const auto &kv : crystal_mean)
    {
        const int crystal_id = kv.first;
        const float cm = kv.second;
        crystal_gain[crystal_id] = (cm > 0) ? (mean_across_crystals / cm) : 1.0f;
    }

    // lambda function Force dead channels to 0 (same as original) IF those crystals exist in mapping
    auto zero_if_present = [&](int crystal_id, int sipm_i)
    {
        auto it = gain_factor.find(crystal_id);
        if (it == gain_factor.end())
            return;
        if (sipm_i < 0 || sipm_i >= LED_SIPMS_PER_CRYSTAL)
            return;
        it->second[sipm_i] = 0.0f;
    };

    zero_if_present(15, 10);
    zero_if_present(21, 1);
    zero_if_present(21, 4);
    zero_if_present(21, 10);
    zero_if_present(21, 11);
    zero_if_present(21, 13);
    zero_if_present(21, 14);
    zero_if_present(21, 15);

    // ---- Output names: unique per run, and compatible with adc_calibration_rf_cc.cxx
    const int run_out = extract_run_number(filename);

    std::string tag = std::string(Form("Run%03d_%.3fV", run_out, voltage));
    if (cc >= 0)
        tag += Form("_CC%d", cc);
    if (rf >= 0)
        tag += Form("_Rf%d", rf);

    std::string pdf_file = std::string(Form("%s/gain_match_%s.pdf", outdir, tag.c_str()));
    std::string root_file = std::string(Form("%s/gain_match_Run%03d.root", outdir, run_out));

    // ---- Save PDF: one page per crystal in mapping (4x4 SiPMs)
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.04);
    text.SetTextFont(42);

    TCanvas *canvas = new TCanvas("gain_matching", "", 900, 700);

    canvas->SaveAs((pdf_file + "[").c_str());

    for (const auto &kv : mapping)
    {
        const int crystal_id = kv.first;

        canvas->Clear();
        canvas->Divide(4, 4);

        for (int sipm_i = 0; sipm_i < LED_SIPMS_PER_CRYSTAL; ++sipm_i)
        {
            canvas->cd(sipm_i + 1);
            TH1F *hist = h_adc[crystal_id][sipm_i];
            if (hist)
                hist->Draw();

            text.DrawLatex(0.15, 0.83, Form("%.3f V", voltage));
            text.DrawLatex(0.15, 0.78, Form("Crystal %d SiPM %d", crystal_id, sipm_i));
            text.DrawLatex(0.15, 0.73, Form("Peak: %.1f", peak[crystal_id][sipm_i]));
            text.DrawLatex(0.15, 0.68, Form("Ch gain: %.2f", gain_factor[crystal_id][sipm_i]));
            text.DrawLatex(0.15, 0.63, Form("Cr gain: %.2f", crystal_gain[crystal_id]));
        }

        canvas->SaveAs(pdf_file.c_str());
    }

    // ---- Summary gain histograms (legacy indexing: crystal*16 + sipm)
    // Keep 400 bins so it matches old plots; fill only what is present.
    TH1F *gain_hist = new TH1F("gain_factors", "Gain Factors;crystal*16+sipm;Gain Factor", 25 * 16, 0, 25 * 16);
    TH1F *crystal_gain_hist = new TH1F("crystal_factor", "Crystal Gain;Crystal;Gain Factor", 25, 0, 25);

    for (int cr = 0; cr < LED_MAX_NUM_CRYSTALS; ++cr)
        crystal_gain_hist->SetBinContent(cr + 1, 1.0f);

    for (const auto &kv : mapping)
    {
        const int crystal_id = kv.first;
        if (crystal_id >= 0 && crystal_id < LED_MAX_NUM_CRYSTALS)
            crystal_gain_hist->SetBinContent(crystal_id + 1, crystal_gain[crystal_id]);

        for (int sipm_i = 0; sipm_i < LED_SIPMS_PER_CRYSTAL; ++sipm_i)
        {
            const int idx = crystal_id * LED_SIPMS_PER_CRYSTAL + sipm_i;
            if (idx >= 0 && idx < LED_MAX_NUM_CRYSTALS * LED_SIPMS_PER_CRYSTAL)
                gain_hist->SetBinContent(idx + 1, gain_factor[crystal_id][sipm_i]);
        }
    }

    canvas->Clear();
    gain_hist->SetMinimum(0);
    gain_hist->SetMaximum(gain_hist->GetMaximum() * 1.2);

    gain_hist->Draw();
    canvas->SaveAs(pdf_file.c_str());

    canvas->Clear();
    crystal_gain_hist->SetMinimum(0);
    crystal_gain_hist->SetMaximum(crystal_gain_hist->GetMaximum() * 1.2);

    crystal_gain_hist->Draw();
    canvas->SaveAs(pdf_file.c_str());

    canvas->Clear();
    TLatex tmeta;
    tmeta.SetNDC();
    tmeta.SetTextFont(42);
    tmeta.SetTextSize(0.04);
    tmeta.DrawLatex(0.12, 0.82, Form("Run: %03d", run_out));
    tmeta.DrawLatex(0.12, 0.74, Form("Voltage label: %.3f V", voltage));
    if (cc >= 0)
        tmeta.DrawLatex(0.12, 0.66, Form("CC gain setting: %d", cc));
    if (rf >= 0)
        tmeta.DrawLatex(0.12, 0.58, Form("Rf setting: %d", rf));
    tmeta.DrawLatex(0.12, 0.50, Form("Hybrid ADC/ToT: %s", use_hybrid_tot ? "enabled" : "disabled"));
    canvas->SaveAs((pdf_file + "]").c_str());

    std::cout << "Saved PDF:  " << pdf_file << "\n";

    // ---- Save ROOT: gain histograms + a TTree with detailed per-channel info
    TFile *out = new TFile(root_file.c_str(), "RECREATE");
    gain_hist->Write();
    crystal_gain_hist->Write();

    TTree *t = new TTree("gain_table", "Gain match results");
    int t_crystal = 0, t_sipm = 0, t_channel = 0;
    float t_peak = 0, t_sigma = 0, t_cr_mean = 0, t_gain = 1, t_cr_gain = 1;

    t->Branch("crystal", &t_crystal, "crystal/I");
    t->Branch("sipm", &t_sipm, "sipm/I");
    t->Branch("channel", &t_channel, "channel/I");
    t->Branch("peak", &t_peak, "peak/F");
    t->Branch("sigma", &t_sigma, "sigma/F");
    t->Branch("crystal_mean", &t_cr_mean, "crystal_mean/F");
    t->Branch("gain", &t_gain, "gain/F");
    t->Branch("crystal_gain", &t_cr_gain, "crystal_gain/F");

    for (const auto &kv : mapping)
    {
        const int crystal_id = kv.first;
        auto chans = get_crystal_channels(mapping, crystal_id, LED_SIPMS_PER_CRYSTAL);

        for (int sipm_i = 0; sipm_i < LED_SIPMS_PER_CRYSTAL; ++sipm_i)
        {
            t_crystal = crystal_id;
            t_sipm = sipm_i;
            t_channel = (sipm_i < (int)chans.size()) ? chans[sipm_i] : -1;

            t_peak = peak[crystal_id][sipm_i];
            t_sigma = sigma[crystal_id][sipm_i];
            t_cr_mean = crystal_mean[crystal_id];
            t_gain = gain_factor[crystal_id][sipm_i];
            t_cr_gain = crystal_gain[crystal_id];

            t->Fill();
        }
    }

    t->Write();

    TParameter<int>("run_number", run_out).Write();
    TParameter<int>("cc_gain", cc).Write();
    TParameter<int>("rf_setting", rf).Write();
    TParameter<float>("label_voltage_V", voltage).Write();

    out->Close();
    delete out;

    std::cout << "Saved ROOT: " << root_file << "\n";

    // ---- Cleanup
    delete canvas;
    delete gain_hist;
    delete crystal_gain_hist;

    for (auto &kv : h_adc)
        for (TH1F *h : kv.second)
            delete h;
    for (auto &kv : h_tot)
        for (TH1F *h : kv.second)
            delete h;

    file->Close();
    delete file;
}

// ------------------------------------------------------------
// Scan like led_scan list (run, voltage)
// (uses data/Run%03d.root)
// ------------------------------------------------------------

// ------------------------------------------------------------
// Scan for the RF/CC campaign
// ------------------------------------------------------------
struct LedGainRunRFCC
{
    int run = -1;
    float voltage = 0.0f;
    int cc = -1;
    int rf = -1;
};

void gain_scan_rf_cc(const char *data_dir = "data",
                     const char *mapping_csv = "eeemcal_desy_dec2025_mapping_v2.csv",
                     const char *outdir = "outputs_rf_cc",
                     bool use_hybrid_tot = false)
{
    std::vector<LedGainRunRFCC> runs = {
        {381, 1.26f, 12, 3},
        {382, 1.26f, 12, 1},
        {383, 1.26f, 12, 2},
        {385, 1.26f, 12, 4},
        {386, 1.26f, 12, 5},
        {387, 1.26f, 12, 6},
        {388, 1.26f, 12, 7},
        {389, 1.26f, 12, 8},
        {390, 1.26f, 12, 9},
        {391, 1.26f, 12, 10},
        {392, 1.26f, 12, 11},
        {394, 1.26f, 12, 12},
        {395, 1.26f, 12, 13},
        {396, 1.26f, 12, 14},
        {397, 1.26f, 12, 15},
        {398, 0.0f, 12, 1},
        {399, 0.0f, 12, 2},
        {400, 0.0f, 12, 3},
        {401, 0.0f, 12, 4},
        {402, 0.0f, 12, 5},
        {404, 0.0f, 12, 6},
        {405, 0.0f, 12, 7},
        {406, 0.0f, 12, 8},
        {407, 0.0f, 12, 9},
        {408, 0.0f, 12, 10},
        {409, 0.0f, 12, 11},
        {410, 0.0f, 12, 12},
        {411, 0.0f, 12, 13},
        {412, 0.0f, 12, 14},
        {413, 0.0f, 12, 15},
        {414, 1.26f, 1, 3},
        {415, 0.0f, 1, 3},
        {416, 1.26f, 2, 3},
        {417, 0.0f, 2, 3},
        {418, 1.26f, 3, 3},
        {419, 0.0f, 3, 3},
        {420, 1.26f, 4, 3},
        {421, 0.0f, 4, 3},
        {422, 1.26f, 5, 3},
        {423, 0.0f, 5, 3},
        {424, 1.26f, 6, 3},
        {429, 0.0f, 6, 3},
        {430, 1.26f, 7, 3},
        {431, 0.0f, 7, 3},
        {432, 1.26f, 8, 3},
        {433, 0.0f, 8, 3},
        {434, 1.26f, 9, 3},
        {435, 0.0f, 9, 3},
        {436, 1.26f, 10, 3},
        {437, 0.0f, 10, 3},
        {438, 1.26f, 11, 3},
        {439, 0.0f, 11, 3},
        {440, 1.26f, 12, 3},
        {441, 0.0f, 12, 3},
        {442, 1.26f, 13, 3},
        {443, 0.0f, 13, 3},
        {444, 1.26f, 14, 3},
        {445, 0.0f, 14, 3},
        {446, 1.26f, 15, 3},
        {447, 0.0f, 15, 3}
    };

    for (const auto &rv : runs)
    {
        const std::string fn = std::string(Form("%s/Run%03d.root", data_dir, rv.run));
        std::cout << "\n--- Gain match RF/CC: " << fn
                  << " @ V=" << rv.voltage
                  << "  CC=" << rv.cc
                  << "  Rf=" << rv.rf
                  << " ---\n";

        gain_match_rf_cc_one(fn.c_str(),
                             rv.voltage,
                             mapping_csv,
                             outdir,
                             use_hybrid_tot,
                             rv.cc,
                             rv.rf);
    }
}
