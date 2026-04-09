// adc_calibration_rf_cc.cxx
//
// Based on adc_calibration.cxx, but focused on RF+CC data with the following differences:
// - Uses common_led.{h,cxx}
// - Uses mapping CSV from generate_mapping.py (FPGA,ASIC,Connector,Crystal,SiPM,Channel)
// - Builds reverse[channel] -> (crystal, sipm)
// - Loops over active channels and routes signals like led_analysis
// - Reads adc/tot/toa as flat buffers using TLeaf::GetLen()
// - Keeps adc_calibration purpose: energy hists, COG, pedestals, energy shares,
//   and produces an "ADC per reference unit" calibration.
// - Takes additional parameters for RF and CC selection, and only processes that subset of events.
// - Skips the LED wavelength scan and related parameters, since those are not relevant for RF+CC.
//
// Reference unit: voltage = 1.27 V (instead of 1 GeV)
// TO RUN automatically on a list of runs/voltages/CC-gains/Rf, use adc_calib_rf_cc_scan() and for cm and not cm:
// root -l -b
// .L common_led.cxx+
// .L adc_calibration_rf_cc.cxx+
// adc_calib_rf_cc_scan("data","eeemcal_desy_dec2025_mapping_v2.csv","outputs_rf_cc","outputs_rf_cc",true,false,false);
//
// Run examples:
//
// root -l -b
// .L common_led.cxx+
// .L adc_calibration_rf_cc.cxx+
// adc_calib_rf_cc_scan("data",
//                      "eeemcal_desy_dec2025_mapping_v2.csv",
//                      "outputs_rf_cc",
//                      "outputs",
//                      true,   // use_hybrid_tot
//                      false,   // use_common_mode
//                      false); // false = use fitted sigma in summary, true = use RMS
//

#include "common_led.h"

#include <map>
#include <array>
#include <vector>
#include <unordered_map>
#include <utility>
#include <string>
#include <iostream>
#include <fstream>
#include <tuple>
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
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TParameter.h>
#include <TError.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TLeaf.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TLegend.h>

static constexpr double WAVELENGTH_NM = 470.0; // for photon energy conversion

// Photon energy helpers (for bookkeeping / optional conversions)
// E_photon [eV] = (h*c)/lambda, with (h*c) = 1239.841984 eV*nm

static constexpr double HC_EV_NM = 1239.841984;
// static constexpr double EV_TO_J = 1.602176634e-19;
static inline double photon_energy_eV(double wavelength_nm)
{
    return (wavelength_nm > 0.0) ? (HC_EV_NM / wavelength_nm) : (HC_EV_NM / WAVELENGTH_NM);
}

// ---------------------- small helpers ----------------------

static void autoset_xrange_1d(TH1 *h, double frac_of_max = 0.01, int pad_bins = 2)
{
    if (!h || h->GetEntries() <= 0)
        return;

    const int nb = h->GetNbinsX();
    const double thr = std::max(1.0, frac_of_max * h->GetMaximum());

    int first = -1;
    int last  = -1;

    for (int b = 1; b <= nb; ++b)
    {
        if (h->GetBinContent(b) > thr)
        {
            first = b;
            break;
        }
    }

    for (int b = nb; b >= 1; --b)
    {
        if (h->GetBinContent(b) > thr)
        {
            last = b;
            break;
        }
    }

    if (first < 0 || last < 0 || first > last)
        return;

    first = std::max(1, first - pad_bins);
    last  = std::min(nb, last + pad_bins);

    const double xmin = h->GetXaxis()->GetBinLowEdge(first);
    const double xmax = h->GetXaxis()->GetBinUpEdge(last);
    h->GetXaxis()->SetRangeUser(xmin, xmax);
}

static void autoset_xy_2d(TH2 *h, double frac_of_max = 0.01, int pad_bins = 1)
{
    if (!h || h->GetEntries() <= 0)
        return;

    const int nx = h->GetNbinsX();
    const int ny = h->GetNbinsY();
    const double thr = std::max(1.0, frac_of_max * h->GetMaximum());

    int firstx = -1, lastx = -1, firsty = -1, lasty = -1;

    for (int bx = 1; bx <= nx; ++bx)
    {
        for (int by = 1; by <= ny; ++by)
        {
            if (h->GetBinContent(bx, by) > thr)
            {
                if (firstx < 0 || bx < firstx) firstx = bx;
                if (lastx  < 0 || bx > lastx ) lastx  = bx;
                if (firsty < 0 || by < firsty) firsty = by;
                if (lasty  < 0 || by > lasty ) lasty  = by;
            }
        }
    }

    if (firstx < 0 || lastx < 0 || firsty < 0 || lasty < 0)
        return;

    firstx = std::max(1, firstx - pad_bins);
    lastx  = std::min(nx, lastx + pad_bins);
    firsty = std::max(1, firsty - pad_bins);
    lasty  = std::min(ny, lasty + pad_bins);

    h->GetXaxis()->SetRangeUser(h->GetXaxis()->GetBinLowEdge(firstx),
                                h->GetXaxis()->GetBinUpEdge(lastx));
    h->GetYaxis()->SetRangeUser(h->GetYaxis()->GetBinLowEdge(firsty),
                                h->GetYaxis()->GetBinUpEdge(lasty));
}

static int32_t get_toa_first_nonzero(uint32_t *toa_values)
{
    for (int i = 0; i < LED_SAMPLES_PER_CHANNEL; ++i)
    {
        if (toa_values[i] > 0)
            return (int32_t)toa_values[i];
    }
    return -1;
}

static bool is_tot_event(uint32_t *tot_values)
{
    for (int i = 0; i < LED_SAMPLES_PER_CHANNEL; ++i)
    {
        if (tot_values[i] > (uint32_t)g_tot_min)
            return true;
    }
    return false;
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

    return point_in_ellipse(x_cog, y_cog, cx, cy, sx, sy);
}

static bool extract_voltage_from_calib_filename(const std::string &name, double &voltage)
{
    const std::string prefix = "adc_to_ref_calibration_";
    const std::string suffix = "V.root";

    if (name.size() <= prefix.size() + suffix.size())
        return false;
    if (name.rfind(prefix, 0) != 0)
        return false;
    if (name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0)
        return false;

    const std::string value = name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());
    try
    {
        voltage = std::stod(value);
        return true;
    }
    catch (...)
    {
        return false;
    }
}

static void make_adc_calibration_summary(const char *outdir)
{
    if (!outdir || !outdir[0])
        return;

    void *dirp = gSystem->OpenDirectory(outdir);
    if (!dirp)
    {
        std::cerr << "Warning: cannot open output directory for summary: " << outdir << "\n";
        return;
    }

    std::vector<std::tuple<double, double, double, double>> rows;

    const char *entry = nullptr;
    while ((entry = gSystem->GetDirEntry(dirp)) != nullptr)
    {
        std::string name(entry);
        double voltage = 0.0;
        if (!extract_voltage_from_calib_filename(name, voltage))
            continue;

        std::string path = std::string(outdir) + "/" + name;
        TFile *in = TFile::Open(path.c_str(), "READ");
        if (!in || in->IsZombie())
        {
            if (in)
            {
                in->Close();
                delete in;
            }
            continue;
        }

        auto *p_mean = dynamic_cast<TParameter<float> *>(in->Get(Form("mean_adc_to_ref_calibration_%.3fV", voltage)));
        auto *p_res = dynamic_cast<TParameter<float> *>(in->Get(Form("central_energy_res_pct_%.3fV", voltage)));
        auto *p_res_err = dynamic_cast<TParameter<float> *>(in->Get(Form("central_energy_res_pct_err_%.3fV", voltage)));

        if (p_mean && p_res)
        {
            const double res_err = p_res_err ? (double)p_res_err->GetVal() : 0.0;
            rows.emplace_back(voltage, (double)p_mean->GetVal(), (double)p_res->GetVal(), res_err);
        }

        in->Close();
        delete in;
    }
    gSystem->FreeDirectory(dirp);

    if (rows.empty())
    {
        std::cerr << "Warning: no adc_to_ref_calibration_*.root files found in " << outdir << "\n";
        return;
    }

    std::sort(rows.begin(), rows.end(),
              [](const auto &a, const auto &b)
              { return std::get<0>(a) < std::get<0>(b); });

    const std::string txt_out = std::string(outdir) + "/adc_to_ref_calibration_summary.txt";
    std::ofstream ofs(txt_out.c_str());
    ofs << "# voltage[V] mean_adc_to_ref_calibration central_energy_res_pct central_energy_res_pct_err\n";
    for (const auto &row : rows)
    {
        ofs << Form("%.3f %.6f %.6f %.6f", std::get<0>(row), std::get<1>(row), std::get<2>(row), std::get<3>(row)) << "\n";
    }
    ofs.close();

    const int n = (int)rows.size();
    TGraphErrors *gr = new TGraphErrors(n);
    gr->SetName("central_energy_res_pct_vs_voltage");
    gr->SetTitle("Central Energy Resolution vs Voltage;Voltage (V);Central Energy Resolution (%)");

    for (int i = 0; i < n; ++i)
    {
        gr->SetPoint(i, std::get<0>(rows[i]), std::get<2>(rows[i]));
        gr->SetPointError(i, 0.0, std::get<3>(rows[i]));
    }

    TCanvas *c = new TCanvas("c_adc_calib_summary", "", 900, 700);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.0);
    gr->SetLineWidth(2);
    c->SetGrid();
    gr->Draw("AP");

    const std::string png_out = std::string(outdir) + "/central_energy_res_pct_vs_voltage.png";
    const std::string pdf_out = std::string(outdir) + "/central_energy_res_pct_vs_voltage.pdf";
    c->SaveAs(png_out.c_str());
    c->SaveAs(pdf_out.c_str());

    const std::string root_out = std::string(outdir) + "/adc_to_ref_calibration_summary.root";
    TFile *fout = TFile::Open(root_out.c_str(), "RECREATE");
    gr->Write();
    fout->Close();
    delete fout;

    delete c;
    delete gr;

    std::cout << "Saved calibration summary: " << txt_out << "\n";
    std::cout << "Saved resolution-vs-voltage plots: " << png_out << " and " << pdf_out << "\n";
    std::cout << "Saved summary ROOT graph: " << root_out << "\n";
}

// ---------------------- main function ----------------------

void adc_calibration_rf_cc_one(const char *filename,
                               float voltage,
                               const char *mapping_csv = "eeemcal_desy_dec2025_mapping_v2.csv",
                               const char *gain_root_file = "outputs/gain_match_1.25V.root",
                               const char *outdir = "outputs",
                               Long64_t max_events = 1000000,
                               float energy_fraction_cut = 0.0f,
                               bool use_hybrid_tot = true,
                               bool use_common_mode = true,
                               float led_wavelength_nm = 470.0f,
                               int cc = -1,
                               int rf = -1)
{
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;
    gSystem->mkdir(outdir, true);

    // Use same signal config style as your other scripts
    g_signal_method = 3;
    g_tot_min = 50;

    const float REF_VOLTAGE = voltage;

    const float LED_WAVELENGTH_NM = (led_wavelength_nm > 0.0f) ? led_wavelength_nm : 470.0f;
    const double EPHOTON_EV = photon_energy_eV((double)LED_WAVELENGTH_NM);
    // const double EPHOTON_J  = EPHOTON_EV * EV_TO_J;

    // Per-crystal common-mode reference channels
    // int common_mode_channels_by_crystal[LED_MAX_NUM_CRYSTALS] = {
    //     224, 296, 188, 386, 152,
    //     512, 548, 242, 134, 404,
    //     206, 368, 332, 458, 8,
    //     440, 116, 350, 62, 98,
    //     422, 170, 260, 314, 278}; // used in DESY

    int common_mode_channels_by_crystal[LED_MAX_NUM_CRYSTALS];
    for (int cr = 0; cr < LED_MAX_NUM_CRYSTALS; ++cr)
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
    auto mapping = read_mapping_csv(mapping_csv, LED_SIPMS_PER_CRYSTAL);
    if (mapping.empty())
    {
        std::cerr << "ERROR: mapping empty. Check: " << mapping_csv << "\n";
        return;
    }

    // reverse[channel] = (crystal, sipm)
    auto reverse = build_reverse_channel_map(mapping, LED_SIPMS_PER_CRYSTAL);

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
        gain_factor = new TH1F("gain_factors_fallback", "gain_factors_fallback", LED_MAX_NUM_CRYSTALS * LED_SIPMS_PER_CRYSTAL, 0, LED_MAX_NUM_CRYSTALS * LED_SIPMS_PER_CRYSTAL);
        for (int i = 1; i <= LED_MAX_NUM_CRYSTALS * LED_SIPMS_PER_CRYSTAL; ++i)
            gain_factor->SetBinContent(i, 1.0);
    }
    if (!crystal_gain_factor)
    {
        crystal_gain_factor = new TH1F("crystal_factor_fallback", "crystal_factor_fallback", LED_MAX_NUM_CRYSTALS, 0, LED_MAX_NUM_CRYSTALS);
        for (int i = 1; i <= LED_MAX_NUM_CRYSTALS; ++i)
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
    if (n_adc <= 0 || (n_adc % LED_SAMPLES_PER_CHANNEL) != 0)
    {
        std::cerr << "ERROR: unexpected adc length " << n_adc
                  << " (not divisible by " << LED_SAMPLES_PER_CHANNEL << ")\n";
        f->Close();
        delete f;
        if (gf)
        {
            gf->Close();
            delete gf;
        }
        return;
    }

    const int n_channels = n_adc / LED_SAMPLES_PER_CHANNEL;

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
    { return &adc_buf[(size_t)ch * LED_SAMPLES_PER_CHANNEL]; };
    auto tot_ptr = [&](int ch) -> uint32_t *
    { return &tot_buf[(size_t)ch * LED_SAMPLES_PER_CHANNEL]; };
    auto toa_ptr = [&](int ch) -> uint32_t *
    { return &toa_buf[(size_t)ch * LED_SAMPLES_PER_CHANNEL]; };

    // Active channels based on mapping, bounded by file channels
    auto active_channels = get_active_channels_from_mapping(mapping, LED_SIPMS_PER_CRYSTAL, LED_MAX_NUM_CRYSTALS, n_channels);

    // ---------------- Histograms (global) ----------------
    TH1 *central_crystal_energy = new TH1F("central_crystal_energy",
                                           Form("Central Crystal Energy (%.3f V);Energy (ADC);Events", voltage),
                                           500, 0, 20000);

    TH1 *central_nine_energy = new TH1F("central_nine_energy",
                                        Form("Central 3x3 Energy (%.3f V);Energy (ADC);Events", voltage),
                                        500, 0, 20000);

    TH1 *total_energy = new TH1F("total_energy",
                                 Form("Total Energy (%.3f V);Energy (ADC);Events", voltage),
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
    E_vs_toa.reserve(LED_SIPMS_PER_CRYSTAL);
    for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; sipm++)
    {
        E_vs_toa.push_back(new TH2F(Form("E_vs_toa_sipm_%d", sipm),
                                    Form("Central Crystal SiPM %d: Energy vs ToA;ToA;Energy (ADC)", sipm),
                                    1024, 0, 1024, 500, 0, 3000));
    }

    // ToA correlations (central crystal only)
    // std::vector<std::vector<TH2 *>> toa_correlations(LED_SIPMS_PER_CRYSTAL);
    // for (int i = 0; i < LED_SIPMS_PER_CRYSTAL; i++)
    // {
    //     toa_correlations[i].reserve(LED_SIPMS_PER_CRYSTAL);
    //     for (int j = 0; j < LED_SIPMS_PER_CRYSTAL; j++)
    //     {
    //         toa_correlations[i].push_back(new TH2F(Form("toa_correlation_sipm_%02d_sipm_%02d", i, j),
    //                                                Form("Central crystal: SiPM %02d vs %02d ToA;SiPM %02d ToA;SiPM %02d ToA", i, j, i, j),
    //                                                1024, 0, 1024, 1024, 0, 1024));
    //     }
    // }

    // ---------------- per-crystal/per-sipm histograms ----------------
    // Allocate for 25 crystals, but only fill when that crystal exists in mapping.
    std::vector<std::array<TH1F *, LED_SIPMS_PER_CRYSTAL>> sipm_energy(LED_MAX_NUM_CRYSTALS);
    std::vector<std::array<TH2F *, LED_SIPMS_PER_CRYSTAL>> sipm_waveform(LED_MAX_NUM_CRYSTALS);
    std::vector<TH1F *> crystal_energy(LED_MAX_NUM_CRYSTALS, nullptr);
    std::vector<TH1F *> crystal_energy_shares(LED_MAX_NUM_CRYSTALS, nullptr);
    std::vector<TH2F *> crystal_common_mode(LED_MAX_NUM_CRYSTALS, nullptr);

    for (int cr = 0; cr < LED_MAX_NUM_CRYSTALS; ++cr)
    {
        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
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
        if ((int)(entry * LED_MAX_NUM_CRYSTALS / nentries) > complete)
        {
            complete = (int)(entry * LED_MAX_NUM_CRYSTALS / nentries);
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
            // use central crystal sipm as reference ToA channel
            int sipm_ref = first_selected_sipm_in_mapping(mapping, central_crystal_id, LED_SIPMS_PER_CRYSTAL);
            if (sipm_ref >= 0)
            {
                int ch0 = get_crystal_channels(mapping, central_crystal_id, LED_SIPMS_PER_CRYSTAL)[sipm_ref];

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

            if (cr < 0 || cr >= LED_MAX_NUM_CRYSTALS || sipm < 0 || sipm >= LED_SIPMS_PER_CRYSTAL)
                continue;
            if (!mapping.count(cr))
                continue; // ignore crystals not present

            if (!use_selected_sipm(cr, sipm))
                continue;

            // gain
            float gain = gain_factor->GetBinContent(cr * LED_SIPMS_PER_CRYSTAL + sipm + 1);

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
                        cm = (cm_adc[0] /*+ cm_adc[1] + cm_adc[2]*/) /* / 2.0f */;
                        cm_cache[cr] = cm;
                        cm_ok[cr] = true;

                        // record CM waveform
                        for (int s = 0; s < LED_SAMPLES_PER_CHANNEL; ++s)
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
            uint32_t adc_corr[LED_SAMPLES_PER_CHANNEL];
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

            for (int s = 0; s < LED_SAMPLES_PER_CHANNEL; ++s)
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
            // pedestals->Fill(ch, adc_used[1]);
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
                    for (int s = 0; s < LED_SAMPLES_PER_CHANNEL; ++s)
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
                    auto chans = get_crystal_channels(mapping, central_crystal_id, LED_SIPMS_PER_CRYSTAL);
                    // for (int other = 0; other < LED_SIPMS_PER_CRYSTAL; ++other)
                    // {
                    //     if (!use_selected_sipm(central_crystal_id, other))
                    //         continue;
                    //     if (other == sipm)
                    //         continue;
                    //     int ch_other = chans[other];
                    //     if (ch_other < 0 || ch_other >= n_channels)
                    //         continue;
                    //     int32_t toa_other = get_toa_first_nonzero(toa_ptr(ch_other));
                    //     if (toa_other < 0)
                    //         continue;
                    //     toa_correlations[sipm][other]->Fill(toa_val, toa_other);
                    // }
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
        // keep &= (!is_tot); // same as your old keep &= (!is_tot_event)
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

        for (int i = 0; i < LED_SIPMS_PER_CRYSTAL; i++)
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
            for (int cr = 0; cr < LED_MAX_NUM_CRYSTALS; cr++)
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

    print_progress_25(LED_MAX_NUM_CRYSTALS);

    std::cout << "Total ToT-tagged events: " << tot_events << " out of " << nentries << "\n";

    // ---------------- calibration: ADC per reference unit (1.27 V or 1.25 V) ----------------
    //   signal_for_ref = mean(crystal_energy) / mean(crystal_energy_share)
    // where "ref" is defined by voltage 1.27 V or 1.25 V.
    int crystal_mapping_for_plots[LED_MAX_NUM_CRYSTALS] = {
        4, 9, 14, 19, 24,
        3, 8, 13, 18, 23,
        2, 7, 12, 17, 22,
        1, 6, 11, 16, 21,
        0, 5, 10, 15, 20}; // crystal 9 excluded from calibration because it has weid behavior in many runs

    std::vector<float> ref_calib;
    float mean_calib = 0.0f;
    int n_used = 0;

    for (int i = 0; i < LED_MAX_NUM_CRYSTALS; i++)
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

    // ---------------- central-crystal energy resolution (at REF_VOLTAGE) ----------------
    double e_mu = 0.0;
    double e_sigma = 0.0;
    double e_res = 0.0; // sigma/mu
    double e_mu_err = 0.0;
    double e_sigma_err = 0.0;
    double e_res_err = 0.0;
    bool e_fit_ok = false;

    if (central_crystal_energy && central_crystal_energy->GetEntries() > 50 && central_crystal_energy->GetRMS() > 0)
    {
        const double xMin = central_crystal_energy->GetXaxis()->GetXmin();
        const double xMax = central_crystal_energy->GetXaxis()->GetXmax();

        // const double mean0 = central_crystal_energy->GetBinCenter(central_crystal_energy->GetMean());
        // const double rms0 = std::max(50.0, (double)central_crystal_energy->GetRMS());

        // double fmin = std::max(xMin, mean0 - 1.5 * rms0);
        // double fmax = std::min(xMax, mean0 + 1.5 * rms0);
        // if (fmax - fmin < 0.03 * (xMax - xMin))
        // {
        //     fmin = std::max(xMin, mean0 - 0.5 * (xMax - xMin));
        //     fmax = std::min(xMax, mean0 + 0.5 * (xMax - xMin));
        // }
        // const double mean0 = central_crystal_energy->GetBinCenter(central_crystal_energy->GetMean()); // central_crystal_energy->GetMaximumBin()
        const double mean0 = central_crystal_energy->GetMean();
        const double rms0 = std::max(1.0, (double)central_crystal_energy->GetRMS());
        const double fmin = std::max(0.0, mean0 - 1.5 * rms0);
        const double fmax = mean0 + 1.5 * rms0;

        // TF1 fE(Form("central_energy_gaus_run%03d_%.3fV", run, REF_VOLTAGE), "gaus", xMin, xMax);
        TF1 fE(Form("central_energy_gaus_run%03d_%.3fV", run, REF_VOLTAGE), "gaus", fmin, fmax);
        fE.SetParameters(central_crystal_energy->GetMean(), mean0, rms0); // (central_crystal_energy->GetMaximum(), mean0, rms0);

        TFitResultPtr rr = central_crystal_energy->Fit(&fE, "RQS"); // stores fit on hist

        if (rr.Get() && rr->Status() == 0)
        {
            e_mu = fE.GetParameter(1);
            e_sigma = std::abs(fE.GetParameter(2));
            e_mu_err = fE.GetParError(1);
            e_sigma_err = fE.GetParError(2);
            e_res = (e_mu != 0.0) ? (e_sigma / e_mu) : 0.0;
            if (e_mu != 0.0 && e_sigma != 0.0)
            {
                const double rel_sigma = e_sigma_err / e_sigma;
                const double rel_mu = e_mu_err / e_mu;
                e_res_err = std::abs(e_res) * std::sqrt(rel_sigma * rel_sigma + rel_mu * rel_mu);
            }
            e_fit_ok = true;
        }
    }

    std::cout << "Central crystal energy peak @ " << REF_VOLTAGE
              << " V: mean=" << e_mu
              << " sigma=" << e_sigma
              << " resolution=" << 100.0 * e_res << " %"
              << (e_fit_ok ? "" : " (fit failed)") << "\n";

    // ---------------- save outputs ----------------
    int run_out = extract_run_number(filename);

    std::string out_tag = std::string(Form("Run%03d_%.3fV", run_out, voltage));
    if (cc >= 0)
        out_tag += Form("_CC%d", cc);
    if (rf >= 0)
        out_tag += Form("_Rf%d", rf);

    std::string pdf = std::string(Form("%s/adc_calibration_%s.pdf", outdir, out_tag.c_str()));
    std::string root_out = std::string(Form("%s/adc_to_ref_calibration_%s.root", outdir, out_tag.c_str()));

    TCanvas *canvas = new TCanvas("adc_calibration_canvas", "", 900, 700);
    canvas->SetRightMargin(0.05);

    canvas->SaveAs((pdf + "[").c_str());

    autoset_xrange_1d(central_crystal_energy, 0.01, 3);
    central_crystal_energy->Draw("hist e");
    if (TF1 *f = central_crystal_energy->GetFunction(Form("central_energy_gaus_run%03d_%.3fV", run_out, REF_VOLTAGE)))
        f->Draw("same");
    canvas->SaveAs(pdf.c_str());

    autoset_xrange_1d(central_nine_energy, 0.01, 3);
    central_nine_energy->Draw("hist e");
    canvas->SaveAs(pdf.c_str());

    autoset_xrange_1d(total_energy, 0.01, 3);
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
    TH1F *pedestals_width = new TH1F("pedestals_width", "Pedestal Width (#sigma) vs Channel;Channel;Width #sigma (ADC)",
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
    for (int i = 0; i < LED_MAX_NUM_CRYSTALS; i++)
    {
        canvas->cd(i + 1);
        int cr = crystal_mapping_for_plots[i];
        autoset_xrange_1d(crystal_energy[cr], 0.01, 2);
        crystal_energy[cr]->Draw("hist e");
    }
    canvas->SaveAs(pdf.c_str());

    // crystal shares + calibration text
    canvas->Clear();
    canvas->Divide(5, 5);
    float mean_calib_check = 0.0f;
    int used_check = 0;

    for (int i = 0; i < LED_MAX_NUM_CRYSTALS; i++)
    {
        canvas->cd(i + 1);
        int cr = crystal_mapping_for_plots[i];

        crystal_energy_shares[cr]->SetLineColor(kRed);
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
            t.DrawLatex(x_coord, 0.82, Form("ADC @ %.3fV: %.1f", REF_VOLTAGE, adc_per_ref));
        }
    }
    if (used_check > 0)
        mean_calib_check /= used_check;

    canvas->SaveAs(pdf.c_str());
    gPad->SetLogx(0);

    // per-crystal pages: sipm energy + waveform + CM waveform
    for (int cr = 0; cr < LED_MAX_NUM_CRYSTALS; cr++)
    {
        if (!mapping.count(cr))
            continue;

        canvas->Clear();
        canvas->Divide(4, 4);
        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; sipm++)
        {
            canvas->cd(sipm + 1);
            autoset_xrange_1d(sipm_energy[cr][sipm], 0.01, 2);
            sipm_energy[cr][sipm]->Draw("hist e");
        }
        canvas->SaveAs(pdf.c_str());

        canvas->Clear();
        canvas->Divide(4, 4);
        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; sipm++)
        {
            canvas->cd(sipm + 1);
            autoset_xy_2d(sipm_waveform[cr][sipm], 0.01, 1);
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
    for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; sipm++)
    {
        canvas->cd(sipm + 1);
        autoset_xy_2d(E_vs_toa[sipm], 0.01, 1);
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
    // TCanvas *canvas2 = new TCanvas("toa_correlations_canvas", "", 8000, 8000);
    // canvas2->Divide(16, 16, 0.0005, 0.0005);
    // for (int i = 0; i < LED_SIPMS_PER_CRYSTAL; i++)
    // {
    //     for (int j = 0; j < LED_SIPMS_PER_CRYSTAL; j++)
    //     {
    //         if (j <= i)
    //         {
    //             canvas2->cd(i * LED_SIPMS_PER_CRYSTAL + j + 1);
    //             toa_correlations[i][j]->Draw("colz");
    //         }
    //     }
    // }
    // std::string corr_png = std::string(Form("%s/adc_toa_correlation_Run%03d.png", outdir, run_out));
    // canvas2->SaveAs(corr_png.c_str());
    // delete canvas2;

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
    t.DrawLatex(0.12, 0.84, Form("Run: %03d   Voltage label: %.3f V", run_out, voltage));
    if (cc >= 0)
        t.DrawLatex(0.12, 0.79, Form("CC gain setting: %d", cc));
    if (rf >= 0)
        t.DrawLatex(0.12, 0.74, Form("Rf setting: %d", rf));
    t.DrawLatex(0.12, 0.69, Form("Reference unit: %.3f V (defined as 1 unit)", REF_VOLTAGE));
    t.DrawLatex(0.12, 0.64, Form("Mean ADC per ref unit: %.3f", mean_calib));
    t.DrawLatex(0.12, 0.59, Form("LED wavelength: %.1f nm  (E_{#gamma}=%.3f eV)", LED_WAVELENGTH_NM, EPHOTON_EV));
    t.DrawLatex(0.12, 0.54, Form("ToT mode used: %s", use_hybrid_tot ? "hybrid ADC/ToT" : "ADC-only"));
    t.DrawLatex(0.12, 0.49, Form("Common-mode subtraction: %s", use_common_mode ? "enabled" : "disabled"));
    t.DrawLatex(0.12, 0.44, Form("Central energy res (sigma/mu): %.2f%%", 100.0 * e_res));

    canvas->SaveAs(pdf.c_str());

    canvas->SaveAs((pdf + "]").c_str());
    std::cout << "Saved PDF: " << pdf << "\n";
    // std::cout << "Saved ToA correlation PNG: " << corr_png << "\n";

    // ROOT output for the reference calibration
    TFile *out = TFile::Open(root_out.c_str(), "RECREATE");
    TParameter<float> p(Form("mean_adc_to_ref_calibration_%.3fV", REF_VOLTAGE), mean_calib);
    p.Write();

    // Resolution parameters for the reference voltage peak
    TParameter<float> p_mu(Form("central_energy_mu_%.3fV", REF_VOLTAGE), (float)e_mu);
    TParameter<float> p_sigma(Form("central_energy_sigma_%.3fV", REF_VOLTAGE), (float)e_sigma);
    TParameter<float> p_res(Form("central_energy_res_%.3fV", REF_VOLTAGE), (float)e_res);
    TParameter<float> p_mu_err(Form("central_energy_mu_err_%.3fV", REF_VOLTAGE), (float)e_mu_err);
    TParameter<float> p_sigma_err(Form("central_energy_sigma_err_%.3fV", REF_VOLTAGE), (float)e_sigma_err);
    TParameter<float> p_res_err(Form("central_energy_res_err_%.3fV", REF_VOLTAGE), (float)e_res_err);
    TParameter<float> p_res_pct(Form("central_energy_res_pct_%.3fV", REF_VOLTAGE), (float)(100.0 * e_res));
    TParameter<float> p_res_pct_err(Form("central_energy_res_pct_err_%.3fV", REF_VOLTAGE), (float)(100.0 * e_res_err));
    TParameter<int> p_fitok(Form("central_energy_fit_ok_%.3fV", REF_VOLTAGE), (int)e_fit_ok);

    p_mu.Write();
    p_sigma.Write();
    p_res.Write();
    p_mu_err.Write();
    p_sigma_err.Write();
    p_res_err.Write();
    p_res_pct.Write();
    p_res_pct_err.Write();
    p_fitok.Write();

    // Bookkeeping: LED wavelength and photon energy (useful if to obtain an optical-energy equivalent)
    TParameter<float> p_led_wavelength_nm("led_wavelength_nm", LED_WAVELENGTH_NM);
    TParameter<double> p_photon_energy_eV("photon_energy_eV", EPHOTON_EV);
    // TParameter<double> p_photon_energy_J("photon_energy_J", EPHOTON_J);
    p_led_wavelength_nm.Write();
    p_photon_energy_eV.Write();
    // p_photon_energy_J.Write();

    TParameter<int>("run_number", run_out).Write();
    TParameter<int>("cc_gain", cc).Write();
    TParameter<int>("rf_setting", rf).Write();
    TParameter<float>("label_voltage_V", voltage).Write();
    TParameter<float>("summary_mean", (float)e_mu).Write();
    TParameter<float>("summary_sigma_fit", (float)e_sigma).Write();
    TParameter<float>("summary_sigma_fit_err", (float)e_sigma_err).Write();
    TParameter<float>("summary_rms_raw", (float)central_crystal_energy->GetRMS()).Write();
    TParameter<float>("summary_rms_raw_err", (float)central_crystal_energy->GetRMSError()).Write();
    TParameter<float>("summary_resolution_pct_fit", (float)(100.0 * e_res)).Write();
    TParameter<float>("summary_resolution_pct_fit_err", (float)(100.0 * e_res_err)).Write();
    {
        const float mean_hist = (float)central_crystal_energy->GetMean();
        const float rms_hist = (float)central_crystal_energy->GetRMS();
        const float res_rms_pct = (mean_hist > 0.0f) ? (100.0f * rms_hist / mean_hist) : 0.0f;
        TParameter<float>("summary_mean_raw", mean_hist).Write();
        TParameter<float>("summary_mean_raw_err", (float)central_crystal_energy->GetMeanError()).Write();
        TParameter<float>("summary_resolution_pct_rms", res_rms_pct).Write();
    }

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

    for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
    {
        E_vs_toa[sipm]->Write();
    }

    // for (int i = 0; i < LED_SIPMS_PER_CRYSTAL; i++)
    // {
    //     for (int j = 0; j < LED_SIPMS_PER_CRYSTAL; j++)
    //     {
    //         toa_correlations[i][j]->Write();
    //     }
    // }

    for (int cr = 0; cr < LED_MAX_NUM_CRYSTALS; ++cr)
    {
        if (!mapping.count(cr))
            continue;
        crystal_energy[cr]->Write();
        crystal_energy_shares[cr]->Write();
        crystal_common_mode[cr]->Write();
        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
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

    // cleanup histograms allocated in this run
    delete central_crystal_energy;
    delete central_nine_energy;
    delete total_energy;
    delete cog_distribution;
    delete pedestals;
    delete pedestals_width;
    delete toa_distribution;
    delete toa_sample;
    delete max_sample_index;
    delete second_max_sample_index;
    delete third_max_sample_index;

    for (auto *h : E_vs_toa)
        delete h;

    // for (int i = 0; i < LED_SIPMS_PER_CRYSTAL; ++i)
    // {
    //     for (int j = 0; j < LED_SIPMS_PER_CRYSTAL; ++j)
    //         delete toa_correlations[i][j];
    // }

    for (int cr = 0; cr < LED_MAX_NUM_CRYSTALS; ++cr)
    {
        delete crystal_energy[cr];
        delete crystal_energy_shares[cr];
        delete crystal_common_mode[cr];

        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
        {
            delete sipm_energy[cr][sipm];
            delete sipm_waveform[cr][sipm];
        }
    }

    delete canvas;
}

// ------------------------------------------------------------
// Scan like led_scan list (run, voltage)
// (uses data/Run%03d.root)
// ------------------------------------------------------------

// ------------------------------------------------------------
// Run-by-run ADC calibration scan for datasets where voltage is
// not unique and extra hardware settings (CC, Rf) changed.
// It keeps the same per-run analysis pages as adc_calibration.cxx,
// but writes unique ROOT/PDF files per run and produces an
// additional summary PDF with:
//   page 1: mean + sigma(or RMS) vs Rf and CC
//   page 2: resolution vs Rf and CC
// ------------------------------------------------------------

struct LedRunCalibRFCC
{
    int run = -1;
    float voltage = 0.0f;
    int cc = -1;
    int rf = -1;
};

struct LedRunSummaryRFCC
{
    int run = -1;
    float voltage = 0.0f;
    int cc = -1;
    int rf = -1;
    double mean = 0.0;
    double emean = 0.0;
    double sigma = 0.0;
    double esigma = 0.0;
    double rms = 0.0;
    double erms = 0.0;
    double res_fit_pct = 0.0;
    double res_fit_pct_err = 0.0;
    double res_rms_pct = 0.0;
};

static std::string make_rf_cc_tag(int run, float voltage, int cc, int rf)
{
    std::string tag = std::string(Form("Run%03d_%.3fV", run, voltage));
    if (cc >= 0)
        tag += Form("_CC%d", cc);
    if (rf >= 0)
        tag += Form("_Rf%d", rf);
    return tag;
}

static bool file_exists_local(const std::string &path)
{
    return !gSystem->AccessPathName(path.c_str());
}

static std::string find_gain_file_for_run(const char *gain_dir, int run, float voltage)
{
    std::vector<std::string> candidates = {
        std::string(Form("%s/gain_match_Run%03d.root", gain_dir, run)),
        std::string(Form("%s/gain_match_run%03d.root", gain_dir, run)),
        std::string(Form("%s/gain_match_Run%03d_%.3fV.root", gain_dir, run, voltage)),
        std::string(Form("%s/gain_match_%.3fV.root", gain_dir, voltage))};

    for (const auto &c : candidates)
        if (file_exists_local(c))
            return c;

    // fall back to the voltage-based path even if it does not exist,
    // so the called function prints a clear message and/or uses unity fallback.
    return std::string(Form("%s/gain_match_%.3fV.root", gain_dir, voltage));
}

static bool read_run_summary_from_root(const std::string &root_path, LedRunSummaryRFCC &row)
{
    TFile *f = TFile::Open(root_path.c_str(), "READ");
    if (!f || f->IsZombie())
    {
        if (f)
        {
            f->Close();
            delete f;
        }
        return false;
    }

    auto getFloat = [&](const char *name, double &out) -> bool
    {
        out = 0.0;
        if (auto *pf = dynamic_cast<TParameter<float> *>(f->Get(name)))
        {
            out = pf->GetVal();
            return true;
        }
        if (auto *pd = dynamic_cast<TParameter<double> *>(f->Get(name)))
        {
            out = pd->GetVal();
            return true;
        }
        return false;
    };

    auto getInt = [&](const char *name, int &out) -> bool
    {
        out = 0;
        if (auto *pi = dynamic_cast<TParameter<int> *>(f->Get(name)))
        {
            out = pi->GetVal();
            return true;
        }
        if (auto *pl = dynamic_cast<TParameter<long> *>(f->Get(name)))
        {
            out = (int)pl->GetVal();
            return true;
        }
        return false;
    };

    bool ok = true;
    ok &= getInt("run_number", row.run);
    ok &= getInt("cc_gain", row.cc);
    ok &= getInt("rf_setting", row.rf);

    double tmp = 0.0;
    ok &= getFloat("label_voltage_V", tmp);
    row.voltage = (float)tmp;

    // Prefer raw mean for consistency with the RMS option; if unavailable use fit mean.
    if (getFloat("summary_mean_raw", row.mean))
        getFloat("summary_mean_raw_err", row.emean);
    else
    {
        ok &= getFloat("summary_mean", row.mean);
        getFloat("summary_sigma_fit_err", row.emean); // not ideal, but keep something if raw mean err absent
    }

    getFloat("summary_sigma_fit", row.sigma);
    getFloat("summary_sigma_fit_err", row.esigma);
    getFloat("summary_rms_raw", row.rms);
    getFloat("summary_rms_raw_err", row.erms);
    getFloat("summary_resolution_pct_fit", row.res_fit_pct);
    getFloat("summary_resolution_pct_fit_err", row.res_fit_pct_err);
    getFloat("summary_resolution_pct_rms", row.res_rms_pct);

    f->Close();
    delete f;
    return ok;
}

static std::vector<float> unique_sorted_voltages_from_rows(const std::vector<LedRunSummaryRFCC> &rows)
{
    std::vector<float> vals;
    for (const auto &r : rows)
        vals.push_back(r.voltage);
    std::sort(vals.begin(), vals.end());
    vals.erase(std::unique(vals.begin(), vals.end(),
                           [](float a, float b)
                           { return std::fabs(a - b) < 1e-6f; }),
               vals.end());
    return vals;
}

static void style_graph_basic(TGraphErrors *g, int marker)
{
    if (!g)
        return;
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(1.1);
    g->SetLineWidth(2);
}

static TLegend *make_basic_legend()
{
    TLegend *leg = new TLegend(0.58, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    return leg;
}

static void fill_graph_from_rows(TGraphErrors *g,
                                 const std::vector<LedRunSummaryRFCC> &rows,
                                 float voltage,
                                 bool scan_rf,
                                 int fixed_value,
                                 bool use_rms,
                                 bool resolution_plot)
{
    if (!g)
        return;

    int ip = 0;
    for (const auto &r : rows)
    {
        if (std::fabs(r.voltage - voltage) > 1e-6f)
            continue;

        if (scan_rf)
        {
            if (r.cc != fixed_value)
                continue;
        }
        else
        {
            if (r.rf != fixed_value)
                continue;
        }

        const double x = scan_rf ? (double)r.rf : (double)r.cc;
        double y = 0.0;
        double ey = 0.0;

        if (resolution_plot)
        {
            if (use_rms)
            {
                y = r.res_rms_pct;
                // simple propagation from RMS/mean if available
                if (r.mean > 0.0 && r.rms > 0.0)
                {
                    double rel2 = 0.0;
                    if (r.erms > 0.0)
                        rel2 += (r.erms / r.rms) * (r.erms / r.rms);
                    if (r.emean > 0.0)
                        rel2 += (r.emean / r.mean) * (r.emean / r.mean);
                    ey = y * std::sqrt(rel2);
                }
            }
            else
            {
                y = r.res_fit_pct;
                ey = r.res_fit_pct_err;
            }
        }
        else
        {
            y = use_rms ? r.rms : r.sigma;
            ey = use_rms ? r.erms : r.esigma;
        }

        const int n = g->GetN();
        g->SetPoint(n, x, y);
        g->SetPointError(n, 0.0, ey);
        ++ip;
    }
}

static void fill_mean_graph_from_rows(TGraphErrors *g,
                                      const std::vector<LedRunSummaryRFCC> &rows,
                                      float voltage,
                                      bool scan_rf,
                                      int fixed_value)
{
    if (!g)
        return;

    for (const auto &r : rows)
    {
        if (std::fabs(r.voltage - voltage) > 1e-6f)
            continue;

        if (scan_rf)
        {
            if (r.cc != fixed_value)
                continue;
        }
        else
        {
            if (r.rf != fixed_value)
                continue;
        }

        const double x = scan_rf ? (double)r.rf : (double)r.cc;
        const int n = g->GetN();
        g->SetPoint(n, x, r.mean);
        g->SetPointError(n, 0.0, r.emean);
    }
}

static void make_adc_rf_cc_summary(const std::vector<LedRunSummaryRFCC> &rows,
                                   const char *outdir,
                                   bool use_rms_for_summary = false)
{
    if (rows.empty())
    {
        std::cerr << "Warning: no rows available for RF/CC summary.\n";
        return;
    }

    const std::vector<float> voltages = unique_sorted_voltages_from_rows(rows);

    for (float voltage_sel : voltages)
    {
        std::vector<LedRunSummaryRFCC> rows_v;
        for (const auto &r : rows)
        {
            if (std::fabs(r.voltage - voltage_sel) < 1e-6f)
                rows_v.push_back(r);
        }

        if (rows_v.empty())
            continue;

        std::map<int, int> cc_counts_for_rf_scan;
        std::map<int, int> rf_counts_for_cc_scan;

        for (const auto &r : rows_v)
        {
            cc_counts_for_rf_scan[r.cc]++;
            rf_counts_for_cc_scan[r.rf]++;
        }

        int ref_cc = rows_v.front().cc;
        int ref_rf = rows_v.front().rf;
        int best = -1;
        for (const auto &kv : cc_counts_for_rf_scan)
            if (kv.second > best)
            {
                best = kv.second;
                ref_cc = kv.first;
            }

        best = -1;
        for (const auto &kv : rf_counts_for_cc_scan)
            if (kv.second > best)
            {
                best = kv.second;
                ref_rf = kv.first;
            }

        std::string pdf = std::string(Form("%s/adc_calibration_rf_cc_summary_%.3fV.pdf", outdir, voltage_sel));
        std::string root_out = std::string(Form("%s/adc_calibration_rf_cc_summary_%.3fV.root", outdir, voltage_sel));

        TCanvas *c = new TCanvas(Form("c_adc_rf_cc_summary_%.3fV", voltage_sel), "", 1500, 900);
        c->SaveAs((pdf + "[").c_str());

        // -------- Page 1: mean and sigma/RMS vs Rf and CC --------
        c->Clear();
        c->Divide(2, 2);

        c->cd(1);
        gPad->SetGrid();
        {
            TGraphErrors *g = new TGraphErrors();
            g->SetName("g_mean_vs_rf");
            g->SetTitle(Form("Central mean vs Rf (V=%.3f, CC=%d);Rf setting;Mean (ADC)", voltage_sel, ref_cc));
            fill_mean_graph_from_rows(g, rows_v, voltage_sel, true, ref_cc);
            style_graph_basic(g, 20);
            if (g->GetN() > 0)
                g->Draw("AP");
            TLatex t;
            t.SetNDC();
            t.SetTextFont(42);
            t.SetTextSize(0.04);
            // t.DrawLatex(0.16, 0.92, Form("Voltage = %.3f V", voltage_sel));
        }

        c->cd(2);
        gPad->SetGrid();
        {
            TGraphErrors *g = new TGraphErrors();
            g->SetName("g_width_vs_rf");
            g->SetTitle(Form("%s vs Rf (V=%.3f, CC=%d);Rf setting;%s (ADC)",
                             use_rms_for_summary ? "Central RMS" : "Central #sigma",
                             voltage_sel, ref_cc,
                             use_rms_for_summary ? "RMS" : "#sigma"));
            fill_graph_from_rows(g, rows_v, voltage_sel, true, ref_cc, use_rms_for_summary, false);
            style_graph_basic(g, 20);
            if (g->GetN() > 0)
                g->Draw("AP");
            TLatex t;
            t.SetNDC();
            t.SetTextFont(42);
            t.SetTextSize(0.04);
            // t.DrawLatex(0.16, 0.92, Form("Voltage = %.3f V", voltage_sel));
        }

        c->cd(3);
        gPad->SetGrid();
        {
            TGraphErrors *g = new TGraphErrors();
            g->SetName("g_mean_vs_cc");
            g->SetTitle(Form("Central mean vs CC (V=%.3f, Rf=%d);CC setting;Mean (ADC)", voltage_sel, ref_rf));
            fill_mean_graph_from_rows(g, rows_v, voltage_sel, false, ref_rf);
            style_graph_basic(g, 20);
            if (g->GetN() > 0)
                g->Draw("AP");
            TLatex t;
            t.SetNDC();
            t.SetTextFont(42);
            t.SetTextSize(0.04);
            // t.DrawLatex(0.16, 0.92, Form("Voltage = %.3f V", voltage_sel));
        }

        c->cd(4);
        gPad->SetGrid();
        {
            TGraphErrors *g = new TGraphErrors();
            g->SetName("g_width_vs_cc");
            g->SetTitle(Form("%s vs CC (V=%.3f, Rf=%d);CC setting;%s (ADC)",
                             use_rms_for_summary ? "Central RMS" : "Central #sigma",
                             voltage_sel, ref_rf,
                             use_rms_for_summary ? "RMS" : "#sigma"));
            fill_graph_from_rows(g, rows_v, voltage_sel, false, ref_rf, use_rms_for_summary, false);
            style_graph_basic(g, 20);
            if (g->GetN() > 0)
                g->Draw("AP");
            TLatex t;
            t.SetNDC();
            t.SetTextFont(42);
            t.SetTextSize(0.04);
            // t.DrawLatex(0.16, 0.92, Form("Voltage = %.3f V", voltage_sel));
        }

        c->SaveAs(pdf.c_str());

        // -------- Page 2: resolution vs Rf and CC --------
        c->Clear();
        c->Divide(2, 2);

        c->cd(1);
        gPad->SetGrid();
        {
            TGraphErrors *g = new TGraphErrors();
            g->SetName("g_res_vs_rf");
            g->SetTitle(Form("Central resolution vs Rf (V=%.3f, CC=%d);Rf setting;Resolution (%%)", voltage_sel, ref_cc));
            fill_graph_from_rows(g, rows_v, voltage_sel, true, ref_cc, use_rms_for_summary, true);
            style_graph_basic(g, 20);
            if (g->GetN() > 0)
                g->Draw("AP");
        }

        c->cd(2);
        gPad->SetGrid();
        {
            TGraphErrors *g = new TGraphErrors();
            g->SetName("g_res_vs_cc");
            g->SetTitle(Form("Central resolution vs CC (V=%.3f, Rf=%d);CC setting;Resolution (%%)", voltage_sel, ref_rf));
            fill_graph_from_rows(g, rows_v, voltage_sel, false, ref_rf, use_rms_for_summary, true);
            style_graph_basic(g, 20);
            if (g->GetN() > 0)
                g->Draw("AP");
        }

        c->cd(3);
        {
            TLatex t;
            t.SetNDC();
            t.SetTextFont(42);
            t.SetTextSize(0.045);
            t.DrawLatex(0.10, 0.86, "Summary configuration");
            t.SetTextSize(0.038);
            t.DrawLatex(0.10, 0.74, Form("Voltage: %.3f V", voltage_sel));
            t.DrawLatex(0.10, 0.66, Form("Width used: %s", use_rms_for_summary ? "RMS" : "Gaussian fit #sigma"));
            t.DrawLatex(0.10, 0.58, Form("Rf scan reference CC: %d", ref_cc));
            t.DrawLatex(0.10, 0.50, Form("CC scan reference Rf: %d", ref_rf));
            t.DrawLatex(0.10, 0.42, Form("Number of processed runs: %d", (int)rows_v.size()));
        }

        c->cd(4);
        {
            TLatex t;
            t.SetNDC();
            t.SetTextFont(42);
            t.SetTextSize(0.040);
            t.DrawLatex(0.10, 0.86, "Interpretation");
            t.SetTextSize(0.035);
            t.DrawLatex(0.10, 0.74, "Mean tracks the signal scale.");
            t.DrawLatex(0.10, 0.66, Form("%s tracks the spread.", use_rms_for_summary ? "RMS" : "#sigma"));
            t.DrawLatex(0.10, 0.58, "Resolution = width / mean.");
            t.DrawLatex(0.10, 0.50, "This PDF contains only one voltage.");
        }

        c->SaveAs(pdf.c_str());
        c->SaveAs((pdf + "]").c_str());
        std::cout << "Saved RF/CC summary PDF: " << pdf << "\n";

        TFile *fout = TFile::Open(root_out.c_str(), "RECREATE");
        if (fout && !fout->IsZombie())
        {
            TTree *t = new TTree("rf_cc_summary", "rf_cc_summary");
            int run = 0, cc = 0, rf = 0;
            float voltage = 0.0f;
            float mean = 0.0f, emean = 0.0f, sigma = 0.0f, esigma = 0.0f, rms = 0.0f, erms = 0.0f;
            float res_fit = 0.0f, res_fit_err = 0.0f, res_rms = 0.0f;

            t->Branch("run", &run, "run/I");
            t->Branch("voltage", &voltage, "voltage/F");
            t->Branch("cc", &cc, "cc/I");
            t->Branch("rf", &rf, "rf/I");
            t->Branch("mean", &mean, "mean/F");
            t->Branch("emean", &emean, "emean/F");
            t->Branch("sigma", &sigma, "sigma/F");
            t->Branch("esigma", &esigma, "esigma/F");
            t->Branch("rms", &rms, "rms/F");
            t->Branch("erms", &erms, "erms/F");
            t->Branch("res_fit_pct", &res_fit, "res_fit_pct/F");
            t->Branch("res_fit_pct_err", &res_fit_err, "res_fit_pct_err/F");
            t->Branch("res_rms_pct", &res_rms, "res_rms_pct/F");

            for (const auto &r : rows_v)
            {
                run = r.run;
                voltage = r.voltage;
                cc = r.cc;
                rf = r.rf;
                mean = (float)r.mean;
                emean = (float)r.emean;
                sigma = (float)r.sigma;
                esigma = (float)r.esigma;
                rms = (float)r.rms;
                erms = (float)r.erms;
                res_fit = (float)r.res_fit_pct;
                res_fit_err = (float)r.res_fit_pct_err;
                res_rms = (float)r.res_rms_pct;
                t->Fill();
            }

            t->Write();
            fout->Close();
            delete fout;
            std::cout << "Saved RF/CC summary ROOT: " << root_out << "\n";
        }

        delete c;
    }
}

void adc_calib_rf_cc_scan(const char *data_dir = "data",
                          const char *mapping_csv = "eeemcal_desy_dec2025_mapping_v2.csv",
                          const char *outdir = "outputs_rf_cc",
                          const char *gain_dir = "outputs",
                          bool use_hybrid_tot = true,
                          bool use_common_mode = true,
                          bool use_rms_for_summary = false)
{
    TH1::AddDirectory(kFALSE);
    gROOT->cd();
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;
    mkdir_p(outdir);

    // IMPORTANT:
    // Use a struct/tuple, not std::vector<std::vector<int,float,...>>.
    std::vector<LedRunCalibRFCC> runs = {
        // {381, 1.26f, 12, 3},
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
        {447, 0.0f, 15, 3}};

    std::vector<LedRunSummaryRFCC> rows;
    rows.reserve(runs.size());

    for (const auto &rv : runs)
    {
        const std::string infile = std::string(Form("%s/Run%03d.root", data_dir, rv.run));
        const std::string gain_file = find_gain_file_for_run(gain_dir, rv.run, rv.voltage);

        std::cout << "\n--- adc calibration RF/CC scan: " << infile
                  << " @ V=" << rv.voltage
                  << "  CC=" << rv.cc
                  << "  Rf=" << rv.rf
                  << "  gain=" << gain_file
                  << " ---\n";

        adc_calibration_rf_cc_one(infile.c_str(),
                                  rv.voltage,
                                  mapping_csv,
                                  gain_file.c_str(),
                                  outdir,
                                  1000000,
                                  0.0f,
                                  use_hybrid_tot,
                                  use_common_mode,
                                  470.0f,
                                  rv.cc,
                                  rv.rf);

        LedRunSummaryRFCC row;
        const std::string root_path = std::string(Form("%s/adc_to_ref_calibration_%s.root",
                                                       outdir,
                                                       make_rf_cc_tag(rv.run, rv.voltage, rv.cc, rv.rf).c_str()));

        if (read_run_summary_from_root(root_path, row))
        {
            rows.push_back(row);
        }
        else
        {
            std::cerr << "Warning: could not read summary back from " << root_path << "\n";
        }
    }

    make_adc_rf_cc_summary(rows, outdir, use_rms_for_summary);
}
