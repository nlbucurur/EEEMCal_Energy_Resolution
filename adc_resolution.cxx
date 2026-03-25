// adc_resolution.cxx
//
// Purpose: Calculate energy resolution for each run and voltage, and plot resolution vs voltage.
// The code reuses parts from adc_calibration.cxx to perform energy analysis and produce the desired plots.
//
// root -l -b
// .L common_led.cxx+
// .L adc_resolution.cxx+
// adc_resolution_scan("data", "eeemcal_desy_dec2025_mapping_v2.csv", "outputs");

#include "common_led.h"
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <istream>
#include <iosfwd>
#include <sstream>

#include <TROOT.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TKey.h>
#include <TError.h>
#include <TSystem.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>


static constexpr int MAX_NUM_CRYSTALS = 25; // Number of crystals
static constexpr int SIPMS_PER_CRYSTAL = 16;
static constexpr int SAMPLES_PER_CHANNEL = 20;
static constexpr int central_crystal_id = 12; // Assuming crystal 12 is the central one (0-indexed)

static bool is_tot_event(uint32_t *tot_values)
{
    for (int i = 0; i < SAMPLES_PER_CHANNEL; ++i)
    {
        if (tot_values[i] > (uint32_t)g_tot_min)
            return true;
    }
    return false;
}

// Fit a single-peak distribution robustly.
// Returns: true on success, and fills mu/sigma and their errors.
static bool fit_peak_robust(TH1 *h, double &mu, double &sigma, double &emu, double &esigma)
{
    mu = sigma = emu = esigma = 0.0;
    if (!h || h->GetEntries() < 50)
        return false;

    const double xMin = h->GetXaxis()->GetXmin();
    const double xMax = h->GetXaxis()->GetXmax();

    const double mean0 = h->GetBinCenter(h->GetMaximumBin());
    const double rms0 = std::max(1e-6, (double)h->GetRMS());

    double fmin = std::max(xMin, mean0 - 1.5 * rms0);
    double fmax = std::min(xMax, mean0 + 1.5 * rms0);
    if (fmax - fmin < 0.07 * (xMax - xMin))
    {
        fmin = std::max(xMin, mean0 - 0.25 * (xMax - xMin));
        fmax = std::min(xMax, mean0 + 0.25 * (xMax - xMin));
    }

    TF1 rough("rough_fit", "gaus", fmin, fmax);
    rough.SetParameters(h->GetMean(), mean0, rms0);
    TFitResultPtr r0 = h->Fit(&rough, "RQS");
    if (!r0.Get() || r0->Status() != 0)
        return false;

    double m = rough.GetParameter(1);
    double s = std::abs(rough.GetParameter(2));
    if (!(s > 0))
        return false;

    // Refined gaussian in m±s
    double gmin = std::max(xMin, m - s);
    double gmax = std::min(xMax, m + s);
    if (gmax - gmin < 0.03 * (xMax - xMin))
    {
        mu = m;
        sigma = s;
        emu = rough.GetParError(1);
        esigma = rough.GetParError(2);
        return true;
    }

    TF1 g2("second_fit", "gaus", gmin, gmax);
    g2.SetParameters(rough.GetParameter(0), m, s);
    TFitResultPtr r1 = h->Fit(&g2, "RQS");
    if (r1.Get() && r1->Status() == 0)
    {
        m = g2.GetParameter(1);
        s = std::abs(g2.GetParameter(2));
    }

    // Final: crystalball in m±s (usually helps tails)
    double cmin = std::max(xMin, m - s);
    double cmax = std::min(xMax, m + s);
    if (cmax - cmin < 0.03 * (xMax - xMin))
    {
        mu = m;
        sigma = s;
        emu = g2.GetParError(1);
        esigma = g2.GetParError(2);
        return true;
    }

    TF1 cb("final_fit", "crystalball", cmin, cmax);
    cb.SetParameters(g2.GetParameter(0), m, s, 1.5, 2.0);
    TFitResultPtr r2 = h->Fit(&cb, "RQS");

    if (r2.Get() && r2->Status() == 0)
    {
        mu = cb.GetParameter(1);
        sigma = std::abs(cb.GetParameter(2));
        emu = cb.GetParError(1);
        esigma = cb.GetParError(2);
        return true;
    }

    // fallback to gaussian
    mu = m;
    sigma = s;
    emu = g2.GetParError(1);
    esigma = g2.GetParError(2);
    return true;
}



// Function to calculate and plot energy resolution
// Function to calculate and plot energy resolution
void adc_resolution_one(const char* filename,
                        float voltage,
                        const char* mapping_csv = "eeemcal_desy_dec2025_mapping_v2.csv",
                        const char* gain_root_file = "outputs/gain_match_1.25V.root",
                        const char* outdir = "outputs",
                        Long64_t max_events = 1000000,
                        bool use_hybrid_tot = true)
{
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;
    gSystem->mkdir(outdir, true);

    // Open the input ROOT file
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open " << filename << "\n";
        return;
    }

    TTree *tree = (TTree*)f->Get("events");
    if (!tree) {
        std::cerr << "ERROR: cannot find TTree 'events' in " << filename << "\n";
        return;
    }

    // Mapping: read the mapping CSV file (similar to adc_calibration.cxx)
    auto mapping = read_mapping_csv(mapping_csv, SIPMS_PER_CRYSTAL);
    if (mapping.empty()) {
        std::cerr << "ERROR: mapping empty. Check: " << mapping_csv << "\n";
        return;
    }

    // Reverse mapping: channel -> (crystal, sipm)
    std::unordered_map<int, std::pair<int, int>> reverse;
    reverse.reserve(mapping.size() * SIPMS_PER_CRYSTAL);
    for (const auto& kv : mapping) {
        int cr = kv.first;
        auto chans = get_crystal_channels(mapping, cr, SIPMS_PER_CRYSTAL);
        for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; ++sipm) {
            int ch = chans[sipm];
            if (ch >= 0)
                reverse[ch] = {cr, sipm};
        }
    }

    // Load gain factors (similar to adc_calibration.cxx)
    TH1 *gain_factor = nullptr;
    TFile *gf = TFile::Open(gain_root_file);
    if (gf && !gf->IsZombie()) {
        gain_factor = (TH1*)gf->Get("gain_factors");
    }

    if (!gain_factor) {
        gain_factor = new TH1F("gain_factors_fallback", "gain_factors_fallback", MAX_NUM_CRYSTALS * SIPMS_PER_CRYSTAL, 0, MAX_NUM_CRYSTALS * SIPMS_PER_CRYSTAL);
        for (int i = 1; i <= MAX_NUM_CRYSTALS * SIPMS_PER_CRYSTAL; ++i) {
            gain_factor->SetBinContent(i, 1.0);
        }
    }

    // Branch leaves for dynamic sizing
    TBranch *br_adc = tree->GetBranch("adc");
    TBranch *br_tot = tree->GetBranch("tot");

    if (!br_adc) {
        std::cerr << "ERROR: missing branch adc\n";
        return;
    }

    TLeaf *leaf_adc = br_adc->GetLeaf("adc");
    if (!leaf_adc) {
        std::cerr << "ERROR: missing leaf adc\n";
        return;
    }

    TLeaf *leaf_tot = (br_tot ? br_tot->GetLeaf("tot") : nullptr);

    const int n_adc = leaf_adc->GetLen();
    const int n_channels = n_adc / SAMPLES_PER_CHANNEL;

    bool have_tot = (leaf_tot && leaf_tot->GetLen() == n_adc);

    std::vector<uint32_t> adc_buf(n_adc);
    std::vector<uint32_t> tot_buf;

    if (have_tot) {
        tot_buf.resize(n_adc);
    }

    tree->SetBranchAddress("adc", adc_buf.data());
    if (have_tot) {
        tree->SetBranchAddress("tot", tot_buf.data());
    }

    // Create adc_ptr and tot_ptr for accessing ADC and ToT values, respectively
    auto adc_ptr = [&](int ch) -> uint32_t* { return &adc_buf[ch * SAMPLES_PER_CHANNEL]; };
    auto tot_ptr = [&](int ch) -> uint32_t* { return &tot_buf[ch * SAMPLES_PER_CHANNEL]; };

    // Active channels based on mapping (bounded by file channels)
    auto active_channels = get_active_channels_from_mapping(mapping, SIPMS_PER_CRYSTAL, MAX_NUM_CRYSTALS, n_channels);

    // Histograms for central crystal energy
    TH1F *central_crystal_energy = new TH1F("central_crystal_energy",
                                            Form("Central Crystal Energy (%.3f V);Energy (ADC);Events", voltage),
                                            500, 0, 75000);

    // Loop over events and accumulate energy for central crystal
    Long64_t nentries = tree->GetEntries();
    if (max_events > 0 && max_events < nentries)
        nentries = max_events;

    for (Long64_t entry = 0; entry < nentries; ++entry) {
        tree->GetEntry(entry);

        // Skip ToT events: Check if it's a ToT event and skip it
        bool is_tot = false;
        if (use_hybrid_tot) {
            is_tot = is_tot_event(tot_ptr);  // Check if it's a ToT event
        }

        if (is_tot) {
            continue;  // Skip ToT events
        }

        // Process signals for the event
        float channel_signal = 0.0f;
        for (int ch : active_channels) {
            if (ch < 0 || ch >= n_channels)
                continue;

            auto it = reverse.find(ch);
            if (it == reverse.end())
                continue;

            int cr = it->second.first;
            int sipm = it->second.second;

            // Gain from the gain factor histogram
            float gain = gain_factor->GetBinContent(cr * SIPMS_PER_CRYSTAL + sipm + 1);

            // Get the ADC value and apply the gain
            uint32_t* adc_values = adc_ptr(ch);
            channel_signal = calculate_signal_adc(adc_values, gain);

            // Accumulate the energy for the central crystal
            if (cr == central_crystal_id) {
                central_crystal_energy->Fill(channel_signal);
            }
        }
    }

    // Fit the central crystal energy histogram to extract resolution
    double mu, sigma, emu, esigma;
    bool fit_success = fit_peak_robust(central_crystal_energy, mu, sigma, emu, esigma);
    double resolution = (fit_success && mu != 0) ? (sigma / mu) : -1.0;

    if (fit_success)
    {
        std::cout << "Energy resolution for run " << filename << " @ " << voltage << " V: " << 100.0 * resolution << "%" << std::endl;
    }
    else
    {
        std::cout << "Fit failed for run " << filename << " @ " << voltage << " V." << std::endl;
    }

    // Save the resolution vs voltage in a new PDF page
    TCanvas *canvas = new TCanvas("adc_resolution_canvas", "", 900, 700);
    canvas->SetRightMargin(0.05);

    // Plot the resolution
    central_crystal_energy->Draw("hist e");
    canvas->SaveAs(Form("%s/adc_resolution_Run%03d_%.3fV.pdf", outdir, 123, voltage)); // Save the resolution plot

    // Save the results in a ROOT file
    TFile *out = TFile::Open(Form("%s/adc_resolution_%.3fV.root", outdir, voltage), "RECREATE");
    central_crystal_energy->Write();
    out->Close();

    delete canvas;
    f->Close();
    delete f;
}

// Function to run resolution calculations for multiple runs/voltages
void adc_resolution_scan(const char* data_dir = "data",
                         const char* mapping_csv = "eeemcal_desy_dec2025_mapping_v2.csv",
                         bool use_hybrid_tot = true)
{
    // List of runs and corresponding voltages
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
        {170, 1.25f},
        {171, 1.259f},
        {172, 1.268f}//,
        // {173, 1.277f},
        // {174, 1.286f},
        // {175, 1.295f},
        // {176, 1.304f},
        // {177, 1.313f},
        // {178, 1.322f},
        // {179, 1.331f},
        // {180, 1.34f}//,
        // {181, 1.349f},
        // {182, 1.358f},
        // {183, 1.367f}
    };

    for (auto& rv : runs) {
        const int run = rv.first;
        const float voltage = rv.second;

        std::string fn = std::string(Form("%s/Run%03d.root", data_dir, run));
        std::cout << "\n--- Calculating resolution for Run " << run << " @ " << voltage << " V ---\n";
        adc_resolution_one(fn.c_str(), voltage, mapping_csv, Form("outputs/gain_match_%.3fV.root", voltage), "outputs", 1000000, use_hybrid_tot);
    }
}