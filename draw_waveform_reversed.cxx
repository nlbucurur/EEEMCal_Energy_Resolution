// draw_waveform.cxx (Option B style: channel-centric + reverse map)
//
// Example:
//   root -l -q 'common_led.cxx+ draw_waveform_reversed.cxx+(23, 1.80)'
//   root -l -q 'common_led.cxx+ draw_waveform_reversed.cxx+(23, 1.80,"data","eeemcal_desy_dec2025_mapping_v2.csv","outputs",16)'
//   root -l -q 'common_led.cxx+ draw_waveform_reversed.cxx+ led_scan()'
// Or interactive:
//   root -l -b
//   .L common_led.cxx+
//   .L draw_waveform_reversed.cxx+
//   led_scan()

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <array>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TSystem.h> // gSystem->mkdir
#include <TPad.h>
#include <TROOT.h>

#include "common_led.h"

void draw_waveform_reversed(int run_number,
                   float voltage,
                   const char *in_dir = "data",
                   const char *mapping_csv = "eeemcal_desy_dec2025_mapping_v2.csv",
                   const char *out_dir = "outputs",
                   int sipms_to_use = g_sipms_to_use)
{
    const int samples_per_channel = LED_SAMPLES_PER_CHANNEL;

    // --- Read mapping
    auto mapping = read_mapping_csv(mapping_csv, sipms_to_use);
    if (mapping.empty())
    {
        std::cerr << "Error: mapping is empty. Check mapping CSV: " << mapping_csv << "\n";
        return;
    }

    // --- Open ROOT file
    TFile *root_file = TFile::Open(Form("%s/Run%03d.root", in_dir, run_number));
    if (!root_file || root_file->IsZombie())
    {
        std::cerr << "Error: could not open input file: "
                  << Form("%s/Run%03d.root", in_dir, run_number) << "\n";
        if (root_file)
        {
            root_file->Close();
            delete root_file;
        }
        return;
    }

    TTree *tree = (TTree *)root_file->Get("events");
    if (!tree)
    {
        std::cerr << "Error: could not find TTree 'events' in file: "
                  << Form("%s/Run%03d.root", in_dir, run_number) << "\n";
        root_file->Close();
        delete root_file;
        return;
    }

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("adc", 1);

    TBranch *br_adc = tree->GetBranch("adc");
    TLeaf *leaf_adc = br_adc ? br_adc->GetLeaf("adc") : nullptr;
    if (!leaf_adc)
    {
        std::cerr << "Error: could not access leaf 'adc'\n";
        root_file->Close();
        delete root_file;
        return;
    }

    const int n_adc = leaf_adc->GetLen();
    if (n_adc <= 0 || (n_adc % samples_per_channel) != 0)
    {
        std::cerr << "Error: unexpected adc length: " << n_adc
                  << " (must be multiple of " << samples_per_channel << ")\n";
        root_file->Close();
        delete root_file;
        return;
    }

    const int nch_from_file = n_adc / samples_per_channel;

    std::vector<uint32_t> adc_buf((size_t)n_adc, 0);
    tree->SetBranchAddress("adc", adc_buf.data());

    auto adc_at = [&](int ch, int t) -> uint32_t
    {
        return adc_buf[(size_t)ch * samples_per_channel + (size_t)t];
    };

    gSystem->mkdir(out_dir, true);

    // ============================================================
    // Option B: build reverse map and active channel list once
    // reverse[ch] = {crystal_id, sipm}
    // ============================================================

    auto reverse = build_reverse_channel_map(mapping, sipms_to_use);

    auto active_channels = get_active_channels_from_mapping(mapping, sipms_to_use, LED_MAX_NUM_CRYSTALS, nch_from_file);
    if (active_channels.empty())
    {
        std::cerr << "Error: active_channels is empty (check mapping / sipms_to_use).\n";
        root_file->Close();
        delete root_file;
        return;
    }

    // ============================================================
    // Histograms per crystal, per sipm (for drawing pages)
    // We'll store as: h2[crystal_id][sipm] = TH2D*
    // ============================================================

    std::map<int, std::vector<TH2D *>> h2; // vector size sipms_to_use

    for (const auto &kv : mapping)
    {
        const int crystal_id = kv.first;
        h2[crystal_id].resize((size_t)sipms_to_use, nullptr);

        for (int sipm = 0; sipm < sipms_to_use; ++sipm)
        {
            std::string hname = Form("crystal_%d_sipm_%d_%.2fV_waveform", crystal_id, sipm, voltage);
            std::string htitle = Form("Run %03d Vol %.2fV | Crystal %d SiPM %d;Time sample;ADC",
                                      run_number, voltage, crystal_id, sipm);

            // Keep your original binning/range
            h2[crystal_id][sipm] = new TH2D(hname.c_str(), htitle.c_str(),
                                           samples_per_channel, 0, samples_per_channel,
                                           1024, 0, 1024);
        }
    }

    // ============================================================
    // Fill ONCE: loop events, then loop active channels
    // ============================================================
    const Long64_t nEntries = tree->GetEntries();

    for (Long64_t entry = 0; entry < nEntries; ++entry)
    {
        tree->GetEntry(entry);

        for (int ch : active_channels)
        {
            if (ch < 0 || ch >= nch_from_file) continue;

            auto it = reverse.find(ch);
            if (it == reverse.end()) continue;

            const int crystal_id = it->second.first;
            const int sipm       = it->second.second;

            // Safety: crystal_id may exist but vector may not (should not happen)
            auto itC = h2.find(crystal_id);
            if (itC == h2.end()) continue;
            if (sipm < 0 || sipm >= (int)itC->second.size()) continue;

            TH2D *hist = itC->second[(size_t)sipm];
            if (!hist) continue;

            for (int t = 0; t < samples_per_channel; ++t)
                hist->Fill(t, adc_at(ch, t));
        }
    }

    // ============================================================
    // Draw PDF: one page per crystal (4x4)
    // Only draw crystals that have at least 1 valid channel in the file
    // ============================================================
    TCanvas *canvas = new TCanvas("waveforms", "Waveforms", 1200, 800);
    std::string output_file = std::string(Form("%s/run%03d_%.2fV_waveform.pdf", out_dir, run_number, voltage));

    canvas->SaveAs((output_file + "[").c_str());

    int pages = 0;

    for (const auto &kv : mapping)
    {
        const int crystal_id = kv.first;

        // Check if this crystal has any channel that exists in this ROOT file
        bool any_valid = false;
        auto chans = get_crystal_channels(mapping, crystal_id, sipms_to_use);
        for (int sipm = 0; sipm < sipms_to_use; ++sipm)
        {
            const int ch = chans[sipm];
            if (ch >= 0 && ch < nch_from_file)
            {
                any_valid = true;
                break;
            }
        }
        if (!any_valid) continue;

        canvas->Clear("D");
        canvas->Divide(4, 4);

        for (int sipm = 0; sipm < sipms_to_use; ++sipm)
        {
            canvas->cd(sipm + 1);
            gPad->Clear();
            gPad->SetLogz();

            TH2D *hist = nullptr;
            auto itC = h2.find(crystal_id);
            if (itC != h2.end() && sipm >= 0 && sipm < (int)itC->second.size())
                hist = itC->second[(size_t)sipm];

            if (hist) hist->Draw("COLZ");
        }

        canvas->SaveAs(output_file.c_str());
        pages++;
    }

    canvas->SaveAs((output_file + "]").c_str());
    std::cout << "Saved " << pages << " page(s) to " << output_file << "\n";

    // --- Cleanup
    for (auto &kv : h2)
        for (TH2D *h : kv.second)
            delete h;

    delete canvas;

    root_file->Close();
    delete root_file;
}

void led_scan()
{
    std::vector<std::pair<int, float>> runs = {
        {23, 0.0f}, {26, 1.2f}, {30, 1.22f}, {33, 1.24f}, {36, 1.25f},
        {39, 1.26f}, {42, 1.27f}, {45, 1.28f}, {48, 1.29f}, {51, 1.30f},
        {54, 1.32f}, {57, 1.33f}, {60, 1.34f}, {63, 1.36f}, {66, 1.37f},
        {69, 1.38f}, {72, 1.40f}, {75, 1.42f}, {78, 1.44f}, {81, 1.46f},
        {84, 1.48f}, {87, 1.50f}, {90, 1.52f}, {93, 1.54f}, {96, 1.56f},
        {99, 1.58f}, {102, 1.60f}, {105, 1.62f}, {108, 1.64f}, {111, 1.66f},
        {114, 1.68f}, {117, 1.70f}, {120, 1.72f}, {123, 1.74f}, {126, 1.76f},
        {129, 1.78f}, {132, 1.80f}, {135, 1.82f}, {138, 1.84f}, {141, 1.86f},
        {144, 1.88f}
    };

    for (auto &rv : runs)
        draw_waveform_reversed(rv.first, rv.second, "data", "eeemcal_desy_dec2025_mapping_v2.csv", "outputs", g_sipms_to_use);
}
