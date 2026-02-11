// draw_waveform.cxx
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <cstdint>
#include <fstream>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TSystem.h> // gSystem->mkdir
#include <TPad.h>

#include "common_led.h"

// LED-friendly version:
// - Uses your mapping CSV (pass mapping_csv)
// - Reads Run%03d.root from a chosen directory (in_dir)
// - Writes a multipage PDF in "outputs/" (out_dir)
// - Still loops over 25 crystals, but safely skips crystals missing from the mapping
// - Safe for files that only have data for 1 crystal
//
// Example:
//   root -l -q 'draw_waveform.cxx+(23, 1.8)'
//   root -l -q 'draw_waveform.cxx+(23, 1.8,".","eeemcal_desy_dec2025_mapping_v2.csv","outputs",16)'
void draw_waveform(int run_number,
                   float voltage,
                   const char *in_dir = "data",
                   const char *mapping_csv = "eeemcal_desy_dec2025_mapping_v2.csv",
                   const char *out_dir = "outputs",
                   int sipms_to_use = g_sipms_to_use)
{
    const int samples_per_channel = 20;

    auto mapping = read_mapping_csv(mapping_csv, sipms_to_use);
    if (mapping.empty())
    {
        std::cerr << "Error: mapping is empty. Check mapping CSV: " << mapping_csv << "\n";
        return;
    }

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

    // --- Create histograms
    std::map<std::string, TH2D *> histograms;

    for (const auto &kv : mapping)
    {
        const int crystal = kv.first;
        const auto &channels = kv.second;
        if (channels.size() < static_cast<size_t>(sipms_to_use))
            continue;

        for (int i = 0; i < sipms_to_use; ++i)
        {
            const std::string hname = Form("crystal_%d_sipm_%d_%.2fV_waveform", crystal, i, voltage);
            const std::string htitle = Form("Run %03d Vol %.2fV | Crystal %d SiPM %d;Time sample;ADC",
                                            run_number, voltage, crystal, i);

            histograms[hname] = new TH2D(hname.c_str(), htitle.c_str(),
                                         samples_per_channel, 0, samples_per_channel,
                                         1024, 0, 1024);
        }
    }

    if (histograms.empty())
    {
        std::cerr << "No histograms created (mapping did not contain enough SiPMs per crystal).\n";
        root_file->Close();
        delete root_file;
        return;
    }

    // --- Fill
    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry)
    {
        tree->GetEntry(entry);

        for (const auto &kv : mapping)
        {
            const int crystal = kv.first;
            const auto &channels = kv.second;
            if (channels.size() < static_cast<size_t>(sipms_to_use))
                continue;

            for (int i = 0; i < sipms_to_use; ++i)
            {
                const int ch = channels[i];
                if (ch < 0 || ch >= nch_from_file)
                    continue;

                const std::string key = Form("crystal_%d_sipm_%d_%.2fV_waveform", crystal, i, voltage);
                auto itH = histograms.find(key);
                if (itH == histograms.end() || !itH->second)
                    continue;

                for (int t = 0; t < samples_per_channel; ++t)
                {
                    itH->second->Fill(t, adc_at(ch, t));
                }
            }
        }
    }

    // --- Draw to PDF (ONE LOOP ONLY)
    TCanvas *canvas = new TCanvas("waveforms", "Waveforms", 1200, 800);
    std::string output_file = std::string(Form("%s/run%03d_%.2fV_waveform.pdf", out_dir, run_number, voltage));

    bool first_page = true;
    int pages = 0;

    canvas->SaveAs((output_file + "[").c_str());

    for (const auto &kv : mapping)
    {
        const int crystal = kv.first;
        const auto &channels = kv.second;
        if (channels.size() < static_cast<size_t>(sipms_to_use))
            continue;

        bool any_valid = false;
        for (int i = 0; i < sipms_to_use; ++i)
        {
            const int ch = channels[i];
            if (ch >= 0 && ch < nch_from_file)
            {
                any_valid = true;
                break;
            }
        }
        if (!any_valid)
            continue;

        canvas->Clear("D");
        canvas->Divide(4, 4);

        for (int i = 0; i < sipms_to_use; ++i)
        {
            canvas->cd(i + 1);
            gPad->Clear();
            gPad->SetLogz();

            const std::string hname = Form("crystal_%d_sipm_%d_%.2fV_waveform", crystal, i, voltage);
            auto itH = histograms.find(hname);
            if (itH != histograms.end() && itH->second)
            {
                itH->second->Draw("COLZ");
            }
        }

        canvas->SaveAs(output_file.c_str());

        pages++;
    }

    canvas->SaveAs((output_file + "]").c_str());
    std::cout << "Saved " << pages << " page(s) to " << output_file << "\n";

    // --- Cleanup
    for (auto &kvh : histograms)
        delete kvh.second;
    delete canvas;

    root_file->Close();
    delete root_file;
}

void led_scan()
{
    // (run, voltage)
    std::vector<std::pair<int, float>> runs = {
        {23, 0.0f},
        {26, 1.2f},
        {30, 1.22f},
        {33, 1.24f},
        {36, 1.25},
        {39, 1.26f},
        {42, 1.27f},
        {45, 1.28f},
        {48, 1.29f},
        {51, 1.3f},
        {54, 1.32f},
        {57, 1.33f},
        {60, 1.34f},
        {63, 1.36f},
        {66, 1.37f},
        {69, 1.38f},
        {72, 1.4f},
        {75, 1.42f},
        {78, 1.44f},
        {81, 1.46f},
        {84, 1.48f},
        {87, 1.5f},
        {90, 1.52f},
        {93, 1.54f},
        {96, 1.56f},
        {99, 1.58f},
        {102, 1.6f},
        {105, 1.62f},
        {108, 1.64f},
        {111, 1.66f},
        {114, 1.68f},
        {117, 1.7f},
        {120, 1.72f},
        {123, 1.74f},
        {126, 1.76f},
        {129, 1.78f},
        {132, 1.8f},
        {135, 1.82f},
        {138, 1.84f},
        {141, 1.86f},
        {144, 1.88f}};

    for (auto &rv : runs)
    {
        int run = rv.first;
        float voltage = rv.second;
        draw_waveform(run, voltage, "data", "eeemcal_desy_dec2025_mapping_v2.csv", "outputs", g_sipms_to_use);
    }
}
