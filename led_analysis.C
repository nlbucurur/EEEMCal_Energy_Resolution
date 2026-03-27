// led_analysis.C
// Run (one line):
//   root -l -q 'common_led.cxx+ led_analysis.C+ led_scan()'
// Or interactive:
//   root -l -b
//   .L common_led.cxx+
//   .L led_analysis.C+
//   led_scan()

#include "common_led.h"

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TH1F.h>

#include <string>
#include <cctype>
#include <iostream>
#include <vector>
#include <utility>
#include <unordered_map>
#include <map>
#include <array>

#include <sys/stat.h>
#include <sys/types.h>

// static constexpr int NUM_SIPMS = 16; // Number of SiPMs to use per crystal (= 16)

void led_analysis_one(const char *filename,
                      float led_voltage,
                      const char *mapping_csv = "eeemcal_desy_dec2025_mapping_v2.csv",
                      const char *outdir = "outputs")
{

    // ===== Configuration =====
    g_signal_method = 3; // 2, 3, 4, 5, 7 (7 = waveform crystal ball fit)
    g_tot_min = 50;      // Minimum ToT value to consider valid

    const int samples_per_channel = LED_SAMPLES_PER_CHANNEL;

    // ===== Run number =====
    int run_number = extract_run_number(filename);
    if (run_number < 0)
    {
        std::cerr << "Warning: could not extract run number from filename, using -1\n";
    }

    // ===== Read mapping =====
    auto mapping = read_mapping_csv(mapping_csv, LED_SIPMS_PER_CRYSTAL);
    if (mapping.empty())
    {
        std::cerr << "Error: mapping empty. Check mapping CSV: " << mapping_csv << "\n";
        return;
    }

    // ===== Open file =====
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Could not open file " << filename << "\n";
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
        std::cerr << "Error: Could not find TTree 'events' in file " << filename << "\n";
        file->Close();
        delete file;
        return;
    }

    // ===== Branches =====
    // Adjust dimensions

    TBranch *br_adc = tree->GetBranch("adc");
    TBranch *br_tot = tree->GetBranch("tot");
    if (!br_adc || !br_tot)
    {
        std::cerr << "Missing adc/tot branches.\n";
        file->Close();
        delete file;
        return;
    }

    TLeaf *leaf_adc = br_adc->GetLeaf("adc");
    TLeaf *leaf_tot = br_tot->GetLeaf("tot");
    if (!leaf_adc || !leaf_tot)
    {
        std::cerr << "Missing adc/tot leaves.\n";
        file->Close();
        delete file;
        return;
    }

    const int n_adc = leaf_adc->GetLen(); // total uint32 elements
    const int n_tot = leaf_tot->GetLen();
    if (n_adc != n_tot || (n_adc % samples_per_channel) != 0)
    {
        std::cerr << "Unexpected adc/tot shape: n_adc=" << n_adc
                  << " n_tot=" << n_tot
                  << " samples_per_channel=" << samples_per_channel << "\n";
        file->Close();
        delete file;
        return;
    }

    const int num_channels = n_adc / samples_per_channel;

    std::vector<uint32_t> adc_buf((size_t)n_adc);
    std::vector<uint32_t> tot_buf((size_t)n_tot);

    tree->SetBranchAddress("adc", adc_buf.data());
    tree->SetBranchAddress("tot", tot_buf.data());

    // ===== Ensure output directory exists =====
    mkdir_p(outdir);

    // reverse[channel] = {crystal_id, sipm_index}
    auto reverse = build_reverse_channel_map(mapping, LED_SIPMS_PER_CRYSTAL);
    // std::unordered_map<int, std::pair<int, int>> reverse;
    // reverse.reserve(mapping.size() * LED_SIPMS_PER_CRYSTAL);

    // for (const auto &kv : mapping)
    // {
    //     int crystal_id = kv.first;
    //     auto chans = get_crystal_channels(mapping, crystal_id, LED_SIPMS_PER_CRYSTAL);

    //     for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
    //     {
    //         int ch = chans[sipm];
    //         if (ch < 0)
    //             continue;

    //         auto it = reverse.find(ch);
    //         if (it != reverse.end())
    //         {
    //             std::cerr << "Duplicate channel in mapping: ch=" << ch
    //                       << " previously (cr=" << it->second.first << ", sipm=" << it->second.second << ")"
    //                       << " now (cr=" << crystal_id << ", sipm=" << sipm << ")\n";
    //         }

    //         reverse[ch] = {crystal_id, sipm};
    //     }
    // }

    // ===== Create histograms for ALL crystals =====
    std::map<int, std::array<TH1F *, LED_SIPMS_PER_CRYSTAL>> h_adc;
    std::map<int, std::array<TH1F *, LED_SIPMS_PER_CRYSTAL>> h_tot;

    for (const auto &kv : mapping)
    {
        const int crystal_id = kv.first;

        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
        {
            h_adc[crystal_id][sipm] = new TH1F(Form("h_adc_run%d_%.3fV_cr%d_sipm%d", run_number, led_voltage, crystal_id, sipm),
                                               Form("ADC | run %d | voltage %.3fV | crystal %d | sipm %d",
                                                    run_number, led_voltage, crystal_id, sipm),
                                               200, 0, 1024);
            h_tot[crystal_id][sipm] = new TH1F(Form("h_tot_run%d_%.3fV_cr%d_sipm%d", run_number, led_voltage, crystal_id, sipm),
                                               Form("ToT | run %d | voltage %.3fV | crystal %d | sipm %d",
                                                    run_number, led_voltage, crystal_id, sipm),
                                               200, 0, 4096);
        }
    }

    // ===== Active channels ======
    // auto active_channels = get_active_channels_from_mapping(mapping, LED_SIPMS_PER_CRYSTAL);
    auto active_channels = get_active_channels_from_mapping(mapping,
                                                            LED_SIPMS_PER_CRYSTAL,
                                                            LED_MAX_NUM_CRYSTALS,
                                                            num_channels);

    // ===== Event loop =====
    const Long64_t nEntries = tree->GetEntries();

    for (Long64_t ev = 0; ev < nEntries; ++ev)
    {
        tree->GetEntry(ev);

        for (int ch : active_channels)
        {
            // safety: file might contain fewer channels than mapping expects
            if (ch < 0 || ch >= num_channels)
                continue;

            // look up which crystal/sipm this channel belongs to
            auto it = reverse.find(ch);
            if (it == reverse.end())
                continue;

            int crystal_id = it->second.first;
            int sipm = it->second.second;

            float sig = 0.0f;
            bool used_tot = false;

            // Pointers to the first sample of this channel
            uint32_t *adc_ptr = &adc_buf[(size_t)ch * samples_per_channel];
            uint32_t *tot_ptr = &tot_buf[(size_t)ch * samples_per_channel];

            calculate_signal_hybrid(adc_ptr, tot_ptr, 1.0f, sig, used_tot);

            if (used_tot)
                h_tot[crystal_id][sipm]->Fill(sig);
            else
                h_adc[crystal_id][sipm]->Fill(sig);
        }
    }

    // ===== Save output =====

    for (const auto &kv : mapping)
    {
        const int crystal_id = kv.first;

        TFile *outFile = new TFile(
            Form("%s/led_output_crystal%d_run%d_%.3fV.root", outdir, crystal_id, run_number, led_voltage),
            "RECREATE");

        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
        {
            h_adc[crystal_id][sipm]->Write();
            h_tot[crystal_id][sipm]->Write();
        }

        outFile->Close();
        delete outFile;

        std::cout << "Saved: "
                  << Form("%s/led_output_crystal%d_run%d_%.3fV.root",
                          outdir, crystal_id, run_number, led_voltage)
                  << "\n";
    }

    // ===== Cleanup histograms =====
    for (auto &kv : h_adc)
        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
            delete kv.second[sipm];

    for (auto &kv : h_tot)
        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
            delete kv.second[sipm];

    // Close input ONCE
    file->Close();
    delete file;
}

void led_scan()
{
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
        {170, 1.25f},
        {171, 1.259f},
        {172, 1.268f},
        {173, 1.277f},
        {174, 1.286f},
        {175, 1.295f},
        {176, 1.304f},
        {177, 1.313f},
        {178, 1.322f},
        {179, 1.331f},
        {180, 1.34f} //,
        // {181, 1.349f},
        // {182, 1.358f},
        // {183, 1.367f}
    };

    for (auto &rv : runs)
    {
        led_analysis_one(Form("data/Run%03d.root", rv.first), rv.second);
    }
}