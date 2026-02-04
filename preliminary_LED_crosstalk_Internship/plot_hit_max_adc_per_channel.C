/** The way to run this code is:
root -l
.L plot_hit_max_adc_per_channel.C
plot_hit_max_adc_per_channel(32); // Just give the run number! **/

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <string>

void plot_hit_max_adc_per_channel(int run_number = 32)
{   
    // Tristan Channels
    std::vector<int> ignored_channels = {8, 17, 26, 35, 44, 53, 62, 71, 80, 89, 98, 107, 116, 125, 134, 143};

    const int num_kcu = 3;
    const int channels_per_kcu = 144;
    const int num_channels = num_kcu * channels_per_kcu;

    // Label maps
    std::vector<int> fpga_ids = {208, 209, 210};
    std::vector<std::string> protoboards = {"pb01", "pb008", "pb06"};

    const TString input_dir = "/home/eic/test_LED/crosstalk_june2025/decoded_data/";
    const TString output_dir = "/home/eic/test_LED/crosstalk_june2025/hit_max_plots/";

    TString run_id = Form("run%03d", run_number);
    TString filename = input_dir + run_id + ".root";

    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie())
    {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    TTree *t = (TTree *)f->Get("events");
    if (!t)
    {
        std::cerr << "Cannot find TTree 'tree' in file." << std::endl;
        return;
    }

    uint hit_max[num_channels];
    uint hit_pedestal[num_channels];
    t->SetBranchAddress("hit_max", hit_max);
    t->SetBranchAddress("hit_pedestal", hit_pedestal);

    std::vector<uint> max_adc(num_channels, 0);
    std::vector<int> max_minus_pedestal(num_channels, 0);

    Long64_t nEntries = t->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        t->GetEntry(i);
        for (int ch = 0; ch < num_channels; ++ch)
        {
            if (hit_max[ch] > max_adc[ch])
            {
                max_adc[ch] = hit_max[ch];
                max_minus_pedestal[ch] = hit_max[ch] - hit_pedestal[ch];
            }
        }
    }

    auto style_histogram = [](TH1 *h, const char *title, const char *ytitle)
    {
        h->SetTitle(title);
        h->GetXaxis()->SetTitle("Local Channel (0 - 143)");
        h->GetYaxis()->SetTitle(ytitle);
        h->GetXaxis()->SetTitleSize(0.045);
        h->GetYaxis()->SetTitleSize(0.045);
        h->GetXaxis()->SetLabelSize(0.04);
        h->GetYaxis()->SetLabelSize(0.04);
        h->GetXaxis()->SetTitleOffset(1.2);
        h->GetYaxis()->SetTitleOffset(1.3);
        h->SetLineWidth(2);
    };

    for (int fpga = 0; fpga < num_kcu; ++fpga)
    {
        TString label = Form("fpga%d_%d_%s", fpga, fpga_ids[fpga], protoboards[fpga].c_str());

        TString out1 = output_dir + run_id + "_" + label + "_max_adc_vs_channel.png";
        TString out2 = output_dir + run_id + "_" + label + "_max_minus_pedestal_vs_channel.png";

        int ch_start = fpga * channels_per_kcu;
        int ch_end = ch_start + channels_per_kcu;

        // Max ADC per channel (per FPGA)
        TCanvas *c1 = new TCanvas(Form("c1_fpga%d", fpga), "Max ADC", 1000, 500);
        c1->SetGrid();
        TH1I *h1 = new TH1I(Form("h1_fpga%d", fpga),
                            Form("Max ADC - FPGA %d;Channel;ADC Max", fpga),
                            channels_per_kcu, 0, channels_per_kcu);
        for (int ch = ch_start; ch < ch_end; ++ch)
        {
            if (std::find(ignored_channels.begin(), ignored_channels.end(), ch % 144) != ignored_channels.end())
                continue;
            h1->SetBinContent(ch - ch_start + 1, max_adc[ch]);
        }
        h1->SetLineColor(kBlue + 1);
        style_histogram(h1, Form("Max ADC - FPGA %d (ID %d, %s)", fpga, fpga_ids[fpga], protoboards[fpga].c_str()), "ADC Max");
        h1->SetMaximum(1024);
        h1->Draw("HIST");
        c1->SaveAs(out1);

        // Max ADC - Pedestal (per FPGA)
        TCanvas *c2 = new TCanvas(Form("c2_fpga%d", fpga), "Max - Pedestal", 1000, 500);
        c2->SetGrid();
        TH1I *h2 = new TH1I(Form("h2_fpga%d", fpga),
                            Form("Max - Pedestal - FPGA %d;Channel;ADC Max - Pedestal", fpga),
                            channels_per_kcu, 0, channels_per_kcu);
        for (int ch = ch_start; ch < ch_end; ++ch)
        {
            h2->SetBinContent(ch - ch_start + 1, max_minus_pedestal[ch]);
        }
        h2->SetLineColor(kRed + 1);
        style_histogram(h2, Form("Max - Pedestal - FPGA %d (ID %d, %s)", fpga, fpga_ids[fpga], protoboards[fpga].c_str()), "ADC Max - Pedestal");
        h2->SetMaximum(1024);
        h2->Draw("HIST");
        c2->SaveAs(out2);
    }
}
