// To run ethe whole code, put
// root -l -q -b 'plot_hit_max_tot_led_per_channel.C++' -e 'run_macro()'

// To run plot_hit_max_tot_led_per_channel(): root -l -q 'plot_hit_max_tot_led_per_channel.C(run_on,run_off, asic = 0 or 1, port = 0, 1, 2 or 3)'
// example: root -l -q -b 'plot_hit_max_tot_led_per_channel.C(321, 345, 0, 1, "008")'

#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <TLatex.h>
#include <TLine.h>
#include <TPaveText.h>

typedef unsigned int uint;

void plot_hit_max_tot_led_per_channel(int run_on = 321, int run_off = 345, int asic = 0, int port = 0, std::string injected_pb = "008")
{
    // gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.03, "X");
    gStyle->SetLabelSize(0.03, "Y");
    gStyle->SetLabelSize(0.03, "Z");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLabelOffset(0.005, "X");
    gStyle->SetLabelOffset(0.005, "Y");

    std::vector<int> ignored_channels = {8, 17, 26, 35, 44, 53, 62, 71, 80, 89, 98, 107, 116, 125, 134, 143};

    const int num_kcu = 3;
    const int channels_per_kcu = 144;
    const int num_channels = num_kcu * channels_per_kcu;
    const int num_channels_per_port = 16;
    int num_ports = 8;
    std::vector<int> port_bins;
    std::vector<int> filled_bins;

    double avg_adc_injected_channels = 0;
    double avg_adc_cross_talk_expected_channels = 0;
    double cross_talk = 0;

    std::vector<int> shihai_channel_map;

    // End channel of each port (in global channels)
    std::vector<int> port_ends_local = {17, 35, 53, 71, 89, 107, 125, 143};
    double y_max = 1024;
    double y_label = y_max * 1.08;

    std::vector<std::string> region_labels = {
        "ASIC0 - Port A", "ASIC0 - Port B", "ASIC0 - Port D", "ASIC0 - Port C",
        "ASIC1 - Port A", "ASIC1 - Port B", "ASIC1 - Port D", "ASIC1 - Port C"};

    for (int ch = 1; ch <= 151; ++ch)
    {
        if (ch == 19 || ch == 38 || ch == 57 || ch == 76 || ch == 95 || ch == 114 || ch == 133)
            continue;

        shihai_channel_map.push_back(ch);
    }

    std::vector<int> fpga_ids = {208, 209, 210};
    std::vector<std::string> protoboards = {"pb01", "pb008", "pb06"};

    const TString input_dir = "/home/eic/test_LED/crosstalk_june2025/decoded_data/";
    const TString output_dir = "/home/eic/test_LED/crosstalk_june2025/hit_max_plots/Led-off-vs-led-on/";

    TString run_on_id = Form("run%03d", run_on);
    TString run_off_id = Form("run%03d", run_off);

    TFile *f_on = TFile::Open(input_dir + run_on_id + ".root");
    TFile *f_off = TFile::Open(input_dir + run_off_id + ".root");

    if (!f_on || f_on->IsZombie() || !f_off || f_off->IsZombie())
    {
        std::cerr << "Cannot open input files." << std::endl;
        return;
    }

    TTree *t_on = (TTree *)f_on->Get("events");
    TTree *t_off = (TTree *)f_off->Get("events");

    if (!t_on || !t_off)
    {
        std::cerr << "Cannot find TTree 'events' in files." << std::endl;
        return;
    }

    uint hit_max_on[num_channels];
    uint hit_max_off[num_channels];
    uint hit_ped_on[num_channels];
    uint hit_ped_off[num_channels];

    uint adc_no_cables[num_channels][20];

    t_on->SetBranchAddress("adc", adc_no_cables);

    t_on->SetBranchAddress("hit_max", hit_max_on);
    t_on->SetBranchAddress("hit_pedestal", hit_ped_on);
    t_off->SetBranchAddress("hit_max", hit_max_off);
    t_off->SetBranchAddress("hit_pedestal", hit_ped_off);

    std::vector<double> sum_adc_on(num_channels, 0);
    std::vector<double> sum_adc2_on(num_channels, 0);
    std::vector<double> max_adc_on(num_channels, 0);

    std::vector<double> sum_adc_off(num_channels, 0);
    std::vector<double> sum_adc2_off(num_channels, 0);
    std::vector<double> max_adc_off(num_channels, 0);

    std::vector<double> sum_adc_minus_ped_on(num_channels, 0);
    std::vector<double> sum_adc_minus_ped2_on(num_channels, 0);
    std::vector<double> max_adc_minus_ped_on(num_channels, 0);

    std::vector<double> sum_adc_minus_ped_off(num_channels, 0);
    std::vector<double> sum_adc_minus_ped2_off(num_channels, 0);
    std::vector<double> max_adc_minus_ped_off(num_channels, 0);

    std::vector<int> counts(num_channels, 0);

    Long64_t nEntries = std::min(t_on->GetEntries(), t_off->GetEntries());
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        t_on->GetEntry(i);
        t_off->GetEntry(i);

        for (int ch = 0; ch < num_channels; ++ch)
        {

            // int int_adc = adc_no_cables[ch][12];

            // std::cout << "phase: " << int_adc << std::endl;

            double adc_minus_ped_on = hit_max_on[ch] - hit_ped_on[ch];
            double adc_minus_ped_off = hit_max_off[ch] - hit_ped_off[ch];

            sum_adc_on[ch] += hit_max_on[ch];
            sum_adc2_on[ch] += hit_max_on[ch] * hit_max_on[ch];

            sum_adc_off[ch] += hit_max_off[ch];
            sum_adc2_off[ch] += hit_max_off[ch] * hit_max_off[ch];

            sum_adc_minus_ped_on[ch] += adc_minus_ped_on;
            sum_adc_minus_ped2_on[ch] += adc_minus_ped_on * adc_minus_ped_on;

            sum_adc_minus_ped_off[ch] += adc_minus_ped_off;
            sum_adc_minus_ped2_off[ch] += adc_minus_ped_off * adc_minus_ped_off;

            counts[ch]++;

            if (hit_max_on[ch] > max_adc_on[ch])
            {
                max_adc_on[ch] = hit_max_on[ch];
                max_adc_minus_ped_on[ch] = adc_minus_ped_on;
            }

            if (hit_max_off[ch] > max_adc_off[ch])
            {
                max_adc_off[ch] = hit_max_off[ch];
                max_adc_minus_ped_off[ch] = adc_minus_ped_off;
            }
        }
    }

    std::vector<double> std_dev_on(num_channels, 0);
    std::vector<double> std_dev_off(num_channels, 0);
    std::vector<double> std_dev_minus_ped_on(num_channels, 0);
    std::vector<double> std_dev_minus_ped_off(num_channels, 0);

    for (int ch = 0; ch < num_channels; ++ch)
    {
        if (counts[ch] > 0)
        {
            // std::cout << counts[ch] << " entries for channel " << ch << std::endl;
            double mean_on = sum_adc_on[ch] / counts[ch];
            double variance_on = (sum_adc2_on[ch] / counts[ch]) - (mean_on * mean_on);
            std_dev_on[ch] = sqrt(variance_on);

            double mean_off = sum_adc_off[ch] / counts[ch];
            double variance_off = (sum_adc2_off[ch] / counts[ch]) - (mean_off * mean_off);
            std_dev_off[ch] = sqrt(variance_off);

            double mean_minus_ped_on = sum_adc_minus_ped_on[ch] / counts[ch];
            double variance_minus_ped_on = (sum_adc_minus_ped2_on[ch] / counts[ch]) - (mean_minus_ped_on * mean_minus_ped_on);
            std_dev_minus_ped_on[ch] = sqrt(variance_minus_ped_on);

            double mean_minus_ped_off = sum_adc_minus_ped_off[ch] / counts[ch];
            double variance_minus_ped_off = (sum_adc_minus_ped2_off[ch] / counts[ch]) - (mean_minus_ped_off * mean_minus_ped_off);
            std_dev_minus_ped_off[ch] = sqrt(variance_minus_ped_off);
        }
    }

    for (int fpga = 0; fpga < num_kcu; ++fpga)
    {
        if (protoboards[fpga].substr(2) != injected_pb)
        {
            std::cout << "Skipping protoboard " << protoboards[fpga] << " (not injected)." << std::endl;
            continue;
        }

        avg_adc_injected_channels = 0;
        avg_adc_cross_talk_expected_channels = 0;
        cross_talk = 0;

        TString label = Form("pb_%s", protoboards[fpga].c_str());

        TString out1 = output_dir + run_on_id + "_vs_" + run_off_id + "_protoboard_" + protoboards[fpga].c_str() + "_max_adc_vs_channel.png";
        TString out2 = output_dir + run_on_id + "_vs_" + run_off_id + "_protoboard_" + protoboards[fpga].c_str() + "_max_minus_pedestal_vs_channel.png";

        int ch_start = fpga * channels_per_kcu;

        TCanvas *c1 = new TCanvas(Form("c1_fpga%d_run%d", fpga, run_on), "Max ADC", 3000, 1200);
        c1->SetGrid();
        c1->SetLeftMargin(0.06);
        c1->SetRightMargin(0.02);
        c1->SetTopMargin(0.32);
        c1->SetBottomMargin(0.15);

        TH1I *h1_on = new TH1I(Form("h1_on_fpga%d", fpga), Form("Max ADC - FPGA %d - LED ON", fpga), channels_per_kcu, 0, channels_per_kcu * 2);
        TH1I *h1_off = new TH1I(Form("h1_off_fpga%d", fpga), Form("Max ADC - FPGA %d - LED OFF", fpga), channels_per_kcu, 0, channels_per_kcu * 2);

        h1_on->GetXaxis()->SetTitle("Channel");
        h1_on->GetYaxis()->SetTitle("Max ADC");
        h1_on->SetTitle(Form("Max ADC per Channel - Protoboard %s", protoboards[fpga].c_str()));
        h1_on->SetTitleSize(0.06);
        h1_on->GetYaxis()->SetTitleOffset(1.3);
        h1_on->GetYaxis()->SetTitleOffset(0.6);

        h1_on->GetXaxis()->SetTitleSize(0.043);
        h1_on->GetYaxis()->SetTitleSize(0.043);

        TGraphErrors *g_on = new TGraphErrors();
        TGraphErrors *g_off = new TGraphErrors();

        for (int local_ch = 0; local_ch < channels_per_kcu; ++local_ch)
        {
            if (std::find(ignored_channels.begin(), ignored_channels.end(), local_ch) != ignored_channels.end())
                continue;
            filled_bins.push_back(local_ch + 1); // ROOT bin index is 1-based

            int global_ch = ch_start + local_ch;
            int bin = local_ch + 1;

            h1_on->SetBinContent(bin, max_adc_on[global_ch]);
            h1_on->SetFillColor(kRed + 1);
            h1_on->SetLineColor(kRed + 1);
            h1_on->SetFillStyle(3004);
            h1_on->GetXaxis()->SetBinLabel(bin, Form("%d", shihai_channel_map[local_ch]));

            h1_off->SetBinContent(bin, max_adc_off[global_ch]);
            h1_off->SetFillColor(kBlue + 1);
            h1_off->SetLineColor(kBlue + 1);
            h1_off->SetFillStyle(3005);
            h1_off->GetXaxis()->SetBinLabel(bin, Form("%d", shihai_channel_map[local_ch]));

            int point_on = g_on->GetN();
            double bin_center = h1_on->GetBinCenter(bin);
            g_on->SetPoint(point_on, bin_center, max_adc_on[global_ch]);
            g_on->SetPointError(point_on, 0, std_dev_on[global_ch]);

            int point_off = g_off->GetN();
            g_off->SetPoint(point_off, bin_center, max_adc_off[global_ch]);
            g_off->SetPointError(point_off, 0, std_dev_off[global_ch]);

            if (protoboards[fpga].substr(2) == injected_pb)
            {
                if (asic == 0)
                {
                    if (port == 0) // port A
                    {
                        if (global_ch < (ch_start + 17))
                            avg_adc_injected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                        else if (global_ch > (ch_start + 17) && global_ch < (ch_start + 35))
                            avg_adc_cross_talk_expected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                    }
                    else if (port == 1) // port B
                    {
                        if (global_ch > (ch_start + 17) && global_ch < (ch_start + 35))
                            avg_adc_injected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                        else if (global_ch < (ch_start + 17))
                            avg_adc_cross_talk_expected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                    }
                    else if (port == 2) // port C
                    {
                        if (global_ch > (ch_start + 53) && global_ch < (ch_start + 71))
                            avg_adc_injected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                        else if (global_ch > (ch_start + 35) && global_ch < (ch_start + 53))
                            avg_adc_cross_talk_expected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                    }

                    else if (port == 3) // port D
                    {
                        if (global_ch > (ch_start + 35) && global_ch < (ch_start + 53))
                            avg_adc_injected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                        else if (global_ch > (ch_start + 53) && global_ch < (ch_start + 71))
                            avg_adc_cross_talk_expected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                    }
                }
                else if (asic == 1)
                {
                    if (port == 0) // port A
                    {
                        if (global_ch > (ch_start + 71) && global_ch < (ch_start + 89))
                            avg_adc_injected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                        else if (global_ch > (ch_start + 89) && global_ch < (ch_start + 107))
                            avg_adc_cross_talk_expected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                    }
                    else if (port == 1) // port B
                    {
                        if (global_ch > (ch_start + 89) && global_ch < (ch_start + 107))
                            avg_adc_injected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                        else if (global_ch > (ch_start + 71) && global_ch < (ch_start + 89))
                            avg_adc_cross_talk_expected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                    }
                    else if (port == 2) // port C
                    {
                        if (global_ch > (ch_start + 125) && global_ch < (ch_start + 143))
                            avg_adc_injected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                        else if (global_ch > (ch_start + 107) && global_ch < (ch_start + 125))
                            avg_adc_cross_talk_expected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                    }

                    else if (port == 3) // port D
                    {
                        if (global_ch > (ch_start + 107) && global_ch < (ch_start + 125))
                            avg_adc_injected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                        else if (global_ch > (ch_start + 125) && global_ch < (ch_start + 143))
                            avg_adc_cross_talk_expected_channels += abs(max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch]);
                    }
                }
            }
        }

        for (int i = 1; i <= num_ports; ++i)
        {
            size_t index = i * num_channels_per_port - 1; // 15, 31, 47...
            if (index < filled_bins.size())
                port_bins.push_back(filled_bins[index]);
        }

        if (protoboards[fpga].substr(2) == injected_pb) // Check if the protoboard matches the injected one
        {
            avg_adc_injected_channels /= num_channels_per_port;
            avg_adc_cross_talk_expected_channels /= num_channels_per_port;
            cross_talk = avg_adc_cross_talk_expected_channels / avg_adc_injected_channels * 100;
        }

        h1_on->SetMaximum(1024);
        h1_on->LabelsOption("v", "X");
        h1_on->Draw("HIST");
        h1_off->Draw("HIST SAME");

        g_on->SetMarkerStyle(20);
        g_on->SetMarkerColor(kBlack);
        g_on->SetLineColor(kBlack);
        g_on->Draw("P SAME");

        g_off->SetMarkerStyle(24);
        g_off->SetMarkerColor(kBlack);
        g_off->SetLineColor(kBlack);
        g_off->Draw("P SAME");

        for (size_t i = 0; i < port_bins.size(); ++i)
        {
            double xline = port_bins[i] * 2 + 1;
            TLine *line = new TLine(xline, 0, xline, y_max);
            line->SetLineColor(kGray + 2);
            line->SetLineStyle(1);
            line->SetLineWidth(2);
            line->Draw("same");

            double xstart = (i == 0) ? 0 : port_bins[i - 1] * 2 + 1;
            double xcenter = (xstart + xline) / 2.0;

            TLatex *label = new TLatex(xcenter, y_label, region_labels[i].c_str());
            label->SetTextAlign(22);
            label->SetTextFont(132);
            label->SetTextSize(0.035);
            label->Draw("same");
        }

        // TLegend *leg1 = new TLegend(0.25, 0.75, 0.45, 0.88);
        TLegend *leg1 = new TLegend(0.06, 0.77, 0.26, 0.93);
        leg1->AddEntry(h1_on, "Max ADC - LED ON", "f");
        leg1->AddEntry(h1_off, "Max ADC - LED OFF", "f");
        leg1->AddEntry(g_on, "Std Dev Bars - LED ON", "p");
        leg1->AddEntry(g_off, "Std Dev Bars - LED OFF", "p");
        leg1->Draw();

        // After drawing everything on c1
        // TPaveText *pt1 = new TPaveText(0.25, 0.70, 0.45, 0.75, "NDC");
        TPaveText *pt1 = new TPaveText(0.40, 0.83, 0.60, 0.88, "NDC");
        pt1->AddText(Form("Events: %lld", nEntries));
        pt1->SetFillColor(kWhite);
        pt1->SetShadowColor(0);
        pt1->SetTextFont(42);
        pt1->SetTextSize(0.035);
        pt1->Draw();

        // c1->SaveAs(out1);
        delete c1;

        TCanvas *c2 = new TCanvas(Form("c1_fpga%d_run%d", fpga, run_on), "Max - Pedestal", 3000, 1200);
        c2->SetGrid();
        c2->SetLeftMargin(0.06);
        c2->SetRightMargin(0.02);
        c2->SetTopMargin(0.32);
        c2->SetBottomMargin(0.15);

        TH1I *h2_on = new TH1I(Form("h2_on_fpga%d", fpga), Form("Max - Pedestal - FPGA %d - LED ON", fpga), channels_per_kcu, 0, channels_per_kcu * 2);
        TH1I *h2_off = new TH1I(Form("h2_off_fpga%d", fpga), Form("Max - Pedestal - FPGA %d - LED OFF", fpga), channels_per_kcu, 0, channels_per_kcu * 2);

        h2_on->GetXaxis()->SetTitle("Channel");
        h2_on->GetYaxis()->SetTitle("Max ADC - Pedestal");
        h2_on->SetTitle(Form("Max ADC minus Pedestal per Channel - Protoboard %s", protoboards[fpga].c_str()));
        h2_on->SetTitleSize(0.06);
        h2_on->GetYaxis()->SetTitleOffset(1.3);
        h2_on->GetYaxis()->SetTitleOffset(0.6);

        h2_on->GetXaxis()->SetTitleSize(0.043);
        h2_on->GetYaxis()->SetTitleSize(0.043);

        TGraphErrors *g2_on = new TGraphErrors();
        TGraphErrors *g2_off = new TGraphErrors();

        for (int local_ch = 0; local_ch < channels_per_kcu; ++local_ch)
        {
            if (std::find(ignored_channels.begin(), ignored_channels.end(), local_ch) != ignored_channels.end())
                continue;

            int global_ch = ch_start + local_ch;
            int bin = local_ch + 1;

            h2_on->SetBinContent(bin, max_adc_minus_ped_on[global_ch]);
            h2_on->SetFillColor(kRed + 1);
            h2_on->SetLineColor(kRed + 1);
            h2_on->SetFillStyle(3004);
            h2_on->GetXaxis()->SetBinLabel(bin, Form("%d", shihai_channel_map[local_ch]));

            h2_off->SetBinContent(bin, max_adc_minus_ped_off[global_ch]);
            h2_off->SetFillColor(kBlue + 1);
            h2_off->SetLineColor(kBlue + 1);
            h2_off->SetFillStyle(3005);
            h2_off->GetXaxis()->SetBinLabel(bin, Form("%d", shihai_channel_map[local_ch]));

            int point_on = g2_on->GetN();
            double bin_center = h2_on->GetBinCenter(bin);
            g2_on->SetPoint(point_on, bin_center, max_adc_minus_ped_on[global_ch]);
            g2_on->SetPointError(point_on, 0, std_dev_minus_ped_on[global_ch]);

            int point_off = g2_off->GetN();
            g2_off->SetPoint(point_off, bin_center, max_adc_minus_ped_off[global_ch]);
            g2_off->SetPointError(point_off, 0, std_dev_minus_ped_off[global_ch]);
        }

        h2_on->SetMaximum(1024);
        h2_on->LabelsOption("v", "X");
        h2_on->Draw("HIST");
        h2_off->Draw("HIST SAME");

        g2_on->SetMarkerStyle(20);
        g2_on->SetMarkerColor(kBlack);
        g2_on->SetLineColor(kBlack);
        g2_on->Draw("P SAME");

        g2_off->SetMarkerStyle(24);
        g2_off->SetMarkerColor(kBlack);
        g2_off->SetLineColor(kBlack);
        g2_off->Draw("P SAME");

        for (size_t i = 0; i < port_bins.size(); ++i)
        {
            double xline = port_bins[i] * 2 + 1;
            TLine *line = new TLine(xline, 0, xline, y_max);
            line->SetLineColor(kGray + 2);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            line->Draw("same");

            double xstart = (i == 0) ? 0 : port_bins[i - 1] * 2 + 1;
            double xcenter = (xstart + xline) / 2.0;

            TLatex *label = new TLatex(xcenter, y_label, region_labels[i].c_str());
            label->SetTextAlign(22);
            label->SetTextSize(0.035);
            label->SetTextFont(132);
            label->Draw("same");
        }

        // TLegend *leg2 = new TLegend(0.25, 0.75, 0.45, 0.88);
        TLegend *leg2 = new TLegend(0.06, 0.77, 0.26, 0.93);
        leg2->AddEntry(h2_on, "Max - Pedestal - LED ON", "f");
        leg2->AddEntry(h2_off, "Max - Pedestal - LED OFF", "f");
        leg2->AddEntry(g2_on, "Std Dev Bars - LED ON", "p");
        leg2->AddEntry(g2_off, "Std Dev Bars - LED OFF", "p");
        leg2->Draw();

        // TPaveText *pt2 = new TPaveText(0.25, 0.70, 0.45, 0.75, "NDC");
        TPaveText *pt2 = new TPaveText(0.40, 0.83, 0.60, 0.88, "NDC");
        pt2->AddText(Form("Events: %lld", nEntries));
        pt2->SetFillColor(0);
        pt2->SetTextFont(42);
        pt2->SetShadowColor(0);
        pt2->SetTextSize(0.035);
        pt2->Draw();

        if (protoboards[fpga].substr(2) == injected_pb)
        {
            TString port_names[] = {"A", "B", "C", "D"};
            TString injection_port = port_names[port];
            TString expected_port;

            if (port == 0)
                expected_port = port_names[1]; // B
            else if (port == 1)
                expected_port = port_names[0]; // A
            else if (port == 2)
                expected_port = port_names[3]; // D
            else if (port == 3)
                expected_port = port_names[2]; // C

            TPaveText *pt_xtalk = new TPaveText(0.35, 0.77, 0.65, 0.82, "NDC");
            // TPaveText *pt_xtalk = new TPaveText(0.77, 0.52, 0.97, 0.63, "NDC");
            pt_xtalk->AddText(Form("Cross talk: %.2f%% in port %s for injection in port %s",
                                   cross_talk, expected_port.Data(), injection_port.Data()));
            pt_xtalk->SetFillColor(0);
            pt_xtalk->SetShadowColor(0);
            pt_xtalk->SetTextFont(42);
            pt_xtalk->SetTextSize(0.035);
            pt_xtalk->Draw();
        }

        // c2->SaveAs(out2);
        delete c2;

        // // Create new canvas for difference plot
        // TCanvas *c3 = new TCanvas(Form("c3_fpga%d", fpga), "Difference (ON - OFF)", 2800, 800);
        // c3->SetGrid();
        // c3->SetLeftMargin(0.06);
        // c3->SetRightMargin(0.02);
        // c3->SetTopMargin(0.26);
        // c3->SetBottomMargin(0.15);

        // // Define histogram and graph for difference
        // TH1I *h3_diff = new TH1I(Form("h3_diff_fpga%d", fpga), Form("Difference (ON - OFF) - FPGA %d", fpga), channels_per_kcu, 0, channels_per_kcu * 2);

        // h3_diff->GetXaxis()->SetTitle("Channel");
        // h3_diff->GetYaxis()->SetTitle("Difference: Max ADC (ON - OFF)");
        // h3_diff->SetTitle(Form("Difference per Channel - Protoboard %s", protoboards[fpga].c_str()));
        // h3_diff->SetTitleSize(0.06);
        // h3_diff->GetYaxis()->SetTitleOffset(1.3);
        // h3_diff->GetYaxis()->SetTitleOffset(0.6);
        // h3_diff->GetXaxis()->SetTitleSize(0.043);
        // h3_diff->GetYaxis()->SetTitleSize(0.043);

        // TGraphErrors *g3_diff = new TGraphErrors();

        // for (int local_ch = 0; local_ch < channels_per_kcu; ++local_ch)
        // {
        //     if (std::find(ignored_channels.begin(), ignored_channels.end(), local_ch) != ignored_channels.end())
        //         continue;

        //     int global_ch = ch_start + local_ch;
        //     int bin = local_ch + 1;

        //     // Calculate difference and combined error
        //     double diff = max_adc_minus_ped_on[global_ch] - max_adc_minus_ped_off[global_ch];
        //     double error_diff = sqrt(pow(std_dev_minus_ped_on[global_ch], 2) + pow(std_dev_minus_ped_off[global_ch], 2));

        //     h3_diff->SetBinContent(bin, diff);
        //     h3_diff->SetFillColor(kRed + 1);
        //     h3_diff->SetLineColor(kRed + 1);
        //     h3_diff->SetFillStyle(3004);
        //     h3_diff->GetXaxis()->SetBinLabel(bin, Form("%d", shihai_channel_map[local_ch]));

        //     int point_diff = g3_diff->GetN();
        //     double bin_center = h3_diff->GetBinCenter(bin);
        //     g3_diff->SetPoint(point_diff, bin_center, diff);
        //     g3_diff->SetPointError(point_diff, 0, error_diff);
        // }

        // h3_diff->SetMaximum(1024);
        // h3_diff->LabelsOption("v", "X");
        // h3_diff->Draw("HIST");

        // g3_diff->SetMarkerStyle(21);
        // g3_diff->SetMarkerColor(kBlack);
        // g3_diff->SetLineColor(kBlack);
        // g3_diff->Draw("P SAME");

        // TLegend *leg3 = new TLegend(0.25, 0.75, 0.45, 0.88);
        // leg3->AddEntry(h3_diff, "Difference (ON - OFF)", "f");
        // leg3->AddEntry(g3_diff, "Std Dev Bars", "p");
        // leg3->Draw();

        // TPaveText *pt3 = new TPaveText(0.25, 0.70, 0.45, 0.75, "NDC");
        // pt3->AddText(Form("Events: %lld", nEntries));
        // pt3->SetFillColor(0);
        // pt3->SetTextFont(42);
        // pt3->SetTextSize(0.035);
        // pt3->Draw();

        // // Save plot
        // TString out3 = output_dir + run_on_id + "_vs_" + run_off_id + "_protoboard_" + protoboards[fpga].c_str() + "_difference_led_on_minus_off.png";
        // c3->SaveAs(out3);
    }

    f_on->Close();
    f_off->Close();
}

struct RunConfig
{
    int run_on;
    int run_off;
    int asic;
    int port;
    std::string pb_id;
};

void run_all_tests()
{
    std::vector<RunConfig> tests = {
        {321, 345, 0, 0, "008"},
        {322, 345, 0, 1, "008"},
        {323, 345, 0, 2, "008"},
        {324, 345, 0, 3, "008"},
        {326, 345, 1, 0, "008"},
        {328, 345, 1, 1, "008"},
        {329, 345, 1, 2, "008"},
        {330, 345, 1, 3, "008"},
        {333, 345, 0, 0, "01"},
        {335, 345, 0, 1, "01"},
        {336, 345, 0, 2, "01"},
        {337, 345, 0, 3, "01"},
        {338, 345, 1, 0, "01"},
        {339, 345, 1, 1, "01"},
        {340, 345, 1, 2, "01"},
        {341, 345, 1, 3, "01"},
        {308, 345, 0, 0, "06"},
        {310, 345, 0, 1, "06"},
        {313, 345, 0, 2, "06"},
        {314, 345, 0, 3, "06"},
        {315, 345, 1, 0, "06"},
        {316, 345, 1, 1, "06"},
        {317, 345, 1, 2, "06"},
        {318, 345, 1, 3, "06"},
    };

    for (const auto &test : tests)
    {
        plot_hit_max_tot_led_per_channel(test.run_on, test.run_off, test.asic, test.port, test.pb_id);
    }
}

// Automatically run when executing the file
void run_macro()
{
    run_all_tests();
}