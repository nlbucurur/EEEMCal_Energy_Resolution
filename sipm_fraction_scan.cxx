// To run
// root -l -b
// .L common_led.cxx+
// .L sipm_fraction_scan.cxx+
// sipm_fraction_scan("data","eeemcal_desy_dec2025_mapping_v2.csv","outputs",12,1000000,false,true) // no common mode
// sipm_fraction_scan("data","eeemcal_desy_dec2025_mapping_v2.csv","outputs",12,1000000,true,true) // yes common mode

#include "common_led.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TParameter.h>
#include <TF1.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TError.h>
#include <TSystem.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TROOT.h>
#include <TPad.h>
#include <TLegend.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{

    static float load_adc_per_refunit_local(const char *adc_calib_root, float ref_voltage)
    {
        if (!adc_calib_root)
            return 0.0f;

        TFile *f = TFile::Open(adc_calib_root, "READ");
        if (!f || f->IsZombie())
        {
            if (f)
            {
                f->Close();
                delete f;
            }
            return 0.0f;
        }

        TParameter<float> *p = (TParameter<float> *)f->Get(Form("mean_adc_to_ref_calibration_%.3fV", ref_voltage));
        float val = p ? p->GetVal() : 0.0f;

        f->Close();
        delete f;
        return val;
    }

    static bool fit_peak_gaus_local(TH1 *h,
                                    double &mu,
                                    double &sigma,
                                    double &emu,
                                    double &esigma,
                                    bool draw_fit = false,
                                    Color_t fit_color = kRed)
    {
        mu = sigma = emu = esigma = 0.0;
        if (!h || h->GetEntries() < 30 || h->GetRMS() <= 0)
            return false;

        const double xMin = h->GetXaxis()->GetXmin();
        const double xMax = h->GetXaxis()->GetXmax();
        const double mean0 = h->GetMean();
        const double rms0 = std::max(1e-6, (double)h->GetRMS());

        double fmin = std::max(xMin, mean0 - 3 * rms0);
        double fmax = std::min(xMax, mean0 + 3 * rms0);
        if (fmax <= fmin)
            return false;

        TF1 fit(Form("fit_%s", h->GetName()), "gaus", fmin, fmax);
        fit.SetParameters(h->GetMaximum(), mean0, rms0);

        TFitResultPtr rr = h->Fit(&fit, "RQS0");
        if (!rr.Get() || rr->Status() != 0)
            return false;

        mu = fit.GetParameter(1);
        sigma = std::abs(fit.GetParameter(2));
        emu = fit.GetParError(1);
        esigma = fit.GetParError(2);

        if (draw_fit)
        {
            fit.SetLineColor(fit_color);
            fit.SetLineWidth(2);
            fit.DrawClone("same");
        }

        return (mu > 0.0 && sigma > 0.0);
    }

    static int find_channel_for_sipm(const std::map<int, std::vector<int>> &mapping,
                                     int crystal,
                                     int sipm)
    {
        auto chans = get_crystal_channels(mapping, crystal, LED_SIPMS_PER_CRYSTAL);
        if (sipm < 0 || sipm >= (int)chans.size())
            return -1;
        return chans[sipm];
    }

    static std::vector<int> default_rank_from_means(const std::array<double, LED_SIPMS_PER_CRYSTAL> &means)
    {
        std::vector<int> order(LED_SIPMS_PER_CRYSTAL);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b)
                  {
        if (means[a] == means[b])
            return a < b;
        return means[a] > means[b]; });
        return order;
    }

    static std::string join_order(const std::vector<int> &order)
    {
        std::string out;
        for (size_t i = 0; i < order.size(); ++i)
        {
            if (i)
                out += ",";
            out += std::to_string(order[i]);
        }
        return out;
    }

} // namespace

void sipm_fraction_and_resolution_one(const char *filename,
                                      float voltage,
                                      const char *mapping_csv = "eeemcal_desy_dec2025_mapping_v2.csv",
                                      const char *gain_root = "outputs/gain_match_1.259V.root",
                                      const char *adc_calib_root = "outputs/adc_to_ref_calibration_1.259V.root",
                                      const char *outdir = "outputs",
                                      int central_crystal = 12,
                                      Long64_t max_events = 1000000,
                                      bool use_common_mode = true,
                                      bool reject_tot_events = true)
{
    TH1::AddDirectory(kFALSE);
    gROOT->cd();

    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;
    mkdir_p(outdir);

    g_signal_method = 3;
    g_tot_min = 50;

    auto mapping = read_mapping_csv(mapping_csv, LED_SIPMS_PER_CRYSTAL);
    if (mapping.empty())
    {
        std::cerr << "Error: empty mapping: " << mapping_csv << "\n";
        return;
    }
    if (!mapping.count(central_crystal))
    {
        std::cerr << "Error: central crystal " << central_crystal << " not found in mapping\n";
        return;
    }

    TH1 *gain_factors = nullptr;
    TH1 *crystal_factor = nullptr;
    TFile *gain_handle = nullptr;
    load_gain_factors(gain_root, gain_factors, crystal_factor, gain_handle);

    const float adc_per_refunit = load_adc_per_refunit_local(adc_calib_root, voltage);
    if (!(adc_per_refunit > 0.0f))
    {
        std::cerr << "Warning: could not load mean_adc_to_ref_calibration_" << voltage
                  << "V from " << adc_calib_root
                  << ". The study will stay in ADC units.\n";
    }

    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie())
    {
        std::cerr << "Error: cannot open " << filename << "\n";
        if (f)
        {
            f->Close();
            delete f;
        }
        if (gain_handle)
        {
            gain_handle->Close();
            delete gain_handle;
        }
        return;
    }

    TTree *tree = (TTree *)f->Get("events");
    if (!tree)
    {
        std::cerr << "Error: TTree 'events' not found in " << filename << "\n";
        f->Close();
        delete f;
        if (gain_handle)
        {
            gain_handle->Close();
            delete gain_handle;
        }
        return;
    }

    TBranch *br_adc = tree->GetBranch("adc");
    TBranch *br_tot = tree->GetBranch("tot");
    if (!br_adc)
    {
        std::cerr << "Error: missing adc branch\n";
        f->Close();
        delete f;
        if (gain_handle)
        {
            gain_handle->Close();
            delete gain_handle;
        }
        return;
    }

    TLeaf *leaf_adc = br_adc->GetLeaf("adc");
    TLeaf *leaf_tot = (br_tot ? br_tot->GetLeaf("tot") : nullptr);
    if (!leaf_adc)
    {
        std::cerr << "Error: missing adc leaf\n";
        f->Close();
        delete f;
        if (gain_handle)
        {
            gain_handle->Close();
            delete gain_handle;
        }
        return;
    }

    const int n_adc = leaf_adc->GetLen();
    if (n_adc <= 0 || (n_adc % LED_SAMPLES_PER_CHANNEL) != 0)
    {
        std::cerr << "Error: unexpected adc length " << n_adc << "\n";
        f->Close();
        delete f;
        if (gain_handle)
        {
            gain_handle->Close();
            delete gain_handle;
        }
        return;
    }

    const int n_channels = n_adc / LED_SAMPLES_PER_CHANNEL;
    const bool have_tot = (leaf_tot && leaf_tot->GetLen() == n_adc);

    std::vector<uint32_t> adc_buf((size_t)n_adc, 0);
    std::vector<uint32_t> tot_buf;
    if (have_tot)
        tot_buf.resize((size_t)n_adc, 0);

    tree->SetBranchAddress("adc", adc_buf.data());
    if (have_tot)
        tree->SetBranchAddress("tot", tot_buf.data());

    auto adc_ptr = [&](int ch) -> uint32_t *
    {
        return &adc_buf[(size_t)ch * LED_SAMPLES_PER_CHANNEL];
    };
    auto tot_ptr = [&](int ch) -> uint32_t *
    {
        return &tot_buf[(size_t)ch * LED_SAMPLES_PER_CHANNEL];
    };

    int common_mode_channels_by_crystal[LED_MAX_NUM_CRYSTALS];
    for (int cr = 0; cr < LED_MAX_NUM_CRYSTALS; ++cr)
        common_mode_channels_by_crystal[cr] = -1;
    common_mode_channels_by_crystal[12] = 44;

    std::array<int, LED_SIPMS_PER_CRYSTAL> chans;
    chans.fill(-1);
    for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
        chans[sipm] = find_channel_for_sipm(mapping, central_crystal, sipm);

    std::array<TH1F *, LED_SIPMS_PER_CRYSTAL> h_share;
    std::array<TH1F *, LED_SIPMS_PER_CRYSTAL> h_signal_adc;
    std::array<TH1F *, LED_SIPMS_PER_CRYSTAL> h_signal_veq;
    for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
    {
        h_share[sipm] = new TH1F(Form("h_share_sipm_%02d", sipm),
                                 Form("SiPM %02d fraction of full 16-SiPM sum;fraction;Events", sipm),
                                 50, 0, 0.2);
        h_signal_adc[sipm] = new TH1F(Form("h_signal_adc_sipm_%02d", sipm),
                                      Form("SiPM %02d signal;ADC;Events", sipm),
                                      500, 0, 5000);
        h_signal_veq[sipm] = new TH1F(Form("h_signal_veq_sipm_%02d", sipm),
                                      Form("SiPM %02d signal;V-equiv;Events", sipm),
                                      500, 0, 5.0);
    }

    TH1F *h_full_adc = new TH1F("h_full16_adc", "Full 16-SiPM central sum;ADC;Events", 600, 0, 40000);
    TH1F *h_full_veq = new TH1F("h_full16_veq", "Full 16-SiPM central sum;V-equiv;Events", 600, 0, 5.0);

    Long64_t nentries = tree->GetEntries();
    if (max_events > 0 && max_events < nentries)
        nentries = max_events;

    Long64_t used_events = 0;
    Long64_t skipped_tot = 0;
    Long64_t skipped_zero = 0;

    std::vector<std::array<float, LED_SIPMS_PER_CRYSTAL>> kept_event_signals;
    kept_event_signals.reserve((size_t)std::min<Long64_t>(nentries, 200000));

    for (Long64_t ev = 0; ev < nentries; ++ev)
    {
        tree->GetEntry(ev);

        float cm = 0.0f;
        bool have_cm = false;
        if (use_common_mode)
        {
            int cm_ch = common_mode_channels_by_crystal[central_crystal];
            if (cm_ch >= 0 && cm_ch < n_channels)
            {
                uint32_t *cm_adc = adc_ptr(cm_ch);
                cm = 0.5f * ((float)cm_adc[0] + (float)cm_adc[1]);
                have_cm = true;
            }
        }

        std::array<float, LED_SIPMS_PER_CRYSTAL> sig_adc_evt{};
        bool event_has_tot = false;
        float sum16_adc = 0.0f;

        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
        {
            const int ch = chans[sipm];
            if (ch < 0 || ch >= n_channels)
                continue;

            uint32_t *adc_raw = adc_ptr(ch);
            uint32_t adc_corr[LED_SAMPLES_PER_CHANNEL];
            uint32_t *adc_used = adc_raw;
            if (have_cm)
            {
                subtract_common_mode(adc_corr, adc_raw, cm);
                adc_used = adc_corr;
            }

            if (have_tot && ::has_tot(tot_ptr(ch)))
            {
                event_has_tot = true;
                if (reject_tot_events)
                    break;
            }

            const float gch = gain_factors ? (float)gain_factors->GetBinContent(central_crystal * LED_SIPMS_PER_CRYSTAL + sipm + 1) : 1.0f;
            float sig = calculate_signal_adc(adc_used, gch);
            sig_adc_evt[sipm] = sig;
            sum16_adc += sig;
        }

        if (event_has_tot && reject_tot_events)
        {
            skipped_tot++;
            continue;
        }

        if (!(sum16_adc > 0.0f))
        {
            skipped_zero++;
            continue;
        }

        const float gcr = crystal_factor ? (float)crystal_factor->GetBinContent(central_crystal + 1) : 1.0f;
        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
            sig_adc_evt[sipm] *= gcr;
        sum16_adc *= gcr;

        h_full_adc->Fill(sum16_adc);
        if (adc_per_refunit > 0.0f)
            h_full_veq->Fill((sum16_adc / adc_per_refunit) * voltage);

        for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
        {
            const float sig = sig_adc_evt[sipm];
            const float frac = sig / sum16_adc;
            h_signal_adc[sipm]->Fill(sig);
            if (adc_per_refunit > 0.0f)
                h_signal_veq[sipm]->Fill((sig / adc_per_refunit) * voltage);
            h_share[sipm]->Fill(frac);
        }

        kept_event_signals.push_back(sig_adc_evt);
        used_events++;
    }

    std::array<double, LED_SIPMS_PER_CRYSTAL> mu_share_fit{};
    std::array<double, LED_SIPMS_PER_CRYSTAL> sigma_share_fit{};
    std::array<double, LED_SIPMS_PER_CRYSTAL> emu_share_fit{};
    std::array<double, LED_SIPMS_PER_CRYSTAL> esigma_share_fit{};
    for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
    {
        double mu = 0.0, sigma = 0.0, emu = 0.0, esigma = 0.0;
        if (fit_peak_gaus_local(h_share[sipm], mu, sigma, emu, esigma))
        {
            mu_share_fit[sipm] = mu;
            sigma_share_fit[sipm] = sigma;
            emu_share_fit[sipm] = emu;
            esigma_share_fit[sipm] = esigma;
        }
    }

    std::vector<int> order = default_rank_from_means(mu_share_fit);

    std::array<TH1F *, LED_SIPMS_PER_CRYSTAL> h_subset_adc;
    std::array<TH1F *, LED_SIPMS_PER_CRYSTAL> h_subset_veq;
    std::array<TH1F *, LED_SIPMS_PER_CRYSTAL> h_subset_frac;
    for (int n = 0; n < LED_SIPMS_PER_CRYSTAL; ++n)
    {
        h_subset_adc[n] = new TH1F(Form("h_subset_adc_top%d", n + 1),
                                   Form("Top-%d cumulative subset;ADC;Events", n + 1),
                                   600, 0, 40000);
        h_subset_veq[n] = new TH1F(Form("h_subset_veq_top%d", n + 1),
                                   Form("Top-%d cumulative subset;V-equiv;Events", n + 1),
                                   600, 0, 2.0);
        h_subset_frac[n] = new TH1F(Form("h_subset_frac_top%d", n + 1),
                                    Form("Top-%d cumulative fraction;fraction;Events", n + 1),
                                    250, 0, 1.05);
    }

    for (const auto &evt : kept_event_signals)
    {
        float full = 0.0f;
        for (int s = 0; s < LED_SIPMS_PER_CRYSTAL; ++s)
            full += evt[s];

        if (!(full > 0.0f))
            continue;

        float acc = 0.0f;
        for (int n = 0; n < LED_SIPMS_PER_CRYSTAL; ++n)
        {
            acc += evt[order[n]];
            h_subset_adc[n]->Fill(acc);
            if (adc_per_refunit > 0.0f)
                h_subset_veq[n]->Fill((acc / adc_per_refunit) * voltage);
            h_subset_frac[n]->Fill(acc / full);
        }
    }

    std::array<double, LED_SIPMS_PER_CRYSTAL> mu_cumshare_fit{};
    std::array<double, LED_SIPMS_PER_CRYSTAL> sigma_cumshare_fit{};
    std::array<double, LED_SIPMS_PER_CRYSTAL> emu_cumshare_fit{};
    std::array<double, LED_SIPMS_PER_CRYSTAL> esigma_cumshare_fit{};
    for (int n = 0; n < LED_SIPMS_PER_CRYSTAL; ++n)
    {
        double mu = 0.0, sigma = 0.0, emu = 0.0, esigma = 0.0;
        if (fit_peak_gaus_local(h_subset_frac[n], mu, sigma, emu, esigma))
        {
            mu_cumshare_fit[n] = mu;
            sigma_cumshare_fit[n] = sigma;
            emu_cumshare_fit[n] = emu;
            esigma_cumshare_fit[n] = esigma;
        }
    }

    TGraphErrors *g_share = new TGraphErrors();
    g_share->SetName("g_mean_share_by_sipm");
    g_share->SetTitle("Fit mean fraction collected by each SiPM;SiPM;Fit mean fraction of full 16-SiPM sum");

    TGraphErrors *g_share_rank = new TGraphErrors();
    g_share_rank->SetName("g_mean_share_ranked");
    g_share_rank->SetTitle("Fit mean fraction ranked;Rank;Fit mean fraction");

    TGraphErrors *g_cumshare = new TGraphErrors();
    g_cumshare->SetName("g_cumulative_mean_share");
    g_cumshare->SetTitle("Fit mean cumulative fraction vs number of SiPMs;Number of SiPMs included;Fit mean collected fraction");

    TGraphErrors *g_res_n = new TGraphErrors();
    g_res_n->SetName("g_resolution_vs_n");
    g_res_n->SetTitle("Resolution vs number of included SiPMs;Number of SiPMs included;Resolution (%)");

    TGraphErrors *g_mu_n = new TGraphErrors();
    g_mu_n->SetName("g_mean_vs_n");
    g_mu_n->SetTitle("Mean reconstructed central signal vs number of included SiPMs;Number of SiPMs included;Mean (V-equiv)");

    for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
    {
        g_share->SetPoint(sipm, sipm, mu_share_fit[sipm]);
        g_share->SetPointError(sipm, 0.0, emu_share_fit[sipm]);

        const int ranked_sipm = order[sipm];
        g_share_rank->SetPoint(sipm, sipm + 1, mu_share_fit[ranked_sipm]);
        g_share_rank->SetPointError(sipm, 0.0, emu_share_fit[ranked_sipm]);

        g_cumshare->SetPoint(sipm, sipm + 1, mu_cumshare_fit[sipm]);
        g_cumshare->SetPointError(sipm, 0.0, emu_cumshare_fit[sipm]);

        double mu = 0.0, sigma = 0.0, emu = 0.0, esigma = 0.0;
        TH1 *href = (adc_per_refunit > 0.0f) ? (TH1 *)h_subset_veq[sipm] : (TH1 *)h_subset_adc[sipm];
        if (fit_peak_gaus_local(href, mu, sigma, emu, esigma))
        {
            const double res = 100.0 * sigma / mu;
            double eres = 0.0;
            if (mu > 0.0 && sigma > 0.0)
            {
                const double rel2 = (emu > 0.0 ? (emu / mu) * (emu / mu) : 0.0) +
                                    (esigma > 0.0 ? (esigma / sigma) * (esigma / sigma) : 0.0);
                eres = res * std::sqrt(rel2);
            }
            g_res_n->SetPoint(sipm, sipm + 1, res);
            g_res_n->SetPointError(sipm, 0.0, eres);
            g_mu_n->SetPoint(sipm, sipm + 1, mu);
            g_mu_n->SetPointError(sipm, 0.0, emu);
        }
    }

    TGraphErrors *g_res_vs_mu = new TGraphErrors();
    g_res_vs_mu->SetName("g_resolution_vs_mu");
    g_res_vs_mu->SetTitle("Resolution vs fitted mean signal;Mean (V-equiv);Resolution (%)");

    TGraphErrors *g_res_sqrt_mu = new TGraphErrors();
    g_res_sqrt_mu->SetName("g_resolution_times_sqrtmu_vs_n");
    g_res_sqrt_mu->SetTitle("Resolution#times#sqrt{mean} vs number of included SiPMs;Number of SiPMs included;Resolution (%)#times#sqrt{mean}");

    int ip_mu = 0;
    for (int i = 0; i < g_mu_n->GetN() && i < g_res_n->GetN(); ++i)
    {
        double n_mu = 0.0, mu = 0.0;
        double n_res = 0.0, res = 0.0;
        g_mu_n->GetPoint(i, n_mu, mu);
        g_res_n->GetPoint(i, n_res, res);

        const double emu = g_mu_n->GetErrorY(i);
        const double eres = g_res_n->GetErrorY(i);

        if (mu <= 0.0 || res <= 0.0)
            continue;

        g_res_vs_mu->SetPoint(ip_mu, mu, res);
        g_res_vs_mu->SetPointError(ip_mu, emu, eres);

        const double y = res * std::sqrt(mu);
        double ey = 0.0;
        if (res > 0.0 && mu > 0.0)
        {
            const double rel2 = (eres > 0.0 ? (eres / res) * (eres / res) : 0.0) +
                                (emu > 0.0 ? 0.25 * (emu / mu) * (emu / mu) : 0.0);
            ey = y * std::sqrt(rel2);
        }

        g_res_sqrt_mu->SetPoint(ip_mu, n_mu, y);
        g_res_sqrt_mu->SetPointError(ip_mu, 0.0, ey);

        ++ip_mu;
    }

    TF1 *f_res_mu = nullptr;
    double R2_res_mu = 0.0;
    double chi2ndf_res_mu = 0.0;

    if (g_res_vs_mu->GetN() >= 2)
    {
        double x0 = 1.0, y0 = 0.0;
        g_res_vs_mu->GetPoint(0, x0, y0);

        double xmin = x0, xmax = x0;
        for (int i = 1; i < g_res_vs_mu->GetN(); ++i)
        {
            double x = 0.0, y = 0.0;
            g_res_vs_mu->GetPoint(i, x, y);
            xmin = std::min(xmin, x);
            xmax = std::max(xmax, x);
        }

        f_res_mu = new TF1("f_res_mu", "sqrt([0]/x + [1]*[1])", xmin, xmax);
        f_res_mu->SetParNames("A", "B");
        f_res_mu->SetParameters(y0 * y0 * x0, 0.0); // simple starting guess
        f_res_mu->SetLineColor(kBlue);
        f_res_mu->SetLineWidth(2);

        g_res_vs_mu->Fit(f_res_mu, "RS");

        if (f_res_mu && g_res_vs_mu->GetN() >= 2)
        {
            const int npts = g_res_vs_mu->GetN();

            double ymean = 0.0;
            for (int i = 0; i < npts; ++i)
            {
                double x = 0.0, y = 0.0;
                g_res_vs_mu->GetPoint(i, x, y);
                ymean += y;
            }
            ymean /= (double)npts;

            double ss_res = 0.0;
            double ss_tot = 0.0;

            for (int i = 0; i < npts; ++i)
            {
                double x = 0.0, y = 0.0;
                g_res_vs_mu->GetPoint(i, x, y);

                const double yfit = f_res_mu->Eval(x);
                ss_res += (y - yfit) * (y - yfit);
                ss_tot += (y - ymean) * (y - ymean);
            }

            if (ss_tot > 0.0)
                R2_res_mu = 1.0 - ss_res / ss_tot;

            // if (f_res_mu->GetNDF() > 0)
            chi2ndf_res_mu = f_res_mu->GetChisquare() / f_res_mu->GetNDF();
        }
    }

    TGraph *g_model_res_vs_n = new TGraph();
    g_model_res_vs_n->SetName("g_model_resolution_vs_n");

    if (f_res_mu)
    {
        for (int i = 0; i < g_mu_n->GetN(); ++i)
        {
            double n = 0.0, mu = 0.0;
            g_mu_n->GetPoint(i, n, mu);
            if (mu <= 0.0)
                continue;

            const double y_model = f_res_mu->Eval(mu);
            g_model_res_vs_n->SetPoint(g_model_res_vs_n->GetN(), n, y_model);
        }

        g_model_res_vs_n->SetLineColor(kBlue);
        g_model_res_vs_n->SetLineWidth(2);
    }

    const int run = extract_run_number(filename);

    const TString run_voltage_label = Form("Run %03d, %.3f V", run, voltage);

    auto draw_run_voltage = [&](double x = 0.4, double y = 0.92, double size = 0.032)
    {
        TLatex t;
        t.SetNDC();
        t.SetTextFont(42);
        t.SetTextSize(size);
        t.DrawLatex(x, y, run_voltage_label);
    };

    const std::string base = Form("%s/sipm_fraction_run%03d_%.3fV", outdir, run, voltage);
    const std::string pdf = base + ".pdf";
    const std::string rootname = base + ".root";
    const std::string txt = base + "_ranking.txt";

    std::ofstream ofs(txt.c_str());
    ofs << "# sipm fit_mean fit_sigma fit_mean_err fit_sigma_err rank\n";
    for (int rank = 0; rank < LED_SIPMS_PER_CRYSTAL; ++rank)
    {
        const int sipm = order[rank];
        ofs << sipm << " " << mu_share_fit[sipm] << " " << sigma_share_fit[sipm] << " "
            << emu_share_fit[sipm] << " " << esigma_share_fit[sipm] << " " << (rank + 1) << "\n";
    }
    ofs << "# order " << join_order(order) << "\n";
    ofs.close();

    TCanvas *c = new TCanvas("c_sipm_fraction", "", 1200, 800);
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.035);

    c->SaveAs((pdf + "[").c_str());

    c->Clear();
    c->Divide(4, 4);
    for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
    {
        c->cd(sipm + 1);
        h_share[sipm]->Draw("hist");
        double mu0 = 0.0, sigma0 = 0.0, emu0 = 0.0, esigma0 = 0.0;
        fit_peak_gaus_local(h_share[sipm], mu0, sigma0, emu0, esigma0, true, kRed);
        lat.DrawLatex(0.5, 0.84, Form("fit mean = %.4f", mu_share_fit[sipm]));
        lat.DrawLatex(0.5, 0.77, Form("fit #sigma = %.4f", sigma_share_fit[sipm]));
        lat.DrawLatex(0.5, 0.70, Form("rank = %d", (int)(std::find(order.begin(), order.end(), sipm) - order.begin()) + 1));
        draw_run_voltage();
    }
    c->SaveAs(pdf.c_str());

    c->Clear();
    g_share->SetMarkerStyle(20);
    c->SetGrid();
    g_share->Draw("AP");
    // draw a horizontal line at 1/16
    TF1 *f_line = new TF1("f_line", "1.0/16.0", 0, LED_SIPMS_PER_CRYSTAL + 1);
    f_line->SetLineColor(kOrange + 7);
    f_line->SetLineStyle(2);
    f_line->Draw("same");
    // add legend entry for the line
    TLegend *leg = new TLegend(0.6, 0.8, 0.9, 0.9);
    leg->AddEntry(f_line, "Ideal equal share (1/16)", "l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();
    draw_run_voltage();
    c->SaveAs(pdf.c_str());

    c->Clear();
    g_share_rank->SetMarkerStyle(20);
    c->SetGrid();
    g_share_rank->Draw("AP");
    f_line->Draw("same");
    leg->Draw();
    draw_run_voltage();
    c->SaveAs(pdf.c_str());

    c->Clear();
    g_cumshare->SetMarkerStyle(20);
    c->SetGrid();
    g_cumshare->Draw("AP");
    draw_run_voltage();
    c->SaveAs(pdf.c_str());

    c->Clear();
    g_res_n->SetMarkerStyle(20);
    c->SetGrid();
    g_res_n->Draw("AP");
    if (g_model_res_vs_n->GetN() > 0)
        g_model_res_vs_n->Draw("L SAME");

    if (f_res_mu)
    {
        lat.SetTextColor(kBlue);
        lat.DrawLatex(0.5, 0.84, "R(mean) = #sqrt{#frac{A}{mean} + B^{2}}");
        lat.DrawLatex(0.5, 0.78, Form("A = %.4g #pm %.2g", f_res_mu->GetParameter(0), f_res_mu->GetParError(0)));
        lat.DrawLatex(0.5, 0.72, Form("B = %.4g #pm %.2g", f_res_mu->GetParameter(1), f_res_mu->GetParError(1)));
        lat.DrawLatex(0.5, 0.66, Form("R^{2} = %.4f", R2_res_mu));
        lat.DrawLatex(0.5, 0.60, Form("#chi^{2}/NDF = %.3f", chi2ndf_res_mu));
        lat.SetTextColor(kBlack);
        draw_run_voltage();
    }
    c->SaveAs(pdf.c_str());

    c->Clear();
    g_res_sqrt_mu->SetMarkerStyle(20);
    g_res_sqrt_mu->Draw("AP");
    // lat.DrawLatex(0.18, 0.84, "If this is roughly flat, the behavior is close to pure statistics");
    draw_run_voltage();
    c->SaveAs(pdf.c_str());

    c->Clear();
    g_mu_n->SetMarkerStyle(20);
    c->SetGrid();
    g_mu_n->Draw("AP");
    draw_run_voltage();
    c->SaveAs(pdf.c_str());

    c->Clear();
    g_res_vs_mu->SetMarkerStyle(20);
    g_res_vs_mu->Draw("AP");

    if (f_res_mu)
        f_res_mu->Draw("SAME");

    if (f_res_mu)
    {
        lat.SetTextColor(kBlue);
        lat.DrawLatex(0.5, 0.84, "R(mean) = #sqrt{#frac{A}{mean} + B^{2}}");
        lat.DrawLatex(0.5, 0.78, Form("A = %.4g #pm %.2g", f_res_mu->GetParameter(0), f_res_mu->GetParError(0)));
        lat.DrawLatex(0.5, 0.72, Form("B = %.4g #pm %.2g", f_res_mu->GetParameter(1), f_res_mu->GetParError(1)));
        lat.SetTextColor(kBlack);
        draw_run_voltage();
    }
    c->SaveAs(pdf.c_str());

    c->Clear();
    c->Divide(4, 4);
    for (int n = 0; n < LED_SIPMS_PER_CRYSTAL; ++n)
    {
        c->cd(n + 1);
        TH1 *href = (adc_per_refunit > 0.0f) ? (TH1 *)h_subset_veq[n] : (TH1 *)h_subset_adc[n];
        href->Draw("hist");

        double mu = 0.0, sigma = 0.0, emu = 0.0, esigma = 0.0;
        if (fit_peak_gaus_local(href, mu, sigma, emu, esigma, true, kRed))
        {
            lat.SetTextColor(kRed);
            lat.DrawLatex(0.5, 0.70, Form("mean = %.3f", mu));
            lat.DrawLatex(0.5, 0.63, Form("#sigma = %.3f", sigma));
            lat.DrawLatex(0.5, 0.56, Form("res = %.2f%%", 100.0 * sigma / mu));
            lat.SetTextColor(kBlack);
            draw_run_voltage();
        }

        lat.DrawLatex(0.5, 0.84, Form("N = %d", n + 1));
        lat.DrawLatex(0.5, 0.77, Form("add SiPM %d", order[n]));
    }
    c->SaveAs(pdf.c_str());

    c->Clear();
    TH1 *hfull = (adc_per_refunit > 0.0f) ? (TH1 *)h_full_veq : (TH1 *)h_full_adc;
    hfull->Draw("hist");
    draw_run_voltage();
    c->SetGrid();

    double mu = 0.0, sigma = 0.0, emu = 0.0, esigma = 0.0;
    if (fit_peak_gaus_local(hfull, mu, sigma, emu, esigma, true, kRed))
    {
        lat.SetTextColor(kRed);
        lat.DrawLatex(0.4, 0.54, Form("mean = %.3f", mu));
        lat.DrawLatex(0.4, 0.49, Form("#sigma = %.3f", sigma));
        lat.DrawLatex(0.4, 0.44, Form("res = %.2f%%", 100.0 * sigma / mu));
        lat.SetTextColor(kBlack);
    }

    lat.DrawLatex(0.4, 0.84, Form("Run %03d", run));
    lat.DrawLatex(0.4, 0.79, Form("Voltage label %.3f V", voltage));
    lat.DrawLatex(0.4, 0.74, Form("Used events %lld", used_events));
    lat.DrawLatex(0.4, 0.69, Form("Skipped ToT %lld", skipped_tot));
    lat.DrawLatex(0.4, 0.64, Form("Skipped zero-sum %lld", skipped_zero));
    lat.DrawLatex(0.4, 0.59, Form("Ranking %s", join_order(order).c_str()));
    draw_run_voltage();
    c->SaveAs(pdf.c_str());

    c->SaveAs((pdf + "]").c_str());

    TFile *fout = TFile::Open(rootname.c_str(), "RECREATE");
    for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
    {
        h_share[sipm]->Write();
        h_signal_adc[sipm]->Write();
        h_signal_veq[sipm]->Write();
        h_subset_adc[sipm]->Write();
        h_subset_veq[sipm]->Write();
        h_subset_frac[sipm]->Write();
    }
    h_full_adc->Write();
    h_full_veq->Write();
    g_share->Write();
    g_share_rank->Write();
    g_cumshare->Write();
    g_res_n->Write();
    g_mu_n->Write();
    g_res_vs_mu->Write();
    g_res_sqrt_mu->Write();
    g_model_res_vs_n->Write();
    if (f_res_mu)
        f_res_mu->Write("f_resolution_vs_mu_model");

    TParameter<int>("run", run).Write();
    TParameter<float>("voltage_label_V", voltage).Write();
    TParameter<int>("central_crystal", central_crystal).Write();
    TParameter<Long64_t>("used_events", used_events).Write();
    TParameter<Long64_t>("skipped_tot", skipped_tot).Write();
    TParameter<Long64_t>("skipped_zero", skipped_zero).Write();
    TParameter<int>("reject_tot_events", (int)reject_tot_events).Write();

    for (int rank = 0; rank < LED_SIPMS_PER_CRYSTAL; ++rank)
    {
        TParameter<int>(Form("rank_%02d_sipm", rank + 1), order[rank]).Write();
        TParameter<float>(Form("sipm_%02d_mean_fraction", rank + 1), (float)mu_share_fit[order[rank]]).Write();
    }

    fout->Close();
    delete fout;

    std::cout << "Saved PDF:  " << pdf << "\n";
    std::cout << "Saved ROOT: " << rootname << "\n";
    std::cout << "Saved TXT:  " << txt << "\n";
    std::cout << "Used events=" << used_events
              << " skipped_tot=" << skipped_tot
              << " skipped_zero=" << skipped_zero << "\n";
    std::cout << "Ranking (largest mean fraction first): " << join_order(order) << "\n";

    delete c;

    for (int sipm = 0; sipm < LED_SIPMS_PER_CRYSTAL; ++sipm)
    {
        delete h_share[sipm];
        delete h_signal_adc[sipm];
        delete h_signal_veq[sipm];
        delete h_subset_adc[sipm];
        delete h_subset_veq[sipm];
        delete h_subset_frac[sipm];
    }
    delete h_full_adc;
    delete h_full_veq;
    delete g_share;
    delete g_share_rank;
    delete g_cumshare;
    delete g_res_n;
    delete g_mu_n;
    delete g_res_vs_mu;
    delete g_res_sqrt_mu;
    delete g_model_res_vs_n;
    delete f_res_mu;

    f->Close();
    delete f;

    if (gain_handle)
    {
        gain_handle->Close();
        delete gain_handle;
    }
    else
    {
        if (gain_factors && std::string(gain_factors->GetName()).find("_unity") != std::string::npos)
            delete gain_factors;
        if (crystal_factor && std::string(crystal_factor->GetName()).find("_unity") != std::string::npos)
            delete crystal_factor;
    }
}

void sipm_fraction_scan(const char *data_dir = "data",
                        const char *mapping_csv = "eeemcal_desy_dec2025_mapping_v2.csv",
                        const char *outdir = "outputs",
                        int central_crystal = 12,
                        Long64_t max_events = 1000000,
                        bool use_common_mode = true,
                        bool reject_tot_events = true)
{
    std::vector<std::pair<int, float>> runs = {
        {170, 1.25f},
        {171, 1.259f},
        {172, 1.268f},
        {173, 1.277f} /*,
         {174, 1.286f},
         {175, 1.295f},
         {176, 1.304f},
         {177, 1.313f},
         {178, 1.322f},
         {179, 1.331f},
         {180, 1.340f}*/
    };

    for (const auto &rv : runs)
    {
        const std::string infile = Form("%s/Run%03d.root", data_dir, rv.first);
        const std::string gain_root = Form("%s/gain_match_%.3fV.root", outdir, rv.second);
        const std::string adc_root = Form("%s/adc_to_ref_calibration_%.3fV.root", outdir, rv.second);

        sipm_fraction_and_resolution_one(infile.c_str(),
                                         rv.second,
                                         mapping_csv,
                                         gain_root.c_str(),
                                         adc_root.c_str(),
                                         outdir,
                                         central_crystal,
                                         max_events,
                                         use_common_mode,
                                         reject_tot_events);
    }
}
