// energy_resolution.cxx (LED/voltage version)
//
// Refactor of Tristan's energy_resolution.cxx to work with the EEEMCal LED scans.
// It follows the same style as led_analysis.C / gain_match.cxx / adc_calibration.cxx / tot_calibration.cxx:
//   - Uses mapping CSV produced by generate_mapping.py
//   - Uses common_led.{h,cxx} for signal extraction
//   - Uses gain_match_XX.XXV.root for per-channel gain factors
//   - Uses adc_to_ref_calibration_1.27V.root to convert ADC-sum to "voltage-equivalent"
//   - Optionally uses tot_calibration_values.root to convert ToT-sum to "voltage-equivalent" for saturated (central) events
//
// Output:
//   outputs/energy_resolution_led.pdf
//   outputs/energy_resolution_led.root
//
// Run examples:
//   root -l -b
//   .L common_led.cxx+
//   .L energy_resolution.cxx+
//   energy_resolution()

//   root -l -q 'common_led.cxx+ energy_resolution.cxx+ energy_resolution_led_scan("data","eeemcal_desy_dec2025_mapping_v2.csv","outputs")'
//
// Notes:
//   - By default, any event with ToT in a NON-central crystal is rejected (same as tot_calibration.cxx).
//     This keeps the reconstruction simple and consistent.
//   - Central-crystal ToT is handled via tot_calibration_values.root (pol2 / pol1). Mixed (ADC+ToT) in the
//     central crystal is allowed; we scale the ToT-sum to 16 channels before conversion, then scale back.
//

#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <istream>
#include <iosfwd>
#include <sstream>

#include "common_led.h"

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

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

static constexpr int SAMPLES_PER_CHANNEL = 20;
static constexpr int SIPMS_PER_CRYSTAL = 16;
static constexpr int MAX_NUM_CRYSTALS = 25;

// -------------------- knobs --------------------
static long g_max_events = 1000000;
static bool g_use_tot_for_central = true;                 // use tot_calibration_values.root for central ToT
static bool g_reject_noncentral_tot = true;               // reject events with ToT in non-central crystals
static bool g_use_totcalib_resolution_for_totonly = true; // if central is ToT-only, take resolution from tot_calibration_values.root

static bool g_force_tot_override_for_voltage = true; // force ToT override for a specific LED voltage label
// Force ToT override for specific LED labels (in volts)
static std::vector<float> g_force_tot_override_voltages = {1.29f};
static float g_vlabel_eps = 5e-3f; // tolerance on voltage label (V)

static bool force_tot_override(float Vlabel)
{
    if (!g_force_tot_override_for_voltage)
        return false;

    for (float v : g_force_tot_override_voltages)
        if (std::fabs(Vlabel - v) < g_vlabel_eps)
            return true;
    return false;
}

// COG ellipse cut (same defaults as tot_calibration.cxx)
static bool g_do_cog_ellipse_cut = true;
static float g_cog_cx = 2.0f;
static float g_cog_cy = 2.0f;
static float g_cog_sx = 0.80f;
static float g_cog_sy = 0.80f;

// Ignore some crystals in sums/COG (e.g. 9 weird, 15 can be problematic)
static std::vector<int> g_skip_crystals = {9};

static inline bool is_skipped_crystal(int cr)
{
    return std::find(g_skip_crystals.begin(), g_skip_crystals.end(), cr) != g_skip_crystals.end();
}

// -------------------- helpers --------------------

static void mkdir_p(const char *dir)
{
    if (!dir)
        return;
    gSystem->mkdir(dir, true);
}

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

// Ellipse cut around (cx,cy) with widths (sx,sy)
static inline bool ellipse_cut(float x, float y)
{
    if (!g_do_cog_ellipse_cut)
        return true;
    const float dx = (x - g_cog_cx) / g_cog_sx;
    const float dy = (y - g_cog_cy) / g_cog_sy;
    return (dx * dx + dy * dy) <= 1.0f;
}

// COG from 5x5 crystal array, returns whether event passes ellipse cut.
static bool calculate_cog(TH2 *distribution, const std::array<float, MAX_NUM_CRYSTALS> &values)
{
    float total_signal = 0.0f;
    for (int i = 0; i < MAX_NUM_CRYSTALS; ++i)
    {
        if (is_skipped_crystal(i))
            continue;
        total_signal += values[i];
    }
    if (total_signal <= 0)
        return false;

    // logarithmic weights (as in Tristan's common.C)
    float x_weighted_sum = 0.0f;
    float y_weighted_sum = 0.0f;
    float total_weight = 0.0f;
    const float w = 4.0f;

    for (int i = 0; i < MAX_NUM_CRYSTALS; ++i)
    {
        if (is_skipped_crystal(i))
            continue;
        const int x = i % 5;
        const int y = i / 5;
        const float sig = values[i];
        if (sig <= 0)
            continue;

        float weight = w + std::log(sig / total_signal);
        if (weight < 0)
            weight = 0;

        total_weight += weight;
        x_weighted_sum += x * weight;
        y_weighted_sum += y * weight;
    }

    if (total_weight <= 0)
        return false;

    const float x_cog = x_weighted_sum / total_weight;
    const float y_cog = y_weighted_sum / total_weight;

    if (distribution)
        distribution->Fill(x_cog, y_cog);

    return ellipse_cut(x_cog, y_cog);
}

static float load_adc_per_refunit(const char *adc_calib_root, float ref_voltage)
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

    // Try exact key first
    std::string exact = std::string(Form("mean_adc_to_ref_calibration_%.2fV", ref_voltage));
    TParameter<float> *p = (TParameter<float> *)f->Get(exact.c_str());

    // Fallback: find any key starting with the prefix
    if (!p)
    {
        TIter nextkey(f->GetListOfKeys());
        TKey *key;
        while ((key = (TKey *)nextkey()))
        {
            std::string name = key->GetName();
            if (name.rfind("mean_adc_to_ref_calibration_", 0) == 0)
            {
                p = (TParameter<float> *)f->Get(name.c_str());
                if (p)
                    break;
            }
        }
    }

    float val = 0.0f;
    if (p)
        val = p->GetVal();

    f->Close();
    delete f;
    return val;
}

struct TotCalibParams
{
    bool ok = false;
    float c0 = 0, c1 = 0, c2 = 0; // pol2
    float a0 = 0, a1 = 0;         // pol1
    float adc_per_refunit = 0;    // stored by tot_calibration
    float ref_voltage_V = 1.27f;  // stored by tot_calibration
};

static TotCalibParams load_tot_calibration(const char *tot_calib_root)
{
    TotCalibParams p;
    if (!tot_calib_root)
        return p;

    TFile *f = TFile::Open(tot_calib_root, "READ");
    if (!f || f->IsZombie())
    {
        if (f)
        {
            f->Close();
            delete f;
        }
        return p;
    }

    auto getPar = [&](const char *name) -> TParameter<float> *
    {
        return (TParameter<float> *)f->Get(name);
    };

    TParameter<float> *c0 = getPar("tot_c0");
    TParameter<float> *c1 = getPar("tot_c1");
    TParameter<float> *c2 = getPar("tot_c2");
    TParameter<float> *a0 = getPar("tot_a0");
    TParameter<float> *a1 = getPar("tot_a1");

    if (c0 && c1 && c2 && a0 && a1)
    {
        p.c0 = c0->GetVal();
        p.c1 = c1->GetVal();
        p.c2 = c2->GetVal();
        p.a0 = a0->GetVal();
        p.a1 = a1->GetVal();
        p.ok = true;
    }

    // optional bookkeeping
    if (auto *ap = (TParameter<float> *)f->Get("adc_per_refunit"))
        p.adc_per_refunit = ap->GetVal();
    if (auto *rv = (TParameter<float> *)f->Get("adc_ref_voltage_V"))
        p.ref_voltage_V = rv->GetVal();

    f->Close();
    delete f;
    return p;
}

// -------------------- optional: use per-run ToT resolution from tot_calibration_values.root --------------------
static inline float tot_to_voltage_equiv(const TotCalibParams &p, float tot_sum_norm);

// tot_calibration.cxx writes (per run):
//   run%03d_tot_sum_mu, run%03d_tot_sum_sigma, run%03d_tot_sum_res_pct (and also run%03d_tot_sum histogram)
// Here we can reuse those for "ToT-only" runs where ADC is saturated.

static bool get_param_as_double(TFile *f, const char *name, double &out)
{
    out = 0.0;
    if (!f || !name)
        return false;

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
}

struct RunTotSumStats
{
    bool ok = false;
    double mu = 0.0;
    double sigma = 0.0;
    double emu = 0.0;
    double esigma = 0.0;
};

static RunTotSumStats load_tot_sum_stats(TFile *f, int run)
{
    RunTotSumStats s;
    if (!f || f->IsZombie())
        return s;

    // Try both padded and unpadded run keys, and both float/double TParameter payloads.
    double mu_val = 0.0, si_val = 0.0;

    bool ok_mu = get_param_as_double(f, Form("run%03d_tot_sum_mu", run), mu_val) ||
                 get_param_as_double(f, Form("run%d_tot_sum_mu", run), mu_val);
    bool ok_si = get_param_as_double(f, Form("run%03d_tot_sum_sigma", run), si_val) ||
                 get_param_as_double(f, Form("run%d_tot_sum_sigma", run), si_val);

    if (!ok_mu || !ok_si)
        return s;

    s.mu = mu_val;
    s.sigma = si_val;
    s.ok = (s.mu > 0.0 && s.sigma > 0.0);

    // try to recover fit errors if the histogram + fitted function were stored
    TH1 *h = dynamic_cast<TH1 *>(f->Get(Form("run%03d_tot_sum", run)));
    if (!h)
        h = dynamic_cast<TH1 *>(f->Get(Form("run%d_tot_sum", run)));

    if (h)
    {
        TF1 *fn = h->GetFunction(Form("tot_sum_gaus_run%03d", run));
        if (!fn)
            fn = h->GetFunction(Form("tot_sum_gaus_run%d", run));

        if (fn)
        {
            s.emu = fn->GetParError(1);
            s.esigma = fn->GetParError(2);
        }
    }
    return s;
}

static inline float tot_to_voltage_equiv(const TotCalibParams &p, float tot_sum_norm)
{
    // Tot calibration is: y = label_minus_missing (treated as volts) as a function of x=ToT sum
    const float y_poly = p.c0 + p.c1 * tot_sum_norm + p.c2 * tot_sum_norm * tot_sum_norm;
    const float y_lin = p.a0 + p.a1 * tot_sum_norm;
    float y = std::max(y_poly, y_lin);
    if (y < 0)
        y = 0;
    return y;
}

// Convert ToT-sum (mu,sigma) -> V-equivalent (mu,sigma) using the calibrated mapping V = f(ToT).
// Use a symmetric finite-difference to handle the piecewise max(pol2,pol1).
static void tot_stats_to_v_equiv(const TotCalibParams &p,
                                 double mu_tot,
                                 double sigma_tot,
                                 double &mu_v,
                                 double &sigma_v)
{
    mu_v = 0.0;
    sigma_v = 0.0;
    if (!p.ok || !(mu_tot > 0.0) || !(sigma_tot > 0.0))
        return;

    auto f = [&](double x) -> double
    {
        if (x < 0)
            x = 0;
        return (double)tot_to_voltage_equiv(p, (float)x);
    };

    mu_v = f(mu_tot);

    const double vp = f(mu_tot + sigma_tot);
    const double vm = f(mu_tot - sigma_tot);
    sigma_v = 0.5 * std::abs(vp - vm);

    // Fallback: if the V(ToT) mapping is locally flat (e.g. piecewise max picks a plateau),
    // preserve the relative width from ToT space so we still get a meaningful resolution point.
    if (!(sigma_v > 0.0) && mu_v > 0.0 && mu_tot > 0.0)
    {
        const double rel = sigma_tot / mu_tot;
        if (rel > 0.0)
            sigma_v = std::abs(mu_v) * rel;
    }

    if (!(mu_v > 0.0) || !(sigma_v > 0.0))
    {
        mu_v = 0.0;
        sigma_v = 0.0;
    }
}
static void load_gain_factors(const char *gain_root,
                              TH1 *&gain_factors,
                              TH1 *&crystal_factor,
                              TFile *&gain_file_handle)
{
    gain_factors = nullptr;
    crystal_factor = nullptr;
    gain_file_handle = nullptr;

    if (gain_root)
    {
        gain_file_handle = TFile::Open(gain_root, "READ");
        if (gain_file_handle && !gain_file_handle->IsZombie())
        {
            gain_factors = (TH1 *)gain_file_handle->Get("gain_factors");
            crystal_factor = (TH1 *)gain_file_handle->Get("crystal_factor");
        }
    }

    // fallbacks
    if (!gain_factors)
    {
        gain_factors = new TH1F("gain_factors_unity", "Gain Factors;crystal*16+sipm;Gain", MAX_NUM_CRYSTALS * SIPMS_PER_CRYSTAL, 0, MAX_NUM_CRYSTALS * SIPMS_PER_CRYSTAL);
        for (int i = 1; i <= gain_factors->GetNbinsX(); ++i)
            gain_factors->SetBinContent(i, 1.0);
    }

    if (!crystal_factor)
    {
        crystal_factor = new TH1F("crystal_factor_unity", "Crystal Factors;crystal;Gain", MAX_NUM_CRYSTALS, 0, MAX_NUM_CRYSTALS);
        for (int i = 1; i <= crystal_factor->GetNbinsX(); ++i)
            crystal_factor->SetBinContent(i, 1.0);
    }
}

static inline float get_channel_gain(TH1 *gain_factors, int crystal, int sipm)
{
    if (!gain_factors)
        return 1.0f;
    const int idx = crystal * SIPMS_PER_CRYSTAL + sipm;
    if (idx < 0)
        return 1.0f;
    return (float)gain_factors->GetBinContent(idx + 1);
}

static inline float get_crystal_gain(TH1 *crystal_factor, int crystal)
{
    if (!crystal_factor)
        return 1.0f;
    return (float)crystal_factor->GetBinContent(crystal + 1);
}

// -------------------- main --------------------

struct RunPoint
{
    int run;
    float voltage;
};

static std::vector<RunPoint> default_led_runs()
{
    // same list as led_analysis.C (Run, voltage)
    return {
        // {23, 0.0f},
        // {26, 1.20f},
        // {30, 1.22f},
        // {33, 1.24f},
        {36, 1.25f},
        {39, 1.26f},
        {42, 1.27f},
        {45, 1.28f},
        {48, 1.29f},
        {51, 1.30f},
        {54, 1.32f},
        {57, 1.33f},
        {60, 1.34f} //,
        // {63, 1.36f},
        // {66, 1.37f},
        // {69, 1.38f},
        // {72, 1.40f},
        // {75, 1.42f},
        // {78, 1.44f},
        // {81, 1.46f},
        // {84, 1.48f},
        // {87, 1.50f},
        // {90, 1.52f},
        // {93, 1.54f},
        // {96, 1.56f},
        // {99, 1.58f},
        // {102, 1.60f},
        // {105, 1.62f},
        // {108, 1.64f},
        // {111, 1.66f},
        // {114, 1.68f},
        // {117, 1.70f},
        // {120, 1.72f},
        // {123, 1.74f},
        // {126, 1.76f},
        // {129, 1.78f},
        // {132, 1.80f},
        // {135, 1.82f},
        // {138, 1.84f},
        // {141, 1.86f},
        // {144, 1.88f}
    };
}

void energy_resolution_led_scan(const char *data_dir = "data",
                                const char *mapping_csv = "eeemcal_desy_dec2025_mapping_v2.csv",
                                const char *outdir = "outputs",
                                const char *adc_calib_root = "outputs/adc_to_ref_calibration_1.27V.root",
                                const char *tot_calib_root = "outputs/tot_calibration_values.root",
                                int central_crystal = 12)
{
    TH1::AddDirectory(kFALSE);
    gROOT->cd();

    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
    gErrorIgnoreLevel = kWarning;

    // signal extraction settings (shared with other macros)
    g_signal_method = 3; // 2,3,4,5,7
    g_tot_min = 50;      // ToT threshold in common_led

    mkdir_p(outdir);

    // ---- mapping
    auto mapping = read_mapping_csv(mapping_csv, SIPMS_PER_CRYSTAL);
    if (mapping.empty())
    {
        std::cerr << "Error: mapping empty. Check mapping CSV: " << mapping_csv << "\n";
        return;
    }

    // reverse[channel] = {crystal,sipm}
    std::unordered_map<int, std::pair<int, int>> reverse;
    reverse.reserve(mapping.size() * SIPMS_PER_CRYSTAL);
    for (const auto &kv : mapping)
    {
        const int cr = kv.first;
        auto chans = get_crystal_channels(mapping, cr, SIPMS_PER_CRYSTAL);
        for (int sipm = 0; sipm < SIPMS_PER_CRYSTAL; ++sipm)
        {
            const int ch = chans[sipm];
            if (ch < 0)
                continue;
            reverse[ch] = {cr, sipm};
        }
    }

    // be generous with max_channels in case mapping uses global channels beyond 0-127
    auto active_channels = get_active_channels_from_mapping(mapping, SIPMS_PER_CRYSTAL, MAX_NUM_CRYSTALS, 576);

    // ---- ToT calibration
    TotCalibParams totP = load_tot_calibration(tot_calib_root);
    if (g_use_tot_for_central && !totP.ok)
    {
        std::cerr << "Warning: could not load ToT calibration from " << tot_calib_root << ". Will fall back to ADC-only.\n";
        g_use_tot_for_central = false;
    }

    // Determine reference voltage and adc_per_refunit
    float refV = totP.ref_voltage_V;
    if (!(refV > 0))
        refV = 1.27f;

    float adc_per_refunit = totP.adc_per_refunit;
    if (!(adc_per_refunit > 0))
    {
        adc_per_refunit = load_adc_per_refunit(adc_calib_root, refV);
    }
    if (!(adc_per_refunit > 0))
    {
        std::cerr << "Error: adc_per_refunit not available. Provide a valid adc_calib_root (from adc_calibration.cxx) or a tot_calibration_values.root that includes it.\n";
        return;
    }

    std::cout << "Using adc_per_refunit=" << adc_per_refunit << " ADC per " << refV << "V\n";
    if (g_use_tot_for_central)
    {
        std::cout << "Using ToT calibration: pol2(" << totP.c0 << "," << totP.c1 << "," << totP.c2 << ") and pol1(" << totP.a0 << "," << totP.a1 << ")\n";
    }

    // Optional: open tot_calibration_values.root to reuse per-run ToT-sum resolution for ToT-only runs
    TFile *f_tot_values = nullptr;
    if (g_use_totcalib_resolution_for_totonly)
    {
        f_tot_values = TFile::Open(tot_calib_root, "READ");
        if (!f_tot_values || f_tot_values->IsZombie())
        {
            std::cerr << "Warning: could not open " << tot_calib_root
                      << " to reuse per-run ToT resolution. Will ignore override.\n";
            if (f_tot_values)
            {
                f_tot_values->Close();
                delete f_tot_values;
                f_tot_values = nullptr;
            }
        }
    }

    // ---- run list
    auto runs = default_led_runs();

    float maxV = 0.0f;
    for (auto &rp : runs)
        maxV = std::max(maxV, rp.voltage);

    // histogram x-range in "voltage-equivalent" units
    const float xMax = std::max(3.0f, 1.6f * maxV);

    // graphs
    TGraphErrors *g_res_c1 = new TGraphErrors();
    g_res_c1->SetName("res_central1");
    g_res_c1->SetTitle("Central crystal resolution;LED voltage label (V);Resolution (%)");

    TGraphErrors *g_res_c9 = new TGraphErrors();
    g_res_c9->SetName("res_central9");
    g_res_c9->SetTitle("Central 3x3 resolution;LED voltage label (V);Resolution (%)");

    TGraphErrors *g_res_tot = new TGraphErrors();
    g_res_tot->SetName("res_total");
    g_res_tot->SetTitle("Total (5x5) resolution;LED voltage label (V);Resolution (%)");

    // mean graphs (sanity)
    TGraphErrors *g_mu_tot = new TGraphErrors();
    g_mu_tot->SetName("mu_total");
    g_mu_tot->SetTitle("Total reconstructed mean;LED voltage label (V);mean (V-equivalent)");

    TGraphErrors *g_mu_c1 = new TGraphErrors();
    g_mu_c1->SetName("mu_central1");
    g_mu_c1->SetTitle("Central reconstructed mean;LED voltage label (V);mean (V-equivalent)");

    TCanvas *canvas = new TCanvas("c_energy_resolution", "", 1100, 700);
    std::string pdf = std::string(Form("%s/energy_resolution_led.pdf", outdir));
    canvas->SaveAs((pdf + "[").c_str());

    // store per-run histograms for the ROOT output
    std::vector<TH1F *> h_c1, h_c9, h_tot, h_c9_rest, h_tot_rest;
    std::vector<TH2F *> h_cog;

    int ipt = 0;
    for (const auto &rp : runs)
    {
        const int run = rp.run;
        const float Vlabel = rp.voltage;

        std::string infile = std::string(Form("%s/Run%03d.root", data_dir, run));

        // gain factors for this voltage
        std::string gain_root = std::string(Form("%s/gain_match_%.2fV.root", outdir, Vlabel));
        TH1 *gain_factors = nullptr;
        TH1 *crystal_factor = nullptr;
        TFile *gain_handle = nullptr;
        load_gain_factors(gain_root.c_str(), gain_factors, crystal_factor, gain_handle);

        // open data
        TFile *f = TFile::Open(infile.c_str());
        if (!f || f->IsZombie())
        {
            std::cerr << "Warning: could not open " << infile << " (skipping)\n";
            if (f)
            {
                f->Close();
                delete f;
            }
            // cleanup gains
            if (gain_handle)
            {
                gain_handle->Close();
                delete gain_handle;
            }
            else
            {
                // created fallbacks
                if (gain_factors && std::string(gain_factors->GetName()).find("_unity") != std::string::npos)
                    delete gain_factors;
                if (crystal_factor && std::string(crystal_factor->GetName()).find("_unity") != std::string::npos)
                    delete crystal_factor;
            }
            continue;
        }

        TTree *t = (TTree *)f->Get("events");
        if (!t)
        {
            std::cerr << "Warning: no TTree 'events' in " << infile << " (skipping)\n";
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
            continue;
        }

        // determine adc/tot shapes
        TBranch *br_adc = t->GetBranch("adc");
        if (!br_adc)
        {
            std::cerr << "Warning: missing branch 'adc' in " << infile << " (skipping)\n";
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
            continue;
        }

        TLeaf *leaf_adc = br_adc->GetLeaf("adc");
        if (!leaf_adc)
        {
            std::cerr << "Warning: missing leaf 'adc' in " << infile << " (skipping)\n";
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
            continue;
        }

        const int n_adc = leaf_adc->GetLen();
        if ((n_adc % SAMPLES_PER_CHANNEL) != 0)
        {
            std::cerr << "Warning: unexpected adc length=" << n_adc << " in " << infile << " (skipping)\n";
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
            continue;
        }

        const int n_channels = n_adc / SAMPLES_PER_CHANNEL;

        // tot branch is optional
        bool have_tot = false;
        TBranch *br_tot = t->GetBranch("tot");
        TLeaf *leaf_tot = nullptr;
        int n_tot = 0;
        if (br_tot)
        {
            leaf_tot = br_tot->GetLeaf("tot");
            if (leaf_tot)
            {
                n_tot = leaf_tot->GetLen();
                if (n_tot == n_adc)
                    have_tot = true;
            }
        }

        std::vector<uint32_t> adc_buf((size_t)n_adc);
        std::vector<uint32_t> tot_buf((size_t)n_tot);
        t->SetBranchAddress("adc", adc_buf.data());
        if (have_tot)
            t->SetBranchAddress("tot", tot_buf.data());

        // histograms for this run
        TH1F *hc1 = new TH1F(Form("run%03d_central1", run), Form("Run %03d | Central crystal (V-equiv);V-equiv;Events", run), 500, 0, xMax);
        TH1F *hc9 = new TH1F(Form("run%03d_central9", run), Form("Run %03d | Central 3x3 (V-equiv);V-equiv;Events", run), 500, 0, xMax);
        TH1F *htot = new TH1F(Form("run%03d_total", run), Form("Run %03d | Total 5x5 (V-equiv);V-equiv;Events", run), 500, 0, xMax);
        TH2F *hcog = new TH2F(Form("run%03d_cog", run), Form("Run %03d | COG;X (crystal index);Y (crystal index)", run), 100, -0.5, 4.5, 100, -0.5, 4.5);

        // For ToT-only runs: keep the "rest-of-sum" (excluding the center) to combine later with ToT resolution
        TH1F *hc9_rest = new TH1F(Form("run%03d_central9_rest", run),
                                  Form("Run %03d | 3x3 (excluding center) (V-equiv);V-equiv;Events", run),
                                  500, 0, xMax);
        TH1F *htot_rest = new TH1F(Form("run%03d_total_rest", run),
                                   Form("Run %03d | Total (excluding center) (V-equiv);V-equiv;Events", run),
                                   500, 0, xMax);

        Long64_t n_central_tot_only = 0; // number of used events where central has ToT on all 16 and no ADC contribution

        // central 3x3 indices around 12
        const int c = central_crystal;
        const std::array<int, 9> idx9 = {c - 6, c - 5, c - 4,
                                         c - 1, c, c + 1,
                                         c + 4, c + 5, c + 6};

        const Long64_t nentries = t->GetEntries();
        const Long64_t nloop = (g_max_events > 0) ? std::min((Long64_t)g_max_events, nentries) : nentries;

        Long64_t skipped_noncentral_tot = 0;
        Long64_t used_events = 0;
        Long64_t skipped_cog = 0;

        for (Long64_t ev = 0; ev < nloop; ++ev)
        {
            t->GetEntry(ev);

            std::array<float, MAX_NUM_CRYSTALS> adc_crystal;
            adc_crystal.fill(0.0f);

            // central ToT bookkeeping
            float central_tot_sum = 0.0f;
            int central_tot_used = 0;

            bool reject = false;

            for (int ch : active_channels)
            {
                if (ch < 0 || ch >= n_channels)
                    continue;

                auto it = reverse.find(ch);
                if (it == reverse.end())
                    continue;

                const int cr = it->second.first;
                const int sipm = it->second.second;
                if (cr < 0 || cr >= MAX_NUM_CRYSTALS || sipm < 0 || sipm >= SIPMS_PER_CRYSTAL)
                    continue;
                if (!mapping.count(cr))
                    continue;
                if (is_skipped_crystal(cr))
                    continue;

                uint32_t *adc_ptr = &adc_buf[(size_t)ch * SAMPLES_PER_CHANNEL];
                uint32_t *tot_ptr = (have_tot ? &tot_buf[(size_t)ch * SAMPLES_PER_CHANNEL] : nullptr);

                const float gch = get_channel_gain(gain_factors, cr, sipm);

                // decide ToT usage
                const bool hasTOT = (have_tot && g_use_tot_for_central && tot_ptr && ::has_tot(tot_ptr));

                if (hasTOT)
                {
                    if (cr != central_crystal)
                    {
                        if (g_reject_noncentral_tot)
                        {
                            reject = true;
                            break;
                        }
                    }
                    else
                    {
                        // central ToT: use max ToT sample (more stable than first)
                        const float raw_tot = (float)get_tot_max(tot_ptr);
                        if (raw_tot > 0)
                        {
                            central_tot_sum += raw_tot * gch;
                            central_tot_used++;
                        }
                    }
                }
                else
                {
                    // ADC
                    const float sig_adc = calculate_signal_adc(adc_ptr, gch);
                    adc_crystal[cr] += sig_adc;
                }
            }

            if (reject)
            {
                skipped_noncentral_tot++;
                continue;
            }

            const bool evt_central_tot_only = (g_use_tot_for_central && (central_tot_used == SIPMS_PER_CRYSTAL) && (adc_crystal[central_crystal] == 0.0f));

            // apply per-crystal factor and convert to V-equivalent
            std::array<float, MAX_NUM_CRYSTALS> v_crystal;
            v_crystal.fill(0.0f);

            for (int cr = 0; cr < MAX_NUM_CRYSTALS; ++cr)
            {
                if (!mapping.count(cr))
                    continue;
                if (is_skipped_crystal(cr))
                    continue;
                const float gcr = get_crystal_gain(crystal_factor, cr);
                const float adc = adc_crystal[cr] * gcr;
                if (adc <= 0)
                    continue;
                v_crystal[cr] = (adc / adc_per_refunit) * refV;
            }

            // central voltage-equivalent includes ADC part + ToT part
            float v_center = v_crystal[central_crystal];
            if (g_use_tot_for_central && totP.ok && central_tot_used > 0)
            {
                const float frac = (float)central_tot_used / (float)SIPMS_PER_CRYSTAL;
                float tot_norm = central_tot_sum;
                if (frac > 0)
                    tot_norm = central_tot_sum / frac; // scale to 16 channels

                float v_norm = tot_to_voltage_equiv(totP, tot_norm);
                float v_tot = v_norm * frac; // scale back

                v_center += v_tot;
                v_crystal[central_crystal] = v_center;
            }

            // COG cut (optional)
            if (!calculate_cog(hcog, v_crystal))
            {
                skipped_cog++;
                continue;
            }

            // build sums
            float v_total = 0.0f;
            for (int cr = 0; cr < MAX_NUM_CRYSTALS; ++cr)
            {
                if (!mapping.count(cr))
                    continue;
                if (is_skipped_crystal(cr))
                    continue;
                v_total += v_crystal[cr];
            }

            float v_3x3 = 0.0f;
            for (int k = 0; k < 9; ++k)
            {
                const int cr = idx9[k];
                if (cr < 0 || cr >= MAX_NUM_CRYSTALS)
                    continue;
                if (!mapping.count(cr))
                    continue;
                if (is_skipped_crystal(cr))
                    continue;
                v_3x3 += v_crystal[cr];
            }

            // "rest" sums (excluding center) for optional ToT-only override
            hc9_rest->Fill(std::max(0.0f, v_3x3 - v_center));
            htot_rest->Fill(std::max(0.0f, v_total - v_center));

            hc1->Fill(v_center);
            hc9->Fill(v_3x3);
            htot->Fill(v_total);

            if (evt_central_tot_only)
                n_central_tot_only++;

            used_events++;
        }

        std::cout << "Run " << run << " (" << Vlabel << " V): used=" << used_events
                  << "  skipped_noncentral_tot=" << skipped_noncentral_tot
                  << "  skipped_cog=" << skipped_cog
                  << "  central_tot_only=" << n_central_tot_only
                  << "  entries_total=" << htot->GetEntries()
                  << "\n";

        // ---- fit + per-run pages
        double mu1, si1, emu1, esi1;
        double mu9, si9, emu9, esi9;
        double mut, sit, emut, esit;

        const bool ok1 = fit_peak_robust(hc1, mu1, si1, emu1, esi1);
        const bool ok9 = fit_peak_robust(hc9, mu9, si9, emu9, esi9);
        const bool okt = fit_peak_robust(htot, mut, sit, emut, esit);

        const bool run_tot_only = (used_events > 0) && ((double)n_central_tot_only / (double)used_events > 0.80);

        const bool force_override = force_tot_override(Vlabel);

        // Values that will be used to fill the resolution graphs (can be overridden for ToT-only runs)
        bool ok1_use = ok1, ok9_use = ok9, okt_use = okt;
        double mu1_use = mu1, si1_use = si1, emu1_use = emu1, esi1_use = esi1;
        double mu9_use = mu9, si9_use = si9, emu9_use = emu9, esi9_use = esi9;
        double mut_use = mut, sit_use = sit, emut_use = emut, esit_use = esit;

        bool used_totcalib_override = false;

        if ((force_override || run_tot_only) && g_use_totcalib_resolution_for_totonly && f_tot_values && totP.ok)
        {
            RunTotSumStats ts = load_tot_sum_stats(f_tot_values, run);
            double muVc = 0.0, siVc = 0.0;
            if (ts.ok)
            {
                tot_stats_to_v_equiv(totP, ts.mu, ts.sigma, muVc, siVc);
                if (muVc > 0.0 && siVc > 0.0)
                {
                    used_totcalib_override = true;

                    ok1_use = true;
                    mu1_use = muVc;
                    si1_use = siVc;
                    emu1_use = 0.0;
                    esi1_use = 0.0;

                    // Combine ToT-based center with measured "rest" (excluding the center).
                    //
                    // Key point for saturated LED scans:
                    //   For ToT-only runs, the non-central contribution can be ~0 event-by-event.
                    //   Then hc9_rest/htot_rest are delta-like at 0, so:
                    //     mean(rest) can legitimately be 0
                    //     RMS(rest) can legitimately be 0
                    // If we require mean(rest)>0, we would DROP those points from the 3x3/5x5 graphs.
                    // So here we always use histogram moments and allow mean(rest)==0.

                    const double muR9 = (hc9_rest ? (double)hc9_rest->GetMean() : 0.0);
                    const double siR9 = (hc9_rest ? (double)hc9_rest->GetRMS() : 0.0);
                    const double muRt = (htot_rest ? (double)htot_rest->GetMean() : 0.0);
                    const double siRt = (htot_rest ? (double)htot_rest->GetRMS() : 0.0);

                    // 3x3 = center(ToT) + rest(ADC)
                    ok9_use = true;
                    mu9_use = muVc + muR9;
                    si9_use = std::sqrt(siVc * siVc + siR9 * siR9);
                    emu9_use = 0.0;
                    esi9_use = 0.0;

                    // 5x5 = center(ToT) + rest(ADC)
                    okt_use = true;
                    mut_use = muVc + muRt;
                    sit_use = std::sqrt(siVc * siVc + siRt * siRt);
                    emut_use = 0.0;
                    esit_use = 0.0;
                }
            }

            if (!ts.ok)
                std::cerr << "WARNING: forced override requested, but ToT stats missing for run "
                          << run << " (V=" << Vlabel << ")\n";
        }

        std::cout << "Run " << run
                  << " V=" << Vlabel
                  << " run_tot_only=" << run_tot_only
                  << " used_totcalib_override=" << used_totcalib_override
                  << " mu1_fit=" << mu1 << " si1_fit=" << si1
                  << " mu1_use=" << mu1_use << " si1_use=" << si1_use
                  << "\n";

        auto fill_graph = [&](TGraphErrors *g, double x, bool ok, double mu, double si, double emu, double esi)
        {
            if (!ok || mu <= 0 || si <= 0)
                return;
            const double res = 100.0 * (si / mu);
            const double rel_mu = (emu > 0) ? (emu / mu) : 0.0;
            const double rel_si = (esi > 0) ? (esi / si) : 0.0;
            const double eres = res * std::sqrt(rel_mu * rel_mu + rel_si * rel_si);
            const int n = g->GetN();
            g->SetPoint(n, x, res);
            g->SetPointError(n, 0.0, eres);
        };

        auto res_pct = [&](double mu, double si)
        {
            return 100.0 * (si / mu);
        };

        std::cout << "Run " << run << " V=" << Vlabel
                  << " -> res_c1_used=" << res_pct(mu1_use, si1_use)
                  << "% (mu=" << mu1_use << ", si=" << si1_use << ")\n";

        fill_graph(g_res_c1, Vlabel, ok1_use, mu1_use, si1_use, emu1_use, esi1_use);
        fill_graph(g_res_c9, Vlabel, ok9_use, mu9_use, si9_use, emu9_use, esi9_use);
        fill_graph(g_res_tot, Vlabel, okt_use, mut_use, sit_use, emut_use, esit_use);

        if (ok1_use && mu1_use > 0)
        {
            const int n = g_mu_c1->GetN();
            g_mu_c1->SetPoint(n, Vlabel, mu1_use);
            g_mu_c1->SetPointError(n, 0.0, emu1_use);
        }

        if (okt_use && mut_use > 0)
        {
            const int n = g_mu_tot->GetN();
            g_mu_tot->SetPoint(n, Vlabel, mut_use);
            g_mu_tot->SetPointError(n, 0.0, emut_use);
        }

        // plot pages
        TLatex text;
        text.SetNDC();
        text.SetTextFont(42);
        text.SetTextSize(0.04);

        canvas->Clear();
        hc1->Draw();
        text.DrawLatex(0.15, 0.86, Form("Run %03d  (%.2f V)", run, Vlabel));
        if (ok1_use)
            text.DrawLatex(0.15, 0.80, Form("Central: mean=%.3f  #sigma=%.3f  res=%.2f%%", mu1_use, si1_use, 100.0 * si1_use / mu1_use));
        if (used_totcalib_override)
            text.DrawLatex(0.15, 0.74, "Central width from tot_calibration_values.root");
        canvas->SaveAs(pdf.c_str());

        canvas->Clear();
        hc9->Draw();
        text.DrawLatex(0.15, 0.86, Form("Run %03d  (%.2f V)", run, Vlabel));
        if (ok9_use)
            text.DrawLatex(0.15, 0.80, Form("3x3: mean=%.3f  #sigma=%.3f  res=%.2f%%", mu9_use, si9_use, 100.0 * si9_use / mu9_use));
        if (used_totcalib_override)
            text.DrawLatex(0.15, 0.74, "Center from tot_calib + rest in quadrature");
        canvas->SaveAs(pdf.c_str());

        canvas->Clear();
        htot->Draw();
        text.DrawLatex(0.15, 0.86, Form("Run %03d  (%.2f V)", run, Vlabel));
        if (okt_use)
            text.DrawLatex(0.15, 0.80, Form("Total: mean=%.3f  #sigma=%.3f  res=%.2f%%", mut_use, sit_use, 100.0 * sit_use / mut_use));
        if (used_totcalib_override)
            text.DrawLatex(0.15, 0.74, "Center from tot_calib + rest in quadrature");
        canvas->SaveAs(pdf.c_str());

        canvas->Clear();
        hcog->Draw("colz");
        if (g_do_cog_ellipse_cut)
        {
            TEllipse *e = new TEllipse(g_cog_cx, g_cog_cy, g_cog_sx, g_cog_sy);
            e->SetFillStyle(0);
            e->SetLineWidth(2);
            e->Draw("same");
        }
        text.DrawLatex(0.15, 0.86, Form("Run %03d  (%.2f V)  COG cut: %s", run, Vlabel, g_do_cog_ellipse_cut ? "ellipse" : "OFF"));
        canvas->SaveAs(pdf.c_str());

        // keep in vectors for output
        h_c1.push_back(hc1);
        h_c9.push_back(hc9);
        h_tot.push_back(htot);
        h_c9_rest.push_back(hc9_rest);
        h_tot_rest.push_back(htot_rest);
        h_cog.push_back(hcog);

        // cleanup input
        f->Close();
        delete f;

        if (gain_handle)
        {
            gain_handle->Close();
            delete gain_handle;
        }

        // Delete only the local unity fallbacks (objects loaded from file are owned by the file/directory)
        if (gain_factors && std::string(gain_factors->GetName()).find("_unity") != std::string::npos)
            delete gain_factors;
        if (crystal_factor && std::string(crystal_factor->GetName()).find("_unity") != std::string::npos)
            delete crystal_factor;

        ipt++;
    }

    // summary pages
    auto draw_res_summary = [&](TGraphErrors *g, const char *label)
    {
        canvas->Clear();
        g->SetMarkerStyle(20);
        g->Draw("AP");
        TLatex t;
        t.SetNDC();
        t.SetTextFont(42);
        t.SetTextSize(0.04);
        t.DrawLatex(0.15, 0.86, label);
        canvas->SaveAs(pdf.c_str());
    };

    draw_res_summary(g_res_c1, "Resolution vs LED voltage (central crystal)");
    draw_res_summary(g_res_c9, "Resolution vs LED voltage (central 3x3)");
    draw_res_summary(g_res_tot, "Resolution vs LED voltage (total 5x5)");

    canvas->Clear();
    g_mu_c1->SetMarkerStyle(20);
    g_mu_c1->Draw("AP");
    {
        TLatex t;
        t.SetNDC();
        t.SetTextFont(42);
        t.SetTextSize(0.04);
        t.DrawLatex(0.15, 0.86, "Central reconstructed mean vs LED voltage (sanity)");
    }
    canvas->SaveAs(pdf.c_str());

    canvas->Clear();
    g_mu_tot->SetMarkerStyle(20);
    g_mu_tot->Draw("AP");
    {
        TLatex t;
        t.SetNDC();
        t.SetTextFont(42);
        t.SetTextSize(0.04);
        t.DrawLatex(0.15, 0.86, "Total reconstructed mean vs LED voltage (sanity)");
    }
    canvas->SaveAs(pdf.c_str());

    canvas->SaveAs((pdf + "]").c_str());
    std::cout << "Saved PDF:  " << pdf << "\n";

    // ROOT output
    std::string root_out = std::string(Form("%s/energy_resolution_led.root", outdir));
    TFile *out = TFile::Open(root_out.c_str(), "RECREATE");

    // store settings
    TParameter<long> p_nev("max_events", g_max_events);
    TParameter<int> p_central("central_crystal", central_crystal);
    TParameter<float> p_adc("adc_per_refunit", adc_per_refunit);
    TParameter<float> p_refV("ref_voltage_V", refV);
    TParameter<int> p_use_tot("use_tot_for_central", (int)g_use_tot_for_central);
    TParameter<int> p_rej_nctot("reject_noncentral_tot", (int)g_reject_noncentral_tot);
    TParameter<int> p_cog("do_cog_ellipse_cut", (int)g_do_cog_ellipse_cut);
    TParameter<int> p_totcal_res("use_totcalib_resolution_for_totonly", (int)g_use_totcalib_resolution_for_totonly);

    p_nev.Write();
    p_central.Write();
    p_adc.Write();
    p_refV.Write();
    p_use_tot.Write();
    p_rej_nctot.Write();
    p_cog.Write();
    p_totcal_res.Write();

    if (totP.ok)
    {
        TParameter<float>("tot_c0", totP.c0).Write();
        TParameter<float>("tot_c1", totP.c1).Write();
        TParameter<float>("tot_c2", totP.c2).Write();
        TParameter<float>("tot_a0", totP.a0).Write();
        TParameter<float>("tot_a1", totP.a1).Write();
    }

    g_res_c1->Write();
    g_res_c9->Write();
    g_res_tot->Write();
    g_mu_c1->Write();
    g_mu_tot->Write();

    for (auto *h : h_c1)
        h->Write();
    for (auto *h : h_c9)
        h->Write();
    for (auto *h : h_tot)
        h->Write();
    for (auto *h : h_c9_rest)
        h->Write();
    for (auto *h : h_tot_rest)
        h->Write();
    for (auto *h : h_cog)
        h->Write();

    out->Close();
    delete out;

    std::cout << "Saved ROOT: " << root_out << "\n";

    // cleanup
    delete canvas;
    delete g_res_c1;
    delete g_res_c9;
    delete g_res_tot;
    delete g_mu_c1;
    delete g_mu_tot;

    for (auto *h : h_c1)
        delete h;
    for (auto *h : h_c9)
        delete h;
    for (auto *h : h_tot)
        delete h;
    for (auto *h : h_c9_rest)
        delete h;
    for (auto *h : h_tot_rest)
        delete h;
    for (auto *h : h_cog)
        delete h;

    if (f_tot_values)
    {
        f_tot_values->Close();
        delete f_tot_values;
    }
}

// Backward-compatible wrapper
void energy_resolution()
{
    energy_resolution_led_scan();
}
