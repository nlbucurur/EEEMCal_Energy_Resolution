// Original code: https://github.com/tlprotzman/eeemcal_desy_dec2025/blob/main/common.cxx

#include "common_led.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <map>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <utility>
#include <string>
#include <cctype>

#include <TH1.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TH2.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TError.h>

// ===== Globals =====

int g_signal_method = 3; // 2, 3, 4, 5, 7 (7 = waveform crystal ball fit)
int g_sipms_to_use = 16; // Number of SiPMs to use per crystal (= 16)
uint32_t g_tot_min = 50; // Minimum ToT value to consider valid

int extract_run_number(const char *filename)
{
    if (!filename)
        return -1;

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

void print_progress_25(int progress)
{
    std::cout << " [";
    for (int i = 0; i < LED_MAX_NUM_CRYSTALS; ++i)
        std::cout << (i < progress ? "*" : " ");
    std::cout << "]\r" << std::flush;
    if (progress >= LED_MAX_NUM_CRYSTALS)
        std::cout << "\n";
}

void build_xy_for_crystal(int crystal_id, float &x, float &y)
{
    const int crystal_mapping[LED_MAX_NUM_CRYSTALS] = {
        4, 9, 14, 19, 24,
        3, 8, 13, 18, 23,
        2, 7, 12, 17, 22,
        1, 6, 11, 16, 21,
        0, 5, 10, 15, 20};

    int found = -1;
    for (int i = 0; i < LED_MAX_NUM_CRYSTALS; ++i)
    {
        if (crystal_mapping[i] == crystal_id)
        {
            found = i;
            break;
        }
    }

    if (found < 0)
    {
        x = -999.0f;
        y = -999.0f;
        return;
    }

    const int row = found / 5;
    const int col = found % 5;
    x = (float)col;
    y = (float)row;
}

bool point_in_ellipse(float x, float y, float cx, float cy, float sx, float sy)
{
    if (sx <= 0 || sy <= 0)
        return true;
    const float dx = (x - cx) / sx;
    const float dy = (y - cy) / sy;
    return (dx * dx + dy * dy) <= 1.0f;
}

void subtract_common_mode(uint32_t out_adc[LED_SAMPLES_PER_CHANNEL],
                          uint32_t *in_adc,
                          float cm)
{
    for (int i = 0; i < LED_SAMPLES_PER_CHANNEL; ++i)
    {
        float v = (float)in_adc[i] - cm;
        if (v < 0.0f)
            v = 0.0f;
        out_adc[i] = (uint32_t)(v + 0.5f);
    }
}

std::unordered_map<int, std::pair<int, int>> build_reverse_channel_map(
    const std::map<int, std::vector<int>> &mapping,
    int sipms_to_use)
{
    std::unordered_map<int, std::pair<int, int>> reverse;
    reverse.reserve(mapping.size() * (size_t)sipms_to_use);

    for (const auto &kv : mapping)
    {
        const int crystal_id = kv.first;
        const auto chans = get_crystal_channels(mapping, crystal_id, sipms_to_use);

        for (int sipm = 0; sipm < sipms_to_use; ++sipm)
        {
            const int ch = chans[sipm];
            if (ch < 0)
                continue;
            reverse[ch] = {crystal_id, sipm};
        }
    }

    return reverse;
}

std::map<int, std::vector<int>> read_mapping_csv(const std::string &filename,
                                                 int expected_sipms_per_crystal)
{
    std::map<int, std::vector<int>> mapping;
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open mapping file " << filename << std::endl;
        return mapping;
    }

    std::string line;
    std::getline(file, line); // skip header

    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }

        std::istringstream iss(line);
        int fpga, asic, connector, crystal, sipm, channel;
        char comma;

        // Parse CSV: FPGA,ASIC,Connector,Crystal,SiPM,Channel
        if (iss >> fpga >> comma >> asic >> comma >> connector >> comma >> crystal >> comma >> sipm >> comma >> channel)
        {
            // Ensure the vector is large enough for all SiPMs (16 per crystal)
            if (sipm < 0)
                continue;

            auto &vec = mapping[crystal];
            if ((int)vec.size() < expected_sipms_per_crystal)
            {
                vec.resize(expected_sipms_per_crystal, -1);
            }
            if (sipm < expected_sipms_per_crystal)
            {
                vec[sipm] = channel;
            }
        }
    }
    return mapping;
}

std::vector<int> get_crystal_channels(const std::map<int, std::vector<int>> &mapping,
                                     int crystal,
                                     int expected_sipms_per_crystal)
{
    std::vector<int> out(expected_sipms_per_crystal, -1);
    auto it = mapping.find(crystal); 
    if (it == mapping.end()) return out;
    
    const auto& v = it->second;
    for (int i = 0; i < expected_sipms_per_crystal && i < (int)v.size(); ++i)
    {
        out[i] = v[i];
    }
    return out;
}

// ===== Signal extraction (ADC) =====
static inline float pedestal_2samples(uint32_t *adc_values)
{
    return (adc_values[0] + adc_values[1] /*+ adc_values[2]*/) / /*3.0f*/2.0f;
}

float calculate_signal_v2(uint32_t *adc_values, float gain)
{
    // Pedestal is the mean of the first two samples
    float ped = pedestal_2samples(adc_values);

    // Signal is the sum of the three greatest values
    float max1 = 0.0f, max2 = 0.0f, max3 = 0.0f;

    for (int i = 5; i < 10; ++i)
    {
        float sample = adc_values[i] - ped;
        if (sample <= 0) continue;

        if (sample > max1)
        {
            max3 = max2;
            max2 = max1;
            max1 = sample;
        } else if (sample > max2)
        {
            max3 = max2;
            max2 = sample;
        } else if (sample > max3)
        {
            max3 = sample;
        }
    }
    return (max1 + max2 + max3) * gain;
}

float calculate_signal_v3(uint32_t *adc_values, float gain)
{
    // Pedestal is the mean of the first two samples
    float ped = pedestal_2samples(adc_values);

    // Signal is the sum of all samples above pedestal
    float signal = 0.0f;
    for (int i = 3; i <= 5; ++i)
    {
        float sample = adc_values[i] - ped;
        if (sample > 0) signal += sample;
    }
    return signal * gain;
}

float calculate_signal_v4(uint32_t *adc_values, float gain)
{
    // Pedestal is the mean of the first two samples
    float ped = pedestal_2samples(adc_values);
    
    // Signal is sample 6 minus pedestal
    float signal = adc_values[6] - ped;
    if (signal < 0)
    {
        signal = 0;
    }
    return signal * gain;
}

float calculate_signal_v5(uint32_t *adc_values, float gain)
{
    // Pedestal is the mean of the first two samples
    float ped = pedestal_2samples(adc_values);

    // Find the sample with the maximum adc
    float max_sample = 0.0f;
    int idx = 0;

    for (int i = 3; i < LED_SAMPLES_PER_CHANNEL; ++i)
    {
        float sample = adc_values[i] - ped;
        if (sample > max_sample)
        {
            max_sample = sample;
            idx = i;
        }
    }
    if (idx < 5 || idx > 16) idx = 6;

    float signal = adc_values[idx - 1] + adc_values[idx] + adc_values[idx + 1] + adc_values[idx + 2] - (4 * ped);

    if (signal < 4 * 4)
    { // pedestal * samples
        signal = 0;
    }
    return signal * gain;
}

// ===== Crystal Ball function for fitting =====
static double crystal_ball_waveform(double *inputs, double *par)
{
    // Parameters
    // alpha: Where the gaussian transitions to the power law tail - fix?
    // n: The exponent of the power law tail - fix?
    // x_bar: The mean of the gaussian - free
    // sigma: The width of the gaussian - fix ?
    // N: The normalization of the gaussian - free
    // B baseline - fix? Maybe offset

    double x      = inputs[0];
    double alpha  = par[0];
    double n      = par[1];
    double x_bar  = par[2];
    double sigma  = par[3];
    double N      = par[4];
    double offset = par[5];

    double A = std::pow(n / std::fabs(alpha), n) * std::exp(-0.5 * alpha * alpha);
    double B = n / std::fabs(alpha) - std::fabs(alpha);
    // std::cout << "A: " << A << std::endl;

    double t = (x - x_bar) / sigma;
    double y = 0.0;

    // std::cout << "alpha: " << alpha << " n: " << n << " x_bar: " << x_bar << " sigma: " << sigma << " N: " << N << " B: " << B << " A: " << A << std::endl;

    if (t < alpha)
    {
        // std::cout << "path a" << std::endl;
        y = std::exp(-0.5 * t * t);
    }
    else
    {
        // std::cout << "path b" << std::endl;
        y = A * std::pow(B + t, -1 * n);
    }
    // std::cout << "x: " << x << " y: " << ret_val << std::endl;
    return N * y + offset;
}

float calculate_signal_v6(uint32_t *adc_values, uint32_t *common_mode, float gain)
{
    // Pedestal is the mean of the first two samples
    float ped = pedestal_2samples(adc_values);
    float common_mode_pedestal = pedestal_2samples(common_mode);

    // Sum samples 6, 7 and 8 after common mode subtraction
    // float sig_0 = adc_values[5] - pedestal;
    // float comm_0 = common_mode[5] - common_mode_pedestal;

    // float sig_a = adc_values[6] - pedestal;
    // float comm_a = common_mode[6] - common_mode_pedestal;

    // float sig_b = adc_values[7] - pedestal;
    // float comm_b = common_mode[7] - common_mode_pedestal;

    // float sig_c = adc_values[8] - pedestal;
    // float comm_c = common_mode[8] - common_mode_pedestal;

    // float sig_d = adc_values[9] - pedestal;
    // float comm_d = common_mode[9] - common_mode_pedestal;

    // float sig = sig_a + sig_b + sig_c + sig_0 + sig_d;
    // float cm = comm_a + comm_b + comm_c + comm_0 + comm_d;
    // if (cm > -12) {  // pedestal rms * samples
    //     cm = 0;
    // } else {
    //     // std::cout << "Common mode avg: " << cm / 3.0f << std::endl;
    // }
    // sig -= cm;
    float sig = 0.0f;
    for (int i = 5; i <= 8; i++)
    {
        float sample = adc_values[i] - ped;
        float comm_sample = common_mode[i] - common_mode_pedestal;
        float corrected_sample = sample - comm_sample;
        if (corrected_sample > 0)
        {
            sig += corrected_sample;
        }
    }

    // Signal is the sum of the three greatest values after common mode subtraction
    float signal = sig;
    if (signal < 0)
    {
        signal = 0;
    }

    return signal * gain;
}

float calculate_signal_v7(uint32_t *adc_values, float gain)
{
    // Fit samples in a histogram (bins correspond to time sample index)

    TH1F temp("temp_waveform", "temp_waveform", LED_SAMPLES_PER_CHANNEL, 0, LED_SAMPLES_PER_CHANNEL);
    for (int i = 3; i < LED_SAMPLES_PER_CHANNEL; ++i)
    {
        float sample = adc_values[i];
        temp.SetBinContent(i, sample); // ROOT histograms start at bin 1
    }
    float ped = pedestal_2samples(adc_values);

    TF1 fit("cb_wave", crystal_ball_waveform, 4, LED_SAMPLES_PER_CHANNEL, 6);
    fit.SetParameters(1.1, 0.4, 6.0, 0.45, 50.0, ped);
    fit.SetParLimits(0, 1.0, 1.2);   // alpha
    fit.SetParLimits(1, 0.2, 0.8);   // n
    fit.SetParLimits(2, 0.5, 8.5);   // x_bar
    fit.SetParLimits(3, 0.25, 0.65); // sigma
    fit.SetParLimits(4, 0, 2000);    // N
    fit.FixParameter(5, ped);        // offset fixed to pedestal

    // Quiet fit
    TFitResultPtr fit_result = temp.Fit(&fit, "RQ0S");
    
    // Signal proxy = N parameter
    float sig = (fit_result.Get() && fit_result->Status()==0) ? fit.GetParameter(4) : 0.0f;

    return sig * gain;
}

float calculate_signal_adc(uint32_t *adc_values, float gain)
{
    if (g_signal_method == 2)
    {
        return calculate_signal_v2(adc_values, gain);
    }
    else if (g_signal_method == 3)
    {
        return calculate_signal_v3(adc_values, gain);
    }
    else if (g_signal_method == 4)
    {
        return calculate_signal_v4(adc_values, gain);
    }
    else if (g_signal_method == 5)
    {
        return calculate_signal_v5(adc_values, gain);
    }
    else if (g_signal_method == 7)
    {
        return calculate_signal_v7(adc_values, gain);
    }
    else
    {
        // Default to v3
        return calculate_signal_v3(adc_values, gain);
    }
}

// ===== ToT helpers =====

bool has_tot(uint32_t *tot_values)
{
    for (int i = 0; i < LED_SAMPLES_PER_CHANNEL; ++i)
    {
        if (tot_values[i] > g_tot_min)
        {
            return true;
        }
    }
    return false;
}

uint32_t get_tot_first(uint32_t *tot_values)
{
    for (int i = 0; i < LED_SAMPLES_PER_CHANNEL; ++i)
    {
        if (tot_values[i] > g_tot_min)
        {
            return tot_values[i];
        }
    }
    return 0;
}

uint32_t get_tot_max(uint32_t *tot_values)
{
    uint32_t max_tot = 0;
    for (int i = 0; i < LED_SAMPLES_PER_CHANNEL; ++i)
    {
        if (tot_values[i] > g_tot_min && tot_values[i] > max_tot)
        {
            max_tot = tot_values[i];
        }
    }
    return max_tot;
}

// Hybrid (ToT preferred if present)
bool calculate_signal_hybrid(uint32_t *adc_values,
                      uint32_t *tot_values,
                      float gain_adc,
                      float &signal_out,
                      bool &used_tot)
{
    if (has_tot(tot_values))
    {
        signal_out = get_tot_first(tot_values) * gain_adc; // Assuming same gain for ToT, adjust if different
        used_tot = true;
        return true;
    }
    else
    {
        signal_out = calculate_signal_adc(adc_values, gain_adc);
        used_tot = false;
        return false;
    }
}

// ===== Peak fitting =====
bool fit_peak_gaus(TH1 *hist, float &peak, float &sigma)
{
    if (!hist || hist->GetEntries() < 10)
    {
        std::cout << "Histogram is empty or null!" << std::endl;
        return false;
    }

    float mean = hist->GetMean();
    float rms = hist->GetRMS();
    if (rms <= 0)
    {
        std::cout << "Histogram has non-positive RMS!" << std::endl;
        return false;
    }

    TF1 fit("gauss_fit", "gaus", mean - 1.5f * rms, mean + 1.5f * rms);

    TFitResultPtr result = hist->Fit(&fit, "RQ0S");

    if (!result.Get() || result->Status() != 0)
    {
        std::cout << "Initial fit failed!" << std::endl;
        return false;
    }
    
    peak = fit.GetParameter(1);
    sigma = fit.GetParameter(2);

    return true;
}

void draw_text_basic(TF1 *fit, const std::string &label1, const std::string &label2, const std::string &label3)
{
    if (!fit) return;

    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.04);
    text.SetTextFont(42);
    text.SetTextAlign(31);

    float peak = fit->GetParameter(1);
    float sigma = fit->GetParameter(2);
    float res   = (peak != 0) ? (sigma/peak)*100.0f : 0.0f;

    // text.DrawLatex(0.93, 0.85, Form("%.01f GeV Electrons", beam_energy));
    // text.DrawLatex(0.93, 0.80, Form("Run %d", run_number));
    text.DrawLatex(0.93, 0.85, label1.c_str());
    if (!label2.empty()) text.DrawLatex(0.93, 0.80, label2.c_str());
    if (!label3.empty()) text.DrawLatex(0.93, 0.75, label3.c_str());

    text.DrawLatex(0.93, 0.70, Form("Peak: %.03f", peak));
    text.DrawLatex(0.93, 0.65, Form("Sigma: %.03f", sigma));
    text.DrawLatex(0.93, 0.60, Form("Resolution: %.2f%%", res));
    text.DrawLatex(0.93, 0.55, Form("Method %d", g_signal_method));
}

std::vector<int> get_active_channels_from_mapping(const std::map<int, std::vector<int>>& mapping,
                                                  int sipms_to_use,
                                                  int n_crystals,
                                                  int max_channels)
{
    std::set<int> unique;
    for (int cr = 0; cr < n_crystals; ++cr) {
        auto it = mapping.find(cr);
        if (it == mapping.end()) continue;
        const auto& v = it->second;
        for (int s = 0; s < sipms_to_use && s < (int)v.size(); ++s) {
            int ch = v[s];
            if (ch >= 0 && ch < max_channels) unique.insert(ch);
        }
    }
    return std::vector<int>(unique.begin(), unique.end());
}

// float sigma_cut = 2;

// int sipms_to_use = 16;

// float center_x = 1.95396;
// float sigma_x = 0.184758 * sigma_cut;
// float center_y = 1.9688;
// float sigma_y = 0.195137 * sigma_cut;

// // float adc_calib = 26704.4;  // 32444.1 Signal_ADC = 1 GeV
// // float adc_calib = 57854.3;  // v3

// int signal_method = 3;

// bool calculate_signal(uint32_t *adc_values, uint32_t *tot_values, float gain, float &signal)
// {
//     if (signal_method == 2)
//     {
//         // Check if there is a ToT value
//         uint32_t tot = 0;
//         for (int i = 0; i < 20; ++i)
//         {
//             if (tot_values[i] > 0)
//             {
//                 tot = tot_values[i];
//                 break;
//             }
//         }
//         // If there is no ToT value
//         if (tot == 0)
//         {
//             signal = calculate_signal(adc_values, gain);
//             return false;
//         }

//         // Otherwise, just return the first ToT value
//         signal = tot * gain;
//         return true;
//     }
//     else if (signal_method == 3)
//     {
//         // Check if there is a ToT value
//         uint32_t tot = 0;
//         for (int i = 0; i < 20; ++i)
//         {
//             if (tot_values[i] > 0)
//             {
//                 tot = tot_values[i];
//                 break;
//             }
//         }
//         // If there is no ToT value
//         if (tot == 0)
//         {
//             signal = calculate_signal_v3(adc_values, gain);
//             return false;
//         }
//         signal = tot * gain;
//         return true;
//     }
//     else if (signal_method == 4)
//     {
//         // Check if there is a ToT value
//         uint32_t tot = 0;
//         for (int i = 0; i < 20; ++i)
//         {
//             if (tot_values[i] > 0)
//             {
//                 tot = tot_values[i];
//                 break;
//             }
//         }
//         // If there is no ToT value
//         if (tot == 0)
//         {
//             signal = calculate_signal_v4(adc_values, gain);
//             return false;
//         }
//         signal = tot * gain;
//         return true;
//     }
//     else if (signal_method == 5)
//     {
//         // Check if there is a ToT value
//         uint32_t tot = 0;
//         for (int i = 0; i < 20; ++i)
//         {
//             if (tot_values[i] > 0)
//             {
//                 tot = tot_values[i];
//                 break;
//             }
//         }
//         // If there is no ToT value
//         if (tot == 0)
//         {
//             signal = calculate_signal_v5(adc_values, gain);
//             return false;
//         }
//         signal = tot * gain;
//         return true;
//     }
//     else
//     {
//         // Default to ADC signal
//         signal = calculate_signal_adc(adc_values, gain);
//         return false;
//     }
// }

// bool is_tot(uint32_t *tot_values)
// {
//     for (int i = 0; i < 20; ++i)
//     {
//         // std::cout << "TOT Value[" << i << "]: " << tot_values[i] << std::endl;
//         if (tot_values[i] > g_tot_min)
//         {
//             return true;
//         }
//     }
//     return false;
// }

// int32_t get_toa(uint32_t *toa_values)
// {
//     for (int i = 0; i < 20; ++i)
//     {
//         if (toa_values[i] > 0)
//         {
//             return toa_values[i];
//         }
//     }
//     return -1;
// }




// bool position_cut(float x, float y)
// {
//     // Cut anything outside of the center
//     if (std::abs(x - center_x) > sigma_x || std::abs(y - center_y) > sigma_y)
//     {
//         return false;
//     }
//     return true;
// }

// bool calculate_cog(TH2 *distribution, float *values)
// {
//     float total_signal = 0;
//     float x_weighted_sum = 0;
//     float y_weighted_sum = 0;
//     float x_cog = 0;
//     float y_cog = 0;
//     float w = 4;

//     for (int i = 0; i < 25; i++)
//     {
//         total_signal += values[i];
//     }

//     float total_weight = 0;
//     for (int i = 0; i < 25; i++)
//     {
//         int x = i % 5;
//         int y = i / 5;
//         float signal = values[i];
//         float weight = w + std::log(signal / total_signal); // Avoid log(0)
//         if (weight < 0)
//         {
//             weight = 0;
//         }
//         total_weight += weight;
//         x_weighted_sum += x * weight;
//         y_weighted_sum += y * weight;
//     }
//     if (total_weight > 0)
//     {
//         x_cog = x_weighted_sum / total_weight;
//         y_cog = y_weighted_sum / total_weight;
//         distribution->Fill(x_cog, y_cog);
//     }
//     return position_cut(x_cog, y_cog);
// }

// void print_progress(int progress)
// {
//     std::cout << " [";
//     for (int i = 0; i < 25; i++)
//     {
//         if (i < progress)
//         {
//             std::cout << "*";
//         }
//         else
//         {
//             std::cout << " ";
//         }
//     }
//     std::cout << "]\r" << std::flush;
// }