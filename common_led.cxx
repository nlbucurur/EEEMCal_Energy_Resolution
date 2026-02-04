// Original code: https://github.com/tlprotzman/eeemcal_desy_dec2025/blob/main/common.cxx

#include <common_led.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <map>
#include <sstream>
#include <vector>

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
            if ((int)vec.size() <= expected_sipms_per_crystal)
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
static inline float pedestal_3samples(uint32_t *adc_values)
{
    return (adc_values[0] + adc_values[1] + adc_values[2]) / 3.0f;
}

float calculate_signal_v2(uint32_t *adc_values, float gain)
{
    // Pedestal is the mean of the first three samples
    float ped = pedestal_3samples(adc_values);

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
    // Pedestal is the mean of the first three samples
    float ped = pedestal_3samples(adc_values);

    // Signal is the sum of all samples above pedestal
    float signal = 0.0f;
    for (int i = 6; i <= 8; ++i)
    {
        float sample = adc_values[i] - ped;
        if (sample > 0) signal += sample;
    }
    return signal * gain;
}

float calculate_signal_v4(uint32_t *adc_values, float gain)
{
    // Pedestal is the mean of the first three samples
    float ped = pedestal_3samples(adc_values);
    
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
    // Pedestal is the mean of the first three samples
    float ped = pedestal_3samples(adc_values);

    // Find the sample with the maximum adc
    float max_sample = 0.0f;
    int idx = 0;

    for (int i = 3; i < 20; ++i)
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
static double crystal_ball(double *inputs, double *par)
{
    // Parameters
    // alpha: Where the gaussian transitions to the power law tail - fix?
    // n: The exponent of the power law tail - fix?
    // x_bar: The mean of the gaussian - free
    // sigma: The width of the gaussian - fix ?
    // N: The normalization of the gaussian - free
    // B baseline - fix?

    double x      = inputs[0];
    double alpha  = par[0];
    double n      = par[1];
    double x_bar  = par[2];
    double sigma  = par[3];
    double N      = par[4];
    double offset = par[5];

    double A = std::pow(n / fabs(alpha), n) * exp(-0.5 * alpha * alpha);
    double B = n / fabs(alpha) - fabs(alpha);
    // std::cout << "A: " << A << std::endl;

    // std::cout << "alpha: " << alpha << " n: " << n << " x_bar: " << x_bar << " sigma: " << sigma << " N: " << N << " B: " << B << " A: " << A << std::endl;

    double ret_val;
    if ((x - x_bar) / sigma < alpha)
    {
        // std::cout << "path a" << std::endl;
        ret_val = exp(-0.5 * (x - x_bar) * (x - x_bar) / (sigma * sigma));
    }
    else
    {
        // std::cout << "path b" << std::endl;
        ret_val = A * pow(B + (x - x_bar) / sigma, -1 * n);
    }
    ret_val = N * ret_val + offset;
    // std::cout << "x: " << x << " y: " << ret_val << std::endl;
    return ret_val;
}

float sigma_cut = 2;

int sipms_to_use = 16;

float center_x = 1.95396;
float sigma_x = 0.184758 * sigma_cut;
float center_y = 1.9688;
float sigma_y = 0.195137 * sigma_cut;

// float adc_calib = 26704.4;  // 32444.1 Signal_ADC = 1 GeV
// float adc_calib = 57854.3;  // v3

int signal_method = 3;









float calculate_signal_v6(uint32_t *adc_values, uint32_t *common_mode, float gain)
{
    // Pedestal is the mean of the first three samples
    float pedestal = (adc_values[0] + adc_values[1] + adc_values[2]) / 3.0f;
    float common_mode_pedestal = (common_mode[0] + common_mode[1] + common_mode[2]) / 3.0f;

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
        float sample = adc_values[i] - pedestal;
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
    TF1 *crystal_ball_fit = new TF1("crystal_ball", crystal_ball, 4, 20, 6);
    crystal_ball_fit->SetParameters(1, 1, 6, 4, 2);
    crystal_ball_fit->SetParLimits(0, 1, 1.2);     // alpha
    crystal_ball_fit->SetParLimits(1, 0.2, 0.6);   // n
    crystal_ball_fit->SetParLimits(2, 0.5, 4.5);   // x_bar
    crystal_ball_fit->SetParLimits(3, 0.25, 0.65); // sigma
    crystal_ball_fit->SetParLimits(4, 0, 2000);    // N
    float pedestal = (adc_values[0] + adc_values[1] + adc_values[2]) / 3.0f;
    crystal_ball_fit->FixParameter(5, pedestal);

    // Fit the samples from 3 to 19
    TH1F *temp_hist = new TH1F("temp_hist", "temp_hist", 20, 0, 20);
    for (int i = 3; i < 20; ++i)
    {
        float sample = adc_values[i];
        temp_hist->SetBinContent(i, sample);
    }
    temp_hist->Fit(crystal_ball_fit, "RQ");
    // std::cout << "Fit parameters: ";
    // for (int i = 0; i < 5; ++i) {
    //     std::cout << crystal_ball_fit->GetParameter(i) << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "Fit chi2/ndf: " << crystal_ball_fit->GetChisquare() / crystal_ball_fit->GetNDF() << std::endl;
    float signal = crystal_ball_fit->GetParameter(4);
    ;
    delete temp_hist;
    delete crystal_ball_fit;
    return signal * gain;
}

float calculate_signal(uint32_t *adc_values, float gain)
{
    if (signal_method == 2)
    {
        return calculate_signal_v2(adc_values, gain);
    }
    else if (signal_method == 3)
    {
        return calculate_signal_v3(adc_values, gain);
    }
    else if (signal_method == 4)
    {
        return calculate_signal_v4(adc_values, gain);
    }
    else if (signal_method == 5)
    {
        return calculate_signal_v5(adc_values, gain);
    }
    else if (signal_method == 7)
    {
        return calculate_signal_v7(adc_values, gain);
    }
    else
    {
        // Default to v2
        return calculate_signal_v2(adc_values, gain);
    }
}

bool calculate_signal(uint32_t *adc_values, uint32_t *tot_values, float gain, float &signal)
{
    if (signal_method == 2)
    {
        // Check if there is a ToT value
        uint32_t tot = 0;
        for (int i = 0; i < 20; ++i)
        {
            if (tot_values[i] > 0)
            {
                tot = tot_values[i];
                break;
            }
        }
        // If there is no ToT value
        if (tot == 0)
        {
            signal = calculate_signal(adc_values, gain);
            return false;
        }

        // Otherwise, just return the first ToT value
        signal = tot * gain;
        return true;
    }
    else if (signal_method == 3)
    {
        // Check if there is a ToT value
        uint32_t tot = 0;
        for (int i = 0; i < 20; ++i)
        {
            if (tot_values[i] > 0)
            {
                tot = tot_values[i];
                break;
            }
        }
        // If there is no ToT value
        if (tot == 0)
        {
            signal = calculate_signal_v3(adc_values, gain);
            return false;
        }
        signal = tot;
        return true;
    }
    else if (signal_method == 4)
    {
        // Check if there is a ToT value
        uint32_t tot = 0;
        for (int i = 0; i < 20; ++i)
        {
            if (tot_values[i] > 0)
            {
                tot = tot_values[i];
                break;
            }
        }
        // If there is no ToT value
        if (tot == 0)
        {
            signal = calculate_signal_v4(adc_values, gain);
            return false;
        }
        signal = tot;
        return true;
    }
    else if (signal_method == 5)
    {
        // Check if there is a ToT value
        uint32_t tot = 0;
        for (int i = 0; i < 20; ++i)
        {
            if (tot_values[i] > 0)
            {
                tot = tot_values[i];
                break;
            }
        }
        // If there is no ToT value
        if (tot == 0)
        {
            signal = calculate_signal_v5(adc_values, gain);
            return false;
        }
        signal = tot;
        return true;
    }
    else
    {
        // Default to ADC signal
        signal = calculate_signal(adc_values, gain);
        return false;
    }
}

bool is_tot(uint32_t *tot_values)
{
    for (int i = 0; i < 20; ++i)
    {
        // std::cout << "TOT Value[" << i << "]: " << tot_values[i] << std::endl;
        if (tot_values[i] > 50)
        {
            return true;
        }
    }
    return false;
}

int32_t get_toa(uint32_t *toa_values)
{
    for (int i = 0; i < 20; ++i)
    {
        if (toa_values[i] > 0)
        {
            return toa_values[i];
        }
    }
    return -1;
}

bool fit_peak(TH1 *hist)
{
    float lower_bound = hist->GetXaxis()->GetXmin();
    float upper_bound = hist->GetXaxis()->GetXmax();
    float mean = hist->GetMean();
    // Check that the mean is within the histogram range
    float rms = hist->GetRMS();
    float fit_min = mean - 1.5 * rms;
    float fit_max = mean + 1.5 * rms;
    TF1 *rough_fit = new TF1("rough_fit", "gaus", fit_min, fit_max);
    // std::cout << "fit 1 range: " << fit_min << " to " << fit_max << std::endl;
    auto result = hist->Fit(rough_fit, "RQS");
    if (result->Status() != 0 || rough_fit->GetParameter(1) < lower_bound || rough_fit->GetParameter(1) > upper_bound)
    {
        std::cout << "failed!" << std::endl;
        return true;
    }
    float peak = rough_fit->GetParameter(1);
    float sigma = rough_fit->GetParameter(2);
    rough_fit->SetLineColor(kBlue);
    rough_fit->SetLineStyle(2);

    TF1 *second_fit = new TF1("second_fit", "gaus", peak - sigma, peak + sigma);
    // std::cout << "fit 2 range: " << peak - sigma << " to " << peak + sigma << std::endl;
    result = hist->Fit(second_fit, "RQS");
    if (result->Status() != 0 || second_fit->GetParameter(1) < lower_bound || second_fit->GetParameter(1) > upper_bound)
    {
        std::cout << "failed!" << std::endl;
        return true;
    }
    peak = second_fit->GetParameter(1);
    sigma = second_fit->GetParameter(2);

    // TF1* final_fit = new TF1("final_fit", "crystalball", peak - sigma, peak + sigma);
    // final_fit->SetParameters(second_fit->GetParameter(0), peak, sigma, 1.5, 2.0);
    TF1 *final_fit = new TF1("final_fit", "gaus", peak - sigma, peak + sigma);
    // std::cout << "fit 3 range: " << peak - sigma << " to " << peak + sigma << std::endl;
    result = hist->Fit(final_fit, "RQS");
    if (result->Status() != 0 || final_fit->GetParameter(1) < lower_bound || final_fit->GetParameter(1) > upper_bound)
    {
        std::cout << "failed!" << std::endl;
        return true;
    }
    return false;
}

void draw_text(TF1 *fit, int run_number, float beam_energy)
{
    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.04);
    text->SetTextFont(42);
    text->SetTextAlign(31);

    float peak = fit->GetParameter(1);
    ;
    float sigma = fit->GetParameter(2);
    float resolution = (sigma / peak) * 100.0f;
    text->DrawLatex(0.93, 0.85, Form("%.01f GeV Electrons", beam_energy));
    text->DrawLatex(0.93, 0.80, Form("Run %d", run_number));
    text->DrawLatex(0.93, 0.75, Form("Peak: %.03f", peak));
    text->DrawLatex(0.93, 0.70, Form("Sigma: %.03f", sigma));
    text->DrawLatex(0.93, 0.65, Form("Resolution: %.2f%%", resolution));
    text->DrawLatex(0.93, 0.60, Form("Signal method %d", signal_method));
}

bool position_cut(float x, float y)
{
    // Cut anything outside of the center
    if (std::abs(x - center_x) > sigma_x || std::abs(y - center_y) > sigma_y)
    {
        return false;
    }
    return true;
}

bool calculate_cog(TH2 *distribution, float *values)
{
    float total_signal = 0;
    float x_weighted_sum = 0;
    float y_weighted_sum = 0;
    float x_cog = 0;
    float y_cog = 0;
    float w = 4;

    for (int i = 0; i < 25; i++)
    {
        total_signal += values[i];
    }

    float total_weight = 0;
    for (int i = 0; i < 25; i++)
    {
        int x = i % 5;
        int y = i / 5;
        float signal = values[i];
        float weight = w + std::log(signal / total_signal); // Avoid log(0)
        if (weight < 0)
        {
            weight = 0;
        }
        total_weight += weight;
        x_weighted_sum += x * weight;
        y_weighted_sum += y * weight;
    }
    if (total_weight > 0)
    {
        x_cog = x_weighted_sum / total_weight;
        y_cog = y_weighted_sum / total_weight;
        distribution->Fill(x_cog, y_cog);
    }
    return position_cut(x_cog, y_cog);
}

void print_progress(int progress)
{
    std::cout << " [";
    for (int i = 0; i < 25; i++)
    {
        if (i < progress)
        {
            std::cout << "*";
        }
        else
        {
            std::cout << " ";
        }
    }
    std::cout << "]\r" << std::flush;
}