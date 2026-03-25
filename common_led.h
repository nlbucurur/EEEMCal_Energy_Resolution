#ifndef COMMON_LED_H
#define COMMON_LED_H

#include <map>
#include <set>
#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <utility>

class TH1;
class TF1;

// ===== Shared geometry / DAQ constants =====
constexpr int LED_SAMPLES_PER_CHANNEL = 20;
constexpr int LED_SIPMS_PER_CRYSTAL = 16;
constexpr int LED_MAX_NUM_CRYSTALS = 25;

// ===== Configuration =====
extern int g_signal_method; // 2, 3, 4, 5, 7 (7 = waveform crystal ball fit)
extern int g_sipms_to_use; // Number of SiPMs to use per crystal (= 16)
extern uint32_t g_tot_min; // Minimum ToT value to consider valid

// ===== Generic helpers shared by the LED analysis macros =====
int extract_run_number(const char* filename);
void print_progress_25(int progress);
void build_xy_for_crystal(int crystal_id, float& x, float& y);
bool point_in_ellipse(float x, float y, float cx, float cy, float sx, float sy);
void subtract_common_mode(uint32_t out_adc[LED_SAMPLES_PER_CHANNEL],
                          uint32_t* in_adc,
                          float cm);
std::unordered_map<int, std::pair<int, int>> build_reverse_channel_map(
    const std::map<int, std::vector<int>>& mapping,
    int sipms_to_use = LED_SIPMS_PER_CRYSTAL
);

// ===== Mapping =====
// csv format: FPGA,ASIC,Connector,Crystal,SiPM,Channel
// Returns: mapping[crystal][sipm] = channel
std::map<int, std::vector<int>> read_mapping_csv(
    const std::string& filename,
    int expected_sipms_per_crystal = LED_SIPMS_PER_CRYSTAL
);

// Convenience: returns channel list for a given crystal (size = expected_sipms_per_crystal)
// If crystal not present, returns vector filled with -1.
std::vector<int> get_crystal_channels(const std::map<int, std::vector<int>>& mapping,
                                     int crystal,
                                     int expected_sipms_per_crystal = LED_SIPMS_PER_CRYSTAL);

std::vector<int> get_active_channels_from_mapping(const std::map<int, std::vector<int>>& mapping,
                                                  int sipms_to_use = LED_SIPMS_PER_CRYSTAL,
                                                  int n_crystals = LED_MAX_NUM_CRYSTALS,
                                                  int max_channels = 128);

// -------------------- Signal extraction --------------------
float calculate_signal_v2(uint32_t* adc_values, float gain);
float calculate_signal_v3(uint32_t* adc_values, float gain);
float calculate_signal_v4(uint32_t* adc_values, float gain);
float calculate_signal_v5(uint32_t* adc_values, float gain);
float calculate_signal_v7(uint32_t* adc_values, float gain); // waveform Crystal Ball fit

float calculate_signal_adc(uint32_t* adc_values, float gain);

// ToT helpers
bool     has_tot(uint32_t* tot_values);
uint32_t get_tot_first(uint32_t* tot_values);
uint32_t get_tot_max(uint32_t* tot_values);

// Hybrid signal: if ToT present -> return ToT value and used_tot = true, else ADC-derived and used_tot = false
bool calculate_signal_hybrid(uint32_t* adc_values,
                      uint32_t* tot_values,
                      float gain_adc,
                      float& signal_out,
                      bool& used_tot);

// ===== Peak fitting (for LED peaks) =====
// Returns true if fit successful, false otherwise. Outputs peak position and sigma.
bool fit_peak_gaus(TH1* hist, float& peak, float& sigma);

// ===== Plot annotations =====
void draw_text_basic(TF1* fit,
                     const std::string& label1,
                     const std::string& label2 = "",
                     const std::string& label3 = "");

#endif // COMMON_LED_H