# eeemcal_desy_dec2025

Example of analysis
```
root -l -b

.L common_led.cxx+
.L gain_match.cxx+
.L adc_calibration.cxx+
.L tot_calibration.cxx+
.L energy_resolution.cxx+

gain_scan("data", "eeemcal_desy_dec2025_mapping_v2.csv", "outputs", true);

adc_calib_scan("data", "eeemcal_desy_dec2025_mapping_v2.csv", true);

tot_calib_scan("data",
               "eeemcal_desy_dec2025_mapping_v2.csv",
               "outputs",
               "outputs/gain_match_1.259V.root",
               "outputs/adc_to_ref_calibration_1.259V.root",
               "outputs/tot_widths.root",
               12);

energy_resolution_led_scan("data",
                           "eeemcal_desy_dec2025_mapping_v2.csv",
                           "outputs",
                           "outputs/adc_to_ref_calibration_1.259V.root",
                           "outputs/tot_calibration_values.root",
                           12);

.q
```