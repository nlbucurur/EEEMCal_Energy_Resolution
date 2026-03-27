# eeemcal_desy_dec2025

## Step-by-step procedure to run

First, set your subset in `common_led.cxx`, for example:

example: `std::set<int> g_selected_central_sipms = {0, 3, 7, 15};`

the central crystal id i 12

### 0. Data

Save your root files in a folder called `data`

### 1. Generate the mapping

Modify `generate_mapping.py` depending on the setup of the detector, it will return a `csv` file containing the channels correpondent to the different crystals involved in the analysis.

```
python3 generate_mapping.py
```

The name of the csv file is `eeemcal_desy_dec2025_mapping_v2.csv`

### 2. Compile and test the common selector

Open ROOT:

```
    root -l -b
```

Load the helpers

```
    .L common_led.cxx+
```

### 3. Optional sanity check with waveforms and generate root files with adc and tot distributions

This is only to visually confirm that the selected subset is what you expect.

```
    .L draw_waveform_reversed.cxx+
    draw_waveform_reversed(171, 1.259, "data", "eeemcal_desy_dec2025_mapping_v2.csv", "outputs", g_sipms_to_use)
```

or if you have a list of runs you want to inspect.

```
    led_scan()
```

to generate the root files

```
    root -l -b
    .L common_led.cxx+
    .L led_analysis.C+
    led_scan()
    .q
```

### 4. Run gain matching

Still in ROOT, or in a fresh ROOT session:

```
    .L common_led.cxx+
    .L gain_match.cxx+
    gain_scan("data", "eeemcal_desy_dec2025_mapping_v2.csv", "outputs", true)
```

This produces the per-voltage gain files such as `outputs/gain_match_1.259V.root`.

### 5. Run ADC calibration

This uses the gain-match outputs and produces the ADC calibration ROOT files.

```
    .L common_led.cxx+
    .L adc_calibration.cxx+
    adc_calib_scan("data", "eeemcal_desy_dec2025_mapping_v2.csv", true)
```

This should produce files like `outputs/adc_to_ref_calibration_1.259V.root`

### 6. Run ToT calibration

It expects the gain file plus the ADC calibration file as inputs.

```
    .L common_led.cxx+
    .L tot_calibration.cxx+
    tot_calib_scan("data",
                "eeemcal_desy_dec2025_mapping_v2.csv",
                "outputs",
                "outputs/gain_match_1.259V.root",
                "outputs/adc_to_ref_calibration_1.259V.root",
                "outputs/tot_widths.root",
                12)
```

This produces the ToT calibration values that `energy_resolution.cxx` uses.

### 7. Run energy resolution

```
    .L common_led.cxx+
    .L energy_resolution.cxx+
    energy_resolution_led_scan("data",
                            "eeemcal_desy_dec2025_mapping_v2.csv",
                            "outputs",
                            "outputs/adc_to_ref_calibration_1.259V.root",
                            "outputs/tot_calibration_values.root",
                            12)
```

This is the last stage in the chain and it uses the gain files, ADC calibration, and ToT calibration together.

## Recommender order in a ROOT session

```
    root -l -b
    .L common_led.cxx+
    .L gain_match.cxx+
    gain_scan("data", "eeemcal_desy_dec2025_mapping_v2.csv", "outputs", true)

    .L adc_calibration.cxx+
    adc_calib_scan("data", "eeemcal_desy_dec2025_mapping_v2.csv", true)

    .L tot_calibration.cxx+
    tot_calib_scan("data", "eeemcal_desy_dec2025_mapping_v2.csv", "outputs", "outputs/gain_match_1.259V.root", "outputs/adc_to_ref_calibration_1.259V.root", "outputs/tot_widths.root", 12)

    .L energy_resolution.cxx+
    energy_resolution_led_scan("data", "eeemcal_desy_dec2025_mapping_v2.csv", "outputs", "outputs/adc_to_ref_calibration_1.259V.root", "outputs/tot_calibration_values.root", 12)

    .q
```