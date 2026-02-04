import pandas as pd
import os

## Original code in here: https://github.com/tlprotzman/eeemcal_fast_offline/blob/main/mapping/generate_mapping.py

# EEEMCal mapping - instead of "layers", we have a single plane, where each crystal is one connector
# FPGA IP | ID
# 208     | 0  # In this occassion, we are going to use only 208 (pb 17 during the beam test it was 210 (2))
# 209     | 1
# 210     | 2
# # 211     | 3
# eeemcal_fpga_map = [1, 2, 1, 2, 1, ## In beam test configuration
#                     3, 3, 1, 0, 2,
#                     1, 2, 2, 3, 0,
#                     3, 0, 2, 0, 0,
#                     2, 1, 1, 2, 1]

eeemcal_fpga_map = [0]

# // ASIC | ID
# // 0    | 0
# // 1    | 1
# eeemcal_asic_map = [1, 0, 0, 1, 0, ## In the beam test configuration
#                     1, 1, 1, 1, 1,
#                     0, 1, 0, 0, 0,
#                     0, 1, 0, 0, 1,
#                     1, 0, 1, 0, 1]

eeemcal_asic_map = [0]

# // Connector | ID
# // A        | 0
# // B        | 1
# // C        | 2
# // D        | 3
# eeemcal_connector_map = [0,  0,  3,  1,  0, ## In the beam test configuration
#                          0,  3,  1,  2,  3,
#                          2,  0,  3,  1,  0,
#                          0,  3,  2,  2,  1,
#                          2,  1,  3,  1,  2]
eeemcal_connector_map = [3]

    

eeemcal_16i_channel_a_map = [2,  6, 11, 15,  0,  4,  9, 13,
                             1,  5, 10, 14,  3,  7, 12, 16]

eeemcal_16i_channel_b_map = [20, 24, 29, 33, 18, 22, 27, 31,
                             19, 23, 28, 32, 21, 25, 30, 34]

eeemcal_16i_channel_c_map = [67, 63, 59, 55, 69, 65, 61, 57,
                             70, 66, 60, 56, 68, 64, 58, 54]

eeemcal_16i_channel_d_map = [50, 46, 40, 36, 52, 48, 42, 38,
                             51, 47, 43, 39, 49, 45, 41, 37]

eeemcal_16i_channel_map = [eeemcal_16i_channel_a_map, eeemcal_16i_channel_b_map, eeemcal_16i_channel_c_map, eeemcal_16i_channel_d_map]


# crystal_ID = [4, 9, 14, 19, 24,
#               3, 8, 13, 18, 23,
#               2, 7, 12, 17, 22,
#               1, 6, 11, 16, 21,
#               0, 5, 10, 15, 20]

crystal_ID = [12]

def make_16i_mapping():
    rows = []
    for sipm in range(16):
        fpga = eeemcal_fpga_map[0]
        asic = eeemcal_asic_map[0]
        connector = eeemcal_connector_map[0]
        channel = eeemcal_16i_channel_map[connector][sipm] + 144 * fpga + 72 * asic
        rows.append({'FPGA': fpga, 'ASIC': asic, 'Connector': connector,
                     'Crystal': crystal_ID[0], 'SiPM': sipm, 'Channel': channel})
        
    df = pd.DataFrame(rows)
    print("cwd:", __import__("os").getcwd())
    print("rows:", len(rows))
    print(df.head())
    out = 'eeemcal_desy_dec2025_mapping_v2.csv'
    df.to_csv(out, index=False)
    print("wrote:", out)

def main():
    
    make_16i_mapping()

if __name__ == '__main__':
    main()
