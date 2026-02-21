"""
Write files:
time.sim
print.prt
codes.bsn

...


Jose P. Teran - Nov 2025
"""


from setup import *
import os
from datetime import datetime
from jopato_fx import *
import argparse

now     = datetime.now()
s_now   = now.strftime("%Y-%m-%d %H:%M:%S")

if all_regions:
    print('\n Writing files for all regions in CoSWAT+ \n\n')

    # with open(all_regions_list_file,'r') as f:
    #     regions_raw = f.readlines()
    #     regions = []
    #     for region in regions_raw:
    #         region_ok = region.split('\n')[0]
    #         regions.append(region_ok)


else:
    print('\n Writing files for selected regions \n\n')
    # regions = selected_regions


# for region in regions:


 # get model setup version

if __name__ == '__main__':

    # argument and regions
    parser = argparse.ArgumentParser(description="scripts for running coswatv2")
    parser.add_argument("--r", help="region", type=str, default='')
    args = parser.parse_args()



    if len(args.r) > 0: region = args.r
    else: quit()

    txtInOut_path = f'{sim_regions_path}/{region}/Scenarios/Default/TxtInOut'
    # time.sim
    print('\n Writing time.sim ... \n\n')
    time_sim_file = f'{txtInOut_path}/time.sim'

    with open(time_sim_file,'r') as f:
        lines = f.readlines()
        cols  = lines[1]

    with open(time_sim_file,'w') as f:
        f.write(f'time.sim written by Jose P Teran - CoSWAT - {s_now}\n')
        f.write(cols)
        yrc_start = time_sim['yrc_start']
        yrc_end   = time_sim['yrc_end']

        f.write(f'       0      {yrc_start}         0        {yrc_end}         0 ')


    # print.prt
    print('\n Writing print.prt ... \n\n')
    print_prt_file = f'{txtInOut_path}/print.prt'
    out = update_print_prt(print_prt_file,print_prt)
    out[0] = f'print.prt written by Jose P Teran - CoSWAT - {s_now}\n'

    with open(print_prt_file,'w') as f:
        f.writelines(out)

    # codes.bsn
    print('\n Writing codes.bsn ... \n\n')
    codes_bsn_file = f'{txtInOut_path}/codes.bsn'

    with open(codes_bsn_file,'r') as f:
        lines = f.readlines()

    cols        = lines[1]
    data        = lines[2]
    col_names   = cols.split()
    vals        = data.split()

    for i, name in enumerate(col_names):
        if name in codes_bsn:
            vals[i] = str(codes_bsn[name])

    prefix = data[:len(data) - len(data.lstrip())]
    new_data_line = prefix + "  ".join(vals) + "\n"
    lines[2] = new_data_line

    lines[0] = f'codes.bsn written by Jose P Teran - CoSWAT - {s_now}\n'

    with open(codes_bsn_file,'w') as f:
        f.writelines(lines)
