"""
Copy weather files from source & QC

...


Jose P. Teran - Nov 2025
"""

from pathlib import Path
from setup import *
import os
import shutil
from datetime import datetime
from jopato_fx import *
import argparse

processes   = nr_processes
now         = datetime.now()
s_now       = now.strftime("%Y-%m-%d %H:%M:%S")

if all_regions:
    print('\n Copying weather for all regions in CoSWAT+  \n\n')

    # with open(all_regions_list_file,'r') as f:
    #     regions_raw = f.readlines()
    #     regions = []
    #     for region in regions_raw:
    #         region_ok = region.split('\n')[0]
    #         regions.append(region_ok)


else:
    print('\n Copying weather for selected regions \n\n')
    # regions = selected_regions



if __name__ == '__main__':

    # argument and regions
    parser = argparse.ArgumentParser(description="scripts for running coswatv2")
    parser.add_argument("--r", help="region", type=str, default='')
    args = parser.parse_args()



    if len(args.r) > 0: region = args.r
    else: quit()

# for region in regions:

    if copy_new_files:
        src = f'{source_dir}/{path_to_files.format(region=region)}'
        dst = f'{sim_regions_path}/{region}/Scenarios/Default/TxtInOut'

        srcDir = Path(src)
        dstDir = Path(dst)

        print(f'\n Copying weather files from source for {region}... \n')
        copyWeatherFiles(srcDir,dstDir)

        print(f'\n Sorting Weather files Lexicographically for {region}... \n')
        sortLexi(dstDir)
    

        

# # if copy_new_files:
#     # Copy
#     print('\n Copying weather files from source ... \n\n')
#     with ThreadPoolExecutor(max_workers=processes) as ex:
#         futures = {
#             ex.submit(submitCopy, src, dst): (src, dst)
#             for src, dst in zip(src_list, dst_list)
#         }

#         for fut in as_completed(futures):
#             src, dst = futures[fut]
#             try:
#                 fut.result()
#             except Exception as e:
#                 print(f"Error copying from {src} to {dst}: {e}")

#     print("\n \n > Done \n")

# # Lexico
#     print('\n Sorting Lexicographically ... \n\n')
#     with ThreadPoolExecutor(max_workers=processes) as ex:
#         futures = {ex.submit(submitLexi, dst): dst for dst in dst_list}

#         for fut in as_completed(futures):
#             dst = futures[fut]
#             try:
#                 fut.result()
#             except Exception as e:
#                 print(f"Error sorting files in {dst}: {e}")

#     print("\n \n > Done \n")
