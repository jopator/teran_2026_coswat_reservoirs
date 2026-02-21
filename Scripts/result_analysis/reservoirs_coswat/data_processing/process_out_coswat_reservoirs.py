"""
Author: Jose P. Teran
Date: 2025-11-27

Description:
Process monthly and yearly reservoir storage, inflow, outflow, evaporation and seepage from CoSWAT+

It will generate pikle files per region for all reservoirs/lakes in the model (not just where observaitions exist)!

The result is a pikle with a dictionary for each reservoir - It is referred to coswat id, but also specifies HydroLakes and GranD id

It takes --r <region> and --v <version> as argument
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import pickle as pkl
import read_swat
import argparse

def process_res(filePath:str, res_id_df: pd.DataFrame):
    file_name = filePath.split("/")[-1]
    print(f"Getting data from {file_name} into a dictionary")
    res_table       = read_swat.swat_table(filePath)
    out_dict = {} # Empty to store per reservoir

    for idx,row in res_id_df.iterrows():
        reservoir_id = row['LakeId']
        hylak_id     = row['Hylak_id']
        grand_id     = row['Grand_id']

        df = res_table.obj_output(reservoir_id,var_list)
        out_dict[reservoir_id] = {  "Hylak_id":int(hylak_id),
                                    "grand_id":int(grand_id),
                                    "coswatOut":df
                                 }

    return out_dict

if __name__ == '__main__':

    os.chdir('/data/brussel/vo/000/bvo00033') # All relative paths will be based on this ! !
    analysis_folder = f'vsc10883/PHD_main/Projects/CoSWAT/Scripts/result_analysis'

    # argument and regions
    parser = argparse.ArgumentParser(description="scripts for running coswatv2")
    parser.add_argument("--r", help="region", type=str, default='')
    parser.add_argument("--v", help="version", type=str, default='')
    args = parser.parse_args()

    region  = args.r
    version = args.v

    if version == "" or region == "":
        print("You did not specify region or version: do --r <region> and --v <version>")
        quit()
    
    print(f"=== Processing region: {region} - CoSWATv{version} ===")



    # Paths and file names
    # ====================

    # Variables
    var_list        = ['precip','evap','seep','flo_stor','flo_in','flo_out','sed_stor','sed_in','sed_out']
    unit_list       = ["m3","m3","m3","m3","m3","m3","ton","ton","ton"]      # Must be aligned to var_list

    # CoSWAT
    modelSetupPath  = f'CoSWAT-Framework/model-setup/CoSWATv{version}/{region}'
    txtInOutPath    = f'{modelSetupPath}/Scenarios/Default/TxtInOut'
    wshedGisPath    = f'{modelSetupPath}/Watershed'

    reservoir_mon   = f'{txtInOutPath}/reservoir_mon.txt'
    reservoir_year  = f'{txtInOutPath}/reservoir_yr.txt'

    reservoir_shp   = f'{wshedGisPath}/Shapes/lakes-grand-ESRI-54003.shp'

    # Output generation
    outFolder       = f'{analysis_folder}/reservoirs_coswat/coswat_outputs/reservoir_storage/CoSWATv{version}/{region}'
    outMonthly      = f'{outFolder}/monthly'
    outYearly       = f'{outFolder}/yearly'

    outFileMonthly  = f'{outMonthly}/reservoir_mon.pkl'
    outFileYearly   = f'{outYearly}/reservoir_yr.pkl'

    # Processing
    #====================
    
    # Create Out folders if they don't exist
    os.makedirs(outMonthly, exist_ok = True)
    os.makedirs(outYearly, exist_ok  = True)


    # Get ids (Coswat, Hylak and GranD)
    res_gdf = gpd.read_file(reservoir_shp)
    res_id_df = res_gdf[['LakeId','Hylak_id','Grand_id']]

    # Get data dictionaries
    monthlyDict = process_res(reservoir_mon,res_id_df)
    yearlyDict  = process_res(reservoir_year,res_id_df)

    # Dump to pickle
    print("> Dumping to pickle files...")
    with open(outFileMonthly,'wb') as f:
        pkl.dump(monthlyDict,f)
    
    with open(outFileYearly,'wb') as f:
        pkl.dump(yearlyDict,f)

    print("Done!")
