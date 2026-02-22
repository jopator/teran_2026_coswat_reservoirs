"""
Author: Jose P. Teran
Date: 2025-11-27

Description:
Process monthly streamflow at selected stations that are actually snapped to a CoSWAT+ Outlet
It will generate pikle with GRDC station, and a dataframe with simulated and reference time-series; will also provide CoSWAT+ channel Id
Will also flag it if it is downstream of a reservoir

It takes --r <region> and --v <version> as argument
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import pickle as pkl
import read_swat
import argparse
import matplotlib.pyplot as plt
import hydroeval as he
import warnings
from pathlib import Path
warnings.simplefilter("ignore", FutureWarning)


if __name__ == '__main__':

    BASE_DIR = Path(__file__).resolve().parents[3]  # root of repo
    os.chdir(BASE_DIR)  # All relative paths will be based on this ! !

    analysis_folder = f'Scripts/result_analysis'

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
    var_list        = ['flo_out']
    unit_list       = ["m3/s"]      # Must be aligned to var_list
    snap_dist       = 10000         # in meters
    grdc_start      = 1970          # Filter for station time coverage
    grdc_end        = 2015          # Filter for station time coverage
    grdc_minYrs     = 5
    
    # CoSWAT
    modelSetupPath  = f'CoSWAT-Framework/model-setup/CoSWATv{version}/{region}'
    modelDataPath   = f'CoSWAT-Framework/model-data/{region}'

    txtInOutPath    = f'{modelSetupPath}/Scenarios/Default/TxtInOut'
    wshedGisPath    = f'{modelSetupPath}/Watershed'

    channelsd_mon   = f'{txtInOutPath}/channel_sd_mon.txt'

    grdc_shp        = f'{modelDataPath}/shapes/grdc_stations-ESRI-54003.gpkg'
    outlet_shp      = f'{wshedGisPath}/Shapes/outlets_sel_snap.shp'
    rivs_shp        = f'{wshedGisPath}/Shapes/rivs1.shp'
    chan_shp        = f'{wshedGisPath}/Shapes/dem-aster-ESRI-54003-lakeBurntchannel.shp' #f'{wshedGisPath}/Shapes/dem-aster-ESRI-54003channel.shp' #f'{wshedGisPath}/Shapes/dem-aster-ESRI-54003channel.shp'  #f'{wshedGisPath}/Shapes/dem-aster-ESRI-54003-lakeBurntchannel.shp'
    reservoir_shp   = f'{wshedGisPath}/Shapes/lakes-grand-ESRI-54003.shp'

    grdc_data_folder = f'{modelDataPath}/observations'

    # Output generation
    outFolder       = f'{analysis_folder}/reservoirs_coswat/coswat_outputs/streamflow/CoSWATv{version}/{region}'
    outMonthly      = f'{outFolder}/monthly'

    out_shp         = f'{outMonthly}/eval_stations.gpkg'
    out_csv         = f'{outMonthly}/eval_stations.csv'
    outFileMonthly  = f'{outMonthly}/channelsd_mon.pkl'


    # Processing
    #====================
    
    # Create Out folders if they don't exist
    os.makedirs(outMonthly, exist_ok = True)

    # > SNAP STATIONS TO OUTLETS <
    # -----------------------------
    print("> Snapping GRDC stations to outlets ...")
    # Snap grdc stations to CoSWAT outlets
    outlets_gdf = gpd.read_file(outlet_shp)
    grdc_gdf    = gpd.read_file(grdc_shp)
    rivs_gdf    = gpd.read_file(rivs_shp)
    res_gdf     = gpd.read_file(reservoir_shp)
    chan_gdf    = gpd.read_file(chan_shp)
    
    # Filter out stations outside of the simulation period or with not enough years
    grdc_gdf = grdc_gdf[(grdc_gdf['t_start']<grdc_end) & (grdc_gdf['t_end']>grdc_start)].copy()
    grdc_gdf = grdc_gdf[(grdc_gdf['d_yrs']>grdc_minYrs) | (grdc_gdf['m_yrs']>grdc_minYrs)].copy()
    

    # Filter out stations for which we have the time-series in model-data
    grdc_ts_list = os.listdir(grdc_data_folder)
    grdc_ts_list = [int(s.replace(".csv", "")) for s in grdc_ts_list]
    grdc_gdf = grdc_gdf[grdc_gdf['grdc_no'].isin(grdc_ts_list)] # Here we only take those and this will make sure for the rest of the code
    
    
    # Match
    matches = gpd.sjoin_nearest(grdc_gdf, outlets_gdf, how='inner', max_distance = snap_dist, distance_col = 'snap_dist')
    matches = matches.sort_values("snap_dist").drop_duplicates("index_right") # Keep closest
    outlets_snapped = outlets_gdf.loc[matches["index_right"]].copy()

    # Get GRDC properties
    matches = matches.set_index("index_right")
    outlets_snapped["grdc_no"] = matches["grdc_no"]
    outlets_snapped["river"]   = matches["river"]
    outlets_snapped["station"] = matches["station"]

    # Save figure of snapping points
    bbox = rivs_gdf.total_bounds   # Network bounds
    minx, miny, maxx, maxy = bbox
    grdc_subset = grdc_gdf.cx[minx:maxx, miny:maxy]
    fig, ax = plt.subplots(figsize=(10,10),dpi=300)
    outlets_gdf.plot(ax=ax,markersize=10,color='gray',alpha=0.5,label='Model Outlets')
    grdc_subset.plot(ax=ax,markersize=10,color='orange',label='GRDC Stations')
    outlets_snapped.plot(ax=ax,markersize=8,color='green',alpha=0.8,label='Snapped Stations')
    rivs_gdf.plot(ax=ax,color='navy',linewidth=0.2,label='Model River Network')
    ax.legend()
    plt.savefig(f'{outMonthly}/outlet_snap_grdc_{region}.jpg')
    
    if version != "1.1.0":

        # > DEFINE RESERVOIR INFLUENCE<
        # -----------------------------
        print("> Asessing relation with reservoirs ...")
        # Get outlets that are or are not dowstream (influenced by a reservoir)
        links = []
        for idx, row in res_gdf.iterrows():
            main_use    = row['MAIN_USE']

            lake_id         = row['LakeId']
            lake_gdf        = res_gdf[res_gdf['LakeId']==lake_id]                         # get lakes
            lakeOutlet_gdf  = chan_gdf[chan_gdf['LakeMain']==lake_id]                     # get main outlet
            lakeOutlet_gdf  = lakeOutlet_gdf.reset_index(drop=True)
            complete        = False

            outletLINKNO    = lakeOutlet_gdf['LINKNO'].iloc[0]
            outletDSLINKNO  = lakeOutlet_gdf['DSLINKNO'].iloc[0]
            outletLENGTH    = lakeOutlet_gdf['Length'].iloc[0]

            total_length = outletLENGTH
            links.append(outletLINKNO)

            currentLINKNO = outletDSLINKNO

            while not complete:
                new_chan_gdf    = chan_gdf[chan_gdf['LINKNO'] == currentLINKNO].reset_index(drop=True)
                new_chan_lenght = new_chan_gdf['Length'].iloc[0]         
                total_length += new_chan_lenght

                currentUSLINKNO1 = new_chan_gdf['USLINKNO1'].iloc[0]
                currentUSLINKNO2 = new_chan_gdf['USLINKNO2'].iloc[0]

                links.append(currentLINKNO)


                currentLINKNO = new_chan_gdf['DSLINKNO'].iloc[0]
                lakeInFlag    = new_chan_gdf['LakeIn'].iloc[0]


                if currentLINKNO == -1:
                    break

                if lakeInFlag != 0:
                    break

        downs_chans_gdf = rivs_gdf[rivs_gdf['LINKNO'].isin(links)]

        # Get channels that will be evaluated (By DSNODEID)
        eval_chans = chan_gdf[chan_gdf['DSNODEID'].isin(outlets_snapped['ID'].to_list())]
        eval_rivs  = rivs_gdf[rivs_gdf['LINKNO'].isin(eval_chans['LINKNO'].to_list())]
        eval_rivs  = pd.merge(eval_rivs,eval_chans[['LINKNO','DSNODEID']],on='LINKNO',how='left') # get ID in eval_rivs

        # Differentiate stations between downstream of a reservoir of not
        down_chans_links = downs_chans_gdf['LINKNO'].to_list()
        eval_rivs['down_res'] = False
        eval_rivs["down_res"] = eval_rivs["LINKNO"].isin(down_chans_links)
        
        eval_rivs       = eval_rivs.reset_index(drop=True)
        eval_rivs_df    = eval_rivs[['LINKNO','Channel','DSNODEID','down_res']].rename(columns={'DSNODEID':'ID'})
        eval_rivs_gdf   = gpd.GeoDataFrame(data=eval_rivs_df,geometry=eval_rivs.geometry)

        outlets_snapped = pd.merge(outlets_snapped,eval_rivs_df,on='ID',how='left')
        outlets_snapped = outlets_snapped.reset_index(drop=True)
        
        outlets_snapped['down_res'] = outlets_snapped['down_res'].fillna(False).astype(bool)

        eval_stats_df   = outlets_snapped[['ID','grdc_no','river','station','LINKNO','Channel','down_res']]
        eval_stats_gdf  = gpd.GeoDataFrame(data=eval_stats_df,geometry=outlets_snapped.geometry)

        # Save figure of potential reservoir influence in streamflow
        fig, ax = plt.subplots(figsize=(10,10),dpi=300)
        downs_chans_gdf.plot(ax=ax,color='salmon',linewidth=0.8,label="River Influenced by Reservoir")
        rivs_gdf[~rivs_gdf['LINKNO'].isin(links)].plot(ax=ax,color='gray',linestyle='-',linewidth=0.4,label="River Not Influenced by Reservoir")
        res_gdf.plot(ax=ax,color='skyblue')
        eval_stats_gdf[~eval_stats_gdf['down_res']].plot(ax=ax,markersize=10,marker='^',color='k',label='Station Upstream of Reservoir')
        eval_stats_gdf[eval_stats_gdf['down_res']].plot(ax=ax,markersize=15,marker='*',color='green',label='Station Downstream of Reservoir')
        ax.collections[-1].set_zorder(10)
        ax.legend()
        plt.savefig(f'{outMonthly}/reservoir_influence_{region}.jpg')
    
    else:
        # Get channels that will be evaluated (By DSNODEID)
        eval_chans      = chan_gdf[chan_gdf['DSNODEID'].isin(outlets_snapped['ID'].to_list())]
        eval_rivs       = rivs_gdf[rivs_gdf['LINKNO'].isin(eval_chans['LINKNO'].to_list())]
        eval_rivs       = pd.merge(eval_rivs,eval_chans[['LINKNO','DSNODEID']],on='LINKNO',how='left') # get ID in eval_rivs
        eval_rivs_df    = eval_rivs[['LINKNO','Channel','DSNODEID']].rename(columns={'DSNODEID':'ID'})
        eval_rivs_gdf   = gpd.GeoDataFrame(data=eval_rivs_df,geometry=eval_rivs.geometry)
        outlets_snapped = pd.merge(outlets_snapped,eval_rivs_df,on='ID',how='left')
        outlets_snapped = outlets_snapped.reset_index(drop=True)
        outlets_snapped['down_res'] = False
        eval_stats_df   = outlets_snapped[['ID','grdc_no','river','station','LINKNO','Channel','down_res']]
        eval_stats_gdf  = gpd.GeoDataFrame(data=eval_stats_df,geometry=outlets_snapped.geometry)
        
        print(eval_stats_df)

    # >PROCESS OUTPUTS<
    
    def calcSats(val_df, sim_col, obs_col):
        # Force numeric
        s = pd.to_numeric(val_df[sim_col], errors="coerce")
        o = pd.to_numeric(val_df[obs_col], errors="coerce")

        # Keep only rows where both are finite numbers
        mask = np.isfinite(s) & np.isfinite(o)
        sim_arr = s[mask].to_numpy(dtype=float)
        obs_arr = o[mask].to_numpy(dtype=float)

        if sim_arr.size == 0:
            raise ValueError(f"No valid overlapping numeric data for {sim_col} and {obs_col}")
            
        kge, r, alpha, beta = he.evaluator(he.kge,sim_arr,obs_arr)
        pbias = he.evaluator(he.pbias,sim_arr,obs_arr)

        return kge[0], pbias[0], r[0], alpha[0], beta[0]


    print("> Processing time series ...")
    # -----------------------------
    streamflow_dict = {}
   
    for idx,row in eval_stats_df.iterrows():
        grdc_no = row['grdc_no']
        channel = row['Channel']
        observed_fn = f'{grdc_data_folder}/{grdc_no}.csv'
        observed_df = pd.read_csv(observed_fn)
        
        if ' Original' in observed_df.columns:
            observed_df = observed_df[['YYYY-MM-DD',' Original']].copy().rename(columns={'YYYY-MM-DD':'date',' Original':'obs'})
        
        elif ' Value' in observed_df.columns:
            observed_df = observed_df[['YYYY-MM-DD',' Value']].copy().rename(columns={'YYYY-MM-DD':'date',' Value':'obs'})
            
        observed_df['date'] = pd.to_datetime(observed_df['date'])
        observed_df = observed_df[['date','obs']].copy()
        sim_df = read_swat.swat_table(channelsd_mon).obj_output(channel,var_list)
        sim_df = sim_df.copy().rename(columns = {'flo_out':'sim'})

        streamflow_dict[grdc_no] ={'sim_data':sim_df,'obs_data':observed_df,'Channel':channel,'ID':row['ID'],'station':row['station'],'down_res':row['down_res']}
    
    # Save dataframe to csv to be used later
    eval_stats_gdf.to_csv(f'{out_shp}')
    eval_stats_df.to_csv(f'{out_csv}')
    
    # Dump to pickle
    print("> Dumping to pickle file...")
    with open(outFileMonthly,'wb') as f:
        pkl.dump(streamflow_dict,f)