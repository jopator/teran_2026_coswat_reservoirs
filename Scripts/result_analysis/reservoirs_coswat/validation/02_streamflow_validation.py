"""
Author: Jose P. Teran
Date: 2025-12-01

Description:
Validate stream flow against reference data (GRDC)

It will generate csv with summary of KGE, PBIAS and R and reservoir influence; this will be used later for a general analysis across all regions
It will generate boxen plots of stats if it exists
It will generate time series plots with simulated vs observed

It takes --r <region> and --v <version> as argument
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import pickle as pkl
import argparse
import matplotlib.pyplot as plt
import hydroeval as he

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


if __name__ == '__main__':

    os.chdir('/media/jopato/jopato_ssd/PHD/PHD_main/Projects/CoSWAT/Paper_Ch2/codeAndDataAvailability/teran_et_al_2026_coswat_reservoirs') # All relative paths will be based on this ! !
    analysis_folder = f'Scripts/result_analysis' #f'vsc10883/PHD_main/Projects/CoSWAT/Scripts/result_analysis'

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
    sim_outputs_dir = f'{analysis_folder}/reservoirs_coswat/coswat_outputs/streamflow/CoSWATv{version}/{region}'
    sim_mon_pkl_fn  = f'{sim_outputs_dir}/monthly/channelsd_mon.pkl'
    stations_csv    = f'{sim_outputs_dir}/monthly/eval_stations.csv'
    
    
    stat_dir        = f'{analysis_folder}/reservoirs_coswat/validation/stats/CoSWATv{version}/{region}'
    plot_dir        = f'{analysis_folder}/reservoirs_coswat/validation/plots/CoSWATv{version}/{region}/streamflow'
    
    os.makedirs(stat_dir,exist_ok=True)
    os.makedirs(plot_dir,exist_ok=True)
    
    # Read pickles
    with open(sim_mon_pkl_fn,'rb') as simMonPkl:
        simMon_dict = pkl.load(simMonPkl)
    
    
    # Read csv with stations
    stations_df = pd.read_csv(stations_csv,index_col=0)
    
    # Start validation
    kge_s_lst, pbias_s_lst, r_s_lst, alpha_s_lst, beta_s_lst = [], [], [], [], []
    grdc_no_lst  = []
    res_inf_list = []
    
    for id, row in stations_df.iterrows():
        grdc_no   = int(row['grdc_no'])
        river     = row['river']
        stat_name = row['station']
        res_inf   = row['down_res']
        obs_df = simMon_dict[grdc_no]['obs_data']
        sim_df = simMon_dict[grdc_no]['sim_data']
        
        # Format date for safety
        obs_df['date'] = pd.to_datetime(obs_df['date'])
        sim_df['date'] = pd.to_datetime(sim_df['date'])

        cutoff = pd.Timestamp("1970-04-01")
        obs_df = obs_df.loc[obs_df["date"] >= cutoff]
        sim_df = sim_df.loc[sim_df["date"] >= cutoff]
        
        # Merge and correct units
        val_df = sim_df.merge(obs_df, on='date', how='left')
        
        # Calculate stats
        try:
            kge_s, pbias_s, r_s, alpha_s, beta_s = calcSats(val_df,"sim","obs")
        except:
            print(f'Station {grdc_no} was not processed')
            continue
        
        kge_s_lst.append(kge_s) 
        pbias_s_lst.append(pbias_s) 
        r_s_lst.append(r_s) 
        alpha_s_lst.append(alpha_s) 
        beta_s_lst.append(beta_s) 
        grdc_no_lst.append(grdc_no)
        res_inf_list.append(res_inf)
        
        # =============== PLOTS =======================================
        # Mask to match with observed period
        mask_obs = val_df['obs'].notna() & np.isfinite(val_df['obs'])

        # If there is at least some observed data
        if mask_obs.any():
            first_date = val_df.loc[mask_obs, 'date'].iloc[0]
            last_date  = val_df.loc[mask_obs, 'date'].iloc[-1]

            val_df_plot = val_df[(val_df['date'] >= first_date) & (val_df['date'] <= last_date)]
            
        else:
            # fall back to full df if no obs at all
            val_df_plot = val_df
            first_date = val_df['date'].iloc[0]
            last_date  = val_df['date'].iloc[-1]

        val_df = val_df_plot.copy()
        
        if res_inf:
            color = 'indigo'
            
        else:
            color = 'navy'
        # Storage plot
        fig, ax = plt.subplots(figsize=(15,8))
        val_df.plot(x='date', y='obs',ax=ax,color='k',linestyle='--', label = 'Reference')
        val_df.plot(x='date', y='sim',ax=ax,color=color, label = 'Simulated')
        
        ax.set_ylabel("Monthly Streamflow (m3/s)")
        ax.set_title(f"GRDC: {grdc_no}, Station: {stat_name}")
        text = (f"KGE = {kge_s:.2f} PBIAS = {pbias_s:.2f}% R = {r_s:.2f} \u03B1 = {alpha_s:.2f} \u03B2 = {beta_s:.2f}")
        plt.text(0.5, -0.12, text, ha='center', va='top', transform=plt.gca().transAxes)
        plt.tight_layout()
        plt.savefig(f'{plot_dir}/monthly_streamflow_grdc-{grdc_no}.jpg',dpi=400)
        plt.close(fig)
    
    # Create stat summary dataframe
    storageStat_df = pd.DataFrame(data={'grdc_no':grdc_no_lst,'res_inf':res_inf_list, 'kge':kge_s_lst, 'pbias':pbias_s_lst,'R':r_s_lst,'alpha':alpha_s_lst,'beta':beta_s_lst})
    storageStat_df.to_csv(f'{stat_dir}/monthly_streamflow_stats_{region}.csv')