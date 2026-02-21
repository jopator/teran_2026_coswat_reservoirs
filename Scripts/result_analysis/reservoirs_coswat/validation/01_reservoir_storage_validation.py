"""
Author: Jose P. Teran
Date: 2025-11-27

Description:
Validate reservoir/lake storage against reference data (Ensemble, and per dataset, if available)

It will generate csv with summary of KGE, PBIAS and R and reservoir properties; this will be used later for a general analysis across all regions
It will generate boxen plots of stats
It will generate time series plots with simulated vs observed (and shading for observed if more than one source)

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
    sim_outputs_dir = f'{analysis_folder}/reservoirs_coswat/coswat_outputs/reservoir_storage/CoSWATv{version}/{region}'
    obs_pkl_fn      = f'{analysis_folder}/reservoirs_coswat/observations/Reservoir_storage_global_data/processed/data/globalReservoirDataAggregated.pkl'

    sim_mon_pkl_fn  = f'{sim_outputs_dir}/monthly/reservoir_mon.pkl'
    sim_yr_pkl_fn   = f'{sim_outputs_dir}/yearly/reservoir_yr.pkl'
    
    modelSetupPath  = f'CoSWAT-Framework/model-setup/CoSWATv{version}/{region}'
    wshedGisPath    = f'{modelSetupPath}/Watershed'
    reservoir_shp   = f'{wshedGisPath}/Shapes/lakes-grand-ESRI-54003.shp'
    
    
    stat_dir        = f'{analysis_folder}/reservoirs_coswat/validation/stats/CoSWATv{version}/{region}'
    plot_dir        = f'{analysis_folder}/reservoirs_coswat/validation/plots/CoSWATv{version}/{region}'
    
    os.makedirs(stat_dir,exist_ok=True)
    os.makedirs(plot_dir,exist_ok=True)
    
    # Read pickles
    with open(obs_pkl_fn,'rb') as obsPkl:
        obs_dict = pkl.load(obsPkl)
        
    with open(sim_mon_pkl_fn,'rb') as simMonPkl:
        simMon_dict = pkl.load(simMonPkl)
        
    with open(sim_yr_pkl_fn,'rb') as simYrPkl:
        simYr_dict = pkl.load(simYrPkl)
        
    
    # Filter observations that correspond to simulated reservoirs
    # ====================
    grand_id_list = []
    grand_dict = {}
    
    for k,item in simMon_dict.items(): 
        grand_id_list.append(int(item['grand_id']))
        grand_dict[item['grand_id']] = int(k)
    
    grand_id_list = [x for x in grand_id_list if x!=0]
    obsRegion_dict = {k: obs_dict[k] for k in grand_id_list if k in obs_dict}
    
    
    
    # Process Monthly
    # ====================
    res_gdf = gpd.read_file(reservoir_shp)
    
    kge_s_lst, pbias_s_lst, r_s_lst, alpha_s_lst, beta_s_lst = [], [], [], [], []
    kge_i_lst, pbias_i_lst, r_i_lst, alpha_i_lst, beta_i_lst = [], [], [], [], []
    kge_o_lst, pbias_o_lst, r_o_lst, alpha_o_lst, beta_o_lst = [], [], [], [], []
    
    grand_id_o_lst = []
    grand_id_i_lst = []
    grand_id_s_lst = []
    
    for grand_id in grand_id_list:
        if grand_id not in obsRegion_dict:
            continue
        obs_df      = obsRegion_dict[grand_id]
        coswatId    = grand_dict[grand_id]
        hylak_id    = simMon_dict[coswatId]['Hylak_id']
        sim_df      = simMon_dict[coswatId]['coswatOut']
        lake_name   = res_gdf.loc[res_gdf['Hylak_id'] == hylak_id, 'Lake_name'].iloc[0]
        dam_name    = res_gdf.loc[res_gdf['Hylak_id'] == hylak_id, 'DAM_NAME'].iloc[0]
        
        # Format date for safety
        obs_df['date'] = pd.to_datetime(obs_df['date'])  + pd.offsets.MonthEnd(0)
        sim_df['date'] = pd.to_datetime(sim_df['date'])
        
        # Merge and correct units
        val_df = sim_df.merge(obs_df, on='date', how='left')
        sec = val_df['date'].dt.days_in_month * 24 * 3600
        
        val_df['flo_in_m3s']    = val_df['flo_in'] / sec
        val_df['flo_stor_km3']  = val_df['flo_stor']/1000000000
        val_df['flo_out_m3s']   = val_df['flo_out'] / sec
        
        val_df['storage_mean'] = pd.to_numeric(val_df['storage_mean'], errors='coerce')
        val_df['storage_min']  = pd.to_numeric(val_df['storage_min'],  errors='coerce')
        val_df['storage_max']  = pd.to_numeric(val_df['storage_max'],  errors='coerce')
        
        # Filter out observations that have a relative error > 50 % comparing long term storage to hydrolakes/granD storage
        lta_storage_obs   = val_df['storage_mean'].mean()
        lta_storage_hylak = float(res_gdf.loc[res_gdf['Hylak_id'] == hylak_id, 'smax'].iloc[0])/1000000000
        
        rel_err = 100 * abs(lta_storage_hylak-lta_storage_obs)/lta_storage_hylak
        
        if abs(rel_err) > 60:
            continue
        
        #Skip lake victoria
        if hylak_id == 16:
            continue
        
        # Calculate stats
        try:
            kge_s, pbias_s, r_s, alpha_s, beta_s = calcSats(val_df,"flo_stor_km3","storage_mean")
        except:
            continue
        
        if kge_s < -5.0 or pbias_s < -200 or pbias_s > 200: #Likely an error here, plots show an unreal difference between sim and obs
            continue
        
        kge_s_lst.append(kge_s) 
        pbias_s_lst.append(pbias_s) 
        r_s_lst.append(r_s) 
        alpha_s_lst.append(alpha_s) 
        beta_s_lst.append(beta_s) 
        grand_id_s_lst.append(grand_id)
        
        if not val_df['inflow_mean'].isna().all():
            kge_i, pbias_i, r_i, alpha_i, beta_i = calcSats(val_df,"flo_in_m3s","inflow_mean")
            kge_i_lst.append(kge_i) 
            pbias_i_lst.append(pbias_i) 
            r_i_lst.append(r_i) 
            alpha_i_lst.append(alpha_i) 
            beta_i_lst.append(beta_i)
            grand_id_i_lst.append(grand_id) 
        
        if not val_df['outflow_mean'].isna().all():   
            kge_o, pbias_o, r_o, alpha_o, beta_o = calcSats(val_df,"flo_out_m3s","outflow_mean")
            kge_o_lst.append(kge_o) 
            pbias_o_lst.append(pbias_o) 
            r_o_lst.append(r_o) 
            alpha_o_lst.append(alpha_o) 
            beta_o_lst.append(beta_o)
            grand_id_o_lst.append(grand_id)  
        
        
        # =============== PLOTS =======================================
        # Mask to match with observed period
        mask_obs = val_df['storage_mean'].notna() & np.isfinite(val_df['storage_mean'])

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
        
        
        # Storage plot
        fig, ax = plt.subplots(figsize=(15,8))
        val_df.plot(x='date', y='storage_mean',ax=ax,color='k',linestyle='--', label = 'Reference')
        val_df.plot(x='date', y='flo_stor_km3',ax=ax,color='navy', label = 'Simulated')
        
        x = val_df['date'].values
        ymin = val_df['storage_min'].to_numpy(dtype=float)
        ymax = val_df['storage_max'].to_numpy(dtype=float)
        
        mask = np.isfinite(ymin) & np.isfinite(ymax)
        ax.fill_between(
            x,
            ymin,
            ymax,
            where=mask,
            alpha=0.2,)
        
        ax.set_ylabel("Monthly Storage (km3)")
        
        if lake_name == None:
            titlename = dam_name
        elif lake_name != None:
            titlename = lake_name
        else:
            titlename = "No name"
        
            
        ax.set_title(f"{titlename}")
        text = (f"KGE = {kge_s:.2f} PBIAS = {pbias_s:.2f}% R = {r_s:.2f} \u03B1 = {alpha_s:.2f} \u03B2 = {beta_s:.2f}")
        plt.text(0.5, -0.12, text, ha='center', va='top', transform=plt.gca().transAxes)
        plt.tight_layout()
        plt.savefig(f'{plot_dir}/monthly_storage_validation_grand-{grand_id}_hylak-{hylak_id}.jpg',dpi=400)
        plt.close(fig)
        
        # Inflow plot
        
        if not val_df['inflow_mean'].isna().all():    
            fig, ax = plt.subplots(figsize=(15,8))
            val_df.plot(x='date', y='inflow_mean',ax=ax,color='k',linestyle='--', label = 'Reference')
            val_df.plot(x='date', y='flo_in_m3s',ax=ax,color='navy',linewidth = 0.5, label = 'Simulated')
            
            x = val_df['date'].values
            ymin = val_df['inflow_min'].to_numpy(dtype=float)
            ymax = val_df['inflow_max'].to_numpy(dtype=float)
            
            mask = np.isfinite(ymin) & np.isfinite(ymax)
            ax.fill_between(
                x,
                ymin,
                ymax,
                where=mask,
                alpha=0.2,)
            
            ax.set_ylabel("Mean monthly inflow (m3/s)")    
            ax.set_title(f"{titlename}")
            text = (f"KGE = {kge_i:.2f} PBIAS = {pbias_i:.2f}% R = {r_i:.2f} \u03B1 = {alpha_i:.2f} \u03B2 = {beta_i:.2f}")
            plt.text(0.5, -0.12, text, ha='center', va='top', transform=plt.gca().transAxes)
            
            plt.tight_layout()
            
            plt.savefig(f'{plot_dir}/monthly_inflow_validation_grand-{grand_id}_hylak-{hylak_id}.jpg',dpi=400)
            plt.close(fig)
        # Outflow plot
        if not val_df['outflow_mean'].isna().all():    
            fig, ax = plt.subplots(figsize=(15,8))
            val_df.plot(x='date', y='outflow_mean',ax=ax,color='k',linestyle='--', label = 'Reference')
            val_df.plot(x='date', y='flo_out_m3s',ax=ax,color='navy',linewidth = 0.5, label = 'Simulated')
            
            x = val_df['date'].values
            ymin = val_df['outflow_min'].to_numpy(dtype=float)
            ymax = val_df['outflow_max'].to_numpy(dtype=float)
            
            mask = np.isfinite(ymin) & np.isfinite(ymax)
            ax.fill_between(
                x,
                ymin,
                ymax,
                where=mask,
                alpha=0.2,)
            
            ax.set_ylabel("Mean monthly outflow (m3/s)")    
            ax.set_title(f"{titlename}")
            text = (f"KGE = {kge_o:.2f} PBIAS = {pbias_o:.2f}% R = {r_o:.2f} \u03B1 = {alpha_o:.2f} \u03B2 = {beta_o:.2f}")
            plt.text(0.5, -0.12, text, ha='center', va='top', transform=plt.gca().transAxes)
            
            plt.tight_layout()
            
            plt.savefig(f'{plot_dir}/monthly_outflow_validation_grand-{grand_id}_hylak-{hylak_id}.jpg',dpi=400)
            plt.close(fig)
        
        # ======================================================
        
    # Create stat summary dataframe
    storageStat_df = pd.DataFrame(data={'grand_id':grand_id_s_lst, 'kge':kge_s_lst, 'pbias':pbias_s_lst,'R':r_s_lst,'alpha':alpha_s_lst,'beta':beta_s_lst})
    storageStat_df.to_csv(f'{stat_dir}/monthly_storage_stats_{region}.csv')
    
    if kge_i_lst:
        inflowStat_df = pd.DataFrame(data={'grand_id':grand_id_i_lst, 'kge':kge_i_lst, 'pbias':pbias_i_lst,'R':r_i_lst,'alpha':alpha_i_lst,'beta':beta_i_lst})
        inflowStat_df.to_csv(f'{stat_dir}/monthly_inflow_stats_{region}.csv')
    
    if kge_o_lst:
        outflowStat_df = pd.DataFrame(data={'grand_id':grand_id_o_lst, 'kge':kge_o_lst, 'pbias':pbias_o_lst,'R':r_o_lst,'alpha':alpha_o_lst,'beta':beta_o_lst})
        outflowStat_df.to_csv(f'{stat_dir}/monthly_outflow_stats_{region}.csv')