"""
Author: Jose P. Teran
Date: 2026-02-10

Description:
Plot selected reservoir time series
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
from scipy.stats import pearsonr
import pickle as pkl
from matplotlib.lines import Line2D

os.chdir('/media/jopato/jopato_ssd/PHD/PHD_main/Projects/CoSWAT/Paper_Ch2/codeAndDataAvailability/teran_et_al_2026_coswat_reservoirs')#'/data/brussel/vo/000/bvo00033') # All relative paths will be based on this ! !

# Paths settings and file names
# ====================
version         = "1.5.0"

analysis_folder = f'Scripts/result_analysis'#f'vsc10883/PHD_main/Projects/CoSWAT/Scripts/result_analysis'
modelSetupPath  = 'CoSWAT-Framework/model-setup/CoSWATv{version}/{region}'
wshedGisPath    = f'{modelSetupPath}/Watershed'
stat_dir        = f'{analysis_folder}/reservoirs_coswat/validation/stats/CoSWATv{version}'
out_dir         = f'{analysis_folder}/reservoirs_coswat/validation/analysis'
stat_plots_dir  = f'{analysis_folder}/reservoirs_coswat/validation/analysis/plots'

os.makedirs(out_dir,exist_ok=True)
os.makedirs(stat_plots_dir,exist_ok=True)

monInflowCsvFn  = 'monthly_inflow_stats_{region}.csv'
monStorageCsvFn = 'monthly_storage_stats_{region}.csv'
monOutflowCsvFn = 'monthly_outflow_stats_{region}.csv'



regions = [
           'africa-nile',
           'africa-orange',
           'america-bravo',
           'america-colorado',
           'america-mississippi',
           'america-parana',
           'asia-mekong',
           'europe-central',
           'europe-west'
           ]

# def month_climatology(df,column):
#     df = df.copy()
#     df["month"] = df["date"].copy().dt.month
#     return df.groupby("month")[column].agg(["mean", "min", "max"]).reset_index()

def month_climatology(df, column):
    df = df.copy()
    df["month"] = df["date"].dt.month

    g = df.groupby("month")[column]
    out = g.agg(["mean", "median", "std"]).reset_index()

    out["min"]  = out["median"] - 1.5 * out["std"]
    out["max"] = out["median"] + 1.5 * out["std"]

    return out[["month", "mean", "min", "max"]]

# Read observations dictionary
obs_fn = f'{analysis_folder}/reservoirs_coswat/observations/Reservoir_storage_global_data/processed/data/globalReservoirDataAggregated.pkl'

with open(obs_fn,'rb') as f:
    obs_dict = pkl.load(f)


#Read simulation performance tables
df_list = []

df_list_inf = []
df_list_out = []
# Read inflow data
for region in regions:
    try:
        df = pd.read_csv(f'{stat_dir}/{region}/{monInflowCsvFn.format(region=region)}',index_col=0)
        df_list_inf.append(df)
        
        df_out = pd.read_csv(f'{stat_dir}/{region}/{monOutflowCsvFn.format(region=region)}',index_col=0)
        df_list_out.append(df_out)
    except:
        continue
    
df_inflow_perf  = pd.concat(df_list_inf)
df_outflow_perf = pd.concat(df_list_out)


for region in regions:
    df = pd.read_csv(f'{stat_dir}/{region}/{monStorageCsvFn.format(region=region)}',index_col=0)
    df_list.append(df)
    
df_storage = pd.concat(df_list)
df_storage = df_storage.reset_index(drop=True)


# Read time series
grand_ids  = [185,870,4478]

sim_df_dict = {}
obs_df_dict = {}

for region in regions:
    pkl_fn = f'{analysis_folder}/reservoirs_coswat/coswat_outputs/reservoir_storage/CoSWATv1.5.0/{region}/monthly/reservoir_mon.pkl'
    with open(pkl_fn,'rb') as f:
        data_dict = pkl.load(f)
        
        
        for reservoir_id,item in data_dict.items():
            for grand_id in grand_ids:
                try:
                    caught = item['grand_id'] == grand_id
                    
                    if caught:
                        sim_df_dict[grand_id] = item['coswatOut']
                        obs_df_dict[grand_id] = obs_dict[grand_id]
                
                except:
                    dumy=0


sim_inflow_clim_dict = {}
sim_outflow_clim_dict = {}
obs_inflow_clim_dict = {}
obs_outflow_clim_dict = {}

for grand_id in grand_ids:
    sim_df_dict[grand_id]['storage_km3'] = sim_df_dict[grand_id].copy()['flo_stor']*1e-9
    sim_df_dict[grand_id]['inflow_m3s']  = sim_df_dict[grand_id].copy()['flo_in']/86400/30
    sim_df_dict[grand_id]['outflow_m3s']  = sim_df_dict[grand_id].copy()['flo_out']/86400/30
    


    obs_df_dict[grand_id].loc[:, "inflow_mean"]  = obs_df_dict[grand_id]["inflow_mean"].clip(lower=0)
    obs_df_dict[grand_id].loc[:, "outflow_mean"] = obs_df_dict[grand_id]["outflow_mean"].clip(lower=0)


    df_sim_clim_inflow = month_climatology(sim_df_dict[grand_id],'inflow_m3s')
    df_obs_clim_inflow = month_climatology(obs_df_dict[grand_id],'inflow_mean')

    df_sim_clim_outflow = month_climatology(sim_df_dict[grand_id],'outflow_m3s')
    df_obs_clim_outflow = month_climatology(obs_df_dict[grand_id],'outflow_mean')

    date_min = max(sim_df_dict[grand_id]['date'].min(),obs_df_dict[grand_id]['date'].min())
    date_max = min(sim_df_dict[grand_id]['date'].max(),obs_df_dict[grand_id]['date'].max())

    
    sim_df_dict[grand_id] = sim_df_dict[grand_id][(sim_df_dict[grand_id]['date']>=date_min) & (sim_df_dict[grand_id]['date']<=date_max)].copy()
    obs_df_dict[grand_id] = obs_df_dict[grand_id][(obs_df_dict[grand_id]['date']>=date_min) & (obs_df_dict[grand_id]['date']<=date_max)].copy()

    sim_inflow_clim_dict[grand_id]  = df_sim_clim_inflow
    sim_outflow_clim_dict[grand_id] = df_sim_clim_outflow
    obs_inflow_clim_dict[grand_id]  = df_obs_clim_inflow
    obs_outflow_clim_dict[grand_id] = df_obs_clim_outflow


fig, axs = plt.subplots(3,3,figsize=(20,10),width_ratios=[1,2,1])
fig.subplots_adjust(hspace=0.6,wspace=0.3) 
month_labels = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
month_labels2 = ["Jan","Mar","May","Jul","Sep","Nov"]

color_obs = 'dimgray'
color_sims = 'firebrick'

for i, grand_id in enumerate(grand_ids):
    ax1 = axs[i,0]
    ax2 = axs[i,1]
    ax3 = axs[i,2]

    
    
    # Storage plot
    #==============
    kge = float(df_storage.loc[df_storage["grand_id"] == grand_id, "kge"].reset_index(drop=True).iloc[0])
    pbias = float(df_storage.loc[df_storage["grand_id"] == grand_id, "pbias"].reset_index(drop=True).iloc[0])
    r = float(df_storage.loc[df_storage["grand_id"] == grand_id, "R"].reset_index(drop=True).iloc[0])
    

    ax2.text(0.02, 0.12, f"KGE = {kge:.2f}, PBIAS = {pbias:.1f} (%), r = {r:.1f}",transform=ax2.transAxes,ha="left",va="top",bbox=dict(facecolor="white", edgecolor="none", alpha=0.8))

    df_sim = sim_df_dict[grand_id]
    df_obs = obs_df_dict[grand_id]
    
    df_sim_clim_inflow = sim_inflow_clim_dict[grand_id]
    df_obs_clim_inflow = obs_inflow_clim_dict[grand_id]
    
    df_sim_clim_outflow = sim_outflow_clim_dict[grand_id]
    df_obs_clim_outflow = obs_outflow_clim_dict[grand_id]
    
    df_sim.plot(x='date',y='storage_km3',ax=ax2,label='Simulated',color=color_sims,linestyle='--',linewidth=1.2)
    df_obs.plot(x='date',y='storage_mean',ax=ax2,label='Observed',color=color_obs,linewidth=1.5)
    ax2.set_xlabel("Date")
    ax2.set_ylabel("Monthly Storage $(km^{3})$")
    ax2.set_ylim(0,max(df_obs['storage_mean'].max(),df_sim['storage_km3'].max())*1.05)
    
    
    # Inflow plot
    #==============
    
    kge     = float(df_inflow_perf.loc[df_inflow_perf["grand_id"] == grand_id, "kge"].reset_index(drop=True).iloc[0])
    pbias   = float(df_inflow_perf.loc[df_inflow_perf["grand_id"] == grand_id, "pbias"].reset_index(drop=True).iloc[0])
    r       = float(df_inflow_perf.loc[df_inflow_perf["grand_id"] == grand_id, "R"].reset_index(drop=True).iloc[0])
    
    ax1.text(0.58, 0.95, f"KGE = {kge:.2f}\nPBIAS = {pbias:.1f} (%)\nr = {r:.1f}",transform=ax1.transAxes,ha="left",va="top",
             bbox=dict(facecolor="white", edgecolor="none", alpha=0.4))
    
    df_sim_clim_inflow.plot(x='month',y='mean',ax=ax1,label='Simulated',color=color_sims,linestyle='--',linewidth=1.2)
    df_obs_clim_inflow.plot(x='month',y='mean',ax=ax1,label='Observed',color=color_obs,linewidth=2)
    ax1.fill_between(df_sim_clim_inflow['month'],df_sim_clim_inflow['min'],df_sim_clim_inflow['max'],alpha=0.1,color=color_sims)

    
    for c in ["month", "min", "max"]:
        df_obs_clim_inflow[c] = pd.to_numeric(df_obs_clim_inflow[c], errors="coerce")
        
    ax1.fill_between(df_obs_clim_inflow['month'],df_obs_clim_inflow['min'],df_obs_clim_inflow['max'],alpha=0.3,color=color_obs)

    
    ax1.set_ylim(0,max(df_obs_clim_inflow['max'].max(),df_sim_clim_inflow['max'].max())*1.05)
    ax1.set_xlim(1,12)
    ax1.set_xlabel("Month of the year")
    ax1.set_ylabel("Monthly Inflow $(m^{3} s^{-1})$")
    ax1.grid(axis='x',linewidth = 0.4)
    ax1.set_xticks(range(1,13,2))
    ax1.set_xticklabels(month_labels2)
    if i == 0:
        ax1.ticklabel_format(axis='y', style='sci', scilimits=(2,2),useMathText=True)
    
    else:
        ax1.ticklabel_format(axis='y', style='sci', scilimits=(3,3),useMathText=True)
    
    # Outflow plot
    #==============
    
    kge     = float(df_outflow_perf.loc[df_outflow_perf["grand_id"] == grand_id, "kge"].reset_index(drop=True).iloc[0])
    pbias   = float(df_outflow_perf.loc[df_outflow_perf["grand_id"] == grand_id, "pbias"].reset_index(drop=True).iloc[0])
    r       = float(df_outflow_perf.loc[df_outflow_perf["grand_id"] == grand_id, "R"].reset_index(drop=True).iloc[0])
    
    ax3.text(0.57, 0.95, f"KGE = {kge:.2f}\nPBIAS = {pbias:.1f} (%)\nr = {r:.1f}",transform=ax3.transAxes,ha="left",va="top",
             bbox=dict(facecolor="white", edgecolor="none", alpha=0.4))    
    
    df_sim_clim_outflow.plot(x='month',y='mean',ax=ax3,label='Simulated',color=color_sims,linestyle='--',linewidth=1.2)
    df_obs_clim_outflow.plot(x='month',y='mean',ax=ax3,label='Observed',color=color_obs,linewidth=2)
    ax3.fill_between(df_sim_clim_outflow['month'],df_sim_clim_outflow['min'],df_sim_clim_outflow['max'],alpha=0.1,color=color_sims)
    
    for c in ["month", "min", "max"]:
        df_obs_clim_outflow[c] = pd.to_numeric(df_obs_clim_outflow[c], errors="coerce")
        
    ax3.fill_between(df_obs_clim_outflow['month'],df_obs_clim_outflow['min'],df_obs_clim_outflow['max'],alpha=0.3,color=color_obs)      
    ax3.set_ylim(0,max(df_obs_clim_outflow['max'].max(),df_sim_clim_outflow['max'].max())*1.05)
    ax3.set_xlim(1,12)
    ax3.set_xlabel("Month of the year")
    ax3.set_ylabel("Monthly Outflow $(m^{3} s^{-1})$")
    ax3.grid(axis='x',linewidth = 0.4)
    ax3.set_xticks(range(1,13,2))
    ax3.set_xticklabels(month_labels2)
    
    if i == 0:
        ax3.ticklabel_format(axis='y', style='sci', scilimits=(2,2),useMathText=True)
    
    else:
        ax3.ticklabel_format(axis='y', style='sci', scilimits=(3,3),useMathText=True)
    
    ax1.get_legend().remove()
    ax2.get_legend().remove()
    ax3.get_legend().remove()
    
    
plt.text(-0.25, 4.40, "a)  Berryessa Lake (Monticello Dam) - Hydroelectricity - Sacramento River Basin",fontsize=12,fontweight='bold',transform=ax1.transAxes)
plt.text(-0.25, 2.80, "b)  Lake Oahe (Oahe Dam) - Flood Control - Mississippi River System",fontsize=12,fontweight='bold',transform=ax1.transAxes)
plt.text(-0.25, 1.20, "c)  Nasser Lake (Aswan Dam) - Irrigation - Nile River Basin",fontsize=12,fontweight='bold',transform=ax1.transAxes)

sim_handle = Line2D([0], [0],color=color_sims,linestyle="--",linewidth=1.5,label="Simulated")
obs_handle = Line2D([0], [0],color=color_obs,linestyle="-",linewidth=1.5,label="Reference")
sim_handle.set_linewidth(2.5)
obs_handle.set_linewidth(2.5)
ax2.legend(handles=[sim_handle, obs_handle],loc="upper center",bbox_to_anchor=(0.5, -0.25),ncol=2,fontsize=12,frameon=False)

plt.savefig(f'{stat_plots_dir}/figure6.jpeg',dpi=300,bbox_inches='tight')
    