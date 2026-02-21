"""
Author: Jose P. Teran
Date: 2025-11-27

Description:

- Analysis of processed outputs for reservoir storage, inflow, and outflow performance
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
from scipy.stats import pearsonr

os.chdir('/media/jopato/jopato_ssd/PHD/PHD_main/Projects/CoSWAT/Paper_Ch2/codeAndDataAvailability/teran_et_al_2026_coswat_reservoirs')#'/data/brussel/vo/000/bvo00033') # All relative paths will be based on this ! !


# Functions
def fmt_lat(lat, pos=None):
    if lat > 0:
        return f"{abs(lat):.0f}°N"
    elif lat < 0:
        return f"{abs(lat):.0f}°S"
    else:
        return "0°"

def fmt_lon(lon, pos=None):
    if lon > 0:
        return f"{abs(lon):.0f}°E"
    elif lon < 0:
        return f"{abs(lon):.0f}°W"
    else:
        return "0°"
    

def write_stats_summary(log, df, label):
    num_res = int(len(df))

    num_kge_min  = int((df["kge"] > -0.41).sum())
    num_kge_ok   = int((df["kge"] > 0.0).sum())
    num_kge_good = int((df["kge"] > 0.3).sum())

    num_r_min    = int((df["R"] > 0.0).sum())
    num_r_ok     = int((df["R"] > 0.2).sum())
    num_r_good   = int((df["R"] > 0.4).sum())

    num_pbias_ok   = int((df["pbias"].abs() <= 50).sum())
    num_pbias_good = int((df["pbias"].abs() <= 25).sum())

    kge_min, kge_max, kge_avg, kge_std = df["kge"].min(), df["kge"].max(), df["kge"].mean(), df["kge"].std()
    r_min,   r_max,   r_avg,   r_std   = df["R"].min(),   df["R"].max(),   df["R"].mean(),   df["R"].std()
    pb_min,  pb_max,  pb_avg,  pb_std  = df["pbias"].min(), df["pbias"].max(), df["pbias"].mean(), df["pbias"].std()

    log.write(f"=== Summary of stats for monthly reservoir {label} ====\n")
    log.write(f"A total of {num_res} were evaluated\n\n")

    log.write("=== Overall statistics ===\n")
    log.write(f"KGE   -> min: {kge_min:.2f}, max: {kge_max:.2f}, mean: {kge_avg:.2f}, std: {kge_std:.2f}\n")
    log.write(f"R     -> min: {r_min:.2f},   max: {r_max:.2f},   mean: {r_avg:.2f},   std: {r_std:.2f}\n")
    log.write(f"PBIAS -> min: {pb_min:.2f}%, max: {pb_max:.2f}%, mean: {pb_avg:.2f}%, std: {pb_std:.2f}%\n\n")

    log.write("=== Performance threshold distribution ===\n")

    p_kge_min  = 100.0 * num_kge_min  / num_res if num_res else 0.0
    p_kge_ok   = 100.0 * num_kge_ok   / num_res if num_res else 0.0
    p_kge_good = 100.0 * num_kge_good / num_res if num_res else 0.0
    p_r_min    = 100.0 * num_r_min    / num_res if num_res else 0.0
    p_r_ok     = 100.0 * num_r_ok     / num_res if num_res else 0.0
    p_r_good   = 100.0 * num_r_good   / num_res if num_res else 0.0
    p_pbias_ok   = 100.0 * num_pbias_ok   / num_res if num_res else 0.0
    p_pbias_good = 100.0 * num_pbias_good / num_res if num_res else 0.0

    log.write(f"KGE > -0.41 : {num_kge_min} ({p_kge_min:.1f}%)\n")
    log.write(f"KGE > 0.00  : {num_kge_ok} ({p_kge_ok:.1f}%)\n")
    log.write(f"KGE > 0.30  : {num_kge_good} ({p_kge_good:.1f}%)\n")
    log.write(f"R > 0.00    : {num_r_min} ({p_r_min:.1f}%)\n")
    log.write(f"R > 0.20    : {num_r_ok} ({p_r_ok:.1f}%)\n")
    log.write(f"R > 0.40    : {num_r_good} ({p_r_good:.1f}%)\n")
    log.write(f"|PBIAS| ≤ 50 : {num_pbias_ok} ({p_pbias_ok:.1f}%)\n")
    log.write(f"|PBIAS| ≤ 25 : {num_pbias_good} ({p_pbias_good:.1f}%)\n\n")



# Paths settigs and file names
# =======================================================================================================
version         = "1.5.0"

analysis_folder = f'Scripts/result_analysis'#f'vsc10883/PHD_main/Projects/CoSWAT/Scripts/result_analysis'
modelSetupPath  = 'CoSWAT-Framework/model-setup/CoSWATv{version}/{region}'
wshedGisPath    = f'{modelSetupPath}/Watershed'
reservoir_shp   = f'{wshedGisPath}/Shapes/lakes-grand-ESRI-54003.shp'
global_shapes   = f'Scripts/result_analysis/reservoirs_coswat/validation/shapefiles'
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


colors = {
    "kge":  "dimgray",
    "R":    "darkslategray",
    "pbias": "midnightblue",
    "alpha": "#4daab7",
    "beta": "#becc27"
}



logOut = f'{out_dir}/outLogStorage.txt'




#===============================================================================
# General performance distribution plot
#===============================================================================
#Aggregate
df_list          = []
df_inf_list      = []
df_outf_list     = []


for region in regions:
    df      = pd.read_csv(f'{stat_dir}/{region}/{monStorageCsvFn.format(region=region)}',index_col=0)   # Storage
    df_list.append(df)
    
    try:
        df_inf  = pd.read_csv(f'{stat_dir}/{region}/{monInflowCsvFn.format(region=region)}',index_col=0)    # Inflow
        df_inf_list.append(df_inf)
    except:
        print(f'{region} does not have inflow data')
        continue
    
    try:
        df_outf = pd.read_csv(f'{stat_dir}/{region}/{monOutflowCsvFn.format(region=region)}',index_col=0)   # Outflow
        df_outf_list.append(df_outf)
    except:
        print(f'{region} does not have inflow data')
        continue
    
df_gral = pd.concat(df_list)
df_gral = df_gral.reset_index(drop=True)

df_gral_inf = pd.concat(df_inf_list)
df_gral_inf = df_gral_inf.reset_index(drop=True)

df_gral_outf = pd.concat(df_outf_list)
df_gral_outf = df_gral_outf.reset_index(drop=True)

df_gral      = df_gral[df_gral['kge']>-5]
df_gral_inf  = df_gral_inf[df_gral_inf['kge']>-5]
df_gral_outf = df_gral_outf[df_gral_outf['kge']>-5]

df_long = df_gral.melt(id_vars="grand_id", value_vars=["kge", "R", "pbias","alpha","beta"], var_name="metric",value_name="value")
df_long_inf = df_gral_inf.melt(id_vars="grand_id",value_vars=["kge", "R", "pbias","alpha","beta"],var_name="metric",value_name="value")
df_long_outf = df_gral_outf.melt(id_vars="grand_id",value_vars=["kge", "R", "pbias","alpha","beta"],var_name="metric",value_name="value")
 


# Distribution
stats = ["kge", "pbias", "R"]

data_sets = [("Storage", df_long),("Inflow", df_long_inf),("Outflow", df_long_outf),]

fig, axs = plt.subplots(3, 3, figsize=(20, 10),sharey=True)
axs = axs.flatten()

for i, (label, df_sel) in enumerate(data_sets):
    for j, stat in enumerate(stats):

        ax = axs[i*3 + j]
        ax.set_xlabel(stat.upper())
        
        if label == "Inflow":
            sns.histplot(data=df_sel[df_sel["metric"] == stat],x="value",kde=True,stat="percent",bins=15,color=colors[stat],ax=ax)
        elif label == "Outflow":
            sns.histplot(data=df_sel[df_sel["metric"] == stat],x="value",kde=True,stat="percent",bins=25,color=colors[stat],ax=ax)           
        else:
            sns.histplot(data=df_sel[df_sel["metric"] == stat],x="value",kde=True,stat="percent",color=colors[stat],ax=ax)
            
            ax.set_ylim(0,25)

        if stat == 'kge':
            if label == 'Storage':
                ax.set_xlim(-2,0.8)
            else:
                ax.set_xlim(-2,0.8)

            ax.axvline(-0.41, color="k", linestyle=":", linewidth=1.5,label="Minimum")
            ax.axvline(0.0, color="k", linestyle="--", linewidth=2,label="Satisfactory")
            ax.axvline(0.3, color="k", linestyle="dashdot", linewidth=1,label="Good")

        if stat == 'R':
            ax.set_xlim(-0.75,1)
            ax.axvline(0.0, color="k", linestyle="--", linewidth=2,label="Satisfactory")
            ax.axvline(0.4, color="k", linestyle="dashdot", linewidth=1,label="Good")
            ax.set_xlabel("r")

        if stat == 'pbias':
            if label == 'Storage':
                ax.set_xlim(-120,100)
            else:
                ax.set_xlim(-160,120)

            ax.axvline(-25, color="k", linestyle="dashdot", linewidth=1.5,label="Good")
            ax.axvline(25, color="k", linestyle="dashdot", linewidth=1.5)
            ax.axvline(-50, color="k", linestyle="--", linewidth=2,label="Satisfactory")
            ax.axvline(50, color="k", linestyle="--", linewidth=2)

        # Y-axis label only on first column
        if j == 0:
            ax.set_ylabel("% of water bodies",fontsize=12)

        # Panel letters per row
        if j == 0:
            ax.text(-0.1, 1.15,  f"{chr(97+i)}) {label}",transform=ax.transAxes,ha="left", va="top",fontsize=12,fontweight="bold")

        ax.legend()


fig.subplots_adjust(hspace=0.40)

row_titles = ["Storage", "Inflow", "Outflow"]
plt.savefig(f'{stat_plots_dir}/figure5.jpg',bbox_inches="tight")
plt.close()


#===============================================================================
# Summary report
#===============================================================================
with open(logOut, "w") as log:
    write_stats_summary(log, df_gral,      "storage")
    write_stats_summary(log, df_gral_inf,  "inflow")
    write_stats_summary(log, df_gral_outf, "outflow")


#===============================================================================
# Performance spatial distribution
#===============================================================================

# Read all reservoir shapes and place into one geodataframe
gdf_list    = []

for region in regions:
    gdf = gpd.read_file(f'{reservoir_shp.format(region=region,version=version)}')
    gdf_list.append(gdf)

gdf_gral = pd.concat(gdf_list)
gdf_gral = gdf_gral.reset_index(drop=True)

gdf_gral = gdf_gral[['LakeId', 'Hylak_id','smax', 'pvol', 'evol', 'parea', 'earea','Country', 'Continent','Lake_type', 'Depth_avg',
       'Res_time', 'Grand_id','elev_masl','maxDepth','MAIN_USE','geometry']]
gdf_gral = gdf_gral.rename(columns={'Grand_id':'grand_id'})

gdf_eval = gdf_gral[gdf_gral['grand_id'].isin(df_gral['grand_id'].to_list())].reset_index(drop=True)
gdf_eval = gdf_eval.copy()
gdf_eval['geometry'] = gdf_eval.geometry.centroid # Convert to point
gdf_eval = pd.merge(gdf_eval,df_gral,on='grand_id',how='left')

# Continents and major basins - rivers - lakes
url = "https://github.com/nvkelso/natural-earth-vector/raw/master/110m_physical/ne_110m_land.shp"
world = gpd.read_file(url)

basins_file = f"{global_shapes}/basins_analysis.gpkg"
basins_gdf = gpd.read_file(basins_file).to_crs("EPSG:4326")

rivers_file = f"{global_shapes}/major_rivers_analysis.gpkg"
rivers_gdf = gpd.read_file(rivers_file,layer='output').to_crs("EPSG:4326")

lakes_file = f"{global_shapes}/all_lakes.gpkg"
lakes_gdf = gpd.read_file(lakes_file).to_crs("EPSG:4326")


# Plot configuration
lat_limits = (-40,58)
lon_limits = (-128,112)
markersize = 15
marker_cat = "s"
marker_kge = "o"
continent_alpha = 0.05
continent_color = "gainsboro"
continent_width = 0.6
regions_alpha = 0.2
regions_width = 0.05

kge_colors = ["#fc8d59","yellowgreen","forestgreen"]
kge_bounds = [-0.42, 0.0, 0.2, 0.4]

river_color = 'navy'
river_alpha = 0.8

lake_color = 'cadetblue'
lake_alpha = 0.8


# Plot execution
# Plot KGE
gdf_eval_plot = gdf_eval.to_crs(world.crs)
gdf_eval_plot['absPbias'] = abs(gdf_eval_plot['pbias'])

fig, axs = plt.subplots(3,1,figsize=(12, 14))

stats = {0:'kge',1:'absPbias',2:'R'}
statName = {0:'KGE',1:'|PBIAS| (%)',2:'r'}


kge_bounds = [-0.4,-0.2,0.0,0.2,0.4]
r_bounds = [0.0,0.1,0.2,0.3,0.4]
pbias_bounds = [0,12.5,25,50]

bounds = {0:kge_bounds,1:pbias_bounds,2:r_bounds}

vmax   = {0:0.4,1:75,2:0.4}
vmin   = {0:-0.4,1:0,2:0}

cmaps  = {0:'RdYlGn',1:'RdYlGn_r',2:'RdYlGn'}
extend = {0:'both',1:'max',2:'both'}
for i,ax in enumerate(axs):
    # Create a ListedColormap from RdYlGn with len(bounds)-1 bins
    cmap_stat = cmaps[i]
    cmap = plt.get_cmap(cmap_stat, len(bounds[i]) +1)
    norm = mcolors.BoundaryNorm(boundaries=bounds[i], ncolors=cmap.N, extend=extend[i])

    world.plot(ax=ax,color=continent_color, alpha=continent_alpha, aspect='auto')
    world.boundary.plot(ax=ax, linewidth=continent_width, color="black", aspect='auto')
    basins_gdf.plot(ax=ax,color='gray',alpha=regions_alpha)
    basins_gdf.boundary.plot(ax=ax,color='k',linewidth=regions_width,alpha=regions_alpha)
    rivers_gdf.plot(ax=ax,color=river_color,linewidth=0.3,alpha=river_alpha,zorder=9)
    lakes_gdf.plot(ax=ax,color=lake_color,alpha=lake_alpha,zorder=8)
    gdf_eval_plot.plot(ax=ax, column=stats[i], cmap=cmap,markersize=markersize, norm=norm,vmin=vmin[i],vmax=vmax[i],marker='s',zorder=10)

    ax.set_xlim(lon_limits)
    ax.set_ylim(lat_limits)

    # Colorbar
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.02, pad=0.02,shrink=0.4)
    cbar.set_label(statName[i])
    cbar.set_ticks(bounds[i])
    cbar.set_ticklabels([f"{b:.1f}" for b in bounds[i]])

    ax.text(0.01, 0.08, f"{chr(97+i)})", transform=ax.transAxes, ha='left',fontweight="bold", va='top',fontsize=12)

    ax.yaxis.set_major_formatter(mticker.FuncFormatter(fmt_lat))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(fmt_lon))

plt.savefig(f'{stat_plots_dir}/figureb1.jpg',bbox_inches="tight")
plt.close()


#===============================================================================
# Performance relationships
#===============================================================================


# Just storage vs reservoir properties
gdf_eval['smax']      = pd.to_numeric(gdf_eval['smax'], errors='coerce')
gdf_eval['parea']     = pd.to_numeric(gdf_eval['parea'], errors='coerce')
gdf_eval['elev_masl'] = pd.to_numeric(gdf_eval['elev_masl'], errors='coerce')
gdf_eval['maxDepth']  = pd.to_numeric(gdf_eval['maxDepth'], errors='coerce')
gdf_eval['Depth_avg'] = pd.to_numeric(gdf_eval['Depth_avg'], errors='coerce')
gdf_eval['Res_time']  = pd.to_numeric(gdf_eval['Res_time'], errors='coerce')
gdf_eval['Max Storage (Km3)'] = gdf_eval['smax']/1e9
gdf_eval['Max Area (Km2)'] = gdf_eval['parea']/1e6
cols = ["Max Storage (Km3)", "Max Area (Km2)", "elev_masl", "maxDepth", "Res_time", "Depth_avg"]
gdf_eval_clean = gdf_eval[gdf_eval['kge']>-0.41]

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
for ax, col in zip(axes.flatten(), cols):
    ax.scatter(gdf_eval_clean[col], gdf_eval_clean["kge"], s=10,marker='^',color='k')
    ax.set_xlabel(col)
    ax.set_ylabel("KGE")
    
    x = gdf_eval_clean[col]
    y = gdf_eval_clean["kge"]
    r, p = pearsonr(x, y)
    
    if col == 'Res_time':
        ax.set_xlim(0,1500)
        
    if col == 'Max Storage (Km3)':
        ax.set_xlim(0,40)
        
    if col == 'Max Area (Km2)':
        ax.set_xlim(0,0.4)
        
    if col == 'Depth_avg':
        ax.set_xlim(0,140)

    ax.text(0.95, 0.05,f"r = {r:.2f}\np = {p:.3f}",transform=ax.transAxes,ha="right",va="bottom",fontsize=10)
    ax.set_ylim(-0.4,1)

plt.tight_layout()
plt.savefig(f'{stat_plots_dir}/scatter_1.jpg',bbox_inches="tight")
plt.close()

# Copies of existing thingies
storage_eval_gdf    = gdf_eval.copy()    
df_inflow_perf      = df_gral_outf.copy()
df_outflow_perf     = df_gral_inf.copy()

# Create unique dataframe to compare
storage_inflow_eval_gdf = storage_eval_gdf[storage_eval_gdf['grand_id'].isin(df_inflow_perf['grand_id'].to_list())]
storage_outflow_eval_gdf = storage_eval_gdf[storage_eval_gdf['grand_id'].isin(df_outflow_perf['grand_id'].to_list())]

df_inflow_perf = df_inflow_perf.rename(columns={'kge':'kge_in','pbias':'pbias_in','R':'R_in','alpha':'alpha_in','beta':'beta_in'})
df_outflow_perf = df_outflow_perf.rename(columns={'kge':'kge_out','pbias':'pbias_out','R':'R_out','alpha':'alpha_out','beta':'beta_out'})

storage_inflow_eval_gdf = pd.merge(storage_inflow_eval_gdf,df_inflow_perf,on='grand_id',how='left')
storage_outflow_eval_gdf = pd.merge(storage_outflow_eval_gdf,df_outflow_perf,on='grand_id',how='left')


# Scatter plots relating inflow-outflow-storage
storage_inflow_eval_gdf = storage_inflow_eval_gdf[storage_inflow_eval_gdf['kge']>-0.41]
storage_inflow_eval_gdf = storage_inflow_eval_gdf[storage_inflow_eval_gdf['kge_in']>-10]

storage_outflow_eval_gdf = storage_outflow_eval_gdf[storage_outflow_eval_gdf['kge']>-0.41]
storage_outflow_eval_gdf = storage_outflow_eval_gdf[storage_outflow_eval_gdf['kge_out']>-10]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].scatter(storage_inflow_eval_gdf["kge_in"], storage_inflow_eval_gdf["kge"], s=10,marker='s',color='k')
axes[0].set_xlabel("Inflow KGE")
axes[0].set_ylabel("Storage KGE")


x = storage_inflow_eval_gdf["kge_in"]
y = storage_inflow_eval_gdf["kge"]
r, p = pearsonr(x, y)
axes[0].text(0.95, 0.05,f"r = {r:.2f}\np = {p:.3f}",transform=axes[0].transAxes,ha="right",va="bottom",fontsize=10)
m, b = np.polyfit(x, y, 1)
axes[0].plot(x, m*x + b, color='gray', linewidth=0.8,linestyle='--')

axes[1].scatter(storage_outflow_eval_gdf["kge_out"], storage_outflow_eval_gdf["kge"], s=10,marker='p',color='k')
axes[1].set_xlabel("Outflow KGE")
axes[1].set_ylabel("Storage KGE")

x = storage_outflow_eval_gdf["kge_out"]
y = storage_outflow_eval_gdf["kge"]
r, p = pearsonr(x, y)
axes[1].text(0.95, 0.05,f"r = {r:.2f}\np = {p:.4f}",transform=axes[1].transAxes,ha="right",va="bottom",fontsize=10)
m, b = np.polyfit(x, y, 1)
axes[1].plot(x, m*x + b, color='gray', linewidth=0.8,linestyle='--')


fig, axes = plt.subplots(1, 3, figsize=(17, 5))

# 1) Storage KGE vs Inflow & Outflow KGE
ax = axes[0]

# Inflow
x_in = storage_inflow_eval_gdf["kge_in"]
y_in = storage_inflow_eval_gdf["kge"]
ax.scatter(x_in, y_in, s=10, marker='s', color='k', label='Inflow')
R, p = pearsonr(x_in, y_in)
m, b = np.polyfit(x_in, y_in, 1)
ax.plot(x_in, m*x_in + b, color='k', linestyle='--', linewidth=1)

# Outflow
x_out = storage_outflow_eval_gdf["kge_out"]
y_out = storage_outflow_eval_gdf["kge"]
ax.scatter(x_out, y_out, s=10, marker='p', color='tab:red', label='Outflow')
R2, p2 = pearsonr(x_out, y_out)
m, b = np.polyfit(x_out, y_out, 1)
ax.plot(x_out, m*x_out + b, color='tab:red', linestyle='--', linewidth=1)
ax.text(0.02, 0.95, f"Inflow\nR = {R:.2f}\np = {p:.3f}",transform=ax.transAxes, va="top", fontsize=9)
ax.text(0.25, 0.95, f"Outflow\nR = {R2:.2f}\np = {p2:.3f}",transform=ax.transAxes, va="top", fontsize=9)
ax.set_xlabel("Inflow / Outflow KGE")
ax.set_ylabel("Storage KGE")
ax.legend(frameon=False)

# 2) Storage KGE vs Inflow R (only)

ax = axes[1]
x_in = storage_inflow_eval_gdf["R_in"]
y_in = storage_inflow_eval_gdf["kge"]
ax.scatter(x_in, y_in, s=10, marker='s', color='k')
R, p = pearsonr(x_in, y_in)
ax.text(0.02, 0.95, f"Inflow\nR = {R:.2f}\np = {p:.3f}",transform=ax.transAxes, va="top", fontsize=9)
m, b = np.polyfit(x_in, y_in, 1)
ax.plot(x_in, m*x_in + b, color='gray', linestyle='--', linewidth=1)
ax.set_xlabel("Inflow R")
ax.set_ylabel("Storage KGE")

# 3) Elevation vs PBIAS (unchanged)

gdf_eval_clean_plot = gdf_eval_clean[gdf_eval_clean['elev_masl'] > 0]

ax = axes[2]
ax.scatter(gdf_eval_clean_plot["elev_masl"],gdf_eval_clean_plot["pbias"],s=10, marker='^', color='k')
ax.set_xlabel("Elevation (masl)")
ax.set_ylabel("Storage PBIAS")
ax.set_xlim(0, 3000)

x = gdf_eval_clean_plot["elev_masl"]
y = gdf_eval_clean_plot["pbias"]
R, p = pearsonr(x, y)

ax.text(0.95, 0.05, f"R = {R:.2f}\np = {p:.4f}",transform=ax.transAxes, ha="right", va="bottom", fontsize=10)

m, b = np.polyfit(x, y, 1)
ax.plot(x, m*x + b, color='gray', linestyle='--', linewidth=0.8)
labels = ['a)', 'b)', 'c)']
for ax, lab in zip(axes, labels):
    ax.text(-0.12, 1.05, lab,transform=ax.transAxes,fontsize=12, fontweight='bold',va='top', ha='left')

plt.tight_layout()
plt.savefig(f'{stat_plots_dir}/figureb2.jpg',bbox_inches="tight")
plt.close()