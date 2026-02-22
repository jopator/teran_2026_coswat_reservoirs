"""
Author: Jose P. Teran
Date: 2025-12-03

Description:
Analysis of reservoir storage performance per region and across regions (monthly and yearly)
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import hydroeval as he
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import matplotlib.lines as mlines
import matplotlib.cm as cm
from matplotlib.patches import Rectangle
import pickle as pkl
import sys
from pathlib import Path

# Function
def summary_row(series):
    s = series.dropna()
    q1 = s.quantile(0.25)
    q3 = s.quantile(0.75)
    return {
        "mean": s.mean(),
        "median": s.median(),
        "min": s.min(),
        "max": s.max(),
        "q1": q1,
        "q3": q3,
        "iqr": q3 - q1,
        "n": len(s)
    }

def lon_formatter(x, pos):
    return f"{abs(int(x))}°{'W' if x < 0 else 'E'}"

def lat_formatter(y, pos):
    return f"{abs(int(y))}°{'S' if y < 0 else 'N'}"


# Paths settings and file names
# ====================
BASE_DIR = Path(__file__).resolve().parents[4]  # root of repo
os.chdir(BASE_DIR)  # All relative paths will be based on this ! !

versionRes         = "1.5.0"
versionNores       = "1.1.0"

analysis_folder = f'Scripts/result_analysis'#f'vsc10883/PHD_main/Projects/CoSWAT/Scripts/result_analysis'
modelSetupPath  = 'CoSWAT-Framework/model-setup/CoSWATv{version}/{region}'
modelDataPath   = 'CoSWAT-Framework/model-data/{region}'
wshedGisPath    = f'{modelSetupPath}/Watershed'
reservoir_shp   = f'{wshedGisPath}/Shapes/lakes-grand-ESRI-54003.shp'
wshed_shp       = f'{wshedGisPath}/Shapes/subsNoLakes.shp'
stat_dir        = f'{analysis_folder}/reservoirs_coswat/validation/stats'
stat_plots_dir  = f'{analysis_folder}/reservoirs_coswat/validation/analysis/plots'
global_shapes   = f'Scripts/result_analysis/reservoirs_coswat/validation/shapefiles'

qperf_generalShp = f'{stat_dir}/streamflowqPerf_general.gpkg'
qperf_generalCsv = f'{stat_dir}/streamflowqPerf_general.csv'
qperf_ResInfCsv  = f'{stat_dir}/streamflowqPerf_ResInfluence.csv'

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
    "kge":  "#d28b84",
    "r":    "#5DA2CD",
    "pbias": "#485b48",
    "alpha": "#4daab7",
    "beta": "#becc27"
}

#======================================================================================================
# KGE Skill score and performance change categorization
# Separate:
# - If skill is reducing,   but the station was bad anyways (kgeNoRes<-0.41)   qSkill_category  = 4
# - If skill is reducing,   but the station was ok          (kgeNoRes>-0.41)   qSkill_category  = 3
# - If skill is increasing, and the station was bad         (kgeNoRes<-0.41)   qSkill_category  = 2
# - If skill is increasing, and the station was ok          (kgeNoRes>-0.41)   qSkill_category  = 1
#======================================================================================================
change_threshold = 0.075       # Change threshold for KGE skill score



# Calculate skills
qperf_df = pd.read_csv(qperf_generalCsv, index_col=0)
# filter out rows where kge (with reservoir) or kge (no reservoir) is NaN
qperf_df = qperf_df.dropna(subset=['kgeqRes', 'kgeqNoRes']).copy()
qperf_df['kgeqSkill']   = (qperf_df['kgeqRes'] - qperf_df['kgeqNoRes'])/(1 - qperf_df['kgeqNoRes'])
qperf_df['kge_change']   = (qperf_df['kgeqRes'] - qperf_df['kgeqNoRes'])

# qperf_df = qperf_df[qperf_df['kgeqSkill']!=0]
qperf_df['kgeqSkill_abs'] = qperf_df['kgeqSkill'].abs()
qperf_df['kge_change_abs'] = qperf_df['kge_change'].abs()
qperf_df_original = qperf_df.copy()

qperf_df = qperf_df[qperf_df['kgeqSkill_abs']>change_threshold]
qperf_df['kgeqSkill_category'] = 0
qperf_df.loc[(qperf_df['kgeqSkill'] < 0) & (qperf_df['kgeqNoRes'] < 0), 'kgeqSkill_category']   = 4
qperf_df.loc[(qperf_df['kgeqSkill'] < 0) & (qperf_df['kgeqNoRes'] >= 0), 'kgeqSkill_category']  = 3
qperf_df.loc[(qperf_df['kgeqSkill'] >= 0) & (qperf_df['kgeqNoRes'] < 0), 'kgeqSkill_category']  = 2
qperf_df.loc[(qperf_df['kgeqSkill'] >= 0) & (qperf_df['kgeqNoRes'] >= 0), 'kgeqSkill_category'] = 1
qperf_df_original.loc[(qperf_df_original['kge_change_abs'] < change_threshold), 'kgeqSkill_category'] = 0
qperf_df_original.loc[(qperf_df_original['kge_change_abs'] == 0), 'kgeqSkill_category'] = -1


log = open(f"{analysis_folder}/reservoirs_coswat/validation/analysis/OutLogStreamflow.txt", "w")
sys.stdout = log


print("\nSummary of stations by KGE skill category:")
print(qperf_df['kgeqSkill_category'].value_counts().sort_index())
print("\nCategory breakdown:")
print(f"Category 4 (skill reducing, station was bad): {(qperf_df['kgeqSkill_category'] == 4).sum()} ({(qperf_df['kgeqSkill_category'] == 4).sum()/len(qperf_df)*100:.1f} %)")
print(f"Category 3 (skill reducing, station was ok): {(qperf_df['kgeqSkill_category'] == 3).sum()} ({(qperf_df['kgeqSkill_category'] == 3).sum()/len(qperf_df)*100:.1f} %)")
print(f"Category 2 (skill increasing, station was bad): {(qperf_df['kgeqSkill_category'] == 2).sum()} ({(qperf_df['kgeqSkill_category'] == 2).sum()/len(qperf_df)*100:.1f} %)")
print(f"Category 1 (skill increasing, station was ok): {(qperf_df['kgeqSkill_category'] == 1).sum()} ({(qperf_df['kgeqSkill_category'] == 1).sum()/len(qperf_df)*100:.1f} %)")
print(f"Category 0 no change - Based on total stations: {(qperf_df_original['kgeqSkill_category'] == 0).sum()} ({(qperf_df_original['kgeqSkill_category'] == 0).sum()/len(qperf_df_original)*100:.1f} %)")
print(f"Category -1 no change - Based on total stations: {(qperf_df_original['kgeqSkill_category'] == -1).sum()} ({(qperf_df_original['kgeqSkill_category'] == -1).sum()/len(qperf_df_original)*100:.1f} %)")

print(f"\nTotal improving: {(qperf_df['kgeqSkill'] > 0).sum()} ({(qperf_df['kgeqSkill'] > 0).sum()/len(qperf_df)*100:.1f} %)")
print(f"Total becoming worse: {(qperf_df['kgeqSkill'] < 0).sum()} ({(qperf_df['kgeqSkill'] < 0).sum()/len(qperf_df)*100:.1f} %)")


print(f"\nStations in category 3 flipping to negative: {qperf_df[(qperf_df['kgeqSkill_category']==3) & (qperf_df['kgeqRes']<0) & (qperf_df['kgeqNoRes']>0)].shape[0]}")
print(f"Stations in category 2 flipping to positive: {qperf_df[(qperf_df['kgeqSkill_category']==2) & (qperf_df['kgeqRes']>0) & (qperf_df['kgeqNoRes']<0)].shape[0]}")
print(f"Category 3 remaining above 0: {qperf_df[(qperf_df['kgeqSkill_category']==3) & (qperf_df['kgeqRes']>0) & (qperf_df['kgeqNoRes']>0)].shape[0]}")
print(f"Category 2 remaining below 0: {qperf_df[(qperf_df['kgeqSkill_category']==2) & (qperf_df['kgeqRes']<0) & (qperf_df['kgeqNoRes']<0)].shape[0]}")
print(f"\nStations in category 3 crossing below -0.41: {qperf_df[(qperf_df['kgeqSkill_category']==3) & (qperf_df['kgeqRes']<-0.41) & (qperf_df['kgeqNoRes']>-0.41)].shape[0]}")
print(f"Stations in category 2 crossing above -0.41: {qperf_df[(qperf_df['kgeqSkill_category']==2) & (qperf_df['kgeqRes']>-0.41) & (qperf_df['kgeqNoRes']<-0.41)].shape[0]}")

print(f"Category 3 remaining above -0.41: {qperf_df[(qperf_df['kgeqSkill_category']==3) & (qperf_df['kgeqRes']>-0.41) & (qperf_df['kgeqNoRes']>-0.41)].shape[0]}")

print(f"Category 2 remaining above -0.41: {qperf_df[(qperf_df['kgeqSkill_category']==2) & (qperf_df['kgeqRes']>-0.41) & (qperf_df['kgeqNoRes']>-0.41)].shape[0]}")
print(f"Category 2 remaining below -0.41: {qperf_df[(qperf_df['kgeqSkill_category']==2) & (qperf_df['kgeqRes']<-0.41) & (qperf_df['kgeqNoRes']<-0.41)].shape[0]}")

#======================================================================================================
# Table report of general performance
#======================================================================================================

print("\n===== PERFORMANCE SUMMARY (HORIZONTAL) =====\n")

header = ("STAT   | CASE     MEAN     MEDIAN     MIN      MAX      Q1       Q3      IQR      N")
print(header)
print("-" * len(header))

for stat in ['kge','pbias','r']:

    r_no = summary_row(qperf_df[f'{stat}qNoRes'])
    r_re = summary_row(qperf_df[f'{stat}qRes'])

    print(
        f"{stat.upper():<6}| Without "
        f"{r_no['mean']:>8.2f}  "
        f"{r_no['median']:>8.2f}  "
        f"{r_no['min']:>8.2f}  "
        f"{r_no['max']:>8.2f}  "
        f"{r_no['q1']:>8.2f}  "
        f"{r_no['q3']:>8.2f}  "
        f"{r_no['iqr']:>8.2f}  "
        f"{r_no['n']:>6d}"
    )

    print(
        f"{stat.upper():<6}| With    "
        f"{r_re['mean']:>8.2f}  "
        f"{r_re['median']:>8.2f}  "
        f"{r_re['min']:>8.2f}  "
        f"{r_re['max']:>8.2f}  "
        f"{r_re['q1']:>8.2f}  "
        f"{r_re['q3']:>8.2f}  "
        f"{r_re['iqr']:>8.2f}  "
        f"{r_re['n']:>6d}"
    )

    print("-" * len(header))

print("\n===== STATIONS PER THRESHOLD CATEGORY =====\n")

kge_cat = (
    pd.DataFrame({
        "without": qperf_df["kgeqNoRes"],
        "with": qperf_df["kgeqRes"]
    })
    .apply(lambda s: pd.Series({
        "min": ((s >= -0.42) & (s <= 0.00)).sum(),
        "satisfactory": ((s >= 0.0) & (s <= 0.40)).sum(),
        "good": (s > 0.40).sum()
    }))
)

pb_cat = (
    pd.DataFrame({
        "without": qperf_df["pbiasqNoRes"],
        "with": qperf_df["pbiasqRes"]
    })
    .apply(lambda s: pd.Series({
        "within_±25": (s.abs() <= 25).sum(),
        "outside_±25": (s.abs() > 25).sum(),
        "within_±50": ((s.abs() <= 50) & (s.abs() >= 25)).sum(),
        "outside_±50": (s.abs() > 50).sum(),
    }))
)

r_cat = (
    pd.DataFrame({
        "without": qperf_df["rqNoRes"],
        "with": qperf_df["rqRes"]
    })
    .apply(lambda s: pd.Series({
        "good_R": (s > 0.40).sum(),
        "poor_R": (s <= 0.40).sum()
    }))
)

print("KGE (stations):")
print(kge_cat, "\n")

print("PBIAS (stations):")
print(pb_cat, "\n")

print("R (stations):")
print(r_cat)


log.close()

#======================================================================================================
# General performance boxplot
#======================================================================================================
fig, axes = plt.subplots(1, 3, figsize=(10, 4))

# dummy lines only for legend (drawn once)
h_sat, = plt.plot([], [], linestyle="--",  color="grey", linewidth=0.8, label="Satisfactory")
h_good, = plt.plot([], [], linestyle="-.", color="grey", linewidth=0.8, label="Good")
h_min,  = plt.plot([], [], linestyle=":",  color="grey", linewidth=0.8, label="Minimum Expected")

for idx, stat in enumerate(['kge','pbias','r']):
    ax = axes[idx]

    data = [qperf_df[f'{stat}qNoRes'].dropna(),qperf_df[f'{stat}qRes'].dropna()]
    ax.boxplot(data,tick_labels=['Without', 'With'],showfliers=False,patch_artist=True,widths=0.40,boxprops=dict(facecolor=colors[stat], alpha=0.4),
        medianprops=dict(color='black', linestyle='-', linewidth=2),whiskerprops=dict(color='black'),capprops=dict(color='black'))

    if stat == 'kge':
        ax.axhline(0.0,   linestyle='--', linewidth=0.8, color='grey')
        ax.axhline(0.4,   linestyle='-.', linewidth=0.8, color='grey')
        ax.axhline(-0.42, linestyle=':',  linewidth=0.8, color='grey')
        ax.set_ylim(-5, 1)

    if stat == 'pbias':
        ax.axhline(25.0,  linestyle='--', linewidth=0.8, color='grey')
        ax.axhline(-25.0, linestyle='--', linewidth=0.8, color='grey')
        ax.set_ylim(-500, 150)

    if stat == 'r':
        ax.axhline(0.4,  linestyle='-.', linewidth=0.8, color='grey')
        ax.set_ylim(0, 1)

    ax.set_xlabel(stat.upper())
    ax.set_ylabel(stat.upper())
    ax.grid(alpha=0.3,axis = 'x')
    label = chr(97 + idx)  # a, b, c
    ax.text(-0.15, 0.95, f"{label})",transform=ax.transAxes,ha="right", va="bottom",fontsize=11, fontweight="bold")

# single common legend at bottom (by line style)
fig.legend(handles=[h_sat, h_good, h_min],loc="lower center",ncol=3,frameon=False)

plt.tight_layout(rect=[0, 0.08, 1, 1])
plt.savefig(f'{stat_plots_dir}/streamflow_per_boxplot.jpg',bbox_inches="tight")
plt.close()



#======================================================================================================
# Spatial distribuion of categories, and box plots per category
#======================================================================================================
qperf_gdf = gpd.read_file(qperf_generalShp)
qperf_gdf = qperf_gdf.dropna(subset=['kgeqRes', 'kgeqNoRes','pbiasqRes','pbiasqNoRes']).copy()
qperf_gdf = pd.merge(qperf_gdf,qperf_df[['grdc_no','kgeqSkill_category','kgeqSkill']],on='grdc_no',how='right')

url = "https://github.com/nvkelso/natural-earth-vector/raw/master/110m_physical/ne_110m_land.shp"
world = gpd.read_file(url)

basins_file = f"{global_shapes}/basins_analysis.gpkg"
basins_gdf  = gpd.read_file(basins_file).to_crs("EPSG:4326")

rivers_file = f"{global_shapes}/major_rivers_analysis.gpkg"
rivers_gdf  = gpd.read_file(rivers_file,layer='output').to_crs("EPSG:4326")

lakes_file  = f"{global_shapes}/all_lakes.gpkg"
lakes_gdf   = gpd.read_file(lakes_file).to_crs("EPSG:4326")


# Plot configuration
# Maps
color_map = {1: "green", 2: "royalblue", 3: "orange", 4: "red"}
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

# Box plots
box_width = 0.3
box_sep   = 0.175
box_alhpa = 0.80
box_whis  = (10,100)


gdf_eval_plot = qperf_gdf.to_crs(world.crs)

fig = plt.figure(figsize=(24, 12))
gs = gridspec.GridSpec(2, 2, width_ratios=[1, 0.42], height_ratios=[1, 1],wspace=0.02,hspace=0.09)

ax1 = fig.add_subplot(gs[0, 0])
ax3 = fig.add_subplot(gs[1, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax4 = fig.add_subplot(gs[1, 1])


# Map 1
gdf_eval_plot["cat_color"] = gdf_eval_plot["kgeqSkill_category"].map(color_map)
world.plot(ax=ax1,color=continent_color, alpha=continent_alpha, aspect='auto')
world.boundary.plot(ax=ax1, linewidth=continent_width, color="black", aspect='auto')
basins_gdf.plot(ax=ax1,color='gray',alpha=regions_alpha)
basins_gdf.boundary.plot(ax=ax1,color='k',linewidth=regions_width,alpha=regions_alpha)
rivers_gdf.plot(ax=ax1,color=river_color,linewidth=0.3,alpha=river_alpha,zorder=9)
lakes_gdf.plot(ax=ax1,color=lake_color,alpha=lake_alpha,zorder=8)
gdf_eval_plot.plot(ax=ax1, color=gdf_eval_plot["cat_color"], markersize=markersize, marker=marker_cat,aspect='auto',zorder=10)


ax1.set_xlim(lon_limits)
ax1.set_ylim(lat_limits)
ax1.xaxis.set_major_formatter(mticker.FuncFormatter(lon_formatter))
ax1.yaxis.set_major_formatter(mticker.FuncFormatter(lat_formatter))

handles = [ mpatches.Patch(color=color_map[1], label="Category 1"),
            mpatches.Patch(color=color_map[2], label="Category 2"),
            mpatches.Patch(color=color_map[3], label="Category 3"),
            mpatches.Patch(color=color_map[4], label="Category 4")]

ax1.legend(handles=handles, loc="lower left", frameon = False,title="Skill Category",fontsize=11.5,title_fontsize=12)
ax1.tick_params(axis='x', labelbottom=False)

# Map 2
bounds  = kge_bounds
cmap    = mcolors.ListedColormap(kge_colors)
cmap.set_under("red")
cmap.set_over("darkgreen")
norm    = mcolors.BoundaryNorm(bounds, cmap.N)

world.plot(ax=ax3, color=continent_color, alpha=continent_alpha, aspect='auto')
world.boundary.plot(ax=ax3, linewidth=continent_width, color="black", aspect='auto')
gdf_eval_plot.plot(ax=ax3,column="kgeqRes",markersize=markersize,marker=marker_kge,cmap=cmap,norm=norm,legend=False,zorder=10)
basins_gdf.plot(ax=ax3,color='gray',alpha=regions_alpha)
basins_gdf.boundary.plot(ax=ax3,color='k',linewidth=regions_width,alpha=regions_alpha)
rivers_gdf.plot(ax=ax3,color=river_color,linewidth=0.3,alpha=river_alpha,zorder=9)
lakes_gdf.plot(ax=ax3,color=lake_color,alpha=lake_alpha,zorder=8)

ax3.set_xlim(lon_limits)
ax3.set_ylim(lat_limits)
ax3.xaxis.set_major_formatter(mticker.FuncFormatter(lon_formatter))
ax3.yaxis.set_major_formatter(mticker.FuncFormatter(lat_formatter))

sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

cax = fig.add_axes([0.18, 0.155, 0.08, 0.0125])

cbar = fig.colorbar(sm,cax=cax,orientation="horizontal",extend="both",shrink=0.15)
cbar.set_label("KGE (With)",fontsize=12)
cbar.set_ticks(kge_bounds)
cbar.ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))
cbar.ax.xaxis.set_label_position('top')
cbar.ax.tick_params(labelsize=11.5)

# Box plot
categories = [4,3,2,1]
positions = []
box_data = []
for i, cat in enumerate(categories, start=1):
    df_cat = qperf_df[(qperf_df['kgeqSkill_category'] == cat) & (qperf_df['kgeqNoRes'] > -4.5)]
    box_data.append(df_cat['kgeqRes'].dropna())
    positions.append(i - box_sep)
    box_data.append(df_cat['kgeqNoRes'].dropna())
    positions.append(i + box_sep)

bp = ax2.boxplot(box_data,positions=positions,widths=box_width,patch_artist=True,manage_ticks=False,vert=False,whis = box_whis,showfliers=False)

for j, box in enumerate(bp['boxes']):
    face = "cadetblue" if j % 2 == 0 else "rosybrown"
    box.set_facecolor(face)
    box.set_edgecolor('black')
    box.set_alpha(box_alhpa)

for median in bp['medians']:
    median.set_color('black')
    median.set_linewidth(1.5)

for whisker in bp['whiskers']:
    whisker.set_color('black')

for cap in bp['caps']:
    cap.set_color('black')

ax2.set_yticks([1, 2, 3, 4])
ax2.set_yticklabels([str(c) for c in categories])
ax2.set_ylabel("Skill Category",fontsize=12)
ax2.set_xlabel("KGE",fontsize=12)
ax2.set_xlim(-2.5, 1)

ax2.axvline(0.0,   linestyle='--', linewidth=1.2, color='grey')
ax2.axvline(0.4,   linestyle='dashdot', linewidth=0.8, color='grey')
ax2.axvline(-0.42, linestyle=':',  linewidth=0.8, color='grey')


h_sat  = mlines.Line2D([], [], color='grey', linestyle='--',      linewidth=1.2, label='Satisfactory')
h_good = mlines.Line2D([], [], color='grey', linestyle='dashdot',  linewidth=0.8, label='Good')
h_min  = mlines.Line2D([], [], color='grey', linestyle=':',        linewidth=0.8, label='Minimum')

handles = [mpatches.Patch(color="cadetblue", label='With'),mpatches.Patch(color="rosybrown", label='Without'),h_sat,h_good,h_min]
ax2.legend(handles=handles, loc='upper left',frameon=False,fontsize=11.5)
ax2.set_box_aspect(0.85)


# Box plot 2
categories = [4,3,2,1]
positions = []
box_data = []
for i, cat in enumerate(categories, start=1):
    df_cat = qperf_df[(qperf_df['kgeqSkill_category'] == cat) & (qperf_df['kgeqNoRes'] > -4.5)]
    box_data.append(df_cat['pbiasqRes'].dropna())
    positions.append(i - box_sep)
    box_data.append(df_cat['pbiasqNoRes'].dropna())
    positions.append(i + box_sep)

bp = ax4.boxplot(box_data,positions=positions,widths=box_width,patch_artist=True,manage_ticks=False,vert=False,whis = box_whis,showfliers=False)

for j, box in enumerate(bp['boxes']):
    face = "cadetblue" if j % 2 == 0 else "rosybrown"
    box.set_facecolor(face)
    box.set_edgecolor('black')
    box.set_alpha(box_alhpa)

for median in bp['medians']:
    median.set_color('black')
    median.set_linewidth(1.5)

for whisker in bp['whiskers']:
    whisker.set_color('black')

for cap in bp['caps']:
    cap.set_color('black')

ax4.set_yticks([1, 2, 3, 4])
ax4.set_yticklabels([str(c) for c in categories])
ax4.set_ylabel("Skill Category",fontsize=12)
ax4.set_xlabel("PBIAS (%)",fontsize=12)
ax4.set_xlim(-200, 100)

ax4.axvline(-25,   linestyle='--', linewidth=1.2, color='grey')
ax4.axvline(25,   linestyle='--', linewidth=1.2, color='grey')
ax4.axvline(50,   linestyle='dashdot', linewidth=0.8, color='grey')
ax4.axvline(-50, linestyle='dashdot',  linewidth=0.8, color='grey')

handles = [mpatches.Patch(color="cadetblue", label='With'),mpatches.Patch(color="rosybrown", label='Without'),h_sat,h_good]
ax4.legend(handles=handles, loc='upper left',frameon=False,fontsize=11.5)
ax4.set_box_aspect(0.85)


fig.canvas.draw()

pad_x = 0.022
pad_x_right = 0.01
pad_y_bottom = 0.05
pad_y_top = 0.02

# Left column
left_axes = [ax1, ax3]
bxs = [ax.get_position() for ax in left_axes]

x0 = min(b.x0 for b in bxs) - pad_x
y0 = min(b.y0 for b in bxs) - pad_y_bottom
x1 = max(b.x1 for b in bxs) + pad_x_right
y1 = max(b.y1 for b in bxs) + pad_y_top

box_left = Rectangle((x0, y0), x1 - x0, y1 - y0,transform=fig.transFigure,facecolor="none",edgecolor="gray",linewidth=0.8)
fig.add_artist(box_left)

# Right column
right_axes = [ax2, ax4]
bxs = [ax.get_position() for ax in right_axes]

x0 = min(b.x0 for b in bxs) - pad_x - 0.004
y0 = min(b.y0 for b in bxs) - pad_y_bottom-0.001
x1 = max(b.x1 for b in bxs) + pad_x_right
y1 = max(b.y1 for b in bxs) + pad_y_top+0.001

box_right = Rectangle((x0, y0), x1 - x0, y1 - y0,transform=fig.transFigure,facecolor="none",edgecolor="gray",linewidth=0.8)
fig.add_artist(box_right)


fig.text(box_left.get_x() + 0.002,  box_left.get_y() + box_left.get_height() - 0.02,  "a)", fontsize=14, fontweight="bold")
fig.text(box_right.get_x() + 0.002, box_right.get_y() + box_right.get_height() - 0.02, "b)", fontsize=14, fontweight="bold")
fig.savefig(f'{stat_plots_dir}/figure8.jpg', dpi=300, bbox_inches="tight", pad_inches=0.02)

plt.close()


#======================================================================================================
# Spatial distribuion of categories, and box plots per category
#======================================================================================================

def metrics(sim, obs):
    sim = np.asarray(sim, dtype=float)
    obs = np.asarray(obs, dtype=float)

    m = np.isfinite(sim) & np.isfinite(obs)
    sim = sim[m]
    obs = obs[m]

    kge   = float(np.asarray(he.evaluator(he.kge,   sim, obs)).ravel()[0])
    pbias = float(np.asarray(he.evaluator(he.pbias, sim, obs)).ravel()[0])
    r     = float(np.corrcoef(sim, obs)[0, 1]) if sim.size > 1 else np.nan

    return kge, pbias, r, sim.size

output_data_folder = f'{analysis_folder}/reservoirs_coswat/coswat_outputs/streamflow'
cat_1_grdc_df = gdf_eval_plot[(gdf_eval_plot['kgeqSkill_category']==1) & (gdf_eval_plot['kgeqNoRes']<0.5)].sort_values(by='kgeqRes')[['grdc_no','region','kgeqRes','kgeqNoRes']].sort_values(by='region').reset_index(drop=True)
cat_2_grdc_df = gdf_eval_plot[(gdf_eval_plot['kgeqSkill_category']==2) & (gdf_eval_plot['kgeqRes']>0) & (gdf_eval_plot['kgeqNoRes']<0)].sort_values(by='kgeqNoRes')[['grdc_no','region','kgeqRes','kgeqNoRes']].sort_values(by='region').reset_index(drop=True)
cat_3_grdc_df = gdf_eval_plot[(gdf_eval_plot['kgeqSkill_category']==3) & (gdf_eval_plot['kgeqRes']<0) & (gdf_eval_plot['kgeqNoRes']>0)].sort_values(by='kgeqNoRes')[['grdc_no','region','kgeqRes','kgeqNoRes']].sort_values(by='region').reset_index(drop=True)


pairs = [
    ("europe-west",         6113050),
    ("america-parana",      3754100),
    ("america-colorado",    4146112),
    ("america-mississippi", 4126350),
    ("europe-west",         6212505),
    ("asia-mekong",         2969096)
    ]

meta = (gdf_eval_plot.set_index("grdc_no")[["kgeqRes", "kgeqNoRes", "kgeqSkill_category","station"]].to_dict("index"))

fig, axs = plt.subplots(2, 3, figsize=(16, 6), sharex=True, sharey=False)
axs = axs.flatten()
plt.subplots_adjust(hspace=0.3,wspace=0.35)
handles = None
labels = None

for i, ax in enumerate(axs):
    region, station = pairs[i]

    data_Res_fn   = f'{output_data_folder}/CoSWATv1.5.0/{region}/monthly/channelsd_mon.pkl'
    data_NoRes_fn = f'{output_data_folder}/CoSWATv1.1.0/{region}/monthly/channelsd_mon.pkl'

    with open(data_Res_fn, "rb") as f:
        data_Res = pkl.load(f)
    with open(data_NoRes_fn, "rb") as f:
        data_NoRes = pkl.load(f)

    station_dataRes   = data_Res[station]
    station_dataNoRes = data_NoRes[station]

    simRes_df   = station_dataRes["sim_data"]
    simNoRes_df = station_dataNoRes["sim_data"]
    obs_df      = station_dataRes["obs_data"]

    simRes_df["date"]   = pd.to_datetime(simRes_df["date"])
    simNoRes_df["date"] = pd.to_datetime(simNoRes_df["date"])
    obs_df["date"]      = pd.to_datetime(obs_df["date"])

    t0 = max(obs_df["date"].min(), pd.Timestamp("1970-04-01"))
    t1 = min(obs_df["date"].max(), pd.Timestamp("2015-12-31"))

    simRes_df   = simRes_df.loc[simRes_df["date"].between(t0, t1)]
    simNoRes_df = simNoRes_df.loc[simNoRes_df["date"].between(t0, t1)]
    obs_df      = obs_df.loc[obs_df["date"].between(t0, t1)]

    clim_res_mean   = simRes_df.groupby(simRes_df["date"].dt.month)["sim"].mean()
    clim_res_std    = simRes_df.groupby(simRes_df["date"].dt.month)["sim"].std()

    clim_nores_mean = simNoRes_df.groupby(simNoRes_df["date"].dt.month)["sim"].mean()
    clim_nores_std  = simNoRes_df.groupby(simNoRes_df["date"].dt.month)["sim"].std()

    clim_obs_mean   = obs_df.groupby(obs_df["date"].dt.month)["obs"].mean()
    clim_obs_std    = obs_df.groupby(obs_df["date"].dt.month)["obs"].std()

    months = clim_obs_mean.index

    clim_res_mean.plot(ax=ax, color="teal", linewidth=1.5, label="With")
    ax.fill_between(
        months,
        clim_res_mean.reindex(months) - 1.5 * clim_res_std.reindex(months),
        clim_res_mean.reindex(months) + 1.5 * clim_res_std.reindex(months),
        color="teal",
        alpha=0.2,
        linewidth=0,
    )

    clim_nores_mean.plot(ax=ax, color="salmon", linewidth=1.5, label="Without")
    ax.fill_between(
        months,
        clim_nores_mean.reindex(months) - 1.5 * clim_nores_std.reindex(months),
        clim_nores_mean.reindex(months) + 1.5 * clim_nores_std.reindex(months),
        color="salmon",
        alpha=0.2,
        linewidth=0,
    )

    clim_obs_mean.plot(ax=ax, color="k", linewidth=2, linestyle="--", label="GRDC data")
    ax.fill_between(
        months,
        clim_obs_mean.reindex(months) - 1.5 * clim_obs_std.reindex(months),
        clim_obs_mean.reindex(months) + 1.5 * clim_obs_std.reindex(months),
        color="k",
        alpha=0.1,
        linewidth=0,
    )

    ymax = max(
        (clim_res_mean.reindex(months) + 1.5 * clim_res_std.reindex(months)).max(),
        (clim_nores_mean.reindex(months) + 1.5 * clim_nores_std.reindex(months)).max(),
        (clim_obs_mean.reindex(months) + 1.5 * clim_obs_std.reindex(months)).max(),
    )
    
    ax.set_ylim(0, 1.15 * ymax)
    ax.set_xlim(1, 12)

    m = meta.get(station, {})
    stat   = m.get("station", "")
    kge_r  = m.get("kgeqRes", np.nan)
    kge_nr = m.get("kgeqNoRes", np.nan)
    cat    = m.get("kgeqSkill_category", "")
    
    try:
        stat_name = stat.split(",")[0].title() + "," + stat.split(",")[1]
        
    except:
        stat_name = stat.title()
    ax.text(
        0.02, 0.98,
        f"KGE With: {kge_r:.2f}\nKGE Without: {kge_nr:.2f}\nSkill Category: {cat}",
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=10,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.4, pad=2.0),
    )
    ax.set_title(
        f"{stat_name} - GRDC: {station} ",
        fontsize=12,
    )

    if handles is None:
        handles, labels = ax.get_legend_handles_labels()

for ax in axs:
    ax.set_ylabel("Streamflow $(m^{3} s^{-1})$",fontsize=12)

for ax in axs[-3:]:
    ax.set_xlabel("Month of the year")

fig.legend(handles,labels,loc="lower center",bbox_to_anchor=(0.5,-0.03),ncol=3,frameon=False)


month_labels = ["Jan","Mar","May",
                "Jul","Sep","Nov"]

for ax in axs:
    ax.set_xticks(range(1,13,2))
    ax.set_xticklabels(month_labels)
    ax.grid(axis='x',linewidth = 0.4)

stat_plots_dir  = f'{analysis_folder}/reservoirs_coswat/validation/analysis/plots'
plt.savefig(f'{stat_plots_dir}/figure9.jpeg',dpi=300,bbox_inches='tight')
