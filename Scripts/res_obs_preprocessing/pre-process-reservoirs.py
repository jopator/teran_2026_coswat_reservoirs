# Imports
import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import pickle as pkl
import os
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parents[2]  # root of repo
os.chdir(BASE_DIR)

# Directories, names, variables
raw_data_dir    = 'Scripts/res_obs_preprocessing/Reservoir_storage_global_data/raw'
grs_data_dir    = 'GRS/data'
resops_data_dir = 'ResOpsUs/data'
yasin_data_dir  = 'yasin/data'
export_folder   = 'Scripts/res_obs_preprocessing/Reservoir_storage_global_data/processed/data'
grand_fn        = 'CoSWAT-Framework/data-preparation/resources/hydro-lakes/GRanD_Version_1_3/GRanD_reservoirs_v1_3.shp'
global_data_dict = {}

os.makedirs(f'{export_folder}',exist_ok=True)

print("> Starting preprocessing of reservoir storage - inflow - outflow data from GRS/ResOpsUS/Yassin2018")
#=================================
#GRS
#=================================
file_list = os.listdir(f'{raw_data_dir}/{grs_data_dir}')
grs_grand_ids = []
for file in file_list:
    fn          = os.path.abspath(f'{raw_data_dir}/{grs_data_dir}/{file}')
    grand_id    = int(file.split('.')[0])
    df          = pd.read_csv(fn)

    df_proc     = df.rename(columns={'Month':'date','Storage(km3)':'grs_storage'})
    df_proc['date'] = pd.to_datetime(df_proc['date'])
    df_proc     = df_proc[['date','grs_storage']]

    grs_grand_ids.append(grand_id)
    global_data_dict[grand_id] = df_proc


#=================================
#ResOpsUS
#=================================

# first, we need to aggregate across sources
os.makedirs(f'{raw_data_dir}/ResOpsUs/data_allSources',exist_ok=True)
file_list = os.listdir(f'{raw_data_dir}/{resops_data_dir}')
file_list.sort()
resops_data_dict = {}
old_grand_id = -9999

for file in file_list:
    fn          = os.path.abspath(f'{raw_data_dir}/{resops_data_dir}/{file}')
    grand_id    = int(file.split('_')[0])
    source      = file.split('_')[-1].split('.')[0]
    df          = pd.read_csv(fn)
    df_proc     = df.rename(columns={'storage':f'storage_{source}','inflow':f'inflow_{source}','outflow':f'outflow_{source}'})
    df_proc['date'] = pd.to_datetime(df_proc['date'])

    df_proc = df_proc[['date',f'storage_{source}',f'inflow_{source}',f'outflow_{source}']]
    

    if grand_id == old_grand_id:
        df_old    = resops_data_dict[grand_id]
        df_merged = pd.merge(df_proc,df_old,on='date',how='outer')
        resops_data_dict[grand_id] = df_merged
    else:
        resops_data_dict[grand_id] = df_proc
    
    old_grand_id = grand_id

# Agreggate across sources
for id,df in resops_data_dict.items():
    for var in ["storage", "inflow", "outflow"]:
        cols = [c for c in df.columns if c.startswith(f"{var}_")]
        if cols:
            df[f"{var}"] = df[cols].mean(axis=1, skipna=True)

    resops_data_dict[id] = df
    df.to_csv(f'{raw_data_dir}/ResOpsUs/data_allSources/{id}.csv',index=None)

file_list = os.listdir(f'{raw_data_dir}/ResOpsUs/data_allSources')
resops_grand_ids = []
for file in file_list:
    fn          = os.path.abspath(f'{raw_data_dir}/ResOpsUs/data_allSources/{file}')
    grand_id    = int(file.split('.')[0])
    df          = pd.read_csv(fn)
    df_proc     = df.rename(columns={'storage':'resops_storage','inflow':'resops_inflow','outflow':'resops_outflow'})
    df_proc['date'] = pd.to_datetime(df_proc['date'])
    df_proc['resops_storage'] = df_proc['resops_storage']/1000

    df_proc     = df_proc[['date','resops_storage','resops_inflow','resops_outflow']]

    df_proc.set_index('date',inplace=True)
    df_proc = df_proc.resample("MS").mean()
    df_proc.reset_index(inplace=True)

    # Read existing data (already has GRS)
    try:
        prev_df   = global_data_dict[grand_id]
        df_merged = pd.merge(df_proc,prev_df,on='date',how='outer')

    except:
        df_merged = df_proc.copy()
        df_merged['grs_storage'] = None


    resops_grand_ids.append(grand_id)
    global_data_dict[grand_id] = df_merged

# Fill reservoirs from GRS that are not in resopsus
for reservoir in grs_grand_ids:
    if reservoir not in resops_grand_ids:
        df = global_data_dict[reservoir]
        df['resops_storage'] = None
        df['resops_inflow']  = None
        df['resops_outflow'] = None
        global_data_dict[reservoir] = df

#=================================
#Yassin (2018)
#=================================
file_list = os.listdir(f'{raw_data_dir}/{yasin_data_dir}')
yasin_grand_ids = []
for file in file_list:
    fn          = os.path.abspath(f'{raw_data_dir}/{yasin_data_dir}/{file}')
    grand_id    = file.split('.')[0].split('_')[-1]
    grand_id    = int(grand_id)
    df          = pd.read_csv(fn)
    df_proc     = df.rename(columns={'stoobs':'yasin_storage','inflow':'yasin_inflow','outflow':'yasin_outflow'})
    df_proc['date'] = pd.to_datetime(df_proc['date'])
    df_proc['yasin_storage'] = df_proc['yasin_storage']/1000000000

    df_proc     = df_proc[['date','yasin_storage','yasin_inflow','yasin_outflow']]

    df_proc.set_index('date',inplace=True)
    df_proc = df_proc.resample("MS").mean()
    df_proc.reset_index(inplace=True)

    # Read existing data (already has GRS)
    try:
        prev_df   = global_data_dict[grand_id]
        df_merged = pd.merge(df_proc,prev_df,on='date',how='outer')

    except:
        df_merged = df_proc.copy()
        df_merged['grs_storage'] = None
        df_merged['resops_storage'] = None
        df_merged['resops_inflow']  = None
        df_merged['resops_outflow'] = None


    yasin_grand_ids.append(grand_id)
    global_data_dict[grand_id] = df_merged


print("> Aggregating into single dictionary...")

#=================================
#Aggregate
#=================================
global_data_dict_agg = {}

for id,df in global_data_dict.items():
    for var in ["storage", "inflow", "outflow"]:
        cols = [c for c in df.columns if c.endswith(f"{var}")]
        if cols:
            df[f"{var}_mean"]   = df[cols].mean(axis=1, skipna=True)
            df[f"{var}_min"]    = df[cols].min(axis=1, skipna=True)
            df[f"{var}_max"]    = df[cols].max(axis=1, skipna=True)

    df = df[['date','storage_mean','inflow_mean','outflow_mean','storage_min','inflow_min','outflow_min','storage_max','inflow_max','outflow_max']]

    global_data_dict_agg[id] = df


# Filter out
# At least 5 years of data and 1 year of continuous data!
def keep_df(df):
    d = df.copy()
    d['date'] = pd.to_datetime(d['date'])
    s = d.loc[d['storage_mean'].notna(), 'date'].sort_values()

    # ≥60 monthly points
    if s.size < 36:
        return False

    # longest run of consecutive months (≥12)
    p = s.dt.to_period('M').drop_duplicates()
    runs = (p != p.shift() + 1).cumsum()
    longest = runs.groupby(runs).size().max() if len(p) else 0

    return longest >= 12

# filter dict
global_data_dict_agg_filt = {k: v for k, v in global_data_dict_agg.items() if keep_df(v)}


# Flag reservoirs for consistent data in storage, inflow and outflow
def has_coverage(df, col, min_points=36, min_run=12):
    s = pd.to_datetime(df.loc[df[col].notna(), 'date'])
    p = s.dt.to_period('M').drop_duplicates().sort_values()
    if p.size < min_points:
        return False
    runs = (p != p.shift() + 1).cumsum()
    longest = runs.groupby(runs).size().max()
    return longest >= min_run

rows = []
for gid, d in global_data_dict_agg_filt.items():
    rows.append({
        'id': gid,
        'storage_cons': True,  # already filtered on storage
        'inflow_cons': has_coverage(d, 'inflow_mean'),
        'outflow_cons': has_coverage(d, 'outflow_mean'),
    })

flags_df = pd.DataFrame(rows).sort_values('id').reset_index(drop=True)

# Remove reservoirs too inconsistent in comparison to GranD storage values -> Using >50% relative error as the cutoff
grand_gdf = gpd.read_file(grand_fn)

ids_to_remove = []

for grand_id, df in global_data_dict_agg_filt.items():

    # long-term mean storage
    storage_lta = df['storage_mean'].mean()

    # capacity from grand_df
    cap_cols = ['CAP_MCM']
    cap_vals = grand_gdf.loc[grand_gdf['GRAND_ID'] == int(grand_id), cap_cols]
    cap_ref = float(cap_vals.max(axis=1).values[0])*0.001
    rel_err = abs(storage_lta - cap_ref) / cap_ref * 100

    if rel_err > 60:
        ids_to_remove.append(grand_id)


for gid in ids_to_remove:
    global_data_dict_agg_filt.pop(gid, None)

# export data
print(f"> Exporting to {export_folder}")
with open(f'{export_folder}/globalReservoirDataAggregated.pkl','wb') as f:
    pkl.dump(global_data_dict_agg_filt, f)