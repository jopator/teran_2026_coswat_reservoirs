"""
Author: Jose P. Teran
Date: 2026-01-18

Description:
Evaluate outputs of reservoir storage for ISIMIP 3A models against reference data

Per region:
    - We read obs and all sims and get KGE, PBIAS, and R

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



os.chdir('/media/jopato/jopato_ssd/PHD/PHD_main/Projects/CoSWAT/Paper_Ch2/codeAndDataAvailability/teran_et_al_2026_coswat_reservoirs')#'/data/brussel/vo/000/bvo00033') # All relative paths will be based on this ! !
analysis_folder = 'Scripts/result_analysis' #f'vsc10883/PHD_main/Projects/CoSWAT/Scripts/result_analysis'

regions = os.listdir('/media/jopato/jopato_ssd/PHD/PHD_main/Projects/CoSWAT/Paper_Ch2/codeAndDataAvailability/teran_et_al_2026_coswat_reservoirs/Scripts/result_analysis/reservoirs_coswat/coswat_outputs/reservoir_storage/CoSWATv1.5.0') #'/data/brussel/vo/000/bvo00033/vsc10883/PHD_main/Projects/CoSWAT/Scripts/result_analysis/reservoirs_coswat/coswat_outputs/reservoir_storage/CoSWATv1.5.0')

for region in regions:

    # Paths and file names
    # ====================

    sim_outputs_dir = f'{analysis_folder}/reservoirs_coswat/validation/isimipModels/processedData/{region}'
    sim_outputs     = f'{sim_outputs_dir}/isimip_storage_{region}.pkl'
    obs_pkl_fn      = f'{analysis_folder}/reservoirs_coswat/observations/Reservoir_storage_global_data/processed/data/globalReservoirDataAggregated.pkl'
    swatOutputs     = f'{analysis_folder}/reservoirs_coswat/coswat_outputs/reservoir_storage/CoSWATv1.5.0/{region}/monthly/reservoir_mon.pkl'

    stat_dir        = f'{analysis_folder}/reservoirs_coswat/validation/isimipModels/stats/{region}'
    plot_dir        = f'{analysis_folder}/reservoirs_coswat/validation/isimipModels/plots/{region}'

    processed_out   = f'{analysis_folder}/reservoirs_coswat/validation/isimipModels/processedData/{region}'



    # Read pickles
    with open(obs_pkl_fn,'rb') as obsPkl:
        obs_dict = pkl.load(obsPkl)
        
    with open(sim_outputs,'rb') as simMonPkl:
        simMon_dict = pkl.load(simMonPkl)

    with open(swatOutputs,'rb') as swatMonPkl:
        swatMon_dict = pkl.load(swatMonPkl)


    # Evaluation
    # ====================
    print(f"... Starting evaluation for {region} \n")

    # Get grand_id list from the model with the highest number of results
    best_model = max(simMon_dict, key=lambda m: len(simMon_dict[m].keys()))
    grand_id_list = [int(x) for x in simMon_dict[best_model].keys()]

    # Filter out observations for this reservoirs
    grand_dict = {}
    obsRegion_dict = {k: obs_dict[k] for k in grand_id_list if k in obs_dict}

    # Get grand_id per lake id in CoSWAT
    grand_dict = {}
    swat_grand_id_list = []
    for k,item in swatMon_dict.items(): 
        swat_grand_id_list.append(int(item['grand_id']))
        grand_dict[item['grand_id']] = int(k)

    # Calculate stats and plot
    isimipColorDict = {'h08':'blueviolet','cwatm':'olivedrab','lpjml5-7-10-fire':'tomato','miroc-integ-land':'yellowgreen','watergap2-2e':'crimson'}
    stats_rows = []

    for grand_id in grand_id_list:
        obs_df = obsRegion_dict[grand_id].copy()
        obs_df["date"] = pd.to_datetime(obs_df["date"])
        swat_df = swatMon_dict[grand_dict[grand_id]]['coswatOut']
        swat_df['flo_stor_km3']  = swat_df['flo_stor']/1000000000

        xmin = max(obs_df["date"].min(), pd.Timestamp("1970-01-01"))
        xmax = min(obs_df["date"].max(), pd.Timestamp("2016-01-01"))

        model_data_dict = {}
        for model in simMon_dict:
            modelName = model.split('_')[0]
            try:
                model_df = simMon_dict[model][float(grand_id)].copy()
                model_df["time"] = pd.to_datetime(model_df["time"].astype(str))
                model_data_dict[modelName] = model_df[["time", "stor_km3"]]
            except:
                continue

        fig, ax = plt.subplots(figsize=(10, 5))

        obs_df.plot(x="date", y="storage_mean",color='k',linewidth=1.5, label="Reference", ax=ax)
        swat_df.plot(x="date", y="flo_stor_km3",color='navy',linestyle='--',linewidth=1.5, label="CoSWAT", ax=ax)

        for modelName, dfm in model_data_dict.items():
            dfm.plot(x="time", y="stor_km3", linewidth = 0.5,color=isimipColorDict[modelName],label=modelName, ax=ax)

        ax.set_ylabel("Storage (km3)")

        ax.set_xlim(xmin, xmax)

        fig.tight_layout()
        fig.savefig(f"{plot_dir}/isimip_storage_{int(grand_id)}.jpg", dpi=300, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)



        # Stats
        obs_stats = obs_df.loc[(obs_df["date"] >= xmin) & (obs_df["date"] <= xmax), ["date", "storage_mean"]].copy()
        obs_stats["ym"] = obs_stats["date"].dt.to_period("M")

        row = {"grand_id": int(grand_id)}

        # CoSWAT
        swat_stats = swat_df.loc[(swat_df["date"] >= xmin) & (swat_df["date"] <= xmax), ["date", "flo_stor_km3"]].copy()
        swat_stats["date"] = pd.to_datetime(swat_stats["date"])
        swat_stats["ym"] = swat_stats["date"].dt.to_period("M")

        m = obs_stats.merge(swat_stats[["ym", "flo_stor_km3"]], on="ym", how="inner")
        if not m.empty:
            kge, pbias, r, _, _ = calcSats(m, "flo_stor_km3", "storage_mean")
            row["kge_CoSWAT"] = kge
            row["pbias_CoSWAT"] = pbias
            row["r_CoSWAT"] = r

        # ISIMIP
        for modelName, dfm in model_data_dict.items():
            dfm_stats = dfm.rename(columns={"time": "date"}).copy()
            dfm_stats["date"] = pd.to_datetime(dfm_stats["date"])
            dfm_stats = dfm_stats.loc[(dfm_stats["date"] >= xmin) & (dfm_stats["date"] <= xmax), ["date", "stor_km3"]]
            dfm_stats["ym"] = dfm_stats["date"].dt.to_period("M")

            m = obs_stats.merge(dfm_stats[["ym", "stor_km3"]], on="ym", how="inner")
            if m.empty:
                continue

            kge, pbias, r, _, _ = calcSats(m, "stor_km3", "storage_mean")
            row[f"kge_{modelName}"] = kge
            row[f"pbias_{modelName}"] = pbias
            row[f"r_{modelName}"] = r

        stats_rows.append(row)


    stats_df = pd.DataFrame(stats_rows)

    cols = ["grand_id"] + \
        sorted([c for c in stats_df.columns if c.startswith("kge_")]) + \
        sorted([c for c in stats_df.columns if c.startswith("pbias_")]) + \
        sorted([c for c in stats_df.columns if c.startswith("r_")])

    stats_df = stats_df.reindex(columns=cols)
    stats_df.to_csv(f"{stat_dir}/storage_stats_by_reservoir.csv", index=False)