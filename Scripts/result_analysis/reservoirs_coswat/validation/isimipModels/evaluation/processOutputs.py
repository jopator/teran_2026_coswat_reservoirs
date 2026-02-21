"""
Author: Jose P. Teran
Date: 2026-01-18

Description:
Process outputs of reservoir storage for ISIMIP 3A models


- We will get the grand_id from the evaluated reservoirs
- Read the shapefile and get the centroid
- Sample ISIMIP models at the centroid
- If more than one reservoir in same pixel (representative water body) fraction storage

"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import pickle as pkl
import argparse
import matplotlib.pyplot as plt
import hydroeval as he
import xarray as xr
import pickle as pkl

regions = os.listdir('/media/jopato/jopato_ssd/PHD/PHD_main/Projects/CoSWAT/Paper_Ch2/codeAndDataAvailability/teran_et_al_2026_coswat_reservoirs/Scripts/result_analysis/reservoirs_coswat/coswat_outputs/reservoir_storage/CoSWATv1.5.0') #'/data/brussel/vo/000/bvo00033/vsc10883/PHD_main/Projects/CoSWAT/Scripts/result_analysis/reservoirs_coswat/coswat_outputs/reservoir_storage/CoSWATv1.5.0')

for region in regions:
    version = "1.5.0"

    print(f"> Processing region {region} ... \n")
    # Paths and file names
    # ====================
    os.chdir('/media/jopato/jopato_ssd/PHD/PHD_main/Projects/CoSWAT/Paper_Ch2/codeAndDataAvailability/teran_et_al_2026_coswat_reservoirs')#'/data/brussel/vo/000/bvo00033') # All relative paths will be based on this ! !
    analysis_folder = 'Scripts/result_analysis' #f'vsc10883/PHD_main/Projects/CoSWAT/Scripts/result_analysis'

    stat_folder     = f'{analysis_folder}/reservoirs_coswat/validation/stats/CoSWATv1.5.0/{region}'


    modelSetupPath  = f'CoSWAT-Framework/model-setup/CoSWATv1.5.0/{region}'
    wshedGisPath    = f'{modelSetupPath}/Watershed'
    reservoir_shp   = f'{wshedGisPath}/Shapes/lakes-grand-ESRI-54003.shp'
    hydroLakes_shp  = 'CoSWAT-Framework/data-preparation/resources/hydro-lakes/HydroLakes/HydroLAKES_polys_v10_fixed_v2.shp'


    modelOutput     = f'{analysis_folder}/reservoirs_coswat/validation/isimipModels/modelOutputs'

    stat_dir        = f'{analysis_folder}/reservoirs_coswat/validation/isimipModels/stats/{region}'
    plot_dir        = f'{analysis_folder}/reservoirs_coswat/validation/isimipModels/plots/{region}'

    processed_out   = f'{analysis_folder}/reservoirs_coswat/validation/isimipModels/processedData/{region}'

    os.makedirs(stat_dir,exist_ok=True)
    os.makedirs(plot_dir,exist_ok=True)
    os.makedirs(processed_out,exist_ok=True)






    # Get regions and reservoir location
    # ====================
    print("... Getting reservoir location and representative WB weight \n")


    swat_stat_df    = pd.read_csv(f'{stat_folder}/monthly_storage_stats_{region}.csv',index_col=0)
    grand_id_list   = swat_stat_df['grand_id'].to_list()
    modelFileList = os.listdir(modelOutput)

    res_gdf     = gpd.read_file(f'{hydroLakes_shp}')
    # print(res_gdf.columns)
    res_df     = res_gdf[res_gdf['Grand_id'].isin(grand_id_list)][['Hylak_id','Grand_id','Vol_total','Pour_lat','Pour_long']]
    res_df     = res_df.rename(columns={'Pour_lat':'lat','Pour_long':'lon'})

    # Weights for representative waterbodies
    res_df = res_df.copy()
    res_df["weight"] = 1.0

    ds0 = xr.open_dataset(os.path.join(modelOutput, modelFileList[0]))

    latg = ds0['lat'].values
    long = ds0['lon'].values

    lon_pts = res_df["lon"].to_numpy(float)
    lat_pts = res_df["lat"].to_numpy(float)


    # Cell id definition
    def nearest_idx_1d(grid, pts):
        grid = np.asarray(grid)
        if grid.ndim != 1:
            raise ValueError("Expected 1D grid for lat/lon.")
        return np.abs(grid[None, :] - pts[:, None]).argmin(axis=1)

    if latg.ndim == 1 and long.ndim == 1:
        ilat = nearest_idx_1d(latg, lat_pts)
        ilon = nearest_idx_1d(long, lon_pts)
        cell_id = ilat.astype(np.int64) * int(len(long)) + ilon.astype(np.int64)

    else:
        # 2D lat/lon grid case (slower but works)
        lat2 = np.asarray(latg)
        lon2 = np.asarray(long)
        flat_lat = lat2.ravel()
        flat_lon = lon2.ravel()

        # lon wrap for 2D as well (if needed)
        if np.nanmin(flat_lon) >= 0 and np.nanmax(flat_lon) > 180:
            lon_pts2 = np.where(res_df["lon"].to_numpy(float) < 0,
                                res_df["lon"].to_numpy(float) + 360,
                                res_df["lon"].to_numpy(float))
        else:
            lon_pts2 = res_df["lon"].to_numpy(float)

        lat_pts2 = res_df["lat"].to_numpy(float)

        # nearest in (lat,lon) euclidean in degrees (good enough for picking a cell)
        cell_flat = np.empty(len(res_df), dtype=np.int64)
        for i in range(len(res_df)):
            d2 = (flat_lat - lat_pts2[i])**2 + (flat_lon - lon_pts2[i])**2
            cell_flat[i] = int(d2.argmin())
        cell_id = cell_flat

    res_df["cell_id"] = cell_id

    res_df["Vol_total"] = res_df["Vol_total"].astype(float)
    smax_sum = res_df.groupby("cell_id")["Vol_total"].transform("sum")
    res_df["weight"] = res_df["Vol_total"] / smax_sum

    # helper: nearest index + safe window slices (works for 1D lat/lon)
    def nearest_i(grid_1d, x):
        g = np.asarray(grid_1d)
        return int(np.abs(g - float(x)).argmin())

    def window_slices(ilat0, ilon0, nlat, nlon, win):
        i0 = max(0, ilat0 - win)
        i1 = min(nlat - 1, ilat0 + win)
        j0 = max(0, ilon0 - win)
        j1 = min(nlon - 1, ilon0 + win)
        return slice(i0, i1 + 1), slice(j0, j1 + 1)


    # Sample isimip models
    # ====================

    model_storage_df_dict = {}
    print("... Sampling ISIMIP models \n")

    for model in modelFileList:
        model_Fn = f'{modelOutput}/{model}'
        modelName = model.split("_")[0]

        print(f"\t \t > Processing {modelName}")
        model_ds = xr.open_dataset(model_Fn)
        R = 6371000.0
        lats = model_ds["lat"]
        lons = model_ds["lon"]

        dlat = np.deg2rad(float(abs(lats.diff("lat").mean())))
        dlon = np.deg2rad(float(abs(lons.diff("lon").mean())))

        cell_area_lat = (R**2 * dlat * dlon * np.cos(np.deg2rad(lats))).rename("cell_area")

        storage_df_dict = {}

        model_ds = model_ds.sel(time=slice("1970-01-01", "2015-12-31"))

        for index, row in res_df.iterrows():
            grand_id = row['Grand_id']
            lat = row['lat']
            lon = row['lon']
            weight  = row['weight']
            da_pt   = model_ds["reservoirstor"].sel(lat=lat, lon=lon, method="nearest")*weight

            area_pt = cell_area_lat.sel(lat=da_pt["lat"], method="nearest")
            ts_km3  = (da_pt / 1000) * area_pt / 1e9

            #  skip all-zero series
            if float(np.abs(ts_km3).sum().item()) == 0:
                print(f"WARNING: all-zero series for Grand_id {grand_id}, skipping")
                continue


            df_ts   = ts_km3.to_dataframe(name="stor_km3").reset_index()
            storage_df_dict[grand_id] = df_ts

        model_storage_df_dict[model] = storage_df_dict


    out_fn = os.path.join(processed_out, f"isimip_storage_{region}.pkl")

    with open(out_fn, "wb") as f:
        pkl.dump(model_storage_df_dict, f)

    print(f"Saved to {out_fn}")