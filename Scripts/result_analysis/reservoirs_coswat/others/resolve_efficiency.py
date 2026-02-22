"""
Author: Jose P. Teran
Date: 2026-01-16

Count reservoirs resolved into network
"""

import os
import geopandas as gpd
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[4]  # root of repo
os.chdir(BASE_DIR)  # All relative paths will be based on this ! !

analysis_folder = f'Scripts/result_analysis'

# Regions and files
version = '1.5.0'
modelDataPath = 'CoSWAT-Framework/model-data'
modelSetupPath  = f'CoSWAT-Framework/model-setup/CoSWATv{version}'

region_list = os.listdir(modelSetupPath)

baseLake_dict       = {}
resolvedLake_dict   = {}

for region in region_list:
    baseLakeFn      = f'{modelDataPath}/{region}/shapes/lakes-grand-ESRI-54003-demAligned-forSnap.shp'
    resolvedLakeFn  = f'{modelSetupPath}/{region}/Watershed/Shapes/lakes-grand-ESRI-54003.shp'
    
    baseLake_gdf        = gpd.read_file(baseLakeFn)
    resolvedLake_gdf    = gpd.read_file(resolvedLakeFn)
    
    baseLake_dict[region]       = baseLake_gdf
    resolvedLake_dict[region]   = resolvedLake_gdf

num_resolved = 0
num_base     = 0
storage_base = 0
storage_resolved = 0
for region in region_list:
    num_base += len(baseLake_dict[region])
    num_resolved += len(resolvedLake_dict[region])
    
    storage_base += baseLake_dict[region]['smax'].astype(float).sum()*1e-9
    storage_resolved += resolvedLake_dict[region]['smax'].astype(float).sum()*1e-9

print(f"Total # of base reservoirs: {num_base}")
print(f"Total # of resolved reservoirs: {num_resolved}")
print(f"Percentage: {num_resolved/num_base*100:.2f}")

print(f"Total storage (km3) of base reservoirs: {storage_base}")
print(f"Total storage (km3) of resolved reservoirs: {storage_resolved}")
print(f"Percentage: {storage_resolved/storage_base*100:.2f}")