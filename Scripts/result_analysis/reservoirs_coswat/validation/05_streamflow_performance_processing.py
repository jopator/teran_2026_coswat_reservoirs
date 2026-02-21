"""
------------------------------------
Author: Jose P. Teran
Date: 2025-12-03
------------------------------------

1. Create general dataframe with performance and scores
2. Connect with reservoir and relate with it's performance
3. Generate a gpkg and csv with all that data
"""

import os
import sys
import numpy as np
import pandas as pd
import geopandas as gpd


os.chdir('/media/jopato/jopato_ssd/PHD/PHD_main/Projects/CoSWAT/Paper_Ch2/codeAndDataAvailability/teran_et_al_2026_coswat_reservoirs')#'/data/brussel/vo/000/bvo00033') # All relative paths will be based on this ! !


# Paths settings and file names
# ====================
versionRes         = "1.5.0"
versionNores       = "1.1.0"

analysis_folder = f'Scripts/result_analysis'#f'vsc10883/PHD_main/Projects/CoSWAT/Scripts/result_analysis'
modelSetupPath  = 'CoSWAT-Framework/model-setup/CoSWATv{version}/{region}'
modelDataPath   = 'CoSWAT-Framework/model-data/{region}'
wshedGisPath    = f'{modelSetupPath}/Watershed'
reservoir_shp   = f'{wshedGisPath}/Shapes/lakes-grand-ESRI-54003.shp'
wshed_shp       = f'{wshedGisPath}/Shapes/subsNoLakes.shp'
stat_dir        = '{analysis_folder}/reservoirs_coswat/validation/stats/CoSWATv{version}'
out_dir         = f'{analysis_folder}/reservoirs_coswat/validation/analysis'
stat_plots_dir  = f'{analysis_folder}/reservoirs_coswat/validation/analysis/plots'
grdc_shp        = f'{modelDataPath}/shapes/grdc_stations-ESRI-54003.gpkg'
outlet_shp      = f'{wshedGisPath}/Shapes/outlets_sel_snap.shp'
rivs_shp        = f'{wshedGisPath}/Shapes/rivs1.shp'
chan_shp        = f'{wshedGisPath}/Shapes/dem-aster-ESRI-54003-lakeBurntchannel.shp' 
outletsSnapCsv  = '{analysis_folder}/reservoirs_coswat/coswat_outputs/streamflow/CoSWATv{version}/{region}/monthly/eval_stations.csv'
monInflowCsvFn  = 'monthly_inflow_stats_{region}.csv'
monStorageCsvFn = 'monthly_storage_stats_{region}.csv'
monOutflowCsvFn = 'monthly_outflow_stats_{region}.csv'


monSflowCsv  = '{region}/monthly_streamflow_stats_{region}.csv'

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

# First, create geodataframe with performance and skill scores
#=============================================================
print('> Merging all performance data into one GPKG and one CSV ...')
# Read all csvs into one
dfRes_list   = []
dfNoRes_list = []

for region in regions:
    statRes_df      = pd.read_csv(f'{stat_dir.format(analysis_folder=analysis_folder,version=versionRes)}/{monSflowCsv.format(region=region)}',index_col=0)
    statNoRes_df    = pd.read_csv(f'{stat_dir.format(analysis_folder=analysis_folder,version=versionNores)}/{monSflowCsv.format(region=region)}',index_col=0)
    
    statRes_df['region']    = region
    statNoRes_df['region']  = region
    
    statRes_df      = statRes_df.rename(columns={'grand_id':'grdc_no'})
    statNoRes_df    = statNoRes_df.rename(columns={'grand_id':'grdc_no'})
    
    dfRes_list.append(statRes_df)
    dfNoRes_list.append(statNoRes_df)
    
df_stats_Res    = pd.concat(dfRes_list).reset_index(drop=True)
df_stats_NoRes  = pd.concat(dfNoRes_list).reset_index(drop=True)
df_stats_NoRes = df_stats_NoRes[df_stats_NoRes['grdc_no'].isin(df_stats_Res['grdc_no'].to_list())]

# Read grdc locations
grdc_points_list = []

for region in regions:
    grdc_gdf = gpd.read_file(f'{grdc_shp}'.format(region=region))
    grdc_gdf = grdc_gdf[grdc_gdf['grdc_no'].isin(df_stats_Res.reset_index()['grdc_no'].to_list())]
    grdc_points_list.append(grdc_gdf)
    
grdc_gdf = pd.concat(grdc_points_list).drop_duplicates().reset_index(drop=True)

# Rename columns to create a general dataset
df_stats_Res    = df_stats_Res[['grdc_no','region','res_inf','kge','R','pbias']].rename(columns = {'kge':'kgeqRes','R':'rqRes','pbias':'pbiasqRes'})
df_stats_NoRes  = df_stats_NoRes[['grdc_no','kge','R','pbias']].rename(columns = {'kge':'kgeqNoRes','R':'rqNoRes','pbias':'pbiasqNoRes'})

df_stats = pd.merge(df_stats_Res,df_stats_NoRes,on='grdc_no',how='left')

# Filter out unrealistic results
df_stats  = df_stats[df_stats['kgeqRes']>-10]

# Get geometries
grdc_gdf  = grdc_gdf[grdc_gdf['grdc_no'].isin(df_stats['grdc_no'].to_list())].reset_index(drop=True)

stats_gdf = pd.merge(grdc_gdf,df_stats,on='grdc_no',how='right').reset_index(drop=True)
stats_gdf = stats_gdf[['grdc_no', 'area', 'altitude', 'region','res_inf', 
                      'kgeqRes', 'rqRes', 'pbiasqRes', 'kgeqNoRes', 'rqNoRes','pbiasqNoRes',
                      'wmo_reg', 'sub_reg', 'river', 'station', 'country', 'lat','long','geometry']]

# Connect to reservoirs
#=============================================================
# First upstream reservoir found will be assigned to the station

print("> Asessing relation with reservoirs ...")
stationsForAll = []
for region in regions:
    res_gdf             = gpd.read_file(reservoir_shp.format(region=region,version=versionRes))
    chan_gdf            = gpd.read_file(chan_shp.format(region=region,version=versionRes))
    rivs_gdf            = gpd.read_file(rivs_shp.format(region=region,version=versionRes))
    outlets_snapped_df  = pd.read_csv(outletsSnapCsv.format(region=region,version=versionRes,analysis_folder=analysis_folder),index_col=0)    

    stationsForThisRegion = []
    for idx, row in res_gdf.iterrows():
        main_use        = row['MAIN_USE']
        lake_id         = row['LakeId']
        lake_gdf        = res_gdf[res_gdf['LakeId']==lake_id]                         # get lakes
        lakeOutlet_gdf  = chan_gdf[chan_gdf['LakeMain']==lake_id]                     # get main outlet
        lakeOutlet_gdf  = lakeOutlet_gdf.reset_index(drop=True)
        complete        = False
        grand_id        = row['Grand_id']
        hylak_id        = row['Hylak_id']
        outletLINKNO    = lakeOutlet_gdf['LINKNO'].iloc[0]
        outletDSLINKNO  = lakeOutlet_gdf['DSLINKNO'].iloc[0]
        outletLENGTH    = lakeOutlet_gdf['Length'].iloc[0]

        total_length = outletLENGTH
        

        currentLINKNO = outletDSLINKNO

        links = []
        links.append(outletLINKNO)
        
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
        
        reservoir_down_chans = chan_gdf[chan_gdf['LINKNO'].isin(links)]
        stationsForThisRes   = outlets_snapped_df[outlets_snapped_df['ID'].isin(reservoir_down_chans['DSNODEID'].to_list())]
        stationsForThisRes = stationsForThisRes.copy()
        stationsForThisRes['grand_id'] = grand_id
        stationsForThisRes = stationsForThisRes.copy()
        stationsForThisRes['hylak_id'] = hylak_id
        stationsForThisRes = stationsForThisRes.copy()
        stationsForThisRes['main_use'] = main_use
        
        
        
        stationsForThisRegion.append(stationsForThisRes)
    
    stationsForThisRegion_df = pd.concat(stationsForThisRegion).reset_index(drop=True)
    stationsForThisRegion_df = stationsForThisRegion_df.drop_duplicates()
    stationsForThisRegion_df = stationsForThisRegion_df.copy()
    stationsForThisRegion_df['region'] = region
        
    stationsForAll.append(stationsForThisRegion_df)


stationsIdentified_df = pd.concat(stationsForAll).reset_index(drop=True)

# We have now all stations and what reservoir are they influenced by (might be more than one)
# Now for those, we will merge the performance data
stationsIdentified_qPerf_df = pd.merge(stationsIdentified_df[['grdc_no','grand_id','hylak_id','main_use']],stats_gdf,on='grdc_no',how='left')

# Now, we will get reservoir performance data and put it in the table
# All csvs into one
df_list_s = []
df_list_i = []
df_list_o = []

for region in regions:
    df_s = pd.read_csv(f'{stat_dir.format(analysis_folder=analysis_folder,version=versionRes)}/{region}/{monStorageCsvFn.format(region=region)}',index_col=0)
    df_list_s.append(df_s[['grand_id','kge','pbias','R']])
    try:
        df_i = pd.read_csv(f'{stat_dir.format(analysis_folder=analysis_folder,version=versionRes)}/{region}/{monInflowCsvFn.format(region=region)}',index_col=0)
        df_o = pd.read_csv(f'{stat_dir.format(analysis_folder=analysis_folder,version=versionRes)}/{region}/{monOutflowCsvFn.format(region=region)}',index_col=0)
        
        df_list_i.append(df_i[['grand_id','kge','pbias','R']])
        df_list_o.append(df_o[['grand_id','kge','pbias','R']])
    except:
        continue
    
df_gral_s = pd.concat(df_list_s).reset_index(drop=True).rename(columns = {'kge':'kgesRes','R':'rsRes','pbias':'pbiassRes'})
df_gral_i = pd.concat(df_list_i).reset_index(drop=True).rename(columns = {'kge':'kgeiRes','R':'riRes','pbias':'pbiasiRes'})
df_gral_o = pd.concat(df_list_o).reset_index(drop=True).rename(columns = {'kge':'kgeoRes','R':'roRes','pbias':'pbiasoRes'})

df_gral_s = df_gral_s.merge(df_gral_i,on='grand_id',how='left')
df_gral_s = df_gral_s.merge(df_gral_o,on='grand_id',how='left')

# Merge alles
stationsIdentified_allPerf_df = pd.merge(stationsIdentified_qPerf_df,df_gral_s,on='grand_id',how='left').reset_index(drop=True)

# Export two tables:
# - One with general streamflow stats (as csv and gpkg)
# - The other with relation to reservoirs
#=========================
print("> Saving files ...")
stats_gdf.to_file(f'{analysis_folder}/reservoirs_coswat/validation/stats/streamflowqPerf_general.gpkg')
stats_gdf.drop(columns=['geometry']).to_csv(f'{analysis_folder}/reservoirs_coswat/validation/stats/streamflowqPerf_general.csv')
stationsIdentified_allPerf_df.to_csv(f'{analysis_folder}/reservoirs_coswat/validation/stats/streamflowqPerf_ResInfluence.csv')