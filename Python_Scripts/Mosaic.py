# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 19:15:47 2023

@author: laura
"""

## Mosaicing Datasets
import os
import xarray as xr
import tl3_analysis_toolbox
import matplotlib.pyplot as plt
import numpy as np
import fiona
import geopandas as gpd
import pyproj
from rasterio.crs import CRS
import rasterio.mask
import rasterio.warp
import gc
import pandas as pd

def write_nc(xds, fn):
    '''
    Write xarray as netcdf
    '''
    comp = dict(zlib=True, complevel=1)
    encoding = {var: comp for var in xds.data_vars}
    xds.to_netcdf(fn, encoding=encoding)
    #xds.to_netcdf(fn) 

def find_feature_index_by_id(path, ID):
    '''
    
    select feature with QGIS and retrieve ID

    '''
    with fiona.open(path, 'r') as src:
        for i, feature in enumerate(src):
            if feature['properties']['id'] == ID:
                index_i = i
                return index_i
            
    return "PolyID is not in Shapefile"

# Check for missing poly_ids in Results folder
def missing_polys(all_shp, path_in):
    poly_lst_shp = set(gpd.read_file(all_shp)['id'].astype(int).astype(str).values)
    poly_lst = set(np.unique([os.listdir(path_in)[d].split('_')[0] for d in range(len(os.listdir(path_in)))]))
    
    result = list(poly_lst_shp - poly_lst)
    return result

def read_tile_lst(IHO_name):
    tile_path = tile_list_path + IHO_name + '.csv'
    tile_df = pd.read_csv(tile_path)
    tile_lst = tile_df.ID.astype(int).astype(str).values
    
    return tile_lst

def short_from_IHO(IHO_name):
    short = IHO_name.replace(" ", "_").replace("(", "").replace(")", "")
    return short

#%%

def mosaic_mk(IHO_name):
    '''
    Function mosaics 1x1° tiles processed with monthly_anomalies.py on the basis of the tile_lists csv dataframes on Drive.
    Names of the csv files refer to the official boundaries defined in the IHO dataset.

    Parameters
    ----------
    IHO_name : TYPE string
        DESCRIPTION. Official name of the IHO basins, and csv files on drive, e.g. "Adriatic Sea"

    Returns
    -------
    None. Writes mosaiced MK and Anomaly (dif) nc files to output directory. Creates one mosaic per month.

    '''
    
    # create one folder for mosaics per study area
    p = path_out + short + '/'
    if not os.path.exists(p):
        os.makedirs(p)
    
    for month in range(1,13):
        
        print('Mosaicing month ' + str(month))
        
        dif_lst = [path_in + poly_id + '_' + str(month).zfill(2) + '_monthly_dif.nc' for poly_id in poly_lst]
        mk_lst = [path_in + poly_id + '_' + str(month).zfill(2) + '_monthly_mk.nc' for poly_id in poly_lst]
        
        mk_lst_sort = []
        dif_lst_sort = []
        
        # Checking if all files exist in input folder
        for i in range(len(dif_lst)):
            if os.path.isfile(dif_lst[i]):
                dif_lst_sort.append(dif_lst[i])
                mk_lst_sort.append(mk_lst[i])
            else:
                print(poly_lst[i] + ' missing in input folder')       
        
        #dif_ds = [xr.open_dataset(dif_lst_sort[i]) for i in range(len(dif_lst_sort))]
        mk_ds = [xr.open_dataset(mk_lst_sort[i]) for i in range(len(mk_lst_sort))]
        
        # Variables to keep in Mann-Kendall .nc files
        mk_vars = ['p', 'slope', 'intercept', 'Tau', 'trend', 'h']
        
        for d in range(len(dif_lst_sort)):
            
            # mask non-significant pixel
            mask = mk_ds[d]['h'] == 1
            #dif_ds[d]['obs_count'] = dif_ds[d].obs_count.where(dif_ds[d]['obs_count'] != 0)
            #dif_ds[d]['sst_dif_max'] = dif_ds[d]['sst_dif_max'] # Possibility to select specific years, if not all years are of interest

            for var in mk_vars:           
                mk_ds[d][var] = mk_ds[d][var].where(mask)

        # Merge MannKendall
        print('Merging ...')
        chunksize=[50,50,50] # wird nicht mehr benötigt
        monthly_mk = reproj.mosaic_vars(mk_ds, mk_vars, chunksize)
        
        # Merge Monthly Anomalies
        #monthly_dif = reproj.mosaic_vars(dif_ds, ['sst_dif_max', 'obs_count'], chunksize)
        
        # Write to netCDF
        print('Write ...')
        outfile = p + str(month).zfill(2) + '_merged_mosaic_mk_' + short +'.nc'        
        write_nc(monthly_mk, outfile)
        
        outfile= p +str(month).zfill(2) + '_merged_mosaic_dif_' + short + '.nc'        
       # write_nc(monthly_dif, outfile)
    
        print('Closing Datasets ...')
        #monthly_dif.close()
        monthly_mk.close()
        gc.collect()
        

#%% Mosaic Anomalies for 2022

def mosaic_anomalies(IHO_name, year):
    '''
    Mosaics anomalies for specific years only. Otherwise same as in mosaic_ds

    Parameters
    ----------
    IHO_name : TYPE string
        DESCRIPTION.
    year : TYPE
        DESCRIPTION.

    Returns
    -------
    None. Writes mosaiced anomalies to output directory. Creates one mosaic .nc per month.

    '''
    # create one folder for mosaics per study area
    p = path_out + short + '/'
    if not os.path.exists(p):
        os.makedirs(p)
    
    for month in range(1,13):
        
        print('Mosaicing month ' + str(month))
        
        dif_lst = [path_in + poly_id + '_' + str(month).zfill(2) + '_monthly_dif.nc' for poly_id in poly_lst]
        
        dif_lst_sort = []
        for i in range(len(dif_lst)):
            if os.path.isfile(dif_lst[i]):
                dif_lst_sort.append(dif_lst[i])
                #mk_lst_sort.append(mk_lst[i])
            else:
                print(poly_lst[i] + ' missing in input folder')       
        
        dif_ds = [xr.open_dataset(dif_lst_sort[i]) for i in range(len(dif_lst_sort))]
        
        dif_ds_2022 = [dif_ds[d].sel(year = year) for d in range(len(dif_ds))]        
        
        for d in range(len(dif_lst_sort)):
            
            dif_ds_2022[d]['obs_count'] = dif_ds_2022[d].obs_count.where(dif_ds_2022[d]['obs_count'] != 0)
            dif_ds_2022[d]['sst_dif_max'] = dif_ds_2022[d]['sst_dif_max']

        print('Merging ...')
        chunksize=[50,50,50]
        
        # Merge Monthly Anomalies
        monthly_dif = reproj.mosaic_vars(dif_ds_2022, ['sst_dif_max', 'obs_count'], chunksize)
        
        # Write to netCDF
        print('Write ...')        
        outfile= p + str(month).zfill(2) + '_merged_mosaic_dif_2022_' + short + '.nc'        
        write_nc(monthly_dif, outfile)
    
        print('Closing Datasets ...')
        monthly_dif.close()
        gc.collect()


#%% Main Workflow
if __name__ == '__main__':
    
    reproj=tl3_analysis_toolbox.reproj()
    tl_crop = tl3_analysis_toolbox.crop()
    
    path_in = 'E:/TIMELINE_SST/OUT/results/V2/'
    path_out = 'E:/TIMELINE_SST/OUT/Mosaics/'
    tile_list_path = "E:/TIMELINE_SST/Tile_Lists/"
    
    shp_path = 'E:/TIMELINE_SST/GIS/sst_analysis_polygons/'
    all_shp = shp_path + 'intersting_sst_analysis.shp'
    
    study_areas = {
        'Skagerrak': 'A',
        'Adriatic Sea': 'B',
        'Aegean Sea': 'C',
        'Balearic (Iberian Sea)': 'D'
        } 
    
    for key in study_areas.keys():
        
        IHO_name = key 
        print('Mosaicing ' + IHO_name)
        short = short_from_IHO(IHO_name)
        poly_lst = read_tile_lst(IHO_name)
        
        p = path_out + short + '/'
        if not os.path.exists(p):
            os.makedirs(p)    
            
        mosaic_mk(IHO_name)
        mosaic_anomalies(IHO_name, 2022)