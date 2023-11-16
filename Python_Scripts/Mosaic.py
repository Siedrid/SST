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

def missing_polys():
    int_sst_analysis = hard_drive + shp_path +'intersting_sst_analysis.shp'
    poly_lst_shp = set(gpd.read_file(int_sst_analysis)['id'].astype(int).astype(str).values)
    poly_lst = set(np.unique([os.listdir(path_in)[d].split('_')[0] for d in range(len(os.listdir(path_in)))]))
    
    result = list(poly_lst_shp - poly_lst)
    return result

#%%
reproj=tl3_analysis_toolbox.reproj()
tl_crop = tl3_analysis_toolbox.crop()
hard_drive = 'E'

path_in = 'E:/TIMELINE_SST/OUT/results/V2/'
plt_path = 'E:/TIMELINE_SST/OUT/Plots/mosaics/'
path_out = 'E:/TIMELINE_SST/OUT/Mosaics/'

shp_path = ':/TIMELINE_SST/GIS/sst_analysis_polygons/'
at_shp = hard_drive + ':/TIMELINE_SST/GIS/sst_analysis_polygons/atlantic_polys_t03.shp'    
all_shp = hard_drive + shp_path + 'intersting_sst_analysis.shp'
gib_shp = hard_drive + ':/TIMELINE_SST/GIS/sst_analysis_polygons/gibralta_polys_t03.shp'
baltic_shp = hard_drive + shp_path + 'baltic_sea.shp'
north_africa_shp = hard_drive + shp_path + 'norther_africa.shp'
GB_shp = hard_drive + shp_path + 'GB.shp'
east_euro_shp = hard_drive + shp_path + 'E-Euro.shp'
italy_shp = hard_drive + shp_path + 'italy.shp'
spain_shp = hard_drive + shp_path + 'spain.shp'#
france_shp = hard_drive + shp_path + 'french_atlantic.shp'
malle_shp = hard_drive + shp_path + 'medi_malle.shp'


# poly_lst = np.unique([os.listdir(path_in)[d].split('_')[0] for d in range(len(os.listdir(path_in)))])
# all polyids in input path

# Funktion, um all_shp zu IHO basin zu croppen --> get poly_ids

def mosaic_ds(shp_path, short):
    poly_lst = gpd.read_file(shp_path)['id'].astype(int).astype(str).values
    
    for month in range(1,13):
        
        print('Mosaicing month ' + str(month))
        
        dif_lst = [path_in + poly_id + '_' + str(month).zfill(2) + '_monthly_dif.nc' for poly_id in poly_lst]
        mk_lst = [path_in + poly_id + '_' + str(month).zfill(2) + '_monthly_mk.nc' for poly_id in poly_lst]
        
        mk_lst_sort = []
        dif_lst_sort = []
        for i in range(len(dif_lst)):
            if os.path.isfile(dif_lst[i]):
                dif_lst_sort.append(dif_lst[i])
                mk_lst_sort.append(mk_lst[i])
            else:
                print(poly_lst[i] + ' missing in input folder')       
        
        dif_ds = [xr.open_dataset(dif_lst_sort[i]) for i in range(len(dif_lst_sort))]
        mk_ds = [xr.open_dataset(mk_lst_sort[i]) for i in range(len(mk_lst_sort))]
        
        mk_vars = ['p', 'slope', 'intercept', 'Tau', 'trend', 'h']
        
        for d in range(len(dif_lst_sort)):
            
            # mask non-significant pixel
            mask = mk_ds[d]['h'] == 1
            dif_ds[d]['obs_count'] = dif_ds[d].obs_count.where(dif_ds[d]['obs_count'] != 0)
            dif_ds[d]['sst_dif_max'] = (dif_ds[d]['sst_dif_max']*1000).astype(int)

            for var in mk_vars:           
                mk_ds[d][var] = mk_ds[d][var].where(mask)

        # Merge MannKendall
        print('Merging ...')
        chunksize=[50,50,50] # wird nicht mehr benötigt
        monthly_mk = reproj.mosaic_vars(mk_ds, mk_vars, chunksize)
        
        # Merge Monthly Anomalies
        monthly_dif = reproj.mosaic_vars(dif_ds, ['sst_dif_max', 'obs_count'], chunksize)
        monthly_dif = xr.merge(dif_ds)
        
        # Write to netCDF
        print('Write ...')
        outfile = path_out+short + '/' + str(month).zfill(2) + '_merged_mosaic_mk_v2_' + short +'.nc'        
        write_nc(monthly_mk, outfile)
        
        outfile=path_out+ short + '/' +str(month).zfill(2) + '_merged_mosaic_dif_v2_' + short + '.nc'        
        write_nc(monthly_dif, outfile)
    
        print('Closing Datasets ...')
        monthly_dif.close()
        monthly_mk.close()
        gc.collect()
        


mosaic_ds(north_africa_shp, 'northern_africa')
mosaic_ds(italy_shp, 'italy')
mosaic_ds(east_euro_shp, 'eastern_europe')
mosaic_ds(baltic_shp, 'baltic')
mosaic_ds(spain_shp, 'spain')
mosaic_ds(GB_shp, 'GB')
mosaic_ds(france_shp, 'french_atlantic')
mosaic_ds(malle_shp, 'malle')

#%% Mosaic Anomalies for 2022

def mosaic_anomalies(shp_path, short, year):
    poly_lst = gpd.read_file(shp_path)['id'].astype(int).astype(str).values
    
    for month in range(1,13):
        
        print('Mosaicing month ' + str(month))
        
        dif_lst = [path_in + poly_id + '_' + str(month).zfill(2) + '_monthly_dif.nc' for poly_id in poly_lst]
        #mk_lst = [path_in + poly_id + '_' + str(month).zfill(2) + '_monthly_mk.nc' for poly_id in poly_lst]
        
        #mk_lst_sort = []
        dif_lst_sort = []
        for i in range(len(dif_lst)):
            if os.path.isfile(dif_lst[i]):
                dif_lst_sort.append(dif_lst[i])
                #mk_lst_sort.append(mk_lst[i])
            else:
                print(poly_lst[i] + ' missing in input folder')       
        
        dif_ds = [xr.open_dataset(dif_lst_sort[i]) for i in range(len(dif_lst_sort))]
        
        dif_ds_2022 = [dif_ds[d].sel(year = year) for d in range(len(dif_ds))]        
        #mk_vars = ['p', 'slope', 'intercept', 'Tau', 'trend', 'h']
        
        for d in range(len(dif_lst_sort)):
            
            # mask not significant pixel
            #mask = mk_ds[d]['h'] == 1
            dif_ds_2022[d]['obs_count'] = dif_ds_2022[d].obs_count.where(dif_ds_2022[d]['obs_count'] != 0)
            dif_ds_2022[d]['sst_dif_max'] = dif_ds_2022[d]['sst_dif_max']
            '''
            for var in mk_vars:           
                mk_ds[d][var] = mk_ds[d][var].where(mask)
            '''
        # Merge MannKendall
        print('Merging ...')
        chunksize=[50,50,50] # wird nicht mehr benötigt
        #monthly_mk = reproj.mosaic_vars(mk_ds, mk_vars, chunksize)
        
        # Merge Monthly Anomalies
        monthly_dif = reproj.mosaic_vars(dif_ds_2022, ['sst_dif_max', 'obs_count'], chunksize)
        #monthly_dif = xr.merge(dif_ds)
        
        # Write to netCDF
        print('Write ...')
        #outfile = path_out+short + '/' + str(month).zfill(2) + '_merged_mosaic_mk_' + short +'.nc'        
        #write_nc(monthly_mk, outfile)
        
        outfile=path_out+ short + '/' +str(month).zfill(2) + '_merged_mosaic_dif_2022_v3' + short + '.nc'        
        write_nc(monthly_dif, outfile)
    
        print('Closing Datasets ...')
        monthly_dif.close()
        #monthly_mk.close()
        gc.collect()

mosaic_anomalies(malle_shp, 'malle', 2022)
mosaic_anomalies(italy_shp, 'italy', 2022)
mosaic_anomalies(east_euro_shp, 'eastern_europe', 2022)
mosaic_anomalies(baltic_shp, 'baltic', 2022)
mosaic_anomalies(north_africa_shp, 'northern_africa', 2022)
mosaic_anomalies(spain_shp, 'spain', 2022)
mosaic_anomalies(GB_shp, 'GB', 2022)
mosaic_anomalies(france_shp, 'french_atlantic', 2022)