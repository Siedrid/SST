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
import pdb
#from fiona.model import to_dict
import fiona
import argparse

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

def mask_ds_with_shp(ds,mask,polygon,crop):
    geoms=rasterio.warp.transform_geom('EPSG:4326','EPSG:3035',[polygon['geometry']])
    for varname in list(ds.keys()):
        ds[varname]=ds[varname]*mask
    if crop==True:
        bounds=rasterio.features.bounds(geoms[0])
        ds=ds.sel(x=slice(bounds[0],bounds[2]),y=slice(bounds[3],bounds[1]))  
    return ds

def count_valid_pixel_study_area(path_in,short,stats):
    shapefile=fiona.open(shp_study_area, "r")
    #polygon=[poly for poly in shapefile if poly['properties']['NAME']==IHO_name][0]
    polygon=[poly for poly in shapefile if IHO_name[0:6] in poly['properties']['NAME']][0]
    
    file_mk=path_in+short+'/'+'01_merged_mosaic_dif_1990_'+ short + '_'+stats+'.nc' 
    monthly_mk=xr.open_dataset(file_mk)

    # Create binary mask from study area shape
    coord_upper_left=(np.array(monthly_mk.x)[0],np.array(monthly_mk.y)[0])
    resolution=1000
    epsg_poly='4326'
    epsg_mask='3035'
    mask_shape=monthly_mk.sst_dif_med.shape
    mask=tl_crop.create_mask_from_shp(polygon,mask_shape,coord_upper_left,resolution,epsg_poly,epsg_mask)
    
    count=(mask==1).sum()
    return pd.DataFrame({'short':[short],'count':[count]})

def count_trends_study_area(study_area_stats,path_mk,path_stats,short,stats):
    stats_trend=pd.DataFrame()
    for m in range(1,13):
        file_mk=path_mk+short+'/'+str(m).zfill(2)+'_merged_mosaic_mk_'+ short + '_'+stats+'.nc' 
        monthly_mk=xr.open_dataset(file_mk)
        count=int(np.isfinite(monthly_mk.p).sum())
        count_all=int(study_area_stats[study_area_stats.short==short]['count'])
        percent=count/count_all*100
        stats_m=pd.DataFrame({'short':[short], 'month':[m],'count':[count],'count_all':[count_all],
                            'percent':[percent]})
        stats_trend=stats_trend.append(stats_m)
    stats_trend.to_csv(path_stats+'stats_trend_'+short+'.csv')


def mosaic_mk(IHO_name,stats):
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
        
        dif_lst = [path_in + poly_id + '_' + str(month).zfill(2) + '_monthly_dif_'+stats+'.nc' for poly_id in poly_lst]
        mk_lst = [path_in + poly_id + '_' + str(month).zfill(2) + '_monthly_mk_'+stats+'.nc' for poly_id in poly_lst]
        
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
        
        # Mask with actual study area
        shapefile=fiona.open(shp_study_area, "r")
        #polygon=[to_dict(poly) for poly in shapefile if to_dict(poly)['properties']['NAME']==IHO_name][0]
        #polygon=[poly for poly in shapefile if poly['properties']['NAME']==IHO_name][0]
        polygon=[poly for poly in shapefile if IHO_name[0:4] in poly['properties']['NAME']][0]
        
        # Create binary mask from study area shape
        coord_upper_left=(np.array(monthly_mk.x)[0],np.array(monthly_mk.y)[0])
        resolution=1000
        epsg_poly='4326'
        epsg_mask='3035'
        mask_shape=monthly_mk.p.shape
        
        mask=tl_crop.create_mask_from_shp(polygon,mask_shape,coord_upper_left,resolution,epsg_poly,epsg_mask)
        monthly_mk=mask_ds_with_shp(monthly_mk, mask, polygon, crop=True)
        
        # Merge Monthly Anomalies
        #monthly_dif = reproj.mosaic_vars(dif_ds, ['sst_dif_max', 'obs_count'], chunksize)
        
        # Write to netCDF
        print('Write ...')
        outfile = p + str(month).zfill(2) + '_merged_mosaic_mk_' + short +'_'+stats+'.nc'          
        write_nc(monthly_mk, outfile)
        
        outfile= p +str(month).zfill(2) + '_merged_mosaic_dif_' + short + '_'+stats+'.nc'        
       # write_nc(monthly_dif, outfile)
    
        print('Closing Datasets ...')
        #monthly_dif.close()
        monthly_mk.close()
        gc.collect()
        

#%% Mosaic Anomalies for 2022

def mosaic_anomalies(IHO_name, year, shp_study_area,stats):
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
    #for month in [3]:    
        print('Mosaicing month-year ' + str(month)+'-'+str(year))
        
        dif_lst = [path_in + poly_id + '_' + str(month).zfill(2) + '_monthly_dif_'+stats+'.nc' for poly_id in poly_lst]
        
        dif_lst_sort = []
        for i in range(len(dif_lst)):
            if os.path.isfile(dif_lst[i]):
                dif_lst_sort.append(dif_lst[i])
                #mk_lst_sort.append(mk_lst[i])
            else:
                print(poly_lst[i] + ' missing in input folder')    
                
        dif_ds = [xr.open_dataset(dif_lst_sort[i]) for i in range(len(dif_lst_sort))]
        dif_ds_2022 = [dif_ds[d].sel(year = year) for d in range(len(dif_ds)) if year in dif_ds[d].year]        
        
        if len(dif_ds_2022)>0:
            for d in range(len(dif_lst_sort)):
                
                dif_ds_2022[d]['obs_count'] = dif_ds_2022[d].obs_count.where(dif_ds_2022[d]['obs_count'] != 0)
                dif_ds_2022[d]['sst_dif_med'] = dif_ds_2022[d]['sst_dif_med']
    
            print('Merging ...')
            chunksize=[50,50,50]
            
            # Merge Monthly Anomalies
            monthly_dif = reproj.mosaic_vars(dif_ds_2022, ['sst_dif_med', 'obs_count'], chunksize)
            
            # Mask with actual study area
            shapefile=fiona.open(shp_study_area, "r")
            #polygon=[to_dict(poly) for poly in shapefile if to_dict(poly)['properties']['NAME']==IHO_name][0]
            
            #polygon=[poly for poly in shapefile if poly['properties']['NAME']==IHO_name][0]
            polygon=[poly for poly in shapefile if IHO_name[0:4] in poly['properties']['NAME']][0]
            
            # Create binary mask from study area shape
            coord_upper_left=(np.array(monthly_dif.x)[0],np.array(monthly_dif.y)[0])
            resolution=1000
            epsg_poly='4326'
            epsg_mask='3035'
            mask_shape=monthly_dif.sst_dif_med.shape
            
            mask=tl_crop.create_mask_from_shp(polygon,mask_shape,coord_upper_left,resolution,epsg_poly,epsg_mask)
            monthly_dif=mask_ds_with_shp(monthly_dif, mask, polygon, crop=True)
           
            
            # Write to netCDF
            print('Write ...')        
            outfile= p + str(month).zfill(2) + '_merged_mosaic_dif_'+str(year)+'_' + short + '_'+stats+'.nc'        
            write_nc(monthly_dif, outfile)
        
            print('Closing Datasets ...')
            monthly_dif.close()
        else:
            print('There are no datasets available for '+str(month)+'-'+str(year))
        gc.collect()

def initargs():
    '''
    Argument definition, also handles user input
    '''
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    #parser.add_argument('--start', help='No help',required=False)
    #parser.add_argument('--stop', help='No help',required=False)
    parser.add_argument('--stats', help='No help',required=False)
    parser.add_argument('--study_area', help='No help',required=False)
    return parser.parse_args()   

#%% Main Workflow
if __name__ == '__main__':
    args=initargs()
    stats=args.stats
    study_area=args.study_area
    
    reproj=tl3_analysis_toolbox.reproj()
    tl_crop = tl3_analysis_toolbox.crop()
    '''
    path_in = 'E:/TIMELINE_SST/OUT/results/V2/'
    path_out = 'E:/TIMELINE_SST/OUT/Mosaics/'
    tile_list_path = "E:/TIMELINE_SST/Tile_Lists/"
    
    shp_path = 'E:/TIMELINE_SST/GIS/sst_analysis_polygons/'
    all_shp = shp_path + 'intersting_sst_analysis.shp'
    
    path_in = 'E:/Publications/SST_analysis/Results/V1/'
    path_out = 'E:/Publications/SST_analysis/Mosaics/New_Cropped/'
    tile_list_path = "E:/Publications/SST_analysis/to_process/"
    
    shp_path = 'E:/SST_Analysis/Shapes/'
    all_shp = shp_path + 'intersting_sst_analysis.shp'
    shp_study_area='E:/Publications/SST_analysis/Final_Study_Areas/final_areas.shp'
    '''
    path_in = '/nfs/IGARSS_2022/Results_Laura/Monthly_Results/New_Test/'
    path_out = '/nfs/IGARSS_2022/Results_Laura/Composites/Median/'
    #path_out = '/nfs/IGARSS_2022/Results_Laura/Composites/MK/'
    tile_list_path = "/nfs/IGARSS_2022/Results_Laura/to_process/"
    path_stats='/nfs/IGARSS_2022/Results_Laura/Stats/'
    
    shp_path = '/nfs/IGARSS_2022/Results_Laura/Shps/'
    all_shp = shp_path + 'intersting_sst_analysis.shp'
    shp_study_area=shp_path+'final_areas.shp'
    '''
    study_areas = {
        'Skagerrak': 'A',
        'Adriatic Sea': 'B',
        'Aegean Sea': 'C',
        'Balearic (Iberian Sea)': 'D'
        } 
    '''
    study_area_stats=pd.DataFrame()
    
    #for key in study_areas.keys():
    #for key in ['Adriatic Sea']:
    IHO_name = study_area 
    print('Mosaicing ' + IHO_name)
    short = short_from_IHO(IHO_name)
    poly_lst = read_tile_lst(IHO_name)
    
    '''
    p = path_out + short + '/'
    if not os.path.exists(p):
        os.makedirs(p)    
        
    mosaic_mk(IHO_name,stats)
   
    for year in range(1990,2023,1):
    #for year in [2022]:
        mosaic_anomalies(IHO_name, year,shp_study_area,stats)
    '''
    
    stats_sa=count_valid_pixel_study_area(path_out,short,stats)
    study_area_stats=study_area_stats.append(stats_sa)
    
    stats_out=path_stats+'study_area_stats_'+short+'.csv'
    study_area_stats.to_csv(stats_out)
    
    study_area_stats=pd.read_csv(stats_out)
    #path_mk=path_out+'MK/'
    count_trends_study_area(study_area_stats, path_out, path_stats,short,stats)
    
        
        
        
        
        
        
        
        