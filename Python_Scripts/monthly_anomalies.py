# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 10:29:16 2023

@author: laura
"""


import tl3_analysis_toolbox
import tl_cci_sst_explorer
#import sst_validation_2021
import xarray as xr
import matplotlib.pyplot as plt
import rioxarray
import rasterio
from rasterio.enums import Resampling
import pdb
import numpy as np
np.set_printoptions(suppress=True)
#from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import linregress
#import mpl_scatter_density
import pandas as pd
from datetime import datetime
import rioxarray # warum nochmal importiert?
from datetime import timedelta
import os
import fiona
import pymannkendall as mk
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse
from datetime import date
import warnings
import gc

hard_drive = 'E'

def write_nc(xds, fn):
    '''
    Write xarray as netcdf
    '''
    comp = dict(zlib=True, complevel=1)
    encoding = {var: comp for var in xds.data_vars}
    xds.to_netcdf(fn, encoding=encoding)
    #xds.to_netcdf(fn) 

def interpolate_along_time(arr,method,limit):
    result = np.zeros_like(arr, dtype=np.float32)
    x=range(arr.shape[1])
    y=range(arr.shape[2])
    for i in x:
        for j in y:
            line_stack = pd.DataFrame(data=arr[:,i,j], dtype=np.float32,columns=['a'])
            line_stack['a'].interpolate(method=method, axis=0, inplace=True, limit=limit)
            if line_stack.a.isna().sum() >0:
                line_stack['a']=-999
            result[:, i, j] = line_stack['a'].values.astype(np.float32)
    return result
    
def ApplyMannKendall(xds,var,spco):
    yco=spco[0]
    xco=spco[1]
    my_func_name='ApplyMannKendall'
    print(datetime.now().strftime("%H:%M:%S")+' '+my_func_name)
    xds_var=xds[var]
    arr=np.array(xds_var)
    
    #arr[arr<260]=np.nan
    #pdb.set_trace()
    #arr=nan_along_time_to_zero(arr) Doesnt work on Geofarm
    #arr=interpolate_nan(arr,method='linear',limit=1)
    
    arr=interpolate_along_time(arr, 'linear', 9)
    
    result = np.apply_along_axis(mk.original_test, 0, arr) # Default alpha of 0.05
    xds_res=xr.Dataset(coords={yco:xds.coords[yco],xco:xds.coords[xco]},
                                      attrs=xds.attrs)
    res_fields=['trend','h','p','z','Tau','s','var_s','slope','intercept']
    i=0
    for res in result:
        field=res_fields[i]
        if field=='trend':
            res[res=='no trend']='0'
            res[res=='increasing']='1'
            res[res=='decreasing']='-1'
        if field=='h':
            res[res=='False']='0'
            res[res=='True']='1'
        res=res.astype(float)
        xds_res[field]=([yco,xco],res)
        i=i+1
    return xds_res

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

def doi_to_month(day_of_year):
    if not 1 <= day_of_year <= 365:
        raise ValueError("Day of year must be between 1 and 365")

    base_date = datetime(datetime.now().year, 1, 1)
    target_date = base_date + timedelta(days=day_of_year - 1)
    
    return target_date.month

#%%

def anomaly_trends():
    print('Calculating Decadal Anomalies')
    hard_drive = 'E'
    path_out= hard_drive + ':/TIMELINE_SST/OUT/decadal_Anomaliesv2/'
    poly_lst = np.unique([os.listdir(path_out)[d].split('_')[0] for d in range(len(os.listdir(path_out)))])

    #shp='D:/TIMELINE_SST/GIS/coast_dk/coast_dk.shp'
    shp = hard_drive + ':/TIMELINE_SST/GIS/sst_analysis_polygons/intersting_sst_analysis.shp'    
    i = find_feature_index_by_id(shp, 3166)

    start_date=date(2000,1,1)
    end_date=date(2000,12,31)
   
    aux=tl3_analysis_toolbox.aux()
    tl_crop=tl3_analysis_toolbox.crop()
    prep=tl3_analysis_toolbox.l3_lst_sst_prep()
    reproj=tl3_analysis_toolbox.reproj()
    
    shapefile=fiona.open(shp, "r")
    
    for i in range(len(shapefile)):
        poly=shapefile[i]
        poly_id= int(poly['properties']['id'])
        
        if str(poly_id) in poly_lst:
            print('Decadal Anomalies for PolyID ' + str(poly_id) + ' already calculated')
        
        else:
            tiles = tl_crop.what_tiles_for_poly(poly)
            print('Processing PolyID: ' + str(poly_id))
            
            temp_res = 'Decades'
            if temp_res == 'Daily':
                date_range = pd.date_range(start_date, end_date, freq = 'd')
            if temp_res == 'Decades':
                daily_range = pd.date_range(start_date, end_date, freq = 'd')
                date_range = daily_range[daily_range.day.isin([8,18,28])] 
                # make sure that each decade file is only processed once
            
            for dt in date_range:
                print('Calculating Anomalies for: ' + str(dt))
                doi=dt.timetuple().tm_yday
                xds_list=[]
                
                for tile in tiles:
                    #xds_mns=xr.Dataset(coords={'x':[],'y':[],'doi':[]})
                    fns = []
                                  
                        
                    for year in range(1990,2023,1):            
        
                        if temp_res == 'Decades':
                            if year > 2018:
                                version = 't01.40'
                            else: 
                                version = 'v01.01'
                            
                            date_f=str(year)+str(dt.month).zfill(2)+str(dt.day).zfill(2)
                            fn = aux.get_l3_file_from_ida(date_f, temp_res, 'SST', version, tile)
                            if fn !=None:
                                fns.append(fn)
                        
                        else:
                         date_f=str(year)+str(dt.month).zfill(2)+str(dt.day).zfill(2)
                         fn=aux.get_l3_file_from_ida(date_f,temp_res,'SST','v01.01',tile)
                         if fn !=None:
                             fns.append(fn)
                
           
                    # Stack all Observations through the time axis
                    xds = xr.open_mfdataset(fns, combine='nested', concat_dim=['t'], 
                                chunks={'x': 2300, 'y': 3250},
                                preprocess=prep.add_date)
                    xds=xds.drop(['lambert_azimuthal_equal_area'])
                    
                    #if np.any(~np.isnan(xds['sst'][0].values)):
                    if xds.x.size > 0:
                        # Crop to Subpolygon Shape (1x1 Grad)
                        crop=True
                        xds=tl_crop.crop_ds_with_shp_tiles(xds,poly,crop,tile)
                        xds_list.append(xds) # append to list
                chunksize=[1000,1000,37] # wird nicht mehr benötigt
                if temp_res == 'daily':
                    var=['sst','view_time'] # keep only certain variables                
                    # calculate mean and daily anomalies
                    xds_merged= reproj.mosaic_vars(xds_list,[var],chunksize)

                    xds_merged['sst_mean']=xds_merged['sst'].mean(dim='t')
                    xds_merged['sst_dif']=xds_merged['sst']-xds_merged['sst_mean']
                    xds_dif=xds_merged[['sst_dif']] 

                    
                if temp_res == 'Decades':
                    var=['qual_max',
                         'sst_max', 'sst_std', 
                         'view_time_max', 'valid_obs_count']
                    
                    xds_merged= reproj.mosaic_vars(xds_list,[var],chunksize)
                    #xds_merged['sst_mean_dec']=xds_merged['sst_mean'].mean(dim='t')
                    xds_merged['sst_max_dec']=xds_merged['sst_max'].mean(dim='t')
                    #xds_merged['sst_min_dec']=xds_merged['sst_min'].mean(dim='t')
                    #xds_merged['sst_med_dec']=xds_merged['sst_median'].mean(dim='t')

                    #xds_merged['sst_dif_mean']=xds_merged['sst_mean']-xds_merged['sst_mean_dec']
                    #xds_merged['sst_dif_min']=xds_merged['sst_min']-xds_merged['sst_min_dec']
                    xds_merged['sst_dif_max']=xds_merged['sst_max']-xds_merged['sst_max_dec']
                    #xds_merged['sst_dif_med']=xds_merged['sst_median']-xds_merged['sst_med_dec']

                    xds_dif=xds_merged[['sst_dif_max']] 
                    
                    # Mask values with unrealistic anomalies
                    xds_dif = xds_dif.where(xds_dif['sst_dif_max'] < 4)
                    xds_dif = xds_dif.where(xds_dif['sst_dif_max'] > -4)
                    
                    # Obs Count
                    #xds_dif['obs_count']=obs_count
                    xds_dif['obs_count'] = xds_merged['valid_obs_count']
                        

                print('Write Dif Netcdf')
                outfile=path_out+str(poly_id)+'_'+str(doi).zfill(3)+'_' + temp_res +'_dif.nc'        
                write_nc(xds_dif, outfile)
                #xds_doi=xr.Dataset(coords={'x':xds_merged.coords['x'],'y':xds_merged.coords['y'], 'doi':doi})
            

                xds_dif.close()
                xds_merged.close()
                xds.close()
                
                gc.collect() # collect garbage
  

def monthly_anomaly_trends():
    print('Calculating monthly anomalies')
    path_in = hard_drive + ':/TIMELINE_SST/OUT/decadal_Anomaliesv2/'
    path_out= hard_drive + ':/TIMELINE_SST/OUT/monthly_Anomaliesv2/'
    
    dec_doi = np.unique([os.listdir(path_in)[d].split('_')[1] for d in range(len(os.listdir(path_in)))])
    poly_lst = np.unique([os.listdir(path_in)[d].split('_')[0] for d in range(len(os.listdir(path_in)))])
    poly_lst_out = np.unique([os.listdir(path_out)[d].split('_')[0] for d in range(len(os.listdir(path_out)))])
    
    aux=tl3_analysis_toolbox.aux()
    tl_crop=tl3_analysis_toolbox.crop()
    prep=tl3_analysis_toolbox.l3_lst_sst_prep()
    reproj=tl3_analysis_toolbox.reproj()
    
    for poly_id in poly_lst:
        if str(poly_id) in poly_lst_out:
            print('Monthly Anomalies for PolyID ' + str(poly_id) + ' already calculated')
            
        else:

            print('Processing PolyID: ' + poly_id)     
            doi_perMonth = [dec_doi[i:i+3] for i in range(0, len(dec_doi), 3)]
            
            for m in range(12):
                fns = [path_in + poly_id + '_' + str(doi_perMonth[m][i]).zfill(3) + '_Decades_dif.nc' for i in range(3)]
                # 0 durch m ersetzen!
                print(fns)
                
                # Stack all Observations through the time axis
                xds = xr.open_mfdataset(fns, combine='nested', concat_dim=['t'], 
                                chunks={'x': 2300, 'y': 3250})        
        
                xds_max = xds[['sst_dif_max']]
                xds_m = xds_max.groupby('t.year').mean('t')
                xds_m['obs_count'] = xds['obs_count'].groupby('t.year').sum()
                #xds_m = xds_m.rio.write_crs(3035)
        
                print('Write Dif Netcdf')
                outfile=path_out+str(poly_id)+'_'+str(m+1).zfill(2)+'_monthly_dif.nc'        
                write_nc(xds_m, outfile)
                        
                # Apply MannKenddal
                spco=['y','x']
                xds_res=ApplyMannKendall(xds_m,'sst_dif_max',spco) 
                #xds_res = xds_res.rio.write_crs(3035) # set spatial reference
                
                outfile=path_out+str(poly_id)+'_'+str(m+1).zfill(2)+'_monthly_mk.nc' 
                write_nc(xds_res, outfile) 
                
                xds.close()
                xds_res.close()
                xds_m.close()
                
                gc.collect() # collect garbage

#%%

if __name__ == '__main__':
    
    anomaly_trends()
    monthly_anomaly_trends()