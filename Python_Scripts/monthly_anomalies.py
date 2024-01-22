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
def get_from_files(path_missing):
    files_missing=[path_missing+f for f in os.listdir(path_missing)]
    big_df=pd.DataFrame()
    for file in files_missing:
        big_df=big_df.append(pd.read_csv(file))
    return np.array(big_df[big_df.existing=='no']['ID']).astype(int).astype(str)

def get_poly_id_from_file(study_area):
    path_missing='/nfs/IGARSS_2022/Results_Laura/to_process/'
    data=pd.read_csv(path_missing+study_area+'.csv')
    return np.array(data['ID']).astype(int)
   
    
  
def mask_out_pf(xds,pf):
    print('Mask out platform '+str(pf))
    xds['pf_max']=((xds['platform'])/10000).astype(int)
    xds=xds.where(xds.pf_max!=pf)
    return xds
    
  
def anomaly_trends(stats,poly_id):
    print('Calculating Decadal Anomalies')
    '''
    hard_drive = 'E'
    path_out= hard_drive + ':/TIMELINE_SST/OUT/decadal_Anomaliesv2/'
    '''
   
    
    
    path_out= '/nfs/IGARSS_2022/Results_Laura/Results/New_Test/'
    poly_lst = np.unique([os.listdir(path_out)[d].split('_')[0] for d in range(len(os.listdir(path_out)))])
    path_missing='/nfs/IGARSS_2022/Results_Laura/to_process/'
    
    #shp='D:/TIMELINE_SST/GIS/coast_dk/coast_dk.shp'
    #shp = hard_drive + ':/TIMELINE_SST/GIS/sst_analysis_polygons/intersting_sst_analysis.shp'    
    shp = '/nfs/IGARSS_2022/Results_Laura/Shps/grid.shp' 
    #i = find_feature_index_by_id(shp, 3166)
   
    start_date=date(2000,1,1)
    end_date=date(2000,12,31)
   
    aux=tl3_analysis_toolbox.aux()
    tl_crop=tl3_analysis_toolbox.crop()
    prep=tl3_analysis_toolbox.l3_lst_sst_prep()
    reproj=tl3_analysis_toolbox.reproj()
    
    shapefile=fiona.open(shp, "r")
    '''
    missing_poly_lst = ['2927', '2767', '3152', '2863', '2913', '2874', '3350', '2814', '2779', '2830', '3656', '2921', '2922', '2924', 
                        '3912', '2825', '2573', '2961', '2778', '3907', '2468', '2875', '2768', '2816', '3669', '4547', '3349', '2928', '2878',
                        '2826', '2879', '2815', '2780', '2864', '2764', '2781', '2829', '2813', '2827', '2777', '2862', '2923']
    '''
    #missing_poly_lst=get_from_files(path_missing)
    #missing_poly_lst=['3608']
    #for i in range(len(shapefile)):
    
    poly=[poly for poly in shapefile if poly['properties']['id']==int(poly_id)][0]
    #poly_id= int(poly['properties']['id'])
    '''
    if str(poly_id) in poly_lst:
        print('Decadal Anomalies for PolyID ' + str(poly_id) + ' already calculated')
    
   
    if str(poly_id) not in missing_poly_lst:
        print('Decadal Anomalies for PolyID ' + str(poly_id) + ' already calculated')
    
    else:
    '''
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
        if len(tiles)>0:
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
                    
                    #xds=tl_crop.crop_ds_with_shp_tiles(xds,poly,crop,tile)
                    # Create binary mask from study area shape
                    coord_upper_left=(np.array(xds.x)[0],np.array(xds.y)[0])
                    resolution=1000
                    epsg_poly='4326'
                    epsg_mask='3035'
                    mask_shape=(xds.y.shape[0],xds.x.shape[0])
                    
                    mask=tl_crop.create_mask_from_shp(poly,mask_shape,coord_upper_left,resolution,epsg_poly,epsg_mask)
                    xds=tl_crop.mask_ds_with_shp(xds, mask, poly, crop=True)
                    
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
                    '''
                    var=['qual_max',
                         'sst_max', 'sst_std', 
                         'view_time_max', 'valid_obs_count','platform']
                    '''
                    var=['qual_'+stats,
                         'sst_'+stats, 'sst_std', 
                         'view_time_'+stats, 'valid_obs_count','platform']
                    xds_merged= reproj.mosaic_vars(xds_list,[var],chunksize)
                    
                    # Mask out NOAA-11
                    #xds_merged=mask_out_pf(xds_merged, 11)
                    # Mask out early and late observation times
                    xds_merged=xds_merged.where((xds_merged['view_time_'+stats]>10)&(xds_merged['view_time_'+stats]<15))
                    
                    #xds_merged['sst_mean_dec']=xds_merged['sst_mean'].mean(dim='t')
                    #xds_merged['sst_max_dec']=xds_merged['sst_max'].mean(dim='t')
                    #xds_merged['sst_min_dec']=xds_merged['sst_min'].mean(dim='t')
                    xds_merged['sst_med_dec']=xds_merged['sst_'+stats].mean(dim='t')

                    #xds_merged['sst_dif_mean']=xds_merged['sst_mean']-xds_merged['sst_mean_dec']
                    #xds_merged['sst_dif_min']=xds_merged['sst_min']-xds_merged['sst_min_dec']
                    #xds_merged['sst_dif_max']=xds_merged['sst_max']-xds_merged['sst_max_dec']
                    xds_merged['sst_dif_med']=xds_merged['sst_'+stats]-xds_merged['sst_med_dec']
                    xds_dif=xds_merged[['sst_dif_med']]
                    #xds_dif=xds_merged[['sst_dif_max']] 
                    
                    # Mask values with unrealistic anomalies
                    #xds_dif = xds_dif.where(xds_dif['sst_dif_max'] < 4)
                    #xds_dif = xds_dif.where(xds_dif['sst_dif_max'] > -4)
                    xds_dif = xds_dif.where(xds_dif['sst_dif_med'] < 4)
                    xds_dif = xds_dif.where(xds_dif['sst_dif_med'] > -4)
                    
                    # Obs Count
                    #xds_dif['obs_count']=obs_count
                    xds_dif['obs_count'] = xds_merged['valid_obs_count']
                        

                print('Write Dif Netcdf')
                #outfile=path_out+str(poly_id)+'_'+str(doi).zfill(3)+'_' + temp_res +'_dif.nc' 
                outfile=path_out+str(poly_id)+'_'+str(doi).zfill(3)+'_' + temp_res +'_dif_'+stats+'.nc' 
                write_nc(xds_dif, outfile)
                #xds_doi=xr.Dataset(coords={'x':xds_merged.coords['x'],'y':xds_merged.coords['y'], 'doi':doi})
            

                xds_dif.close()
                xds_merged.close()
                xds.close()
                
                gc.collect() # collect garbage
        else:
            print('No tiles found for PolyID ' + str(poly_id))
  

def monthly_anomaly_trends(stats,poly_id):
    print('Calculating monthly anomalies')
    '''
    path_in = hard_drive + ':/TIMELINE_SST/OUT/decadal_Anomaliesv2/'
    path_out= hard_drive + ':/TIMELINE_SST/OUT/monthly_Anomaliesv2/'
    '''
    path_in = '/nfs/IGARSS_2022/Results_Laura/Results/New_Test/'
    path_out= '/nfs/IGARSS_2022/Results_Laura/Monthly_Results/New_Test/'
    dec_doi = np.unique([os.listdir(path_in)[d].split('_')[1] for d in range(len(os.listdir(path_in)))])
    poly_lst = np.unique([os.listdir(path_in)[d].split('_')[0] for d in range(len(os.listdir(path_in)))])
    poly_lst_out = np.unique([os.listdir(path_out)[d].split('_')[0] for d in range(len(os.listdir(path_out)))])
    poly_lst_out=[]
    '''
    aux=tl3_analysis_toolbox.aux()
    tl_crop=tl3_analysis_toolbox.crop()
    prep=tl3_analysis_toolbox.l3_lst_sst_prep()
    reproj=tl3_analysis_toolbox.reproj()
    poly_lst=['3608']
    '''
        
    if str(poly_id) in poly_lst_out:
        print('Monthly Anomalies for PolyID ' + str(poly_id) + ' already calculated')
        
    else:

        print('Processing PolyID: ' + str(poly_id))     
        doi_perMonth = [dec_doi[i:i+3] for i in range(0, len(dec_doi), 3)]
        
        for m in range(12):
            #fns = [path_in + poly_id + '_' + str(doi_perMonth[m][i]).zfill(3) + '_Decades_dif.nc' for i in range(3)]
            fns = [path_in + str(poly_id) + '_' + str(doi_perMonth[m][i]).zfill(3) + '_Decades_dif_'+stats+'.nc' for i in range(3)]
            # 0 durch m ersetzen!
            print(fns)
            
            # Stack all Observations through the time axis
            xds = xr.open_mfdataset(fns, combine='nested', concat_dim=['t'], 
                            chunks={'x': 2300, 'y': 3250})        
    
            #xds_max = xds[['sst_dif_max']]
            xds_max = xds[['sst_dif_med']]
            xds_m = xds_max.groupby('t.year').mean('t')
            xds_m['obs_count'] = xds['obs_count'].groupby('t.year').sum()
            #xds_m = xds_m.rio.write_crs(3035)
    
            print('Write Dif Netcdf')
            outfile=path_out+str(poly_id)+'_'+str(m+1).zfill(2)+'_monthly_dif_'+stats+'.nc'        
            write_nc(xds_m, outfile)
                    
            # Apply MannKenddal
            spco=['y','x']
            #xds_res=ApplyMannKendall(xds_m,'sst_dif_max',spco) 
            xds_res=ApplyMannKendall(xds_m,'sst_dif_med',spco) 
            #xds_res = xds_res.rio.write_crs(3035) # set spatial reference
            
            outfile=path_out+str(poly_id)+'_'+str(m+1).zfill(2)+'_monthly_mk_'+stats+'.nc' 
            write_nc(xds_res, outfile) 
            
            xds.close()
            xds_res.close()
            xds_m.close()
            
            gc.collect() # collect garbage

#%%

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

if __name__ == '__main__':
    args=initargs()
    stats=args.stats
    study_area=args.study_area
    poly_ids=get_poly_id_from_file(study_area)
    poly_ids=['2875']
    for poly_id in poly_ids:
        anomaly_trends(stats,poly_id)
        monthly_anomaly_trends(stats,poly_id)