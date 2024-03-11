# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 13:26:54 2024

@author: rein_pp
"""

import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import fiona
import tl3_analysis_toolbox
import rasterio
from datetime import date
import pdb
from scipy.optimize import leastsq
import warnings
warnings.filterwarnings("ignore")

def mask_ds_with_shp(ds,mask,polygon,crop):
    geoms=rasterio.warp.transform_geom('EPSG:4326','EPSG:3035',[polygon['geometry']])
    for varname in list(ds.keys()):
        ds[varname]=ds[varname]*mask
    if crop==True:
        bounds=rasterio.features.bounds(geoms[0])
        ds=ds.sel(x=slice(bounds[0],bounds[2]),y=slice(bounds[3],bounds[1]))  
    return ds
'''
path_stats='/nfs/IGARSS_2022/Results_Laura/Stats/'

study_areas = {
    'Adriatic Sea': 't04',
    'Aegean Sea': 't04',
    'Balearic (Iberian Sea)': 't03'
    } 

for IHO_name in study_areas.keys():
    print('Processing '+IHO_name)
    print('------------------------------------------------')

    start_date=date(2000,1,1)
    end_date=date(2000,12,31)
    daily_range = pd.date_range(start_date, end_date, freq = 'd')
    date_range = daily_range[daily_range.day.isin([8,18,28])] 
    
    
    tl_crop = tl3_analysis_toolbox.crop()
    aux=tl3_analysis_toolbox.aux()
    
    tile=study_areas[IHO_name]
    
    shp_study_area='/nfs/IGARSS_2022/Results_Laura/Shps/final_areas.shp'
    shapefile=fiona.open(shp_study_area, "r")
    polygon=[poly for poly in shapefile if poly['properties']['NAME']==IHO_name][0]
    #tiles = tl_crop.what_tiles_for_poly(polygon)
    
    file='/nfs/IGARSS_2022/Results_Laura/latlon/'+'lonlat_l3_'+tile+'.nc'
    ds=xr.open_dataset(file)
    
    coord_upper_left=(np.array(ds.x)[0],np.array(ds.y)[0])
    resolution=1000
    epsg_poly='4326'
    epsg_mask='3035'
    mask_shape=ds.lat.shape
    
    mask=tl_crop.create_mask_from_shp(polygon,mask_shape,coord_upper_left,resolution,epsg_poly,epsg_mask)
    
    sst_maximum=[]
    sst_median=[]
    sst_min=[]
    
    months=[]
    years=[]
    days=[]
    hours=[]
    
    for year in range(1990,2023,1): 
        print('Process year '+str(year))
        for dt in date_range:
            print('process_date '+str(dt))
            if year > 2018:
                version = 't01.40'
            else: 
                version = 'v01.01'
            date_f=str(year)+str(dt.month).zfill(2)+str(dt.day).zfill(2)
            file = aux.get_l3_file_from_ida(date_f, 'Decades', 'SST', version, tile)
            if file !=None:
                ds=xr.open_dataset(file)
        
        
                #ds=mask_ds_with_shp(ds, mask, polygon, crop=True)
        
        
                sst_max=np.array(ds.sst_max)
                vt_max=np.array(ds.view_time_max)
        
                sst_max=sst_max*mask
        
                valid=np.isfinite(sst_max)
                sst_max=sst_max[valid]
                vt_max=vt_max[valid]
                
                for i in range(1,24):
                    sst_hour=sst_max[(vt_max>i)&(vt_max<i+1)]
                    if len(sst_hour)>3:
                        sst_median.append(np.quantile(sst_hour,0.5))
                        sst_maximum.append(np.quantile(sst_hour,0.9))
                        sst_min.append(np.quantile(sst_hour,0.1))
                        
                        months.append(dt.month)
                        days.append(dt.day)
                        years.append(year)
                        hours.append(i)
            
    df=pd.DataFrame({'sst_median':sst_median,'sst_maximum':sst_maximum,'sst_min':sst_min,
                     'month':months,'day':days,'year':years,'hour':hours})        
    outfile=path_stats+IHO_name+'_observation_time_stats.csv'
    df.to_csv(outfile)
'''

        
def get_residuals(vars,t,anom):
    a=vars[0]
    b=vars[1]
    sst_mod=dtc(t,a,b)
    return abs(anom-sst_mod)

def dtc(t,a,b):
    return a*np.cos((t*12.29)+4)+b




'''
for i in range(1,13):
    df_sub=df_Adriatic_stats[df_Adriatic_stats.month==i]
    df_sub['dif']=abs(df_sub.sst_median-np.median(df_sub.sst_median))
    df_sub=df_sub[df_sub.dif<3]
    fig=plt.figure()
    ax=fig.add_subplot(111)
    df_sub.boxplot(column=['sst_median'],by='hour',ax=ax)
    fig.suptitle(str(i).zfill(2))
    fig.savefig(path+'obs_time'+str(i).zfill(2)+'.png')

for i in range(1,13):
    df_sub=df_Adriatic_stats[df_Adriatic_stats.month==i]
    df_sub['dif']=abs(df_sub.sst_median-np.median(df_sub.sst_median))
    df_sub=df_sub[df_sub.dif<3]
    medians=[]
    for hour in range(8,20):
        df_hour=df_sub[df_sub.hour==hour]
        #print(np.nanmedian(df_hour.sst_median))
        medians.append(np.nanmedian(df_hour.sst_median))
    
    day_anom=np.array(medians)-np.nanmedian(medians)
    plt.plot(range(8,20),day_anom)
    plt.title(str(i))
    plt.show()
'''    
path='E:/Publications/SST_analysis/daytime_correction/'
study_areas=['Adriatic_Sea','Aegean_Sea','Balearic_Iberian_Sea']

for study_area in study_areas:
    df_Adriatic_stats=pd.read_csv(path+study_area+'_observation_time_stats.csv')
    
    a_s=[]
    b_s=[]
    c_s=[]
    d_s=[]
    i_s=[]
    for i in range(1,13):
        df_sub=df_Adriatic_stats[df_Adriatic_stats.month==i]
        df_sub=df_sub[(df_sub.hour>4)&(df_sub.hour<25)]
        df_sub['dif']=abs(df_sub.sst_median-np.median(df_sub.sst_median))
        df_sub=df_sub[df_sub.dif<3]
        medians=[]
        hours=np.array(range(1,25))
        for hour in hours:
            df_hour=df_sub[df_sub.hour==hour]
            #print(np.nanmedian(df_hour.sst_median))
            medians.append(np.nanmedian(df_hour.sst_median))
        anom=np.array(medians)-np.nanmedian(medians)
        
        valid=np.isfinite(anom)
        hours=hours[valid]
        anom=anom[valid]
        
        vars=[1,1]
        out=leastsq(get_residuals,vars,args=(hours,anom))
        
        
        a=out[0][0]
        b=out[0][1]
        #print(b)
        #print(a)
        #a=1
        #b=12.3
        #c=0
        #d=3
        
        a_s.append(abs(a))
        b_s.append(b)
        i_s.append(i)
    
        
        cycle=dtc(hours, abs(a), b)
        #plt.plot(anom)
        
        fig=plt.figure()
        ax=fig.add_subplot(111)
        df_sub.boxplot(column=['sst_median'],by='hour',ax=ax)
        ax.plot(cycle+np.nanmedian(medians))
        fig.suptitle(str(i).zfill(2))
        fig.savefig(path+study_area+'_obs_time'+str(i).zfill(2)+'.png')
        
        
    
    params=pd.DataFrame({'month':i_s,'a':a_s,'b':b_s})
    params.to_csv(path+study_area+'_daytime_params.csv')
