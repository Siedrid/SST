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

def mask_ds_with_shp(ds,mask,polygon,crop):
    geoms=rasterio.warp.transform_geom('EPSG:4326','EPSG:3035',[polygon['geometry']])
    for varname in list(ds.keys()):
        ds[varname]=ds[varname]*mask
    if crop==True:
        bounds=rasterio.features.bounds(geoms[0])
        ds=ds.sel(x=slice(bounds[0],bounds[2]),y=slice(bounds[3],bounds[1]))  
    return ds

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


        



path='E:/Publications/SST_analysis/Test/dif/'
df_Adriatic_stats=pd.read_csv(path+'Adriatic Sea_observation_time_stats.csv')
for i in range(1,13):
    df_sub=df_Adriatic_stats[df_Adriatic_stats.month==i]
    fig=plt.figure()
    ax=fig.add_subplot(111)
    df_sub.boxplot(column=['sst_maximum'],by='hour',ax=ax)
    fig.suptitle(str(i).zfill(2))
    fig.savefig(path+'obs_time'+str(i).zfill(2)+'.png')



df_Adriatic_Sea=pd.read_csv(path+'Adriatic_SeaCCI_TL_difference.csv')



'''
path='E:/Publications/SST_analysis/Test/dif/'

df_Skagerrak=pd.read_csv(path+'SkagerrakCCI_TL_difference.csv')
df_Adriatic_Sea=pd.read_csv(path+'Adriatic_SeaCCI_TL_difference.csv')
df_Balearic_Iberian_Sea=pd.read_csv(path+'Balearic_Iberian_SeaCCI_TL_difference.csv')
df_Aegean_Sea=pd.read_csv(path+'Aegean_SeaCCI_TL_difference.csv')

file=path+'TL-L3-SST-AVHRR_NOAA-199408_decade3-fv0101_t04.nc'
'''





#plt.plot(sst_list)
#plt.boxplot(sst_list)