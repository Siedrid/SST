# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 19:02:56 2023

@author: laura
"""
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
import rioxarray
import gc
from datetime import datetime
import pandas as pd
from sklearn.linear_model import LinearRegression

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

def read_tile_lst(IHO_name):
    tile_path = tile_list_path + IHO_name + '.csv'
    tile_df = pd.read_csv(tile_path)
    tile_lst = tile_df.ID.astype(int).astype(str).values
    
    return tile_lst

def short_from_IHO(IHO_name):
    short = IHO_name.replace(" ", "_").replace("(", "").replace(")", "")
    return short

# Check for missing poly_ids in Results folder
def missing_polys(all_shp, path_in):
    poly_lst_shp = set(gpd.read_file(all_shp)['id'].astype(int).astype(str).values)
    poly_lst = set(np.unique([os.listdir(path_in)[d].split('_')[0] for d in range(len(os.listdir(path_in)))]))
    
    result = list(poly_lst_shp - poly_lst)
    return result

def get_cci_shp(all_shp, IHO_name, crop2IHO = False):
    
    poly_lst = read_tile_lst(IHO_name)    
    
    gdf = gpd.read_file(all_shp)
    gdf['id'] = gdf['id'].astype(int).astype(str)
    
    shapefile = gdf[gdf['id'].isin(poly_lst)]
    
    dis_shp = shapefile.dissolve()
    geometry = dis_shp.geometry[0]
    
    if crop2IHO == True:
        ocean_shp = gpd.read_file(ocean_path)
        sel_ocean = ocean_shp.loc[ocean_shp['NAME'] == IHO_name]
        ov = gpd.overlay(sel_ocean,dis_shp, how='intersection')
        geometry = ov.geometry[0]
    
    return geometry, poly_lst

def write_sst_ts(df, short):
    sorted_df = df.sort_values(by='Date')

    outfile = path_out + short + '/' + short + '_sst_cci-v2.csv'
    sorted_df.to_csv(outfile, index = False)
    
    return sorted_df

def load_sst_ts(short):

    df = pd.read_csv(path_out + short + '/' + short + '_sst_cci-2022.csv')
    df['datetime']=df.Date.apply(lambda x: datetime(int(x[0:4]),int(x[5:7]),int(x[8:10])))
    df['med_sst'] = df.med_sst.interpolate()
    df = df.set_index(df.datetime)
    return df

def get_site_letter(IHO_name):
    if IHO_name not in study_areas.keys():
        print("IHO_name not in study_areas")
    else:
        letter = study_areas[IHO_name]
        return letter

#%%

def sst_timeseries(IHO_name):
    
    geometry, poly_lst = get_cci_shp(all_shp, IHO_name, crop2IHO=False)
   
    dates = []
    med_sst = []
    obs_count = []    
    med_sst_cci = []
    
    for month in range(1,13):
        
        print('Mosaicing month ' + str(month))
        
        dif_lst = [path_in + poly_id + '_' + str(month).zfill(2) + '_monthly_dif.nc' for poly_id in poly_lst]        
        dif_lst_sort = []
        
        for i in range(len(dif_lst)):
            if os.path.isfile(dif_lst[i]):
                dif_lst_sort.append(dif_lst[i])
            else:
                print(poly_lst[i] + ' missing in input folder')       
        
        dif_ds = [xr.open_dataset(dif_lst_sort[i]) for i in range(len(dif_lst_sort))]

        for y in range(1990, 2023):
            
            # Median per Polygon
            med = [dif_ds[d].sst_dif_max.sel(year = y).median().values for d in range(len(dif_lst_sort)) if y in dif_ds[d].year]
            # Median over all Polygons
            med_sst.append(np.nanmedian(med))
            
            # Sum observations
            obs = [dif_ds[d].obs_count.sel(year = y).sum().values for d in range(len(dif_lst_sort)) if y in dif_ds[d].year]
            obs_count.append(sum(obs))
            
            dates.append(datetime.strptime(str(y)+ '-' + str(month).zfill(2)+'-01', "%Y-%m-%d"))

            if y < 2017:                
                print('Compositing '+str(y)+ '-' + str(month).zfill(2))
                fs_cci=[path_cci + f for f in os.listdir(path_cci) if (f[4:6]==str(month).zfill(2)) and (f[0:4]==str(y)) and ('CDR2.1_anomaly-v02.0' in f)]
                ds_cci=xr.open_mfdataset(fs_cci,chunks={'x': 2300, 'y': 3250},preprocess=prep.add_date)
                
                xds_cci=xr.Dataset(coords={'lat':ds_cci.coords['lat'],'lon':ds_cci.coords['lon']},
                                       attrs=ds_cci.attrs)
                
                # Calculate monthly median CCI composite from daily data
                xds_cci['sst_cci_anom']=ds_cci['analysed_sst_anomaly'].mean(dim='time')
                xds_cci.coords['t'] = y
        
                xds_cci = xds_cci.rio.write_crs("epsg:4326", inplace = True)
                print('Crop to Polygon')
                xds_cci_cropped = xds_cci.rio.clip([geometry])
                med_sst_cci.append(np.nanmedian(xds_cci_cropped.sst_cci_anom))
                                
            if y >= 2017:
                med_sst_cci.append(np.nan)
                
    df = pd.DataFrame({'Date': dates, 'tl_sst': med_sst, 'cci_sst': med_sst_cci, 'obs': obs_count})
    
    return df
        
#%% load csv Data

def get_linear_trend(data):
    
    # TIMELINE
    arr=np.array(data['roll_mean_tl'].fillna(method='bfill').fillna(method='ffill'))
    rows=np.arange(len(data))
    model=LinearRegression()
    model.fit(rows.reshape(-1,1),arr.reshape(-1,1))
    linear=model.predict(rows.reshape(-1,1))
    data['lin_tr_tl']=linear.reshape(len(linear))
    yr_tr_tl=float(model.coef_*12)
    
    # CCI
    arr=np.array(data['roll_mean_cci'].fillna(method='bfill').fillna(method='ffill'))
    rows=np.arange(len(data))
    model=LinearRegression()
    model.fit(rows.reshape(-1,1),arr.reshape(-1,1))
    linear=model.predict(rows.reshape(-1,1))
    data['lin_tr_cci']=linear.reshape(len(linear))
    yr_tr_cci=float(model.coef_*12)
    
    return data, yr_tr_tl, yr_tr_cci

# Plot anomalies
def plot_anomalies(df, site_letter):
    
    df = load_sst_ts(short)
    
    # Calculate Rolling Mean
    df['roll_mean_tl'] = df['med_sst'].rolling(12,min_periods=1).mean()
    df['roll_mean_cci'] = df['cci_sst'].rolling(12,min_periods=1).mean()

    df, yr_tr_tl, yr_tr_cci = get_linear_trend(df)

    fig, ax = plt.subplots() # figsize = ()
    fig.set_figheight(10)
    fig.set_figwidth(30) 
    
    ax.axhline(y = 0, color = 'black')
    ax.plot(df.datetime, df.roll_mean_tl,color='black',linewidth=6,label='TIMELINE 1-year moving mean')
    ax.plot(df.datetime, df.roll_mean_cci,color='green',linewidth=6,label='CCI 1-year moving mean')
    ax.plot(df.datetime, df.lin_tr_tl, color='black',linewidth=2,linestyle='--', label = 'TL Linear Trend')
    ax.plot(df.datetime, df.lin_tr_cci, color = 'grey', linewidth=2, linestyle = '--', label = 'CCI Linear Trend')
    
    fontsize=40
    ax.tick_params(labelsize=fontsize)
    ax.set_xlabel('Year',fontsize=fontsize)
    ax.set_ylabel('SST anomalies [K]',fontsize=fontsize)
    ax.set_ylim(-3,3)
    ax.legend(fontsize = fontsize -10, markerscale = 1.2, loc = 'lower right')
    ax.set_title('SST anomalies for ' + site_letter, fontsize = fontsize)
    ax.grid()
    fig.tight_layout()
    fig.savefig(plt_out + short + '_CCI_TL_anomalies.png')
    

# Plot TL anomalies only
def plot_TL_anomalies(df, site_letter):
    
    # Calculate Rolling Mean
    df['roll_mean_tl'] = df['med_sst'].rolling(12,min_periods=1).mean()

    df, yr_tr_tl, yr_tr_cci = get_linear_trend(df)

    fig, ax = plt.subplots() # figsize = ()
    fig.set_figheight(10)
    fig.set_figwidth(30)
    
    data_pos=df[df['med_sst']>0]
    data_neg=df[df['med_sst']<0]    
    
    ax.bar(data_pos.datetime,data_pos['med_sst'],width=30,color='red',alpha=0.7,label='TIMELINE anomalies (positive)')
    ax.bar(data_neg.datetime,data_neg['med_sst'],width=30,color='blue',alpha=0.7,label='TIMELINE anomalies (negative)')
    
    ax.plot(df.datetime, df.roll_mean_tl,color='black',linewidth=6,label='TIMELINE 1-year moving mean')
    ax.plot(df.datetime, df.lin_tr_tl, color='black',linewidth=2,linestyle='--')
    ax.text(df.datetime[10], 2.5, 'Linear Trend ('+str(round(yr_tr_tl,2))+' K/a)', fontsize = 40)
    
    fontsize=40
    ax.tick_params(labelsize=fontsize)
    ax.set_xlabel('Year',fontsize=fontsize)
    ax.set_ylabel('SST anomalies [K]',fontsize=fontsize)
    ax.set_ylim(-3,3)
    ax.legend(fontsize = fontsize -10, markerscale = 1.2, loc = 'lower right')
    ax.set_title('SST anomalies for ' + site_letter, fontsize = fontsize)
    ax.grid()
    fig.tight_layout()
    fig.savefig(plt_out + short + '_TL_anomalies.png')

#%%% Main Workflow
if __name__ == '__main__':
    
    reproj=tl3_analysis_toolbox.reproj()
    tl_crop = tl3_analysis_toolbox.crop()
    prep=tl3_analysis_toolbox.l3_lst_sst_prep()
    
    path_in = 'E:/TIMELINE_SST/OUT/results/V2/'
    path_cci= 'E:/TIMELINE_SST/CCI/'
    plt_path = 'E:/TIMELINE_SST/OUT/Plots/mosaics/'
    path_out = 'E:/TIMELINE_SST/OUT/Mosaics/Validation/'
    
    tile_list_path = "E:/TIMELINE_SST/Tile_Lists/"
    
    ocean_path = "E:/TIMELINE_SST/GIS/World_Seas_IHO_v3/"
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
        short = short_from_IHO(IHO_name)
        
        plt_out = path_out + short + '/'
        if not os.path.exists(plt_out):
            os.makedirs(plt_out)    
        
        df = sst_timeseries(IHO_name)
        sorted_df = write_sst_ts(df, short)
        
        plot_anomalies(sorted_df, get_site_letter(IHO_name))
        plot_TL_anomalies(sorted_df, get_site_letter(IHO_name))


  