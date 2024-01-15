# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 09:38:18 2022

@author: rein_pp
"""
from datetime import date, datetime, timedelta
import os
from osgeo import gdal
import pdb
import xarray as xr
import numpy as np
from rioxarray.merge import merge_datasets
import rasterio
from rasterio.transform import from_origin
import pandas as pd
import dask
import matplotlib.dates as mdates
from scipy.optimize import curve_fit
from sklearn import preprocessing
import math
#from sunrise import sun
import rioxarray
from rasterio.enums import Resampling


class aux:
    #def get_l3_file_from_ida(self,date,temporal,product,version,tile):
    def get_l3_file_from_ida(self,date,temporal,product,version,tile):    
        path='/Timeline/IDA-mount/archive/AVHRR/L3/'
        v_string=version[0:3]+version[4:6]
        if temporal=='Daily' or temporal=='daily':
            temporal_str='daily'
            filename='TL-L3-'+product+'-AVHRR_NOAA-'+date+'_'+temporal_str+'-f'+v_string+'_'+tile+'.nc'
        
        if temporal=='Decades' or  temporal=='decadal':
            if int(date[6:8])<10:
              temporal_str='decade1' 
            elif int(date[6:8])>19:
              temporal_str='decade3'
            else:
              temporal_str='decade2'  
            
            filename='TL-L3-'+product+'-AVHRR_NOAA-'+date[0:6]+'_'+temporal_str+'-f'+v_string+'_'+tile+'.nc'
        
        if temporal=='Monthly' or temporal=="monthly":
            temporal_str='monthly'
            filename='TL-L3-'+product+'-AVHRR_NOAA-'+date[0:6]+'_'+temporal_str+'-f'+v_string+'_'+tile+'.nc'
            if product=='NDVI':
                temporal=temporal_str
        file=os.path.join(path, product, version,date[0:4],date[4:6],temporal,filename)
        if os.path.isfile(file):
            return file
        else:
            print('Warning: '+date+' '+product+' '+temporal+' does not exist.')
            return None
        
    def xds_to_df(self,xds,varlist):
        df=pd.DataFrame()
        for var in varlist:
            arr=np.array(xds[var]).reshape(np.array(xds[var]).size)
            df[var]=arr
        return df
        
class reproj:
    def to_tl_l3_tiles(file,resampleAlg,outdir):
        my_func_name='to_tl_l3_tiles'
        print(datetime.now().strftime("%H:%M:%S")+' '+my_func_name)
        '''
        Warps given gdal dataset to resolution and extents of TL Tiles and writes output in 
        given dir.
        '''
        '''
        tile_exts={
            "t01":(900000,3200000,4150000,5500000),
            "t02":(4150000,3200000,7400000,5500000),
            "t03":(900000,900000,4150000,3200000),
            "t04":(4150000,900000,7400000,3200000)
            }
        '''
        # Contrary to TIMELINE NC gdal.warp interpretes coordinates as pixel border not center.
        # Therefore coordinates need to be moved about 500m
        tile_exts={
            "t01":(899500,3200500,4149500,5500500),
            "t02":(4149500,3200500,7399500,5500500),
            "t03":(899500,900500,4149500,3200500),
            "t04":(4149500,900500,7399500,3200500)
            }
              
        ds = gdal.Open(file)
        '''
        # 2017 metadatum and projection is missing
        d =ds.GetMetadata_Dict()
        #if not 'AREA_OR_POINT' in d:
        ds.SetMetadata({'AREA_OR_POINT': 'Area'})
        ds.SetProjection('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,'+
                         'AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],'+
                         'UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",'+
                         'NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]')
        #pdb.set_trace()
        '''
        dsReprj = gdal.Warp("", ds, format="vrt", dstSRS = "EPSG:3035")
        
        fn=os.path.basename(file)
        tiles=['t01','t02','t03','t04']
        for tile in tiles:
            outfile=outdir+'/'+os.path.splitext(fn)[0]+'_'+tile+'.tif'
            pdb.set_trace()
            dsRes = gdal.Warp(outfile, dsReprj, xRes = 1000, yRes = 1000, 
                              resampleAlg=resampleAlg,outputBounds=tile_exts[tile])
            dsRes = None
            
        ds = dsReprj = None
        

    def mosaic_ds(self,xds_list,chunksize):
        my_func_name='mosaic_ds'
        print(datetime.now().strftime("%H:%M:%S")+' '+my_func_name)
        xds_new=[]
        for xds in xds_list:
            xds=xds.rio.write_crs(3035)
            #xds=xds.drop('lambert_azimuthal_equal_area')
            #xds=xds.chunk({'x': chunksize[0], 'y': chunksize[1],'t':chunksize[2]})
            xds=xds.chunk({'x': chunksize[0], 'y': chunksize[1]})
            xds.coords['x']=xds.coords['x'].astype('int64')
            xds.coords['y']=xds.coords['y'].astype('int64')
            xds_new.append(xds)
        #merged=merge_datasets(xds_new)
        with dask.config.set(**{'array.slicing.split_large_chunks': True}):
            merged=xr.merge(xds_new)
        merged=merged.reindex(y=list(reversed(merged.y)))
        return merged
    
    def mosaic_vars(self,xds_list,vars,chunksize):
        #mosaic_vars(self,xds_list,vars,chunksize):
        my_func_name='mosaic_vars'
        print(datetime.now().strftime("%H:%M:%S")+' '+my_func_name)
        xds_new=[]
        for var in vars:           
            for xds in xds_list:
                xds=xds.rio.write_crs(3035)
                #xds=xds.drop('lambert_azimuthal_equal_area')
                #xds=xds.chunk({'x': chunksize[0], 'y': chunksize[1],'t':chunksize[2]})
                xds.coords['x']=xds.coords['x'].astype('int64')
                xds.coords['y']=xds.coords['y'].astype('int64')
                xds_new.append(xds[var])
        #pdb.set_trace()
        merged=xr.merge(xds_new)
        merged=merged.reindex(y=list(reversed(merged.y)))
        return merged
    
    def coord_4326_to_3035(self,lat,lon):
        geom={'type':'MultiPoint',
              'coordinates':[(lon,lat),(lon,lat)]}
        trans=rasterio.warp.transform_geom('EPSG:4326','EPSG:3035',geom)
        lon=trans['coordinates'][0][0]
        lat=trans['coordinates'][0][1]
        return lat,lon
    
    def warp_togeth_cci_tl(self,xds_tl,xds_cci,prep):
        # Reproject TIMELINE dataset to WGS84
        xds_tl==prep.int64_to_int32(xds_tl) # GDAL cant handle int32 /:E
        xds_lonlat = xds_tl.rio.reproject("EPSG:4326")
        xmin=float(min(xds_lonlat['x']))
        xmax=float(max(xds_lonlat['x']))
        ymin=float(min(xds_lonlat['y']))
        ymax=float(max(xds_lonlat['y']))
        # Crop CCI dataset to TIMELINE extent
        xds_cci=xds_cci.sel(lon=slice(xmin,xmax),lat=slice(ymin,ymax)) 
        xds_cci=xds_cci.reindex(lat=list(reversed(xds_cci.lat)))
        # Resample CCI Dataset to TIMELINE resoulution
        if 'lat_bnds' in xds_cci.keys():
            xds_cci=xds_cci.drop('lat_bnds')
            xds_cci=xds_cci.drop('lon_bnds')
            xds_cci=xds_cci.drop('time_bnds')
        xds_cci=xds_cci.rio.write_crs("EPSG:4326")
        if 'sst_dtime' in xds_cci.keys(): # Rio cannot deal with timedelta variables in the L3 products
            xds_cci=xds_cci.drop(['sst_dtime','sst_depth_dtime'])
        xds_cci=xds_cci.rio.reproject(xds_cci.rio.crs,shape=(xds_lonlat.rio.height,
                                                          xds_lonlat.rio.width),
                                      resampling=Resampling.bilinear)
        #xds_merged=xds_tl.merge(xds_cci)
        for var in list(xds_cci.variables):
            if len(xds_cci[var].shape)==3:
                xds_lonlat[var]=(["y","x"],np.array(xds_cci[var])[0,:,:])
            if len(xds_cci[var].shape)==2:
                xds_lonlat[var]=(["y","x"],np.array(xds_cci[var]))
        return xds_lonlat
    
        
        

class l3_lst_sst_prep:

    def add_date(self,ds):
        #add_date(self,ds): 
        if 'platform' in list(ds.coords.keys()):
           ds=ds.reset_coords('platform')
        starttime=ds.attrs['time_coverage_start']
        date=datetime(int(starttime[0:4]),int(starttime[4:6]),int(starttime[6:8]))
        ds.coords['t']=date
        return ds
    
    def add_date_and_pf(self,ds):
        #add_date(self,ds): 
        if 'platform' in list(ds.coords.keys()):
           ds=ds.reset_coords('platform')
        starttime=ds.attrs['time_coverage_start']
        date=datetime(int(starttime[0:4]),int(starttime[4:6]),int(starttime[6:8]))
        ds.coords['t']=date
        ds.coords['pf1_day']=ds.lct1_day.attrs['platform']
        ds.coords['pf2_day']=ds.lct2_day.attrs['platform']
        ds.coords['pf3_day']=ds.lct3_day.attrs['platform']
        ds.coords['pf1_night']=ds.lct1_night.attrs['platform']
        ds.coords['pf2_night']=ds.lct2_night.attrs['platform']
        ds.coords['pf3_night']=ds.lct3_night.attrs['platform']
        return ds
    
    def l3_lst_sst_prepare(self,xds):
        my_func_name='l3_lst_sst_prep'
        print(datetime.now().strftime("%H:%M:%S")+' '+my_func_name)
        xds['pf_min']=xds['platform']%100
        xds['pf_median']=np.round((xds['platform']%10000)/100)
        xds['pf_median']=xds.where(xds['pf_median'] != 99)['pf_median']
        xds['pf_max']=np.round((xds['platform'])/10000)
        
       
        xds['view_day_min']=xds['view_day']%1000
        xds['view_day_median']=np.round((xds['view_day']%1000000)/1000)
        xds['view_day_median']=xds.where(xds['view_day_median'] != 999)['view_day_median']  
        xds['view_day_max']=np.round((xds['view_day'])/1000000)
        
        xds['view_time_median']=xds.where(xds['view_time_median'] != 99)['view_time_median']
        xds=xds.drop(['platform','view_day','lambert_azimuthal_equal_area'])
        return xds
    
    def split_pf(self,pf):
        # Splits platform variable into min, max, median
        pf_min=pf%100
        pf_median=((pf%10000)/100).astype(int)
        pf_max=((pf)/10000).astype(int)
        return pf_min,pf_median,pf_max
    
    def split_view_day(self,view_day):
        # Splits platform variable into min, max, median
        view_day_min=view_day%1000
        view_day_median=np.round((view_day%1000000)/1000)
        view_day_max=np.round((view_day)/1000000)
        return view_day_min,view_day_median,view_day_max
  
        
        
    
    def unpackbits(self,x, num_bits):
        if np.issubdtype(x.dtype, np.floating):
            raise ValueError("numpy data type needs to be int-like")
        xshape = list(x.shape)
        x = x.reshape([-1, 1])
        mask = 2**np.arange(num_bits, dtype=x.dtype).reshape([1, num_bits])
        #mask[0]=mask[0][::-1] #change endian
        return (x & mask).astype(bool).astype(int).reshape(xshape + [num_bits])

    def flag_bit_pos(xds,qual,pos,val):
        # returns new dataset with masked variables according to flags
        xds_new=xds.copy(deep=True)
        bit_arr=unpackbits(qual, 16)
        mask=np.zeros(qual.shape)
       
        mask[np.where((bit_arr[:,:,pos[0]]==val[0])&(bit_arr[:,:,pos[1]]==val[1]))]=1
        #pdb.set_trace()
        varlist=list(xds.keys())
        varlist.remove('lambert_azimuthal_equal_area')
        for var in varlist:
            xds_new[var] = xr.where(mask==1, np.nan, xds_new[var])
        return xds_new

    def filter_quality(xds,parameter,qual_param,level):
        # level 1: Just filter out worst quality (value==11)
        # level 2: Filter out worst, and reduced quality (value==11,10)
        # level 3: Only keep highest quality. Filter (value==11,10,01)
        if parameter=='sat_zenith':
            pos=(0,1)
        if parameter=='sat_azimuth':
            pos=(2,3)
        if parameter=='sun_zenith':
            pos=(4,5)
        if parameter=='sun_azimuth':
            pos=(6,7)
        if parameter=='tcwv':
            pos=(8,9)
        if parameter=='obs_count':
            pos=(10,11)
        if parameter=='qual_input':
            pos=(12,13)
                   
        qual=np.array(xds[qual_param])
        qual=qual.reshape(2300,3250)
        qual[~np.isfinite(qual)]=0
        qual=qual.astype(int)
       
        if level==1:
           xds_filtered=flag_bit_pos(xds,qual,pos,(1,1))
        if level==2:
           xds_filtered=flag_bit_pos(xds,qual,pos,(1,1))
           xds_filtered=flag_bit_pos(xds_filtered,qual,pos,(1,0))
        if level==3:
           xds_filtered=flag_bit_pos(xds,qual,pos,(1,1))
           xds_filtered=flag_bit_pos(xds_filtered,qual,pos,(1,0))
           xds_filtered=flag_bit_pos(xds_filtered,qual,pos,(0,1))
        return xds_filtered
    
    def get_qual_sz_tcwv(self,xds,var,prep):
        qual=np.array(xds[var])
        qual[~np.isfinite(qual)]=0
        qual=qual.astype('int32')
        bit_arr=prep.unpackbits(qual, 16)
        if len(bit_arr.shape)==4: # If xds is multitemporal
            tcwv_qual=bit_arr[:,:,:,8]+bit_arr[:,:,:,9]*2
            xds['tcwv_qual']=(['t','y','x'],tcwv_qual)
            sz_qual=bit_arr[:,:,:,0]+bit_arr[:,:,:,1]*2
            xds['sz_qual']=(['t','y','x'],sz_qual)
        if len(bit_arr.shape)==3: # If xds is monotemporal
            tcwv_qual=bit_arr[:,:,8]+bit_arr[:,:,9]*2
            xds['tcwv_qual']=(['y','x'],tcwv_qual)
            sz_qual=bit_arr[:,:,0]+bit_arr[:,:,1]*2
            xds['sz_qual']=(['y','x'],sz_qual)
        return xds

    def mask_tcwv_sz_qual(self,xds):
        xds=xds.where(xds['tcwv_qual']<3)
        xds=xds.where(xds['sz_qual']<3)
        return xds
    
    
    def int64_to_int32(self,xds):
        for i in xds.keys():
            if xds[i].dtype=='int64':
                xds[i]=xds[i].astype('int32')
        for i in xds.coords.keys():
            if xds.coords[i].dtype=='int64':
                xds.coords[i]=xds.coords[i].astype('int32')
        return xds
    
class crop:
    
    def create_mask_from_shp(self,polygon,mask_shape,coord_upper_left,resolution,epsg_poly,epsg_mask):
        geoms=rasterio.warp.transform_geom('EPSG:'+epsg_poly,'EPSG:'+epsg_mask,[polygon['geometry']])
        transform=from_origin(coord_upper_left[0],coord_upper_left[1] , resolution, resolution) 
        arr = np.full((mask_shape[0], mask_shape[1]),fill_value=1, dtype="uint8") # Mask template 
        # Create in memory rasterio ds
        memfile=rasterio.io.MemoryFile()
        dataset=memfile.open(
            driver='GTiff',
            height=arr.shape[0],
            width=arr.shape[1],
            count=1,
            dtype=arr.dtype,
            transform=transform,
            crs=rasterio.crs.CRS.from_epsg(int(epsg_mask))
        )
        dataset.write(arr,1)
        output, m_transform = rasterio.mask.mask(dataset, geoms)
        dataset.close()
        output=output.astype(float)
        output[output==0]=np.nan
        return output.squeeze(0)
    
    def create_mask_from_shp_tiles(self,geoms,ds):
        my_func_name='create_mask_from_shp_tiles'
        print(datetime.now().strftime("%H:%M:%S")+' '+my_func_name)
        tile_corn={
            "t01":(899500,5500500),
            "t02":(4149500,5500500),
            "t03":(899500,3200500),
            "t04":(4149500,3200500)
            }
        pdb.set_trace()
        transform=from_origin(tile_corn[tl][0],tile_corn[tl][1] , 1000, 1000) # TL Coordinates upper left and pixel size (meters)
        arr = np.full((2300, 3250),fill_value=1, dtype="uint8") # Mask template shape of the TL extent
        # Create in memory rasterio ds
        memfile=rasterio.io.MemoryFile()
        dataset=memfile.open(
            driver='GTiff',
            height=arr.shape[0],
            width=arr.shape[1],
            count=1,
            dtype=arr.dtype,
            transform=transform,
            crs=rasterio.crs.CRS.from_epsg(3035)
        )
        dataset.write(arr,1)
        output, m_transform = rasterio.mask.mask(dataset, geoms)
        dataset.close()
        output=output.astype(float)
        output[output==0]=np.nan
        return output.squeeze(0)
    
    def crop_ds_with_shp_tiles(self,ds,shp,crop,tl):
        my_func_name='crop_ds_with_shp_tiles'
        print(datetime.now().strftime("%H:%M:%S")+' '+my_func_name)
        geoms=rasterio.warp.transform_geom('EPSG:4326','EPSG:3035',[shp['geometry']])
        pdb.set_trace()
        mask=self.create_mask_from_shp_tiles(geoms, ds)
        for varname in list(ds.keys()):
            ds[varname]=ds[varname]*mask
        if crop==True:
            bounds=rasterio.features.bounds(geoms[0])
            ds=ds.sel(x=slice(bounds[0],bounds[2]),y=slice(bounds[3],bounds[1]))  
        return ds
    '''
    def create_mask_from_shp(self,geoms):   
        my_func_name='create_mask_from_shp'
        print(datetime.now().strftime("%H:%M:%S")+' '+my_func_name)
        transform=from_origin(900000,5500000 , 1000, 1000) # TL Coordinates upper left and pixel size (meters)
        arr = np.full((4600,6500),fill_value=1) # Mask template shape of the TL extent
        # Create in memory rasterio ds
        memfile=rasterio.io.MemoryFile()
        dataset=memfile.open(
            driver='GTiff',
            height=arr.shape[0],
            width=arr.shape[1],
            count=1,
            dtype=arr.dtype,
            transform=transform,
            crs=rasterio.crs.CRS.from_epsg(3035)
        )
        dataset.write(arr,1)
        output, m_transform = rasterio.mask.mask(dataset, geoms)
        dataset.close()
        return output.squeeze(0)
    '''
    def crop_ds_with_shp(self,ds,shp,crop):
        my_func_name='crop_ds_with_shp'
        print(datetime.now().strftime("%H:%M:%S")+' '+my_func_name)
        geoms=rasterio.warp.transform_geom('EPSG:4326','EPSG:3035',[shp['geometry']])
        mask=self.create_mask_from_shp(geoms)
        for varname in list(ds.keys()):
            ds[varname]=ds[varname]*mask
        if crop==True:
            bounds=rasterio.features.bounds(geoms[0])
            ds=ds.sel(x=slice(bounds[0],bounds[2]),y=slice(bounds[3],bounds[1]))  
        return ds
    
    def create_mask_from_shp_cci(self,geoms):   
        my_func_name='create_mask_from_shp_cci'
        print(datetime.now().strftime("%H:%M:%S")+' '+my_func_name)
        transform=from_origin(-180,90 , 0.05, 0.05) # CCI Coordinates upper left and pixel size (°)
        arr = np.full((3600,7200),fill_value=1,dtype='int32') # Mask template shape of the CCI extent
        # Create in memory rasterio ds
        memfile=rasterio.io.MemoryFile()
        dataset=memfile.open(
            driver='GTiff',
            height=arr.shape[0],
            width=arr.shape[1],
            count=1,
            dtype=arr.dtype,
            transform=transform,
            crs=rasterio.crs.CRS.from_epsg(3035)
        )
        dataset.write(arr,1)
        output, m_transform = rasterio.mask.mask(dataset, geoms)
        dataset.close()
        return output.squeeze(0)
    
    def crop_ds_with_shp_cci(self,ds,shp,crop):
        my_func_name='crop_ds_with_shp_cci'
        print(datetime.now().strftime("%H:%M:%S")+' '+my_func_name)
        geoms=rasterio.warp.transform_geom('EPSG:4326','EPSG:4326',[shp['geometry']])
        mask=self.create_mask_from_shp_cci(geoms)
        for varname in list(ds.keys()):
            if len(ds[varname].shape)==3:
                if ds[varname][0,:,:].shape==mask.shape:
                    ds[varname]=ds[varname]*mask
        if crop==True:
            bounds=rasterio.features.bounds(geoms[0])
            ds=ds.sel(lon=slice(bounds[0],bounds[2]),lat=slice(bounds[3],bounds[1]))  
        return ds
    
    def what_tiles_for_poly(self,poly):
        my_func_name='what_tiles_for_poly'
        print(datetime.now().strftime("%H:%M:%S")+' '+my_func_name)
        geoms=rasterio.warp.transform_geom('EPSG:4326','EPSG:3035',[poly['geometry']])
        tiles=['t01','t02','t03','t04']
        nec_tiles=[]
        for tile in tiles:
            output=self.create_mask_from_shp_tiles(geoms,tile)
            if np.nanmax(output)==1:
                nec_tiles.append(tile)
        return nec_tiles
    
    def get_vars_by_coords(self,xds,variables,lat,lon):
        xds_pnt=xds.sel(x=lon,y=lat, method='nearest') 
        stats=[]
        for var in variables:
            stats.append(xds_pnt[var].values)
        df=pd.DataFrame(data=np.transpose(stats),columns=variables)
        #df=pd.DataFrame(data=np.transpose(stats),columns=variables)
        return df
    
    def get_vars_by_coords_buffer(self,xds,variables,lat,lon,buffer):
        my_func_name='get_vars_by_coords_buffer'
        print(datetime.now().strftime("%H:%M:%S")+' '+my_func_name)
        xds_pnt=xds.where((xds.x>(lon-buffer*1000))&(xds.x<(lon+buffer*1000))&
                  (xds.y>(lat-buffer*1000))&(xds.y<(lat+buffer*1000)),drop=True)
        stats=[]
        varnames=[]
        for var in variables:
            values=xds_pnt[var].values
            var_mean=np.nanmean(values,axis=(1,2))
            var_std=np.nanstd(values,axis=(1,2))
            var_val=(np.count_nonzero(~np.isnan(values),axis=(1,2))/float(values.shape[1]*values.shape[2]))*100
            stats.append(var_mean)
            stats.append(var_std)
            stats.append(var_val)
            varnames.append(var)
            varnames.append(var+'_std')
            varnames.append(var+'_percval')
            #stats.append(xds_pnt[var].values)
        df=pd.DataFrame(data=np.transpose(stats),columns=varnames)
        #df=pd.DataFrame(data=np.transpose(stats),columns=variables)
        return df
        
    def coord_in_tile(self,lat,lon,tl):
        tile_x={
            "t01":(900000,4150000),
            "t02":(4150000,7400000),
            "t03":(900000,4150000),
            "t04":(4150000,7400000)
            }
        tile_y={
            "t01":(3200000,5500000),
            "t02":(3200000,5500000),
            "t03":(900000,3200000),
            "t04":(900000,3200000)
            }
        if tile_x[tl][0] <= lon <= tile_x[tl][1] and tile_y[tl][0] <= lat <= tile_y[tl][1]:
            return True
        else:
            return False
        
class lst_od_correction:
    def daily_anomalies(self,data,var):
        aggr=data.groupby(['day',"month"])[var].agg(lambda x: np.nanmean(x)).rename(var+'_agg')
        data=pd.merge(data,aggr,how='inner', on=['day',"month"])
        data[var+'_anom']=data[var]-data[var+'_agg']
        return data
    
    def second_dgr_fit(self,dates,data):
        dates=mdates.date2num(dates)
        fit_param=np.polyfit(dates,data,2)
        fit=fit_param[0]*(dates**2)+fit_param[1]*dates+fit_param[2]
        return fit       
    
    def julien_c0(self,X, a, b, c):
        dates,sza_fit=X
        dates=mdates.date2num(dates)
        return a+b*dates+c*sza_fit
    
    def julien_c1(self,sza_fit, a,b):
        return a+b*sza_fit
    
    def dtm_corr_lst(self,data,t1,dtm_params):
        # data must contain lst, view time and doy
        data=data[np.isfinite(data['view_time'])&np.isfinite(data['lst'])]
        data=pd.merge(data,dtm_params,how='inner', on=['doy'])
        t2=np.array(data['view_time'])
        lst=np.array(data['lst'])
        Ta=np.array(data['Ta'])
        om=np.array(data['omega'])
        T0=np.array(data['T0'])
        tm=np.array(data['tm'])
        data['lst_dtc']=lst+Ta*np.cos((math.pi/om)*(t1-tm))-Ta*np.cos((math.pi/om)*(t2-tm))
        return data
    
    def get_vars_by_coords_seviri(self,xds,variables,lat,lon):    
        xds_pnt=xds.sel(lon=lon,lat=lat, method='nearest')  
        stats=[]
        for var in variables:
            stats.append(xds_pnt[var].values)
        variables.append('time')
        stats.append(xds_pnt['time'].values)    
        df=pd.DataFrame(data=np.transpose(stats),columns=variables)
        return df
    
    def get_T0_Ta_tm(self,df,doy,ts):
        df_sub=df[df['doys']==doy]
        if len(df_sub[np.isfinite(df_sub['lst'])]) >0:      
            tm=float(df_sub[df_sub['lst']==max(df_sub['lst'])]['times'])
            Ta=max(df_sub['lst'])-min(df_sub['lst'])
            T0=df_sub[(abs(df_sub['times']-ts))==(min(abs(df_sub['times']-ts)))]['lst']
            if len(T0)>1:
                T0=T0.iloc[0]
            elif len(T0)<1:
                T0=np.nan
            return float(T0), Ta, tm
        else:
            return np.nan,np.nan,np.nan
        
    def get_sunrise_daytime_length(self,sunrise,sunset):
        sunrise=sunrise.hour+sunrise.minute/60
        sunset=sunset.hour+sunset.minute/60
        omega=sunset-sunrise
        return sunrise, omega
        
    def dtm_from_seviri(self,file,lat,lon):
        ds=xr.open_dataset(file,decode_timedelta=False)
        vals=self.get_vars_by_coords_seviri(ds,['mast','yast','theta'],lat,lon)
        # Compute continues lst with parameters
        datetimes=[]
        for date in self.datespan(datetime(2000,1,1), datetime(2000,12,31),delta=timedelta(hours=0.5)):
            datetimes.append(date)
        doys=[]
        times=[]
        for doy in range(-79,286): # days from march equinox
            for time in np.arange(0,24,0.5):
                doys.append(doy)
                times.append(time)
        df=pd.DataFrame(list(zip(datetimes,doys,times)),columns =['datetimes','doys','times'])
        df['lst']=df.apply(lambda row: self.daily_annual_lst_from_seviri(vals,row['doys'],row['times']),axis=1)
        df['doys']=df['doys']+80 # Change doy to normal again
        s = sun(lat=lat,long=lon)
        dtm_params=pd.DataFrame()
        for doy in df['doys'].unique(): 
            # Compute sunrise and daytime length
            #pdb.set_trace()
            dt=datetime(2000, 1, 1) + timedelta(int(doy) - 1)
            sunrise=s.sunrise(when=dt)
            sunset=s.sunset(when=dt)
            sunrise,omega=self.get_sunrise_daytime_length(sunrise,sunset)
            #pdb.set_trace()
            # Get other parameters
            T0, Ta, tm=self.get_T0_Ta_tm(df, doy,sunrise) 
            dat=pd.DataFrame([[doy,T0,Ta,tm,omega]],columns=['doy','T0','Ta','tm','omega'])
            dtm_params=dtm_params.append(dat)
        return dtm_params    
        
    
    def datespan(self, startDate, endDate, delta):
        currentDate = startDate
        while currentDate < endDate:
            yield currentDate
            currentDate += delta
        
    def daily_annual_lst_from_seviri(self,acp,d,time):
        mast=float(acp[acp['time']==time].mast)
        theta=float(acp[acp['time']==time].theta)
        yast=float(acp[acp['time']==time].yast) 
        lst=mast+yast*math.sin((2*math.pi/365)*(d+theta))
        return lst
    
    
    
   
        
        
    

        
        
        
        
        
        
        
        
        
