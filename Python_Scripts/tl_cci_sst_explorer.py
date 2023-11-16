# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 11:28:38 2022

@author: rein_pp
"""
import tl3_analysis_toolbox
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
import rioxarray
from datetime import timedelta
import os
import fiona
import pymannkendall as mk
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse
from datetime import date
import warnings


#import seaborn as sns
def daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + timedelta(n)
        
        
def interpolate_3darr(arr,axis,method,limit):
    result = np.zeros_like(arr, dtype=np.float32)
    for i in range(arr.shape[axis]):
        if axis==0:
            line_stack = pd.DataFrame(data=arr[i,:,:], dtype=np.float32)
        if axis==1:
            line_stack = pd.DataFrame(data=arr[:,i,:], dtype=np.float32)
        if axis==2:
            line_stack = pd.DataFrame(data=arr[:,:,i], dtype=np.float32)
        line_stack.interpolate(method=method, axis=0, inplace=True, limit=limit)
        pdb.set_trace()
        if axis==0:
            result[i, :, :] = line_stack.values.astype(np.float32)
        if axis==1:
            result[:, i, :] = line_stack.values.astype(np.float32)
        if axis==2:
            result[:, :, i] = line_stack.values.astype(np.float32)
    return result

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
    
        
def interpolate_nan(arr, method="linear", limit=3): # Von Christina
    """return array interpolated along time-axis to fill missing values"""
#    result = np.zeros_like(arr, dtype=np.int16)
    result = np.zeros_like(arr, dtype=np.float32)

#    for i in range(arr.shape[1]):
#        # slice along y axis, interpolate with pandas wrapper to interp1d
    for i in range(arr.shape[0]):
        # slice along y axis, interpolate with pandas wrapper to interp1d ---> # a-axis???
#        line_stack = pd.DataFrame(data=arr[:,i,:], dtype=np.float32)
        line_stack = pd.DataFrame(data=arr[i,:,:], dtype=np.float32)
#        line_stack.replace(to_replace=-37268, value=np.NaN, inplace=True)
        line_stack.interpolate(method=method, axis=0, inplace=True, limit=limit)
        line_stack.replace(to_replace=np.NaN, value=-37268, inplace=True)
#        result[:, i, :] = line_stack.values.astype(np.int16)
        result[i, :, :] = line_stack.values.astype(np.float32)
    return result

def interp_along_axis(y, x, newx, axis, inverse=False, method='linear'):
    """ Interpolate vertical profiles, e.g. of atmospheric variables
    using vectorized numpy operations

    This function assumes that the x-xoordinate increases monotonically

    ps:
    * Updated to work with irregularly spaced x-coordinate.
    * Updated to work with irregularly spaced newx-coordinate
    * Updated to easily inverse the direction of the x-coordinate
    * Updated to fill with nans outside extrapolation range
    * Updated to include a linear interpolation method as well
        (it was initially written for a cubic function)

    Peter Kalverla
    March 2018

    --------------------
    More info:
    Algorithm from: http://www.paulinternet.nl/?page=bicubic
    It approximates y = f(x) = ax^3 + bx^2 + cx + d
    where y may be an ndarray input vector
    Returns f(newx)

    The algorithm uses the derivative f'(x) = 3ax^2 + 2bx + c
    and uses the fact that:
    f(0) = d
    f(1) = a + b + c + d
    f'(0) = c
    f'(1) = 3a + 2b + c

    Rewriting this yields expressions for a, b, c, d:
    a = 2f(0) - 2f(1) + f'(0) + f'(1)
    b = -3f(0) + 3f(1) - 2f'(0) - f'(1)
    c = f'(0)
    d = f(0)

    These can be evaluated at two neighbouring points in x and
    as such constitute the piecewise cubic interpolator.
    """

    # View of x and y with axis as first dimension
    if inverse:
        _x = np.moveaxis(x, axis, 0)[::-1, ...]
        _y = np.moveaxis(y, axis, 0)[::-1, ...]
        _newx = np.moveaxis(newx, axis, 0)[::-1, ...]
    else:
        _y = np.moveaxis(y, axis, 0)
        _x = np.moveaxis(x, axis, 0)
        _newx = np.moveaxis(newx, axis, 0)

    # Sanity checks
    if np.any(_newx[0] < _x[0]) or np.any(_newx[-1] > _x[-1]):
        # raise ValueError('This function cannot extrapolate')
        warnings.warn("Some values are outside the interpolation range. "
                      "These will be filled with NaN")
    if np.any(np.diff(_x, axis=0) < 0):
        raise ValueError('x should increase monotonically')
    if np.any(np.diff(_newx, axis=0) < 0):
        raise ValueError('newx should increase monotonically')

    # Cubic interpolation needs the gradient of y in addition to its values
    if method == 'cubic':
        # For now, simply use a numpy function to get the derivatives
        # This produces the largest memory overhead of the function and
        # could alternatively be done in passing.
        ydx = np.gradient(_y, axis=0, edge_order=2)

    # This will later be concatenated with a dynamic '0th' index
    ind = [i for i in np.indices(_y.shape[1:])]

    # Allocate the output array
    original_dims = _y.shape
    newdims = list(original_dims)
    newdims[0] = len(_newx)
    newy = np.zeros(newdims)

    # set initial bounds
    i_lower = np.zeros(_x.shape[1:], dtype=int)
    i_upper = np.ones(_x.shape[1:], dtype=int)
    x_lower = _x[0, ...]
    x_upper = _x[1, ...]

    for i, xi in enumerate(_newx):
        # Start at the 'bottom' of the array and work upwards
        # This only works if x and newx increase monotonically

        # Update bounds where necessary and possible
        needs_update = (xi > x_upper) & (i_upper+1<len(_x))
        # print x_upper.max(), np.any(needs_update)
        while np.any(needs_update):
            i_lower = np.where(needs_update, i_lower+1, i_lower)
            i_upper = i_lower + 1
            x_lower = _x[[i_lower]+ind]
            x_upper = _x[[i_upper]+ind]

            # Check again
            needs_update = (xi > x_upper) & (i_upper+1<len(_x))

        # Express the position of xi relative to its neighbours
        xj = (xi-x_lower)/(x_upper - x_lower)

        # Determine where there is a valid interpolation range
        within_bounds = (_x[0, ...] < xi) & (xi < _x[-1, ...])

        if method == 'linear':
            f0, f1 = _y[[i_lower]+ind], _y[[i_upper]+ind]
            a = f1 - f0
            b = f0

            newy[i, ...] = np.where(within_bounds, a*xj+b, np.nan)

        elif method=='cubic':
            f0, f1 = _y[[i_lower]+ind], _y[[i_upper]+ind]
            df0, df1 = ydx[[i_lower]+ind], ydx[[i_upper]+ind]

            a = 2*f0 - 2*f1 + df0 + df1
            b = -3*f0 + 3*f1 - 2*df0 - df1
            c = df0
            d = f0

            newy[i, ...] = np.where(within_bounds, a*xj**3 + b*xj**2 + c*xj + d, np.nan)

        else:
            raise ValueError("invalid interpolation method"
                             "(choose 'linear' or 'cubic')")

    if inverse:
        newy = newy[::-1, ...]

    return np.moveaxis(newy, 0, axis)

def interpolate_nan2(arr, method="linear", limit=3): # Von Christina
    """return array interpolated along time-axis to fill missing values"""
#    result = np.zeros_like(arr, dtype=np.int16)
    result = np.zeros_like(arr, dtype=np.float32)

#    for i in range(arr.shape[1]):
#        # slice along y axis, interpolate with pandas wrapper to interp1d
    for i in range(arr.shape[1]):
        # slice along y axis, interpolate with pandas wrapper to interp1d ---> # a-axis???
#        line_stack = pd.DataFrame(data=arr[:,i,:], dtype=np.float32)
        line_stack = pd.DataFrame(data=arr[:,1,:], dtype=np.float32)
        pdb.set_trace()
#        line_stack.replace(to_replace=-37268, value=np.NaN, inplace=True)
        line_stack.interpolate(method=method, axis=0, inplace=True, limit=limit)
        #line_stack.replace(to_replace=np.NaN, value=-37268, inplace=True)
#        result[:, i, :] = line_stack.values.astype(np.int16)
        result[i, :, :] = line_stack.values.astype(np.float32)
    return result

def nan_along_time_to_zero(arr):
    my_func_name='nan_along_time_to_zero'
    print(datetime.now().strftime("%H:%M:%S")+' '+my_func_name)
    # Set all pixel with complete NaN timeseries to 0
    x=np.isnan(np.nanmean(arr,0))
    x_3d=np.broadcast_to(x,arr.shape)
    arr[x_3d]=0
    return arr
       
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


def write_nc(xds, fn):
    '''
    Write xarray as netcdf
    '''
    comp = dict(zlib=True, complevel=1)
    encoding = {var: comp for var in xds.data_vars}
    xds.to_netcdf(fn, encoding=encoding)
    #xds.to_netcdf(fn) 
    


def plt_image(arr,vmin,vmax,title,fig,pltnr):
    fontsize=60
    cmap = mpl.cm.Spectral_r
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    #fig=plt.figure()
    #fig.set_figwidth(24)
    #fig.set_figheight(12)
    ax=fig.add_subplot(pltnr)
    im=ax.imshow(arr,cmap=cmap,norm=norm)
    ax.set_axis_off()
    ax.set_title(title,fontsize=fontsize)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar=fig.colorbar(im,cax=cax)
    cbar.ax.tick_params(axis='both',labelsize=fontsize)
    return fig
    
def mk_maps(xds,product,path_plot):
    var_dict={'trend':[-1,1],
              'p':[0,0.4],
              'Tau':[0,0.4],
              'slope':[-0.1,0.1],
              'obs_count':[0,35]
        }
    mask=np.array(xds['slope'])==0 # Land pixels
    fig=plt.figure()
    fig.set_figwidth(72)
    fig.set_figheight(24)
    pltnrs=[231,232,233,234,235,236]
    i=0
    #for var in var_dict.keys(): #var_dict.keys()
    for var in ['trend','p','Tau','slope']:
        vmin=var_dict[var][0]
        vmax=var_dict[var][1] 
        arr=np.array(xds[var])
        arr=arr.astype(float)
        arr[mask]=np.nan
        pltnr=pltnrs[i]
        fig=plt_image(arr,vmin,vmax,var,fig,pltnr)
        i=i+1
    
    
    # Plot slope only for significant areas
    mask=np.array(xds['trend'])==0
    vmin=var_dict['slope'][0]
    vmax=var_dict['slope'][1]
    arr=np.array(xds['slope'])
    arr[mask]=np.nan
    pltnr=pltnrs[i]
    fig=plt_image(arr,vmin,vmax,'slope_masked',fig,pltnr)     
    outfile=path_plot+product+'.png'
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close()
    
    

        
def confidence_plot_sst_dif(data,title,ylabel,outfile):
    fontsize=15
    fig=plt.figure()
    fig.set_figwidth(20)
    fig.set_figheight(4)
    ax= fig.add_subplot(111)
    # Median plot
    p=sns.lineplot(data=data, x='date',y='q2',ax=ax,marker='o',color='blue',
                   legend='brief',label='Median',ms=1)
    # confidence intervals
    f1=ax.fill_between(x=data['date'],
                    y1=data['q1'],
                    y2=data['q3'], alpha=0.2,color='blue')    
    # legends
    l1=ax.legend(loc='upper left', bbox_to_anchor=(0.8, 0.5, 0.5, 0.5),frameon=False)
    l2=ax.legend(handles=[f1], labels=['Quantile 0.05-0.95'],loc='upper left', 
                                       bbox_to_anchor=(0.8, 0.42, 0.5, 0.5),frameon=False)
    # axes
    start=datetime(1982,1,1)
    stop=datetime(2017,1,1)
    ax.set_ylim(-10,10)
    ax.set_xlim(start,stop)
    ax.add_artist(l1)
    ax.add_artist(l2)
    ax.set_title(title,fontsize=fontsize)
    ax.set_xlabel("Time",fontsize=fontsize)
    ax.set_ylabel(ylabel,fontsize=fontsize)
    # zero line
    ax.plot([start,stop],[0,0],color='black')     
    fig.savefig(outfile)
    
def confidence_plot_sst(data_tl,data_cci,title,ylabel,outfile):
    fontsize=15
    fig=plt.figure()
    fig.set_figwidth(20)
    fig.set_figheight(4)
    ax= fig.add_subplot(111)
    # Median plot
    p=sns.lineplot(data=data_cci, x='date',y='q2',ax=ax,marker='o',color='red',
                   legend='brief',label='CCI Median',ms=4)
    # Median plot
    p=sns.lineplot(data=data_tl, x='date',y='q2',ax=ax,marker='o',color='blue',
                   legend='brief',label='TIMELINE Median',ms=4)  
    # confidence intervals
    f1=ax.fill_between(x=data_cci['date'],
                    y1=data_cci['q1'],
                    y2=data_cci['q3'], alpha=0.4,color='red')    
    f2=ax.fill_between(x=data_tl['date'],
                    y1=data_tl['q1'],
                    y2=data_tl['q3'], alpha=0.4,color='blue')    
    # legends
    l1=ax.legend(loc='upper left', bbox_to_anchor=(0.83, 0.5, 0.5, 0.5),frameon=False)
    l2=ax.legend(handles=[f1,f2], labels=['CCI Quantile 0.05-0.95','TIMELINE Quantile 0.05-0.95'],loc='upper left', 
                                       bbox_to_anchor=(0.83, 0.35, 0.5, 0.5),frameon=False)
    # axes
    start=datetime(1982,1,1)
    stop=datetime(2017,1,1)
    ax.set_ylim(270,305)
    ax.set_xlim(start,stop)
    ax.add_artist(l1)
    ax.add_artist(l2)
    ax.set_title(title,fontsize=fontsize)
    ax.set_xlabel("Time",fontsize=fontsize)
    ax.set_ylabel(ylabel,fontsize=fontsize)
    # zero line
    ax.plot([start,stop],[0,0],color='black')     
    fig.savefig(outfile)
    
def confidence_plot_sst_dif_filtered(data_list,title,ylabel,outfile):
    fontsize=15
    fig=plt.figure()
    fig.set_figwidth(20)
    fig.set_figheight(4)
    ax= fig.add_subplot(111)
    # Median plot
    p=sns.lineplot(data=data_list[0], x='date',y='q2',ax=ax,marker='o',color='blue',
                   legend='brief',label='CCI Median',ms=4)
    # Median plot
    p=sns.lineplot(data=data_list[1], x='date',y='q2',ax=ax,marker='o',color='green',
                   legend='brief',label='TIMELINE Median',ms=4)  
    # Median plot
    p=sns.lineplot(data=data_list[2], x='date',y='q2',ax=ax,marker='o',color='orange',
                   legend='brief',label='CCI Median',ms=4)
    # Median plot
    p=sns.lineplot(data=data_list[3], x='date',y='q2',ax=ax,marker='o',color='red',
                   legend='brief',label='TIMELINE Median',ms=4)  
    
    # confidence intervals
    f1=ax.fill_between(x=data_list[0]['date'],
                    y1=data_list[0]['q1'],
                    y2=data_list[0]['q3'], alpha=0.4,color='blue')    
    f2=ax.fill_between(x=data_list[1]['date'],
                    y1=data_list[1]['q1'],
                    y2=data_list[1]['q3'], alpha=0.4,color='green')    
    f3=ax.fill_between(x=data_list[2]['date'],
                    y1=data_list[2]['q1'],
                    y2=data_list[2]['q3'], alpha=0.4,color='orange')    
    if len(data_list[3]>0):
        f4=ax.fill_between(x=data_list[3]['date'],
                        y1=data_list[3]['q1'],
                        y2=data_list[3]['q3'], alpha=0.4,color='red')    
    # legends
    l1=ax.legend(loc='upper left', bbox_to_anchor=(0.83, 0.5, 0.5, 0.5),frameon=False)
    l2=ax.legend(handles=[f1,f2], labels=['CCI Quantile 0.05-0.95','TIMELINE Quantile 0.05-0.95'],loc='upper left', 
                                       bbox_to_anchor=(0.83, 0.35, 0.5, 0.5),frameon=False)
    # axes
    start=datetime(1982,1,1)
    stop=datetime(2017,1,1)
    ax.set_ylim(-5,5)
    ax.set_xlim(start,stop)
    ax.add_artist(l1)
    ax.add_artist(l2)
    ax.set_title(title,fontsize=fontsize)
    ax.set_xlabel("Time",fontsize=fontsize)
    ax.set_ylabel(ylabel,fontsize=fontsize)
    # zero line
    ax.plot([start,stop],[0,0],color='black')     
    fig.savefig(outfile)
    
def calc_quantiles(path,files,var):
    q1=[]
    q2=[]
    q3=[]
    date=[]
    for file in files:
        dat=pd.read_csv(path+file)
        dat=dat[np.isfinite(dat.sst_dif)]
        if len(dat)>0:
            q1.append(np.quantile(dat[var],0.05))
            q2.append(np.quantile(dat[var],0.5))
            q3.append(np.quantile(dat[var],0.95))
            if file[0:6]=='Baltic':
                date.append(datetime(int(file[11:15]),int(file[16:18]),int(file[19:21])))
            if file[0:5]=='North':
                date.append(datetime(int(file[10:14]),int(file[15:17]),int(file[18:20])))
    data=pd.DataFrame({'date':date,'q1':q1,'q2':q2,'q3':q3})
    return data
    
def add_sst_qual(data):
    data['sst_qual']=np.nan
    data.loc[(data['sst_cci']>275)&(data['sst_cci']<280),'sst_qual']=0
    data.loc[(data['sst_cci']>280)&(data['sst_cci']<285),'sst_qual']=1
    data.loc[(data['sst_cci']>285)&(data['sst_cci']<290),'sst_qual']=2
    data.loc[(data['sst_cci']>290)&(data['sst_cci']<295),'sst_qual']=3
    return data
                                                   

def calc_quantiles_filtered(path,files,var,var_filt):
    data_list=[]
    for i in range(4):
        q1=[]
        q2=[]
        q3=[]
        date=[]
        for file in files:
            dat=pd.read_csv(path+file)
            dat=dat[np.isfinite(dat.sst_dif)]
            dat=add_sst_qual(dat)
            dat=dat[dat[var_filt]==i]
            if len(dat)>0:
                q1.append(np.quantile(dat[var],0.05))
                q2.append(np.quantile(dat[var],0.5))
                q3.append(np.quantile(dat[var],0.95))
                if file[0:6]=='Baltic':
                    date.append(datetime(int(file[11:15]),int(file[16:18]),int(file[19:21])))
                if file[0:5]=='North':
                    date.append(datetime(int(file[10:14]),int(file[15:17]),int(file[18:20])))
        data=pd.DataFrame({'date':date,'q1':q1,'q2':q2,'q3':q3})
        data_list.append(data)
    return data_list

def write_masked_slope(xds,p,outfile):
    xds=xds.where(xds.p<p)
    write_nc(xds, outfile)
    
    
def analyze_mk():
    path_mk='E:/SST_Analysis/MannKendall/'
    path_plot='E:/SST_Analysis/Plots/'
    path_lax='E:/Conferences_Presentations/LAX_Seminar_202210/Data/'
    suffixes=['max_no_filter','max_dif_filter','Stdev_test']
    suffixes=['_test2023']
    for suffix in suffixes:
        for month in [7,8]: #range(1,13)
            m_str=str(month).zfill(2)           
            '''
            # Baltic Sea
            # CCI
            product='Baltic Sea_'+m_str+'_sst_cci_MannKendall'
            file_mk=path_mk+product+'.nc'
            xds_mk_cci=xr.open_dataset(file_mk)
            mk_maps(xds_mk_cci,product,path_plot)
            outfile=path_lax+product+'_slope.nc'
            write_masked_slope(xds_mk_cci, 0.1, outfile)
            
            # Timeline
            product='Baltic Sea_'+m_str+'_sst_max'+suffix+'_MannKendall'
            file_mk=path_mk+product+'.nc'
            xds_mk_tl=xr.open_dataset(file_mk)
            mk_maps(xds_mk_tl,product,path_plot)
            outfile=path_lax+product+'_slope.nc'
            write_masked_slope(xds_mk_tl, 0.1, outfile)
           
            '''
            # North Sea
            # CCI
            product='North Sea_'+m_str+'_sst_cci_MannKendall'
            file_mk=path_mk+product+'.nc'
            xds_mk_cci=xr.open_dataset(file_mk)
            mk_maps(xds_mk_cci,product,path_plot)
            outfile=path_lax+product+'_slope.nc'
            write_masked_slope(xds_mk_cci, 0.1, outfile)
           
            # Timeline
            product='North Sea_'+m_str+'_sst_max'+suffix+'_MannKendall'
            file_mk=path_mk+product+'.nc'
            xds_mk_tl=xr.open_dataset(file_mk)
            mk_maps(xds_mk_tl,product,path_plot)    
            outfile=path_lax+product+'_slope.nc'
            write_masked_slope(xds_mk_tl, 0.1, outfile)

def analyze_val():
    path_val='E:/SST_Analysis/Validationfiles/'
    path_plot='E:/SST_Analysis/Plots/'
    path_res='E:/SST_Analysis/Results/'
    
    # Baltic Sea
    files=[file for file in os.listdir(path_val) if file[0:6]=='Baltic']
    
    res_file=path_res+'baltic_dif_quantiles.csv'
    var='sst_dif'
    #data=calc_quantiles(path_val,files,var)
    #data.to_csv(res_file)
    
    data=pd.read_csv(res_file)
    data['date']=data.date.apply(lambda x:datetime(int(x[0:4]),int(x[5:7]),int(x[8:10])))
    outfile=path_plot+'baltic_confidence_dif.png'
    confidence_plot_sst_dif(data, 'Baltic Sea SST difference', 'TIMELINE - CCI SST[K]', outfile)
    
    res_file=path_res+'baltic_cci_quantiles.csv'
    var='sst_cci'
    #data=calc_quantiles(path_val,files,var)
    #data.to_csv(res_file)
    
    res_file=path_res+'baltic_tl_quantiles.csv'
    var='sst_max'
    #data=calc_quantiles(path_val,files,var)
    #data.to_csv(res_file)
    
    res_file_cci=path_res+'baltic_cci_quantiles.csv'
    data_cci=pd.read_csv(res_file_cci)
    data_cci['date']=data_cci.date.apply(lambda x:datetime(int(x[0:4]),int(x[5:7]),int(x[8:10])))
    res_file_tl=path_res+'baltic_tl_quantiles.csv'
    data_tl=pd.read_csv(res_file_tl)
    data_tl['date']=data_tl.date.apply(lambda x:datetime(int(x[0:4]),int(x[5:7]),int(x[8:10])))  
    
    outfile=path_plot+'baltic_confidence_sst.png'
    confidence_plot_sst(data_tl, data_cci,'Baltic Sea TIMELINE and CCI SST', 'SST[K]', outfile)
   
    # Filter sst by variables
    var='sst_dif'
    var_filts=['tcwv_qual','sz_qual','sst_qual']
    for var_filt in var_filts:
        data_list=calc_quantiles_filtered(path_val,files,var,var_filt)
        res_file=path_res+'baltic_dif_quantiles_'+var_filt+'.csv'
        for i in range(4):
            data=data_list[i]
            res_file=path_res+'baltic_dif_quantiles_'+var_filt+'_'+str(i)+'.csv'
            data.to_csv(res_file)
    
    var_filts=['tcwv_qual','sz_qual','sst_qual']
    for var_filt in var_filts:
        data_list=[]
        for i in range(4):
            res_file=path_res+'baltic_dif_quantiles_'+var_filt+'_'+str(i)+'.csv'
            data=pd.read_csv(res_file)
            if len(data)>0:
                data['date']=data.date.apply(lambda x:datetime(int(x[0:4]),int(x[5:7]),int(x[8:10])))
                outfile=path_plot+'baltic_confidence_dif'+var_filt+'_'+str(i)+'.png'
                confidence_plot_sst_dif(data, 'Baltic Sea SST difference '+var_filt+'_'+str(i), 'TIMELINE - CCI SST[K]', outfile)
    
    # North Sea
    files=[file for file in os.listdir(path_val) if file[0:5]=='North']
    
    res_file=path_res+'north_dif_quantiles.csv'
    var='sst_dif'
    #data=calc_quantiles(path_val,files,var)
    #data.to_csv(res_file)
    
    data=pd.read_csv(res_file)
    data['date']=data.date.apply(lambda x:datetime(int(x[0:4]),int(x[5:7]),int(x[8:10])))
    outfile=path_plot+'north_confidence_dif.png'
    confidence_plot_sst_dif(data, 'North Sea SST difference', 'TIMELINE - CCI SST[K]', outfile)
    
    res_file=path_res+'north_cci_quantiles.csv'
    var='sst_cci'
    #data=calc_quantiles(path_val,files,var)
    #data.to_csv(res_file)
    
    res_file=path_res+'north_tl_quantiles.csv'
    var='sst_max'
    #data=calc_quantiles(path_val,files,var)
    #data.to_csv(res_file)
    
    res_file_cci=path_res+'north_cci_quantiles.csv'
    data_cci=pd.read_csv(res_file_cci)
    data_cci['date']=data_cci.date.apply(lambda x:datetime(int(x[0:4]),int(x[5:7]),int(x[8:10])))
    res_file_tl=path_res+'north_tl_quantiles.csv'
    data_tl=pd.read_csv(res_file_tl)
    data_tl['date']=data_tl.date.apply(lambda x:datetime(int(x[0:4]),int(x[5:7]),int(x[8:10])))  
    
    outfile=path_plot+'north_confidence_sst.png'
    confidence_plot_sst(data_tl, data_cci,'North Sea TIMELINE and CCI SST', 'SST[K]', outfile)
    
def filter_stdev(xds,tem_var): 
    print(datetime.now().strftime("%H:%M:%S")+'Filter Temperatures with STD')
    # Filter out unrealistic Temperature values 
    var=np.array(xds[tem_var])
    var[var<260]=np.NAN
    var[var>340]=np.NAN
    xds[tem_var]=(['t','y','x'],var)
        
    # Filter all sst 
    var_std=np.array(xds[tem_var].std(dim='t'))
    var_mean=np.array(xds[tem_var].mean(dim='t'))
    obs_count=np.array(xds[tem_var].count(dim='t'))
    xds.coords['var_std']=(["y","x"],var_std)
    xds.coords['var_mean']=(["y","x"],var_mean)
    xds.coords['obs_count']=(["y","x"],obs_count)
    xds['var_difmean']=abs(xds[tem_var]-xds.coords['var_mean'])# Calculate difference from mean
    xds['var_difmean']=xds.var_difmean.where(xds.coords['obs_count']>4,0) #Too less values for STD filtering
    xds[tem_var]=xds[tem_var].where(xds.var_difmean < xds.coords['var_std']*1.5) # Remove outliers (Dif > 1.5*STD)
    xds=xds.drop(['var_difmean','var_std','obs_count','var_mean']) # Drop temporal variables  
    return xds
    
def initargs():
    '''
    Argument definition, also handles user input
    '''
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--filter', help='Dif or Stdev',required=False)
    parser.add_argument('--suffix', help='Dif or Stdev',required=False)
    return parser.parse_args()



def main_kendall():
    path_out='/nfs/IGARSS_2022/MannKendall/'
    shp='/nfs/IGARSS_2022/Shps/coast_dk.shp'
    
    aux=tl3_analysis_toolbox.aux()
    tl_crop=tl3_analysis_toolbox.crop()
    prep=tl3_analysis_toolbox.l3_lst_sst_prep()
    reproj=tl3_analysis_toolbox.reproj()
    
    month=4
    
    
    
    
    shapefile=fiona.open(shp, "r")
    for poly in [shapefile[0]]:
        monthly_xds=xr.Dataset(coords={'x':[],'y':[],'t':[]})
        tiles = tl_crop.what_tiles_for_poly(poly)
        for year in range(1990,2018,1):
            print(year)
            xds_list=[]
            for tile in tiles:
                fs_month=[]
                for day in [1,10,20]:    
                    date_f=str(year)+str(month).zfill(2)+str(day).zfill(2)
                    fn=aux.get_l3_file_from_ida(date_f,'Decades','SST','v01.01',tile)
                    fs_month.append(fn)
                
                xds = xr.open_mfdataset(fs_month, combine='nested', concat_dim=['t'], 
                                         chunks={'x': 2300, 'y': 3250},
                                         preprocess=prep.add_date)
                xds=xds.drop(['lambert_azimuthal_equal_area'])
                crop=True
                xds=tl_crop.crop_ds_with_shp_tiles(xds,poly,crop,tile)
                xds_list.append(xds)
                
            chunksize=[1000,1000,37]
            var=['sst_max','view_time_max','view_day']
            xds_merged= reproj.mosaic_vars(xds_list,[var],chunksize)
            # Compute view day max
            view_day_max=(np.array(xds_merged.view_day)/1000000).astype(int)
            view_day_max=view_day_max.astype(float)
            view_day_max[view_day_max<0]=np.nan
            xds_merged['view_day_max']=(['t','y','x'],view_day_max)
            
            yearly_xds=xr.Dataset(coords={'x':xds_merged.coords['x'],'y':xds_merged.coords['y'],'t':year})         
            yearly_xds['sst_max']=(['y','x'],xds_merged['sst_max'].max(dim='t').data) 
            doy_max=np.array(xds_merged['view_day_max'].where(xds_merged['sst_max']==xds_merged['sst_max'].max(dim='t')).max(dim='t'))
            yearly_xds['doy_max']=(['y','x'],doy_max)
            monthly_xds=xr.concat([monthly_xds,yearly_xds],dim='t')
        poly_id=poly['properties']['id']
        
        outfile=path_out+str(poly_id)+'_'+str(month).zfill(2)+'_monthly_sst.nc'
        write_nc(monthly_xds, outfile)
        spco=['y','x']
        xds_res=ApplyMannKendall(monthly_xds,'sst_max',spco) 
        outfile=path_out+str(poly_id)+'_'+str(month).zfill(2)+'_mk.nc'
        write_nc(xds_res, outfile) 
                 
def anomaly_trends():
    path_out='/nfs/IGARSS_2022/MannKendall/'
    shp='/nfs/IGARSS_2022/Shps/coast_dk.shp'
    start_date=date(2000,8,1)
    end_date=date(2000,8,4)
   
    aux=tl3_analysis_toolbox.aux()
    tl_crop=tl3_analysis_toolbox.crop()
    prep=tl3_analysis_toolbox.l3_lst_sst_prep()
    reproj=tl3_analysis_toolbox.reproj()
    
    shapefile=fiona.open(shp, "r")
    poly=shapefile[0]
    poly_id=poly['properties']['id']
    tiles = tl_crop.what_tiles_for_poly(poly)
    xds_mns=xr.Dataset(coords={'x':[],'y':[],'doi':[]})
    for dt in daterange(start_date, end_date):
        print(dt)
        doi=dt.timetuple().tm_yday
        xds_list=[]
        for tile in tiles:
            fns=[]
            for year in range(1990,2018,1):
                date_f=str(year)+str(dt.month).zfill(2)+str(dt.day).zfill(2)
                fn=aux.get_l3_file_from_ida(date_f,'Daily','SST','v01.01',tile)
                if fn !=None:
                    fns.append(fn)
            xds = xr.open_mfdataset(fns, combine='nested', concat_dim=['t'], 
                                     chunks={'x': 2300, 'y': 3250},
                                     preprocess=prep.add_date)
            xds=xds.drop(['lambert_azimuthal_equal_area'])
            crop=True
            xds=tl_crop.crop_ds_with_shp_tiles(xds,poly,crop,tile)
            xds_list.append(xds)
        chunksize=[1000,1000,37]
        var=['sst','view_time']
        xds_merged= reproj.mosaic_vars(xds_list,[var],chunksize)
        xds_merged['sst_mean']=xds_merged['sst'].mean(dim='t')
        xds_merged['sst_dif']=xds_merged['sst']-xds_merged['sst_mean']
        xds_dif=xds_merged[['sst_dif']]     
        print('Write Dif Netcdf')
        outfile=path_out+str(poly_id)+'_'+str(doi).zfill(3)+'_dif.nc'        
        write_nc(xds_dif, outfile)
        xds_doi=xr.Dataset(coords={'x':xds_merged.coords['x'],'y':xds_merged.coords['y'],
                                   'doi':doi})
        xds_mns=xr.concat([xds_mns,xds_doi],dim='doi')
        # Apply MannKenddal
        spco=['y','x']
        xds_res=ApplyMannKendall(xds_dif,'sst_dif',spco) 
        outfile=path_out+str(poly_id)+'_'+str(doi).zfill(3)+'_mk.nc' 
        write_nc(xds_res, outfile) 
        
   
    outfile=(path_out+str(start_date.month).zfill(2)+str(start_date.day).zfill(2)+'_' 
    +str(end_date.month).zfill(2)+str(end_date.day).zfill(2)+'_mns.nc')
    write_nc(xds_mns, outfile)




    

def main():
    # Set paths and load functions
    path='/nfs/IGARSS_2022/Composites/max/'
    path_mk='/nfs/IGARSS_2022/MannKendall/'
    path_val='/nfs/IGARSS_2022/Validationfiles/'
    shp='/nfs/IGARSS_2022/Shps/North_Sea.shp'
    aux=tl3_analysis_toolbox.aux()
    prep=tl3_analysis_toolbox.l3_lst_sst_prep()
    reproj=tl3_analysis_toolbox.reproj()
    tl_crop=tl3_analysis_toolbox.crop()
    
    # Get arguments
    args = initargs()
    filt=args.filter
    suffix=args.suffix
    
    # Find and read Composites 
    for month in range(3,10,1): #range(1,13,1)
        print('Start processing month '+str(month))
        print('Find and read Composites')
        #fs_cci=[path+ f for f in os.listdir(path) if str(month).zfill(2)+'_' in f]
        fs_cci=[path+ f for f in os.listdir(path) if (str(month).zfill(2)+'_' in f)& (int(f[0:4])>1989)]
        xds_cci=xr.open_mfdataset(fs_cci,chunks={'x': 2300, 'y': 3250},combine='nested',concat_dim=['t'])
        xds_cci=xds_cci.rename({'x': 'lon','y': 'lat'})
        xds_cci=xds_cci.sortby('t')
        # Crop composites
        print('Crop composites')
        shapefile=fiona.open(shp, "r")
        poly=shapefile[0]
        name=poly['properties']['NAME']
        crop=True
        xds_cropped=tl_crop.crop_ds_with_shp_cci(xds_cci, poly, crop)
        
        if filt=='Dif':
            xds_cropped=xds_cropped.where(abs(xds_cropped.sst_dif)<5)
        if filt=='Stdev':
            xds_cropped=filter_stdev(xds_cropped, 'sst_max')
        obs_count=np.array(xds_cropped['sst_max'].count(dim='t'))
        
        '''
        # Create csv file for each composite
        print('Create csv file for each composite')
        for day in xds_cropped.groupby('t'):
            xds_sub=day[1]
            date=str(day[0])[0:10]
            print('Process Date '+date)
            varlist=['sst_max', 'qual_max', 'tcwv_qual', 'sz_qual', 'sst_cci', 'sst_dif']
            df=aux.xds_to_df(xds_sub,varlist)
            df=df[np.isfinite(df['sst_dif'])]
            outfile=path_val+name+'_'+date+'_validationfile.csv'
            df.to_csv(outfile)
        '''
        # Compute Mann Kendall and write them to nc files
        print('Compute Mann Kendall and write them to nc files')
        #varlist=['sst_max','sst_cci']
        varlist=['sst_max']
        for var in varlist:
            print('Compute Mann Kendall for '+var)
            xds_res=ApplyMannKendall(xds_cropped,var)
            xds_res['obs_count']=(['lat','lon'],obs_count)
            fn=path_mk+name+'_'+str(month).zfill(2)+'_'+var+'_'+suffix+'_MannKendall.nc'
            write_nc(xds_res, fn)

if __name__ == '__main__':
    anomaly_trends()
    #main_kendall()
    #main()
    #analyze_mk()
    #analyze_val()
    #file_mk='E:/SST_Analysis/Test/North Sea_07_sst_max_Dif_test_MannKendall.nc'
    #xds_mk_tl=xr.open_dataset(file_mk)
    #mk_maps(xds_mk_tl,'Baltic_tl_','E:/SST_Analysis/Test/')



'''
import glob

path='E:/SST_Analysis/MannKenndal_Grid/'
grid_id='3103'
file_mk=glob.glob(path+'*'+grid_id+'*mk*.nc')
file_sst=glob.glob(path+'*'+grid_id+'*_monthly_sst*.nc')
xds_mk=xr.open_dataset(file_mk[0])
xds_sst=xr.open_dataset(file_sst[0])

xds_sst['p']=(['y','x'],xds_mk.p.data)

xds_bigp=xds_sst.where(xds_sst.p>=0.1)

big_arr=[]
for year in xds_bigp.groupby('t'):
    xds_year=year[1]
    arr=np.array(xds_year.sst_max)
    arr=arr[np.isfinite(arr)]
    big_arr.append(arr)
   
plt.boxplot(big_arr)
plt.show()

xds_smallp=xds_sst.where(xds_sst.p<0.1)

big_arr=[]
for year in xds_smallp.groupby('t'):
    xds_year=year[1]
    arr=np.array(xds_year.sst_max)
    arr=arr[np.isfinite(arr)]
    big_arr.append(arr)
   
plt.boxplot(big_arr)
plt.show()

big_arr=[]
for year in xds_sst.groupby('t'):
    xds_year=year[1]
    arr=np.array(xds_year.sst_max)
    arr=arr[np.isfinite(arr)]
    big_arr.append(arr)
plt.boxplot(big_arr)
plt.show()

big_arr=[]
for year in xds_sst.groupby('t'):
    xds_year=year[1]
    arr=np.array(xds_year.doy_max)
    arr=arr[np.isfinite(arr)]
    big_arr.append(arr)
plt.boxplot(big_arr)
plt.show()

sst_max=np.array(xds_sst.sst_max)
doy_max=np.array(xds_sst.doy_max)

sst_max=sst_max[np.isfinite(sst_max)]
doy_max=doy_max[np.isfinite(doy_max)]

doy_arr=[]
for doy in np.unique(doy_max):
    sst_doy=sst_max[doy_max==doy]
    doy_arr.append(sst_doy)
    
plt.boxplot(doy_arr)
plt.plot([0,30],[280,280])
plt.show()    
'''
'''    
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d as scipy1d

# toy coordinates and data
nx, ny, nz = 25, 30, 10
x = np.arange(nx)
y = np.arange(ny)
z = np.tile(np.arange(nz), (nx,ny,1)) + np.random.randn(nx, ny, nz)*.1
testdata = np.random.randn(nx,ny,nz) # x,y,z

# Desired z-coordinates (must be between bounds of z)
znew = np.tile(np.linspace(2,nz-2,50), (nx,ny,1)) + np.random.randn(nx, ny, 50)*0.01

# Inverse the coordinates for testing
z = z[..., ::-1]
znew = znew[..., ::-1]

# Now use own routine
ynew = interp_along_axis(testdata, z, znew, axis=2, inverse=True)

# Check some random profiles
for i in range(5):
    randx = np.random.randint(nx)
    randy = np.random.randint(ny)

    checkfunc = scipy1d(z[randx, randy], testdata[randx,randy], kind='cubic')
    checkdata = checkfunc(znew)

    fig, ax = plt.subplots()
    ax.plot(testdata[randx, randy], z[randx, randy], 'x', label='original data')
    ax.plot(checkdata[randx, randy], znew[randx, randy], label='scipy')
    ax.plot(ynew[randx, randy], znew[randx, randy], '--', label='Peter')
    ax.legend()
    plt.show()
'''
'''
def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]   

arr=np.random.randint(1, 10, size=(4, 4,4)).astype(float)
arr[1:3,1:3,1:3]=np.nan
'''



