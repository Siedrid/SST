# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 15:37:06 2023

@author: rein_pp
"""
import fiona
from shapely.geometry import shape
import os
import pandas as pd
import pdb


path_seas='E:/Conferences_Presentations/Strukturkommision_2022/Folien/\
World_Seas_IHO_v3/World_Seas_IHO_v3/'
path_grid='E:/l3_lst_orbit_drift_corr/julien2022_correction/Data/grid/'
path_results='E:/Publications/SST_analysis/Results/V1/'
path_out='E:/Publications/SST_analysis/to_process/'
#pdb.set_trace()
breakpoint()
sea_file=path_seas+'World_Seas_IHO_v3.shp'
grid_file=path_grid+'grid.shp'
avail_list=[float(file[0:4]) for file in os.listdir(path_results)]


grid_shp=fiona.open(grid_file)
sea_shp=fiona.open(sea_file)
study_areas=[['Balearic (Iberian Sea)'],['Adriatic Sea'],['Skagerrak', 'Kattegat',
              'Baltic Sea', 'North Sea'], ['Aegean Sea' , 'Sea of Marmara']]
output=dict()

for study_area in study_areas:
    id_list=[]
    existing=[]
    left=[]
    top=[]
    right=[]
    bottom=[]
    for sea_id in study_area:
        print(sea_id)
        sea_poly=[poly for poly in sea_shp if poly['properties']['NAME']==sea_id][0]
        sea_geom=shape(sea_poly['geometry'])
        for grid_poly in grid_shp:
            grid_geom=shape(grid_poly['geometry'])
            if grid_geom.intersection(sea_geom).area >0:
                grid_id=grid_poly['properties']['id']
                id_list.append(grid_id)
                left.append(grid_poly['properties']['left'])
                top.append(grid_poly['properties']['top'])
                right.append(grid_poly['properties']['right'])
                bottom.append(grid_poly['properties']['bottom'])
                if grid_id in avail_list:
                    existing.append('yes')
                else:
                    existing.append('no')
    df_out=pd.DataFrame({'ID':id_list,'existing':existing,'left':left,'right':right,
                         'bottom':bottom,'top':top})
    if study_area[0] == 'Skagerrak':
        df_out=df_out[(df_out.left>5.4)&(df_out.right<18.6)&(df_out.bottom>52.4)&(df_out.top<62.6)]
    output[study_area[0]]=df_out
    

    
for key in output.keys():
    outfile=path_out+key+'.csv'
    output[key].to_csv(outfile)            
        
        
        
