# Initialise
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import stat
import math
import geopandas as gpd
import netCDF4
import datetime
import itertools
import shapely 
from shapely.geometry import LineString, shape
import folium
from colormap import rgb2hex
from folium.plugins import FloatImage
from scipy import interpolate
import alphashape
import descartes
import pyproj
import xarray as xr
import json

transformer = \
    pyproj.Transformer.from_crs(pyproj.CRS("EPSG:32760"),pyproj.CRS("EPSG:4326")) 

import os
import rpy2
os.environ['R_HOME'] = '/lib/R'

def efth_stat_calc(df):
    '''
        Function for grouping each dataframe into freq bins and calc stats
    '''
    bin_stats_dict = {}
    for lower,upper in zip(freq_vals_bin_edges[:-1],freq_vals_bin_edges[1:]):
        df_freq_bin = df[(df.freq>=lower)&(df.freq<upper)]
        bin_max = np.max(df_freq_bin.efth)
        bin_mean = np.mean(df_freq_bin.efth)

        bin_stats_dict.update({
            tuple((lower,upper)):{
                'max':bin_max,
                'mean':bin_mean
            }
        })

    return(bin_stats_dict)

#########################################################################################################
####### Load Superpoint Data ############################################################################
#########################################################################################################

# Load the data
ds_nanumaga = netCDF4.Dataset('SuperPoint_Nanumaga.nc')
ds_nanumea = netCDF4.Dataset('SuperPoint_Nanumea.nc')

var_dict = {}

for ds,atoll in zip([ds_nanumaga,ds_nanumea],['Nanumaga','Nanumea']):
    # Extract the variables
    efth = np.array(ds.variables['efth'])
    time = np.array(ds.variables['time'])
    dirr = np.array(ds.variables['dir'])
    freq = np.array(ds.variables['freq'])
    wdir = np.array(ds.variables['Wdir'])
    wspd = np.array(ds.variables['Wspeed'])

    # Adjust time to be datetime (need to confirm the start time)
    time_start = datetime.datetime(1980,1,1,0,0)
    time = [(time_start+datetime.timedelta(hours=x)) for x in time]
    
    var_dict.update({
        atoll:{
            'efth':efth,
            'time':time,
            'dirr':dirr,
            'freq':freq,
            'wdir':wdir,
            'wspd':wspd
        }
    })
    
#########################################################################################################
####### Load Shoreline Data #############################################################################
#########################################################################################################

# Load the data using geopandas

proxies = [
    r'TOB',
    r'VL',
    r'WM'
]

atolls = [
    'Nanumea',
    'Nanumaga'
]

locations_dict = {
    'Nanumea':[-5.667723, 176.094928],
    'Nanumaga':[-6.287944, 176.321295]
}

combinations = list(itertools.product(atolls,proxies))

# Define the years the shoreline change datafile is for. This is atoll specific
years_dict = {
    'Nanumea':'1971_2020',
    'Nanumaga':'2003_2020'
}

geopandas_dict = {}

for combination in combinations:
    atoll = combination[0]
    proxy = combination[1]
    year = years_dict[combination[0]]
    geopandas_dict.update({
        (atoll,proxy):gpd.read_file('../D9_Tuvalu_shoreline//Shoreline_Definition_Shp/{}/{}_{}_{}.shp'.format(proxy,atoll,proxy,year))
    })


# converting the geopandas dataframes into a useable pd.dataframe
    
waves_per_atoll_year_dict = {}

for combination in combinations:
    atoll = combination[0]
    proxy = combination[1]
#     year = years_dict[combination[0]]
    
    gdf_test = geopandas_dict[atoll,proxy].copy()

    # Correct some typos
    gdf_test.loc[gdf_test.layer=='Nanumea_WM_2015_','layer'] = 'Nanumea_WM_2015'

    gdf_test['layer'] = gdf_test['layer'].fillna(value=1)
    gdf_test['id'] = gdf_test['id'].fillna(value=1)
    gdf_test.dropna(inplace=True)

    years_list = []
    gdf_right_years_dict = {}

    for layer,group in gdf_test.groupby('layer'):
        if (proxy=='VL')&(atoll!='Nanumaga'):
            year = int(layer.split('_')[1])
        else:
            year = int(layer.split('_')[-1])
        years_list.append(year)

        gdf_right_years_dict.update({
            year:group
        })

    gdf_test = pd.concat(gdf_right_years_dict)
    # gdf_test = gdf_test.rename_axis('year')
    gdf_test.reset_index(drop=False,inplace=True)

    # Convert polygons to linestrings if there are any
    if type(gdf_test.loc[0,'geometry'])==shapely.geometry.Polygon:
        gdf_test['geometry'] = [x.boundary for x in gdf_test.geometry]
        
    # Interpolate the linestrings so that they are all of the same length (x)
    x = 100 #length of interpolation
    for i in np.arange(0,len(gdf_test),1):
        gdf_test.loc[i,'geometry'] = \
            shapely.geometry.linestring.LineString(
                [gdf_test.loc[i,'geometry'].interpolate((j/x), normalized=True) for j in range(1, x)]
            )

    dict_of_df_xy = {}

    for idx,row in gdf_test.iterrows():
        linestring = row.geometry
        XY_list = []
        
        if type(linestring)==shapely.geometry.linestring.LineString:
            XY_list = XY_list+[(x,y) for x,y in linestring.coords]
        elif type(linestring)==shapely.geometry.multilinestring.MultiLineString:
            XY_list = XY_list+[(x,y) for x,y in linestring[0].coords]
                
        df_xy = pd.DataFrame(XY_list)
        df_xy.columns = ['lon','lat']
        df_xy['id'] = int(row.id)
        df_xy['year'] = int(row.level_0)
#         print(row.id)
        
        dict_of_df_xy.update({
            idx:df_xy
        })

    df_xy = pd.concat(dict_of_df_xy)

    df_xy = df_xy.reset_index(drop=True)

    df_xy[('lat,lon')] = [transformer.transform(x,y) for x,y in zip(df_xy.lon,df_xy.lat)]
    df_xy['lon'] = [x[1] for x in df_xy[('lat,lon')]]
    df_xy['lat'] = [x[0] for x in df_xy[('lat,lon')]]

    # Get all the years
    years = np.sort(np.unique(df_xy.year))    
    years_beginning = years[:-1]
    years_end = years[1:]

    # Loop over the years, defining the beginning and the end year
    for year_beginning,year_end in zip(years_beginning,years_end):

        beginning = datetime.datetime(year_beginning,1,1,0,0)
        end = datetime.datetime(year_end,12,31,23,59)
        beginning_idx = min(range(len(time)), key=lambda i: abs(time[i]-beginning))
        end_idx = min(range(len(time)), key=lambda i: abs(time[i]-end))
        
        if (atoll,proxy,year_beginning,year_end)==('Nanumea', 'TOB', 1971, 1984):
            continue

        #################################################
        ## get the wave energy for the right time, place
        #################################################

        xr_atoll = xr.DataArray(data=var_dict[atoll]['efth'],coords=[var_dict[atoll]['dirr'],var_dict[atoll]['freq'],var_dict[atoll]['time']],
                                dims=['dirr','freq','time'])
        xr_atoll = xr_atoll[:,:,beginning_idx:end_idx]
        df_xr_atoll = xr_atoll.to_dataframe('efth').reset_index(drop=False)
        
        print(df_xr_atoll.head)
        
        # Normalise all values back to 0->360 degress
        df_xr_atoll.dirr= df_xr_atoll.dirr-np.min(df_xr_atoll.dirr)

        # Create bins by quadrants
        df_xr_NE = df_xr_atoll[df_xr_atoll.dirr<=90]
        df_xr_ES = df_xr_atoll[(df_xr_atoll.dirr>90)&(df_xr_atoll.dirr<=180)]
        df_xr_SW = df_xr_atoll[(df_xr_atoll.dirr>180)&(df_xr_atoll.dirr<=270)]
        df_xr_WN = df_xr_atoll[(df_xr_atoll.dirr>270)&(df_xr_atoll.dirr<=360)]

        # Bin the frequency bins
        freq_vals = np.unique(df_xr_atoll.freq)
        freq_vals_bin_edges = freq_vals[::int(round(len(freq_vals)/2.5,0))]
        freq_vals_bin_edges = list(freq_vals_bin_edges)+[np.max(df_xr_ES.freq)]
        
        print(freq_vals_bin_edges)

        # Get stats of efth
        NE_bin_stats_dict = efth_stat_calc(df_xr_NE)
        ES_bin_stats_dict = efth_stat_calc(df_xr_ES)
        SW_bin_stats_dict = efth_stat_calc(df_xr_SW)
        WN_bin_stats_dict = efth_stat_calc(df_xr_WN)
        
        
        waves_per_atoll_year_dict.update({
            tuple((atoll,proxy,year_beginning,year_end)):{
                'NE':NE_bin_stats_dict,
                'ES':ES_bin_stats_dict,
                'SW':SW_bin_stats_dict,
                'WN':WN_bin_stats_dict
            }
        })
        print((atoll,proxy,year_beginning,year_end))
        
# Save the data
with open('data.json', 'w') as fp:
    json.dump(waves_per_atoll_year_dict, fp)
    
        
