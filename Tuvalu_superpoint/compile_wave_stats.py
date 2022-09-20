## Load packages

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
from joblib import Parallel, delayed
from collections import ChainMap
from random import sample
import os
import rpy2
import json
import warnings
warnings.filterwarnings("ignore")

transformer = \
    pyproj.Transformer.from_crs(pyproj.CRS("EPSG:32760"),pyproj.CRS("EPSG:4326")) 

## Define functions and variables

rho = 1025
g = 9.81
depth = 500 # Should be changed depending on the location
# For now we start with the deep water approximation

# Wave number k (number of waves per metre) #pg 212
def calculate_k(T,d,g,rho):
    #set error tolerance
    tol = 0.0000000001
    omega=(2*np.pi)/float(T)
    #first guess k based on deep water approximation
    k1=omega**2/g
    error = 100
    cnt = 0
    while error > tol:
        #iterate until error < tolerance
        k2=k1-((g*k1*np.tanh(k1*d)-omega**2)/(g*np.tanh(k1*d)+(g*k1*d)/np.cosh(k1*d)**2))
        error = k2 - k1;
        k1=k2;
        cnt = cnt+1
    k = k1
    return(k)

# Phase speed C0
def calculate_C0(g,k,d):
    C0 = np.sqrt((g/k) * np.tanh(k * d)) # pg 199, 202
    return(C0)

# Group velocity Cg
def calculate_Cg(C0,k,d):
    Cg = C0 * 1/2 * (1 + (2*k*d)/np.sinh(2*k*d)) # pg 199
    return(Cg)

def calc_power(p,timeframe):
    if timeframe=='month':
        year = p[0]
        month = p[1]

        time_idxs = [idx for idx,y,m in zip(range(len(times)),pd.DatetimeIndex(times).year,pd.DatetimeIndex(times).month) \
         if (y==year)&(m==month)]
    elif timeframe=='year':
        year = p
        time_idxs = [idx for idx,y in zip(range(len(times)),pd.DatetimeIndex(times).year) \
         if (y==year)]
    elif timeframe=='monthly':
        month = p
        time_idxs = [idx for idx,m in zip(range(len(times)),pd.DatetimeIndex(times).month) \
         if (m==month)]
        
    # Create a dataframe for the variables
    df = var_dict['Nanumaga'][:,:,time_idxs].to_dataframe('E').reset_index()
    
    # Calculate C?
    df['T'] = 1.0/df.freq# which is think is the period
    df['k'] = [calculate_k(T,depth,g,rho) for T in df['T']]
    df['C0'] = [calculate_C0(g,k,depth) for k in df['k']]
    df['Cg'] = [calculate_Cg(C0,k,depth) for C0,k in zip(df['C0'],df['k'])]

    # for each freq-dir combo, calculate the energy coming from x and y
    df['E_CgX'] = df.Cg*df.E*df.freq*-np.sin(np.deg2rad(df.dirr))
    df['E_CgY'] = df.Cg*df.E*df.freq*-np.cos(np.deg2rad(df.dirr))

    # # THis is how you calc power for each direction
    # df_per_dir = df.groupby('dir').sum()
    # df_per_dir['E_CgXY'] = (df_per_dir.E_CgX**2+df_per_dir.E_CgY**2)**0.5

    # Calculate power across all dir-freq combos for this timestep
    PX_sum = df['E_CgX'].sum()
    PY_sum = df['E_CgY'].sum()
    P_sum = np.sqrt(PX_sum**2+PY_sum**2)**.5
    
    PX_count = df['E_CgX'].count()
    PY_count = df['E_CgY'].count()
    P_count = np.sqrt(PX_count**2+PY_count**2)**.5
    
    PX_std = df['E_CgX'].std()
    PY_std = df['E_CgY'].std()
    P_std = np.sqrt(PX_std**2+PY_std**2)**.5
    
    results = {
                "sum":P_sum,
                "count":P_count,
                "std":P_std
            }
    
    # power is in W/m... but what is m? Is it /metre of coastline or W/metre perpendicular to the coastline?
    
    if timeframe=='month':
        return({
            tuple((year,month)):results
        })
    elif timeframe=='year':
        return({
            year:results
        })
    elif timeframe=='monthly':
        return({
            month:results
        })
    
    
## Load superpoint data

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
    
    xr_atoll = xr.DataArray(data=efth,coords=[dirr,freq,time],
                            dims=['dirr','freq','time'])
    
    var_dict.update({
        atoll:xr_atoll
    })
    print("loaded "+str(atoll))
    
    ## Load seasonal wave power

    # Get a list of all the times
    times = list(np.array(var_dict[atoll].time))

    # Find all the months and years
    months = np.arange(1,13,1)

    # Create an empty dictionary for the power at each timestep
    power_dict = {}

    # Calculate the power for each month-year combo
    # for p in itertools.product(years,months)
    power_dicts = Parallel(n_jobs=5)(delayed(calc_power)(p,'monthly') for p in months)
    power_dict = dict(ChainMap(*power_dicts))

    df_wave_power = pd.DataFrame.from_dict(power_dict,orient='index')
    df_wave_power.reset_index(inplace=True)
    df_wave_power = df_wave_power.rename(columns={'index':'month',0:'power'})
    df_wave_power = df_wave_power.sort_values('month').reset_index(drop=True)
    df_wave_power = df_wave_power.rename(columns={'index':'month'})

    df_wave_power.to_csv("{}_wave_power_monthly.csv".format(atoll),index=False)
    print('Finished monthly')

    ######## For each month over the timeseries

    # Get a list of all the times
    times = list(np.array(var_dict[atoll].time))

    # Find all the months and years
    years = np.unique(pd.DatetimeIndex(times).year)
    months = np.arange(1,13,1)

    # Create an empty dictionary for the power at each timestep
    power_dict = {}

    # Calculate the power for each month-year combo
    # for p in itertools.product(years,months)
    power_dicts = Parallel(n_jobs=12)(delayed(calc_power)(p,'month') for p in list(itertools.product(years,months)))
    power_dict = dict(ChainMap(*power_dicts))

    df_wave_power = pd.DataFrame.from_dict(power_dict,orient='index')
    df_wave_power.reset_index(inplace=True)
    try:
        df_wave_power['year'] = [x[0] for x in df_wave_power['index']]
    except:
        from IPython import embed
        embed()
    df_wave_power['month'] = [x[1] for x in df_wave_power['index']]
    df_wave_power = df_wave_power.drop('index',axis=1).rename(columns={0:'power'})
    df_wave_power = df_wave_power.sort_values(['year','month']).reset_index(drop=True)
    df_wave_power = df_wave_power.drop('month',axis=1).reset_index()
    df_wave_power = df_wave_power.rename(columns={'index':'months'})

    df_wave_power.to_csv("{}_wave_power_per_month.csv".format(atoll),index=False)
    print('Finished per month')

    #### Per annum

    # Get a list of all the times
    times = list(np.array(var_dict[atoll].time))

    # Find all the months and years
    years = np.unique(pd.DatetimeIndex(times).year)

    # Create an empty dictionary for the power at each timestep
    power_dict = {}

    # Calculate the power for each month-year combo
    # for p in itertools.product(years,months)
    power_dicts = Parallel(n_jobs=12)(delayed(calc_power)(p,'year') for p in years)
    power_dict = dict(ChainMap(*power_dicts))

    df_wave_power = pd.DataFrame.from_dict(power_dict,orient='index')
    df_wave_power.reset_index(inplace=True)
    df_wave_power = df_wave_power.rename(columns={'index':'year',0:'power'})
    df_wave_power = df_wave_power.sort_values('year').reset_index(drop=True)
    df_wave_power = df_wave_power.rename(columns={'index':'months'})

    df_wave_power.to_csv("{}_wave_power_per_annum.csv".format(atoll),index=False)

    print('Finished per annum')