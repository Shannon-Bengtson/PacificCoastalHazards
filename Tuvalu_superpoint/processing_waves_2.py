import pandas as pd
import numpy as np
import netCDF4
import datetime
import itertools
import xarray as xr
from joblib import Parallel, delayed
from collections import ChainMap
import json 

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
    
def averaging_waves(netcdf_file_name,atoll):

    # Load MEI data
    df_MEI = pd.DataFrame.from_dict(json.load(open('../D5_ENSO/MEI_preprocessed.json')),orient='index').T
    df_MEI = df_MEI.reset_index().melt('index').rename(columns={'index':'year','variable':'month','value':'MEI'})
    df_MEI['year'] = df_MEI.year.astype(int)
    df_MEI['month'] = df_MEI.month.astype(int)
    df_MEI['date'] = [pd.to_datetime("01/{}/{}".format(month,year)) for month,year in zip(df_MEI.month,df_MEI.year)]
    
    # Load the data
    ds = netCDF4.Dataset(netcdf_file_name)

    # Extract the variables
    time = np.array(ds.variables['time'])
    
    # Adjust time to be datetime (need to confirm the start time)
    time_start = datetime.datetime(1980,1,1,0,0)
    time = [(time_start+datetime.timedelta(hours=x)) for x in np.array(ds.variables['time'])]

    xr_test = xr.DataArray(data=np.array(ds.variables['efth']),
                           coords=[np.array(ds.variables['dir']),
                                   np.array(ds.variables['freq']),
                                   time],
                           dims=['dirr','freq','time'])
    
    T = 1/np.array(ds.variables['freq'])
    K = [calculate_k(t,depth,g,rho) for t in T]
    C0 = [calculate_C0(g,k,depth) for k in K]
    Cg = [calculate_Cg(c0,k,depth) for c0,k in zip(C0,K)]

    xr_test['freq1'] = xr_test.freq.copy()
    xr_test = xr_test.assign_coords({'freq':np.array(Cg)})
    xr_test = xr_test.rename({'freq':'Cg'})    
    
    df_clusters = pd.read_csv('cluster_averages.csv')
    for index,row in df_clusters.iterrows():
        
        if row.left==1:
            xr_left = xr_test.where(xr_test.dirr>(row.shoreline_direction-90))
            xr_left = xr_left.where(xr_left.dirr<(row.shoreline_direction-30))
        else:
            xr_left = xr_test.copy()

        if row.mid==1:
            xr_mid = xr_test.where(xr_test.dirr>(row.shoreline_direction-30))
            xr_mid = xr_mid.where(xr_mid.dirr<(row.shoreline_direction+30))
        else:
            xr_mid = xr_test.copy()
        
        if row.right==1:
            xr_right = xr_test.where(xr_test.dirr>(row.shoreline_direction+30))
            xr_right = xr_right.where(xr_right.dirr<(row.shoreline_direction+90))
        else:
            xr_right = xr_test.copy()
            
        xr_left = xr_left.fillna(0)
        xr_mid = xr_mid.fillna(0)
        xr_right = xr_right.fillna(0)
    
#         from IPython import embed
#         embed()

        for direction,xr_direction in zip(['left','mid','right'],[xr_left,xr_mid,xr_right]):
        
            print('preprocessing this combo')
            print(index,direction)
            xr_E_CgX = xr_direction*xr_direction.Cg*xr_direction.freq1*-np.sin(np.deg2rad(xr_direction.dirr))
            xr_E_CgY = xr_direction*xr_direction.Cg*xr_direction.freq1*-np.cos(np.deg2rad(xr_direction.dirr))

        #     # Get a list of all the times
        #     times = [str(x) for x in list(np.array(xr_E_CgX.time))]

        #     with open(r'time_{}.txt'.format(atoll), 'w') as fp:
        #         for time in times:
        #             fp.write(time+'\n')

        #     xr_E_CgX.to_csv('xr_E_CgX_{}.csv'.format(atoll),index=True)
        #     xr_E_CgY.to_csv('xr_E_CgY_{}.csv'.format(atoll),index=True)

            xr_E_CgX.to_netcdf(path='xr_E_CgX_{}_{}_{}.nc'.format(atoll,index,direction),engine='netcdf4')
            xr_E_CgY.to_netcdf(path='xr_E_CgY_{}_{}_{}.nc'.format(atoll,index,direction),engine='netcdf4')

        
def calc_wave_power(year_1,year_2,json_file_name,atoll,clust,direction):

#     xr_E_CgX = pd.read_csv('xr_E_CgX_{}.csv'.format(atoll),index=True)
#     xr_E_CgY = pd.read_csv('xr_E_CgY_{}.csv'.format(atoll),index=True)

    xr_E_CgX = xr.open_dataset('xr_E_CgX_{}_{}_{}.nc'.format(atoll,clust,direction),engine='netcdf4')
    xr_E_CgY = xr.open_dataset('xr_E_CgY_{}_{}_{}.nc'.format(atoll,clust,direction),engine='netcdf4')
    
    print(year_1)
    print(clust)
    print(direction)

    xr_nanumaga_Y = xr_E_CgY.copy()
    xr_nanumaga_X = xr_E_CgX.copy()

    time_filter = (xr_nanumaga_Y.time>pd.to_datetime(year_1,format='%d/%m/%Y'))&(xr_nanumaga_Y.time<pd.to_datetime(year_2,format='%d/%m/%Y'))

    xr_nanumaga_Y = xr_nanumaga_Y.where(time_filter,drop=True)   
    xr_nanumaga_X = xr_nanumaga_X.where(time_filter,drop=True)
    
    xr_nanumaga_X = xr_nanumaga_X.rename({'__xarray_dataarray_variable__':'P'})
    xr_nanumaga_Y = xr_nanumaga_Y.rename({'__xarray_dataarray_variable__':'P'})

    df_nanumaga_X = xr_nanumaga_X.to_dataframe().drop('freq1',axis=1)
    df_nanumaga_Y = xr_nanumaga_Y.to_dataframe().drop('freq1',axis=1)

    P = (df_nanumaga_X.P.sum()**2+df_nanumaga_Y.P.sum()**2)**0.5
    P_mean = P/(pd.to_datetime(year_2)-pd.to_datetime(year_1)).days
    P_median = np.median(P)
    P_stdev = np.std(P)

#     # Get mean MEI for this period
#     MEI_values = df_MEI[(df_MEI.date>=year_1)&(df_MEI.date<=year_2)].MEI.values
#     mean_MEI = np.mean([float(x) for x in MEI_values])

    P_dirr = ((df_nanumaga_X.groupby('dirr').sum().P**2+df_nanumaga_Y.groupby('dirr').sum().P**2)**0.5).to_dict()

    P_per_hours = (df_nanumaga_X.P**2+df_nanumaga_Y.P**2)**0.5

    def x_size_event(percent_low,percent_high):

        minn = (np.max(P_per_hours)-np.min(P_per_hours))*percent_low/100
        maxx = (np.max(P_per_hours)-np.min(P_per_hours))*percent_high/100

        events = P_per_hours[(P_per_hours>=minn)&(P_per_hours<maxx)]

        events_dict = pd.DataFrame(events).reset_index()[['time','P']].set_index('time').to_dict()
        events_dict = {str(k):v for k,v in events_dict['P'].items()}

        return(events_dict)

    percentages = [50,75,80,85,90,100]
    events_dict = {str((percent_low,percent_high)):x_size_event(percent_low,percent_high) for percent_low,percent_high in zip(percentages[:-1],percentages[1:])}

    total_wave_power = {
            'total':P,
            'daily_mean':P_mean,
            'median_wave':P_median,
            'std_wave':P_stdev,
            'P_per_dirr':P_dirr,
#             'mean_MEI':mean_MEI,
            'events':events_dict
        }

    return(total_wave_power)






