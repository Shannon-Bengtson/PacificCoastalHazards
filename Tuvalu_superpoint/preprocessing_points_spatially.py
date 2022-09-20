import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import stat
from datetime import datetime
import itertools
from netCDF4 import Dataset
import sys
import scipy as sp
import scipy.io as sio
import xarray as xr
import math

def preprocessing_points_spatially():
    '''
    
    '''
    # Load wave data
    waves_dict = sio.loadmat('/src/Dataset/D8_tarawa_inundation/Waves_around_tarawa.mat')

    # Load wind data and do some basic preprocessing
    winds_dict = sio.loadmat('/src/Dataset/D8_tarawa_inundation/ERA5_Winds_Tarawa.mat')
    winds_dict = {x[0]:(x[1],x[2]) for x in winds_dict['data']}
    winds_array = [winds_dict[time[0]] if time[0] in winds_dict else winds_dict[min(winds_dict.keys(), key=lambda k: abs(k-time[0]))]\
             for time in waves_dict['Timeo']]
    time_array = [time[0] if time[0] in winds_dict else winds_dict[min(winds_dict.keys(), key=lambda k: abs(k-time[0]))]\
             for time in waves_dict['Timeo']]## new 
    
    # Load the inundation data
    inundation_dict = sio.loadmat('/src/Dataset/D8_tarawa_inundation/Historical_inundation_Risk.mat')

    # Load the lagoon and ocean side profiles 
    df_lagoon_profiles = pd.read_csv(
        '/src/Dataset/D8_tarawa_inundation/Profiles_definition_inner_lagoon_xyxy.txt',delim_whitespace=True,header=None)
    df_lagoon_profiles.columns = ['shore_long','shore_lat','reef_long','reef_lat','reef_depth']
    df_ocean_profiles = pd.read_csv(
        '/src/Dataset/D8_tarawa_inundation/Profiles_definition_outer_reef_xyxy.txt',delim_whitespace=True,header=None)
    df_ocean_profiles.columns = ['shore_long','shore_lat','reef_long','reef_lat','reef_width','reef_depth','forereef_slope']

    # Calculate the lat-long gradient between the reef and the shore
    df_ocean_profiles['m'] = (df_ocean_profiles.reef_lat-df_ocean_profiles.shore_lat)/(df_ocean_profiles.reef_long-df_ocean_profiles.shore_long)

    # Calculate the direction that the shoreline faces using the lat-long gradient
    phi_list = []
    for index,row in df_ocean_profiles.iterrows():
        if (row.reef_lat>row.shore_lat)&(row.reef_long>row.shore_long):
            phi = 90-180*math.atan(row.m)/math.pi
        elif (row.reef_lat<row.shore_lat)&(row.reef_long>row.shore_long):
            phi = 90+180*math.atan(-row.m)/math.pi
        elif (row.reef_lat<row.shore_lat)&(row.reef_long<row.shore_long):
            phi = 270-180*math.atan(row.m)/math.pi
        elif (row.reef_lat>row.shore_lat)&(row.reef_long<row.shore_long):
            phi = 270+180*math.atan(-row.m)/math.pi
        else:
            print('Halp')

        phi_list.append(phi)

    # add gradient to dictionary
    df_ocean_profiles['shore_dir'] = phi_list

    # Create dictionaries to put data into. Each location will have one key, and one combined dataframe of all the variables for the model in it
    lagoon_twl_and_wave_dict = {}
    ocean_twl_and_wave_dict = {}

    # Loop over each position, split into ocean side and into shoreside, create df, and add to dictionary
    for position,position_index in zip(inundation_dict['Ptos'],np.arange(0,len(inundation_dict['Ptos']),1)):

        waves_subset_dict = {key:value[:,position_index] for key,value in waves_dict.items() if key in ['Diro','Hso','Tmo']}
        inundation_subset_dict = {key:value[:,position_index] for key,value in inundation_dict.items() if key in \
                 ['TWL']}
        tide_dict = {'Tide':inundation_dict['Tide'].reshape(inundation_dict['Tide'].shape[0],)}
        sla_dict = {'MSL':inundation_dict['sladac2'].reshape(inundation_dict['sladac2'].shape[0],)}
        winds_dict = {'wind_u':[x[0] for x in winds_array],'wind_v':[x[1] for x in winds_array]}
        time_dict = {'time':time_array} # new

        BN_vars_dict = {**inundation_subset_dict,**waves_subset_dict,**tide_dict,**sla_dict,**winds_dict,**time_dict}

        BN_vars_dict['Hs_offshore'] = BN_vars_dict.pop('Hso')
        BN_vars_dict['Tm_offshore'] = BN_vars_dict.pop('Tmo')
        BN_vars_dict['Dir_offshore'] = BN_vars_dict.pop('Diro') 

        if not df_lagoon_profiles[df_lagoon_profiles.shore_long==position[0]].empty:

            df_position_characteristics = df_lagoon_profiles[df_lagoon_profiles.shore_long==position[0]].reset_index(drop=True)
            pos_dict = {'long':[df_position_characteristics.loc[0,'shore_long']]*len(tide_dict['Tide']),
                             'lat':[df_position_characteristics.loc[0,'shore_lat']]*len(tide_dict['Tide'])}

            BN_vars_position_dict = {**BN_vars_dict,**pos_dict}

            lagoon_twl_and_wave_dict.update({
                tuple(position):pd.DataFrame(BN_vars_position_dict).sample(1000)
            })

        elif not df_ocean_profiles[df_ocean_profiles.shore_long==position[0]].empty:

            # Add reef characteristicsdf_lagoon_profiles
            df_position_characteristics = df_ocean_profiles[df_ocean_profiles.shore_long==position[0]].reset_index(drop=True)
            reef_width_dict = {'reef_width':[df_position_characteristics.loc[0,'reef_width']]*len(tide_dict['Tide'])}
            reef_depth_dict = {'reef_depth':[df_position_characteristics.loc[0,'reef_depth']]*len(tide_dict['Tide'])}
            forereef_slope_dict = {'forereef_slope':[df_position_characteristics.loc[0,'forereef_slope']]*len(tide_dict['Tide'])}
            pos_dict = {'long':[df_position_characteristics.loc[0,'shore_long']]*len(tide_dict['Tide']),
                             'lat':[df_position_characteristics.loc[0,'shore_lat']]*len(tide_dict['Tide']),
                             'shore_dir':[df_position_characteristics.loc[0,'shore_dir']]*len(tide_dict['Tide'])}

            BN_vars_position_dict = {**BN_vars_dict,**reef_width_dict,**reef_depth_dict,**forereef_slope_dict, **pos_dict}

            ocean_twl_and_wave_dict.update({
                tuple(position):pd.DataFrame(BN_vars_position_dict).sample(1000)
            })
        else:
            print('Error! Neither in lagoon nor ocean?')

    # create two combined dataframes from the two dictionaries
    df_lagoon = pd.concat(lagoon_twl_and_wave_dict).dropna()
    df_ocean = pd.concat(ocean_twl_and_wave_dict).dropna()
    
    # Reset in the dataframe indexes
    df_lagoon.reset_index(drop=True,inplace=True)
    df_ocean.reset_index(drop=True,inplace=True)

    # Save the dataframes
    return(df_ocean,df_lagoon)
