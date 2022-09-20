# Load required modules
import pandas as pd 
from datetime import datetime
import numpy as np
import json

# def BN_Antonio_preprocessing_lagoon(df_Pto_lagoon):
#     '''
#         Function used for preprocessing the data provided by Antonio (DS2)
#         Reads in the files (lagoon side) and formates the data
#         Adds MEI data (ENSO indicator) to the file (DS5)
        
#         Returns:
#          - Two pandas.DataFrames; lagoon and ocean side data
#     '''

#     # Remove the tide from the total water level in both files
#     # (since there is a very obvious relationship between them)
#     df_Pto_lagoon['WL_wave_comp'] = df_Pto_lagoon['TWL']-df_Pto_lagoon['Tide']-df_Pto_lagoon['MSL'] 
    
#     # Get variables for SM1
#     wind_u_array = np.array(df_Pto_lagoon['wind_u'])
#     wind_v_array = np.array(df_Pto_lagoon['wind_v'])
#     Hs_offshore_array = np.array(df_Pto_lagoon['Hs_offshore'])
#     Tm_offshore_array = np.array(df_Pto_lagoon['Tm_offshore'])
#     Dir_offshore_array = np.array(df_Pto_lagoon['Dir_offshore'])
#     Tide_array = np.array(df_Pto_lagoon['Tide'])
#     TWL_array = np.array(df_Pto_lagoon['TWL'])
#     WL_wave_comp = np.array(df_Pto_lagoon['WL_wave_comp'])
#     MSL_array = np.array(df_Pto_lagoon['MSL'])
    
#     # Put all variables for BN into a dictionary
#     variables_dict = {
#         'wind_u':wind_u_array,
#         'wind_v':wind_v_array,
#         'Hs_offshore':Hs_offshore_array,
#         'Tm_offshore':Tm_offshore_array,
#         'Dir_offshore':Dir_offshore_array,
# #         'Tide':Tide_array,
# #         'TWL':TWL_array,
#         'WL_wave_comp':WL_wave_comp,
# #         'MSL':MSL_array,
#     }
    
#     return(df_Pto_lagoon,variables_dict)

# def BN_Antonio_preprocessing_ocean(df_Pto_ocean):
#     '''
#         Function used for preprocessing the data provided by Antonio (DS2)
#         Reads in the files (ocean side) and formates the data
#         Adds MEI data (ENSO indicator) to the file (DS5)
        
#         Returns:
#          - Two pandas.DataFrames; lagoon and ocean side data
#     '''

#     # Remove the tide from the total water level in both files
#     # (since there is a very obvious relationship between them)
#     df_Pto_ocean['WL_wave_comp'] = df_Pto_ocean['TWL']-df_Pto_ocean['Tide']-df_Pto_ocean['MSL'] 

# #     ### Add DS5 (ENSO)
# #     # Load MEI (ENSO indicator)
# #     data_location = '/mnt/c/Users/shannonb/Documents/Model_and_data/Dataset/D5_ENSO/'
# #     file_name = "MEI_preprocessed"
# #     with open("{}{}.json".format(data_location,file_name), 'r') as fp:
# #         MEI_dict = json.load(fp)
# #     # Add MEI to dataframe
# #     df_Pto_ocean['MEI'] = [float(MEI_dict[str(x.month)][str(x.year)]) for x in df_Pto_ocean.Time]
    
#     TWL_array = np.array(df_Pto_ocean['TWL'])
#     MSL_array = np.array(df_Pto_ocean['MSL'])
#     Tide_array = np.array(df_Pto_ocean['Tide'])
#     wind_u_array = np.array(df_Pto_ocean['wind_u'])
#     wind_v_array = np.array(df_Pto_ocean['wind_v'])
#     Hs_offshore_array = np.array(df_Pto_ocean['Hs_offshore'])
#     Tm_offshore_array = np.array(df_Pto_ocean['Tm_offshore'])
#     Dir_offshore_array = np.array(df_Pto_ocean['Dir_offshore'])
#     reef_width_array = np.array(df_Pto_ocean['reef_width'])
#     reef_depth_array = np.array(df_Pto_ocean['reef_depth'])
#     forereef_slope_array = np.array(df_Pto_ocean['forereef_slope'])
#     WL_wave_comp = np.array(df_Pto_ocean['WL_wave_comp'])
#     shore_dir_array = np.array(df_Pto_ocean['shore_dir'])
    
#     # Put all variables for BN into a dictionary
#     variables_dict = {
# #         'TWL':TWL_array,
#         'MSL':MSL_array,
#         'Tide':Tide_array,
# #         'Wind_u':wind_u_array,
# #         'Wind_v':wind_v_array,
#         'Hs_offshore':Hs_offshore_array,
#         'Tm_offshore':Tm_offshore_array,
#         'Dir_offshore':Dir_offshore_array,
#         'WL_wave_comp':WL_wave_comp,
#         'reef_width':reef_width_array,
#         'reef_depth':reef_depth_array,
#         'forereef_slope':forereef_slope_array,
#         'shore_dir':shore_dir_array
#     }
#     return(df_Pto_ocean,variables_dict)



# # Load required modules
# import pandas as pd 
# from datetime import datetime
# import numpy as np
# import json

def BN_Antonio_preprocessing_lagoon(df_Pto_lagoon):
    '''
        Function used for preprocessing the data provided by Antonio (DS2)
        Reads in the files (lagoon side) and formates the data
        Adds MEI data (ENSO indicator) to the file (DS5)
        
        Returns:
         - Two pandas.DataFrames; lagoon and ocean side data
    '''

    # Remove the tide from the total water level in both files
    # (since there is a very obvious relationship between them)
    df_Pto_lagoon['TWL_less_Tide'] = df_Pto_lagoon['Tide']-df_Pto_lagoon['TWL']
    
    # Get variables for SM1
    wind_u_array = np.array(df_Pto_lagoon['wind_u'])
    wind_v_array = np.array(df_Pto_lagoon['wind_v'])
    Hs_offshore_array = np.array(df_Pto_lagoon['Hs_offshore'])
    Tm_offshore_array = np.array(df_Pto_lagoon['Tm_offshore'])
    Dir_offshore_array = np.array(df_Pto_lagoon['Dir_offshore'])
    Tide_array = np.array(df_Pto_lagoon['Tide'])
    TWL_array = np.array(df_Pto_lagoon['TWL'])
    TWL_less_Tide_array = np.array(df_Pto_lagoon['TWL_less_Tide'])
    MSL_array = np.array(df_Pto_lagoon['MSL'])
    
    # Put all variables for BN into a dictionary
    variables_dict = {
        'wind_u':wind_u_array,
        'wind_v':wind_v_array,
        'Hs_offshore':Hs_offshore_array,
        'Tm_offshore':Tm_offshore_array,
        'Dir_offshore':Dir_offshore_array,
        'Tide':Tide_array,
        'TWL':TWL_array,
        'TWL_less_Tide':TWL_less_Tide_array,
        'MSL':MSL_array,
    }
    
    return(df_Pto_lagoon,variables_dict)

def BN_Antonio_preprocessing_ocean(df_Pto_ocean):
    '''
        Function used for preprocessing the data provided by Antonio (DS2)
        Reads in the files (ocean side) and formates the data
        Adds MEI data (ENSO indicator) to the file (DS5)
        
        Returns:
         - Two pandas.DataFrames; lagoon and ocean side data
    '''

    # Remove the tide from the total water level in both files
    # (since there is a very obvious relationship between them)
    df_Pto_ocean['TWL_less_Tide'] = df_Pto_ocean['Tide']-df_Pto_ocean['TWL']

#     ### Add DS5 (ENSO)
#     # Load MEI (ENSO indicator)
#     data_location = '/mnt/c/Users/shannonb/Documents/Model_and_data/Dataset/D5_ENSO/'
#     file_name = "MEI_preprocessed"
#     with open("{}{}.json".format(data_location,file_name), 'r') as fp:
#         MEI_dict = json.load(fp)
#     # Add MEI to dataframe
#     df_Pto_ocean['MEI'] = [float(MEI_dict[str(x.month)][str(x.year)]) for x in df_Pto_ocean.Time]
    
    TWL_array = np.array(df_Pto_ocean['TWL'])
    MSL_array = np.array(df_Pto_ocean['MSL'])
    Tide_array = np.array(df_Pto_ocean['Tide'])
    wind_u_array = np.array(df_Pto_ocean['wind_u'])
    wind_v_array = np.array(df_Pto_ocean['wind_v'])
    Hs_offshore_array = np.array(df_Pto_ocean['Hs_offshore'])
    Tm_offshore_array = np.array(df_Pto_ocean['Tm_offshore'])
    Dir_offshore_array = np.array(df_Pto_ocean['Dir_offshore'])
    reef_width_array = np.array(df_Pto_ocean['reef_width'])
    reef_depth_array = np.array(df_Pto_ocean['reef_depth'])
    forereef_slope_array = np.array(df_Pto_ocean['forereef_slope'])
    TWL_less_Tide_array = np.array(df_Pto_ocean['TWL_less_Tide'])
    shore_dir_array = np.array(df_Pto_ocean['shore_dir'])
    
    # Put all variables for BN into a dictionary
    variables_dict = {
        'TWL':TWL_array,
        'MSL':MSL_array,
        'Tide':Tide_array,
#         'Wind_u':wind_u_array,
#         'Wind_v':wind_v_array,
        'Hs_offshore':Hs_offshore_array,
        'Tm_offshore':Tm_offshore_array,
        'Dir_offshore':Dir_offshore_array,
        'TWL_less_Tide':TWL_less_Tide_array,
        'reef_width':reef_width_array,
        'reef_depth':reef_depth_array,
        'forereef_slope':forereef_slope_array,
        'shore_dir':shore_dir_array
    }
    return(df_Pto_ocean,variables_dict)

