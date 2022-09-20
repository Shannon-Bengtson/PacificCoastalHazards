import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from processing_waves_2 import averaging_waves,calc_wave_power
import json
import datetime

atolls = ['Nanumea']

for atoll in atolls:

    dates = list(pd.read_csv('{}_dates.tsv'.format(atoll),sep='\t'))
    dates = \
        [datetime.datetime(int(x.split('/')[2]),int(x.split('/')[1]),int(x.split('/')[0])) for x in dates][:-1]
    
    averaging_waves('SuperPoint_{}.nc'.format(atoll),atoll)

    for clust in [1,2,3,4,5]:
        for direction in ['left','mid','right']:
            total_wave_power_list = list(map(calc_wave_power,
                                             dates[:-1],
                                             dates[1:],
                                            ['digitised_averages_{}'.format(atoll)]*len(dates[1:]),
                                            [atoll]*len(dates[1:]),
                                            [clust]*len(dates[1:]),
                                            [direction]*len(dates[1:])))

            total_wave_power = {str((year_1,year_2)):x for x,year_1,year_2 in zip(total_wave_power_list,dates[:-1],dates[1:])}

            with open('processed_waves/{}_{}_{}.json'.format(atoll,clust,direction), 'w') as fp:
                json.dump(total_wave_power, fp)