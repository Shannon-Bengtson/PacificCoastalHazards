import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4
import requests
from datetime import date, timedelta
import xarray as xr
import os
from threading import Thread
from time import sleep, perf_counter

start_date = date(1990,1,15) 
end_date = date(2021,1,15)
delta = end_date-start_date

def task(day):
    print(f'Starting the task {day}...')
    
    file = f"https://www.star.nesdis.noaa.gov/pub/sod/mecb/crw/data/5km/v3.1_op/nc/v1.0/daily/ssta/{day.strftime('%Y')}/ct5km_ssta_v3.1_{day.strftime('%Y')}{day.strftime('%m')}{day.strftime('%d')}.nc"

    os.system(f'wget {file}')

    print('done')
    
    
start_time = perf_counter()

######

# create and start threads
threads = []
for i in range(delta.days + 1):
    
    day = start_date + timedelta(days=i)
    
    t = Thread(target=task, args=(day,))
    
    sleep(2)
    threads.append(t)
    t.start()

# wait for the threads to complete
for t in threads:
    t.join()

end_time = perf_counter()
