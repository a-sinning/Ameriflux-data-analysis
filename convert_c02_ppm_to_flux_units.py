#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#convert c02 to flux units from Ameriflux data


# In[4]:


# Import Libraries
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
from scipy import stats
import csv
import math


# In[5]:


directory=r'/Users/Amanda1/Downloads'
filename=r'/Field Trip Data (1).csv'
file=directory+filename

#print(type(pd.read_csv(file)))
data = np.array(pd.read_csv(file))
data[data==-9999]=np.nan
data2 =pd.read_csv(file)
data2[data2==-9999]=np.nan


# In[13]:


#timeseries that I need to make: GPP, respiration
Time = data2['TIME']
#GPP=carbon flux + respiration
ppm= data2['CO2'] #overall carbon flux (ppm) (need to transform this into a flux)
#convert C02 to flux in micromol/m^2 s .. use IGL and conversion factors for time and surface area
#solve for n/v in n/v=Rd/T, which is mol/m^3
#define a function to convert ppm to desired flux units, that takes the 3 variables needed for conversion
def ppm_to_umol_per_m2_per_s(ppm, pressure, temperature):
    R = 8.314 #ideal gas constant
    pressure_pa = pressure * 1000 # Convert pressure to Pascals
    temperature_k = temperature + 273.15  # Convert temperature to Kelvin
    
    # Calculate the conversion factor
    cf = (10**6 / 22.4) * (pressure_pa / R / temperature_k) * (1 / 3600) * (1 / 10**6)
    
    # Convert CO2 concentration from ppm to micromol/m^3
    conc_umol_per_m3 = ppm * 10**-6 * (pressure / 101.3) * (temperature + 273.15) / 298.15
    
    # Convert micromol/m^3 to micromol/m^2/s
    conc_umol_per_m2_per_s = conc_umol_per_m3 * cf
    
    return conc_umol_per_m2_per_s #carbon flux returned in proper flux units


# In[ ]:


# Read in the orig CSV with CO2 concentrations in ppm
df = pd.read_csv(file, usecols=['CO2'])

# Set the atmospheric pressure and temperature
pressure = 101.0  # in kPa ...was a fair weather day, will say was around 1010 hPa
temperature = 16  # temp was around 16 deg C the morning we took the measurements

# Convert CO2 concentrations from ppm to micromol/m^2/s by applying function above
df['umol_per_m2_per_s'] = df['CO2'].apply(lambda x: ppm_to_umol_per_m2_per_s(x, pressure, temperature))

# Print the resulting Dataframe and save new data
print(df)
df.to_csv('c02_fluxunits.csv', index=False)

#respiration is the times when we were using the dark chamber, so avg just those values for calculating the respiration

