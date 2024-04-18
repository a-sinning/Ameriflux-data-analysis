#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Code to calculate stability functions and wind speeds using Monin-Obuhov similarity theory
#Thanks to Dr. Camilo Rey-Sanchez for providing parts of this code and helping with the code

import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
from scipy import stats
import csv


# In[2]:


directory=r'/Users/Amanda1/Downloads'
filename=r'/AMF_US-xGR_BASE-BADM_6-5 (1)/AMF_US-xGR_BASE_HH_6-5.csv'
file=directory+filename

data = np.array(pd.read_csv(file))
data[data==-9999]=np.nan
data2 =pd.read_csv(file, engine='python', skiprows= [i for i in range(1,78817)], skipfooter=14581)
#skipping rows of the .csv file before July 2021 data and after July 2021, adjust based on what month is of interest
#July 2021: i for i in range(1,78817)], skipfooter=14581
#August 2021:    (tropical storm Fred month, but will want to look at 2 weeks before and 2 weeks after Fred)
data2[data2==-9999]=np.nan

## Get date
time_start=data[:,0]

date_c_all=[]
for i in np.arange(1,len(time_start)):
    ts=str(data[i,0])
    year=ts[0:4]
    month=ts[4:6]
    day=ts[6:8]
    if day=='00': 
        continue
    hour=ts[8:10]
    minute=ts[10:12]
    print(year)
    print(month)
    print(day)
    print(ts)
    print(time_start)
    date_c = datetime(int(year),int(month),int(day),int(hour),int(minute))   
    date_c_all.append(date_c)
    #if year != '2021' and month != '07':    #add this to only look at the month of July, 2021 (for milestone 3)
        #continue
    
SH= data2['SH']         #pull SH (W/m^2) from the data
TA= data2['TA_1_1_1'] #Air temperature, pull from the data (TA, in deg C), 
TA_k= TA+273.15 #convert air temp to kelvin
pa= data2['PA']    #pressure in kPa, pull from data (PA, in kPa, convert to hPa below)
#Tbar=      #averaged air temp, convert to Kelvin
timestamp = data2['TIMESTAMP_START']
timestamp=timestamp.to_numpy()
ustar=data2['USTAR']
windspeed=data2['WS_1_1_1']


# In[3]:


#Parameters
h=30; #canopy height (m)
z_tower=45.40; # tower height (m)
z=100 #height want to calculate winds at (m)
d=19.8; # displacement height (m)
zo=3; # roughness length for momentum (m)
zoh=3; #roughness length for heat (m)
g= 9.8 #gravitational constant (m/s)
cp = 1.004 #J/k*g, heat capacity
rho = 1200 #g/m^3, air density (1.2 kg/m^3)

k=0.4;  # von Karman Constant
#ustar=0.36; #friction velocity m/s

p= pa*10 #pressure in hPa
thetabar= TA_k*(1000/1000)**0.286; # mean Potential temperature (K) (calc from Tbar) #assume p=1000 (surface)
thetaw= (SH/cp)/rho ; # k m s-1, potential temp flux, calculate from SH flux in kinematic form
thetastar=-(thetaw/ustar);

L= (-ustar**3)/(k*(g/thetabar)*thetaw); # Obukhov length (m)
zL=(z_tower-d)/L;  # Stability parameter
zL= zL.to_numpy() #convert to numpy array, works better for if statement
L= L.to_numpy()
thetastar=thetastar.to_numpy()
thetabar=thetabar.to_numpy()
ustar=ustar.to_numpy()
#print('Obukhov length', L)
#print('MOST parameter', zL)

type(thetastar)
print(np.shape(TA))
#ind_pos = [82832, 82845, 4090] #find the values of array at different observation times
#print (L[ind_pos])

print(type(ustar))


# In[4]:


#Parameters
bool_arr= zL >0
#if zL>0:      #stable conditions
if np.any(bool_arr):   
    zt=0;
    phim=1+5*zL; # Eq 4.37 Lee. Stability coefficient for wind shear
    phih=1+5*zL; # Eq 4.37 Lee. Stability coefficient for eddy diffusivity of any scalar(k). Assumed to equal to diffusivity for heat (van Ulden,1978)
    psim=-5 * zL;
    psih=psim;
        
else:          #unstable conditions
    zt=(1 - 16 * zL)**-0.25;
    phim=(1-16*zL)**-0.25; # Eq 4.35 Lee. Stability correction function for wind shear
    phih=(1-16*zL)**-0.5; # Eq 4.36 Lee. Stability correction function for heat or any scalar(k). Assumed to equal to diffusivity for heat (van Ulden,1978)
    psim=np.log((1 + zt**2)/2)*(((1 + zt)/2)**2) - (2 * np.arctan(zt)) + (np.pi/2);
    psih=2*np.log((1+zt**2)/2);
        

Km = k * ustar * z/phim;# Eddy diffusivity for momentum
Kh = k * ustar * z/phih;# Eddy diffusicity for heat


print(['Eddy diffuvity for momentum (K_m) equals', str(Km), 'm^2 s^{-1})'])
print(['Eddy diffuvity for heat (K_h) equals', str(Kh), 'm^2 s^{-1})'])

#ubar=ustar/k*((np.log((z-d)/zo))-psim);# eq 4.40 Lee
#print(['Wind speed at height', str(z), ' m, equals ', str(ubar), ' m s^{-1})'])

# Temperature profile
thetadiff=thetastar/k*(np.log((z-d)/zoh)-psih);# eq 4.41 Lee
print(['Temperature difference between the height z and the effective surface (z-d) equals', str(thetadiff), ' K'])
#print(type(thetadiff))
#print(type(thetabar))
thetasurf=thetabar-thetadiff;
print(['Temperature at the effective surface (z-d) equals', str(thetasurf), 'K'])


# In[5]:


#print(np.nanmax(ubar)) #find the time with maximum wind speed, excluding nans
#print(np.nanargmax(ubar)) #find the index at which ubar is maximum, so can find the value of other variables at that index
print(np.nanmax(windspeed), 'max windspeed in m/s')
print(np.nanargmax(windspeed),'time index where windspeed is max, WS measured at z=45.89m at site')
#max_wind= time
#print(type(ubar))
#print(ubar)
halfhour_maxwind= timestamp[940]
print(halfhour_maxwind) #half hour with the max wind speed, index 940 for July 2021
print(ustar[940], 'ustar value at time index of maxwind speed')
print(p[940])
print(psim[940],'psim value at time index of maxwind speed')


# In[13]:


#loop to calculate the wind profile for the time step of highest wind speed  (index 940 for July 2021)
#just for my specific half hour timestep index [940], calculate wind from z=1m to z=100m to make profile
n_z=100 #create array of heights/vertical levels (from 1 to 100)
u=np.zeros(n_z)
z_levels=np.arange(1,101)
for i_z in range (22,n_z):     #loop calculates wind speed (in m/s) for each height between 22m and 100m (nan or negative below 19m)
    u[i_z] = 0.48/k*((np.log((z_levels[i_z]-d)/zo))-(-0.34)) #ubar loop calculation
    #print(u[i_z])

#make plot of this wind profile
plt.figure(figsize=(10,8))
plt.plot(u,z_levels,'.')
plt.xlabel('windspeed (m/s)')
plt.ylabel('z (m)')
plt.ylim(20,105) #shift the wind profile up the magnitude of the displacement height (19.8) so that not plotting winds below there since not accurate)
plt.title('Wind Profile US-xGR: NEON Great Smoky Mtns. National Park, July ..... 2021') #this timestep is July 19, 2021 at noon
plt.grid()
#plt.legend()


# In[12]:


#wind profile using power law, for same time index
#WS measured at z=45.89m at site, WS at this height is 3.89 m/s (use as reference height)
n_z=100
z_levels=np.arange(1,101)
z_r=45.89 #reference height (m)
u_r=3.89 #wind speed at reference height (m/s)
u_pl=np.zeros(n_z) #power law wind profile
m=0.28 #use the power law exponent for the "forest" surface type
#use for loop to calculate wind speed at each height from 0-100m
for i_z in range (1,n_z): 
    u_pl[i_z]=u_r*((z_levels[i_z])/(z_r-d))**m #needed to adjust z_r to z_r -d to incorporate the displacement height, since wiond profile looks funny at lower heights if not incorporate displacement height
    #print(u_pl[i_z])
    
plt.figure(figsize=(10,8))
plt.plot(u_pl,z_levels,'.')
plt.ylim(20,105) #shift the wind profile up the magnitude of the displacement height (19.8) so that not plotting winds below there since not accurate)
plt.xlabel('windspeed (m/s)')
plt.ylabel('z (m)')
plt.title('Wind Profile US-xGR: NEON Great Smoky Mtns. National Park, July 19, 2021 12Z') #figure out which day/hr this is 
plt.grid()


# In[15]:


#comparison plot between MOST profile and power law profile (for the specific half hour of high wind speed)
plt.figure(figsize=(10,8))
plt.plot(u_pl,z_levels,label='Power Law',color='skyblue')
plt.plot(u,z_levels,label='MOST', color='forestgreen')
plt.plot(3.89,45.89,'x',label='Obs wind speed',color='black')
plt.xlabel('windspeed (m/s)')
plt.ylabel('z (m)')
plt.ylim(20,105) #shift the wind profile up the magnitude of the displacement height (19.8) so that not plotting winds below there since not accurate)
plt.title('Wind Profile Comparison for US-xGR: NEON Great Smoky Mtns. National Park, July 19, 2021 12Z') #figure out which day/hr this is 
plt.grid()
plt.legend()


# In[24]:


#HISTOGRAMS FOR STABLE AND UNSTABLE WIND SPEEDS AT 100m

stable=[]
unstable=[]

for i in range (1,len(ustar)):
    if zL[i]>0:      #stable conditions   
        zt=0;
        phim=1+5*zL[i]; # Eq 4.37 Lee. Stability coefficient for wind shear
        phih=1+5*zL[i]; # Eq 4.37 Lee. Stability coefficient for eddy diffusivity of any scalar(k). Assumed to equal to diffusivity for heat (van Ulden,1978)
        psim=-5 * zL[i];
        psih=psim;
        ubar= ustar[i]/k*(np.log((z-d)/zo)-psim)
        stable.append(ubar)
        
    else:          #unstable conditions
        zt=(1 - 16 * zL[i])**-0.25;
        phim=(1-16*zL[i])**-0.25; # Eq 4.35 Lee. Stability correction function for wind shear
        phih=(1-16*zL[i])**-0.5; # Eq 4.36 Lee. Stability correction function for heat or any scalar(k). Assumed to equal to diffusivity for heat (van Ulden,1978)
        psim=np.log((1 + zt**2)/2)*(((1 + zt)/2)**2) - (2 * np.arctan(zt)) + (np.pi/2);
        psih=2*np.log((1+zt**2)/2);
        ubar= ustar[i]/k*(np.log((z-d)/zo)-psim);
        unstable.append(ubar)

        
Km = k * ustar[i] * z/phim;# Eddy diffusivity for momentum
Kh = k * ustar[i] * z/phih;# Eddy diffusicity for heat

#wind profile
ubar=ustar[i]/k*((np.log((z-d)/zo))-(psim))
ubar=UBAR[i]

#temp profile
thetadiff=thetastar[i]/k*(np.log((z-d)/zoh)-psih)
theta100=thetabar[i]-thetadiff;
theta100=THETA[i]

#plt.subplot(1,2,1)
plt.figure(figsize=(10,8))
plt.hist(stable,bins=np.arange(0,18,1))
plt.title('Wind Speed Observations at 100m for July 2021 - STABLE Conditions')
plt.xlabel('Wind speed (m/s)')
plt.ylabel('# of Observations')

#plt.subplot(1,2,2)
plt.figure(figsize=(10,8))
plt.hist(unstable,bins=np.arange(0,18,1), color='indianred')
plt.title('Wind Speed Observations at 100m for July 2021 - UNSTABLE Conditions')
plt.xlabel('Wind speed (m/s)')
plt.ylabel('# of Observations')


# In[ ]:




