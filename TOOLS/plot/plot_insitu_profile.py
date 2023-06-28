#!/usr/bin/python
#################
# Python script used for plotting maps
# at specific times/dates on specific
# vaiables.
#################

import sys, getopt
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
import numpy as np
from mpl_toolkits.basemap import Basemap
import time
import scipy.io as spio
#from colormaps import cmaps
from netCDF4 import Dataset as NetCDFFile
from plot_lib import *
import cmocean


if __name__ == "__main__":
    coupled = ''
    figure_name = 'Profile_sst+chl_temp_Arkona.png'
    year_p = 2015
    month_p = [2]
    lat_stat = 54.8830
    lon_stat = 13.8660
    dist = 10           # Grid point distance to check around central location 
    plot_free = 1 
    plot_fcst = 1 
    plot_ana = 1   
    #path_data = '/scratch/usr/hbxsovli/SEAMLESS/validation_data/ICES_2015_Bottle_LR_CTD/ICES_2015_Bottle_LR_CTD.xlsx' # path to validation data .xlsx file 
    path_data = '/scratch/usr/hbxsovli/SEAMLESS/validation_data/ICES_2015_HR_CTD/ICES_2015_HR_CTD.xlsx'
    z_mean=0     #0(no vertical mean)
    z_integral=0 #(0 for no vertical integral; 1 for vertical integral)
    z1=0         # upper Z boundary (for z_mean and z_integral)
    z2=5         # lower Z boundary (for z_mean and z_integral)
    save = 0

    domain = 'no'
    varnum = 2	 	# Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP
    plot_save = 1
    ######################################################

    varstr, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()
    
    varnum_data2 = [2, 3, 21] # variables contained in HR_CTD data set
    hour_p = [0, 24] # valid hours 

    data_temp = []
    depth_temp = []
    minute_temp = []
    hour_temp = []
    day_temp = []
    month_temp = []
    lat_temp = []
    lon_temp = []

    # read data 
    year, month, day, hour, minute, lon, lat, depth, data, data_qv = read_insitu_vali_ICES(varnum, path_data)
     
    print 'Data source 1 total length = ', len(data) 
    if len(year)>0:
      for r in range(len(data)):
        if type(data[r]) != str and year[r]==year_p and (month[r] in month_p) and (hour[r] in hour_p):
          data_temp.append(data[r])
          lat_temp.append(lat[r])
          lon_temp.append(lon[r])
          depth_temp.append(depth[r])
          day_temp.append(day[r])
          month_temp.append(month[r])
   
    # find closest lat/lon 
    latdiff = np.array(lat_temp) - lat_stat
    londiff = np.array(lon_temp) - lon_stat
    lat_idx = np.argmin(abs(latdiff))
    lon_idx = np.argmin(abs(londiff))
 
    data_final = []
    depth_final = []
    minute_final = []
    hour_final = []
    day_final = []
    month_final = []
    lat_final = []
    lon_final = []

    if lat_idx != lon_idx: 
      print 'ERROR: lat and lon idx is not the same, this scenario still needs to be programmed'
    else: 
      lat_min = lat_temp[lat_idx]
      lon_min = lon_temp[lon_idx]
      print 'Closest model lat/lon = ', lat_min, lon_min
      print 'At lat/lon index = ', lat_idx, lon_idx
      for i in range(len(data_temp)):
        if abs(lat_temp[i]-lat_min) < 0.1 and abs(lon_temp[i]-lon_min) < 0.1: 
          lat_final.append(lat_temp[i])
          lon_final.append(lon_temp[i])
	  data_final.append(data_temp[i])
	  day_final.append(day_temp[i])
	  month_final.append(month_temp[i])
          depth_final.append(depth_temp[i])

    print 'Length of data_final = ', len(data_final) 

    grid_area = 'coarse'
    time_stamp = '00'
    coupled = 'weak'
    depth = ':'

    # read model data
    # free run
    if plot_free == 1:  
      DA_switch = 0
      free_temp, free_depth_temp = read_latlon_day(varnum, grid_area, year_p, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2, month_final[0], \
                        day_final[0],lat_final[0], lon_final[0], dist) # function
      free_data = []
      free_depth = []

      for i in range(len(free_temp)):
        if free_temp[i]!=0:
          free_data.append(free_temp[i])
          free_depth.append(free_depth_temp[i])

    # Forecast
    if plot_fcst==1: 
      DA_switch = 1
      fcst_temp, fcst_depth_temp = read_latlon_day(varnum, grid_area, year_p, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2, month_final[0], \
                        day_final[0],lat_final[0], lon_final[0], dist) # function
      fcst_data = []
      fcst_depth = []

      for i in range(len(fcst_temp)):
        if fcst_temp[i]!=0:
          fcst_data.append(fcst_temp[i])
          fcst_depth.append(fcst_depth_temp[i])

    # Analysis
    if plot_ana==1: 
      DA_switch = 2
      ana_temp, ana_depth_temp = read_latlon_day(varnum, grid_area, year_p, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2, month_final[0], \
                        day_final[0],lat_final[0], lon_final[0], dist) # function
      ana_data = []
      ana_depth = []

      for i in range(len(ana_temp)):
        if ana_temp[i]!=0:
          ana_data.append(ana_temp[i])
          ana_depth.append(ana_depth_temp[i])

    
    # plot 
    fig, ax = plt.subplots(figsize=(8,6))
    plt.plot(data_final, depth_final, c='r', label='Insitu')
    if plot_free == 1: 
      plt.plot(free_data, free_depth, c='b', label='Free')
    if plot_fcst == 1:
      plt.plot(fcst_data, fcst_depth, c='#bbbb00', label='Forecast')
    if plot_ana == 1: 
      plt.plot(ana_data, ana_depth, c='k', label='Analysis')
       
    plt.xlabel(var_unit)
    plt.ylabel('Depth [m]')
    plt.legend()
    ax.invert_yaxis()
    plt.title(Variable + ' Profile at location lat/lon = ' + str(round(lat_final[0],2)) +'/'+ str(round(lon_final[0],2)) + ' on ' + \
		str(day_final[0]) + '-' + str(month_final[0])+'-'+str(year_p))
    if plot_save == 1: 
      plt.savefig(figure_name)
    plt.show()
