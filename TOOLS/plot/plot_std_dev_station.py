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
import matplotlib.dates as mdates
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
    depth  = '0'	# Depth (Arkona bottom = 21)
    year   = 2015       # Year
    month  = '03'       # Month has to be string and two digits i.e '05' and '10'
    day    = '03'	# Day has to be double digits
    ampm   = '00'	# am or pm (string)
    dist = 10           # Grid point distance to check around central location 
    station = 'ArkonaWR' # Set station name according to the XLSX file
    istation = 0        # (0) use 'station', (1) Arkona, (2) Baltic Proper, (3) Bothnian Sea

    ref_free = 1        # Difference to free std dev
    rel_conc = 1        # Divided by concentration
    minmax = [0,0]
    save = 1

    assim1 = 'chl_weakly_sstL3_strongly_N30'
    #assim1 = 'free_N30'

    tstart = '20150201'
    tend = '20151231'
    
    z_mean=0     #0(no vertical mean)
    z_integral=0 #(0 for no vertical integral; 1 for vertical integral)
    z1=0         # upper Z boundary (for z_mean and z_integral)
    z2=5         # lower Z boundary (for z_mean and z_integral)

    varnum_range = [27] #range(1,25)	
			# Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP


    ######################################################
    varstr_list = []
    mean_list = []

    for varnum in varnum_range: 
      if varnum == 27: 
        z_integral = 1
      else:
        z_integral = 0
      # minmax 
      if ref_free == 0: 
        if varnum == 2: 
          minmax = [0, 0.5]
      if ref_free == 1: 
        if varnum == 2: 
          minmax = [0,0] 
      date_range = generate_date_range(tstart, tend)
      varstr, mat_var, Variable, var_unit = var_names(varnum)# loop over each time step 
      varstr_list.append(varstr)
      DA_type = 0
      loopcnt = 0
      data_list = []
      conc_list = []
      assim_list = []
      free_list = []

      months = [1]
      station, station_idx, station_coords, noba =  read_station_idx_obs(varnum, year, months, istation, station, dist)
      model_idx = get_model_idx(station_coords[0], station_coords[1], dist)

      # Determine tick localizations and labels
      freeday = [datetime.strptime(date_str, '%Y%m%d').timetuple().tm_yday for date_str in date_range]
      months_free = range(2,13)
      tickloc = []
      tickloc_min = []
      loc = freeday[0]
      tickloc.append(loc)
      months_str = []
      for i in range(len(months_free)):
        _, mstr, _, days_in_mon = month_names(i+months_free[0])
        loc = loc + days_in_mon
        tickloc.append(loc)
        tickloc_min.append(tickloc[i]+days_in_mon/2)
        months_str.append(mstr)
 
      for d in date_range:
        print d
        path = '/scratch/projects/hbk00095/exp/' + assim1 + '/DA/'
        # read variance  
        assim1_data = read_latlon_day_path(varnum, 'coarse', d[0:4], '00', depth, DA_type, coupled, z_mean, z_integral, z1, z2, d[4:6], \
                      d[6:8], model_idx, dist, path, 1) # function
        if rel_conc:
          conc_data = read_latlon_day_path(varnum, 'coarse', d[0:4], '00', depth, DA_type, coupled, z_mean, z_integral, z1, z2, d[4:6], \
                      d[6:8], model_idx, dist, path, 0) # function
          conc_list.append(conc_data)
        assim1_data = np.sqrt(assim1_data)
        if ref_free == 1: 
          path_free = '/scratch/projects/hbk00095/exp/free_N30/DA/'
          free_data = read_latlon_day_path(varnum, 'coarse', d[0:4], '00', depth, DA_type, coupled, z_mean, z_integral, z1, z2, d[4:6], \
                      d[6:8], model_idx, dist, path_free, 1) # function
          free_data = np.sqrt(free_data)
          print assim1_data, free_data
          assim_list.append(assim1_data)
          free_list.append(free_data)
  	  assim1_data = assim1_data - free_data
          if rel_conc: 
	    assim1_data = assim1_data/np.abs(conc_data)
        data_list.append(assim1_data)

      mean_list.append(np.nanmean(data_list))

      #### Plotting
      if ref_free == 0:
        fname= str(varstr)+'_stddev_'+station+'_'+assim1+'.png'
      elif ref_free == 1:
        if rel_conc == 0: 
          fname= str(varstr)+'_stddev_'+station+'_'+assim1+'-free.png'
	elif rel_conc == 1:
          if depth != '0':
            fname= str(varstr)+'_stddev_'+station+'_'+depth+'_'+assim1+'-free-div-conc.png'
          else:
            fname= str(varstr)+'_stddev_'+station+'_'+assim1+'-free-div-conc.png'
      print ('File '+fname)

      fig, ax = plt.subplots(figsize=(8,3))
      plt.plot_date(freeday, data_list, '-', label='CHLw + SSTs') 
      plt.axhline(0, color='k', linestyle='dashed', linewidth=1.5)
      ax.xaxis.set_major_locator(ticker.FixedLocator(tickloc))
      ax.xaxis.set_minor_locator(ticker.FixedLocator(tickloc_min))
      ax.xaxis.set_major_formatter(ticker.NullFormatter())
      ax.tick_params(axis='x', which='minor', tick1On=False, tick2On=False)
      ax.xaxis.set_minor_formatter(ticker.FixedFormatter(months_str))
      plt.legend()
      if station == 'ArkonaWR' and depth == '21':
        plt.title('Bottom '+Variable)
      else:
        plt.title(Variable)
      plt.ylabel(var_unit)

      if save == 1:
        plt.savefig(fname, dpi=300)

    print 'assim std', assim_list
    print 'free std', free_list 
    print 'concentration', conc_list

    'MEAN STD DEV'
    for t in range(0,len(mean_list)):
      print varstr_list[t], mean_list[t]
