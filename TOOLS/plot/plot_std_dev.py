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
    depth  = '0'	# Depth
    year   = 2015       # Year
    month  = '03'       # Month has to be string and two digits i.e '05' and '10'
    day    = '03'	# Day has to be double digits
    ampm   = '00'	# am or pm (string)
    dist = 10           # Grid point distance to check around central location 
    station = 'ArkonaWR' # Set station name according to the XLSX file
    istation = 0        # (0) use 'station', (1) Arkona, (2) Baltic Proper, (3) Bothnian Sea

    ref_free = 1
    rel_conc = 1
    minmax = [0,0]
    save = 1

    assim1 = 'chl_weakly_sstL3_strongly_N30'
    #assim1 = 'free_N30'

    tstart = '20150201'
    tend = '20150331'
    
    z_mean=0     #0(no vertical mean)
    z_integral=0 #(0 for no vertical integral; 1 for vertical integral)
    z1=0         # upper Z boundary (for z_mean and z_integral)
    z2=5         # lower Z boundary (for z_mean and z_integral)

    varnum = 7	 	# Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP


    ######################################################

    # minmax 
    if ref_free == 0: 
      if varnum == 2: 
        minmax = [0, 0.5]
    if ref_free == 1: 
      if varnum == 2: 
        minmax = [-0.125,0.125] 
      elif varnum == 24: 
        minmax = [-0.7, 0.7]
      elif varnum == 10: 
        minmax = [-0.2, 0.2]
      elif varnum == 11: 
        minmax = [-0.25, 0.25]
      elif varnum == 12: 
        minmax = [-0.08, 0.08]
      elif varnum == 6: 
        minmax = [-0.03, 0.03]
      elif varnum == 7: 
        minmax = [-300, 300]
        minmax = [0,0]
      elif varnum == 8:
        minmax = [-500, 500] 
      elif varnum == 23: 
        minmax = [-0.04, 0.04] 
      elif varnum == 27: 
        minmax = [-50, 50]
      else: 
        minmax = [1,1]

    date_range = generate_date_range(tstart, tend)
    varstr, mat_var, Variable, var_unit = var_names(varnum)# loop over each time step 
    DA_type = 0
    loopcnt = 0
    data_list = np.array([])
    free_list = np.array([])
    conc_list = np.array([])

    #strcmap = 'RdGy'
    strcmap = 'plasma'

    for d in date_range:
      print d
      path = '/scratch/projects/hbk00095/exp/' + assim1 + '/DA/'
      # read variance  
      data_coarse, coarse_keys = read_model_var_path(varstr, 'coarse', d[0:4], d[4:6], d[6:8], ampm, depth, DA_type, coupled, z_mean, z_integral, z1, z2, path)
      data_coarse = np.flip(np.array(data_coarse), 0)
      if rel_conc == 1:
        conc_coarse, _ = read_model_path(varstr, 'coarse', d[0:4], d[4:6], d[6:8], ampm, depth, DA_type, coupled, z_mean, z_integral, z1, z2, path)
        conc_coarse = np.flip(np.array(data_coarse), 0)
      if loopcnt == 0: 
        fill_value = data_coarse[0,0]
        mask = (data_coarse == fill_value)
      data_coarse = np.where(mask, np.nan, data_coarse) 
      if rel_conc == 1: 
        conc_coarse = np.where(mask, np.nan, data_coarse)
      # from variance to standard deviation
      data_coarse = np.sqrt(data_coarse)
      if ref_free == 1: 
        path_free = '/scratch/projects/hbk00095/exp/free_N30/DA/'
        free_coarse, _ = read_model_var_path(varstr, 'coarse', d[0:4], d[4:6], d[6:8], ampm, depth, DA_type, coupled, z_mean, z_integral, z1, z2, path_free)
        free_coarse = np.flip(np.array(free_coarse), 0)
        free_coarse = np.where(mask, np.nan, free_coarse) 
        # from variance to standard deviation
        free_coarse = np.sqrt(free_coarse)
        # calc difference to free
        data_coarse = data_coarse - free_coarse
        if rel_conc == 1: 
          data_coarse = np.divide(data_coarse, conc_coarse)
      if loopcnt > 1: 
        data_coarse = data_coarse.reshape(data_coarse.shape[0], data_coarse.shape[1], 1)
        data_list = np.concatenate((data_list, data_coarse), axis=2)
      elif loopcnt == 1: 
        data_list = np.stack((data_list, data_coarse), axis=2)
      else: 
        data_list = data_coarse
      if loopcnt == 0:
        # set coordinates
        lat = coarse_keys['lat'][:]
        lon = coarse_keys['lon'][:]

        lat1 = lat[:,0]
        lon1 = lon[0,:]
        for i in range(len(lat[:,0])):
          lat1[i] = np.max(lat[i,:])
        for j in range(len(lat[0,:])):
          lon1[j] = np.max(lon[:,j])
        lon1[-1] = lon1[-2] + lon1[1]-lon1[0]
        # Flip arrays for plotting
        lat1 = np.flip(lat1)
      loopcnt += 1

    # calculate mean per time step over grid 
    timeseries = np.nanmean(data_list, axis=0)
    print timeseries.shape
    timeseries = np.nanmean(timeseries, axis=0)
    print 'TEST', timeseries.shape

    ### Calculate mean per grid cell over timespan
    std_dev = np.nanmean(data_list, axis=2)
    print('Mean = ', np.nanmean(std_dev))
    
    #### Plotting
    if ref_free == 0:
      fname= str(mat_var)+'_stddev_field_'+assim1+'.png'
    elif ref_free == 1:
      fname= str(mat_var)+'_stddev_field_'+assim1+'-free.png'
    print ('File '+fname)

    domain = 'ba'
    plotlog = 0
    plotcb = 1
    if ref_free == 0:
      title = assim1 
    elif ref_free == 1: 
      title = assim1 + ' difference to free' 
    
    # plot map 
    plot_map(std_dev, lat1, lon1, varnum, domain, \
	strcmap, minmax, save, title, fname, plotlog, plotcb)

    # plot timeseries 
    # Determine tick localizations and labels
    if ref_free == 1 and rel_conc == 1: 
      fname= str(mat_var)+'_stddev_field_timeseries_'+assim1+'-free-div_conc.png'

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

    fig, ax = plt.subplots(figsize=(8,3))
    plt.plot_date(freeday, timeseries, '-')
    plt.axhline(0, color='k', linestyle='dashed', linewidth=1.5)
    ax.xaxis.set_major_locator(ticker.FixedLocator(tickloc))
    ax.xaxis.set_minor_locator(ticker.FixedLocator(tickloc_min))
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.tick_params(axis='x', which='minor', tick1On=False, tick2On=False)
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(months_str))
    plt.legend()
    plt.title(Variable)
    plt.ylabel(var_unit)
    if save == 1: 
      plt.savefig(fname, dpi=300)
