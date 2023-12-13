#!/usr/bin/python
#################
# Script to plot maps of satellite data.
# Both SST and CHL are supported
#################

import sys, getopt
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import scale
import matplotlib.colors as colors
import numpy as np
from mpl_toolkits.basemap import Basemap
import time
import scipy.io as spio
from netCDF4 import Dataset as NetCDFFile
from plot_lib import *
import cmocean

if __name__ == "__main__":
    # Pick a date:
    title = 'Satellite observations'
    tstart = '20150701' # start day (in format: 'yyyymmdd')
    tend = '20150731'   # end day (in format: 'yyyymmdd')
    ampm   = 'am'	# am or pm (string)
    domain = 'ba'       # Domain to plot: 'no' for both domain or 'ku' for fine only
    plotcb = 1            # Whether to plot the colorbar
    save   = 1          # (1) save file
    varnum = 2 	# Variable number from the var_names routine: 2-SST, 24-CHL
    ssttype='L3S'        # Choose SST observation type: 'L4' or 'L3S'
    chltype='ba'        # Choose CHL observation data set: 'ba' or 'no'

    if varnum==2:
        minmax = [16, 16]        # max/max plotted values - set min=max for automatic
        plotlog=0
        if tstart[4:6]=='02':
          minmax = [-1, 6.0]
        elif tstart[4:6]=='03':
	  minmax = [-1, 6.0]
	elif tstart[4:6]=='04':
	  minmax = [-1, 10]
	elif tstart[4:6] == '05':
	  minmax = [-1, 14]
	elif tstart[4:6] == '07':
	  minmax = [11, 21]
	elif tstart[4:6] == '09':
	  minmax = [12, 18.5]
	elif tstart[4:6] == '11':
	  minmax = [0.5, 13]
    else:
        minmax = [0.1, 10.0]        # max/max plotted values - set min=max for automatic
        plotlog=1
            

    ####################################################
    # Ignore these variables
    var, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()
    ######################################################

    if varnum == 2: #TEMP
      strcmap = 'coolwarm'   # colormap
    elif varnum==10:
        strcmap = 'viridis'
    elif varnum==11:
        strcmap = 'inferno'
    elif varnum==12:
        strcmap = cmocean.cm.deep
    elif varnum==21:     # OXY
        strcmap = cmocean.cm.haline
    elif varnum==23:     # pH
        strcmap = 'gist_ncar_r'
    elif varnum==24:     # Chlorophyll
        strcmap = cmocean.cm.thermal
    elif varnum==25:   # TE
        strcmap = cmocean.cm.dense
    elif varnum==26:   # PFT
        strcmap = cmocean.cm.rain
    else:
        strcmap = 'gist_ncar_r'

    if varnum==2:
        if ssttype=='L3S':
            datatype = 'L3S'  # 'multi', 'SLSTRA', 'METOPB'
        else:
            datatype = 'REP_L4'  # 'multi', 'SLSTRA', 'METOPB'
    else:
        if chltype=='ba':
            datatype = 'ba_MY'  # 'multi', 'SLSTRA', 'METOPB'
        else:
            datatype = 'no_MY'  # 'multi', 'SLSTRA', 'METOPB'
 
    date_range = generate_date_range(tstart, tend)

    data_list = []

    for d in date_range:
      # Read observations
      if varnum==2:
        data, lat, lon = get_cmems_sst(d[0:4], d[4:6], d[6:8], ampm, datatype, ssttype)
#        if ssttype=='L3S':
        data = np.flip(data,0)
        lat = np.flip(lat,0)
      else:
        log = 0
        data, lat, lon = get_cmems_chl(d[0:4], d[4:6], d[6:8], ampm, datatype, log)
        if datatype=='ba_MY':
            data = np.flip(data,0)
            lat = np.flip(lat,0)
      data_list.append(data)

    data_mean = np.ma.mean(data_list, axis=0)
    lons, lats = np.meshgrid(lon.compressed(),lat.compressed())
    out_lat = np.array(lats) > 63
    out_lon = np.array(lons) < 15
    out_bal = out_lat & out_lon
    data_mean[out_bal] = np.nan

    #### Plotting


    # Set file name
    if varnum==2:
        fname= 'Obs_'+str(MAT_VAR)+'-'+ssttype+'_'+tstart+'-'+tend+'_'+domain+'.png'
    else:
        fname= 'Obs_'+str(MAT_VAR)+'-'+chltype+'_'+tstart+'-'+tend+'_'+domain+'.png'

    print('File: '+ fname)
    title = 'Satellite observations'
    plot_map(data_mean, lat, lon, varnum, domain, \
             strcmap, minmax, save, title, fname, plotlog, plotcb)

