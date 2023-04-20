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
    year   = 2015       # for year 2012 we plot NOAA SST data, for >=2017 Sentinel 3a
    month  = '03'       # Month has to be string and two digits i.e '05' and '10'
    day    = 14	        # Day 
    ampm   = 'am'	# am or pm (string)
    domain = 'no'       # Domain to plot: 'no' for both domain or 'ku' for fine only
    plotcb=1            # Whether to plot the colorbar
    save   = 0          # (1) save file
    varnum = 24 	# Variable number from the var_names routine: 2-SST, 24-CHL
    ssttype='L4'        # Choose SST observation type: 'L4' or 'L3S'
    chltype='ba'        # Choose CHL observation data set: 'ba' or 'no'

    if varnum==2:
        minmax = [0, 10]        # max/max plotted values - set min=max for automatic
        plotlog=0
    else:
        minmax = [0.1, 10.0]        # max/max plotted values - set min=max for automatic
        plotlog=1
            

    ####################################################
    # Ignore these variables
    var, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()
    ######################################################

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
 

    # Read observations
    if varnum==2:
        data, lat, lon = get_cmems_sst(year, month, day, ampm, datatype, ssttype)
#        if ssttype=='L3S':
        data = np.flip(data,0)
        lat = np.flip(lat,0)
    else:
        log = 0
        data, lat, lon = get_cmems_chl(year, month, day, ampm, datatype, log)
        if datatype=='ba_MY':
            data = np.flip(data,0)
            lat = np.flip(lat,0)


    #### Plotting

    # Set plot title
    if day<10:
        nullstr = "0"
    else:
        nullstr = ""
    title='Satellite '+str(Variable)+': '+str(year)+'-'+str(month)+'-'+nullstr+str(day)

    # Set file name
    if varnum==2:
        fname= 'Obs_'+str(MAT_VAR)+'-'+ssttype+'_'+str(year)+str(month)+nullstr+str(day)+'_'+domain+'.png'
    else:
        fname= 'Obs_'+str(MAT_VAR)+'-'+chltype+'_'+str(year)+str(month)+nullstr+str(day)+'_'+domain+'.png'

    print('File: '+ fname)

    strcmap = 'coolwarm'   # colormap
    strcmap = 'jet'   # colormap
    if varnum==2:
        strcmap = cmocean.cm.thermal
        strcmap = 'coolwarm'   # colormap
    else:
        strcmap = cmocean.cm.algae
        strcmap = 'viridis'
#        strcmap = cmocean.cm.thermal
#    strcmap = 'gist_earth'

    plot_map(data, lat, lon, varnum, domain, \
             strcmap, minmax, save, title, fname, plotlog, plotcb)

