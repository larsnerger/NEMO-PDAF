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
    month  = '04'       # Month has to be string and two digits i.e '05' and '10'
    day    = 5	        # Day 
    ampm   = 'am'	# am or pm (string)
    domain = 'no'       # Domain to plot: 'no' for both domain or 'ku' for fine only
    plotcb=1            # Whether to plot the colorbar
    save   = 0          # (1) save file
    varnum = 2 	# Variable number from the var_names routine: 2-SST, 24-CHL
    ssttype = 'L3S'

    if varnum==2:
        minmax = [3, 8]        # max/max plotted values - set min=max for automatic
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
       datatype = 'ba_MY'  # 'multi', 'SLSTRA', 'METOPB'
 

    # Read observations
    if varnum==2:
        data, lat, lon = get_cmems_sst(year, month, day, ampm, datatype, ssttype)
#        if ssttype=='L3S':
        data = np.flip(data,0)
        lat = np.flip(lat,0)
    else:
        log = 0
        data, lat, lon = get_cmems_chl(year, month, day, ampm, datatype, log)

    # To determine indices for some coordinate
#    for i in range(len(lat)):
#        if lat[i]>54.8 and lat[i]<55.0:
#            print i, lat[i]
#    for i in range(len(lon)):
#        if lon[i]>13.8 and lon[i]<14.0:
#            print i, lon[i]

#    cnt = 0
#    for j in range(250,271):
#        for i in range(138,159):
#         if data.mask[i,j]==False:
#            cnt = cnt + 1
#    print "Number of observations around Arkona: ", cnt

    #### Plotting

    # Set plot title
    if day<10:
        nullstr = "0"
    else:
        nullstr = ""
    title='Satellite '+str(Variable)+': '+str(year)+'-'+str(month)+'-'+nullstr+str(day)

    # Set file name
    fname= 'Obs_'+str(MAT_VAR)+'_'+str(year)+str(month)+nullstr+str(day)+'_'+domain+'.png'
    print('File: '+ fname)

    strcmap = 'coolwarm'   # colormap
    strcmap = 'jet'   # colormap
    strcmap = cmocean.cm.thermal
#    strcmap = cmocean.cm.algae
#    strcmap = 'gist_earth'

    plot_map(data, lat, lon, varnum, domain, \
             strcmap, minmax, save, title, fname, plotlog, plotcb)

