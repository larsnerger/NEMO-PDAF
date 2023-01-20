#!/usr/bin/python
#################
# Script to plot maps of satellite data.
# Both SST and CHL are supported
#################

import sys, getopt
import matplotlib.pyplot as plt
from matplotlib import ticker
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
    day    = 10	        # Day 
    forana = 'a'
    member = 0          # Ensemble member to plot (0=ensemble mean)
    ampm = 'am'
    domain = 'ba'       # Domain to plot: 'no' for both domain or 'ku' for fine only
    minmax = [16, 16]        # max/max plotted values - set min=max for automatic
    minmax = [0, 2]        # max/max plotted values - set min=max for automatic
    save   = 0          # (1) save file
    varnum = 24 	# Variable number from the var_names routine: 2-SST, 24-CHL

    ####################################################
    # Ignore these variables
    var, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()
    ######################################################

    datatype = 'lestkf_N30'

    if member==0:
        path_exp = '/scratch/usr/hzfblner/SEAMLESS/run/ens_ERGOM_DA/001/output/DA_mobs'
    else:
        path_exp = '/home/hzfblner/SEAMLESS/rmse/chl_mobs_lestkf_201504_01-15'
 

    # Read model on observation grid
    log = 0
    data_c, lat_c, lon_c = get_mobs_cmems_chl(year, month, day, forana, path_exp, datatype, log, member)

    #### Plotting

    # Flip field and latitude for plotting
    data_c = np.flip(data_c,0)
    lat_c = np.flip(lat_c)

    # Set plot title
    if day<10:
        nullstr = "0"
    else:
        nullstr = ""
    title='CMEMS '+datatype+' '+str(MAT_VAR)+' : '+str(year)+'-'+str(month)+'-'+nullstr+str(day)

    # Set file name
    fname= 'Obs_'+str(MAT_VAR)+'_plot_'+str(year)+'-'+str(month)+str(day)+'_'+domain+'.png'
    print 'File: ', fname

    strcmap = 'coolwarm'   # colormap
    strcmap = 'jet'   # colormap
    strcmap = cmocean.cm.thermal
    strcmap = cmocean.cm.algae
#    strcmap = 'gist_earth'

    plotlog = 0
    plotcb=1
    plot_map(data_c, lat_c, lon_c, varnum, domain, \
             strcmap, minmax, save, title, fname, plotlog, plotcb)
