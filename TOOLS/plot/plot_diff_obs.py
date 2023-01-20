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
    month  = '03'       # Month has to be string and two digits i.e '05' and '10'
    day    = 1	        # Day 
    ampm   = 'am'	# am or pm (string)
    forana = 'f'
    domain = 'no'       # Domain to plot: 'no' for both domain or 'ku' for fine only
    minmax = [16, 16]        # max/max plotted values - set min=max for automatic
   # minmax = [-1, 1]        # max/max plotted values - set min=max for automatic
    save   = 0          # (1) save file
    varnum = 24 		# Variable number from the var_names routine: 2-SST, 24-CHL
    log = 0

    ####################################################
    # Ignore these variables
    var, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()
    ######################################################

    if varnum==2:
       datatype = 'REP_L4'  # 'multi', 'SLSTRA', 'METOPB'
    else:
       datatype = 'hnk_ska'
       datatype = 'lestkf_N30'

    path_exp = '/scratch/projects/hbk00095/exp/LESTKF/DA_mobs'
 

    # Read SST observations
    if varnum==2:
        data_c, lat_c, lon_c = get_cmems_sst(year, month, day, ampm, datatype)
    else:
#        log = 0
        data_m, lat_c, lon_c = get_mobs_cmems_chl(year, month, day, forana, path_exp, datatype, log, 0)

    if varnum==2:
       datatype = 'REP_L4'  # 'multi', 'SLSTRA', 'METOPB'
    else:
       datatype = 'ba_MY'  # 'multi', 'SLSTRA', 'METOPB'


    # Read SST observations
    if varnum==2:
        data_c, lat_c, lon_c = get_cmems_sst(year, month, day, ampm, datatype)
    else:
#        log = 0
        data_o, lat_o, lon_o = get_cmems_chl(year, month, day, ampm, datatype, log)


    data_c = data_m - data_o

    meanerr = np.ma.mean(data_c)
    print 'mean error', meanerr

    #### Plotting

    # Flip field and latitude for plotting
    data_c = np.flip(data_c,0)
    lat_c = np.flip(lat_c)

    # Set plot title
    if day<10:
        nullstr = "0"
    else:
        nullstr = ""
    title='Model-obs '+datatype+' '+str(MAT_VAR)+' : '+str(year)+'-'+str(month)+'-'+nullstr+str(day)

    # Set file name
    fname= 'figs/Obs_'+str(MAT_VAR)+'_plot_'+str(year)+'-'+str(month)+str(day)+'_'+domain+'.png'

    strcmap = 'coolwarm'   # colormap
    strcmap = 'jet'   # colormap
    strcmap = cmocean.cm.thermal
    if varnum==2:
        strcmap = 'RdBu'
    else:
        strcmap = 'PuOr'
#    strcmap = 'gist_earth'

    plotlog = 0
    plotcb=1
    plot_map(data_c, lat_c, lon_c, varnum, domain, \
             strcmap, minmax, save, title, fname, plotlog, plotcb)
