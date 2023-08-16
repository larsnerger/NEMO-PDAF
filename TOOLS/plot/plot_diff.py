#!/usr/bin/python
#################
# Script to plot the differences
#    analysis - forecast
# or
#    analysis-free 
# or
#    forecast-free
#################

import sys, getopt
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
import numpy as np
from mpl_toolkits.basemap import Basemap
import time
import scipy.io as spio
#from colormaps import cmaps
from netCDF4 import Dataset as NetCDFFile
from plot_lib import *
import cmocean


if __name__ == "__main__":
    title = 'CHL weakly + SST strongly'
    depth  = '0'	# Depth level
    year   = 2015       # Year
    month  = '05'       # Month has to be string and two digits i.e '05' and '10'
    day    = '29'	      # Day has to be double digits
    ampm   = '00'	     # am or pm (string)  - not relevant for NEMO
    assim  = 'Fcst'     # Fcst (background/forecast), Ana (analysis)
    assim2  = 'Free'     # Free (freerun), Fcst (background/forecast), Ana (analysis)
    domain = 'ba'       # Domain to plot: 'no' for both domains or 'ku' for fine only
    minmax = [-20.0,20.0]      # max/max plotted values - set min=max for automatic
    save = 1
    plotcb = 0         # plot colorbar (0-no, 1-yes)
    varnum = 21	 	# Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP

    ######################################################

    asml = 1
    coupled = 'weakly' #irrelevant variable but needs to be specified
    full_path = get_exp_path(coupled, '2015', '2', '01', '00', asml)
    model_run = extractBetween(full_path, 'exp/', '/DA/')
    print 'Model run = ', model_run

    if varnum == 2:
      minmax = [-1.5, 1.5]
    elif varnum == 21:
      minmax = [-20, 20]

    varstr, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()

    if assim == 'Free': DAtype = 0
    if assim == 'Fcst': DAtype = 1
    if assim == 'Ana': DAtype = 2

    if coupled=='':
        dotcoupled=coupled
    else:
        dotcoupled='.'+coupled

    # Read Forecast
    data_coarse1, coarse_keys = read_model(varstr, 'coarse', year, month, day, ampm, depth, DAtype, dotcoupled, \
                                           0, 0, 0, 0)

    # Read analysis
    if assim2 == 'Free': DAtype2 = 0
    if assim2 == 'Fcst': DAtype2 = 1
    if assim2 == 'Ana': DAtype2 = 2
    data_coarse2, coarse_keys2 = read_model(varstr, 'coarse', year, month, day, ampm, depth, DAtype2, dotcoupled,\
                                            0, 0, 0, 0)

    data_coarse = data_coarse1 - data_coarse2

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

    lat1 = np.flip(lat1)
    data_coarse = np.flip(data_coarse,0)

    #### Plotting

    # Set plot title
    if 'title' not in locals():
      if assim2=='Ana' or assim2=='Fcst':
         title='Assimiliation increment for '+str(Variable)+' on '+str(year)+'-'+str(month)+'-'+str(day)
      elif assim2=='Free' and assim=='Fcst':
         title='Difference forecast-free for '+str(Variable)+' on '+str(year)+'-'+str(month)+'-'+str(day)
      elif assim2=='Free' and assim=='Ana':
         title='Difference analysis-free for '+str(Variable)+' on '+str(year)+'-'+str(month)+'-'+str(day)
    # Set file name
    fname= str(MAT_VAR)+'_'+model_run+'_'+assim+'-'+assim2+'_'+str(year)+str(month)+str(day)+\
			'_'+domain+'.png'
    print 'filename', fname
    
    strcmap = 'coolwarm'   # colormap
    if varnum == 2:
      strcmap='seismic'
    elif varnum == 21:
      strcmap = 'cmo.delta'

    plotlog = 0
    plot_map(data_coarse, lat1, lon1, varnum, domain, \
             strcmap, minmax, save, title, fname, plotlog, plotcb)

