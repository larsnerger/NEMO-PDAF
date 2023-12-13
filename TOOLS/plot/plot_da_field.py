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
import numpy as np
from mpl_toolkits.basemap import Basemap
import time
import scipy.io as spio
#from colormaps import cmaps
from netCDF4 import Dataset as NetCDFFile
from plot_lib import *
import cmocean


if __name__ == "__main__":
    title = 'CHL strongly + SST strongly'
    coupled = ''
    depth  = '0'	# Depth
    year   = 2015       # Year
    month  = '05'       # Month has to be string and two digits i.e '05' and '10'
    day    = '29'	# Day has to be double digits
    ampm   = '00'	# am or pm (string)
    assim  = 'Fcst'     # Free (freerun), Fcst (background/forecast), Ana (analysis)
    domain = 'ba'       # Domain to plot: 'no' for both domains or 'ba' for Baltic, 'Ar' for Arkona
    exp = 'lestkf'
    plotcb = 0          # Whether to show the colorbar
    z_mean=0     #0(no vertical mean)
    z_integral=0 #(0 for no vertical integral; 1 for vertical integral)
    z1=0         # upper Z boundary (for z_mean and z_integral)
    z2=5         # lower Z boundary (for z_mean and z_integral)
    save = 1

    varnum = 21  	# Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP

    if varnum==2:
        if month == '05' and day == '29':
          minmax = [-1, 17.0]
        elif month == '05' and day == '01':
          minmax = [-1, 13.0]
        elif month == '04' and day == '01':
          minmax = [-1, 7]
        elif month == '03' and day == '03':
          minmax = [-1, 6.0]
        plotlog = 0        
    elif varnum==10:
        minmax = [0.1, 10.0]      # DIA
        plotlog = 1
    elif varnum==11:
        minmax = [0.01, 10.0]     # FLA
        plotlog = 1
    elif varnum==12:
        minmax = [0.00001, 0.1]   # CYA
        plotlog = 1
    elif varnum==21:
        minmax = [300.0, 500.0]      # OXY
        plotlog = 0
    elif varnum==23:
        minmax = [6.8 ,  8.2]      # pH
        plotlog = 0
    elif varnum==24:
        minmax = [0.1, 10.0]      # CHL
        plotlog = 1
    elif varnum==25:
	if month == '05':
          minmax = [0, 6]           # TE
	else:
	  minmax = [0,4]
        plotlog = 0
    elif varnum==26:
        minmax = [0, 1]           # PFT
        plotlog = 0
    else:
        minmax = [1, 1]
        plotlog = 0

    ######################################################

    if assim == 'Free': asml = 0
    else: asml = 1
    full_path = get_exp_path(coupled, '2015', '2', '01', '00', asml)
    if 'runfolder' in full_path: 
      model_run = extractBetween(full_path, '_N30/', '/DA/')
    else: 
      model_run = extractBetween(full_path, 'exp/', '/DA/')
    print 'Model run = ', model_run

    varstr, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()

    if assim == 'Free': DAtype = 0
    if assim == 'Fcst': DAtype = 1
    if assim == 'Ana': DAtype = 2

    if coupled=='':
        dotcoupled=coupled
    else:
        dotcoupled='.'+coupled
    data_coarse, coarse_keys = read_model(varstr, 'coarse', year, month, day, ampm, depth, DAtype, dotcoupled, z_mean, z_integral, z1, z2)

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

    # To determine indices for some coordinate
#    for i in range(len(lat1)):
#        if lat1[i]>61.6 and lat1[i]<61.8:
#            print i, lat1[i]
#    for i in range(len(lon1)):
#        if lon1[i]>19.4 and lon1[i]<19.6:
#            print i, lon1[i]

    # Flip arrays for plotting
    lat1 = np.flip(lat1)
    data_coarse = np.flip(data_coarse,0)


    #### Plotting

    if assim=='Ana':
        assimstr = 'analysis'
    elif assim=='Fcst':
        assimstr = 'forecast'
    else:
        assimstr = 'free run'


    print "Prepare for plotting"
    # Set plot title
    #title=str(Variable)+' - '+assimstr+' on '+str(year)+'-'+str(month)+'-'+str(day)

    # Set file name
    if assim=='Free':
      fname= str(MAT_VAR)+'_'+str(assim)+'_'+domain+'_'+\
               str(year)+str(month)+str(day)+'.png'
      fname = varstr+'_'+str(assim)+'_'+str(year)+'-'+str(month)+'-'+str(day)+'.png'
    else:
      fname= str(MAT_VAR)+'_'+exp+'_'+str(assim)+'_'+domain+'_'+\
               str(year)+str(month)+str(day)+'.png'
      fname = varstr+'_'+model_run+'_'+str(assim)+'_'+str(year)+'-'+str(month)+'-'+str(day)+'.png'
    
    if varnum==2:     # TEMP
        strcmap = 'coolwarm' #cmocean.cm.thermal
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

    
    #title = str(Variable)+' on '+str(year)+'-'+str(month)+'-'+str(day)
    #title = 'SST_DA'
    # plot_map(data_coarse, data_fine, lat1, lon1, lat_f, lon_f, varnum, domain, \
    #        strcmap, minmax, save, title, fname)
    plot_map(data_coarse, lat1, lon1, varnum, domain, \
             strcmap, minmax, save, title, fname, plotlog, plotcb)

