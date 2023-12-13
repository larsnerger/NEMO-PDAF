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
from datetime import datetime 


if __name__ == "__main__":
    title = 'PH'
    tstart = '20151101' # start day (in format: 'yyyymmdd')
    tend = '20151130'   # end day (in format: 'yyyymmdd')
    coupled = ''
    depth  = '0'	# Depth
    ampm   = '00'	# am or pm (string)
    domain = 'ba'       # Domain to plot: 'no' for both domains or 'ba' for Baltic, 'Ar' for Arkona
    exp = 'lestkf'
    plotcb = 1          # Whether to show the colorbar
    z_mean=0     #0(no vertical mean)
    z_integral=0 #(0 for no vertical integral; 1 for vertical integral)
    z1=0         # upper Z boundary (for z_mean and z_integral)
    z2=5         # lower Z boundary (for z_mean and z_integral)
    save = 1

    varnum_range = [23] 	# Variable number from the var_names routine.
                                # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                                #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                                #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                                #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP

    for varnum in varnum_range:

    ######################################################
      print tstart[5]
      if varnum == 24: 
	if tstart[5] == '4':
	  minmax = [-10,10]
	elif tstart[5] == '3':
  	  minmax = [-10, 10]
	elif tstart[5] == '2':
	  minmax = [-3.5, 3.5]
	elif tstart[5] == '5':
	  minmax = [1,1]
        elif tstart[5] == '7':
          minmax = [-13, 13]
        elif tstart[5] == '9': 
          minmax = [-8, 8]
        elif tstart[4:6] == '11':
	  minmax = [-5, 5]
        else: 
          minmax = [1,1]
      elif varnum == 2: 
        if tstart[5] == '7' or tstart[4:6] == '11':
          minmax = [-1.5,1.5]
        else:
	  minmax = [-1.25, 1.25]
      elif varnum == 21: 
	minmax = [-20, 20]
        minmax = [-75, 75]
        minmax = [-2000, 2000]
      elif varnum == 27: 
        minmax = [1,1]
      elif varnum == 26: 
        if tstart[4:6] == '07':
          minmax = [-0.75, 0.75]
        else: 
          minmax = [1,1]
      elif varnum == 25:
        minmax = [-3, 3]
      elif varnum == 23:
        if tstart[4:6] == '07':
          minmax = [-0.55, 0.55]
	elif tstart[4:6] == '09':
	  minmax = [-0.2, 0.2]
	elif tstart[4:6] == '11':
	  minmax = [-0.1, 0.1]

      if 'minmax' not in locals():
        minmax = [1,1]
      date_range = generate_date_range(tstart, tend)
      
      loopcnt = 0
      data_list = []
      free_list = []


      for d in date_range:
        full_path = get_exp_path(coupled, '2015', '2', '01', '00', 1)
        if 'runfolder' in full_path: 
          model_run = extractBetween(full_path, '_N30/', '/DA/')
        else: 
          model_run = extractBetween(full_path, 'exp/', '/DA/')
        print 'Model run = ', model_run

        varstr, mat_var, Variable, var_unit = var_names(varnum)
        MAT_VAR = mat_var.upper()

        if coupled=='':
          dotcoupled=coupled
        else:
          dotcoupled='.'+coupled
        data_coarse, coarse_keys = read_model(varstr, 'coarse', d[0:4], d[4:6], d[6:8], ampm, depth, 1, dotcoupled, z_mean, z_integral, z1, z2)
        free_coarse, free_keys = read_model(varstr, 'coarse', d[0:4], d[4:6], d[6:8], ampm, depth, 0, dotcoupled, z_mean, z_integral, z1, z2)
        # flip array for plotting
        data_coarse = np.flip(data_coarse,0)
        free_coarse = np.flip(free_coarse,0)
        data_list.append(data_coarse)
        free_list.append(free_coarse)
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

      #### Plotting
      data_mean = np.ma.mean(data_list, axis=0)
      free_mean = np.ma.mean(free_list, axis=0)

      data_mean = data_mean - free_mean

      print "Prepare for plotting"
      # Set plot title
      #title=str(Variable)+' - '+assimstr+' on '+str(year)+'-'+str(month)+'-'+str(day)

      # Set file name
      fname = varstr+'_diff_free_'+model_run+'_'+tstart+'-'+tend+'.png'
    
      if varnum==2:     # TEMP
        strcmap = 'seismic' #cmocean.cm.thermal
      elif varnum==10:
        strcmap = 'viridis'
      elif varnum==11:
        strcmap = 'inferno'
      elif varnum==12:
        strcmap = cmocean.cm.deep
      elif varnum==21:     # OXY
        strcmap = 'cmo.delta'
      elif varnum==23:     # pH
        strcmap = 'PiYG'
      elif varnum==24:     # Chlorophyll
        strcmap = 'cmo.curl'
      elif varnum==25:   # TE
        strcmap = cmocean.cm.tarn
      elif varnum==26:   # PFT
        strcmap = cmocean.cm.diff
      else:
        strcmap = 'gist_ncar_r'

      if minmax[0] == minmax[1]:
        print 'Inside loop bcs', minmax[0], minmax[1]
	minmaxabs = max(abs(np.min(data_mean)), abs(np.max(data_mean)))
   	minmax = [-minmaxabs, minmaxabs] 
      print 'minmax = ', minmax
      #title = str(Variable)+' on '+str(year)+'-'+str(month)+'-'+str(day)
      #title = 'SST_DA'
      # plot_map(data_coarse, data_fine, lat1, lon1, lat_f, lon_f, varnum, domain, \
      #        strcmap, minmax, save, title, fname)
      plotlog = 0
      plot_map(data_mean, lat1, lon1, varnum, domain, \
             strcmap, minmax, save, title, fname, plotlog, plotcb)

