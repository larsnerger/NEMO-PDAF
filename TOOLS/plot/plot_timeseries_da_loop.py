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
  exp = 'LESTKF'
  
  plot_fcst = 1      # 0 - analysis, 1 - forecast
  plot_ana = 0       # if only assim1 is to be plotted (all others are out-commented), then one can additionally plot analysis to forecast
  assim1 = 'sstL3_strongly_N30'    # DA-run1 - exp name as in project folder
  assim2 = 'chl_strongly_N30'    # DA-run2 (can be out-commented)
  assim3 = 'sst_weakly_N30'    # DA-run3 (can be out-commented)
  assim4 = 'chl_weakly_N30'    # DA-run4 (can be out-commented) 

  plotlog = 0
  plotcb = 1          # Whether to show the colorbar

  months_free = [1, 2, 3, 4, 5] #, 2, 3, 4, 5]
  months_obs = [1, 2, 3, 4, 5] #, 3] #, 4, 5]
  months_da = [2, 3, 4, 5] #, 2, 3] #3, 4, 5]
    
  z_mean=0     #0(no vertical mean)
  z_integral=0 #(0 for no vertical integral; 1 for vertical integral)
  z1=0         # upper Z boundary (for z_mean and z_integral)
  z2=5         # lower Z boundary (for z_mean and z_integral)
  save = 1


#  varnum = 11	 	# Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP

  for varnum in range(2,25):


    ######################################################

    varstr, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()

    if coupled=='':
        dotcoupled=coupled
    else:
        dotcoupled='.'+coupled

    DAtype=0
    freemean, freeday = read_station_series(varnum, 'coarse', year, ampm, depth, DAtype, dotcoupled, z_mean, z_integral, z1, z2, months_free, istation, station, dist)

    if plot_fcst == 1: 
      DAtype=1
    else: 
      DAtype = 0

    basedir_pro = '/scratch/projects/hbk00095/exp/'
    basedir_work = '/scratch/usr/hbxsovli/SEAMLESS/runfolders/'
    if varnum == 25 or varnum == 26:
      path1 = basedir+assim1+'/Post_DA/'
    else:
      #path1 = basedir_pro + assim1 + '/DA/'
      path1 = basedir_work + assim1 + '/' + assim1 + '/DA/'
    assimmean1, assimday1 = read_station_series_path(varnum, 'coarse', year, ampm, depth, DAtype, dotcoupled, z_mean, z_integral, z1, z2, months_da, istation, station, dist, path1)
    assimmean1 = np.array(assimmean1)
    if 'assim2' in locals(): 
      print '\n \n \n ----------------------------------------------------------------------------------------------------'
      if varnum == 26 or varnum == 25:
        path2 = basedir+assim2+'/Post_DA/'
      else:
        #path2 = basedir_pro + assim2 + '/DA/'
        path2 = basedir_work + assim2 + '/' + assim2 + '/DA/'
      assimmean2, assimday2 = read_station_series_path(varnum, 'coarse', year, ampm, depth, DAtype, dotcoupled, z_mean, z_integral, z1, z2, months_da, istation, station, dist, path2)
      assimmean2 = np.array(assimmean2)
    else: 
      print '\n \n \n ----------------------------------------------------------------------------------------------------'
      print 'NO'
    if 'assim3' in locals(): 
      if varnum == 26 or varnum == 25:
        path3 = basedir+assim3+'/Post_DA/'
      else:
        path3 = basedir_pro + assim3 + '/DA/'
        #path3 = basedir_work + assim3 + '/' + assim3 + '/DA/'
      assimmean3, assimday3 = read_station_series_path(varnum, 'coarse', year, ampm, depth, DAtype, dotcoupled, z_mean, z_integral, z1, z2, months_da, istation, station, dist, path3)
      assimmean3 = np.array(assimmean3)
    if 'assim4' in locals(): 
      if varnum == 26 or varnum == 25:
        path4 = basedir+assim4+'/Post_DA/'
      else:
        path4 = basedir_pro + assim4 + '/DA/'
        #path4 = basedir_work + assim4 + '/' + assim4 + '/DA/'
      assimmean4, assimday4 = read_station_series_path(varnum, 'coarse', year, ampm, depth, DAtype, dotcoupled, z_mean, z_integral, z1, z2, months_da, istation, station, dist, path4)
      assimmean4 = np.array(assimmean4)



    if varnum==2 or varnum==24:
        obsmean, obsday = read_station_series_obs(varnum, year, months_obs, istation, station, dist)
        print 'size obsmean=', size(obsmean)
	print 'size obsday=', size(obsday)

        ### RMSE
        # make sure the same time span is covered to calculate RMSE
        # convert everything into numpy arrays so one can use masks
        obsday = np.array(obsday)
        freeday = np.array(freeday)
        obsmean = np.array(obsmean)
        freemean = np.array(freemean)
        # assim & free 
        obs_mask = np.in1d(obsday, np.array(assimday1))
        free_mask = np.in1d(freeday, np.array(assimday1))
        freeday_temp = freeday[free_mask]
        freemean_temp = freemean[free_mask]
        obsday_temp = obsday[obs_mask]
        obsmean_temp = obsmean[obs_mask]
        # obs & model 
        model_mask = np.in1d(freeday_temp, obsday_temp)
        rmse_free = calc_rmse(obsmean_temp, freemean_temp[model_mask])
        rmse_assim1 = calc_rmse(obsmean_temp, assimmean1[model_mask]) 
        print '\n ---------- RMSE -----------'
        print 'Free: ', rmse_free
        print assim1, ': ', rmse_assim1
        if 'assim2' in locals():
          rmse_assim2 = calc_rmse(obsmean_temp, assimmean2[model_mask])
          print assim2, ': ', rmse_assim2
        if 'assim3' in locals(): 
          rmse_assim3 = calc_rmse(obsmean_temp, assimmean3[model_mask])
          print assim3, ': ', rmse_assim3
        if 'assim4' in locals(): 
          rmse_assim4 = calc_rmse(obsmean_temp, assimmean4[model_mask])
          print assim4, ': ', rmse_assim4


    #### Plotting

    # Determine tick localizations and labels
    tickloc = []
    tickloc_min = []
    loc = freeday[0]
    tickloc.append(loc)
    months_str = []
    for i in range(len(months_free)):
        _, _, mstr, days_in_mon = month_names(i+months_free[0])
        loc = loc + days_in_mon
        tickloc.append(loc)
        tickloc_min.append(tickloc[i]+days_in_mon/2)
        months_str.append(mstr)

    if station == 'DarsserS':
      station = 'Darsser Sill' #english name of the station

    fname= str(MAT_VAR)+'_series_'+exp+'_'+station.replace(' ','_')+'_'+str(year)+'_'+str(months_free[0])+'-'+str(months_free[-1])+'.png'
    print ('File '+fname)

    if 1==1:
        fig, ax = plt.subplots(figsize=(8,3))
        plt.plot(freeday, freemean, 'b',label='FREE')
        #plt.plot(assimday1, assimmean1,'#bbbb00',label=assim1)
        plt.plot(assimday1, assimmean1,'#bbbb00',label='SST strongly')
        if 'assim2' in locals(): 
          #plt.plot(assimday2, assimmean2,'g',label=assim2)
          plt.plot(assimday2, assimmean2,'g',label='CHL strongly')
        if 'assim3' in locals():
          #plt.plot(assimday3, assimmean3,'k',label=assim3)
          plt.plot(assimday3, assimmean3,'k',label='SST weakly')
        if 'assim4' in locals():
          #plt.plot(assimday4, assimmean4,'k',label=assim4)
          plt.plot(assimday4, assimmean4,'y',label='CHL weakly')
        if plot_ana==1:
            plt.plot(anaday, anamean,'k',label='analysis')
        if varnum==2 or varnum==24:
            plt.plot(obsday, obsmean,'r+',label='obs')
        plt.ylabel(Variable+' ('+var_unit+')')
        plt.xlim(tickloc[0],tickloc[-1])
        if varnum == 11 or varnum == 13 or varnum == 9 or varnum == 8:
	  plt.ylim(ymin=0.0)
        ax.xaxis.set_major_locator(ticker.FixedLocator(tickloc))
        ax.xaxis.set_minor_locator(ticker.FixedLocator(tickloc_min))
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.tick_params(axis='x', which='minor', tick1On=False, tick2On=False)
        ax.xaxis.set_minor_formatter(ticker.FixedFormatter(months_str))        
        plt.legend()
#        plt.title(Variable+' at '+station+': '+exp)
        plt.title(Variable, fontsize=16)

        if save==1:
 	    plt.savefig(fname, dpi=300)
        else:
 	    plt.show()


