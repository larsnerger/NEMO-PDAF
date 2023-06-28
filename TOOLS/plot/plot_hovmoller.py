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
import matplotlib.colors as mcolors
import numpy as np
from mpl_toolkits.basemap import Basemap
import time
import scipy.io as spio
#from colormaps import cmaps
from netCDF4 import Dataset as NetCDFFile
from plot_lib import *
import cmocean


if __name__ == "__main__":
    runname = 'chl+sstL3'
    coupled = ''
    t_start = '20150201'
    t_end = '20150531'
    dist = 10           # Grid point distance to check around central location 
    station = 'ArkonaWR' # Set station name according to the XLSX file
    istation = 0        # (0) use 'station', (1) Arkona, (2) Baltic Proper, (3) Bothnian Sea
    plotlog = 0
    plotcb = 1          # Whether to show the colorbar
    max_depth = 45      # entire depth ranges till ~700m but most places dont have data after around maybe 70m
    #months_da = [1, 2, 3] #3, 4, 5]
    plot_ffa = 'fcst' # 'free', 'fcst' or 'ana'
    z_mean=0     #0(no vertical mean)
    z_integral=0 #(0 for no vertical integral; 1 for vertical integral)
    z1=0         # upper Z boundary (for z_mean and z_integral)
    z2=5         # lower Z boundary (for z_mean and z_integral)
    save = 1
    #minmax = [0.8, 2] # can be out-commented if cbar axis should be set automatically

    varnum = 24	 	# Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP

    ######################################################
    
    # colormap based on variable 

    if varnum==2:     # Temperature
        strcmap = 'coolwarm' 
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

    year = t_start[0:4]
    months_free = range(int(t_start[4:6]), int(t_end[4:6])+1)
    num_days = 0
    for i in months_free:
      _, _, _, days_temp = month_names(i)
      num_days += days_temp
    months_da = months_free
    day = t_start[6:8]
    print 'Starting ', year, t_start[4:6], t_start[6:8]
    ampm = '0'
    print 'Total no. days = ', str(num_days)
    varstr, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()
    
    if plot_ffa != 'free':
      figure_name = 'Hovmoeller_'+runname+'_'+plot_ffa+'_'+varstr+'_'+station+'.png'
    else:
      figure_name = 'Hovmoeller_'+runname+'_'+varstr+'_'+station+'.png'

    if coupled=='':
        dotcoupled=coupled
    else:
        dotcoupled='.'+coupled

    depth = ':' #entire column 

    if plot_ffa == 'free':
      DAtype=0
    elif plot_ffa == 'fcst':
      DAtype = 1
    elif plot_ffa == 'ana':
      DAtype = 2

    data_mean, data_day, dep = read_station_series(varnum, 'coarse', year, ampm, ':', DAtype, dotcoupled, z_mean, z_integral, z1, z2, months_free, istation, station, dist)

    dep = dep[0]

    for depth_cut in range(0,len(dep)):
      if dep[depth_cut]>max_depth:
        break

    print 'Cut depth off after layer no.', str(depth_cut)
    dep = dep[0:depth_cut]
    data_mean = [row[0:depth_cut] for row in data_mean] # delete all data at depths > max_depth

    # Determine tick localizations and labels
    tickloc = []
    tickloc_min = []
    loc = data_day[0]
    tickloc.append(loc)
    months_str = []
    for i in range(len(months_free)):
        _, _, mstr, days_in_mon = month_names(i+months_free[0])
    #    print 'CAUTION hardcoded in plot_hovmoeller days_in_mon=14 and mstr'
    #	days_in_mon = 14
        if t_start=='20150501':
          mstr = 'May 1-14'
        elif t_start == '20150301':
	  mstr = 'March 1-14'
        loc = loc + days_in_mon
        tickloc.append(loc)
        tickloc_min.append(tickloc[i]+days_in_mon/2)
        months_str.append(mstr)
    
    print 'Min max data values = ', np.min(data_mean), np.max(data_mean)
    #### Plotting
    yticks = [0, 10, 20, 30, 40]
    print 'data_day'
    print data_day
    fig, ax = plt.subplots(figsize=(8,3))
    X, Y = np.meshgrid(data_day, dep)
    if plotlog == 1:
      im = ax.pcolormesh(X,Y,transpose(data_mean), cmap=strcmap, norm=mcolors.LogNorm(vmin=0.1, vmax=10))
    else:
      if 'minmax' in locals():
        im = ax.pcolormesh(X,Y,transpose(data_mean), cmap=strcmap, vmin=minmax[0], vmax=minmax[1])
      else:
	im = ax.pcolormesh(X,Y,transpose(data_mean), cmap=strcmap)
    if plotcb == 1: 
      cbar = fig.colorbar(im)
      cbar.set_label(var_unit)
    levels1 = [1.2, 2.0, 5.0]
    #levels1 = [6.0, 6.5, 7.0, 7.5, 8.0, 8.5]
    #levels1 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]#, 7.0, 8.0]
    #levels1 = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    #levels1 = [6.0, 7.0, 8.0]
    #levels1 = [3.0, 3.5, 4.0, 4.5, 5.0]
    #levels2 = [6.5, 7.5, 8.5]
    contour1 = ax.contour(X, Y, transpose(data_mean), levels1, colors='black') # add contour lines
    plt.clabel(contour1, inline=True, fontsize=10, fmt='%.1f')
    #contour2 = ax.contour(X, Y, transpose(data_mean), levels2, colors='black') # add contour lines
    #plt.clabel(contour2, inline=True, fontsize=10, fmt='%.2f')
    if plot_ffa == 'free':
      ax.set_title(Variable + ' for Free run')
    elif plot_ffa == 'fcst':
      ax.set_title(Variable + ' for DA - Forecast')
    elif plot_ffa == 'ana':
      ax.set_title(Variable + ' for DA - Analysis')
    ax.set_ylabel('depth [m]')
    ax.set_yticks(yticks)
    plt.xlim(tickloc[0],tickloc[-1])
    ax.xaxis.set_major_locator(ticker.FixedLocator(tickloc))
    ax.xaxis.set_minor_locator(ticker.FixedLocator(tickloc_min))
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.tick_params(axis='x', which='minor', tick1On=False, tick2On=False)
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(months_str))
    ax.invert_yaxis() # flip y-axis
    if save == 1: 
      print 'Saving figure'
      plt.savefig(figure_name)
    plt.show()
