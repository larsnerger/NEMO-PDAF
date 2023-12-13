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
    ens = range(21,28)
    #ens = [10, 11, 14, 15, 16, 17] #range(1,7) # which ensemble state 
    #ens = [18, 19, 20, 21, 26, 27]
    #ens = [28, 29, 30]
    tstart = '20150301' # start day (in format: 'yyyymmdd')
    tend = '20150331'   # end day (in format: 'yyyymmdd')
    coupled = ''
    depth  = '0'	# Depth
    ampm   = '00'	# am or pm (string)
    assim  = 'Fcst'     # Free (freerun), Fcst (background/forecast), Ana (analysis)
    domain = 'ba'       # Domain to plot: 'no' for both domains or 'ba' for Baltic, 'Ar' for Arkona
    exp = 'lestkf'
    plotcb = 1          # Whether to show the colorbar
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
        if tstart[4:6] == '02':
          minmax = [-1, 6.0]
        elif tstart[4:6] == '03':
          minmax = [-1, 6.0]
        elif tstart[4:6] == '04':
          minmax = [-1, 10.0]
        elif tstart[4:6] == '05':
          minmax = [-1, 14.0]
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
        minmax = [-3.5, 3.5]      # OXY
        plotlog = 0
    elif varnum==23:
        minmax = [6.8 ,  8.2]      # pH
        plotlog = 0
    elif varnum==24:
        minmax = [0.1, 10.0]      # CHL
        plotlog = 1
    elif varnum==25:
        if tstart[4:6] == '05':
          minmax = [0, 6]           # TE
        else:
          minmax = [0,4]
          plotlog = 0
    elif varnum==26:
        minmax = [0, 1]           # PFT
        plotlog = 0
    else:
    #    minmax = [1, 1]
        plotlog = 0

    ######################################################

    
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
        strcmap = 'gist_ncar_r'
    elif varnum==24:     # Chlorophyll
        strcmap = 'cmo.curl'
    elif varnum==25:   # TE
        strcmap = cmocean.cm.dense
    elif varnum==26:   # PFT
        strcmap = cmocean.cm.rain
    date_range = generate_date_range(tstart, tend)

    loopcnt = 0
    data_list = []

    if len(ens) % 6 != 0:
      print '!!! ERROR needs to be divisible by 6 !!!'  

    ens_cnt = 0 
    
    for d in date_range: 
      varstr, mat_var, Variable, var_unit = var_names(varnum)
      MAT_VAR = mat_var.upper()

      if assim == 'Free': DAtype = 0
      if assim == 'Fcst': DAtype = 1
      if assim == 'Ana': DAtype = 2

      if coupled=='':
        dotcoupled=coupled
      else:
        dotcoupled='.'+coupled
      exp = '/scratch-emmy/projects/hbk00095/exp/free_scaleparams/data/'
      path = exp + '001_NORDIC_1d_ERGOM_T_' + d + '-' + d + '.nc' #default
      #data_coarse, coarse_keys = read_model(varstr, 'coarse', d[0:4], d[4:6], d[6:8], ampm, depth, DAtype, dotcoupled, z_mean, z_integral, z1, z2)
      data_coarse, coarse_keys = read_model_path(varstr, 'coarse', d[0:4], d[4:6], d[6:8], ampm, depth, DAtype, dotcoupled, z_mean, z_integral, z1, z2, path)
      # flip array for plotting
      data_coarse = np.flip(data_coarse,0)
      data_list.append(data_coarse)

    default_mean = np.ma.mean(data_list, axis=0)

    for e in ens:
      data_list = []
      print 'Start ens = ', e
      for d in date_range:
        if e < 10: 
          ensstr = '00'+str(e)
        else: 
          ensstr = '0'+str(e)
        exp = '/scratch-emmy/projects/hbk00095/exp/free_scaleparams/data/'
        path = exp + ensstr + '_NORDIC_1d_ERGOM_T_' + d + '-' + d + '.nc'
        #data_coarse, coarse_keys = read_model(varstr, 'coarse', d[0:4], d[4:6], d[6:8], ampm, depth, DAtype, dotcoupled, z_mean, z_integral, z1, z2)
        data_coarse, coarse_keys = read_model_path(varstr, 'coarse', d[0:4], d[4:6], d[6:8], ampm, depth, DAtype, dotcoupled, z_mean, z_integral, z1, z2, path)
        # flip array for plotting
        data_coarse = np.flip(data_coarse,0)
        data_list.append(data_coarse)
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
      data_mean = np.ma.mean(data_list, axis=0) - default_mean

      if assim=='Ana':
        assimstr = 'analysis'
      elif assim=='Fcst':
        assimstr = 'forecast'
      else:
        assimstr = 'free run'

      if ens_cnt % 6 == 0: # new subplot  
        fname = varstr+'_ens'+ensstr+'_'+tstart+'-'+tend+'.png'
        # Create a figure with 3x2 subplots
        varstr, mat_var, Variable, var_unit = var_names(varnum)

        if domain=='no':
          llon=lon_c[0]
          llat=lat_c[-1]
          ulon=lon_c[-1]
          ulat=lat_c[0]
          res = 5.0  # resolution for meridians & parallels
          bmres = 'l'
        elif domain=='ba':
          llon = 9.25
          ulon = 30.375
          llat = 53.229
          ulat = 66.0 #65.875
          res = 3.0
          bmres = 'i'
        else:
          llon = 9.25
          ulon = 18.0
          llat = 53.229
          ulat = 60.0 #65.875
          res = 2.0
          bmres = 'i'

        fig, axes = plt.subplots(2,3,figsize=(12,12))
      if ens_cnt == 0:
        ax = axes[0,0]
      elif ens_cnt == 1: 
        ax = axes[0,1]
      elif ens_cnt == 2:
        ax = axes[0,2]
      elif ens_cnt == 3: 
        ax = axes[1,0]
      elif ens_cnt == 4: 
        ax = axes[1,1]
      elif ens_cnt == 5: 
        ax = axes[1,2]
      m = Basemap(projection='merc',
                llcrnrlon=llon,
                llcrnrlat=llat,
                urcrnrlon=ulon,
                urcrnrlat=ulat,
                lon_0=(llon + ulon)/2.0,
                lat_0=(ulat + llat)/2.0,
                resolution=bmres, area_thresh = 10, ax=ax)

      # Map Land Details
      m.drawcoastlines()                      # Coastlines
      m.drawcountries()                       # Country Borders
      m.drawlsmask(land_color='Linen')       # No geographic detail
      #m.shadedrelief()                        # Low geographic detail
      #m.bluemarble()                         # High geographic detail

      # draw parallels.
      parallels = np.arange(0.,360.,res)
      m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
      # draw meridians
      meridians = np.arange(0.,360.,res)
      m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

      if 'minmax' not in locals():
        vmin = np.min(data_mean)
        vmax = np.max(data_mean)
      else:
        vmin = minmax[0]
        vmax = minmax[1]
      freeclim = 0
      print 'Limit color range to: ', vmin, 'to', vmax
      cmap = plt.cm.get_cmap(strcmap)
    #    cmap1 = plt.cm.get_cmap('terrain_r')
      clev = np.linspace(vmin, vmax, 100)

        # Data - Coarse
      lons,lats = np.meshgrid(lon1, lat1)
      x, y = m(lons,lats)
      
      if vmin < 0.0:
          if freeclim==1:
            cs1 = m.contourf(x, y, data_mean, cmap=cmap, vmin=vmin, vmax=vmax, \
                             levels=clev)#, norm=MidpointNormalize(midpoint=0.0))
          else:
            cs1 = m.pcolor(x, y, data_mean, cmap=cmap, vmin=vmin, vmax=vmax)#, \
                           #norm=MidpointNormalize(midpoint=0.0))
      else:
#        if freeclim==1:
#            cs1 = m.contourf(x, y, data_coarse, cmap=cmap, vmin=vmin, vmax=vmax, levels=clev)
#        else:
           if plotlog==1:
             print 'Plot log values'
             cs1 = m.pcolor(x, y, data_mean, cmap=cmap, vmin=vmin, vmax=vmax, norm=colors.LogNorm(vmin=vmin, vmax=vmax))
           else:
             cs1 = m.pcolor(x, y, data_mean, cmap=cmap, vmin=vmin, vmax=vmax)

      if ens_cnt % 6 == 0:
        if plotcb==1:
        #    tick_locator = ticker.MaxNLocator(nbins=5)
          cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
          cbar = plt.colorbar(cs1, cax = cbar_ax)
        #    cbar.locator = tick_locator
        #    cbar.update_ticks()
          cbar.set_label(var_unit,fontsize=18)
          cbar.ax.tick_params(labelsize=16)
      
      title = str(e)
      ax.set_title(title, fontsize=16)
 
      ens_cnt += 1
      if ens_cnt == 6:
        ens_cnt = 0
        if save==1:
            plt.savefig(fname, dpi=200)
        else:
            plt.show()

	ens_cnt = 0
