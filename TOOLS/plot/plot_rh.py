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
    year   = 2015       # Year
    month  = '04'       # Month has to be string and two digits i.e '05' and '10'
    firstday = 1        # Initial day 
    ndays  = 1          # Number of days to process
    dimens = 30         # Ensemble size
    dist = 10           # Grid point distance to check around central location 
    forana = 'a'        # (a) analysis, (f) forecast
    istation = 1        # (1) Arkona, (2) Baltic Proper, (3) Bothnian Sea
    varnum = 4 		# Variable number from the var_names routine: 1-SST, 4-CHL
    save   = 0          # (1) save file

    ####################################################
    # Ignore these variables
    var, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()
    ######################################################


#    path_exp = '/scratch/usr/hzfblner/SEAMLESS/run/ens_ERGOM_DA/001/output/data_mobs'
    path_exp = '/scratch/projects/hbk00095/exp/LESTKF/data_mobs'
    datatype = 'lestkf_N30'
    datatype_model = 'lestkf_N30'
    exp = 'LESTKF'


    if istation == 1:
       station = 'Arkona'
       station_lon = 260   # 13.87 deg E
       station_lat = 148   # 54.88 deg N
    elif istation == 2:
       station = 'Baltic Proper'
       station_lon = 520 # 18.5 deg E
       station_lat = 247 # 56. 0 deg N
    else:
       station= 'Bothnian Sea'
       station_lon = 580 # 19.5 deg E
       station_lat = 770 # 61.7 deg N
  

    
    lon_min = station_lon-dist
    lon_max = station_lon+dist+1
    lat_min = station_lat-dist
    lat_max = station_lat+dist+1
    range_lon = lon_max-lon_min
    range_lat = lat_max-lat_min

    hist = np.zeros([dimens+1])
    cnt = 0

    print 'EXP: ', path_exp

    for iday in range(ndays):

        day = iday+firstday

        print 'Process day ', day
        # Read observations
        log = 0
        datatype_obs= 'ba_MY'
        data_o, lat_o, lon_o = get_cmems_chl(year, month, day, 'am', datatype_obs, log)

        if day==firstday:
            print 'check observations within ', dist, ' grid points around location ', station
            print 'longitude range: ', lon_o[lon_min], lon_o[lon_max]
            print 'longitude range: ', lat_o[lat_min], lat_o[lat_max]
        
        obsval= np.zeros((range_lon, range_lat))
        obsval = data_o[lat_min:lat_max, lon_min:lon_max]

        # Read observed model ensemble
        mobsval = np.zeros((range_lon, range_lat, dimens))

        for iens in range(dimens):
            print 'ENS: ', iens, ' forana', forana
            log = 0
            data_m, lat_c, lon_c = get_mobs_cmems_chl(year, month, day, forana, path_exp, datatype_model, log, iens+1)
            mobsval[:,:,iens] = data_m[lat_min:lat_max, lon_min:lon_max]
            
        # Compute rank histogram
        cnt0=0
        for i in range(len(obsval[:,0])):
            for j in range(len(obsval[0,:])):
                oval = obsval[i,j]
                movals = mobsval[i,j,:]
                movals = np.sort(movals)

                rank = 0
                if obsval.mask[i,j]==False:
                    cnt = cnt + 1
                    for iens in range(dimens):
                        if movals[iens] < oval:
                            cnt0 = cnt0 + 1
                            rank = rank + 1

                if cnt0>0:
                    hist[rank] = hist[rank]+1

    	print 'accumulated observation points: ', cnt

    print 'Number of observation points: ', cnt
    print 'Histogram: ', hist


    xval = np.zeros([dimens+1])
    for i in range(dimens+1):
        xval[i] = i

    if day<10:
        nullstr='0'
    else:
        nullstr=''

    fname = 'hist_'+exp+'-'+str(forana)+'_'+str(year)+str(month)+nullstr+str(day)+'.png'
    print ('File '+fname)

    # Plot histogram
    lastday = firstday+ndays-1
    if month=='04':
    	mstr = 'April '+str(year)
    else:
    	mstr = 'May '+str(year)
    if forana=='a':
        foranastr='analysis'
    else:	
        foranastr='forecast'

    if 1==1:
       fig = plt.figure(figsize=(6,3))
       plt.bar(xval,hist/cnt)
       if ndays==1:
           plt.title(station+': '+exp+' '+foranastr+' - '+str(firstday)+' '+mstr)
       else:
           plt.title(station+': '+exp+' '+foranastr+' - '+mstr+' days '+str(firstday)+'-'+str(lastday))

       if save==1:
	    plt.savefig(fname, dpi=300)
       else:
	    plt.show()
    
    if 1==2:
       fig = plt.figure
       for iens in range(dimens):
          p_mobsval = mobsval[:,:,iens]
          plt.plot(p_mobsval.reshape(range_lon*range_lat),'^')
       plt.plot(obsval.reshape(range_lon*range_lat),'ko')
       plt.show()


