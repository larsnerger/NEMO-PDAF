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

    months_free = [1, 2, 3] #, 2, 3, 4, 5]
    #months_obs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    months_obs = [1, 2, 3] #, 3] #, 4, 5]
    months_da = [2, 3] #3, 4, 5]
    
    z_mean=0     #0(no vertical mean)
    z_integral=0 #(0 for no vertical integral; 1 for vertical integral)
    z1=0         # upper Z boundary (for z_mean and z_integral)
    z2=5         # lower Z boundary (for z_mean and z_integral)
    save = 0


    varnum = 2	 	# Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP

    if varnum==24 or varnum==10 or varnum==11 or varnum==12:
        plotana = 1
    else:
        plotana = 0


    ######################################################

    varstr, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()

    if coupled=='':
        dotcoupled=coupled
    else:
        dotcoupled='.'+coupled

    DAtype=0
    freemean, freeday = read_station_series(varnum, 'coarse', year, ampm, depth, DAtype, dotcoupled, z_mean, z_integral, z1, z2, months_free, istation, station, dist)

    DAtype=1
    fcstmean, fcstday = read_station_series(varnum, 'coarse', year, ampm, depth, DAtype, dotcoupled, z_mean, z_integral, z1, z2, months_da, istation, station, dist)

    if plotana==1:
        DAtype=2
        anamean, anaday = read_station_series(varnum, 'coarse', year, ampm, depth, DAtype, dotcoupled, z_mean, z_integral, z1, z2, months_da, istation, station, dist)

    if varnum==2 or varnum==24:
        obsmean, obsday = read_station_series_obs(varnum, year, months_obs, istation, station, dist)

    #print freemean[0:10]
    #print fcstmean[0:10]
    #print anamean[0:10]

    #### Plotting

    # Determine tick locations and labels
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


    fname= str(MAT_VAR)+'_series_'+exp+'_'+station+'_'+str(year)+'_'+str(months_free[0])+'-'+str(months_free[-1])+'.png'
    print ('File '+fname)

    if 1==1:
        fig, ax = plt.subplots(figsize=(8,3))
        plt.plot(freeday, freemean, 'b',label='free')
        plt.plot(fcstday, fcstmean,'#bbbb00',label='forecast')
        if plotana==1:
            plt.plot(anaday, anamean,'k',label='analysis')
        if varnum==2 or varnum==24:
            plt.plot(obsday, obsmean,'r+',label='obs')
        plt.ylabel(Variable+' ('+var_unit+')')
        plt.xlim(tickloc[0],tickloc[-1])
        ax.xaxis.set_major_locator(ticker.FixedLocator(tickloc))
        ax.xaxis.set_minor_locator(ticker.FixedLocator(tickloc_min))
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.tick_params(axis='x', which='minor', tick1On=False, tick2On=False)
        ax.xaxis.set_minor_formatter(ticker.FixedFormatter(months_str))        
        plt.legend()
        plt.title(Variable+' at '+station+': '+exp)


        if save==1:
 	    plt.savefig(fname, dpi=300)
        else:
 	    plt.show()


