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
import pickle
from datetime import timedelta, datetime


if __name__ == "__main__":
    coupled = ''
    fname = 'OlandSodraBS.png'
    year_p = 2015
    lat = 56.17
    lon = 16.64
    plot_save = 1
    #####################################################
    domain = 'ba'

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

    fig = plt.figure(figsize=(8,8))
    m = Basemap(projection='merc',
                llcrnrlon=llon,
                llcrnrlat=llat,
                urcrnrlon=ulon,
                urcrnrlat=ulat,
                lon_0=(llon + ulon)/2.0,
                lat_0=(ulat + llat)/2.0,
                resolution=bmres, area_thresh = 10)
    # Map Land Details
    m.drawcoastlines()                      # Coastlines
    m.drawcountries()                       # Country Borders
    m.drawlsmask(land_color='Linen')       # No geographic detail
    #m.shadedrelief()                        # Low geographic detail
    #m.bluemarble()                         # High geographic detail

    # draw parallels.
    parallels = np.arange(0.,360.,res)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14)
    # draw meridians
    meridians = np.arange(0.,360.,res)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14)


    x, y = m(lon,lat)

#    cs = m.scatter(x, y, c=data, cmap=cmap, vmin=vmin, vmax=vmax, marker="s", s=90, edgecolor=[0.0,0.0,0.0], linewidth=0.5, \
#            norm=MidpointNormalize(midpoint=0.0))
    cs = m.scatter(x, y, facecolor = 'red',marker="s", s=90, edgecolor=[0.0,0.0,0.0], linewidth=0.5)

    if plot_save==1:
       plt.savefig(fname, dpi=fig.dpi)

