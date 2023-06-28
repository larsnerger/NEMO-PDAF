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


if __name__ == "__main__":
    coupled = ''
    figure_name = 'insitu_cmems_all.png'
    year_p = 2015
    month_p = [2,3,4,5]
    lat_stat = 57.5
    lon_stat = 17.5
    dist = 10           # Grid point distance to check around central location 
    plot_free = 1 
    plot_fcst = 1 
    plot_ana = 1   
    data_path = '/scratch/usr/hbxsovli/SEAMLESS/validation_data/INSITU_BAL_PHYBGCWAV_DISCRETE_MYNRT_013_032'
    z_mean=0     #0(no vertical mean)
    z_integral=0 #(0 for no vertical integral; 1 for vertical integral)
    z1=0         # upper Z boundary (for z_mean and z_integral)
    z2=5         # lower Z boundary (for z_mean and z_integral)

    domain = 'no'
    varnum = 2	 	# Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP
    plot_save = 1
    ######################################################

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
        strcmap = 'viridis'
#        strcmap = cmocean.cm.thermal
    elif varnum==25:   # TE
        strcmap = cmocean.cm.dense
    elif varnum==26:   # PFT
        strcmap = cmocean.cm.rain
    else:
        strcmap = 'gist_ncar_r'


    varstr, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()

    file_data = 'insitu_'+varstr+'_cmems_2015.pkl'
    print 'Read ', file_data
    with open(file_data, 'rb') as f:
      dlist = pickle.load(f) 
      print 'Done reading'

    file_coords = 'unique_coords_'+varstr+'_cmems.pkl'
    print 'Read ', file_coords
    with open(file_coords, 'rb') as f:
      uni_coords = pickle.load(f)
      print 'Done reading'


    time = dlist['time']
    lat = dlist['lat'] 
    lon = dlist['lon']
    depth = dlist['depth']
    data = dlist['data']
    data_qc = dlist['data_qc']
   
    data_final = []
    data_qc_final = []
    time_final = []
    depth_final = []
    
    print 'Original data length = ', len(data)

    print '-----------------------------------'
    print 'Unique coords = ', len(uni_coords)
    
    month_mask = np.array([date.month in month_p for date in time])
    depth_mask = np.array([d<1 for d in depth]) 
    combined_mask = month_mask & depth_mask
    lat_sel = lat[combined_mask]
    lon_sel = lon[combined_mask]
    time_sel = time[combined_mask]
    depth_sel = depth[combined_mask]
    data_sel = data[combined_mask]

    print 'Length after time and depth selection = ', len(data_sel)
    
    uni_coords_sel, uni_idx = unique_tuples(zip(lat_sel,lon_sel),0.1)
    print uni_idx
    print uni_coords_sel
    print 'Unique coords in ', month_p ,' = ', len(uni_coords_sel)

    lat_sel = lat_sel[uni_idx]
    lon_sel = lon_sel[uni_idx]
    data_sel = data_sel[uni_idx]
    
    
    print '#Final data points', len(lat_sel)

    if domain=='no':
        llon=-4.0417
        llat=48.575
        ulon=30.375
        ulat=65.875
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

    m = Basemap(projection='merc', \
            llcrnrlat=50, urcrnrlat=65, \
            llcrnrlon=0, urcrnrlon=30, \
            lat_ts=20, \
            resolution=bmres, area_thresh=10)
    # Draw coastlines, countries, and states
    m.drawcoastlines()
    m.drawcountries()
    m.drawlsmask(land_color='Linen')       # No geographic detail
    #m.shadedrelief()                        # Low geographic detail
    #m.bluemarble()                         # High geographic detail

    # draw parallels.
    parallels = np.arange(0.,360.,res)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14)
    # draw meridians
    meridians = np.arange(0.,360.,res)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14)

    print 'min max colorbar=', min(data_sel), max(data_sel)
    #clev = np.linspace(min(data_temp), max(data_temp), 100)

# Convert latitudes and longitudes to map coordinates
    x, y = m(lon_sel, lat_sel)

# Scatter plot the data on the basemap
    m.scatter(x, y, c=data_sel, cmap=strcmap, vmin=min(data_sel), vmax=max(data_sel))

# Add a colorbar
    cbar = plt.colorbar()
    cbar.set_label('Data')

    plt.show()

    if plot_save == 1:
      plt.savefig(figure_name)
 
