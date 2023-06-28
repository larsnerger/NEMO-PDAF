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
    figure_name = 'insitu_cmems_ts.png'
    year_p = 2015
    lat_stat = 54.883
    lon_stat = 13.8660
    depth_plot = 1 # 0-surface, 1-bottom
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
    #####################################################

    file_data = 'insitu_DOX1_cmems_2015.pkl'
    print 'Read ', file_data
    with open(file_data, 'rb') as f:
      dlist_oxy1 = pickle.load(f)
      print 'Done reading'


    file_data = 'insitu_DOXY_cmems_2015.pkl'
    print 'Read ', file_data
    with open(file_data, 'rb') as f:
      dlist_doxy = pickle.load(f)
      print 'Done reading'

    keys_list = dlist_doxy.keys()
    d_merged = {}

    for key in keys_list:
      if key in dlist_oxy1 and key in dlist_doxy:
        d_merged[key] = np.concatenate((dlist_oxy1[key], dlist_doxy[key]))
      elif key in dlist_oxy1:
        d_merged[key] = dlist_oxy1[key]
      elif key in dlist_doxy:
        d_merged[key] = dlist_doxy[key]

    #all_depths = d_merged['depth']
    #dsources = d_merged['data_source']
    #unique_depths, indices = np.unique(all_depths, return_index = True)
    #print 'unique_depths=', unique_depths
    #for i in range(0, len(indices)):
    #  print  'dsources[indices]', dsources[indices[i]]

    file_data = 'insitu_OXY_cmems_2015.pkl'
    with open(file_data, 'wb') as f:
      pickle.dump(d_merged,f)
      print 'Done saving'

