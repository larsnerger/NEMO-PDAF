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
    save_pkl = 1
    varnum = 2	 	# Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP
    ######################################################


    varstr, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()

    file_data = 'insitu_'+varstr+'_cmems_2015.pkl'
    print 'Read ', file_data
    with open(file_data, 'rb') as f:
      dlist = pickle.load(f) 
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
    
    print 'Original data len'
    print 'time', len(time)
    print 'lat', len(lat)
    print 'lon', len(lon) 
    print 'data', len(data)
    print 'depth', len(depth)
  
    # get unique coordinates
    print 'Start finding unique coordinates'
    uni_coord, uni_coord_idx = unique_tuples(zip(lat, lon), 0.05) 
    print 'Finished'
    print 'lat = ',[item[0] for item in uni_coord]
    print 'lon = ', [item[1] for item in uni_coord]
   
    # MAP COORDINATES TO MODEL GRID
    #print len(uni_coord) 	
    #for i in range(0,len(uni_coord)):
    #  model_idx = get_model_idx(uni_coord[i][0], uni_coord[i][1], 1)
    #  print 'Length of model_idx = ', len(model_idx)

    if save_pkl == 1: 
      file_coords = 'unique_coords_'+varstr+'_cmems.pkl'
      print 'Save to', file_coords
      with open(file_coords, 'wb') as f:
        dlist = pickle.dump(uni_coord,f)
        print 'Done saving'

