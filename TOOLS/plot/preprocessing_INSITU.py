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
    data_path = '/scratch/usr/hbxsovli/SEAMLESS/validation_data/INSITU_BAL_PHYBGCWAV_DISCRETE_MYNRT_013_032'
    save = 1
    varnum = 21	 	# Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP
    ######################################################

    varstr, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()

    # get cmems insitu data in months 2,3,4,5 of 2015
    time, lon, lat, depth, data, data_qc, data_source = read_insitu_vali_CMEMS(varnum, data_path)    
    
    print 'Done reading all insitu data'
    print 'Length = ', len(data)
    print 'Save time, lon, lat, depth, data, data_qc, data_source'

    data_dict = {'time': time, 'lon': lon, 'lat':lat, 'depth':depth, 'data':data, 'data_qc':data_qc, 'data_source':data_source}

    print 'Original data len'
    print 'time', len(time)
    print 'lat', len(lat)
    print 'lon', len(lon)
    print 'data', len(data)
    print 'data_qc', len(data_qc)
    print 'depth', len(depth)
    print 'data_source', len(data_source)


    file_data = 'insitu_'+varstr+'_cmems_2015.pkl'
    file_data = 'insitu_DOX1_cmems_2015.pkl'
    if save == 1: 
      with open(file_data, 'wb') as f:
        pickle.dump(data_dict,f)
        print 'Done saving'


