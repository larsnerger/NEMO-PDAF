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
from datetime import datetime


if __name__ == "__main__":
    coupled = ''
    path_data1 = '/scratch/usr/hbxsovli/SEAMLESS/validation_data/ICES_2015_Bottle_LR_CTD/ICES_2015_Bottle_LR_CTD.xlsx' # path to validation data .xlsx file
    path_data2 = '/scratch/usr/hbxsovli/SEAMLESS/validation_data/ICES_2015_HR_CTD/ICES_2015_HR_CTD.xlsx'
    save = 1
    varnum = 2	 	# Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP
    map2model = 0 
    ######################################################

    varstr, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()

    varnum_data2 = [2, 3, 21] # variables contained in HR_CTD data set
    # READ DATA1 
    year1, month1, day1, hour1, minute1, lon1, lat1, depth1, data1, data_qv1 = read_insitu_vali_ICES(varnum, path_data1)
    
    time1 = [datetime(year, month, day) for year, month, day in zip(year1, month1, day1)]
    #time1 = []
    #for i in range(0,len(year1)):
    #  print i,'/',len(range(0,len(year1)))
    #  time1.append(datetime(year1[i], month1[i], day1[i], hour1[i], minute1[i]))
    #   time1.append(datetime(year1[i], month1[i], day1[i])

    # save insitu data in file
    #data_dict1 = {'year':year1, 'month':month1, 'day':day1, 'hour':hour1, 'minute':minute1, 'lon':lon1, 'lat':lat1, 'depth':depth1, 'data':data1, 'data_qv':data_qv1}
    data_dict1 = {'time':time1, 'lon':lon1, 'lat':lat1, 'depth':depth1, 'data':data1, 'data_qv':data_qv1}
   
    # MAP TO MODEL COORDINATES
    # read model grid to check that correct coordinates are chosen 
    if map2model == 1: 
      model_file = '/scratch/projects/hbk00095/exp/free_N30/DA/state_20150201.nc'
      ncfile = NetCDFFile(model_file) 
      mlon = ncfile.variables['nav_lon'][:]
      mlat = ncfile.variables['nav_lat'][:]
      lat_vec = mlat[:,0]
      lon_vec = mlon[0,:]

      empty_modelgrid = [[[] for _ in range(len(mlon[0]))] for _ in range(len(mlon))]
      model_lats = []
      model_lons = []

      print 'Map to model coordinates'
      for s in range(0,len(time1)):
        print s,'/',len(range(0,len(time1)))
        for i in range(len(mlat[:,0])):
          lat_vec[i] = np.max(mlat[i,:])
        for j in range(len(mlat[0,:])):
          lon_vec[j] = np.max(mlon[:,j])
        lon_vec[-1] = lon_vec[-2] + lon_vec[1]-lon_vec[0]

        latdiff = lat_vec - lat1[s]
        londiff = lon_vec - lon1[s]
        model_idx = []
        model_idx.append(np.argmin(abs(latdiff)))
        model_idx.append(np.argmin(abs(londiff)))
        empty_modelgrid[model_idx[0]][model_idx[1]].append(s)
        model_lats.append(model_idx[0])
        model_lons.append(model_idx[1])

      print 'Find unique model coordinates'
      uni_modelcoords, _ = unique_tuples(zip(model_lats, model_lons))
      print uni_modelcoords

    file_data1 = 'ICES_LR_'+varstr+'_2015.pkl'
    print 'Save ', file_data1
    with open(file_data1, 'wb') as f:
      pickle.dump(data_dict1, f)
      print 'Done saving data_dict1'

    
   
    if varnum in varnum_data2:
      year2, month2, day2, hour2, minute2, lon2, lat2, depth2, data2, data_qv2 = read_insitu_vali_ICES(varnum, path_data2)
      # Save excel data in external .pkl file
      time2 = [datetime(year, month, day) for year, month, day in zip(year2, month2, day2)]
      data_dict2 = {'time':time2, 'lon':lon2, 'lat':lat2, 'depth':depth2, 'data':data2, 'data_qv':data_qv2}
      file_data2 = 'ICES_HR_'+varstr+'_2015.pkl'
      with open(file_data2, 'wb') as f:
        pickle.dump(data_dict2, f)
        print 'Done saving data_dict2'

