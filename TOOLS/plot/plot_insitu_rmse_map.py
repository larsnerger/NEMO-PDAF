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
import xlrd
import re
from collections import OrderedDict
import local_settings

if __name__ == "__main__":
    t_start = datetime(2015,3,1)
    t_end = datetime(2015,5, 31)
    coupled = ''
    strcmap = 'cmo.amp' # colormap, temp-cmo.solar, 
    depth_plot = 0 # 0-surface (only possible option atm) 
    dist = 10           # Grid point distance to check around central location 
    plot_ffa = 0 # 0 - free, 1 - forecast, 2 - analysis
    z_mean=0     # 0 - no vertical mean
    z_integral=0 #(0 for no vertical integral; 1 for vertical integral)
    z1=0         # upper Z boundary (for z_mean and z_integral)
    z2=5         # lower Z boundary (for z_mean and z_integral)

    domain = 'ba'
    varnum = 2	 	# Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP
    plot_save = 0
    plot_show = 0 
    plot_numobs = 0 # plot the scatter plot with number of observations
    minmax = [0, 1.1] # set equal min/max for automatic minmax
    ######################################################

    # get name of run 
    if plot_ffa == 0: asml = 0
    else: asml = 1
    full_path = get_exp_path(coupled, '2015', '2', '01', '00', asml)
    model_run = extractBetween(full_path, 'exp/', '/DA/') 
    print 'Model run = ', model_run

    # create datetime array for entire time span specified by t_start - t_end
    num_days = (t_end - t_start).days + 1
    date_array = [t_start + timedelta(days=i) for i in range(num_days)]

    varstr, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()

    # Read .xlsx: all stations with unique averaged coords
    xlsx_path = '/scratch/usr/hbxsovli/SEAMLESS/validation_data/Station_info.xlsx'
    wb = xlrd.open_workbook(xlsx_path)
    sheet = wb.sheet_by_index(0)
    num_rows = sheet.nrows
    column_values = [sheet.col_values(col, start_rowx=0, end_rowx=num_rows) for col in range(3)]
    stat_name = column_values[0][1:]
    stat_name = [str(element) for element in stat_name]
    stat_name = [s.rstrip('.0') for s in stat_name]
    stat_lat = column_values[1][1:]
    stat_lon = column_values[2][1:]

    # Read data measured by insitu stations
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
    data_source = dlist['data_source']

    data_final = np.array([])
    data_qc_final = np.array([])
    time_final = np.array([])
    depth_final = np.array([])
    data_source_final = []
    stat_coord_final = []

    # select only the data from the mooring stations 
    substring = '_MO_'
    mo_mask = [substring in element for element in data_source] 
    time = time[mo_mask]
    lat = lat[mo_mask]
    lon = lon[mo_mask]
    depth = depth[mo_mask]
    data = data[mo_mask]
    data_qc = data_qc[mo_mask]
    data_source = data_source[mo_mask] 

    stat_name_final = []
    lat_final = []
    lon_final = []
    rmse_final = []
    numobs_final = []
    uni_stations = list(set(data_source))

    grid_area = 'coarse'
    time_stamp = '00'
    if plot_ffa == 0: #free 
      DA_switch = 0
    elif plot_ffa == 1: #forecast
      DA_switch = 1
    elif plot_ffa == 2: #analysis
      DA_switch = 2

    coupled = 'weak'
    if depth_plot == 0: # surface
      depth_model = '0'
    elif depth_plot == 1: # bottom
      depth_model = '20' #~40m
    
    # CALC RMSE
    # loop over each station   
    for s in range(0,len(uni_stations)):
      # OBSERVATIONS
      print uni_stations[s]
      #station_type = extractBefore(uni_stations[s], '/') 
      uni_stat_short = extractBetween(uni_stations[s],'_MO_','.nc')
      #uni_stat_short = extractBetween(uni_stations[s], station_type+"_", ".nc")
      idx = [index for index, string in enumerate(stat_name) if uni_stat_short==string]
      print 'idx=', idx
      lon_mask = abs(lon-stat_lon[idx[0]])<0.2
      lat_mask = abs(lat - stat_lat[idx[0]])<0.2
      model_idx = get_model_idx(stat_lat[idx[0]], stat_lon[idx[0]],dist)
      if depth_plot == 0: #surface
        depth_mask = depth<=2
      coord_mask = lat_mask & lon_mask & depth_mask
      time_mask = [t_start <= date <= t_end for date in time]
      #station_mask = [index for index, string in enumerate(data_station) if uni_stat_short[0]==string]
      #final_mask = coord_mask & station_mask 
      final_mask = coord_mask & time_mask
      if sum(final_mask)>0: 
        stat_name_final.append(uni_stat_short[0])
        lat_final.append(stat_lat[idx[0]])
        lon_final.append(stat_lon[idx[0]])
        # loop over all time steps and average 
        time_sel = time[final_mask]
        data_sel = data[final_mask]
        uni_time = list(OrderedDict.fromkeys(map(lambda dt: dt.date(), time_sel)))
        obs_temp = []
        model_temp = []
        numobs_final.append(len(uni_time))
        for t in range(0,len(uni_time)):   
          t_mask = [dt.date() == uni_time[t] for dt in time_sel]
          obs_temp.append(np.mean(data_sel[t_mask]))
          # MODEL 
          model_val = read_latlon_day(varnum, grid_area,uni_time[t]. year, time_stamp, depth_model, DA_switch, coupled, z_mean, \
		z_integral, z1, z2, uni_time[t].month, uni_time[t].day, model_idx, dist)
          model_temp.append(model_val)
        if check_list_type(model_temp) == 2:  
	  # if model_temp is list of lists, convert to list 
	  model_temp = [item for sublist in model_temp for item in sublist]
	# RMSE 
        #print 'OBS_TEMP = ', obs_temp
        #print 'MODEL_TEMP = ', model_temp
	rmse_final.append(calc_rmse(obs_temp, model_temp))
 	#rmse_final.append(np.mean(obs_temp))

    print 'rmse_final = ', size(rmse_final) 
    print 'numobs_final = ', size(numobs_final)
    print 'lat_final = ', size(lat_final)
    print 'lon_final = ', size(lon_final)
    # PLOTTING
    print 'MEAN RMSE = ', np.mean(rmse_final)
    # RMSE
    print '\n Test RMSE_final = ', rmse_final
    fname = varstr+'_insitu_rmse_map_'+model_run+'.png'
    cbar_label = 'RMSE ['+var_unit+']'
    title_str = 'RMSE for '+Variable+' ('+model_run+')'
    plot_scatter_map(lat_final, lon_final, np.array(rmse_final), strcmap, domain, minmax, title_str, cbar_label, plot_show, plot_save, fname)
    # number of observations
    fname = varstr+'_insitu_numobs_map_'+model_run+'.png'
    cbar_label = 'Number of obs'
    title_str = 'Number of observations for '+Variable+' RMSE calculation'
    if plot_numobs==1:
      minmax = [0,0]  
      strcmap = 'bone_r'
      plot_scatter_map(lat_final, lon_final, np.array(numobs_final), strcmap, domain, minmax, title_str, cbar_label, plot_show, plot_save, fname)
    
