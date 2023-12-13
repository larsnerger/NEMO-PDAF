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
    year_p = 2015
    lat_stat = 54.7 # DarsserS: 54.70, Ostergarnsholm: 57.42
    lon_stat = 12.7 # DarsserS: 12.70, Osstergarnsholm: 18.99
    depth_plot = 2 # 0 - surface, 2 - DarsserS CHL, 0 - Ostergarnsholm OXY or any depth at the specific station (check Station_info.xlsx for available depths)
    dist = 10           # Grid point distance to check around central location 
    plot_fcst = 1 # 0 - analysis, 1 - forecast
    plot_ana = 0   # if only assim1 is to be plotted (all others are out-commented), then one can additionally plot analysis to forecast
    data_path = '/scratch/usr/hbxsovli/SEAMLESS/validation_data/INSITU_BAL_PHYBGCWAV_DISCRETE_MYNRT_013_032'
    ens_range = [10, 1]
    z_mean=0     #0(no vertical mean)
    z_integral=0 #(0 for no vertical integral; 1 for vertical integral)
    z1=0         # upper Z boundary (for z_mean and z_integral)
    z2=5         # lower Z boundary (for z_mean and z_integral)

    domain = 'no'
    varnum = 24	 	# Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4, 
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2, 
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP
    plot_save = 1
    ######################################################

    # Define colormap 
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

    # 
    file_data = 'insitu_'+varstr+'_cmems_2015_feb-may.pkl'
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
  
    # MODEL
    grid_area = 'coarse'
    time_stamp = '00'
    DA_switch = 0
    coupled = 'weak'
    if depth_plot == 0: # surface
      depth_model = '0'
    else: # other depth
      depth_model = str(get_model_depth_idx(depth_plot))

    model_idx = get_model_idx(lat_stat, lon_stat, dist)
    free_final = np.array([])
    parm_final = np.array([])

    one_hour = timedelta(hours=1)
    model_idx = get_model_idx(lat_stat, lon_stat, dist)

    idx_lat = [i for i, la in enumerate(lat) if (abs(la - lat_stat) <= 0.1)]
    idx_lon = [i for i, lo in enumerate(lon) if (abs(lo - lon_stat) <= 0.1)]
    if depth_plot == 0: #surface 
      #idx_depth = [i for i, de in enumerate(depth) if (de <= 1 and de >= 0)]
      idx_depth = [i for i, de in enumerate(depth) if (de <= 1.5 and de >= 0)]
    else:
      idx_depth = [i for i, de in enumerate(depth) if abs(de-depth_plot)<0.001]
    idx_coord = np.intersect1d(idx_lat, idx_lon)
    idx_coord = np.intersect1d(idx_coord, idx_depth)
    start_date = min(time[idx_coord])
    end_date = max(time[idx_coord])
    date_array = [datetime(2015,1,1) + timedelta(days=i) for i in range((datetime(2015,5,31)-datetime(2015,1,1)).days + 1)]
    months_plot = list({dt.month for dt in date_array})

    print start_date, end_date
    
    for t in date_array:
      idx_t = [i for i, dt in enumerate(time) if abs(dt - t) <= one_hour]
      idx_final = np.intersect1d(idx_t, idx_coord)
      if len(idx_final)>0:
        data_final = np.append(data_final, np.mean(data[idx_final]))
        data_source_final.append([data_source[idx_final]])
        stat_coord_final.append(zip(lat[idx_final],lon[idx_final]))
      else:      
        data_final = np.append(data_final, np.nan)
        data_source_final.append([])

    DA_switch = 0

      #basedir = '/scratch/usr/hbxsovli/SEAMLESS/runfolders/'
    if plot_fcst == 1:
      DA_switch = 1 
    elif plot_fcst == 0: 
      DA_switch = 2

    data_source_final = [str(item) for sublist in data_source_final for array in sublist for item in array]

    uni_datasource = list(set(data_source_final)) 
    print 'Unique data source(s) = ', uni_datasource
    station_name = extractBetween(uni_datasource[0],'MO_','.nc')
    stat_coord_final = [item for sublist in stat_coord_final for item in sublist]
    uni_statcoord = unique_tuples(stat_coord_final)
    print 'Unique station lat/lon = ', uni_statcoord
 
    # Determine tick localizations and labels
    tickloc = []
    tickloc_min = []
    loc = date_array[0].timetuple().tm_yday
    plot_day = range(loc,date_array[-1].timetuple().tm_yday+1)
    tickloc.append(loc)
    months_str = []
    
    for i in range(len(months_plot)):
        _, _, mstr, days_in_mon = month_names(i+months_plot[0])
        loc = loc + days_in_mon
        tickloc.append(loc)
        tickloc_min.append(tickloc[i]+days_in_mon/2)
        months_str.append(mstr)


    # Calculate RMSE 
    
    if station_name == 'DarsserS':
      station_name = 'Darss Sill'
    # PLOTTING
    fname = varstr+'_'+station_name.replace(' ','_')+'_'+str(depth_plot)+'m.png'
    print 'fname', fname
    fig, ax = plt.subplots(figsize=(8,3))
    for e in ens_range:
      if e < 10:
	ensstr = '00'+str(e)
      else:
	ensstr = '0'+str(e)
      parm_final = np.array([])
      for t in date_array:
	if t.month < 5:
	  if t.month < 10:
	    mstr = '0'+str(t.month)
	  else:
	    mstr = str(t.month)
  	  if t.day < 10:	
	    dstr = '0'+str(t.day)
	  else:
	    dstr = str(t.day)
	  datestr = str(t.year)+mstr+dstr
          if plot_fcst == 1:
            DA_switch = 1
          elif plot_fcst == 0:
            DA_switch = 2
          if t.month == 1:
            parm_data = np.nan
          else:
            exp_path = '/scratch-emmy/projects/hbk00095/exp/free_scaleparams/data/'
            path = exp_path + ensstr + '_NORDIC_1d_ERGOM_T_' + datestr + '-' + datestr + '.nc'
            parm_data = read_latlon_day_path(varnum, grid_area, '2015', time_stamp, depth_model, DA_switch, coupled, z_mean, z_integral, z1, z2, t.month, \
                        t.day, model_idx, dist, path) # function
          parm_final = np.append(parm_final, parm_data)
      if e == 1:
        plt.plot(plot_day[0:len(parm_final)], parm_final, 'k', label='default')
      else:
        plt.plot(plot_day[0:len(parm_final)], parm_final, label=str(e))
    plt.plot(plot_day[0:len(parm_final)], data_final[0:len(parm_final)], 'm.', label='insitu obs')
    plt.ylabel(Variable+' ('+var_unit+')')
    print 'TEST tickloc', tickloc[0], tickloc[-1]
    plt.xlim(tickloc[0],tickloc[-1])
    ax.xaxis.set_major_locator(ticker.FixedLocator(tickloc))
    ax.xaxis.set_minor_locator(ticker.FixedLocator(tickloc_min))
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.tick_params(axis='x', which='minor', tick1On=False, tick2On=False)
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(months_str))   
    plt.legend()
    if depth_plot == 0:
      plt.title(Variable+' at '+station_name+' (surface)')
    else: 
      plt.title(Variable+' at '+station_name+' (depth='+str(depth_plot)+'m)')
    #  plt.title(Variable+' at depth='+str(depth_plot)+'m', fontsize=16)
     
    
    if plot_save == 1: 
      plt.savefig(fname, dpi=300)
    #plt.show()
