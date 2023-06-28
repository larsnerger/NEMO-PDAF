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
from matplotlib.patches import Patch
import cmocean
import pickle

if __name__ == "__main__":
    name_save = 'Test_all_TEM_stations_surface_April_2hours.png'
    coupled = ''
    year_p = 2015
    month_p = [4]
    #day_p = range(1,29)           
    dist = 10           # Grid point distance to check around central location 
    plotlog = 1 
    plotcb = 1          # Whether to show the colorbar
    strcmap = 'viridis'   
    path_data1 = '/scratch/usr/hbxsovli/SEAMLESS/validation_data/ICES_2015_Bottle_LR_CTD/ICES_2015_Bottle_LR_CTD.xlsx' # path to validation data .xlsx file 
    path_data2 = '/scratch/usr/hbxsovli/SEAMLESS/validation_data/ICES_2015_HR_CTD/ICES_2015_HR_CTD.xlsx'
    dist = 10 # distance to station  
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
    var_plot = 3        # To plot data or data information: 1 = data (varnum value), 2 = month, 3 = day, 4 = hour
    save_plot = 1
 
    ######################################################

    varstr, mat_var, Variable, var_unit = var_names(varnum)
    MAT_VAR = mat_var.upper()
    
    varnum_data2 = [2, 3, 21] # variables contained in HR_CTD data set
    
    hour_p = [0,1, 2, 3, 21, 22, 23, 24] # valid hours 
   # hour_p = range(0,25)
    if 'day_p' not in locals():
      day_p = range(1,32)

    data_final1 = []
    lat_final1 = []
    lon_final1 = []
    depth_final1 = []
    day_final1 = []
    month_final1 = []
    year_final1 = []
    hour_final1 = []

    data_final2 = []
    year_final2 = []
    lat_final2 = []
    lon_final2 = []
    depth_final2 = []
    day_final2 = []
    month_final2 = []
    hour_final2 = []

    # read data 
    year1, month1, day1, hour1, minute1, lon1, lat1, depth1, data1, data_qv1 = read_insitu_vali_ICES(varnum, path_data1)
    print 'TEST ', len(year1)     
    # save insitu data in file 
    data_dict1 = {'year':year1, 'month':month1, 'day':day1, 'hour':hour1, 'minute':minute1, 'lon':lon1, 'lat':lat1, 'depth':depth1, 'data':data1, 'data_qv':data_qv1}
    file_data1 = 'ICES_LR_'+varstr+'_2015.pkl'

    print 'Data source 1 total length = ', len(data1) 
    if len(year1)>0:
      print 'Source 1:'
      for r in range(len(data1)):
	if depth1[r]<3 and month1[r] in month_p and hour1[r]!=str and day1[r]>0 and hour1[r] in hour_p:
        #if abs(lat1[r]-54.883)<0.1 and abs(lon1[r]-13.866)<0.1 and month1[r]<6 and depth1[r]<3:
        #if (day1[r] in day_p) and type(data1[r]) != str and year1[r]==year_p and (month1[r] in month_p) and (hour1[r] in hour_p) and \
#			(abs(lat1[r]-54.883)<0.1 and abs(lon1[r]-13.866)<0.1):
	  if var_plot == 1: 
	    data_final1.append(data1[r])
	  elif var_plot == 2: 
	    data_final1.append(month1[r])
	  elif var_plot == 3: 
	    data_final1.append(day1[r])
	  elif var_plot == 4: 
	    data_final1.append(hour1[r])
	  else: 
	    print 'ERROR: need to specify a var_plot between 1-4'
          lat_final1.append(lat1[r])
          lon_final1.append(lon1[r])
          depth_final1.append(depth1[r])
          day_final1.append(day1[r])
          month_final1.append(month1[r])
          hour_final1.append(hour1[r])
          year_final1.append(year1[r])

    len_data1 = len(data_final1)
    print 'Data source 1 after selection = ', len_data1 
   
    # get unique lat/lons
    unique_stat1, unique_idx1 = unique_tuples(zip(lat_final1, lon_final1))
    print '--------- All available stations: DATA SET 1 lat/lons ------------'
    print unique_stat1

    #print '+++++++++++++ All time steps ++++++++++++++++'
    #print zip(year_final1, month_final1, day_final1, hour_final1)

    if varnum in varnum_data2:
      year2, month2, day2, hour2, minute2, lon2, lat2, depth2, data2, data_qv2 = read_insitu_vali_ICES(varnum, path_data2)
      # Save excel data in external .pkl file
      data_dict2 = {'year': year2, 'month':month2, 'day': day2, 'hour':hour2, 'minute':minute2, 'lon':lon2, 'lat':lat2, 'depth':depth2, 'data':data2, 'data_qv':data_qv2}
      file_data2 = 'ICES_HR_'+varstr+'_2015.pkl'

      print 'Data source 2 total length = ', len(data2) 
      if len(year2) > 0:
        print 'Source 2:'
        for r in range(len(data2)):
	  if depth2[r]<3 and month2[r] in month_p and day2[r]>0 and hour2[r]!=str and hour2[r] in hour_p:
          #if abs(lat2[r]-54.883)<0.1 and abs(lon2[r]-13.866)<0.1 and month2[r]<6 and depth2[r]<3:
          #if (day2[r] in day_p) and type(data2[r]) != str and year2[r]==year_p and (month2[r] in month_p) and (hour2[r] in hour_p) and depth2[r]<=1:
	    if var_plot == 1: 
	      data_final1.append(data2[r])
	    elif var_plot == 2: 
	      data_final1.append(month2[r])
	    elif var_plot == 3: 
	      data_final1.append(day2[r])
	    elif var_plot == 4: 
	      data_final1.append(hour2[r])
	    else: 
              print 'ERROR: need to specify a var_plot between 1-4'
	    lat_final1.append(lat2[r])
            lon_final1.append(lon2[r])
            depth_final1.append(depth2[r])
            day_final1.append(day2[r])
            month_final1.append(month2[r])
	    hour_final1.append(hour2[r])
            year_final1.append(year2[r])

      #len_data2 = len(data_final2)
      #'Data source 2 after selection length = ', len_data2
      #minmax = [min(min(data_final1), min(data_final2)), max(max(data_final1), max(data_final2))]
      #print 'Length of final data after selection = ', len_data1+len_data2
      #unique_stat2, unique_idx2 = unique_tuples(zip(lat_final2, lon_final2), 0.01)
      #print '--------- All available stations: DATA SET 2 lat/lons ------------'
      #print unique_stat2

      #print '+++++++++++++ All time steps ++++++++++++++++'
      #print zip(year_final2, month_final2, day_final2, hour_final2)
    else: 
      minmax = [min(data_final1), max(data_final1)]
      print 'Length of final data after selection = ', len_data1

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

    #print 'min max colorbar=', minmax
    #clev = np.linspace(min(data_temp), max(data_temp), 100)

    # Convert latitudes and longitudes to map coordinates
    print 'test lon1 lat1 = ', len(lon_final1), len(lat_final1) 
    print 'data = ', data_final1
    print 'len data', len(data_final1)
    x1, y1 = m(lon_final1, lat_final1)
    m.scatter(x1, y1, c=data_final1, cmap=strcmap, marker='s')#, vmin=min(data_final1), vmax=max(data_final1))
    
    if len(lon_final2)>0:
      x2, y2 = m(lon_final2, lat_final2)
      m.scatter(x2, y2, c=data_final2, cmap=strcmap, marker='o', vmin=min(data_final2), vmax=max(data_final2))
# Scatter plot the data on the basemap

    # add the legend to the plot a colorbar
    cbar = plt.colorbar()
    if var_plot == 1: 
      plt.title('Data for '+Variable)
      cbar.set_label('Data')
    elif var_plot == 2: 
      plt.title('Month '+ str(month_p) +'for '+Variable)
      cbar.set_label('Month')
    elif var_plot == 3: 
      plt.title('Day in month '+str(month_p[0])+ ' for '+Variable)
      cbar.set_label('Day')
    elif var_plot == 4: 
      plt.title('Hour on day '+str(day_p)+' in month '+str(month_p)+' for '+Variable)
      cbar.set_label('Hour')

    print 'Data final 1'
    print data_final1
    #print 'Data final 2'
    #print data_final2
# Show the basemap
    #plot_map(np.array(data_temp), np.array(lat_temp), np.array(lon_temp), varnum, 'no', 'jet', minmax, 0, 'Test', 'test.png', 0, 1)
    if save_plot == 1:
      plt.savefig(name_save)

    plt.show()
