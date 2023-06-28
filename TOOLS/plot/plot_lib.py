#!/usr/bin/python

import sys, getopt
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
import numpy as np
from mpl_toolkits.basemap import Basemap
import time
import scipy.io as spio
from netCDF4 import Dataset as NetCDFFile
from local_settings import *
from scipy import stats
import xlrd
import numpy.ma as ma
import os 
import datetime

def stations(station):

    print 'Read coordinates for station:', station

    # Open Excel file
    wb = xlrd.open_workbook("/home/hzfblner/SEAMLESS/plot/InfoWas_insitu_Stationen.xlsx")
    sheet = wb.sheet_by_index(0)

    # Read index ranges for Baltic and North Sea from Excel file
    bs_start = int(sheet.cell_value(1, 3))
    bs_end = int(sheet.cell_value(1, 4))
    ns_start = int(sheet.cell_value(2, 3))
    ns_end = int(sheet.cell_value(2, 4))
    print 'Ranges', bs_start, bs_end, ns_start, ns_end

    # Read coordinates
    for i in range(bs_start-1, bs_end):
        if sheet.cell_value(i, 1)==station:
            lon = sheet.cell_value(i, 3)
            lat = sheet.cell_value(i, 2)
            noba = 'ba'

    for i in range(ns_start-1, ns_end):
        if sheet.cell_value(i, 1)==station:
            lon = sheet.cell_value(i, 3)
            lat = sheet.cell_value(i, 2)
            noba = 'no'
    
    # DarsserS is not in .xlsx file above
    if station == 'DarsserS':
      lat = 54.70 
      lon = 12.70 
      noba = 'ba'
            
    return lon, lat, noba

def station_info(istation, mod_or_obs):

    station_idx = []
    station_coords = []

    if istation == 1:
        station = 'Arkona'
        station_coords.append(54.88)
        station_coords.append(13.87)
        noba = 'ba'
    elif istation == 2:
        station = 'Baltic_Proper'
        station_coords.append(56.0)
        station_coords.append(18.5)
        noba = 'ba'
    elif istation == 3:
        station= 'Bothnian_Sea'
        station_coords.append(61.7)
        station_coords.append(19.5)
        noba = 'ba'
    elif istation == 4:
        station= 'FINO1WR'
        station_coords.append(54.00)
        station_coords.append(6.58)
        noba = 'no'
    else:
        station= 'Sleipner-A'
        station_coords.append(58.37)
        station_coords.append(1.91)
        noba = 'no'
                    
    return station, station_coords, noba

def month_names(mon):
    # Define Month strings and number of days of month
    month_nam  = ['jan','feb','mar','apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct','nov','dec']
    month_NAM  = ['Jan.','Feb.','Mar.','Apr.', 'May', 'Jun.', 'Jul.', 'Aug.', 'Sep.', 'Oct.','Nov.','Dec.']
    month_Nfull  = ['January','February','March','April', 'May', 'June', 'July', 'August', \
                    'Sept.', 'Oct.','Nov.','Dec.']
#                    'September', 'October','November','December']
    month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    return month_nam[mon-1], month_NAM[mon-1], month_Nfull[mon-1], month_days[mon-1]

def day_of_year(month, day):
    # Set day of year

    DoY = 1
    for i in range(month-1):
        _, _, _, days_in_mon = month_names(i+1)
        DoY = DoY + days_in_mon

    return DoY

def var_names(X):
    # A: variable name in file
    # B: variable name in in situ data file
    # C: Full variable string
    # unit: Unit of variable
    if X == 1:
        A = 'SSH_inst'
        B = 'z'
        C = 'SSH'
        unit = 'm'
    elif X == 2:
        A = 'votemper'
        B = 'tem'
        C = 'Temperature'
        unit = '$^o$C'
    elif X == 3:
        A = 'vosaline'
        B = 'sal'
        C = 'Salinity'
        unit = 'PSU'
    elif X == 4:
        A = 'uos'
        B = 'uvel'
        C = 'uvel'
        unit = 'm/s'
    elif X == 5:
        A = 'vos'
        B = 'vvel'
        C = 'vvel'
        unit = 'm/s'
    elif X == 6:
        A = 'NH4'
        B = 'amm'
        C = 'Ammonium'
        unit = 'mmol N/m$^3$'
    elif X == 7:
        A = 'NO3'
        B = 'nit'
        C = 'Nitrate'
        unit = 'mmol N/m$^3$'
    elif X == 8:
        A = 'PO4'
        B = 'pho'
        C = 'Phosphate'
        unit = 'mmol P/m$^3$'
    elif X == 9:
        A = 'SIL'
        B = 'sil'
        C = 'Silicate'
        unit = 'mmol Si/m$^3$'
    elif X == 10:
        A = 'DIA'
        B = 'DIA'
        C = 'Diatoms'
        unit = 'mmol N/m$^3$'
    elif X == 11:
        A = 'FLA'
        B = 'FLA'
        C = 'Flagellates'
        unit = 'mmol F/m$^3$'
    elif X == 12:
        A = 'CYA'
        B = 'CYA'
        C = 'Cyanobacteria'
        unit = 'mmol C/m$^3$'
    elif X == 13:
        A = 'MEZ'
        B = 'MEZ'
        C = 'Meso zooplankton'
        unit = 'mmol M/m$^3$'
    elif X == 14:
        A = 'MIZ'
        B = 'MIZ'
        C = 'Micro zooplankton'
        unit = 'mmol M/m$^3$'
    elif X == 15:
        A = 'DET'
        B = 'POC'
        C = 'Particulate Organic Carbon'
        unit = 'mmol/m$^3$'
    elif X == 16:
        A = 'DETs'
        B = 'DETs'
        C = 'Detritus Silica'
        unit = 'mmol/m$^3$'
    elif X == 17:
        A = 'FE'
        B = 'FE'
        C = 'Iron'
        unit = 'mmol/m$^3$'
    elif X == 18:
        A = 'LDON'
        B = 'LDON'
        C = 'lDON'
        unit = 'mmol/m$^3$'
    elif X == 19:
        A = 'DIC'
        B = 'DIC'
        C = 'Dissolved Inorganic Carbon'
        unit = 'mmol/m$^3$'
    elif X == 20:
        A = 'ALK'
        B = 'ALK'
        C = 'Total Alkalinity'
        unit = 'mmol/m$^3$'
    elif X == 21:
        A = 'OXY'
        B = 'OXY'
        C = 'Dissolved Oxygen'
        unit = 'mmol/m$^3$'
    elif X == 22:
        #A = 'xpco2'
	A = 'PCO2'
        B = 'PCO2'
        C = 'Partial Pressure of Carbon Dioxide'
        unit = 'micro atm'
    elif X == 23:
        #A = 'xph'
	A = 'PH'
        B = 'pH'
        C = 'pH'
        unit = ''
    elif X == 24:
        A = 'CHL'
        B = 'chl'
        C = 'Chlorophyll'
        unit = 'mg Chl/m$^3$'
    elif X == 25:
        A = 'TE'
        B = 'TE'
        C = 'Trophic Efficiency'
        unit = 'ratio'
    elif X == 26:
        A = 'PFT'
        B = 'PFT'
        C = 'PFT Ratio'
        unit = 'ratio'
    elif X == 27:
        A = 'PP'
        B = 'PP'
        C = 'Primary Production'
        unit = 'mmolC/m$^2$/d'

    return A, B, C, unit

class MidpointNormalize(colors.Normalize):
        ## This class normalises the colorbar so white is at zero
        def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
                self.midpoint = midpoint
                colors.Normalize.__init__(self, vmin, vmax, clip)

        def __call__(self, value, clip=None):
                x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
                return np.ma.masked_array(np.interp(value, x, y))

def loadmat(filename):
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict        

def _todict(matobj):
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

def read_mat(myfile, length, depth, variable):
    # Reads .mat files	
    # returns values and lon/lat as indices
    mat_aux = loadmat(myfile)

    Value=[]
    time=[]
    dep=[]
    lat=[]
    lon=[]

    mat = mat_aux['RefData']['data']['variable']
    if variable!='secchi_depth_of_sea_water':
        if variable!='mole_concentration_of_oxygen_in_seawater':
            for lo in range(llon,ulon): #=1:412:
                for la in range(347): #1:347:
                    for d in range(depth+1): #1:1%25:
                        for t in range(length):
                            if mat[t,d,la,lo] > -100 and mat[t,d,la,lo] < 1e6:
                                Value.append(mat[t,d,la,lo])
                                time.append(t)
                                dep.append(d)
                                lat.append(la)
                                lon.append(lo)
        else:
            # Exclude outliers for oxygen
            for lo in range(llon,ulon): #=1:412:
                for la in range(347): #1:347:
                    for d in range(depth+1): #1:1%25:
                        for t in range(length):
                            if mat[t,d,la,lo] > 240 and mat[t,d,la,lo] < 350:
                                Value.append(mat[t,d,la,lo])
                                time.append(t)
                                dep.append(d)
                                lat.append(la)
                                lon.append(lo)
    else:
        d=0
        for lo in range(412): #=1:412:
            for la in range(347): #1:347:
                for t in range(length):
                    if mat[t,la,lo] > -100 and mat[t,la,lo] < 1e6:
                        Value.append(mat[t,la,lo])
                        time.append(t)
                        dep.append(d)
                        lat.append(la)
                        lon.append(lo)
        
    lat = np.array(lat)
    lon = np.array(lon)
    Value = np.array(Value)
    time = np.array(time)
    dep = np.array(dep)

    return Value, lat, lon, time 

def read_mat_sea(myfile, length, depth, sea, limits, variable):
    # Reads .mat files	
    # returns values and lon/lat as indices
    mat_aux = loadmat(myfile)

    Value=[]
    time=[]
    dep=[]
    lat=[]
    lon=[]

    if sea=='all':
        # North and Baltic Seas
        llon = 1
        ulon = 412
    elif sea=='North':
        # North Sea (west of 10deg East) only
        llon = 1
        ulon = 167
    elif sea=='Baltic':
        # Baltic Sea (east of 10deg East) only
        llon = 168
        ulon = 412
        

    mat = mat_aux['RefData']['data']['variable']
    if variable!='secchi_depth_of_sea_water':
        if variable!='mole_concentration_of_oxygen_in_seawater':
            for lo in range(llon,ulon): #=1:412:
                for la in range(347): #1:347:
                    for d in range(depth+1): #1:1%25:
                        for t in range(length):
                            if mat[t,d,la,lo] > -100 and mat[t,d,la,lo] < 1e6:
                                Value.append(mat[t,d,la,lo])
                                time.append(t)
                                dep.append(d)
                                lat.append(la)
                                lon.append(lo)
        else:
            #exclude outliers for oxygen
            for lo in range(llon,ulon): #=1:412:
                for la in range(347): #1:347:
                    for d in range(depth+1): #1:1%25:
                        for t in range(length):
                            if mat[t,d,la,lo] > limits[0] and mat[t,d,la,lo] < limits[1]:
#                            if mat[t,d,la,lo] > -100 and mat[t,d,la,lo] < 1e6:
#                            if mat[t,d,la,lo] > 240 and mat[t,d,la,lo] < 350:
                                Value.append(mat[t,d,la,lo])
                                time.append(t)
                                dep.append(d)
                                lat.append(la)
                                lon.append(lo)
    else:
        d=0
        for lo in range(llon,ulon): #=1:412:
            for la in range(347): #1:347:
                for t in range(length):
                    if mat[t,la,lo] > -100 and mat[t,la,lo] < 1e6:
                            Value.append(mat[t,la,lo])
                            time.append(t)
                            dep.append(d)
                            lat.append(la)
                            lon.append(lo)

    lat = np.array(lat)
    lon = np.array(lon)
    Value = np.array(Value)
    time = np.array(time)
    dep = np.array(dep)

    return Value, lat, lon, time 

def read_mat_coords(myfile, length, depth, varnum):
    # Reads .mat files
    # returns values and lon/lat as coordinates
    mat_aux = loadmat(myfile)

    Value=[]
    time=[]
    dep=[]
    lat=[]
    lon=[]
    mat = mat_aux['RefData']['data']['variable']
    mat_lat = mat_aux['RefData']['data']['lat']
    mat_lon = mat_aux['RefData']['data']['lon']

    if varnum!=12:
        if varnum!=9:
            for lo in range(412): #=1:412:
                for la in range(347): #1:347:
                    for d in range(depth+1): #1:1%25:
                        for t in range(length):
                            if mat[t,d,la,lo] > -100 and mat[t,d,la,lo] < 1e6:
                                Value.append(mat[t,d,la,lo])
                                time.append(t)
                                dep.append(d)
                                lat.append(mat_lat[la])
                                lon.append(mat_lon[lo])
        else:
            # exclude outliers for oxygen
            for lo in range(412): #=1:412:
                for la in range(347): #1:347:
                    for d in range(depth+1): #1:1%25:
                        for t in range(length):
#                            if mat[t,d,la,lo] > -100 and mat[t,d,la,lo] < 1e6:
                            if mat[t,d,la,lo] > 240 and mat[t,d,la,lo] < 350:
                                Value.append(mat[t,d,la,lo])
                                time.append(t)
                                dep.append(d)
                                lat.append(mat_lat[la])
                                lon.append(mat_lon[lo])
    else:
        d = 0
        for lo in range(412): #=1:412:
            for la in range(347): #1:347:
                for t in range(length):
                    if mat[t,la,lo] > -100 and mat[t,la,lo] < 1e6:
                            Value.append(mat[t,la,lo])
                            time.append(t)
                            dep.append(d)
                            lat.append(mat_lat[la])
                            lon.append(mat_lon[lo])

        
    lat = np.array(lat)
    lon = np.array(lon)
    Value = np.array(Value)
    time = np.array(time)
    dep = np.array(dep)

    return Value, lat, lon, time 


def get_sst(month, day, time):

    sst_fine, sst_coarse, _,_,_,_ = get_noaa_sst(month, day, time)

    return sst_fine, sst_coarse


def get_sat_sst(year, month, day, time):

#    if year<2017:
    data_f, data_c, lat_f, lat_c, lon_f, lon_c = get_noaa_sst(month, day, time)
#    else:
#        data_f, data_c, lat_f, lat_c, lon_f, lon_c = get_s3a_sst(month, day, time)

    return data_f, data_c, lat_f, lat_c, lon_f, lon_c


def get_noaa_sst(month, day, time):

    from netCDF4 import Dataset as NetCDFFile

    basedir = base_dir()

    if day==1:
        print 'Read NOAA satellite SST, day ', day

    if time == "am":
        time_aux = "00"
    elif time == "pm":
        time_aux = "12"

    if day < 10:
        day = '0'+str(day) 

    argv = '2012'+str(month)+str(day)+'_'+str(time_aux)
    sst_coarse   = NetCDFFile(basedir+'/inputs_2012/SST_Data/NO/NO_unzipped/'+str(argv)+'_no.nc')
    sst_coarse_keys   = sst_coarse.variables.keys()

    sst_fine   = NetCDFFile(basedir+'/inputs_2012/SST_Data/KU/KU_unzipped/'+str(argv)+'_ku.nc')
    sst_fine_keys   = sst_fine.variables.keys()

    coarse_keys = {}
    for i in range(len(sst_coarse_keys)):
        coarse_keys[str(sst_coarse_keys[i])+'_coarse'] = sst_coarse.variables[str(sst_coarse_keys[i])][:]
    sst_coarse = coarse_keys['mcsst_coarse'][:,:]
    lat_c = coarse_keys['lat_coarse'][:]
    lon_c = coarse_keys['lon_coarse'][:]

    fine_keys = {}
    for i in range(len(sst_fine_keys)):
        fine_keys[str(sst_fine_keys[i])+'_fine'] = sst_fine.variables[str(sst_fine_keys[i])][:]
    sst_fine = fine_keys['mcsst_fine'][:,:]
    lat_f = fine_keys['lat_fine'][:]
    lon_f = fine_keys['lon_fine'][:]

    return sst_fine, sst_coarse, lat_f, lat_c, lon_f, lon_c


def get_s3a_sst(month, day, time):
    from netCDF4 import Dataset as NetCDFFile

    basedir = base_dir()
    
    if day==1:
        print 'Read Sentinel-3 SST'
    
    if time == "am":
        time_aux = "00"
    elif time == "pm":
        time_aux = "12"

    if day < 10:
        day = '0'+str(day) 

    argv = '2017'+str(month)+str(day)+'_'+str(time_aux)
    sst_coarse   = NetCDFFile(basedir+'/inputs_2017/SST_Data/NO_unzipped/S3A_'+str(argv)+'_sst_no.nc')
    sst_coarse_keys   = sst_coarse.variables.keys()

    sst_fine   = NetCDFFile(basedir+'/inputs_2017/SST_Data/KU_unzipped/S3A_'+str(argv)+'_sst_ku.nc')
    sst_fine_keys   = sst_fine.variables.keys()

    coarse_keys = {}
    for i in range(len(sst_coarse_keys)):
        coarse_keys[str(sst_coarse_keys[i])] = sst_coarse.variables[str(sst_coarse_keys[i])][:]
    sst_c = coarse_keys['sea_surface_temperature'][:,:]
    lat_c = coarse_keys['lat'][:]
    lon_c = coarse_keys['lon'][:]

    fine_keys = {}
    for i in range(len(sst_fine_keys)):
        fine_keys[str(sst_fine_keys[i])] = sst_fine.variables[str(sst_fine_keys[i])][:]
    sst_f = fine_keys['sea_surface_temperature'][:,:]
    lat_f = fine_keys['lat'][:]
    lon_f = fine_keys['lon'][:]

    return sst_f, sst_c, lat_f, lat_c, lon_f, lon_c


def get_cmems_chl(year, month, day, time, datatype, log):

    from netCDF4 import Dataset as NetCDFFile

    basedir = base_dir()

    if day==1:
        print 'Read CMEMS satellite CHL, day ', day
        
    if day < 10:
        day = '0'+str(day)

    if datatype=='ba_MY':
        # Baltic Sea
        
        print '/scratch/usr/hzfblner/SEAMLESS/observations/CHL_BA_2015/chl_'+datatype+'_'+str(year)+str(month)+'.nc'
        ncid   = NetCDFFile('/scratch/usr/hzfblner/SEAMLESS/observations/CHL_BA_2015/chl_'+datatype+'_'+str(year)+str(month)+'.nc')

        chl = ncid.variables['CHL'][int(day)-1,:,:]
        lat = ncid.variables['lat'][:]
        lon = ncid.variables['lon'][:]

    else:

        # North Sea
        
        print '/scratch/usr/hzfblner/SEAMLESS/observations/CHL_NO_2015/chl_'+datatype+'_'+str(year)+str(month)+'.nc'
        ncid   = NetCDFFile('/scratch/usr/hzfblner/SEAMLESS/observations/CHL_NO_2015/chl_'+datatype+'_'+str(year)+str(month)+'.nc')

        chl = ncid.variables['CHL'][int(day)-1,:,:]
        lat = ncid.variables['latitude'][:]
        lon = ncid.variables['longitude'][:]

    if log==1:
        chl = np.log(chl)


    return chl, lat, lon


def get_mobs_cmems_chl(year, month, day, forana, path, datatype, chltype, log, ens):

    basedir = base_dir()

#    if day==1:
#        print 'Read observed model CHL, day ', day

    if day < 10:
        day = '0'+str(day)

    if ens==0:
       print path+'/chl_'+chltype+'_mobs_'+datatype+'_'+str(year)+str(month)+'.nc'
       ncid   = NetCDFFile(path+'/chl_'+chltype+'_mobs_'+datatype+'_'+str(year)+str(month)+'.nc')
    else:
        if ens < 10:
            ensstr = '00'+str(ens)
        else: 
            ensstr = '0'+str(ens)
        print path+'/chl_'+chltype+'_mobs_'+datatype+'_'+str(year)+str(month)+'.nc'
        ncid   = NetCDFFile(path+'/chl_'+chltype+'_mobs_'+datatype+'_'+str(year)+str(month)+'_'+ensstr+'.nc')

    chl = ncid.variables['CHL_'+forana][int(day)-1,:,:]
    if datatype=='no_multi':
        lat = ncid.variables['latitude'][:]
        lon = ncid.variables['longitude'][:]
    else:
        lat = ncid.variables['lat'][:]
        lon = ncid.variables['lon'][:]

    if log==1:
        chl = np.log(chl)

    return chl, lat, lon



def get_cmems_sst(year, month, day, time, datatype, ssttype):

    basedir = base_dir()

    if day==1:
        print 'Read CMEMS satellite SST, day ', day

    if time == "am":
        time_aux = "00"
    elif time == "pm":
        time_aux = "12"

    if day < 10:
        day = '0'+str(day)

    if ssttype=='L3S':
        print '/scratch/usr/hzfblner/SEAMLESS/observations/SST_L3S_2015/sst_'+datatype+'_'+str(year)+str(month)+'.nc'
        ncid   = NetCDFFile('/scratch/usr/hzfblner/SEAMLESS/observations/SST_L3S_2015/sst_'+datatype+'_'+str(year)+str(month)+'.nc')
        sst_coarse = ncid.variables['sea_surface_temperature'][int(day)-1,:,:]
    else:
        print '/scratch/usr/hzfblner/SEAMLESS/observations/SST_2015/sst_'+datatype+'_'+str(year)+str(month)+'.nc'
        ncid   = NetCDFFile('/scratch/usr/hzfblner/SEAMLESS/observations/SST_2015/sst_'+datatype+'_'+str(year)+str(month)+'.nc')
        sst_coarse = ncid.variables['analysed_sst'][int(day)-1,:,:]

    lat = ncid.variables['lat'][:]
    lon = ncid.variables['lon'][:]

    # Convert to degree Celsius
    sst_coarse = sst_coarse - 273.15

    return sst_coarse, lat, lon




def get_mobs_cmems_sst(year, month, day, forana, path, datatype, ens, ssttype):

    basedir = base_dir()

#    if day==1:
#        print 'Read CMEMS satellite SST, day ', day

#    if time == "am":
#        time_aux = "00"
#    elif time == "pm":
#        time_aux = "12"

    if day < 10:
        day = '0'+str(day)

    if ens==0:
        if ssttype=='L3S':
            print path+'/sst_L3S_mobs_'+datatype+'_'+str(year)+str(month)+'.nc'
            ncid   = NetCDFFile(path+'/sst_L3S_mobs_'+datatype+'_'+str(year)+str(month)+'.nc')
        else:
            print path+'/sst_mobs_'+datatype+'_'+str(year)+str(month)+'.nc'
            ncid   = NetCDFFile(path+'/sst_L4_mobs_'+datatype+'_'+str(year)+str(month)+'.nc')
    else:
      if ens < 10:
        ensstr = '00'+str(ens)
      else: 
        ensstr = '0'+str(ens)
      ncid   = NetCDFFile(path+'/sst_mobs_'+datatype+'_'+str(year)+str(month)+'_'+ensstr+'.nc')

#    if day < 10:
#        day = '0'+str(day) 
#    print '/scratch/usr/hzfblner/SEAMLESS/observations/SST_2015/sst_'+datatype+'_'+str(year)+str(month)+'.nc'
#    ncid   = NetCDFFile('/scratch/usr/hzfblner/SEAMLESS/observations/SST_2015/sst_'+datatype+'_'+str(year)+str(month)+'.nc')

    sst = ncid.variables['SST_'+forana][int(day)-1,:,:]
    lat = ncid.variables['lat'][:]
    lon = ncid.variables['lon'][:]

    # Convert to degree Celsius
    sst = sst

    return sst, lat, lon



def get_mean_sat_sst(year, month, day, time):

    # Function to compute the sst composite over one month

    _,_,_,ndays = month_names(int(month))
    print 'Month:', int(month), 'days', ndays

    ampm = 'am'
    iday = 1
    for i in range(ndays*2):

        print 'day', iday, ampm

        # Read satellite data
        if year<2017:
            sst_f, sst_c, lat_f, lat_c, lon_f, lon_c \
                    = get_noaa_sst(month, iday, ampm)
        else:
            sst_f, sst_c, lat_f, lat_c, lon_f, lon_c \
                    = get_s3a_sst(month, iday, ampm)

        if ampm == "am":
            time_aux = "00"
            ampm = "pm"
        elif ampm == "pm":
            time_aux = "12"
            ampm = "am"
            iday = iday+1

        if i==0:
            # Initialize arrays and counter at first step
            msst_f = np.ma.array(np.zeros(shape=(len(lat_f),len(lon_f))),mask=(np.zeros(shape=(len(lat_f),len(lon_f)))))
            msst_c = np.ma.array(np.zeros(shape=(len(lat_c),len(lon_c))),mask=(np.zeros(shape=(len(lat_c),len(lon_c)))))
            cntsst_f = np.zeros(shape=(len(lat_f),len(lon_f)))
            cntsst_c = np.zeros(shape=(len(lat_c),len(lon_c)))

        # Increment sums and counters
        for m in range(len(lat_f)):
            for n in range(len(lon_f)):
                if sst_f[m, n]>-100 and  sst_f[m, n] < 1e6:
                    msst_f[m,n] = msst_f[m,n] + sst_f[m,n]
                    cntsst_f[m,n] = cntsst_f[m,n]+1

        for m in range(len(lat_c)):
            for n in range(len(lon_c)):
                if sst_c[m, n]>-100 and  sst_c[m, n] < 1e6:
                    msst_c[m,n] = msst_c[m,n] + sst_c[m,n]
                    cntsst_c[m,n] = cntsst_c[m,n]+1

    for m in range(len(lat_f)):
        for n in range(len(lon_f)):
            if cntsst_f[m, n]>0:
                msst_f[m,n] = msst_f[m,n]/cntsst_f[m,n]
                msst_f.mask[m,n]=False
            else:
                msst_f.mask[m,n]=True

    for m in range(len(lat_c)):
        for n in range(len(lon_c)):
            if cntsst_c[m, n]>0:
                msst_c[m,n] = msst_c[m,n]/cntsst_c[m,n]
                msst_c.mask[m,n]=False
            else:
                msst_c.mask[m,n]=True

    return msst_f, msst_c, cntsst_f, cntsst_c, lat_f, lat_c, lon_f, lon_c


def load_data(arg1):
    ## Takes arg1 (variable) and finds the corresponding .nc file

    letter = 'state'
    # Temperature, Salt
    if arg1 == 'temp' or arg1 == 'salt':
        letter = 't'

    elif arg1 == 'elev' or arg1 == 'h' or arg1 == 'wvel':
        letter = 'z'
    # Eastward/Northward Wind Velocity
    elif arg1 == 'wuvel' or arg1 == 'wvvel':
        letter = 'm'
    # Eastward/Northward Sea Water Velocity
    elif arg1 == 'uvel' or arg1 == 'vvel':
        letter = 'c'
    # Casus, Sea Ice Thickness/Area fraction, Snow Thickness
    # Surface Temp Snow Ice, Eastward/Northward Ice Velocity
    elif arg1 == 'casus' or arg1 == 'ice_thk' or arg1 == 'ice_cov' \
      or arg1 == 'snow_thk' or arg1 == 'tsnei' or arg1 == 'ice_uvel' \
      or arg1 == 'ice_vvel':
        letter = 'i'
    # Secchi Depth of Sea Water
    elif arg1 == 'secchi_depth_of_sea_water':
        letter = 'secchidat'
    # Photosynthetically Active Radiation
    elif arg1 == 'photosynthetically_active_radiation':
        letter = 'light'
    # Mass Concentration of Chlorophyll in Sea Water
    elif arg1 == 'mass_concentration_of_chlorophyll_in_sea_water':
        letter = 'chldat'
    # Mole Concentration of ... in the Sea Water
    elif arg1 == 'mole_concentration_of_ammonium_in_seawater' or \
         arg1 == 'mole_concentration_of_nitrate_in_seawater' or \
         arg1 == 'mole_concentration_of_phosphate_in_seawater' or \
         arg1 == 'mole_concentration_of_silicate_in_seawater' or \
         arg1 == 'mole_concentration_of_diatoms_expressed_as_nitrogen_in_seawater' or \
         arg1 == 'mole_concentration_of_flagellates_expressed_as_nitrogen_in_seawater' or \
         arg1 == 'mole_concentration_of_cyanobacteria_expressed_as_nitrogen_in_seawater' or \
         arg1 == 'mole_concentration_of_zooplankton_expressed_as_nitrogen_in_seawater' or \
         arg1 == 'mole_concentration_of_protozooplankton_expressed_as_nitrogen_in_seawater' or \
         arg1 == 'mole_concentration_of_detritus_expressed_as_nitrogen_in_seawater' or \
         arg1 == 'mole_concentration_of_detritus_expressed_as_silicate_in_seawater' or \
         arg1 == 'mole_concentration_of_labile_dissolved_organic_nitrogen_expressed_as_nitrogen_in_seawater' or \
         arg1 == 'mole_concentration_of_dissolved_molecular_oxygen_in_sea_water':
        letter = 'biodat'
    # Mole Concentration of Nitrogen/Silicate in the Sediment
    elif arg1 == 'mole_concentration_of_nitrogen_in_the_sediment' or \
         arg1 == 'mole_concentration_of_silicate_in_the_sediment':
        letter = 'bendat'
    elif arg1 == 'kpar':
        letter = 'kpardat'
    elif arg1 == 'elev':
        letter = 'z' 

    else:
        print 'variable not linked'

    return letter

def param_and_vars(variable, coarse_keys, fine_keys, depth, DA_step):
    ## Takes the coarse/fine arg1 data from all the .nc file data

    numdim = coarse_keys[str(variable)+'_coarse'].ndim
    if numdim == 4:
        data_coarse = coarse_keys[str(variable)+'_coarse'][DA_step,depth,:,:]
        data_fine = fine_keys[str(variable)+'_fine'][DA_step,depth,:,:]
    elif numdim == 3:
        data_coarse = coarse_keys[str(variable)+'_coarse'][depth,:,:]
        data_fine = fine_keys[str(variable)+'_fine'][depth,:,:]

    return data_coarse, data_fine

def remove_errors(data_coarse, data_fine, arg1):
    ## This can be used to remove outliers of certain variables if needed. 
    ## So far, only one has been found in temp
    if arg1 == 'temp':
        import heapq
        A = np.where(data_fine < 0.0)[0][0]
        B = np.where(data_fine < 0.0)[1][0]
        data_fine[A,B] = heapq.nsmallest(2, data_fine.compressed())[-1]
        data_coarse[A,B] = heapq.nsmallest(2, data_fine.compressed())[-1]

    return data_coarse, data_fine



def read_station_idx(varnum, grid_area, year, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2, months, istation, station, dist):

    # Get indices in the model grid of a station according to its coordinates

    # Get variable names
    varstr, mat_var, variable, var_unit = var_names(varnum)
    print 'istation', istation
    # Station station name and indices
    if istation==0:
        lon_stat, lat_stat, _ = stations(station)
        station_coords = []
        station_coords.append(lat_stat)
        station_coords.append(lon_stat)
    else:
        station, station_coords, _ = station_info(istation, 'mod')

    if months[0]<10:
      month = '0'+str(months[0])
    else:
      month = months[0]
   
    day = 1
    _, mod_keys = read_model(varstr, grid_area, year, month, day, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2)

    lat = mod_keys['lat'][:]
    lon = mod_keys['lon'][:]

    lat1 = lat[:,0]
    lon1 = lon[0,:]
    for i in range(len(lat[:,0])):
       lat1[i] = np.max(lat[i,:])
    for j in range(len(lat[0,:])):
       lon1[j] = np.max(lon[:,j])
    lon1[-1] = lon1[-2] + lon1[1]-lon1[0]

    latdiff = lat1 - station_coords[0]
    londiff = lon1 - station_coords[1]
    station_idx = []
    station_idx.append(np.argmin(abs(latdiff)))
    station_idx.append(np.argmin(abs(londiff)))

    print 'Station:', station
    print '   model grid: coords lat/lon', station_coords[0], station_coords[1]
    print '   model grid: index lat/lon ', station_idx[0], station_idx[1]
    print '   model grid: lat/lon', lat1[station_idx], lon1[station_idx]

    return station, station_idx, station_coords



def read_station_idx_obs(varnum, year, months, istation, station, dist):

    # Get indices in the observation grid of a station according to its coordinates

    # Station station name and indices
    if istation==0:
        lon_stat, lat_stat, noba = stations(station)
        station_coords = []
        station_coords.append(lat_stat)
        station_coords.append(lon_stat)
    else:
        station, station_coords, noba = station_info(istation, 'mod')

    if months[0]<10:
        month = '0'+str(months[0])
    else:
        month = str(months[0])

    day = 1
    if varnum==24:
        # Chlorophyll
        log = 0
        datatype_obs= noba+'_MY'
        _, lat_o, lon_o = get_cmems_chl(year, month, day, 'am', datatype_obs, log)
    else:
        # SST
        datatype_obs= 'REP_L4'
	datatype_obs = 'L3S'
        _, lat_o, lon_o = get_cmems_sst(year, month, day, 'am', datatype_obs, noba)


    latdiff = lat_o - station_coords[0]
    londiff = lon_o - station_coords[1]
    station_idx = []
    station_idx.append(np.argmin(abs(latdiff)))
    station_idx.append(np.argmin(abs(londiff)))

    print 'Station:', station
    print '   obs grid: coords lat/lon', station_coords[0], station_coords[1]
    print '   obs grid: index lat/lon ', station_idx[0], station_idx[1]

    return station, station_idx, station_coords, noba



def read_station_series(varnum, grid_area, year, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2, months, istation, strstation, dist):

    # Get day of year
    DoY = day_of_year(months[0], 1)
    print 'DoY ', DoY

    # Get variable names
    varstr, mat_var, variable, var_unit = var_names(varnum)
    print 'istation1', istation
    # Station station name and indices
    station, station_idx, station_coords =  read_station_idx(varnum, grid_area, year, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2, months, istation, strstation, dist)

    lon_min = station_idx[1]-dist
    lon_max = station_idx[1]+dist+1
    lat_min = station_idx[0]-dist
    lat_max = station_idx[0]+dist+1
    range_lon = lon_max-lon_min
    range_lat = lat_max-lat_min
    range_lon = lon_max-lon_min
    range_lat = lat_max-lat_min

    modmean = []
    modday = []

    alldays = DoY

    for imonth in range(len(months)):

      _, _, _, ndays = month_names(months[imonth])

      if months[imonth]<10:
        month = '0'+str(months[imonth])
      else:
        month = str(months[imonth])

      for iday in range(ndays):

        day = iday+1

        data_mod, mod_keys = read_model(varstr, grid_area, year, month, day, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2)

        # set coordinates
        lat = mod_keys['lat'][:]
        lon = mod_keys['lon'][:]
	dep_temp = mod_keys['lev'][:]
	if iday == 0:
	  dep = dep_temp 
	else: 
	  dep = np.vstack((dep,dep_temp))

        if day==1:
            print 'check observations within ', dist, ' grid points around location ', station
            print 'longitude range: ', lon[lon_min], lon[lon_max]
            print 'longitude range: ', lat[lat_min], lat[lat_max]

	if depth == ':':
          modval= np.zeros((len(dep_temp), range_lon, range_lat))
          modval = data_mod[:,lat_min:lat_max, lon_min:lon_max]	
          mean_mobs= np.zeros(len(dep_temp))
          cnt_mean = np.zeros(len(dep_temp))
          for d in range(len(modval)):
            for i in range(len(modval[0])):
              for j in range(len(modval[0][0])):
                if modval[d,i,j]>-999.0:
                    mean_mobs[d] = mean_mobs[d] + modval[d,i,j]
                    cnt_mean[d] = cnt_mean[d] + 1
            if cnt_mean[d]>0:
              mean_mobs[d] = mean_mobs[d] / cnt_mean[d]
	else:
          modval= np.zeros((range_lon, range_lat))
          modval = data_mod[lat_min:lat_max, lon_min:lon_max]
          mean_mobs= 0.0
          cnt_mean = 0
          for i in range(len(modval[:,0])):
            for j in range(len(modval[0,:])):
                if modval[i,j]>-999.0:
                    mean_mobs = mean_mobs + modval[i,j]
                    cnt_mean = cnt_mean + 1
          if cnt_mean>0:
            mean_mobs = mean_mobs / cnt_mean

        if varnum==10:
            if mean_mobs<0.0:
                mean_mobs = 1.0e-05
        modmean.append(mean_mobs)
        modday.append(alldays)

        alldays = alldays+1
    if depth == ':':
      return modmean, modday, dep 
    else: 
      return modmean, modday

def read_station_series_path(varnum, grid_area, year, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2, months, istation, strstation, dist, path):

    # Get day of year
    DoY = day_of_year(months[0], 1)
    print 'DoY ', DoY

    # Get variable names
    varstr, mat_var, variable, var_unit = var_names(varnum)
    print 'istation1', istation
    # Station station name and indices
    station, station_idx, station_coords =  read_station_idx(varnum, grid_area, year, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2, months, istation, strstation, dist)

    lon_min = station_idx[1]-dist
    lon_max = station_idx[1]+dist+1
    lat_min = station_idx[0]-dist
    lat_max = station_idx[0]+dist+1
    range_lon = lon_max-lon_min
    range_lat = lat_max-lat_min
    range_lon = lon_max-lon_min
    range_lat = lat_max-lat_min

    modmean = []
    modday = []

    alldays = DoY

    for imonth in range(len(months)):

      _, _, _, ndays = month_names(months[imonth])
      if months[imonth]<10:
        month = '0'+str(months[imonth])
      else:
        month = str(months[imonth])

      for iday in range(ndays):

        day = iday+1
        data_mod, mod_keys = read_model_path(varstr, grid_area, year, month, day, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2, path)

        # set coordinates
        lat = mod_keys['lat'][:]
        lon = mod_keys['lon'][:]

        dep_temp = mod_keys['lev'][:]
        if iday == 0:
          dep = dep_temp
        else:
          dep = np.vstack((dep,dep_temp))

        if day==1:
            print 'check observations within ', dist, ' grid points around location ', station
            print 'longitude range: ', lon[lon_min], lon[lon_max]
            print 'latitude range: ', lat[lat_min], lat[lat_max]

        if depth==':':
          modval= np.zeros((len(dep_temp), range_lon, range_lat))
          modval = data_mod[:,lat_min:lat_max, lon_min:lon_max]
        else:
          modval= np.zeros((range_lon, range_lat))
          modval = data_mod[lat_min:lat_max, lon_min:lon_max]

        if depth==':':
          mean_mobs= np.zeros(len(dep_temp))
          cnt_mean = np.zeros(len(dep_temp))
        else:
          mean_mobs = 0.0
          cnt_mean = 0

        if depth==':':
          for d in range(len(modval)):
            for i in range(len(modval[0])):
              for j in range(len(modval[0][0])):
                if modval[d,i,j]>-999.0:
                    mean_mobs[d] = mean_mobs[d] + modval[d,i,j]
                    cnt_mean[d] = cnt_mean[d] + 1
            if cnt_mean[d]>0:
              mean_mobs[d] = mean_mobs[d] / cnt_mean[d]
        else:
          for i in range(len(modval[:,0])):
            for j in range(len(modval[0,:])):
                if modval[i,j]>-999.0:
                    mean_mobs = mean_mobs + modval[i,j]
                    cnt_mean = cnt_mean + 1
          if cnt_mean>0:
            mean_mobs = mean_mobs / cnt_mean

        if varnum==10:
            if mean_mobs<0.0:
                mean_mobs = 1.0e-05

        modmean.append(mean_mobs)
        modday.append(alldays)

        alldays = alldays+1
    if depth==':':
      return modmean, modday, dep
    else:
      return modmean, modday

def read_latlon_series(varnum, grid_area, year, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2, months, lat_stat, lon_stat, dist):

    # Get day of year
    DoY = day_of_year(months[0], 1)
    print 'DoY ', DoY

    # Get variable names
    varstr, mat_var, variable, var_unit = var_names(varnum)
    model_idx = get_model_idx(lat_stat, lon_stat, dist)
    lon_min = model_idx[1]-dist
    lon_max = model_idx[1]+dist+1
    lat_min = model_idx[0]-dist
    lat_max = model_idx[0]+dist+1
    range_lon = lon_max-lon_min
    range_lat = lat_max-lat_min

    modmean = []
    modday = []

    alldays = DoY

    for imonth in range(len(months)):
      _, _, _, ndays = month_names(months[imonth])
      if months[imonth]<10:
        month = '0'+str(months[imonth])
      else:
        month = str(months[imonth])

      for iday in range(ndays):

        day = iday+1

        data_mod, mod_keys = read_model(varstr, grid_area, year, month, day, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2)

        # set coordinates
        lat = mod_keys['lat'][:]
        lon = mod_keys['lon'][:]

        dep_temp = mod_keys['lev'][:]
        if iday == 0:
          dep = dep_temp
        else:
          dep = np.vstack((dep,dep_temp))

        if day==1:
            print 'check observations within ', dist, ' grid points around location lat/lon = ', lat_stat, lon_stat
            print 'longitude range: ', lon[lat_min, lon_min], lon[lat_max, lon_max]
            print 'latitude range: ', lat[lat_min, lon_min], lat[lat_max, lon_max]
            if isinstance(lon[lat_min, lon_min], np.ma.core.MaskedConstant):
              print lon[lat_min:lat_max, lon_min:lon_max]
              print 'longitude range first real value: ', np.ma.minimum.reduce(lon[lat_min:lat_max, lon_min:lon_max]), np.max(lon[lat_min:lat_max, lon_min:lon_max])
            if isinstance(lat[lat_min, lon_min], np.ma.core.MaskedConstant):
              print lat[lat_min:lat_max, lon_min:lon_max]
              print 'latitude range first real value: ', np.ma.min(lat[lat_min:lat_max, lon_min:lon_max]), np.max(lat[lat_min:lat_max, lon_min:lon_max])
        if depth==':':
          modval= np.zeros((len(dep_temp), range_lat, range_lon))
          modval = data_mod[:,lat_min:lat_max, lon_min:lon_max]
        else:
          modval= np.zeros((range_lat, range_lon))
          modval = data_mod[lat_min:lat_max, lon_min:lon_max]

        if depth==':':
          mean_mobs= np.zeros(len(dep_temp))
          cnt_mean = np.zeros(len(dep_temp))
        else:
          mean_mobs = 0.0
          cnt_mean = 0

        if depth==':':
          for d in range(len(modval)):
            for i in range(len(modval[0])):
              for j in range(len(modval[0][0])):
                if modval[d,i,j]>-999.0:
                    mean_mobs[d] = mean_mobs[d] + modval[d,i,j]
                    cnt_mean[d] = cnt_mean[d] + 1
            if cnt_mean[d]>0:
              mean_mobs[d] = mean_mobs[d] / cnt_mean[d]
        else:
          for i in range(len(modval[:,0])):
            for j in range(len(modval[0,:])):
                if modval[i,j]>-999.0:
                    mean_mobs = mean_mobs + modval[i,j]
                    cnt_mean = cnt_mean + 1
          if cnt_mean>0:
            mean_mobs = mean_mobs / cnt_mean

        if varnum==10:
            if mean_mobs<0.0:
                mean_mobs = 1.0e-05

        modmean.append(mean_mobs)
        modday.append(alldays)

        alldays = alldays+1
    if depth==':':
      return modmean, modday, dep
    else:
      return modmean, modday

def read_latlon_day(varnum, grid_area, year, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2, month, day, model_idx, dist):

    # if month and day input is number convert to string
    if not isinstance(month, str):
      if month<10:
        month = '0'+str(month)
      else:
        month = str(month)
    if not isinstance(day, str):
      if day<10:
        day = '0'+str(day)
    else:
        day = str(day)

    # Get variable names
    varstr, mat_var, variable, var_unit = var_names(varnum)
    #model_idx = get_model_idx(lat_stat, lon_stat, dist)
    lon_min = model_idx[1]-dist
    lon_max = model_idx[1]+dist+1
    lat_min = model_idx[0]-dist
    lat_max = model_idx[0]+dist+1
    range_lon = lon_max-lon_min
    range_lat = lat_max-lat_min

    modmean = []
    modday = []

    data_mod, mod_keys = read_model(varstr, grid_area, year, month, day, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2)

    # set coordinates
    lat = mod_keys['lat'][:]
    lon = mod_keys['lon'][:]
    dep = mod_keys['lev'][:]

    if day==1:
        print 'check observations within ', dist, ' grid points around location lat/lon = ', lat_stat, lon_stat
        print 'longitude range: ', lon[lat_min, lon_min], lon[lat_max, lon_max]
        print 'latitude range: ', lat[lat_min, lon_min], lat[lat_max, lon_max]
        if isinstance(lon[lat_min, lon_min], np.ma.core.MaskedConstant):
          print lon[lat_min:lat_max, lon_min:lon_max]
          print 'longitude range first real value: ', np.ma.minimum.reduce(lon[lat_min:lat_max, lon_min:lon_max]), np.max(lon[lat_min:lat_max, lon_min:lon_max])
        if isinstance(lat[lat_min, lon_min], np.ma.core.MaskedConstant):
          print lat[lat_min:lat_max, lon_min:lon_max]
          print 'latitude range first real value: ', np.ma.min(lat[lat_min:lat_max, lon_min:lon_max]), np.max(lat[lat_min:lat_max, lon_min:lon_max])
    if depth==':':
      modval= np.zeros((len(dep), range_lat, range_lon))
      modval = data_mod[:,lat_min:lat_max, lon_min:lon_max]
    else:
      modval= np.zeros((range_lat, range_lon))
      modval = data_mod[lat_min:lat_max, lon_min:lon_max]

    if depth==':':
      mean_mobs= np.zeros(len(dep))
      cnt_mean = np.zeros(len(dep))
    else:
      mean_mobs = 0.0
      cnt_mean = 0

    if depth==':':
      for d in range(len(modval)):
        for i in range(len(modval[0])):
          for j in range(len(modval[0][0])):
            if modval[d,i,j]>-999.0:
                mean_mobs[d] = mean_mobs[d] + modval[d,i,j]
                cnt_mean[d] = cnt_mean[d] + 1
        if cnt_mean[d]>0:
          mean_mobs[d] = mean_mobs[d] / cnt_mean[d]
    else:
      for i in range(len(modval[:,0])):
        for j in range(len(modval[0,:])):
           if modval[i,j]>-999.0:
                mean_mobs = mean_mobs + modval[i,j]
                cnt_mean = cnt_mean + 1
      if cnt_mean>0:
        mean_mobs = mean_mobs / cnt_mean

    if varnum==10:
      if mean_mobs<0.0:
        mean_mobs = 1.0e-05

    modmean.append(mean_mobs)

    if depth==':':
      return modmean[0], dep
    else:
      return modmean

def read_latlon_day_path(varnum, grid_area, year, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2, month, day, model_idx, dist, path):

    # if month and day input is number convert to string
    if not isinstance(month, str):
      if month<10:
        month = '0'+str(month)
      else:
        month = str(month)
    if not isinstance(day, str):
      if day<10:
        day = '0'+str(day)
    else:
        day = str(day)

    # Get variable names
    varstr, mat_var, variable, var_unit = var_names(varnum)
    #model_idx = get_model_idx(lat_stat, lon_stat, dist)
    lon_min = model_idx[1]-dist
    lon_max = model_idx[1]+dist+1
    lat_min = model_idx[0]-dist
    lat_max = model_idx[0]+dist+1
    range_lon = lon_max-lon_min
    range_lat = lat_max-lat_min

    modmean = []
    modday = []

    data_mod, mod_keys = read_model_path(varstr, grid_area, year, month, day, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2, path)

    # set coordinates
    lat = mod_keys['lat'][:]
    lon = mod_keys['lon'][:]
    dep = mod_keys['lev'][:]

    if day==1:
        print 'check observations within ', dist, ' grid points around location lat/lon = ', lat_stat, lon_stat
        print 'longitude range: ', lon[lat_min, lon_min], lon[lat_max, lon_max]
        print 'latitude range: ', lat[lat_min, lon_min], lat[lat_max, lon_max]
        if isinstance(lon[lat_min, lon_min], np.ma.core.MaskedConstant):
          print lon[lat_min:lat_max, lon_min:lon_max]
          print 'longitude range first real value: ', np.ma.minimum.reduce(lon[lat_min:lat_max, lon_min:lon_max]), np.max(lon[lat_min:lat_max, lon_min:lon_max])
        if isinstance(lat[lat_min, lon_min], np.ma.core.MaskedConstant):
          print lat[lat_min:lat_max, lon_min:lon_max]
          print 'latitude range first real value: ', np.ma.min(lat[lat_min:lat_max, lon_min:lon_max]), np.max(lat[lat_min:lat_max, lon_min:lon_max])
    if depth==':':
      modval= np.zeros((len(dep), range_lat, range_lon))
      modval = data_mod[:,lat_min:lat_max, lon_min:lon_max]
    else:
      modval= np.zeros((range_lat, range_lon))
      modval = data_mod[lat_min:lat_max, lon_min:lon_max]

    if depth==':':
      mean_mobs= np.zeros(len(dep))
      cnt_mean = np.zeros(len(dep))
    else:
      mean_mobs = 0.0
      cnt_mean = 0

    if depth==':':
      for d in range(len(modval)):
        for i in range(len(modval[0])):
          for j in range(len(modval[0][0])):
            if modval[d,i,j]>-999.0:
                mean_mobs[d] = mean_mobs[d] + modval[d,i,j]
                cnt_mean[d] = cnt_mean[d] + 1
        if cnt_mean[d]>0:
          mean_mobs[d] = mean_mobs[d] / cnt_mean[d]
    else:
      for i in range(len(modval[:,0])):
        for j in range(len(modval[0,:])):
           if modval[i,j]>-999.0:
                mean_mobs = mean_mobs + modval[i,j]
                cnt_mean = cnt_mean + 1
      if cnt_mean>0:
        mean_mobs = mean_mobs / cnt_mean

    if varnum==10:
      if mean_mobs<0.0:
        mean_mobs = 1.0e-05

    modmean.append(mean_mobs)

    if depth==':':
      return modmean[0], dep
    else:
      return modmean


def read_station_series_obs(varnum, year, months, istation, strstation, dist):

    DoY = day_of_year(months[0], 1)
    print 'DoY ', DoY
 
    # Station station name and indices
    station, station_idx, station_coords, noba =  read_station_idx_obs(varnum, year, months, istation, strstation, dist)

    lon_min = station_idx[1]-dist
    lon_max = station_idx[1]+dist+1
    lat_min = station_idx[0]-dist
    lat_max = station_idx[0]+dist+1
    range_lon = lon_max-lon_min
    range_lat = lat_max-lat_min
    range_lon = lon_max-lon_min
    range_lat = lat_max-lat_min

    obsmean = []
    obsday = []

    alldays = DoY
    
    for imonth in range(len(months)):
        
      _, _, _, ndays = month_names(months[imonth])

      if months[imonth]<10:
        month = '0'+str(months[imonth])
      else:
        month = str(months[imonth])

      for iday in range(ndays):

        day = iday+1

        if varnum==24:
            # Chlorophyll
            log = 0
            datatype_obs= noba+'_MY'
            data_o, lat_o, lon_o = get_cmems_chl(year, month, day, 'am', datatype_obs, log)
        else:
            # SST
            datatype_obs= 'REP_L4'
	    datatype_obs = 'L3S'
            data_o, lat_o, lon_o = get_cmems_sst(year, month, day, 'am', datatype_obs, noba)

        if day==1:
            print 'check observations within ', dist, ' grid points around location ', station
            print 'longitude range: ', lon_o[lon_min], lon_o[lon_max]
            print 'longitude range: ', lat_o[lat_min], lat_o[lat_max]

        obsval= np.zeros((range_lon, range_lat))
        obsval = data_o[lat_min:lat_max, lon_min:lon_max]

        #print 'modval: ', modval


        mean_obs= 0.0
        cnt_mean = 0
        for i in range(len(obsval[:,0])):
            for j in range(len(obsval[0,:])):
                if obsval[i,j]>-999.0:
                    mean_obs = mean_obs + obsval[i,j]
                    cnt_mean = cnt_mean + 1
        if cnt_mean>0:
            mean_obs = mean_obs / cnt_mean

        if cnt_mean>0:
            obsmean.append(mean_obs)
            obsday.append(alldays)

        alldays = alldays+1

    return obsmean, obsday



def read_model(variable, grid_area, year, month, day, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2):

    if day < 10:
        day = '0'+str(day)
    else:
        day = str(day)

    letter = 'state' #load_data(variable)
    print 'Read variable:', variable

    if depth == ':':
      print 'depth: entire column'
    elif depth != '0':
      print 'Depth level '+depth

    # Set Scale factor
    scale = 1.0

    # Get paths for assimilation experiment and free run
    #if day=='01':
    if DA_switch==1 or DA_switch==2 or DA_switch==11 or DA_switch==12:
       path = get_exp_path(coupled, year, month, day, time_stamp, 1)
    elif DA_switch==0:
       path = get_exp_path(coupled, year, month, day, time_stamp, 0)

    if variable == 'PFT' or variable == 'TE':
      path = path.replace('DA','Post_DA')
    # Set file name including path
    timestr = str(year)+str(month)+str(day)

    print path+str(letter)+'_'+timestr+'.nc'
    if DA_switch == 1 or DA_switch==2 or DA_switch==11 or DA_switch==12 : # Background or analysis
        nc_coarse = NetCDFFile(path+str(letter)+'_'+timestr+'.nc')
    elif DA_switch == 0: # Free Ens
        nc_coarse = NetCDFFile(path+str(letter)+'_'+timestr+'.nc')

    # Set file index to read
    #    0: forecast, 1: analysis
    if DA_switch == 1: #Background
        DA_step =0
    elif DA_switch == 2: #Analysis
        DA_step =1
    elif DA_switch == 0: #Free Ens
        DA_step = 0
    print 'DA_step', DA_step
    nc_coarse_keys = ['nav_lon', 'nav_lat', 'nav_lev', str(variable)]

    # Coarse Data
    coarse_keys = {}
    data_coarse = {}
    if grid_area == 'all' or grid_area == 'coarse':
        coarse_keys['lon'] = nc_coarse.variables['nav_lon'][:]
        coarse_keys['lat'] = nc_coarse.variables['nav_lat'][:]
	coarse_keys['lev'] = nc_coarse.variables['nav_lev'][:]
        numdim =  nc_coarse.variables[str(variable)].ndim
        if numdim == 4:
            if z_mean == 0 and z_integral == 0:
	        if depth == ':':
                  data_coarse = nc_coarse.variables[str(variable)][DA_step,:,:,:]
		else:
                  data_coarse = nc_coarse.variables[str(variable)][DA_step,depth,:,:]
            elif z_mean == 1 and z_integral ==0:
                data_coarse_pre = nc_coarse.variables[str(variable)][DA_step,:,:,:]
                data_coarse = vertical_mean_4D(data_coarse_pre, z1, z2)
            elif z_mean == 0 and z_integral ==1:
                data_coarse_pre = nc_coarse.variables[str(variable)][DA_step,:,:,:]
                data_coarse = vertical_integral_4D(data_coarse_pre, z1, z2)
            else:
                print 'Wrong z_mean or z_integral'
                exit()
        elif numdim == 3:
            if str(variable)=='SST_inst':
                data_coarse = nc_coarse.variables[str(variable)][DA_step,:,:]
            else:
		if depth == ':':
                  data_coarse = nc_coarse.variables[str(variable)][:,:,:]
		else:
                  data_coarse = nc_coarse.variables[str(variable)][depth,:,:]

        # Rescale data
        data_coarse = data_coarse * scale

    return data_coarse, coarse_keys

def read_model_path(variable, grid_area, year, month, day, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2, path):

    if day < 10:
        day = '0'+str(day)
    else:
        day = str(day)

    letter = 'state' #load_data(variable)
    print 'Read variable:', variable

    if depth != '0':
      if depth == ':':
        print 'depth: entire column'
      else:
        print 'Depth level '+depth

    # Set Scale factor
    scale = 1.0

    # Set file name including path
    timestr = str(year)+str(month)+str(day)

    print path+str(letter)+'_'+timestr+'.nc'
    if DA_switch == 1 or DA_switch==2 or DA_switch==11 or DA_switch==12 : # Background or analysis
        nc_coarse = NetCDFFile(path+str(letter)+'_'+timestr+'.nc')
    elif DA_switch == 0: # Free Ens
        nc_coarse = NetCDFFile(path+str(letter)+'_'+timestr+'.nc')

    # Set file index to read
    #    0: forecast, 1: analysis
    if DA_switch == 1: #Background
        DA_step =0
    elif DA_switch == 2: #Analysis
        DA_step =1
    elif DA_switch == 0: #Free Ens
        DA_step = 0
    print 'DA_step', DA_step
    nc_coarse_keys = ['nav_lon', 'nav_lat', 'nav_lev', str(variable)]

    # Coarse Data
    coarse_keys = {}
    data_coarse = {}
    if grid_area == 'all' or grid_area == 'coarse':
        coarse_keys['lon'] = nc_coarse.variables['nav_lon'][:]
        coarse_keys['lat'] = nc_coarse.variables['nav_lat'][:]
        coarse_keys['lev'] = nc_coarse.variables['nav_lev'][:]
        numdim =  nc_coarse.variables[str(variable)].ndim
        if numdim == 4:
            if z_mean == 0 and z_integral == 0:
                if depth==':':
                  data_coarse = nc_coarse.variables[str(variable)][DA_step,:,:,:]
                else:
                  data_coarse = nc_coarse.variables[str(variable)][DA_step,depth,:,:]
            elif z_mean == 1 and z_integral ==0:
                data_coarse_pre = nc_coarse.variables[str(variable)][DA_step,:,:,:]
                data_coarse = vertical_mean_4D(data_coarse_pre, z1, z2)
            elif z_mean == 0 and z_integral ==1:
                data_coarse_pre = nc_coarse.variables[str(variable)][DA_step,:,:,:]
                data_coarse = vertical_integral_4D(data_coarse_pre, z1, z2)
            else:
                print 'Wrong z_mean or z_integral'
                exit()
        elif numdim == 3:
            if str(variable)=='SST_inst':
                data_coarse = nc_coarse.variables[str(variable)][DA_step,:,:]
            else:
                if depth==':':
                  data_coarse = nc_coarse.variables[str(variable)][:,:,:]
                else:
                  data_coarse = nc_coarse.variables[str(variable)][depth,:,:]


        # Rescale data
        data_coarse = data_coarse * scale

    return data_coarse, coarse_keys

def read_model_var(variable, grid_area, year, month, day, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2):

    if day < 10:
        day = '0'+str(day)
    else:
        day = str(day)

    letter = 'variance' #load_data(variable)
    print 'Read variable:', variable

    if depth != '0':
      print 'Depth level '+depth

    # Set Scale factor
    scale = 1.0

    # Get paths for assimilation experiment and free run
    #if day=='01':
    if DA_switch==1 or DA_switch==2 or DA_switch==11 or DA_switch==12:
       path = get_exp_path(coupled, year, month, day, time_stamp, 1)
    elif DA_switch==0:
       path = get_exp_path(coupled, year, month, day, time_stamp, 0)

    # Set file name including path
    timestr = str(year)+str(month)+str(day)

    print path+str(letter)+'_'+timestr+'.nc'
    if DA_switch == 1 or DA_switch==2 or DA_switch==11 or DA_switch==12 : # Background or analysis
        nc_coarse = NetCDFFile(path+str(letter)+'_'+timestr+'.nc')
    elif DA_switch == 0: # Free Ens
        nc_coarse = NetCDFFile(path+str(letter)+'_'+timestr+'.nc')

    # Set file index to read
    #    0: forecast, 1: analysis
    if DA_switch == 1: #Background
        DA_step =0
    elif DA_switch == 2: #Analysis
        DA_step =1
    elif DA_switch == 0: #Free Ens
        DA_step = 0
    print 'DA_step', DA_step
    nc_coarse_keys = ['nav_lon', 'nav_lat', str(variable)]

    # Coarse Data
    coarse_keys = {}
    data_coarse = {}
    if grid_area == 'all' or grid_area == 'coarse':
        coarse_keys['lon'] = nc_coarse.variables['nav_lon'][:]
        coarse_keys['lat'] = nc_coarse.variables['nav_lat'][:]

        numdim =  nc_coarse.variables[str(variable)].ndim
        if numdim == 4:
            if z_mean == 0 and z_integral == 0:
                data_coarse = nc_coarse.variables[str(variable)][DA_step,depth,:,:]
            elif z_mean == 1 and z_integral ==0:
                data_coarse_pre = nc_coarse.variables[str(variable)][DA_step,:,:,:]
                data_coarse = vertical_mean_4D(data_coarse_pre, z1, z2)
            elif z_mean == 0 and z_integral ==1:
                data_coarse_pre = nc_coarse.variables[str(variable)][DA_step,:,:,:]
                data_coarse = vertical_integral_4D(data_coarse_pre, z1, z2)
            else:
                print 'Wrong z_mean or z_integral'
                exit()
        elif numdim == 3:
            if str(variable)=='SST_inst':
                data_coarse = nc_coarse.variables[str(variable)][DA_step,:,:]
            else:
                data_coarse = nc_coarse.variables[str(variable)][depth,:,:]

        # Rescale data
        data_coarse = data_coarse * scale

    return data_coarse, coarse_keys


# Read from file with hourly output
def read_model_hr(variable, grid_area, year, month, day, time_stamp, hour, depth, DA_switch, coupled):

    ## Main code
    if day < 10:
        day = '0'+str(day)	
    else:
        day = str(day)

    letter = load_data(variable)	
#    print 'var:', variable
    # Set Scale factor
#    if letter=='secchidat' or letter=='light' or letter=='biodat' or letter=='bendat' \
#           or letter=='kpardat' or letter=='chldat':
#        if variable=='mole_concentration_of_dissolved_molecular_oxygen_in_sea_water':
#        scale=1.0
#        else:
#           scale = 1.0/1000.0
#    else:
    scale = 1.0

    # Get paths for assimilation experiment and free run
    path = get_exp_path(coupled, year, month, day, time_stamp, 1)
    if day=='01':
        if DA_switch==1 or DA_switch==2:
            print path
    path_free = get_exp_path(coupled, year, month, day, time_stamp, 0)
    if day=='01' and DA_switch==0:
        print path_free

    # Set file name including path
    timestr = str(year)+str(month)+str(day)+str(time_stamp)

    if DA_switch == 1 or DA_switch==2 : # Background or analysis
        if grid_area == 'all':
            nc_coarse = NetCDFFile(path+str(letter)+'_pdaf_coarse_'+timestr+'.nc')
            nc_fine   = NetCDFFile(path+str(letter)+'_pdaf_fine_'+timestr+'.nc')
        elif grid_area == 'fine':
            nc_fine   = NetCDFFile(path+str(letter)+'_pdaf_fine_'+timestr+'.nc')
            nc_coarse = nc_fine
        elif grid_area == 'coarse':
            nc_coarse = NetCDFFile(path+str(letter)+'_pdaf_coarse_'+timestr+'.nc')
            nc_fine = nc_coarse
    elif DA_switch == 0: # Free Ens
        if grid_area == 'all':
            nc_coarse = NetCDFFile(path_free+str(letter)+'_file_coarse_'+timestr+'.nc')
            nc_fine   = NetCDFFile(path_free+str(letter)+'_file_fine_'+timestr+'.nc')
        elif grid_area == 'fine':
            nc_fine   = NetCDFFile(path_free+str(letter)+'_file_fine_'+timestr+'.nc')
            nc_coarse = nc_fine
        elif grid_area == 'coarse':
            nc_coarse = NetCDFFile(path_free+str(letter)+'_file_coarse_'+timestr+'.nc')
            nc_fine = nc_coarse

    # Set file index to read
    if DA_switch == 1: #Background
        # If the output only contains fields (write_variances=.false.):
        #    0: current forecast, 1: analysis, 2: next 12h forecast
        # the 12h forecast is identical to the current forecast in the next output file
        # If write_variances=.true. is set, the output contains variances and fields. Then
        #    0,2,4: variance forcast, analysis, next 12h forecast;
        #    1,3,5: fields: forecast, analysis, next 12h forecast
        DA_step = 0
    elif DA_switch == 2: #Analysis
        DA_step = 1
    elif DA_switch == 0: #Free Ens
        DA_step = hour-1

    nc_coarse_keys = ['lon', 'lat', str(variable)]
    nc_fine_keys = ['lon', 'lat', str(variable)]

    # Coarse Data
    coarse_keys = {}
    data_coarse = {}
    if grid_area == 'all' or grid_area == 'coarse':
        coarse_keys['lon_coarse'] = nc_coarse.variables['lon'][:]
        coarse_keys['lat_coarse'] = nc_coarse.variables['lat'][:]
         
        numdim =  nc_coarse.variables[str(variable)].ndim
        if numdim == 4:
            data_coarse = nc_coarse.variables[str(variable)][DA_step,depth,:,:]
        elif numdim == 3:
            data_coarse = nc_coarse.variables[str(variable)][DA_step,:,:]

        # Rescale data
        data_coarse = data_coarse * scale

    # Fine Data
    fine_keys = {}
    data_fine = {}

    if grid_area == 'all' or grid_area == 'fine':
        fine_keys['lon_fine'] = nc_fine.variables['lon'][:]
        fine_keys['lat_fine'] = nc_fine.variables['lat'][:]

        numdim =  nc_fine.variables[str(variable)].ndim
        if numdim == 4:
            data_fine = nc_fine.variables[str(variable)][DA_step,depth,:,:]
        elif numdim == 3:
            data_fine = nc_fine.variables[str(variable)][DA_step,:,:]

        # Rescale data
        data_fine = data_fine * scale

    return data_fine, data_coarse, coarse_keys, fine_keys


def read_model_point(variable, grid_area, year, month, day, time_stamp, depth, ilon, ilat,DA_switch, coupled):

    ## Main code
    if day < 10:
        day = '0'+str(day)	
    else:
        day = str(day)

    letter = load_data(variable)	

    # Set Scale factor
#    if letter=='secchidat' or letter=='light' or letter=='biodat' or letter=='bendat' \
#           or letter=='kpardat' or letter=='chldat':
#        scale = 1.0/1000.0
#    else:
    scale = 1.0

    # Get paths for assimilation experiment and free run
    path = get_exp_path(coupled, year, month, day, time_stamp, 1)
    if day=='01':
        if DA_switch==1 or DA_switch==2:
            print path
    path_free = get_exp_path(coupled, year, month, day, time_stamp, 0)
    if day=='01' and DA_switch==0:
        print path_free

    # Set file name including path
    timestr = str(year)+str(month)+str(day)+str(time_stamp)
    if DA_switch == 1 or DA_switch==2 : # Background or analysis
        if grid_area == 'all':
            nc_coarse = NetCDFFile(path+str(letter)+'_pdaf_coarse_'+timestr+'.nc4')
            nc_fine   = NetCDFFile(path+str(letter)+'_pdaf_fine_'+timestr+'.nc4')
        elif grid_area == 'fine':
            nc_fine   = NetCDFFile(path+str(letter)+'_pdaf_fine_'+timestr+'.nc4')
            nc_coarse = nc_fine
#            print path+str(letter)+'_pdaf_fine_'+timestr+'.nc4'
        elif grid_area == 'coarse':
            nc_coarse = NetCDFFile(path+str(letter)+'_pdaf_coarse_'+timestr+'.nc4')
            nc_fine = nc_coarse
#            print path+str(letter)+'_pdaf_coarse_'+timestr+'.nc4'
    elif DA_switch == 0: # Free Ens
        if grid_area == 'all':
            nc_coarse = NetCDFFile(path_free+str(letter)+'_pdaf_coarse_'+timestr+'.nc4')
            nc_fine   = NetCDFFile(path_free+str(letter)+'_pdaf_fine_'+timestr+'.nc4')
        elif grid_area == 'fine':
            nc_fine   = NetCDFFile(path_free+str(letter)+'_pdaf_fine_'+timestr+'.nc4')
            nc_coarse = nc_fine
#            print path_free+str(letter)+'_pdaf_fine_'+timestr+'.nc4'
        elif grid_area == 'coarse':
            nc_coarse = NetCDFFile(path_free+str(letter)+'_pdaf_coarse_'+timestr+'.nc4')
            nc_fine = nc_coarse

    # Set file index to read
    if DA_switch == 1: #Background
        # If the output only contains fields (write_variances=.false.):
        #    0: current forecast, 1: analysis, 2: next 12h forecast
        # the 12h forecast is identical to the current forecast in the next output file
        # If write_variances=.true. is set, the output contains variances and fields. Then
        #    0,2,4: variance forcast, analysis, next 12h forecast;
        #    1,3,5: fields: forecast, analysis, next 12h forecast
        DA_step = 0
    elif DA_switch == 2: #Analysis
        DA_step = 1
    elif DA_switch == 0: #Free Ens
        DA_step = 0

    nc_coarse_keys = ['lon', 'lat', str(variable)]
    nc_fine_keys = ['lon', 'lat', str(variable)]

    # Coarse Data
    coarse_keys = {}
    data_coarse = {}
    if grid_area == 'all' or grid_area == 'coarse':
        coarse_keys['lon'] = nc_coarse.variables['lon'][:]
        coarse_keys['lat'] = nc_coarse.variables['lat'][:]
         
        numdim =  nc_coarse.variables[str(variable)].ndim
        if numdim == 4:
            data_coarse = nc_coarse.variables[str(variable)][DA_step,depth,ilat,ilon]
        elif numdim == 3:
            data_coarse = nc_coarse.variables[str(variable)][depth,ilat,ilon]

        # Rescale data
        data_coarse = data_coarse * scale

    # Fine Data
    fine_keys = {}
    data_fine = {}

    if grid_area == 'all' or grid_area == 'fine':
        fine_keys['lon'] = nc_fine.variables['lon'][:]
        fine_keys['lat'] = nc_fine.variables['lat'][:]

        numdim =  nc_fine.variables[str(variable)].ndim
        if numdim == 4:
            data_fine = nc_fine.variables[str(variable)][DA_step,depth,ilat,ilon]
        elif numdim == 3:
            data_fine = nc_fine.variables[str(variable)][depth,ilat,ilon]

        # Rescale data
        data_fine = data_fine * scale

    if grid_area=='fine':
        return data_fine, fine_keys
    else:
        return data_coarse, coarse_keys

def read_model_profile(variable, grid_area, year, month, day, time_stamp, ilon, ilat,DA_switch, coupled):

    ## Main code
    if day < 10:
        day = '0'+str(day)	
    else:
        day = str(day)

    letter = load_data(variable)	

    # Set Scale factor
#    if letter=='secchidat' or letter=='light' or letter=='biodat' or letter=='bendat' \
#           or letter=='kpardat' or letter=='chldat':
#        scale = 1.0/1000.0
#    else:
    scale = 1.0

    # Get paths for assimilation experiment and free run
    path = get_exp_path(coupled, year, month, day, time_stamp, 1)
    if day=='01':
        if DA_switch==1 or DA_switch==2:
            print path
    path_free = get_exp_path(coupled, year, month, day, time_stamp, 0)
    if day=='01' and DA_switch==0:
        print path_free

    # Set file name including path
    timestr = str(year)+str(month)+str(day)+str(time_stamp)
    if DA_switch == 1 or DA_switch==2 : # Background or analysis
        if grid_area == 'all':
            nc_coarse = NetCDFFile(path+str(letter)+'_pdaf_coarse_'+timestr+'.nc4')
            nc_fine   = NetCDFFile(path+str(letter)+'_pdaf_fine_'+timestr+'.nc4')
        elif grid_area == 'fine':
            nc_fine   = NetCDFFile(path+str(letter)+'_pdaf_fine_'+timestr+'.nc4')
            nc_coarse = nc_fine
#            print path+str(letter)+'_pdaf_fine_'+timestr+'.nc4'
        elif grid_area == 'coarse':
            nc_coarse = NetCDFFile(path+str(letter)+'_pdaf_coarse_'+timestr+'.nc4')
            nc_fine = nc_coarse
#            print path+str(letter)+'_pdaf_coarse_'+timestr+'.nc4'
    elif DA_switch == 0: # Free Ens
        if grid_area == 'all':
            nc_coarse = NetCDFFile(path_free+str(letter)+'_pdaf_coarse_'+timestr+'.nc4')
            nc_fine   = NetCDFFile(path_free+str(letter)+'_pdaf_fine_'+timestr+'.nc4')
        elif grid_area == 'fine':
            nc_fine   = NetCDFFile(path_free+str(letter)+'_pdaf_fine_'+timestr+'.nc4')
            nc_coarse = nc_fine
#            print path_free+str(letter)+'_pdaf_fine_'+timestr+'.nc4'
        elif grid_area == 'coarse':
            nc_coarse = NetCDFFile(path_free+str(letter)+'_pdaf_coarse_'+timestr+'.nc4')
            nc_fine = nc_coarse

    nc_depth = NetCDFFile('/home/ollie/lnerger/work/BSH/MERAMO/assim_weak_f0.95_N40/strong-output-201204.h_in_state/'+'z_pdaf_coarse_'+timestr+'.nc4')

    # Set file index to read
    if DA_switch == 1: #Background
        # If the output only contains fields (write_variances=.false.):
        #    0: current forecast, 1: analysis, 2: next 12h forecast
        # the 12h forecast is identical to the current forecast in the next output file
        # If write_variances=.true. is set, the output contains variances and fields. Then
        #    0,2,4: variance forcast, analysis, next 12h forecast;
        #    1,3,5: fields: forecast, analysis, next 12h forecast
        DA_step = 0
    elif DA_switch == 2: #Analysis
        DA_step = 1
    elif DA_switch == 0: #Free Ens
        DA_step = 0

    nc_coarse_keys = ['lon', 'lat', 'h', str(variable)]
    nc_fine_keys = ['lon', 'lat', str(variable)]

    # Coarse Data
    coarse_keys = {}
    data_coarse = {}
    if grid_area == 'all' or grid_area == 'coarse':
        coarse_keys['lon'] = nc_coarse.variables['lon'][:]
        coarse_keys['lat'] = nc_coarse.variables['lat'][:]
        coarse_keys['h'] = nc_depth.variables['h'][DA_step,:,ilat,ilon]
         
        numdim =  nc_coarse.variables[str(variable)].ndim
        if numdim == 4:
            data_coarse = nc_coarse.variables[str(variable)][DA_step,:,ilat,ilon]
        elif numdim == 3:
            data_coarse = nc_coarse.variables[str(variable)][depth,ilat,ilon]

        # Rescale data
        data_coarse = data_coarse * scale

    # Fine Data
    fine_keys = {}
    data_fine = {}

    if grid_area == 'all' or grid_area == 'fine':
        fine_keys['lon'] = nc_fine.variables['lon'][:]
        fine_keys['lat'] = nc_fine.variables['lat'][:]

        numdim =  nc_fine.variables[str(variable)].ndim
        if numdim == 4:
            data_fine = nc_fine.variables[str(variable)][DA_step,depth,ilat,ilon]
        elif numdim == 3:
            data_fine = nc_fine.variables[str(variable)][depth,ilat,ilon]

        # Rescale data
        data_fine = data_fine * scale

    if grid_area=='fine':
        return data_fine, fine_keys
    else:
        return data_coarse, coarse_keys

def read_mld(variable, grid_area, year, month, day, time_stamp, depth, DA_switch, coupled, limit):

    ## Main code
    if day < 10:
        day = '0'+str(day)	
    else:
        day = str(day)

    letter = load_data(variable)	

    # Get paths for assimilation experiment and free run
    path = get_exp_path(coupled, year, month, day, time_stamp, 1)
    if day=='01':
        if DA_switch==1 or DA_switch==2:
            print path
    path_free = get_exp_path(coupled, year, month, day, time_stamp, 0)
    if day=='01' and DA_switch==0:
        print path_free

    # Set file name including path
    timestr = str(year)+str(month)+str(day)+str(time_stamp)
    if DA_switch == 1 or DA_switch==2 : # Background or analysis
        if grid_area == 'all':
            nc_coarse = NetCDFFile(path+str(letter)+'_pdaf_coarse_'+timestr+'.nc4')
            nc_fine   = NetCDFFile(path+str(letter)+'_pdaf_fine_'+timestr+'.nc4')
        elif grid_area == 'fine':
            nc_fine   = NetCDFFile(path+str(letter)+'_pdaf_fine_'+timestr+'.nc4')
            nc_coarse = nc_fine
        elif grid_area == 'coarse':
            nc_coarse = NetCDFFile(path+str(letter)+'_pdaf_coarse_'+timestr+'.nc4')
            nc_fine = nc_coarse
    elif DA_switch == 0: # Free Ens
        if grid_area == 'all':
            nc_coarse = NetCDFFile(path_free+str(letter)+'_pdaf_coarse_'+timestr+'.nc4')
            nc_fine   = NetCDFFile(path_free+str(letter)+'_pdaf_fine_'+timestr+'.nc4')
        elif grid_area == 'fine':
            nc_fine   = NetCDFFile(path_free+str(letter)+'_pdaf_fine_'+timestr+'.nc4')
            nc_coarse = nc_fine
        elif grid_area == 'coarse':
            nc_coarse = NetCDFFile(path_free+str(letter)+'_pdaf_coarse_'+timestr+'.nc4')
            nc_fine = nc_coarse

    # Set file index to read
    if DA_switch == 1: #Background
        # If the output only contains fields (write_variances=.false.):
        #    0: current forecast, 1: analysis, 2: next 12h forecast
        # the 12h forecast is identical to the current forecast in the next output file
        # If write_variances=.true. is set, the output contains variances and fields. Then
        #    0,2,4: variance forcast, analysis, next 12h forecast;
        #    1,3,5: fields: forecast, analysis, next 12h forecast
        DA_step = 0
    elif DA_switch == 2: #Analysis
        DA_step = 1
    elif DA_switch == 0: #Free Ens
        DA_step = 0

    nc_coarse_keys = ['lon', 'lat', 'depth',str(variable)]
    nc_fine_keys = ['lon', 'lat', 'depth', str(variable)]

    # Coarse Data
    coarse_keys = {}
    data_coarse = {}
    depth_key_c = {}
    temp = {}
    if grid_area == 'all' or grid_area == 'coarse':
        coarse_keys['lon_coarse'] = nc_coarse.variables['lon'][:]
        coarse_keys['lat_coarse'] = nc_coarse.variables['lat'][:]
        depth_key_c['depth'] = nc_coarse.variables['depth'][:]
         
        numdim =  nc_coarse.variables[str(variable)].ndim
        if numdim == 4:
            temp = nc_coarse.variables[str(variable)][DA_step,:,:,:]
            data_coarse = nc_coarse.variables[str(variable)][DA_step,0,:,:]
        elif numdim == 3:
            temp = nc_coarse.variables[str(variable)][:,:,:]

        nlat = len(temp[1,:,1])
        nlon = len(temp[1,1,:])
        ndep = len(temp[:,1,1])
        print 'Coarse dimensions:',ndep, nlat, nlon

        for ilat in range(nlat):
            for ilon in range(nlon):
                stemp = temp[0,ilat,ilon]
                if stemp>-2.0 or stemp<40.0:
                    data_coarse[ilat, ilon] = depth_key_c['depth'][0]
                    for idep in range(1,ndep):
                        if np.abs(stemp - temp[idep,ilat,ilon]) <= limit:
                            data_coarse[ilat, ilon] = depth_key_c['depth'][idep]

    # Fine Data
    fine_keys = {}
    data_fine = {}
    depth_key_f = {}
    temp_f = {}
    if grid_area == 'all' or grid_area == 'fine':
        fine_keys['lon_fine'] = nc_fine.variables['lon'][:]
        fine_keys['lat_fine'] = nc_fine.variables['lat'][:]
        depth_key_f['depth'] = nc_fine.variables['depth'][:]

        numdim =  nc_fine.variables[str(variable)].ndim
        if numdim == 4:
            temp_f = nc_fine.variables[str(variable)][DA_step,:,:,:]
            data_fine = nc_fine.variables[str(variable)][DA_step,0,:,:]
        elif numdim == 3:
            temp_f = nc_fine.variables[str(variable)][:,:,:]

        nlat = len(temp_f[1,:,1])
        nlon = len(temp_f[1,1,:])
        ndep = len(temp_f[:,1,1])
        print 'Fine dimensions:',ndep, nlat, nlon

        for ilat in range(nlat):
            for ilon in range(nlon):
                stemp = temp_f[0,ilat,ilon]
                if stemp>-2.0 or stemp<40.0:
                    data_fine[ilat, ilon] = depth_key_f['depth'][0]
                    for idep in range(ndep):
                        if np.abs(stemp - temp_f[idep,ilat,ilon]) <= limit:
                            data_fine[ilat, ilon] = depth_key_f['depth'][idep]

    return data_fine, data_coarse, coarse_keys, fine_keys


def rmse(predictions, targets):
    if np.nansum(predictions) == np.nan or np.nansum(targets) == np.nan:
	rmse_val = np.nan
	return rmse_val
    predictions = predictions[~np.isnan(predictions)]
    targets = targets[~np.isnan(targets)]
    differences = predictions - targets                       #the DIFFERENCEs.
    differences_squared = differences ** 2                    #the SQUAREs of ^
    mean_of_differences_squared = differences_squared.mean()  #the MEAN of ^
    rmse_val = np.sqrt(mean_of_differences_squared)           #ROOT of ^

    return rmse_val     #get the ^

def rmse_bio(predictions, targets, inflate):
    if inflate == 'on':
        targets = targets*1000.0
    predictions = predictions[~np.isnan(predictions)]
    targets = targets[~np.isnan(targets)]
    differences = predictions - targets                       #the DIFFERENCEs.
    differences_squared = differences ** 2                    #the SQUAREs of ^
    mean_of_differences_squared = differences_squared.mean()  #the MEAN of ^
    rmse_val = np.sqrt(mean_of_differences_squared)           #ROOT of ^

    return rmse_val     #get the ^


def ssum_diff(predictions, targets, checknan):
    ssum = 0.0
    if checknan==1:
        predictions = predictions[~np.isnan(predictions)]
        targets = targets[~np.isnan(targets)]

    for i in range(0,len(predictions)):
        diff = predictions[i] - targets[i]                       #the DIFFERENCEs.
        diff_squared = diff ** 2                    #the SQUAREs of ^
        ssum = ssum + diff_squared

    return ssum   


def ssum_sst_diff(predictions, targets, checknan):
    ssum = 0.0
#    if checknan==1:
#        predictions = predictions[~np.isnan(predictions)]
#        targets = targets[~np.isnan(targets)]

    cntvals = 0
    for i in range(0,len(predictions)):
        if ~np.isnan(predictions[i]) & ~np.isnan(targets[i]):
            diff = predictions[i] - targets[i]                       #the DIFFERENCEs.
            diff_squared = diff ** 2                    #the SQUAREs of ^
            ssum = ssum + diff_squared
            cntvals = cntvals + 1

    return ssum , cntvals 

def sum_diff(predictions, targets, checknan):
    sumall = 0.0
    if checknan==1:
        predictions = predictions[~np.isnan(predictions)]
        targets = targets[~np.isnan(targets)]

    for i in range(0,len(predictions)):
        diff = predictions[i] - targets[i]                       #the DIFFERENCEs.
        sumall = sumall + diff

    return sumall

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def rms_time_loop(length, daytime, depth, variable, grid_area, arg1, coupled, bio_switch, mat_data, mat_lat, \
               mat_lon, mat_time, SST_switch, year, month_aux, rms_or_bias):

    # Loop over days to compute RMS error with respect to in situ validation data

    lat_offset = 94   # Off set in latitude from end of array
    lat_inset =  158  # End if fine grid in latitude from end of array
    lon_offset = 123 # off set to fine grid
    lon_inset = 228  # end of fine grid

    k = 1 # day
    total_vali_points = 0
    sst_vali_points = 0
    day_vali_points = []
    sum_squared = 0
    bias_all = 0
    sst_sum_squared = 0
    print 'coupled', coupled

    val_insitu = []
    val_model  = []

    for j in range(1,length+1):
        mat = []
        lat = []
        lon = []
        day_vali_points.append(0)
        if daytime == 'both':
            if (j % 2 == 0): #even
                data_fine, data_coarse, coarse_keys, fine_keys \
                           = read_model(variable, grid_area, year, month_aux, k, 'pm', str(depth),arg1, coupled)
                if SST_switch == 'on':
                    SST_fine, SST_coarse = get_sst(month_aux, k, 'pm')
            else:
                data_fine, data_coarse, coarse_keys, fine_keys = \
                           read_model(variable, grid_area, year, month_aux, k, 'am', str(depth),arg1, coupled)
                if SST_switch == 'on':
                    SST_fine, SST_coarse = get_sst(month_aux, k, 'am')
            for i in range(len(mat_data)):
                if (k-1) == mat_time[i,]:
                    mat.append(mat_data[i,])
                    lat.append(mat_lat[i,])
                    lon.append(mat_lon[i,])
        elif daytime == 'am':
            data_fine, data_coarse, coarse_keys, fine_keys = \
                       read_model(variable, grid_area, year, month_aux, j, 'am', str(depth),arg1, coupled)
            if SST_switch == 'on':
                SST_fine, SST_coarse = get_sst(month_aux, j, 'am')
            for i in range(len(mat_data)):
                if (j-1) == mat_time[i,]:
                    mat.append(mat_data[i,])
                    lat.append(mat_lat[i,])
                    lon.append(mat_lon[i,])
        elif daytime == 'pm':
            data_fine, data_coarse, coarse_keys, fine_keys = \
                       read_model(variable, grid_area, year, month_aux, j, 'pm', str(depth),arg1, coupled)
            if SST_switch == 'on':
                SST_fine, SST_coarse = get_sst(month_aux, j, 'pm')
            for i in range(len(mat_data)):
                if (j-1) == mat_time[i,]:
                    mat.append(mat_data[i,])
                    lat.append(mat_lat[i,])
                    lon.append(mat_lon[i,])



        if len(mat) == 0:
            if daytime == 'both' and (j % 2 == 0):
                k = k + 1
            continue

        mat = np.array(mat)
        mat = mat.reshape((len(mat),1))
        data = np.empty(shape=(len(mat),1))
        sst = np.empty(shape=(len(mat),1))
        cnt = 0
        for m, n in zip(lat, lon):
            if grid_area == 'coarse':
                if data_coarse[-(m+1), (n)] > -100 and data_coarse[-(m+1), (n)] < 1e6:
                    if np.isnan(data_coarse[-(m+1), (n)]) ==0:
                        total_vali_points = total_vali_points + 1
                        day_vali_points[j-1] = day_vali_points[j-1] + 1
                        data[cnt,0] = data_coarse[-(m+1), (n)]
                        if SST_switch == 'on':
                            sst[cnt,0] = SST_coarse[-(m+1), (n)]
                        else:
                            sst[cnt,0] = np.nan
                    else:
                        data[cnt,0] = np.nan
                else:
                    data[cnt,0] = np.nan
                    sst[cnt,0] = np.nan
            elif grid_area == 'fine':
                if m < lat_inset and m > lat_offset and n > lon_offset and n < lon_inset:
                    if np.isnan(data_fine[-6*(m-lat_offset)-2, 6*(n-lon_offset)-3]) == 0:
                        total_vali_points = total_vali_points + 1
                        day_vali_points[j-1] = day_vali_points[j-1] + 1
                        data[cnt,0] = data_fine[-6*(m-lat_offset)-2, 6*(n-lon_offset)-3]
                        if SST_switch == 'on':
                            sst[cnt,0] = SST_fine[-6*(m-lat_offset)-2, 6*(n-lon_offset)-3]
                        else:
                            sst[cnt,0] = np.nan
                    else:
                        data[cnt,0] = np.nan
                else:
                    data[cnt,0] = np.nan
                    sst[cnt,0] = np.nan
            elif grid_area == 'all':
                if m < lat_inset and m > lat_offset and n > lon_offset and n < lon_inset:
                    if np.isnan(data_fine[-6*(m-lat_offset)-2, 6*(n-lon_offset)-3]) == 0:
                        total_vali_points = total_vali_points + 1
                        day_vali_points[j-1] = day_vali_points[j-1] + 1
                        data[cnt,0] = data_fine[-6*(m-lat_offset)-2, 6*(n-lon_offset)-3]
                        if SST_switch == 'on':
                            sst[cnt,0] = SST_fine[-6*(m-lat_offset)-2, 6*(n-lon_offset)-3]
                        else:
                            sst[cnt,0] = np.nan
                    else:
                        data[cnt,0] = np.nan
                        sst[cnt,0] = np.nan
                elif data_coarse[-(m+1), (n)] > -100 and data_coarse[-(m+1), (n)] < 1e6:
                    if np.isnan(data_coarse[-(m+1), (n)]) ==0:
                        total_vali_points = total_vali_points + 1
                        day_vali_points[j-1] = day_vali_points[j-1] + 1
                        data[cnt,0] = data_coarse[-(m+1), (n)]
                        if SST_switch == 'on':
                            sst[cnt,0] = SST_coarse[(m+1), (n)]
                        else:
                            sst[cnt,0] = np.nan
                    else:
                        data[cnt,0] = no.nan
                        sst[cnt,0] = np.nan
                else:
                    data[cnt,0] = np.nan
                    sst[cnt,0] = np.nan
            cnt = cnt + 1

        for i in range(len(lat)):
            if SST_switch == 'on':
                if np.isnan(data[i,0]) == 1 or np.isnan(sst[i,0]) == 1:
                    mat[i,0] = np.nan
                    data[i,0] = np.nan
                    sst[i,0] = np.nan
            if SST_switch == 'off':
                if np.isnan(data[i,0]) == 1:
                    mat[i,0] = np.nan
                    data[i,0] = np.nan
                    sst[i,0] = np.nan  

        # Collect data arrays
        for i in range(len(mat[:,0])):
            val_insitu.append(mat[i,0])
            val_model.append(data[i,0])

        ssum_day = ssum_diff(data, mat, 1)
        
        sum_day = sum_diff(data, mat, 1)
        sum_squared = sum_squared + ssum_day
        bias_all = bias_all + sum_day
        if SST_switch == 'on':
            ssum_day, ndata = ssum_sst_diff(data, mat, 1 )
            sst_vali_points = sst_vali_points + ndata
            sst_sum_squred = sst_sum_squared + ssum_day

        if daytime == 'both' and (j % 2 == 0):
            k = k + 1
        if daytime =='am' or daytime=='pm':
            k = k + 1

    print total_vali_points, sst_vali_points
    if total_vali_points>0:
        rmse_mon = np.sqrt(sum_squared / total_vali_points)
        if SST_switch == 'on':
            rmse_mon_sst = np.sqrt(sst_sum_squared / sst_vali_points)
        else:
            rmse_mon_sst = 0
        bias_mon = bias_all / total_vali_points
        bias_mon_sst = bias_all / total_vali_points
    else:
        rmse_mon = 0
        rmse_mon_sst = 0
        bias_mon = 0
        bias_mon_sst = 0

    if SST_switch=='on':
        total_vali_points = sst_vali_points
    if (rms_or_bias==0) :
        return rmse_mon, rmse_mon_sst, total_vali_points, val_insitu, val_model
    else :
        return bias_mon, bias_mon_sst, total_vali_points, val_insitu, val_model


def get_rmse(length, days, daytime, mat_file, depth, variable, grid_area, coupled, year, month, \
                   SST_switch, bio_switch, forana, sea):

    print 'Sea region: ', sea

    # Set flag to compute bias
    rms_or_bias = 0

    print ' -Read Mat'
    mat_data, mat_lat, mat_lon, mat_time = read_mat_sea(mat_file, days, depth, sea, variable)
    
    if forana=='for':
        print ' -DA Background'
        rmse_mon, rmse_mon_sst, val_points, _, _ \
                  = rms_time_loop(length, daytime, depth, variable, grid_area, 1, coupled, bio_switch, \
                               mat_data, mat_lat, mat_lon, mat_time, SST_switch, year, month, rms_or_bias)
        print 'val_points for', val_points
    elif forana=='ana':
        print ' -DA Analysis'
        rmse_mon, rmse_mon_sst, val_points, _, _ \
                  = rms_time_loop(length, daytime, depth, variable, grid_area, 2, coupled, bio_switch, \
                               mat_data, mat_lat, mat_lon, mat_time, SST_switch, year, month, rms_or_bias)
        print 'val_points ana', val_points
    elif forana=='fre':
        print ' -Free Ensemble'
        rmse_mon, rmse_mon_sst, val_points, _, _ \
                  = rms_time_loop(length, daytime, depth, variable, grid_area, 0, coupled, bio_switch, \
                               mat_data, mat_lat, mat_lon, mat_time, SST_switch, year, month, rms_or_bias)
        print 'val_points free', val_points

    return rmse_mon, rmse_mon_sst, val_points, mat_data, mat_lon, mat_lat


def get_rmse_sst_vali(length, days, daytime, mat_file, depth, variable, grid_area, coupled, year, month, \
                   SST_switch, bio_switch, forana, sea):

    print 'Sea region: ', sea

    # Set flag to compute bias
    rms_or_bias = 0

    print ' -Read Mat'
    mat_data, mat_lat, mat_lon, mat_time = read_mat_sea(mat_file, days, depth, sea, variable)
    
    if forana=='for':
        print ' -DA Background'
        rmse_mon, rmse_mon_sst, val_points, _, _ \
                  = rms_sst_time_loop(length, daytime, depth, variable, grid_area, 1, coupled, bio_switch, \
                               mat_data, mat_lat, mat_lon, mat_time, SST_switch, year, month, rms_or_bias)
        print 'val_points for', val_points
    elif forana=='ana':
        print ' -DA Analysis'
        rmse_mon, rmse_mon_sst, val_points, _, _ \
                  = rms_sst_time_loop(length, daytime, depth, variable, grid_area, 2, coupled, bio_switch, \
                               mat_data, mat_lat, mat_lon, mat_time, SST_switch, year, month, rms_or_bias)
        print 'val_points ana', val_points
    elif forana=='fre':
        print ' -Free Ensemble'
        rmse_mon, rmse_mon_sst, val_points, _, _ \
                  = rms_sst_time_loop(length, daytime, depth, variable, grid_area, 0, coupled, bio_switch, \
                               mat_data, mat_lat, mat_lon, mat_time, SST_switch, year, month, rms_or_bias)
        print 'val_points free', val_points

    return rmse_mon, rmse_mon_sst, val_points, mat_data, mat_lon, mat_lat



def get_bias_vali(length, days, daytime, mat_file, depth, variable, grid_area, coupled, year, month, \
                   SST_switch, bio_switch, forana, sea):

    print 'Sea region: ', sea

    # Set flag to compute bias
    rms_or_bias = 1

    print ' -Read Mat'
    mat_data, mat_lat, mat_lon, mat_time = read_mat_sea(mat_file, days, depth, sea, variable)
    
    if forana=='for':
        print ' -DA Background'
        bias_mon, bias_mon_sst, val_points, _, _ \
                  = rms_time_loop(length, daytime, depth, variable, grid_area, 1, coupled, bio_switch, \
                               mat_data, mat_lat, mat_lon, mat_time, SST_switch, year, month, rms_or_bias)
        print 'val_points for', val_points
    elif forana=='ana':
        print ' -DA Analysis'
        bias_mon, bias_mon_sst, val_points, _, _ \
                  = rms_time_loop(length, daytime, depth, variable, grid_area, 2, coupled, bio_switch, \
                               mat_data, mat_lat, mat_lon, mat_time, SST_switch, year, month, rms_or_bias)
        print 'val_points ana', val_points
    elif forana=='fre':
        print ' -Free Ensemble'
        bias_mon, bias_mon_sst, val_points, _, _ \
                  = rms_time_loop(length, daytime, depth, variable, grid_area, 0, coupled, bio_switch, \
                               mat_data, mat_lat, mat_lon, mat_time, SST_switch, year, month, rms_or_bias)
        print 'val_points free', val_points

    return bias_mon, bias_mon_sst, val_points, mat_data, mat_lon, mat_lat


def get_rms_sst3(length, year, days, daytime, depth, variable, grid_area, coupled, month_aux):

    # Compute daily RMS errors for SST over time
    # For the years before 2017 NOAA satellite data is read, from year 2017 we read Sentinel 3A data

    rmse_sst_vals_fr = np.empty(shape=(length,))
    rmse_sst_vals_fo = np.empty(shape=(length,))
    rmse_sst_vals_a  = np.empty(shape=(length,))
    sst_vals = np.empty(shape=(length,))

    k = 1 # day
    total_vali_points = 0
    for j in range(1,length+1):

        # Read model data and SST observations
        if daytime == 'both':
            if (j % 2 == 0): #even
                # Free run
                data_f_fr, data_c_fr, _, _ \
                    = read_model(variable, grid_area, year, month_aux, k, 'pm', str(depth), 0, coupled)
                data_f_fo, data_c_fo, _, _ \
                    = read_model(variable, grid_area, year, month_aux, k, 'pm', str(depth), 1, coupled)
                data_f_a, data_c_a, keys_c, keys_f \
                    = read_model(variable, grid_area, year, month_aux, k, 'pm', str(depth), 2, coupled)
                if year < 2017:
                    SST_fine, SST_coarse = get_sst(month_aux, k, 'pm')
                else:
                    SST_fine, SST_coarse, _,_,_,_ = get_s3a_sst(month_aux, k, 'pm')
            else:
                data_f_fr, data_c_fr, _, _ \
                    = read_model(variable, grid_area, year, month_aux, k, 'am', str(depth), 0, coupled)
                data_f_fo, data_c_fo, _, _ \
                    = read_model(variable, grid_area, year, month_aux, k, 'am', str(depth), 1, coupled)
                data_f_a, data_c_a, keys_c, keys_f \
                    = read_model(variable, grid_area, year, month_aux, k, 'am', str(depth), 2, coupled)
                if year < 2017:
                    SST_fine, SST_coarse = get_sst(month_aux, k, 'am')
                else:
                    SST_fine, SST_coarse, _,_,_,_ = get_s3a_sst(month_aux, k, 'am')
        elif daytime == 'am':
            # Free
            data_f_fr, data_c_fr, _, _ \
                = read_model(variable, grid_area, year, month_aux, j, 'am', str(depth),0, coupled)
            # Forecast
            data_f_fo, data_c_fo, _, _ \
                = read_model(variable, grid_area, year, month_aux, j, 'am', str(depth),1, coupled)
            # Analysis
            data_f_a, data_c_a, keys_c, keys_f \
                = read_model(variable, grid_area, year, month_aux, j, 'am', str(depth),2, coupled)
            if year < 2017:
                SST_fine, SST_coarse = get_sst(month_aux, j, 'am')
            else:
                SST_fine, SST_coarse, _,_,_,_ = get_s3a_sst(month_aux, j, 'am')
        elif daytime == 'pm':
            # Free
            data_f_fr, data_c_fr, _, _ \
                = read_model(variable, grid_area, year, month_aux, j, 'pm', str(depth),0, coupled)
            # Forecast
            data_f_fo, data_c_fo, _, _ \
                = read_model(variable, grid_area, year, month_aux, j, 'pm', str(depth),1, coupled)
            # Analysis
            data_f_a, data_c_a, keys_c, keys_f \
                = read_model(variable, grid_area, year, month_aux, j, 'pm', str(depth),2, coupled)
            if year < 2017:
                SST_fine, SST_coarse = get_sst(month_aux, j, 'pm')
            else:
                SST_fine, SST_coarse, _,_,_,_ = get_s3a_sst(month_aux, j, 'pm')

        data_fr = []
        data_fo = []
        data_a = []
        sst = []
        cnt_day = 0
        
        # Get collocated data for fine grid
        if grid_area == 'fine' or grid_area == 'all':
            for m in range(1,len(keys_f['lat_fine'])):
                for n in range(1,len(keys_f['lon_fine'])):
                    if data_f_fr[m, n] > -100 and data_f_fr[m, n] < 1e6 \
                           and SST_fine[m, n]>-100 and  SST_fine[m, n] < 1e6:
                        total_vali_points = total_vali_points + 1
                        data_fr.append(data_f_fr[(m), (n)])
                        data_fo.append(data_f_fo[(m), (n)])
                        data_a.append(data_f_a[(m), (n)])
                        sst.append(SST_fine[(m), (n)])
                        cnt_day = cnt_day + 1

        # Get collocated data for coarse grid
        if grid_area == 'coarse' or grid_area == 'all':
            for m in range(1,len(keys_c['lat_coarse'])):
                for n in range(1,len(keys_c['lon_coarse'])):
                    if data_c_fr[m, n] > -100 and data_c_fr[m, n] < 1e6 \
                           and SST_coarse[m, n]>-100 and  SST_coarse[m, n] < 1e6:
                        total_vali_points = total_vali_points + 1
                        data_fr.append(data_c_fr[(m), (n)])
                        data_fo.append(data_c_fo[(m), (n)])
                        data_a.append(data_c_a[(m), (n)])
                        sst.append(SST_coarse[(m), (n)])
                        cnt_day = cnt_day + 1


        # Compute RMS error
        ssum_sst_fr = ssum_diff(data_fr, sst, 0)
        ssum_sst_fo = ssum_diff(data_fo, sst, 0)
        ssum_sst_a  = ssum_diff(data_a, sst, 0)

        if cnt_day>0:
            rmse_sst_vals_fr[j-1] = np.sqrt(ssum_sst_fr / cnt_day)
            rmse_sst_vals_fo[j-1] = np.sqrt(ssum_sst_fo / cnt_day)
            rmse_sst_vals_a[j-1]  = np.sqrt(ssum_sst_a / cnt_day)
        else:
            rmse_sst_vals_fr[j-1] = np.nan
            rmse_sst_vals_fo[j-1] = np.nan
            rmse_sst_vals_a[j-1] = np.nan

        print 'step',j, 'day', k, ' collocation points: ', cnt_day, 'RMS fr,for,ana: ',rmse_sst_vals_fr[j-1], rmse_sst_vals_fo[j-1], rmse_sst_vals_a[j-1]

        if daytime == 'both' and (j % 2 == 0):
            k = k + 1
        if daytime == 'am' or daytime == 'pm':
            k=k+1

    print 'val_points monthly avg.', total_vali_points / length

    return rmse_sst_vals_fr, rmse_sst_vals_fo, rmse_sst_vals_a, total_vali_points


def plot_rms_time_series(length, free, fcst, ana, title, fname, save, m_Names, m_days):
    # Plot a time series of RMS errors for forecast, analysis, and free run
    plt.figure
    plt.title(str(title), fontsize=16)
    plt.plot(ana, c='g', label='Analysis') # plotting by columns
    plt.plot(fcst, c='b', label='Forecast')
    plt.plot(free, c='k', label='Free')
#    plt.plot(free, c='k', linewidth=2.0, label='Modell')
#    plt.plot(ana, c='g', linewidth=2.0, label='Assim.') # plotting by columns
    tickloc,ticknam = plt.xticks()
    plt.xticks(m_days,m_Names, fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylabel('RMSe ($^oC$)', fontsize=16)
    plt.xlim([0, len(ana)])
    plt.ylim([0.2, 2.2])
    plt.legend(loc='upper left')
#    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    if save==1:
#	    plt.savefig('figs/'+fname, bbox_extra_artists=(lgd,), bbox_inches='tight')
	    plt.savefig(fname)
    else:
	    plt.show()


def plot_2bars(A, ens, my_xticks, val, grid, var, var_unit, fname, save, sea, rms_or_bias):
    x = np.arange(len(A))
    print x
    print A
    print val
    print ens
    fig, ax1 = plt.subplots()
    # Bar
    bar_width = 0.25
    opacity = 0.8
    ax1.bar(x*0.75+0*bar_width,A,bar_width,alpha=opacity,color='r',label='DA')
    ax1.bar(x*0.75+1*bar_width,ens,bar_width,alpha=opacity,color='b',label='Free Ensemble')  
    ax1.set_xticks(x*0.75 + bar_width)

    ax1.set_xlabel('2012', fontsize=18)
    ax1.set_xticklabels(my_xticks, fontsize=18)
    plt.xlim(0, len(x)*2*bar_width+bar_width*(len(x)-1))# + 1.5*bar_width)
    if rms_or_bias==0:
        ax1.set_ylabel('RMSE ['+var_unit+']', fontsize=18)
    else:
        ax1.set_ylabel('bias ['+var_unit+']', fontsize=18)
    for item in (ax1.get_yticklabels()):
        item.set_fontsize(18)

    ax2 = ax1.twiny()
    ax2.set_xlabel('validation points', fontsize=12)
    ax2.set_xticks(0.75*x + bar_width)
    ax2.set_xticklabels(val, fontsize=18)
    plt.xlim(0, len(x)*2*bar_width+bar_width*(len(x)-1))# + 1.5*bar_width)

    fig.tight_layout()
    lgd = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=18)
    if sea=='North':
        strsea = ' in North Sea'
    elif sea=='Baltic':
        strsea = ' in Baltic Sea'
    elif sea=='all':
        strsea = ''
    if rms_or_bias==0:
        rmsstr = 'RMS errors'
    else:
        rmsstr = 'Bias'
    figure_title = str(var)+": "+rmsstr+" on "+str(grid)+" grid"+strsea
    plt.title(figure_title, y=1.07, fontsize=18)
    if save==1:
	    plt.savefig(str(fname), bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
	    plt.show()


def plot_3bars(rmsw, rmss, rmsf, my_xticks, val, grid, var, fname, save, sea):
    x = np.arange(len(rmsw))
    print rmsw
    print rmss
    print rmsf
    fig, ax1 = plt.subplots()
    # Bar
    bar_width = 0.25
    opacity = 0.8
    ax1.bar(x+0*bar_width,rmsw,bar_width,alpha=opacity,color='r',label='DA Weakly')
    ax1.bar(x+1*bar_width,rmss,bar_width,alpha=opacity,color='m',label='DA Strongly')
    ax1.bar(x+2*bar_width,rmsf,bar_width,alpha=opacity,color='b',label='Free Ensemble')  
    ax1.set_xticks(x + 1.5*bar_width)
    
    ax1.set_xlabel('2012')
    ax1.set_xticklabels(my_xticks)
    plt.xlim(0, len(x)*3*bar_width+bar_width*(len(x)-1))# + 1.5*bar_width)
    ax1.set_ylabel('RMSE', fontsize=16)
    ax1.ylabel(fonsize=16)
    
    ax2 = ax1.twiny()
    ax2.set_xlabel('validation points',fontsize=16)
    ax2.set_xticks(x + 1.5*bar_width)
    ax2.set_xticklabels(val)
    #plt.xlim(0, len(x) + 1.5*bar_width)
    plt.xlim(0, len(x)*3*bar_width+bar_width*(len(x)-1))# + 1.5*bar_width)
    
    fig.tight_layout()
    lgd = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if sea=='North':
        strsea = ' in North Sea'
    elif sea=='Baltic':
        strsea = ' in Baltic Sea'
    elif sea=='all':
        strsea = ''     
    figure_title = str(var)+": RMS errors on "+str(grid)+" grid"+strsea
    plt.title(figure_title, y=1.1, fonztsize=16)
    if save==1:
	    plt.savefig(str(fname), bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
	    plt.show()


def plot_map(data_coarse, lat_c, lon_c, varnum, domain, \
             strcmap, minmax, save, title, fname, plotlog, plotcb):
    # Scrip to plot a field as a map

    varstr, mat_var, Variable, var_unit = var_names(varnum)

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

    # Parameters/Variables for plots
    vmin = np.nanmin(data_coarse) #; fine_min = np.min(data_fine)
    vmax = np.nanmax(data_coarse) #; fine_max = np.max(data_fine)
    freeclim=1

    print 'Min/Max data values: ', vmin, ', ', vmax
    if minmax[0]!=minmax[1]:
        vmin = minmax[0]
        vmax = minmax[1]
        freeclim=0
        print 'Limit color range to: ', minmax[0], 'to', minmax[1]

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

    cmap = plt.cm.get_cmap(strcmap)
#    cmap1 = plt.cm.get_cmap('terrain_r')
    clev = np.linspace(vmin, vmax, 100)

    # Data - Coarse
    lons,lats = np.meshgrid(lon_c, lat_c)
    x, y = m(lons,lats)
    if vmin < 0.0:
        if freeclim==1:
            cs1 = m.contourf(x, y, data_coarse, cmap=cmap, vmin=vmin, vmax=vmax, \
                             levels=clev, norm=MidpointNormalize(midpoint=0.0))
        else:
            cs1 = m.pcolor(x, y, data_coarse, cmap=cmap, vmin=vmin, vmax=vmax, \
                           norm=MidpointNormalize(midpoint=0.0))
    else:
#        if freeclim==1:
#            cs1 = m.contourf(x, y, data_coarse, cmap=cmap, vmin=vmin, vmax=vmax, levels=clev)
#        else:
       if plotlog==1:
           print 'Plot log values'
           cs1 = m.pcolor(x, y, data_coarse, cmap=cmap, vmin=vmin, vmax=vmax, norm=colors.LogNorm(vmin=vmin, vmax=vmax))
       else: 
           cs1 = m.pcolor(x, y, data_coarse, cmap=cmap, vmin=vmin, vmax=vmax)

    # Save/Show Map
    plt.title(title,fontsize=16)
    if plotcb==1:
        #    tick_locator = ticker.MaxNLocator(nbins=5)
        cbar = m.colorbar(cs1)
        #    cbar.locator = tick_locator
        #    cbar.update_ticks()
        cbar.set_label(var_unit,fontsize=18)
        cbar.ax.tick_params(labelsize=16) 

    if save==1:
	    plt.savefig(fname, dpi=200)
    else:
	    plt.show()


def plot_map_baltic(data_coarse, data_fine, lat_c, lon_c, lat_f, lon_f, varnum, domain, \
             strcmap, minmax, save, title, fname):

    # Scrip to plot a field as a map

    varstr, mat_var, Variable, var_unit = var_names(varnum)

    if domain=='no':
            llon=9.0 #lon_c[0]
            llat=53.5 #lat_c[-1]
            ulon=lon_c[-1]
            ulat=lat_c[0]
            res = 5.0  # resolution for meridians & parallels
            bmres = 'l'
    else:
            llon=lon_f[0]
            llat=lat_f[-1]
            ulon=lon_f[-1]
            ulat=lat_f[0]
            res = 1.0  # resolution for meridians & parallels
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

    # Parameters/Variables for plots
    coarse_min = np.min(data_coarse); fine_min = np.min(data_fine)
    coarse_max = np.max(data_coarse); fine_max = np.max(data_fine)
    if domain=='no':
            vmin = min(coarse_min, fine_min)
            vmax = max(coarse_max, fine_max)
            freeclim=1
    else:
	    vmin = fine_min
            vmax = fine_max
            freeclim=1
    print 'Min/Max data values: ', vmin, ', ', vmax
    if minmax[0]!=minmax[1]:
        vmin = minmax[0]
        vmax = minmax[1]
        freeclim=0
        print 'Limit color range to: ', minmax[0], 'to', minmax[1]

    cmap = plt.cm.get_cmap(strcmap)
#    cmap1 = plt.cm.get_cmap('terrain_r')
    clev = np.linspace(vmin, vmax, 100)

    # Data - Coarse
    if domain=='no':
        lons,lats = np.meshgrid(lon_c, lat_c)
        x, y = m(lons,lats)
        if vmin < 0.0:
            if freeclim==1:
                cs1 = m.contourf(x, y, data_coarse, cmap=cmap, vmin=vmin, vmax=vmax, \
                                 levels=clev, norm=MidpointNormalize(midpoint=0.0))
            else:
                cs1 = m.pcolor(x, y, data_coarse, cmap=cmap, vmin=vmin, vmax=vmax, \
                               norm=MidpointNormalize(midpoint=0.0))
        else:
            if freeclim==1:
                cs1 = m.contourf(x, y, data_coarse, cmap=cmap, vmin=vmin, vmax=vmax, levels=clev)
            else:
                cs1 = m.pcolor(x, y, data_coarse, cmap=cmap, vmin=vmin, vmax=vmax)
		    
    # Data - Fine
    lons,lats = np.meshgrid(lon_f, lat_f)
    x, y = m(lons,lats)
    if 1==1:   # Use this to exclude the fine mesh when plotting the coarse
        if vmin < 0.0:
            if freeclim==1:
                cs2 = m.contourf(x, y, data_fine, cmap=cmap, vmin=vmin, vmax=vmax, \
                                 levels=clev, norm=MidpointNormalize(midpoint=0.0))
            else:
                cs2 = m.pcolor(x, y, data_fine, cmap=cmap, vmin=vmin, vmax=vmax, \
                               norm=MidpointNormalize(midpoint=0.0))
        else:
            if freeclim==1:
                cs2 = m.contourf(x, y, data_fine, cmap=cmap, vmin=vmin, vmax=vmax, levels=clev)
            else:
                cs2 = m.pcolor(x, y, data_fine, cmap=cmap, vmin=vmin, vmax=vmax)

    # Save/Show Map
    plt.title(title,fontsize=16)
    tick_locator = ticker.MaxNLocator(nbins=5)
    if domain=='no':
            cbar = m.colorbar(cs1)
    else:
            cbar = m.colorbar(cs2)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label(var_unit,fontsize=18)
    cbar.ax.tick_params(labelsize=16) 

    if save==1:
	    plt.savefig(fname, dpi=fig.dpi)
    else:
	    plt.show()

def plot_map_dataandvali(data_coarse, data_fine, lat_c, lon_c, lat_f, lon_f, varnum, domain, \
             strcmap, minmax, save, title, fname, lat_v, lon_v, data_v, lines):

    # Scrip to plot a field as a map

    varstr, mat_var, Variable, var_unit = var_names(varnum)

    if domain=='no':
            llon=lon_c[0]
            llat=lat_c[-1]
            ulon=lon_c[-1]
            ulat=lat_c[0]
            res = 5.0  # resolution for meridians & parallels
            bmres = 'l'
    else:
            llon=lon_f[0]
            llat=lat_f[-1]
            ulon=lon_f[-1]
            ulat=lat_f[0]
            res = 1.0  # resolution for meridians & parallels
            bmres = 'i'
    if domain=='no':
        llon = -4.0417
        ulon = 30.375
        llat = 48.575
        ulat = 65.875
        res = 5.0
        bmres = 'l'
    else:
        llon = 6.1736
        ulon = 14.9097
        llat = 53.229
        ulat = 56.4458
        res = 1.0
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

    # Parameters/Variables for plots
    coarse_min = np.min(data_coarse); fine_min = np.min(data_fine)
    coarse_max = np.max(data_coarse); fine_max = np.max(data_fine)
    if domain=='no':
            vmin = min(coarse_min, fine_min)
            vmax = max(coarse_max, fine_max)
            freeclim=1
    else:
	    vmin = fine_min
            vmax = fine_max
            freeclim=1
    print 'Min/Max data values: ', vmin, ', ', vmax
    print 'Min/Max in situ values: ', min(data_v), ', ', max(data_v)
    if minmax[0]!=minmax[1]:
        vmin = minmax[0]
        vmax = minmax[1]
        freeclim=0
        print 'Limit color range to: ', minmax[0], 'to', minmax[1]

    cmap = plt.cm.get_cmap(strcmap)
#    cmap1 = plt.cm.get_cmap('terrain_r')
    clev = np.linspace(vmin, vmax, 100)

    # Data - Coarse
    if domain=='no':
        lons,lats = np.meshgrid(lon_c, lat_c)
        x, y = m(lons,lats)
        if vmin < 0.0:
            if freeclim==1:
                cs1 = m.contourf(x, y, data_coarse, cmap=cmap, vmin=vmin, vmax=vmax, \
                                 levels=clev, norm=MidpointNormalize(midpoint=0.0))
            else:
                cs1 = m.pcolor(x, y, data_coarse, cmap=cmap, vmin=vmin, vmax=vmax, \
                               norm=MidpointNormalize(midpoint=0.0))
        else:
            if freeclim==1:
                cs1 = m.contourf(x, y, data_coarse, cmap=cmap, vmin=vmin, vmax=vmax, levels=clev)
            else:
                cs1 = m.pcolor(x, y, data_coarse, cmap=cmap, vmin=vmin, vmax=vmax)
		    
    # Data - Fine
    lons,lats = np.meshgrid(lon_f, lat_f)
    x, y = m(lons,lats)
    if 1==1:   # Use this to exclude the fine mesh when plotting the coarse
        if vmin < 0.0:
            if freeclim==1:
                cs2 = m.contourf(x, y, data_fine, cmap=cmap, vmin=vmin, vmax=vmax, \
                                 levels=clev, norm=MidpointNormalize(midpoint=0.0))
            else:
                cs2 = m.pcolor(x, y, data_fine, cmap=cmap, vmin=vmin, vmax=vmax, \
                               norm=MidpointNormalize(midpoint=0.0))
        else:
            if freeclim==1:
                cs2 = m.contourf(x, y, data_fine, cmap=cmap, vmin=vmin, vmax=vmax, levels=clev)
            else:
                cs2 = m.pcolor(x, y, data_fine, cmap=cmap, vmin=vmin, vmax=vmax)


    lons,lats = np.meshgrid(lon_v[:], lat_v[:])
    x, y = m(lon_v,lat_v)

#    s = m.scatter(x, y, marker='s', c=data, edgecolor='none',cmap='gist_rainbow')
    if lines==1:
        s = m.scatter(x, y, c=data_v, cmap=cmap, vmin=vmin, vmax=vmax, marker="s", s=30, edgecolor=[0.0,0.0,0.0], linewidth=0.5)
    else:
        s = m.scatter(x, y, marker='s', c=data_v, edgecolor='none',cmap=cmap, vmin=vmin, vmax=vmax)
    plt.clim(vmin, vmax)

    # Save/Show Map
    plt.title(title,fontsize=16)
    tick_locator = ticker.MaxNLocator(nbins=5)

    
    if domain=='no':
            cbar = m.colorbar(cs1)
    else:
            cbar = m.colorbar(cs2)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label(var_unit,fontsize=18)
    cbar.ax.tick_params(labelsize=16) 

    if save==1:
	    plt.savefig(fname, dpi=fig.dpi)
    else:
	    plt.show()

# Scatter plots like Himansu's paper
def plot_scatter(aux, insitu, title, varnum, ptype, coupling, fname, save):
    # Plots scatter

    varstr, mat_var, Variable, var_unit = var_names(varnum)

    plt.figure()
    plt.title(title, fontsize=16)
    plt.scatter(aux, insitu, marker='o',edgecolor='none',color=[0,0,0])
    rmse_t = rmse(aux, insitu); rmse_temp = round(rmse_t, 2)
    bias_t = bias(aux, insitu); bias_temp = round(bias_t, 2)
    _, pval = stats.ttest_ind(aux, insitu); ttest_temp = round(pval, 3)
    _, ppval = stats.ttest_rel(aux, insitu); pttest_temp = round(ppval, 3)
    if ptype<3:
        print mat_var.upper(), ': rmse, bias, pval, ppval: ', rmse_t, bias_t, pval, ppval
    else:
        print mat_var.upper(), ': ppval, pval, rmse, bias: ', ppval, pval, rmse_t, bias_t
    corr_temp = np.corrcoef(aux,insitu)[1][0]; corr_temp = round(corr_temp, 2)
    scatter_plot_with_correlation_line(aux, insitu)
    
    m, b = np.polyfit(aux, insitu, 1)

    # Labels etc
    y_lim = plt.ylim()[1]; x_lim = plt.xlim()[0]+1
    xmin=plt.xlim()[0]; xmax=plt.xlim()[1]
    ymin=plt.ylim()[0]; ymax=plt.ylim()[1]
    t_loc1 = ymax-(ymax-ymin)*0.05
    t_loc2 = ymax-(ymax-ymin)*0.1
    t_loc3 = ymax-(ymax-ymin)*0.15
    t_loc4 = ymax-(ymax-ymin)*0.2
    t_loc5 = ymax-(ymax-ymin)*0.25
    t_loc6 = ymax-(ymax-ymin)*0.3
    xpos = xmin + (xmax-xmin)*0.02
    plt.text(xpos, t_loc1, r'RMSE = '+str(rmse_temp),    {'color': 'k', 'fontsize': 12})
    plt.text(xpos, t_loc2, r'bias = '+str(bias_temp),    {'color': 'k', 'fontsize': 12})
    plt.text(xpos, t_loc3, r'r = '+str(corr_temp),    {'color': 'k', 'fontsize': 12})
    plt.text(xpos, t_loc4, r'n = '+str(len(aux)), {'color': 'k', 'fontsize': 12})
    if ptype<3:
        if ttest_temp>0.05:
            plt.text(xpos, t_loc5, r'pval = '+str(ttest_temp), {'color': 'k', 'fontsize': 12})
        else:
            plt.text(xpos, t_loc5, r'pval = '+str(ttest_temp)+'*', {'color': 'k', 'fontsize': 12})
    else:
        # Show values for paied t-test when assessbias bias or RMSE
        if ttest_temp>0.05:
            plt.text(xpos, t_loc5, r'pval = '+str(pttest_temp), {'color': 'k', 'fontsize': 12})
        else:
            plt.text(xpos, t_loc5, r'pval = '+str(pttest_temp)+'*', {'color': 'k', 'fontsize': 12})
    plt.text(xpos, t_loc6, r'm,b = '+str(round(m,2))+', '+str(round(b,2)), {'color': 'k', 'fontsize': 12})
    if ptype<3:
        xstr = 'in situ'
    else:
        xstr = 'FREE'
    if ptype==1:
        ystr = 'FREE'
    elif ptype>1:
        if coupling=='weak':
            ystr = 'WEAK'
        else:
            ystr = 'STRONG'
    plt.xlabel(xstr+' ('+var_unit+')')
    plt.ylabel(ystr+' ('+var_unit+')')

    
    if save==1:
        #plt.savefig(fname, bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.savefig(fname)
    elif save==0:
        plt.show()

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def bias(predictions, targets):
    forecast_errors = [targets[i]-predictions[i] for i in range(len(targets))]
    b = sum(forecast_errors) * 1.0/len(targets)
    return b

def scatter_plot_with_correlation_line(x, y):
    plt.gca()
    m, b = np.polyfit(x, y, 1)
    print 'm, b: ', m, b
    X_plot = np.linspace(plt.xlim()[0],plt.xlim()[1],2)
    limx = plt.xlim()
    limy = plt.ylim()
#    plt.plot(X_plot, m*X_plot + b, ls='--', c='b')
    minx = min(limx)
    maxx = max(limx)
    miny = min(limy)
    maxy = max(limy)
    minxy = min(minx, miny)
    maxxy = max(maxx, maxy)
    xmin = minxy-(maxxy-minxy)*0.1
    xmax = maxxy+(maxxy-minxy)*0.1
    plt.plot([xmin,xmax],[m*xmin+b,m*xmax+b], ls="--", c="k")
    plt.plot([xmin,xmax], [xmin,xmax], ls="-", c="k")
    plt.xlim(min(limx[0],limy[0]), max(limx[1],limy[1]))
    plt.ylim(min(limx[0],limy[0]), max(limx[1],limy[1]))
#    plt.xlim(230,330)
#    plt.ylim(230,330)

def read_ssh_track(ncfile, offset, tracktime):
    
    ssh_track = NetCDFFile(ncfile)
    ssh_keys = ssh_track.variables.keys()

    ssh_lon = ssh_track.variables['longitude'][:]
    ssh_lat = ssh_track.variables['latitude'][:]
    ssh_vals = ssh_track.variables['sla_filtered'][:]
    mdt_vals = ssh_track.variables['mdt'][:]
    tide_vals = ssh_track.variables['ocean_tide'][:]
    time_vals = ssh_track.variables['time'][:]

    ssh_vals = ssh_vals + mdt_vals + tide_vals - offset

    # For plotting the time of the tracks
    if tracktime==1:
        time_null = np.floor(time_vals[-1])
        time_vals = (time_vals - time_null)
        time_vals = time_vals *24.0 
        ssh_vals = time_vals

    for i in (range(len(ssh_lon))):
        if (ssh_lon[i]>=180.0):
            ssh_lon[i] = ssh_lon[i] - 360.0

    return ssh_vals, ssh_lat, ssh_lon


def vertical_mean_4D(data_coarse_pre, z1, z2):

    REFfile=NetCDFFile('/home/hbkycsun/plot/cellthc.nc')
    deptht=REFfile.variables['deptht'][:]

    cell=list(range(len(deptht)))
    cell[0]=(deptht[0]+deptht[1])/2
    cell[-1]=(deptht[-1]-deptht[-2])/2
    for i in (range(len(deptht))[1:-1]):
        cell[i]=(deptht[i+1]-deptht[i-1])/2

    if z1<=deptht[0]:
        id1=0
    else:
        for i in (range(len(deptht))[1:]):
            if deptht[i]>=z1 and deptht[i-1]<=z1:
                id1=i
    if z2>deptht[-1]:
        print 'Wrong z2'
        exit()
    elif z2==deptht[-1]:
        id2=len(deptht)-1
    else:
        for i in (range(len(deptht))):
            if deptht[i]<=z2 and deptht[i+1]>=z1:
                id2=i
	print 'Z id1=',id1,'Z id2=',id2

    if id1==id2:
        data_coarse=data_coarse_pre[id1,:,:]
    else:
        data_coarse_pre1=data_coarse_pre[id1:id2+1,:,:]
        print data_coarse_pre1.shape
        data_coarse_pre2=np.array(data_coarse_pre1)
        data_coarse_pre3=ma.masked_where(data_coarse_pre2>1000000, data_coarse_pre2)
        print data_coarse_pre3[:,100,100]
        data_coarse_pre4=np.ma.MaskedArray(data_coarse_pre3, mask=np.isnan(data_coarse_pre3))
        print data_coarse_pre4[:,100,100]
        data_coarse=np.ma.average(data_coarse_pre4, axis=0, weights=cell[id1:id2+1])
        #print cell[id1:id2+1]
        #data_coarse=np.average(data_coarse_pre3, axis=0, weights=cell[id1:id2+1])
        #for i in (range(id2-id1+1)):
        #    print i,data_coarse_pre3[i,100,100],cell[id1+i]
        #    data_coarse_pre3[i,:,:]=cell[id1+i]*data_coarse_pre3[i,:,:]
        #    print i,data_coarse_pre3[i,100,100]
        #data_coarse_pre2=np.array(data_coarse_pre1)
        #fill_value = 0
        #mask = np.isclose(data_coarse_pre2, fill_value)
        #data_coarse_pre2[mask] = np.nan
        #data_coarse=data_coarse_pre2.nanmean(axis=0)
        #data_coarse_pre4=np.nansum(data_coarse_pre3,axis=0)
        #print data_coarse_pre4[100,100]
        #data_coarse=data_coarse_pre4/np.sum(cell[id1:id2+1])
        #print data_coarse[100,100]
        #data_coarse_pre3=ma.masked_equal(data_coarse_pre2, fill_value)
        #data_coarse_pre4=data_coarse_pre3.mean(axis=0)
        #fill_value = 3.6135494e+20
        #mask = np.isclose(data_coarse, fill_value)
        #data_coarse[mask] = np.nan
        #data_coarse=ma.masked_values(data_coarse_pre4, fill_value)
        #data_coarse=ma.masked_where(data_coarse_pre4>1000000, data_coarse_pre4)
	print data_coarse[100,100]
    return data_coarse

def vertical_integral_4D(data_coarse_pre, z1, z2):

    REFfile=NetCDFFile('/home/hbkycsun/plot/cellthc.nc')
    deptht=REFfile.variables['deptht'][:]

    cell=list(range(len(deptht)))
    cell[0]=(deptht[0]+deptht[1])/2
    cell[-1]=(deptht[-1]-deptht[-2])/2
    for i in (range(len(deptht))[1:-1]):
        cell[i]=(deptht[i+1]-deptht[i-1])/2

    if z1<=deptht[0]:
        id1=0
    else:
        for i in (range(len(deptht))[1:]):
            if deptht[i]>=z1 and deptht[i-1]<=z1:
                id1=i
    if z2>deptht[-1]:
        print 'Wrong z2'
        exit()
    elif z2==deptht[-1]:
        id2=len(deptht)-1
    else:
        for i in (range(len(deptht))):
            if deptht[i]<=z2 and deptht[i+1]>=z1:
                id2=i
	print 'Z id1=',id1,'Z id2=',id2

    if id1==id2:
        data_coarse=data_coarse_pre[id1,:,:]
    else:
        data_coarse_pre1=data_coarse_pre[id1:id2+1,:,:]
        print data_coarse_pre1.shape
        data_coarse_pre2=np.array(data_coarse_pre1)
        data_coarse_pre3=ma.masked_where(data_coarse_pre2>1000000, data_coarse_pre2)
        print data_coarse_pre3[:,100,100]
        data_coarse_pre4=np.ma.MaskedArray(data_coarse_pre3, mask=np.isnan(data_coarse_pre3))
        print data_coarse_pre4[:,100,100]
        for i in (range(id2-id1+1)):
            print i,data_coarse_pre4[i,100,100],cell[id1+i]
            data_coarse_pre4[i,:,:]=cell[id1+i]*data_coarse_pre4[i,:,:]
            print i,data_coarse_pre4[i,100,100]
        data_coarse=np.nansum(data_coarse_pre4,axis=0)
	print data_coarse[100,100]
    return data_coarse

def get_model_idx(lat_stat, lon_stat, dist):
    # get model indices in dist 10 to input lat&lon

    # irrelevant variables that need to be specified to read model and obtain model lat lon
    varstr, _, _, _ = var_names(2)
    grid_area = 'coarse'
    year = 2015
    month = '02'
    day = '01'
    time_stamp = '00'
    depth = '0'
    DA_switch = 0
    coupled = 'weak'
    z_mean = 0
    z_integral = 0
    z1 = 0
    z2 = 0

    _, mod_keys = read_model(varstr, grid_area, year, month, day, time_stamp, depth, DA_switch, coupled, z_mean, z_integral, z1, z2)

    lat = mod_keys['lat'][:]
    lon = mod_keys['lon'][:]

    lat1 = lat[:,0]
    lon1 = lon[0,:]
    for i in range(len(lat[:,0])):
       lat1[i] = np.max(lat[i,:])
    for j in range(len(lat[0,:])):
       lon1[j] = np.max(lon[:,j])
    lon1[-1] = lon1[-2] + lon1[1]-lon1[0]

    latdiff = lat1 - lat_stat
    londiff = lon1 - lon_stat
    model_idx = []
    model_idx.append(np.argmin(abs(latdiff)))
    model_idx.append(np.argmin(abs(londiff)))

    return model_idx

def size(array):
    # function that returns the dimensions of a list (or list of lists etc.)
    dims = []
    while hasattr(array, '__iter__'):
        dims.append(len(array))
        array = array[0]
    return tuple(dims)

def transpose(matrix):
    return [list(row) for row in zip(*matrix)]

def read_insitu_vali_ICES(var_num, path_data):
    # read specific variable + other info from ICES insitu data

    print 'Start reading ', path_data

    # open Excel file
    wb = xlrd.open_workbook(path_data)
    worksheet = wb.sheet_by_index(0)

    headers = worksheet.row_values(0)
    data_rows = []
    for row_index in range(1,worksheet.nrows):
      row = worksheet.row(row_index)
      data_rows.append(row)

    # read relevant information
    cruise = read_var_insitu(0, data_rows)
    station = read_var_insitu(1, data_rows)
    type_station = read_var_insitu(2, data_rows)
    year = read_var_insitu(3,data_rows)
    month = read_var_insitu(4, data_rows)
    day = read_var_insitu(5, data_rows)
    hour = read_var_insitu(6, data_rows)
    minute = read_var_insitu(7, data_rows)
    lon = read_var_insitu(8, data_rows) # [deg east]
    lat = read_var_insitu(9, data_rows) # [deg north]
    bot_depth = read_var_insitu(10, data_rows) # [m]
    if "ICES_2015_Bottle_LR_CTD" in path_data:
      secchi_depth = read_var_insitu(12, data_rows) # [m]
      depth = read_var_insitu(14, data_rows) # [m]
    elif "ICES_2015_HR_CTD" in path_data:
      depth = read_var_insitu(12, data_rows) #[m]

    # read specified obs
    # Variable number from the var_names routine.
                        # Quick ref: 1=z, 2=TEM, 3=SAL, 4=uvel, 5=vvel, 6=NH4, 7=NO3, 8=PO4,
                        #            9=SIL, 10=DIA, 11=FLA, 12=CYA, 13=MEZ, 14=MIZ, 15=DET,
                        #            16=DETs, 17=FE, 18=LDON, 19=DIC, 20=ALK, 21=OXY, 22=pCO2,
                        #            23=PH, 24=CHL, 25=TE, 26=PFT, 27=PP
    if "ICES_2015_Bottle_LR_CTD" in path_data:
      if var_num == 2: #TEM
        idx = 18 # [deg C]
        idx_qv = 19
      elif var_num == 3:
        idx = 20 # practical salinity [dmnless?]
        idx_qv = 21
      elif var_num == 6: # Ammonium Nitrogen (NH4-N) [umol/l]
        idx = 34
        idx_qv = 35
        print 'WARNING: This is Ammonium Nitrogen (NH4-N) and not Ammonium (NH4)'
      elif var_num == 9:
        idx = 28 # Silicate Silicon (SiO4-Si) [umol/l]
        idx_qv = 29
      elif var_num == 20: # alkalinity
        idx = 42
        idx_qv = 43
      elif var_num == 21: # oxygen
        idx = 22 # dissolved oxygen [ml/l]
        idx_qv = 23
      elif var_num == 23: # PH
        idx = 40 # Hydrogen Ion Concentration (pH) [pH]
        idx_qv = 41
      elif var_num == 24: # CHL
        idx = 44 # Chlorophyll a [ug/l]
        idx_qv = 45
      else:
        print 'ERROR: this variable is not contained in the validation data set'
    elif "ICES_2015_HR_CTD" in path_data:
      print 'var_num == ', var_num
      if var_num == 2: # TEM
        print 'inside correct: TEM'
        idx = 16
        idx_qv = 17
      elif var_num == 3: # SAL
        idx = 18
        idx_qv = 19
      elif var_num == 21: #OXY
        idx = 20
        idx = 21
      else:
        print 'test ERROR: this variable is not contained in the validation data set'

    data = read_var_insitu(idx, data_rows)
    data_qv = read_var_insitu(idx_qv, data_rows)

    return year, month, day, hour, minute, lon, lat, depth, data, data_qv

def read_insitu_vali_CMEMS(varnum, data_path):
  # Read insitu data for specific variable
  if varnum == 1: #z
    varstr = 'SLEV' # water surface height above a specific datum
  elif varnum == 2: #tem
    varstr = 'TEMP'
    #varstr = 'TEMP_ADJUSTED' # sea temperature adjusted
    #varstr = 'SSJT' # sea temperature from TSG
  elif varnum == 3: # sal
    varstr = 'PSAL' # practical salinity
    #varstr = 'PSAD_ADJUSTED' # practical salinity adjusted
  elif varnum == 6: #NH4
    print 'CAUTION: Insitu data is not NH4 but NH4-N (Ammonium)'
    varstr = 'AMON'# NH4-N - Ammonium !!!
    #varstr = 'AMOW' # NH4-N !!!
  elif varnum == 7: # NO3
    print 'CAUTION: Insitu data is not NO3 but NO3-N'
    varstr = 'NTRA' # NO3-N - Nitrate!!
    #varstr = 'NTAW' # NO3-N !!
    #varstr = 'NTAW_ADJUSTED' # NO3-N adjusted
  elif varnum == 8: #PO4
    print 'CAUTION: Insitu data is not PO4 but PO4-P'
    varstr = 'PHOS'
    #varstr = 'PHOW'
  elif varnum == 9: #SIL
    varstr = 'SLCA' # Silicate (SIO4-SI)
    #varstr = 'SLCW' # Silicate (SIO4-SI)
  elif varnum == 20: #ALK
    varstr = 'ALKY' # Total alkalinity
    #varstr = 'ALKW' # Total alkalinity
    #varstr = 'ALK3' # Total alkalinity
  elif varnum == 21: #OXY
    #varstr = 'DOXY' # dissolved oxygen
    varstr = 'DOX1' # dissolved oxygen
    #varstr = 'DOX2' # dissolved oxygen
    #varstr = 'DOX2_ADJUSTED' ' dissolved oxygen adjusted
    #varstr = 'OSAT' # oxygen saturation
  elif varnum == 23: # PH
    varstr = 'PHPH' # Ph
  elif varnum == 24: # CHL
    #varstr = 'CPHL' # Chlorophyll-a
    #varstr = 'CH1P' # Chlorophyll-a (less divinyl chlorophyll-a)
    #varstr = 'Chlorophyll-a adjusted' # Chlorophyll-a adjusted
    #varstr = 'FLU2' # Chlorophyll-a fluorescence
    varstr = 'CHLT' # Total chlorophyll

  time_final = []
  data_final = []
  data_qc_final = []
  lat_final = []
  lon_final = []
  depth_final = []
  depth_qc_final = []
  data_source_final = np.array([])

  print 'READING VARIABLE = ', varstr
  for dirpath, dirnames, filenames in os.walk(data_path):
    for filename in filenames:
      if filename.endswith('.nc'):
        file_path = os.path.join(dirpath, filename)
        folder_path = '/scratch/usr/hbxsovli/SEAMLESS/validation_data/INSITU_BAL_PHYBGCWAV_DISCRETE_MYNRT_013_032/'
        idx_split = file_path.index(folder_path)
        data_source = file_path[idx_split + len(folder_path):]
        #file_path = '/scratch/usr/hbxsovli/SEAMLESS/validation_data/INSITU_BAL_PHYBGCWAV_DISCRETE_MYNRT_013_032/MO/BO_TS_MO_BothnianSea.nc'
        nc = NetCDFFile(file_path, 'r')
        if varstr in nc.variables and 'DEPH' in nc.variables and 'SYKE.nc' not in file_path:
          data = np.array(nc.variables[varstr][:])
          data_qc = np.array(nc.variables[varstr+'_QC'][:])
          lat = np.array(nc.variables['LATITUDE'][:])
          lon = np.array(nc.variables['LONGITUDE'][:])
          time = np.array(nc.variables['TIME'][:])
          depth = np.array(nc.variables['DEPH'][:])
          ref_date = datetime.datetime(1950,1,1,0,0,0)
          deltas = [datetime.timedelta(days=f) for f in time]
          time_dt = ref_date + np.array(deltas)
          # select only relevent timespan
          notnan_mask = ~np.isnan(data) & (data>-999).astype(bool)
          year_mask = np.array([date.year == 2015 for date in time_dt])
          month_mask = np.array([date.month in [2,3,4,5] for date in time_dt])
          lat_mask = (lat>50).astype(bool) & (lat<66).astype(bool)
          lon_mask = (lon>-4).astype(bool) & (lon<31).astype(bool)
          hour_mask = np.array([date.hour in [22,23,0,24,1,2] for date in time_dt])
          if len(data.shape)>1:
            combined_mask = year_mask & month_mask & hour_mask & lat_mask & lon_mask
            if np.any(combined_mask):
              depth_temp = depth[combined_mask,:]
              data_temp = data[combined_mask,:]
              data_qc_temp = data_qc[combined_mask,:]
              lat_temp = lat[combined_mask]
              lon_temp = lon[combined_mask]
              time_temp = time_dt[combined_mask]
              if len(data_temp)==len(lat_temp):
                data_len = len(data_temp[0])
              elif len(data_temp[0])==len(lat_temp):
                data_len = len(data_temp)
              lat_matrix = np.transpose(np.tile(lat_temp,(data_len,1)))
              lon_matrix = np.transpose(np.tile(lon_temp,(data_len,1)))
              time_matrix = np.transpose(np.tile(time_temp,(data_len,1)))
              notnan_mask = ~np.isnan(data_temp) & ~np.isnan(depth_temp) & (data_temp>-999).astype(bool)
              if np.any(notnan_mask):
                print file_path
                print len(time_matrix[notnan_mask])
                time_final = np.append(time_final, time_matrix[notnan_mask])
                lat_final = np.append(lat_final, lat_matrix[notnan_mask])
                lon_final = np.append(lon_final, lon_matrix[notnan_mask])
                if varstr == 'DOX1':
                  data_temp = 1/0.022391*data_temp[notnan_mask]#*0.001
                  print data_temp
                else: data_temp = data_temp[notnan_mask]
                data_final = np.append(data_final, data_temp)
                data_qc_final = np.append(data_qc_final, data_qc_temp[notnan_mask])
                depth_final = np.append(depth_final, depth_temp[notnan_mask])
                data_source_final = np.append(data_source_final, np.transpose([str(data_source)]*len(depth_temp[notnan_mask])))
          else:
            notnan_mask = ~np.isnan(data)
            combined_mask= year_mask & month_mask & hour_mask & notnan_mask & lat_mask & lon_mask
            if np.any(combined_mask):
              print file_path
              print len(time_dt[combined_mask])
              time_final = np.append(time_final, time_dt[combined_mask])
              if varstr == 'DOX1':
                data_temp = 1/0.022391*data[combined_mask]#*0.001
                print data_temp
              else: data_temp = data[combined_mask]
              data_final = np.append(data_final, data_temp)
              data_qc_final = np.append(data_qc_final, data_qc[combined_mask])
              lat_final = np.append(lat_final, lat[combined_mask])
              lon_final = np.append(lon_final, lon[combined_mask])
              depth_final = np.append(depth_final, depth[combined_mask])
              data_source_final.append([str(data_source)]*len(depth[combined_mask]))
              if len(data_final) != len(lat_final):
                print 'PROBLEM LENGTH NOT THE SAME'
                break
              else:
                print len(lat_final), len(data_final)
          del data, data_qc, lat, lon, time, depth, deltas, time_dt, year_mask, month_mask, hour_mask, combined_mask, data_source
        #else:
        #  print 'Variable not found in', file_path

  return time_final, lon_final, lat_final, depth_final, data_final, data_qc_final, data_source_final

def unique_list(lst, tol):
    unique_lst = []
    for x in lst:
        for uniq in unique_lst:
            if abs(x-uniq) < tol:
                break
        else:
            unique_lst.append(x)
    return unique_lst

def unique_tuples(a, tol=0):
    unique_list = []
    unique_indices = []
    count = 0
    for i, x in enumerate(a):
    #    print count,'/',len(a)
        count = count + 1
        found = False
        for j, y in enumerate(unique_list):
            if len(x) != len(y):
                continue
            match = all(abs(xi - yi) <= tol for xi, yi in zip(x, y))
            if match:
                found = True
                break
        if not found:
            unique_list.append(x)
            unique_indices.append(i)
    return unique_list, unique_indices


def unique_tuples2(a, tol=0):
    unique_dict = {}
    unique_indices = []
    for i, x in enumerate(a):
        key = tuple(x)
        if key not in unique_dict:
            unique_dict[key] = i
            unique_indices.append(i)
        else:
            existing_index = unique_dict[key]
            existing_tuple = a[existing_index]
            if len(existing_tuple) == len(x) and all(abs(xi - yi) <= tol for xi, yi in zip(x, existing_tuple)):
                continue
            unique_dict[key] = i
            unique_indices.append(i)
    unique_list = [a[i] for i in unique_indices]
    return unique_list, unique_indices

