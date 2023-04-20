#!/usr/bin/python
# Functions for local settings:
# base_dir()
#     Set path to root directory of HBM experiments
# get_exp_path()
#     Set full path to experiment

# Define path to base directory for HBM experiments
def base_dir():

#    basedir = '/scratch/usr/hzfblner/SEAMLESS/run'
    basedir = '/scratch/projects/hbk00095/exp'

    return basedir


# define full path for output directory for current month
def get_exp_path(coupled, year, month, day, ampm, asml):

    data_loc = base_dir()

    YYYYMM = str(year)+str(month)+str(day)+str(ampm)

    if asml==1:
        # Path to assimilation experiment
#        path = data_loc+'/free_N30/DA/'
        path = '/scratch/usr/hbknerge/SEAMLESS/run/DA-SST-CHL/sst-chl_Tonly_N30/DA/'
        #path = '/scratch/usr/hbknerge/SEAMLESS/run/DA-SST-CHL/free_initMarch1_N30/DA/'
        #path = '/scratch/usr/hbxsovli/SEAMLESS/runfolders/DA_chl/exp.chl_30ens/DA/'
    else:
        # path to free running experiment
        path = data_loc+'/free_N30/DA/'
        #path = '/scratch/usr/hbknerge/SEAMLESS/run/DA-SST-CHL/sst-chl_Tonly_N30/DA/'


    return path
