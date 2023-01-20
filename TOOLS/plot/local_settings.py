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
        #path = data_loc+'/ens_NORDIC_ndays/output_N20_DA_Smago/DA/'
        #path = data_loc+'/ens_ERGOM_ndays/output_N4_bio/DA/'
        #path = data_loc+'/ens_ERGOM_ndays/001/' #output/DA/'
        #path = data_loc+'/ens_ERGOM_ndays/output_N30_log_chl30k/DA/'
        path = data_loc+'/LESTKF/DA/'
        #path = data_loc+'/Hybrid/DA/'
    else:
#        data_loc = '/scratch/usr/hbkycsun/data'
        # path to free running experiment
#        path = data_loc+'/FreeRun/Post_DA/'
        path = data_loc+'/free_N30/DA/'
#        path = data_loc+'/Spinup/DA/'

    return path