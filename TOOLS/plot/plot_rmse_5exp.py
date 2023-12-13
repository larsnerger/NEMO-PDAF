import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as bp
from plot_lib import *

# CONFIGURATION PART
fname = 'rmse_one-year_sst_ba.png'
nc_folder = '../rmse/' # if not in a different path just just empty string ('')
varnum = 2 #2-TEM, 24-CHL
allmonths = range(1,13)
exp_free = 'free'
exp1 = 'chl_weakly_sstL3_strongly'
#exp2 = 'sstL3_strongly'
#exp3 = 'chl_weakly_sstL3_strongly'
#exp4 = 'chl_strongly_sstL3_strongly'
#exp5 = 'chl_weakly_sstL3_weakly' #with dashed line
plot_ana = 1 # whether to plot analysis step
operational = 1 # 1 for operational, 0 otherwise
save = 1 #whether to save plot 
plot_leg = 1 #whether to plot legend in plot 
# END CONFIGURATION

if varnum == 24:
  logrms = '_log'
  rmsfile = 'rms_chl_ba_'
  plot_title = 'RMSe CHL log'
elif varnum == 2:
  logrms = ''
  rmsfile = 'rms_sst_L3S_ba_'
  plot_title = 'RMSe SST'

file_free = nc_folder + rmsfile + exp_free + '_N30' + logrms + '.nc'
if 'exp1' in locals():
  file1 = nc_folder + rmsfile + exp1 + '_N30' + logrms + '.nc'
  print file1
  ncfile1=nc.Dataset(file1)
  rmse1 = ncfile1['rmse']
 
if 'exp2' in locals():
  file2 = nc_folder + rmsfile + exp2 + '_N30' + logrms + '.nc'
  print file2
  ncfile2=nc.Dataset(file2)
  rmse2 = ncfile2['rmse']

if 'exp3' in locals():
  file3 = nc_folder + rmsfile + exp3 + '_N30' + logrms + '.nc'
  print file3
  ncfile3=nc.Dataset(file3)
  rmse3 = ncfile3['rmse']

if 'exp4' in locals():
  file4 = nc_folder + rmsfile + exp4 + '_N30' + logrms + '.nc'
  print file4
  ncfile4=nc.Dataset(file4);
  rmse4 = ncfile4['rmse']

if 'exp5' in locals():
  file5 = nc_folder + rmsfile + exp5 + '_N30' + logrms + '.nc'
  print file5
  ncfile5=nc.Dataset(file5);
  rmse5 = ncfile5['rmse']

print file_free
print(plot_title)
# read rmse netcdfs
ncfile_free=nc.Dataset(file_free);
rmse_free = ncfile_free['rmse']

# PLOTTING
# Determine tick localizations and labels 
allday_free = ncfile_free['doy']
allday_assim = ncfile1['doy']
tickloc = []
tickloc_min = []
if varnum == 24:
  loc = allday_free[0]
elif varnum == 2:
  #loc = allday_free[30]
  loc = allday_free[0]
tickloc.append(loc)
months_str = []
for i in range(len(allmonths)):
    _, mstr, _, days_in_mon = month_names(i+allmonths[0])
    loc = loc + days_in_mon
    tickloc.append(loc)
    tickloc_min.append(tickloc[i]+days_in_mon/2)
    months_str.append(mstr)

if varnum == 2:
  print len(rmse_free[0,:])
  #print len(rmse_free[0,30:-1])
elif varnum == 24:
  print len(rmse_free[0,:])
print len(rmse1[0,:])

#fig, ax = plt.subplots(figsize=(8,3))
fig, ax = plt.subplots(figsize=(8,3))
if operational==1:
   if varnum == 24:
     plt.plot(allday_free, rmse_free[0,:], 'b', label='FREE')
   elif varnum == 2:
     plt.plot(allday_free, rmse_free[0,:], 'b', label='FREE')
   if plot_ana == 1:
     #plt.plot(allday_assim, rmse1[1,:], '#bbbb00',linestyle='--', label='CHL strongly')
     plt.plot(allday_assim, rmse1[1,:], 'grey', label='CHL w + SST s Ana')
     #plt.plot(allday_assim, rmse2[1,:], 'g--', label='SST strongly')
     #plt.plot(allday_assim, rmse3[1,:], 'k--', label='CHL weakly + SST strongly')
     #plt.plot(allday_assim, rmse_free[0,30:-1], 'b', label='FREE')
   if 'exp1' in locals():
     plt.plot(allday_assim, rmse1[0,:], color='k', label='CHL w + SST s Fcst')
   if 'exp2' in locals():
     plt.plot(allday_assim, rmse2[0,:], 'g', label='SST s')
   if 'exp3' in locals():
     plt.plot(allday_assim, rmse3[0,:], 'k', label='CHL w + SST s')
   if 'exp4' in locals():
     plt.plot(allday_assim, rmse4[0,:], 'darkviolet', label='CHL s + SST s')
   if 'exp5' in locals(): 
     plt.plot(allday_assim, rmse5[0,:], 'saddlebrown', linestyle = '--', label='CHL w + SST w')
else:
   plt.plot(allday_free, rmse_free[1,:], 'b', label='FREE')
   plt.plot(allday_assim, rmse1[1,:], '#bbbb00', label='CHL strongly')
   plt.plot(allday_assim, rmse2[1,:], 'g', label='SST strongly')
   plt.plot(allday_assim, rmse3[1,:], 'k', label='CHL weakly + SST strongly')
   if plot_ana == 1:
     plt.plot(allday_assim, rmse1[2,:], '#bbbb00', linestyle='--', label='CHL strongly')
     plt.plot(allday_assim, rmse2[2,:], 'g--', label='SST strongly')
     plt.plot(allday_assim, rmse3[2,:], 'k--', label='CHL weakly + SST strongly')

plt.title(plot_title)
if varnum == 2:
  plt.ylabel('RMSe ($^oC$)', fontsize=14)
if varnum == 24: 
  plt.ylabel('log RMSe', fontsize=10)
if plot_leg == 1:
  plt.legend(loc='upper right')
plt.xlim(tickloc[0],tickloc[-1])
plt.ylim(ymin=0.0)
ax.xaxis.set_major_locator(ticker.FixedLocator(tickloc))
ax.xaxis.set_minor_locator(ticker.FixedLocator(tickloc_min))
ax.xaxis.set_major_formatter(ticker.NullFormatter())
ax.tick_params(axis='x', which='minor', tick1On=False, tick2On=False)
ax.xaxis.set_minor_formatter(ticker.FixedFormatter(months_str))

# print RMSE mean to screen 
print '----- MEAN RMSE ------'
print 'Free = ', np.mean(rmse_free[0,:])
#print 'Free = ', np.mean(rmse_free[0,30:-1])
print '_________ Forecast __________'
if 'exp1' in locals():
  print exp1,' = ', np.mean(rmse1[0,:])
if 'exp2' in locals():
  print exp2, ' = ', np.mean(rmse2[0,:])
if 'exp3' in locals():
  print exp3, ' = ', np.mean(rmse3[0,:])
if 'exp4' in locals():
  print exp4, ' = ', np.mean(rmse4[0,:])
if 'exp5' in locals():
  print exp5, ' = ', np.mean(rmse5[0,:])
if plot_ana == 1: 
  print '_________ Analysis __________'
  print 'chl_N30 = ', np.mean(rmse1[1,:])
  #print 'sstL3_N30 = ', np.mean(rmse2[1,:])
  #print 'chl+sstL3_N30 = ', np.mean(rmse3[1,:])

#plt.show()
if save == 1: 
  plt.savefig(fname, dpi=200)

