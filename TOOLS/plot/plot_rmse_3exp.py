import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as bp
from plot_lib import *

# CONFIGURATION PART
fname = 'rmse_log_chl_ba.png'
plot_title = 'RMSe CHL (log)'
nc_folder = '../rmse/netcdfs/' # if not in a different path just just empty string ('')
file_free = nc_folder+'rms_chl_ba_free_N30_log.nc';
file1 = nc_folder+'rms_chl_ba_chl_N30_log.nc'; #chl_N30
file2 = nc_folder+'rms_chl_ba_sstL3_N30_log.nc'; #sstL3
file3 = nc_folder+'rms_chl_ba_chl+sstL3_N30_log.nc' #chl+sstL3
plot_ana = 0 # whether to plot analysis step
operational = 1 # 1 for operational, 0 otherwise

# END CONFIGURATION

# read rmse netcdfs
ncfile_free=nc.Dataset(file_free);
ncfile1=nc.Dataset(file1);
ncfile2=nc.Dataset(file2);
ncfile3=nc.Dataset(file3);

rmse1 = ncfile1['rmse']
rmse2 = ncfile2['rmse']
rmse3 = ncfile3['rmse']
rmse_free = ncfile_free['rmse']

# PLOTTING
# Determine tick localizations and labels 
allmonths = [2, 3, 4, 5]
allday_free = ncfile_free['doy']
allday_assim = ncfile1['doy']
tickloc = []
tickloc_min = []
loc = allday_free[0]
tickloc.append(loc)
months_str = []
for i in range(len(allmonths)):
    _, _, mstr, days_in_mon = month_names(i+allmonths[0])
    loc = loc + days_in_mon
    tickloc.append(loc)
    tickloc_min.append(tickloc[i]+days_in_mon/2)
    months_str.append(mstr)


fig, ax = plt.subplots(figsize=(8,3))
if operational==1:
   plt.plot(allday_free, rmse_free[0,:], 'b', label='FREE')
   plt.plot(allday_assim, rmse1[0,:], '#bbbb00', label='CHL_DA')
   plt.plot(allday_assim, rmse2[0,:], 'g', label='SST_DA')
   plt.plot(allday_assim, rmse3[0,:], 'k', label='CHL+SST_DA')
   if plot_ana == 1:
     plt.plot(allday_assim, rmse1[1,:], '#bbbb00',linestyle='--', label='CHL_DA')
     plt.plot(allday_assim, rmse2[1,:], 'g--', label='SST_DA')
     plt.plot(allday_assim, rmse3[1,:], 'k--', label='CHL+SST_DA')
else:
   plt.plot(allday_free, rmse_free[1,:], 'b', label='FREE')
   plt.plot(allday_assim, rmse1[1,:], '#bbbb00', label='CHL_DA')
   plt.plot(allday_assim, rmse2[1,:], 'g', label='SST_DA')
   plt.plot(allday_assim, rmse3[1,:], 'k', label='CHL+SST_DA')
   if plot_ana == 1:
     plt.plot(allday_assim, rmse1[2,:], '#bbbb00', linestyle='--', label='CHL_DA')
     plt.plot(allday_assim, rmse2[2,:], 'g--', label='SST_DA')
     plt.plot(allday_assim, rmse3[2,:], 'k--', label='CHL+SST_DA')

plt.title(plot_title)
#plt.ylabel('RMSe ($^oC$)', fontsize=16)
plt.ylabel('log RMSe')
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
print 'Free = ', np.mean(rmse_free)
print '_________ Forecast __________'
print 'chl_N30 = ', np.mean(rmse1[0,:])
print 'sstL3_N30 = ', np.mean(rmse2[0,:])
print 'chl+sstL3_N30 = ', np.mean(rmse3[0,:])
if plot_ana == 1: 
  print '_________ Analysis __________'
  print 'chl_N30 = ', np.mean(rmse1[1,:])
  print 'sstL3_N30 = ', np.mean(rmse2[1,:])
  print 'chl+sstL3_N30 = ', np.mean(rmse3[1,:])

#plt.show()
plt.savefig(fname, dpi=200)
