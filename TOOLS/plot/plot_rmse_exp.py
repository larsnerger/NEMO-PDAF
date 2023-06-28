import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as bp
from plot_lib import *

# CONFIGURATION PART
fname = 'rmse_sst_ba.png'
plot_title = 'RMSe SST'
nc_folder = '../rmse/netcdfs/' # if not in a different path just just empty string ('')
file_free = nc_folder+'rms_sst_L3S_ba_free_N30.nc';
file1 = nc_folder+'rms_sst_L3S_ba_sstL3_N30.nc'; #chl_N30
plot_ana = 1 # whether to plot analysis step
operational = 1 # 1 for operational, 0 otherwise

# END CONFIGURATION

# read rmse netcdfs
ncfile_free=nc.Dataset(file_free);
ncfile1=nc.Dataset(file1);

rmse1 = ncfile1['rmse']
rmse_free = ncfile_free['rmse']

# PLOTTING
# Determine tick localizations and labels 
allmonths = [1, 2, 3, 4, 5]
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
   plt.plot(allday_assim, rmse1[0,:], 'g', label='Forecast SST_DA')
   if plot_ana == 1:
     plt.plot(allday_assim, rmse1[1,:], 'g',linestyle='--', label='Analysis SST_DA')
else:
   plt.plot(allday_free, rmse_free[1,:], 'b', label='FREE')
   plt.plot(allday_assim, rmse1[1,:], '#bbbb00', label='Forecast SST_DA')
   if plot_ana == 1:
     plt.plot(allday_assim, rmse1[2,:], '#bbbb00', linestyle='--', label='Analysis SST_DA')

plt.title(plot_title)
plt.ylabel('RMSe ($^oC$)')
#plt.ylabel('RMSe (mg Chl/m$^3$)', fontsize=16)
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
print 'sst_N30 = ', np.mean(rmse1[0,:])
if plot_ana == 1: 
  print '_________ Analysis __________'
  print 'sst_N30 = ', np.mean(rmse1[1,:])

#plt.show()
plt.savefig(fname, dpi=200)

