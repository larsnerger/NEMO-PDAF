import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
from plot_lib import *

# Select variable to plot
var = 'chl_ba'
var = 'chl_no'
var = 'sst_L4'


experiment = 'sst-chl_Tonly_N30'

path_free = '/scratch/projects/hbk00095/exp/free_N30/rmse'
path_da = '/scratch/projects/hbk00095/exp/'+experiment+'/rmse'

file_free = 'rms_'+var+'_free_N30.nc';
file_da = 'rms_'+var+'_'+experiment+'.nc';

print 'Free: ', file_free 
print 'DA:   ', file_da

save = 0


##########################

ncfile_free=nc.Dataset(path_free+'/'+file_free);
ncfile_da=nc.Dataset(path_da+'/'+file_da);

rmse_da = ncfile_da['rmse']
doy_da = ncfile_da['doy']
rmse_free = ncfile_free['rmse']
doy_free = ncfile_free['doy']
days_free = ncfile_free['days']


# Determine tick locations and labels

tickloc = []
tickloc_min = []
loc = doy_free[0]
tickloc.append(loc)
months_str = []
for i in range(len(days_free)):
   _, _, mstr, days_in_mon = month_names(i+1)
   loc = loc + days_in_mon
   tickloc.append(loc)
   tickloc_min.append(tickloc[i]+days_in_mon/2)
   months_str.append(mstr)

# Plotting

fig, ax=plt.subplots(figsize=(8,5))
plt.plot(doy_free[:], rmse_free[0,:], 'k', label='Free')
plt.plot(doy_da[:], rmse_da[0,:], 'b', label='Forecast')
plt.plot(doy_da[:], rmse_da[1,:], 'r', label='Analysis')

plt.title(file_da, fontsize=16)
if var=='sst' or var=='sst_L3S':
   plt.ylabel('RMSe ($^oC$)', fontsize=16)
else:
   plt.ylabel('RMSe (mg Chl/m$^3$)', fontsize=16)
plt.xlim(tickloc[0],tickloc[-1])
ax.xaxis.set_major_locator(ticker.FixedLocator(tickloc))
ax.xaxis.set_minor_locator(ticker.FixedLocator(tickloc_min))
ax.xaxis.set_major_formatter(ticker.NullFormatter())
ax.tick_params(axis='x', which='minor', tick1On=False, tick2On=False)
ax.xaxis.set_minor_formatter(ticker.FixedFormatter(months_str))        
plt.legend(loc='upper right',fontsize=16)

if save==1:
  fname=file+'.png'
  print 'Save to', fname
  plt.savefig(fname)
else:
  plt.show()

