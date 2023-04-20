import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
from plot_lib import *

# Select variable to plot
var = 'chl_ba'
var = 'chl_no'
var1 = 'sst_L3S'
var2 = 'sst_L4'
var = 'sst'

experiment1 = 'sst-chl_Tonly_N30'
experiment2 = 'sst-chl_Tonly_N30'

fname = 'rmse_'+var1+'_'+var2+'_'+experiment1+'_'+experiment1

path_free = '/scratch/projects/hbk00095/exp/free_N30/rmse'
path_da1 = '/scratch/projects/hbk00095/exp/'+experiment1+'/rmse'
path_da2 = '/scratch/projects/hbk00095/exp/'+experiment1+'/rmse'

file_free = 'rms_'+var2+'_free_N30.nc';
file_da1 = 'rms_'+var1+'_'+experiment1+'.nc';
file_da2 = 'rms_'+var2+'_'+experiment2+'.nc';

print 'Free: ', file_free 
print 'DA1:  ', file_da1
print 'DA2:  ', file_da2


save = 0


##########################

ncfile_free = nc.Dataset(path_free+'/'+file_free);
ncfile_da1 = nc.Dataset(path_da1+'/'+file_da1);
ncfile_da2 = nc.Dataset(path_da2+'/'+file_da2);

rmse_da1 = ncfile_da1['rmse']
doy_da1 = ncfile_da1['doy']
rmse_da2 = ncfile_da2['rmse']
doy_da2 = ncfile_da2['doy']
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

fig, ax=plt.subplots(figsize=(10,6))
plt.plot(doy_free[:], rmse_free[0,:], 'k', label='Free')
plt.plot(doy_da1[:], rmse_da1[0,:], 'b', label='Fcst '+var1)
plt.plot(doy_da1[:], rmse_da1[1,:], 'r', label='Ana '+var1)
plt.plot(doy_da2[:], rmse_da2[0,:], 'b--', label='Fcst '+var2)
plt.plot(doy_da2[:], rmse_da2[1,:], 'r--', label='Ana '+var2)

plt.title('1: '+experiment1+' 2:'+experiment2, fontsize=16)
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

if save>0:
  if save==1:
     fname=fname+'.png'
     print 'Save to', fname
     plt.savefig(fname)
  else:
     fname=fname+'_prt.png'
     print 'Save to', fname
     plt.savefig(fname, dpi=300)
else:
  plt.show()

