import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as bp


file_free = 'rms_out_N20_free';
file = 'rms_out_N20_f09_SST';
file = 'rms_out_N20_f09_SSH+SST_varscale04';

operational=0 # 1 for operational, 0 otherwise

save = 1

ncfile_free=nc.Dataset(file_free+'.nc');
ncfile=nc.Dataset(file+'.nc');

rmse = ncfile['rmse']
rmse_free = ncfile_free['rmse']

fig=plt.figure(figsize=(12,6))
if operational==1:
   plt.plot(rmse_free[0,:], c='k', label='Free')
   plt.plot(rmse[2,:], c='r', label='Forecast')
   plt.plot(rmse[1,:], c='b', label='Analysis')
else:
   plt.plot(rmse_free[1,:], c='k', label='Free')
   plt.plot(rmse[1,:], c='b', label='Forecast')
   plt.plot(rmse[2,:], c='r', label='Analysis')

plt.title(file, fontsize=20)
plt.ylabel('RMSe ($^oC$)', fontsize=16)
plt.xlabel('day', fontsize=16)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.legend(loc='lower left',fontsize=16)

if save==1:
  fname=file+'.png'
  plt.savefig(fname)
else:
  plt.show()

