import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as bp


file_free = 'rms_out_N20_free.nc';
file1 = 'rms_out_N20_f09.nc';
file2 = 'rms_out_N20_f09_SSH+SST_varscale04.nc';

operational=0 # 1 for operational, 0 otherwise


ncfile_free=nc.Dataset(file_free);
ncfile1=nc.Dataset(file1);
ncfile2=nc.Dataset(file2);

rmse1 = ncfile1['rmse']
rmse2 = ncfile2['rmse']
rmse_free = ncfile_free['rmse']

fig=plt.figure(figsize=(12,6))
if operational==1:
   plt.plot(rmse_free[0,:], 'k', label='Free')
   plt.plot(rmse1[0,:], 'r', label='Forecast 1')
   plt.plot(rmse1[1,:], 'b', label='Analysis 1')
   plt.plot(rmse2[0,:], 'r--', label='Forecast 2')
   plt.plot(rmse2[1,:], 'b--', label='Analysis 2')
else:
   plt.plot(rmse_free[1,:], 'k', label='Free')
   plt.plot(rmse1[1,:], 'b', label='Forecast 1')
   plt.plot(rmse1[2,:], 'r', label='Analysis 1')
   plt.plot(rmse2[1,:], 'b--', label='Forecast 2')
   plt.plot(rmse2[2,:], 'r--', label='Analysis 2')

plt.title('1: '+file1+' 2:'+file2, fontsize=16)
plt.ylabel('RMSe ($^oC$)', fontsize=16)
plt.legend(loc='lower left')

plt.show()

