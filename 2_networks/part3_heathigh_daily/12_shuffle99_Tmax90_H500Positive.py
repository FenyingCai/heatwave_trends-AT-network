import numpy as np
import xarray as xr

##--##--##--##--  读取 15 天滑动气温 10%, 90% 阈值 (JJA)   --##--##--##--##
ds00 = xr.open_dataset("/public/home/fcai/extreme1_AT/NC2_1980_2010/1_threshold/era5_Tmax_th90_JJA_1981_2010.nc")
t_th90 = ds00.air_th90.loc[:,85:10,:][:,::-1,:]
ds01 = xr.open_dataset("/public/home/fcai/extreme1_AT/NC2_1980_2010/1_threshold/era5_H500_std_JJA_1981_2010.nc")
h_std = ds01.h_std.loc[:,85:10,:][:,::-1,:]
ds02 = xr.open_dataset("/public/home/fcai/extreme1_AT/NC2_1980_2010/1_threshold/era5_H500_clm_JJA_1981_2010.nc")
h_clm = ds02.h_clm.loc[:,85:10,:][:,::-1,:]


##--##--##--##--  读取 ERA5 2米气温数据  --##--##--##--##
ds00 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/t2m_max/t2m.2021-05.daily.nc")
t0 = ds00.t2m.loc[:,85:10,:][:,::-1,:]; ds00.close()
lat = t0.latitude; lon = t0.longitude
print(np.shape(t0))
t = np.zeros((44, 92, np.shape(t0)[1],np.shape(t0)[2]))




##--##--##--##--  Tmax, 90th  --##--##--##--##
for iyear in range(44):
  print('iyear = ',iyear+1979)
  ds2 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/t2m_max/t2m."+str(iyear+1979)+"-06.daily.nc")
  t[iyear,0:30,:,:] = ds2.t2m.loc[:,85:10,:][:,::-1,:]
  ds3 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/t2m_max/t2m."+str(iyear+1979)+"-07.daily.nc")
  t[iyear,30:61,:,:] = ds3.t2m.loc[:,85:10,:][:,::-1,:]
  ds4 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/t2m_max/t2m."+str(iyear+1979)+"-08.daily.nc")
  t[iyear,61:92,:,:] = ds4.t2m.loc[:,85:10,:][:,::-1,:]
  ds2.close(); ds3.close(); ds4.close()

t_True = np.zeros((44, 92, np.shape(t0)[1],np.shape(t0)[2]))
for iyear in range(44):
  t_True[iyear,:,:,:] = np.array(t[iyear,:,:,:]) >= np.array(t_th90)
del t_th90, t  #; print(np.sum(t_True)); print(np.size(t_True))

t_True2 = t_True.copy()
t_True2[:,2:-2,:,:][t_True[:,:-4,:,:]+t_True[:,1:-3,:,:]+t_True[:,3:-1,:,:]+t_True[:,4:,:,:]<=1] = 0.0
t_True2[:,2:-2,:,:][t_True[:,1:-3,:,:]+t_True[:,3:-1,:,:]==0] = 0.0
t_True2[:,2:-2,:,:][(t_True[:,1:-3,:,:]+t_True[:,3:-1,:,:]==1)&(t_True[:,3:-1,:,:]+t_True[:,4:,:,:]<=1)&(t_True[:,:-4,:,:]+t_True[:,1:-3,:,:]<=1)] = 0.0

t_True2[:,0,:,:][t_True[:,1,:,:]+t_True[:,2,:,:]<=1] = 0.0
t_True2[:,1,:,:][t_True[:,0,:,:]+t_True[:,2,:,:]+t_True[:,3,:,:]<=1] = 0.0
t_True2[:,1,:,:][(t_True[:,0,:,:]==1)&(t_True[:,2,:,:]==0)&(t_True[:,3,:,:]==1)] = 0.0
t_True2[:,-1,:,:][t_True[:,-2,:,:]+t_True[:,-3,:,:]<=1] = 0.0
t_True2[:,-2,:,:][t_True[:,-1,:,:]+t_True[:,-3,:,:]+t_True[:,-4,:,:]<=1] = 0.0
t_True2[:,-2,:,:][(t_True[:,-1,:,:]==1)&(t_True[:,-3,:,:]==0)&(t_True[:,-4,:,:]==1)] = 0.0
del t_True





##--##--##--##--  H500, >1std  --##--##--##--##
h = np.zeros((44, 92, np.shape(t0)[1],np.shape(t0)[2]))
for iyear in range(44):
  print('iyear = ',iyear+1979)
  ds2 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/zg/zg."+str(iyear+1979)+"-06.daily.nc")
  h[iyear,0:30,:,:] = ds2.z.loc[:,500,85:10,:][:,::-1,:]
  ds3 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/zg/zg."+str(iyear+1979)+"-07.daily.nc")
  h[iyear,30:61,:,:] = ds3.z.loc[:,500,85:10,:][:,::-1,:]
  ds4 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/zg/zg."+str(iyear+1979)+"-08.daily.nc")
  h[iyear,61:92,:,:] = ds4.z.loc[:,500,85:10,:][:,::-1,:]
  ds2.close(); ds3.close(); ds4.close()

h_Positive = np.zeros((44, 92, np.shape(t0)[1],np.shape(t0)[2]))
for iyear in range(44):
  h_ano = np.array(h[iyear,:,:,:]) - np.array(h_clm)
  h_Positive[iyear,:,:,:] = np.array(h_ano) >= np.array(h_std); del h_ano
del h, h_clm  #; print(np.sum(t_True)); print(np.size(t_True))








##--##--##--##--  把不通过 99% 显著性检验的 links 去掉  --##--##--##--##
t_True_days = np.sum(np.sum(t_True2, axis=1), axis=0)
print(" Heat days  ",np.min(t_True_days), np.max(t_True_days))
h_Positive_days = np.sum(np.sum(h_Positive, axis=1), axis=0)
print(" H Positive days  ",np.min(h_Positive_days), np.max(h_Positive_days))
del t_True2, h_Positive

f0 = xr.open_dataset("/public/home/fcai/extreme1_AT/NC2_1980_2010/1_threshold/syn_shuffle_significant_test1200.nc")
matrix99 = f0.matrix99
matrix4D_99 = np.zeros((np.shape(t0)[1],np.shape(t0)[2], np.shape(t0)[1],np.shape(t0)[2]))




for blat in range(np.shape(t0)[1]):
  print('blat = ',blat)
  for blon in range(np.shape(t0)[2]):
    b_days = int(h_Positive_days[blat,blon])
    
    for a_daysT in range(0,800):
      concurrent_days = matrix99[a_daysT, a_daysT]
      matrix4D_99[:,:,blat,blon][np.where(t_True_days[:,:]==a_daysT)] = concurrent_days
      del concurrent_days
    del b_days

print(" matrix4D_99  ",np.min(matrix4D_99), np.max(matrix4D_99))




##--##--##--##--  同步网络矩阵, 存入nc文件  --##--##--##--##
matrix4D_99_array0 = xr.DataArray(data=matrix4D_99, dims=['lat1', 'lon1', 'lat2', 'lon2'],
                                coords={'lat1':lat.data,'lon1':lon.data,
                                        'lat2':lat.data,'lon2':lon.data})
ds0 = xr.Dataset(data_vars=dict(matrix4D_99=matrix4D_99_array0))
ds0.to_netcdf("/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/part3_heathigh_daily/shuffle99_Tmax90_H500Positive_lag0.nc")
ds0.close()




