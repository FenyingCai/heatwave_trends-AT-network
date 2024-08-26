import numpy as np
import xarray as xr

##--##--##--##--  读取 15 天滑动气温 10%, 90% 阈值 (JJA)   --##--##--##--##
ds00 = xr.open_dataset("/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/1_threshold/era5_Tmax_th90_JJA_1981_2010.nc")
t_th90 = ds00.air_th90.loc[:,85:10,:][:,::-1,:]
ds01 = xr.open_dataset("/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/1_threshold/era5_H500_std_JJA_1981_2010.nc")
h_std = ds01.h_std.loc[:,85:10,:][:,::-1,:]
ds02 = xr.open_dataset("/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/1_threshold/era5_H500_clm_JJA_1981_2010.nc")
h_clm = ds02.h_clm.loc[:,85:10,:][:,::-1,:]


##--##--##--##--  读取 ERA5 2米气温数据  --##--##--##--##
ds00 = xr.open_dataset("/home/ys17-23/Extension/ERA5/ERA5-daily/surface/t2m_max/t2m.2021-05.daily.nc")
t0 = ds00.t2m.loc[:,85:10,:][:,::-1,:]; ds00.close()
lat = t0.latitude; lon = t0.longitude
print(np.shape(t0))




##--##--##--##--                   --##--##--##--##
##--##--##--##--  1979-2022, 40年  --##--##--##--##
##--##--##--##--                   --##--##--##--##

for iyear in range(44):
  print("                   ")
  print('iyear = ',iyear+1979)
  print("                   ")
  
  
  ##--##--##--##--  Tmax, 90th  --##--##--##--##
  t = np.zeros((92, np.shape(t0)[1],np.shape(t0)[2]))
  ds2 = xr.open_dataset("/home/ys17-23/Extension/ERA5/ERA5-daily/surface/t2m_max/t2m."+str(iyear+1979)+"-06.daily.nc")
  t[0:30,:,:] = ds2.t2m.loc[:,85:10,:][:,::-1,:]
  ds3 = xr.open_dataset("/home/ys17-23/Extension/ERA5/ERA5-daily/surface/t2m_max/t2m."+str(iyear+1979)+"-07.daily.nc")
  t[30:61,:,:] = ds3.t2m.loc[:,85:10,:][:,::-1,:]
  ds4 = xr.open_dataset("/home/ys17-23/Extension/ERA5/ERA5-daily/surface/t2m_max/t2m."+str(iyear+1979)+"-08.daily.nc")
  t[61:92,:,:] = ds4.t2m.loc[:,85:10,:][:,::-1,:]
  ds2.close(); ds3.close(); ds4.close()
  
  t_True = np.zeros((92, np.shape(t0)[1],np.shape(t0)[2]))
  t_True[:,:,:] = np.array(t) >= np.array(t_th90)
  del t; print("Heat  days ",np.sum(t_True), np.size(t_True))

  ##--##--##--##--  热浪事件：连续三天超过90%阈值； 偏冷事件：或低于10%阈值  --##--##--##--##
  t_True2 = t_True.copy()
  t_True2[2:-2,:,:][t_True[:-4,:,:]+t_True[1:-3,:,:]+t_True[3:-1,:,:]+t_True[4:,:,:]<=1] = 0.0
  t_True2[2:-2,:,:][t_True[1:-3,:,:]+t_True[3:-1,:,:]==0] = 0.0
  t_True2[2:-2,:,:][(t_True[1:-3,:,:]+t_True[3:-1,:,:]==1)&(t_True[3:-1,:,:]+t_True[4:,:,:]<=1)&(t_True[:-4,:,:]+t_True[1:-3,:,:]<=1)] = 0.0

  t_True2[0,:,:][t_True[1,:,:]+t_True[2,:,:]<=1] = 0.0
  t_True2[1,:,:][t_True[0,:,:]+t_True[2,:,:]+t_True[3,:,:]<=1] = 0.0
  t_True2[1,:,:][(t_True[0,:,:]==1)&(t_True[2,:,:]==0)&(t_True[3,:,:]==1)] = 0.0
  t_True2[-1,:,:][t_True[-2,:,:]+t_True[-3,:,:]<=1] = 0.0
  t_True2[-2,:,:][t_True[-1,:,:]+t_True[-3,:,:]+t_True[-4,:,:]<=1] = 0.0
  t_True2[-2,:,:][(t_True[-1,:,:]==1)&(t_True[-3,:,:]==0)&(t_True[-4,:,:]==1)] = 0.0
  del t_True
  




  
  ##--##--##--##--  H500, >1std  --##--##--##--##
  h = np.zeros((92, np.shape(t0)[1],np.shape(t0)[2]))
  ds2 = xr.open_dataset("/home/ys17-23/Extension/ERA5/ERA5-daily/pressure/zg/zg."+str(iyear+1979)+"-06.daily.nc")
  h[0:30,:,:] = ds2.z.loc[:,500,85:10,:][:,::-1,:]
  ds3 = xr.open_dataset("/home/ys17-23/Extension/ERA5/ERA5-daily/pressure/zg/zg."+str(iyear+1979)+"-07.daily.nc")
  h[30:61,:,:] = ds3.z.loc[:,500,85:10,:][:,::-1,:]
  ds4 = xr.open_dataset("/home/ys17-23/Extension/ERA5/ERA5-daily/pressure/zg/zg."+str(iyear+1979)+"-08.daily.nc")
  h[61:92,:,:] = ds4.z.loc[:,500,85:10,:][:,::-1,:]
  ds2.close(); ds3.close(); ds4.close()
  h_ano = np.array(h) - np.array(h_clm); del h

  h_Positive = np.zeros((92, np.shape(t0)[1],np.shape(t0)[2]))
  h_Positive[:,:,:] = np.array(h_ano) >= np.array(h_std); del h_ano
  print(" High  pressure  days ",np.sum(h_Positive), np.size(h_Positive))
  




  ##--##--##--##--  热浪事件同步网络  --##--##--##--##
  hw_networks0 = np.zeros((np.shape(t0)[1],np.shape(t0)[2], np.shape(t0)[1],np.shape(t0)[2]))
  # print(hw_networks0.size * hw_networks0.itemsize / 1073741824.0)

  for alat in range(np.shape(t0)[1]):
    print('alat = ',alat)
    for alon in range(np.shape(t0)[2]):

      t_True_3D = np.tile(t_True2[:,alat,alon], (np.shape(t0)[1],np.shape(t0)[2],1)).transpose(2,0,1)   ## day|92* lat|66* lon|360
      
      concurrent_days = np.sum(t_True_3D[:,:,:]+h_Positive[:,:,:]==2,axis=0)
      hw_networks0[alat,alon,:,:] = concurrent_days

      del  t_True_3D, concurrent_days
  del t_True2, h_Positive
  print("network ", np.max(hw_networks0), np.min(hw_networks0))




  ##--##--##--##--  同步网络矩阵, 存入nc文件  --##--##--##--##
  hw_networks_array0= xr.DataArray(data=hw_networks0.data, dims=['lat1', 'lon1', 'lat2', 'lon2'],
                                  coords={'lat1':lat.data,'lon1':lon.data,
                                          'lat2':lat.data,'lon2':lon.data})
  ds0 = xr.Dataset(data_vars=dict(networks0=hw_networks_array0))
  ds0.to_netcdf("/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/2_networks/part3_heathigh_daily/Networks_Tmax90_H500Positive_lag0_"+str(iyear+1979)+".nc")
  ds0.close()
  del hw_networks0, hw_networks_array0




