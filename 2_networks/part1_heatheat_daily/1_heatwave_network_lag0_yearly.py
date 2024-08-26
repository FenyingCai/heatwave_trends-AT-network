import numpy as np
import xarray as xr

##--##--##--##--  读取 15 天滑动气温 10%, 90% 阈值 (JJA)   --##--##--##--##
ds00 = xr.open_dataset("/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/1_threshold/era5_Tmax_th90_JJA_1981_2010.nc")
t_th90 = ds00.air_th90.loc[:,85:10,:][:,::-1,:]


##--##--##--##--  读取 ERA5 2米气温数据  --##--##--##--##
ds00 = xr.open_dataset("/home/ys17-23/Extension/ERA5/ERA5-daily/surface/t2m_max/t2m.2021-05.daily.nc")
t0 = ds00.t2m.loc[:,85:10,:][:,::-1,:]; ds00.close()
lat = t0.latitude; lon = t0.longitude
print(np.shape(t0))


##--##--##--##--                   --##--##--##--##
##--##--##--##--  1979-2022, 40年  --##--##--##--##
##--##--##--##--                   --##--##--##--##
for iyear in range(44):
  print('iyear = ',iyear+1979)
  
  t = np.zeros((92, np.shape(t0)[1],np.shape(t0)[2]))
  ds2 = xr.open_dataset("/home/ys17-23/Extension/ERA5/ERA5-daily/surface/t2m_max/t2m."+str(iyear+1979)+"-06.daily.nc")
  t[0:30,:,:] = ds2.t2m.loc[:,85:10,:][:,::-1,:]
  ds3 = xr.open_dataset("/home/ys17-23/Extension/ERA5/ERA5-daily/surface/t2m_max/t2m."+str(iyear+1979)+"-07.daily.nc")
  t[30:61,:,:] = ds3.t2m.loc[:,85:10,:][:,::-1,:]
  ds4 = xr.open_dataset("/home/ys17-23/Extension/ERA5/ERA5-daily/surface/t2m_max/t2m."+str(iyear+1979)+"-08.daily.nc")
  t[61:92,:,:] = ds4.t2m.loc[:,85:10,:][:,::-1,:]
  ds2.close(); ds3.close(); ds4.close()



  ##--##--##--##--  对于某一时刻，气温是否超过90%阈值, 或低于10%阈值  --##--##--##--##
  t_True = np.zeros((92, np.shape(t0)[1],np.shape(t0)[2]))
  t_False = np.zeros((92, np.shape(t0)[1],np.shape(t0)[2]))
  t_True[:,:,:] = np.array(t) >= np.array(t_th90)
  del t; print(np.sum(t_True), np.size(t_True))


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
  
  del t_True; print(np.sum(t_True2), np.size(t_True2))



  ##--##--##--##--  热浪事件同步网络  --##--##--##--##
  hw_networks0 = np.zeros((np.shape(t0)[1],np.shape(t0)[2], np.shape(t0)[1],np.shape(t0)[2]))
  # print(hw_networks0.size * hw_networks0.itemsize / 1073741824.0)

  for alat in range(np.shape(t0)[1]):
    print('alat = ',alat)
    for alon in range(np.shape(t0)[2]):
      t_True_3D = np.tile(t_True2[:,alat,alon], (np.shape(t0)[1],np.shape(t0)[2],1)).transpose(2,0,1)  ## day|92* lat|66* lon|360
      # if np.sum(np.sum(np.sum(t_True_3D[5,:,:])))>0: print(np.sum(np.sum(np.sum(t_True_3D[5,:,:]))))
      
      concurrent_days = np.sum(t_True_3D[:,:,:]+t_True2[:,:,:]==2,axis=0)
      hw_networks0[alat,alon,:,:] = concurrent_days
      # hw_networks0[:,:,alat,alon] = concurrent_days
      del concurrent_days, t_True_3D
  del t_True2
  print(" Network ", np.nanmax(hw_networks0), np.nanmin(hw_networks0))

  

  ##--##--##--##--  同步网络矩阵, 存入nc文件  --##--##--##--##
  hw_networks_array0= xr.DataArray(data=hw_networks0.data, dims=['lat1', 'lon1', 'lat2', 'lon2'],
                                  coords={'lat1':lat.data,'lon1':lon.data,
                                          'lat2':lat.data,'lon2':lon.data})
  ds0 = xr.Dataset(data_vars=dict(networks0=hw_networks_array0))
  ds0.to_netcdf("/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/2_networks/part1_heatheat_daily/Networks_Tmax90_lag0_part1_heatheat_"+str(iyear+1979)+".nc")
  ds0.close()
  del hw_networks0, hw_networks_array0




