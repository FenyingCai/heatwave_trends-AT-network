
import numpy as np
import xarray as xr

##--##--##--##--  read Tmax 90th threshold   --##--##--##--##
ds00 = xr.open_dataset("/public/home/fcai/extreme1_AT/NC2_1980_2010/1_threshold/era5_Tmax_th90_JJA_1981_2010.nc")
t_th90 = ds00.air_th90.loc[:,90:0,:][:,::-1,:]
ds01 = xr.open_dataset("/public/home/fcai/data0/topography/topo_adapt_era5_air2m.nc")
topo = ds01.topo.loc[0:90,:]; ds01.close()


##--##--##--##--  read Tmax  --##--##--##--##
ds00 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/t2m_max/t2m.2021-05.daily.nc")
t0 = ds00.t2m.loc[:,90:0,:][:,::-1,:]; ds00.close()
lat = t0.latitude;  lon = t0.longitude;  print(np.shape(t0))
t = np.zeros((44, 92, np.shape(t0)[1],np.shape(t0)[2]))

for iyear in range(44):
  # print('iyear = ',iyear+1979)
  ds2 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/t2m_max/t2m."+str(iyear+1979)+"-06.daily.nc")
  t[iyear,0:30,:,:] = ds2.t2m.loc[:,90:0,:][:,::-1,:]
  ds3 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/t2m_max/t2m."+str(iyear+1979)+"-07.daily.nc")
  t[iyear,30:61,:,:] = ds3.t2m.loc[:,90:0,:][:,::-1,:]
  ds4 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/t2m_max/t2m."+str(iyear+1979)+"-08.daily.nc")
  t[iyear,61:92,:,:] = ds4.t2m.loc[:,90:0,:][:,::-1,:]
  ds2.close(); ds3.close(); ds4.close()



##--##--##--##--  Identify heatwaves (>90th, >=3days)  --##--##--##--##
t_True = np.zeros((44, 92, np.shape(t0)[1],np.shape(t0)[2]))
for iyear in range(44):
  t_True[iyear,:,:,:] = np.array(t[iyear,:,:,:]) >= np.array(t_th90)
print(np.sum(t_True), np.size(t_True))
del t, t_th90


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
del t_True; print(np.sum(t_True2), np.size(t_True2))



##--##--##--##--  mask ocean  --##--##--##--##
topo_4D = np.zeros_like(t_True2)
for iyear in range(44):
  # print(' iyear = ',iyear)
  for iday in range(92):
    topo_4D[iyear,iday,:,:] = topo.loc[0:90].copy()

t_True2[topo_4D<=0.0] = 0.0
del topo_4D
print(np.sum(t_True2), np.size(t_True2))


##--##--##--##--  regional heatwave grid cells: PC (year|44 × day|92)  --##--##--##--##
western_lon = np.zeros(36,dtype=int)
western_lon[0:10] = np.linspace(350,359,10)
western_lon[10:36] = np.linspace(0,25,26)
print(western_lon)
PC_WE = np.sum(np.sum(np.array(t_True2[:,:,35:55+1,western_lon]),axis=3),axis=2)
print('PC_WE', np.shape(PC_WE), np.max(PC_WE), np.min(PC_WE), np.nanpercentile(PC_WE, 90))



##--##--##--##--  heatwave grid cells > 90th  --##--##--##--##
PC_index = np.array(np.where(PC_WE>=np.nanpercentile(PC_WE, 90)))
print('   90% percentile,   ', np.nanpercentile(PC_WE, 90))
print('   Days,   ', np.shape(PC_index))




##--##--##--##--##--##--##--##--                       --##--##--##--##--##--##--##--##
##--##--##--##--##--##--##--##--   read data (single)   --##--##--##--##--##--##--##--##
t = np.zeros((44, 92, np.shape(t0)[1],np.shape(t0)[2]))
u = np.zeros((44, 92, np.shape(t0)[1],np.shape(t0)[2]))
v = np.zeros((44, 92, np.shape(t0)[1],np.shape(t0)[2]))
h = np.zeros((44, 92, np.shape(t0)[1],np.shape(t0)[2]))

for iyear in range(44):
  # print('iyear = ',iyear+1979)
  ds2 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/t2m/t2m."+str(iyear+1979)+"-06.daily.nc")
  t[iyear,0:30,:,:] = ds2.t2m.loc[:,90:0,:][:,::-1,:]
  ds3 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/t2m/t2m."+str(iyear+1979)+"-07.daily.nc")
  t[iyear,30:61,:,:] = ds3.t2m.loc[:,90:0,:][:,::-1,:]
  ds4 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/t2m/t2m."+str(iyear+1979)+"-08.daily.nc")
  t[iyear,61:92,:,:] = ds4.t2m.loc[:,90:0,:][:,::-1,:]
  ds2.close(); ds3.close(); ds4.close()
  
  ds2 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/uwind/uwind."+str(iyear+1979)+"-06.daily.nc")
  u[iyear,0:30,:,:]  = ds2.u.loc[:,500,90:0,:][:,::-1,:]
  ds3 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/uwind/uwind."+str(iyear+1979)+"-07.daily.nc")
  u[iyear,30:61,:,:] = ds3.u.loc[:,500,90:0,:][:,::-1,:]
  ds4 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/uwind/uwind."+str(iyear+1979)+"-08.daily.nc")
  u[iyear,61:92,:,:] = ds4.u.loc[:,500,90:0,:][:,::-1,:]
  ds2.close(); ds3.close(); ds4.close()
  
  ds5 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/vwind/vwind."+str(iyear+1979)+"-06.daily.nc")
  v[iyear,0:30,:,:]  = ds5.v.loc[:,500,90:0,:][:,::-1,:]
  ds6 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/vwind/vwind."+str(iyear+1979)+"-07.daily.nc")
  v[iyear,30:61,:,:] = ds6.v.loc[:,500,90:0,:][:,::-1,:]
  ds7 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/vwind/vwind."+str(iyear+1979)+"-08.daily.nc")
  v[iyear,61:92,:,:] = ds7.v.loc[:,500,90:0,:][:,::-1,:]
  ds5.close(); ds6.close(); ds7.close()

  ds5 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/zg/zg."+str(iyear+1979)+"-06.daily.nc")
  h[iyear,0:30,:,:]  = ds5.z.loc[:,500,90:0,:][:,::-1,:]
  ds6 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/zg/zg."+str(iyear+1979)+"-07.daily.nc")
  h[iyear,30:61,:,:] = ds6.z.loc[:,500,90:0,:][:,::-1,:]
  ds7 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/zg/zg."+str(iyear+1979)+"-08.daily.nc")
  h[iyear,61:92,:,:] = ds7.z.loc[:,500,90:0,:][:,::-1,:]
  ds5.close(); ds6.close(); ds7.close()




##--##--##--##--   Composite   --##--##--##--##
t_WE = np.zeros((np.shape(PC_index)[1], np.shape(t0)[1],np.shape(t0)[2]))
u_WE = np.zeros((np.shape(PC_index)[1], np.shape(t0)[1],np.shape(t0)[2]))
v_WE = np.zeros((np.shape(PC_index)[1], np.shape(t0)[1],np.shape(t0)[2]))
h_WE = np.zeros((np.shape(PC_index)[1], np.shape(t0)[1],np.shape(t0)[2]))
for i in range(np.shape(PC_index)[1]):
  t_WE[i,:,:] = t[PC_index[0,i],PC_index[1,i],:,:].copy()
  u_WE[i,:,:] = u[PC_index[0,i],PC_index[1,i],:,:].copy()
  v_WE[i,:,:] = v[PC_index[0,i],PC_index[1,i],:,:].copy()
  h_WE[i,:,:] = h[PC_index[0,i],PC_index[1,i],:,:].copy()

day1 = np.linspace(1,np.shape(PC_index)[1],np.shape(PC_index)[1])
t_WE = xr.DataArray(data=t_WE.data, dims=['day','lat','lon'], coords={'day':day1, 'lat':lat.data, 'lon':lon.data})
u_WE = xr.DataArray(data=u_WE.data, dims=['day','lat','lon'], coords={'day':day1, 'lat':lat.data, 'lon':lon.data})
v_WE = xr.DataArray(data=v_WE.data, dims=['day','lat','lon'], coords={'day':day1, 'lat':lat.data, 'lon':lon.data})
h_WE = xr.DataArray(data=h_WE.data, dims=['day','lat','lon'], coords={'day':day1, 'lat':lat.data, 'lon':lon.data})





##--##--##--##--  All days --##--##--##--##
t_4D = np.zeros((44*92, np.shape(t0)[1],np.shape(t0)[2]))
u_4D = np.zeros((44*92, np.shape(t0)[1],np.shape(t0)[2]))
v_4D = np.zeros((44*92, np.shape(t0)[1],np.shape(t0)[2]))
h_4D = np.zeros((44*92, np.shape(t0)[1],np.shape(t0)[2]))
for iyear in range(44):
  t_4D[iyear*92:iyear*92+92,:,:] = t[iyear,:,:,:].copy()
  u_4D[iyear*92:iyear*92+92,:,:] = u[iyear,:,:,:].copy()
  v_4D[iyear*92:iyear*92+92,:,:] = v[iyear,:,:,:].copy()
  h_4D[iyear*92:iyear*92+92,:,:] = h[iyear,:,:,:].copy()
t_4D = xr.DataArray(data=t_4D.data, dims=['day2','lat','lon'], coords={'day2':np.linspace(1,92*44,92*44), 'lat':lat.data, 'lon':lon.data})
u_4D = xr.DataArray(data=u_4D.data, dims=['day2','lat','lon'], coords={'day2':np.linspace(1,92*44,92*44), 'lat':lat.data, 'lon':lon.data})
v_4D = xr.DataArray(data=v_4D.data, dims=['day2','lat','lon'], coords={'day2':np.linspace(1,92*44,92*44), 'lat':lat.data, 'lon':lon.data})
h_4D = xr.DataArray(data=h_4D.data, dims=['day2','lat','lon'], coords={'day2':np.linspace(1,92*44,92*44), 'lat':lat.data, 'lon':lon.data})





##--##--##--##--   Composited differences   --##--##--##--##
t_diff = np.mean(t_WE.data, axis=0) - np.mean(t_4D.data, axis=0)
u_diff = np.mean(u_WE.data, axis=0) - np.mean(u_4D.data, axis=0)
v_diff = np.mean(v_WE.data, axis=0) - np.mean(v_4D.data, axis=0)
h_diff = np.mean(h_WE.data, axis=0) - np.mean(h_4D.data, axis=0)
print(" t_diff max ", np.nanmax(t_diff), np.nanmin(t_diff))

t_diff = xr.DataArray(data=t_diff.data, dims=['lat','lon'], coords={'lat':lat.data, 'lon':lon.data})
u_diff = xr.DataArray(data=u_diff.data, dims=['lat','lon'], coords={'lat':lat.data, 'lon':lon.data})
v_diff = xr.DataArray(data=v_diff.data, dims=['lat','lon'], coords={'lat':lat.data, 'lon':lon.data})
h_diff = xr.DataArray(data=h_diff.data, dims=['lat','lon'], coords={'lat':lat.data, 'lon':lon.data})




ds2 = xr.Dataset(data_vars=dict(t_WE=t_WE, u_WE=u_WE, v_WE=v_WE, h_WE=h_WE,\
                                t=t_4D, u=u_4D, v=v_4D, h=h_4D,\
                                t_diff=t_diff, u_diff=u_diff, v_diff=v_diff, h_diff=h_diff))
ds2.to_netcdf("/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/case_western_Europe/era5_WE_composite_T2m_UVH500.nc")
ds2.close()
del t_WE, u_WE, v_WE, h_WE, t_4D, u_4D, v_4D, h_4D










##--##--##--##--##--##--##--##--                        --##--##--##--##--##--##--##--##
##--##--##--##--##--##--##--##--   read data (Levels)   --##--##--##--##--##--##--##--##
u = np.zeros((44, 92, 5, np.shape(t0)[1],np.shape(t0)[2]))
v = np.zeros((44, 92, 5, np.shape(t0)[1],np.shape(t0)[2]))
h = np.zeros((44, 92, 5, np.shape(t0)[1],np.shape(t0)[2]))

for iyear in range(44):
  # print('iyear = ',iyear+1979)
  ds2 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/uwind/uwind."+str(iyear+1979)+"-06.daily.nc")
  u[iyear,0:30,:,:]  = ds2.u.loc[:,[600,500,400,300,200],90:0,:][:,:,::-1,:]
  ds3 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/uwind/uwind."+str(iyear+1979)+"-07.daily.nc")
  u[iyear,30:61,:,:] = ds3.u.loc[:,[600,500,400,300,200],90:0,:][:,:,::-1,:]
  ds4 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/uwind/uwind."+str(iyear+1979)+"-08.daily.nc")
  u[iyear,61:92,:,:] = ds4.u.loc[:,[600,500,400,300,200],90:0,:][:,:,::-1,:]
  ds2.close(); ds3.close(); ds4.close()
  
  ds5 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/vwind/vwind."+str(iyear+1979)+"-06.daily.nc")
  v[iyear,0:30,:,:]  = ds5.v.loc[:,[600,500,400,300,200],90:0,:][:,:,::-1,:]
  ds6 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/vwind/vwind."+str(iyear+1979)+"-07.daily.nc")
  v[iyear,30:61,:,:] = ds6.v.loc[:,[600,500,400,300,200],90:0,:][:,:,::-1,:]
  ds7 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/vwind/vwind."+str(iyear+1979)+"-08.daily.nc")
  v[iyear,61:92,:,:] = ds7.v.loc[:,[600,500,400,300,200],90:0,:][:,:,::-1,:]
  ds5.close(); ds6.close(); ds7.close()

  ds5 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/zg/zg."+str(iyear+1979)+"-06.daily.nc")
  h[iyear,0:30,:,:]  = ds5.z.loc[:,[600,500,400,300,200],90:0,:][:,:,::-1,:]
  ds6 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/zg/zg."+str(iyear+1979)+"-07.daily.nc")
  h[iyear,30:61,:,:] = ds6.z.loc[:,[600,500,400,300,200],90:0,:][:,:,::-1,:]
  ds7 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/zg/zg."+str(iyear+1979)+"-08.daily.nc")
  h[iyear,61:92,:,:] = ds7.z.loc[:,[600,500,400,300,200],90:0,:][:,:,::-1,:]
  ds5.close(); ds6.close(); ds7.close()




##--##--##--##--   Composite   --##--##--##--##
u_WE = np.zeros((np.shape(PC_index)[1], 5, np.shape(t0)[1],np.shape(t0)[2]))
v_WE = np.zeros((np.shape(PC_index)[1], 5, np.shape(t0)[1],np.shape(t0)[2]))
h_WE = np.zeros((np.shape(PC_index)[1], 5, np.shape(t0)[1],np.shape(t0)[2]))
for i in range(np.shape(PC_index)[1]):
  u_WE[i,:,:,:] = u[PC_index[0,i],PC_index[1,i],:,:,:].copy()
  v_WE[i,:,:,:] = v[PC_index[0,i],PC_index[1,i],:,:,:].copy()
  h_WE[i,:,:,:] = h[PC_index[0,i],PC_index[1,i],:,:,:].copy()

level = np.array([600,500,400,300,200])
day1 = np.linspace(1,np.shape(PC_index)[1],np.shape(PC_index)[1])
u_WE = xr.DataArray(data=u_WE.data, dims=['day','level','lat','lon'], coords={'day':day1, 'level':level, 'lat':lat.data, 'lon':lon.data})
v_WE = xr.DataArray(data=v_WE.data, dims=['day','level','lat','lon'], coords={'day':day1, 'level':level, 'lat':lat.data, 'lon':lon.data})
h_WE = xr.DataArray(data=h_WE.data, dims=['day','level','lat','lon'], coords={'day':day1, 'level':level, 'lat':lat.data, 'lon':lon.data})




##--##--##--##--  All days  --##--##--##--##
u_4D = np.zeros((44*92, 5, np.shape(t0)[1],np.shape(t0)[2]))
v_4D = np.zeros((44*92, 5, np.shape(t0)[1],np.shape(t0)[2]))
h_4D = np.zeros((44*92, 5, np.shape(t0)[1],np.shape(t0)[2]))
for iyear in range(44):
  u_4D[iyear*92:iyear*92+92,:,:,:] = u[iyear,:,:,:,:].copy()
  v_4D[iyear*92:iyear*92+92,:,:,:] = v[iyear,:,:,:,:].copy()
  h_4D[iyear*92:iyear*92+92,:,:,:] = h[iyear,:,:,:,:].copy()
u_4D = xr.DataArray(data=u_4D.data, dims=['day2','level','lat','lon'], coords={'day2':np.linspace(1,92*44,92*44), 'level':level, 'lat':lat.data, 'lon':lon.data})
v_4D = xr.DataArray(data=v_4D.data, dims=['day2','level','lat','lon'], coords={'day2':np.linspace(1,92*44,92*44), 'level':level, 'lat':lat.data, 'lon':lon.data})
h_4D = xr.DataArray(data=h_4D.data, dims=['day2','level','lat','lon'], coords={'day2':np.linspace(1,92*44,92*44), 'level':level, 'lat':lat.data, 'lon':lon.data})


##--##--##--##--   Composited differences   --##--##--##--##
u_diff = np.mean(u_WE.data, axis=0) - np.mean(u_4D.data, axis=0)
v_diff = np.mean(v_WE.data, axis=0) - np.mean(v_4D.data, axis=0)
h_diff = np.mean(h_WE.data, axis=0) - np.mean(h_4D.data, axis=0)

u_diff = xr.DataArray(data=u_diff.data, dims=['level','lat','lon'], coords={'level':level, 'lat':lat.data, 'lon':lon.data})
v_diff = xr.DataArray(data=v_diff.data, dims=['level','lat','lon'], coords={'level':level, 'lat':lat.data, 'lon':lon.data})
h_diff = xr.DataArray(data=h_diff.data, dims=['level','lat','lon'], coords={'level':level, 'lat':lat.data, 'lon':lon.data})



ds2 = xr.Dataset(data_vars=dict(u_WE=u_WE, v_WE=v_WE, h_WE=h_WE,\
                                u=u_4D, v=v_4D, h=h_4D,\
                                u_diff=u_diff, v_diff=v_diff, h_diff=h_diff))
ds2.to_netcdf("/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/case_western_Europe/era5_WE_composite_UVH_5levels.nc")
ds2.close()









