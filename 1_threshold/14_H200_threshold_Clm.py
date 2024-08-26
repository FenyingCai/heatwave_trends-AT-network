import numpy as np
import xarray as xr

##--##--##--##--  读取 ERA5 2米气温数据  --##--##--##--##
ds0 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/zg/zg.2021-05.daily.nc")
h0 = ds0.z.loc[:,200, 90:0,:]
print(np.shape(h0))
h = np.zeros((30, 92, np.shape(h0)[1],np.shape(h0)[2]))

for iyear in range(30):   ## 1981-2010
  print('iyear = ',iyear+1981)

  ds2 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/zg/zg."+str(iyear+1981)+"-06.daily.nc")
  h[iyear,0:30,:,:] = ds2.z.loc[:,200, 90:0,:]
  ds3 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/zg/zg."+str(iyear+1981)+"-07.daily.nc")
  h[iyear,30:61,:,:] = ds3.z.loc[:,200, 90:0,:]
  ds4 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/zg/zg."+str(iyear+1981)+"-08.daily.nc")
  h[iyear,61:92,:,:] = ds4.z.loc[:,200, 90:0,:]
  ds2.close(); ds3.close(); ds4.close()

##--##--##--##--  计算夏天一倍标准差阈值  --##--##--##--##
h_clm = np.nanmean(h, axis=0); del h
np.shape(h_clm)



##--##--##--##--  写进nc文件内  --##--##--##--##
h_clm = xr.DataArray(data=h_clm.data, dims=['day','lat','lon'],
                      coords={'day':np.arange(1,93,1), 'lat':h0.latitude.data, 'lon':h0.longitude.data})
ds2 = xr.Dataset(data_vars=dict(h_clm=h_clm))
ds2.to_netcdf("/public/home/fcai/extreme1_AT/NC2_1980_2010/1_threshold/era5_H200_clm_JJA_1981_2010.nc")
ds2.close()


