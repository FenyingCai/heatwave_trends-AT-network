import numpy as np
import xarray as xr

##--##--##--##--  read lat, lon  --##--##--##--##
ds00 = xr.open_dataset("/public/home/fcai/data-observation/ERA5_daily/t2m_max/t2m.2021-05.daily.nc")
t0 = ds00.t2m.loc[:,90:0,:][:,::-1,:]; ds00.close()
lat = t0.latitude; lon = t0.longitude
print(lat)
print(np.shape(t0))
distance = np.zeros((np.shape(t0)[1],np.shape(t0)[2], np.shape(t0)[1],np.shape(t0)[2]))



##--##--##--##--  distance between i and j  --##--##--##--##
lat_2D = np.zeros((np.shape(t0)[1],np.shape(t0)[2]))
lon_2D = np.zeros((np.shape(t0)[1],np.shape(t0)[2]))
for i in range(np.shape(t0)[1]):
  lon_2D[i,:] = lon[:]
for j in range(np.shape(t0)[2]):
  lat_2D[:,j] = lat[:]


for i in range(np.shape(t0)[1]):
  if i%10 == 0: print("i = ",i)
  for j in range(np.shape(t0)[2]):
    y_r = 6371.0 * (lat_2D[:,:] - lat_2D[i,j])*3.1416/180.0    
    x_r1 = 6371.0 * np.abs(lon_2D[:,:] - lon_2D[i,j])*3.1416/180.0 * np.cos((lat_2D[:,:]+lat_2D[i,j])/2.0 *3.1416/180.0)
    x_r2 = 6371.0 * (360 - np.abs(lon_2D[:,:]-lon_2D[i,j]))*3.1416/180.0 * np.cos((lat_2D[:,:]+lat_2D[i,j])/2.0 *3.1416/180.0)
    x_r = np.minimum(x_r1, x_r2)
    r1 = np.sqrt(y_r*y_r + x_r*x_r)
    distance[i,j,:,:] = r1
    del y_r, r1
print("%e"%np.nansum(distance))


print(np.max(x_r1), np.min(x_r1))
print(np.max(x_r2), np.min(x_r2))
print(np.max(x_r), np.min(x_r))
print(np.max(distance), np.min(distance))



##--##--##--##--  to nc file  --##--##--##--##
distance_array0= xr.DataArray(data=distance, dims=['lat1', 'lon1', 'lat2', 'lon2'],
                                coords={'lat1':lat.data,
                                        'lon1':lon.data,
                                        'lat2':lat.data,
                                        'lon2':lon.data})
ds0 = xr.Dataset(data_vars=dict(distance=distance_array0))
ds0.to_netcdf("/public/home/fcai/extreme1_AT/NC2_1980_2010/5_response/11_distance_test/distance_i_j.nc")
ds0.close()




