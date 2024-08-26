load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin
;;********************************** read lat, lon ************************************
f00 = addfile("data/topo_adapt_era5_air2m.nc","r")
topo0 = f00->topo({10:85},{0:360})
lat = topo0&lat
lon = topo0&lon
dim_0 = dimsizes(topo0)



;;-----------------------   read H500 high-pressure, low-pressure   ----------------------;;
f2 = addfile("data/hgt_ano_44years_JJAyearly.nc","r")
h500Positive_yearly = f2->hgt_Positive_intensity(:,{500},{10:85},:)
h500Positive_yearly@_FillValue = 1.0e20
h500Negative_yearly = f2->hgt_Negative_intensity(:,{500},{10:85},:)
h500Negative_yearly@_FillValue = 1.0e20
printVarSummary(h500Positive_yearly)


;;----------------------- remove zonal mean  ----------------------;;
h500Positive_yearly = dim_rmvmean_n_Wrap(h500Positive_yearly, 2)
h500Negative_yearly = dim_rmvmean_n_Wrap(h500Negative_yearly, 2)

h500_yearly = h500Positive_yearly - h500Negative_yearly
h500_yearly = dim_rmvmean_n_Wrap(h500_yearly, 2)




;;-----------------------  long-term trend  ----------------------;;
trend_Positive = regCoef_n(fspan(1,44,44), h500Positive_yearly, 0,0) / 9.8 * 10.0 ; * 1.0e-2
copy_VarCoords(h500Positive_yearly(0,:,:), trend_Positive)
printMinMax(trend_Positive, 1)
cor_High = escorc_n(fspan(1,44,44), h500Positive_yearly, 0,0)
copy_VarCoords(h500Positive_yearly(0,:,:), cor_High)

cor_High_Posit = where(cor_High.ge.0.0, cor_High, 0.0)
cor_High_Negat = where(cor_High.le.0.0, cor_High, 0.0)
copy_VarCoords(h500Positive_yearly(0,:,:), cor_High_Posit)
copy_VarCoords(h500Positive_yearly(0,:,:), cor_High_Negat)




trend_Negative = regCoef_n(fspan(1,44,44), h500Negative_yearly, 0,0) / 9.8 * 10.0 ; * 1.0e-2
copy_VarCoords(h500Negative_yearly(0,:,:), trend_Negative)
printMinMax(trend_Negative, 1)
cor_Low = escorc_n(fspan(1,44,44), h500Negative_yearly, 0,0)
copy_VarCoords(h500Negative_yearly(0,:,:), cor_Low)

cor_Low_Posit = where(cor_Low.ge.0.0, cor_Low, 0.0)
cor_Low_Negat = where(cor_Low.le.0.0, cor_Low, 0.0)
copy_VarCoords(h500Positive_yearly(0,:,:), cor_Low_Posit)
copy_VarCoords(h500Positive_yearly(0,:,:), cor_Low_Negat)




trend_ = regCoef_n(fspan(1,44,44), h500_yearly, 0,0) ; * 1.0e-2
copy_VarCoords(h500Negative_yearly(0,:,:), trend_)
cor_ = escorc_n(fspan(1,44,44), h500_yearly, 0,0)
copy_VarCoords(h500Negative_yearly(0,:,:), cor_)

cor_Posit = where(cor_.ge.0.0, cor_, 0.0)
cor_Negat = where(cor_.le.0.0, cor_, 0.0)
copy_VarCoords(h500Positive_yearly(0,:,:), cor_Posit)
copy_VarCoords(h500Positive_yearly(0,:,:), cor_Negat)


print("------------ trend -----------")
printMinMax(trend_Positive, 1)
printMinMax(trend_Negative, 1)
printMinMax(trend_, 1)









;;-----------------------                                                  ----------------------;;
;;-----------------------   most-connected points (Correlation networks)   ----------------------;;
f0 = addfile("data/distance_i_j.nc","r")
distance = f0->distance(10:85,:,10:85,:)

f1 = addfile("data/era5_Tmax_heatwave_intensity_yearly.nc","r")
hwi = doubletofloat(f1->heatwave_intensity(:,10:85,:))
hwi@_FillValue = 1.0e20



f5 = addfile("data/CorNetorks_Tmax90th_H500Positive.nc","r")
cor_4D_Positive = f5->cor_4D
dim_2 = dimsizes(cor_4D_Positive(0,0,:,:))
f5 = addfile("data/CorNetorks_Tmax90th_H500Negative.nc","r")
cor_4D_Negative = f5->cor_4D
print(dim_0)
print(dim_2)

WE_points_High1 = new((/dim_0(0),dim_0(1)/), float)
WE_points_Low1 = new((/dim_0(0),dim_0(1)/), float)
WE_points_High1 = 0.0
WE_points_Low1 = 0.0

do ilat = 35-10,55-10-1   ;; western Europe
print(" ilat = "+ilat)
do jlon = 0,25-1
if(.not.ismissing(dim_avg(dim_avg(cor_4D_Positive(ilat,jlon,:,:)))))then

cor_4D_Positive_new = where(distance(ilat,jlon,:,:).le.1000.0, cor_4D_Positive(ilat,jlon,:,:), 0.0)
cor_4D_Positive_new_1D = ndtooned(cor_4D_Positive_new)
indices  = ind_resolve(maxind(cor_4D_Positive_new_1D), dim_2)
WE_points_High1(indices(0,0),indices(0,1)) = WE_points_High1(indices(0,0),indices(0,1)) + 1.0
delete(([/cor_4D_Positive_new,cor_4D_Positive_new_1D,indices/]))

cor_4D_Negative_new = where(distance(ilat,jlon,:,:).gt.1500.0 .and. distance(ilat,jlon,:,:).le.4000.0, cor_4D_Negative(ilat,jlon,:,:), 0.0)
cor_4D_Negative_new_1D = ndtooned(cor_4D_Negative_new)
indices  = ind_resolve(maxind(cor_4D_Negative_new_1D), dim_2)
WE_points_Low1(indices(0,0),indices(0,1)) = WE_points_Low1(indices(0,0),indices(0,1)) + 1.0
delete(([/cor_4D_Negative_new,cor_4D_Negative_new_1D,indices/]))

end if
end do
end do
do ilat = 35-10,55-10-1   ;; western Europe
print(" ilat = "+ilat)
do jlon = 350,360-1
if(.not.ismissing(dim_avg(dim_avg(cor_4D_Positive(ilat,jlon,:,:)))))then

cor_4D_Positive_new = where(distance(ilat,jlon,:,:).le.1000.0, cor_4D_Positive(ilat,jlon,:,:), 0.0)
cor_4D_Positive_new_1D = ndtooned(cor_4D_Positive_new)
indices  = ind_resolve(maxind(cor_4D_Positive_new_1D), dim_2)
WE_points_High1(indices(0,0),indices(0,1)) = WE_points_High1(indices(0,0),indices(0,1)) + 1.0
delete(([/cor_4D_Positive_new,cor_4D_Positive_new_1D,indices/]))

cor_4D_Negative_new = where(distance(ilat,jlon,:,:).gt.1500.0 .and. distance(ilat,jlon,:,:).le.4000.0, cor_4D_Negative(ilat,jlon,:,:), 0.0)
cor_4D_Negative_new_1D = ndtooned(cor_4D_Negative_new)
indices  = ind_resolve(maxind(cor_4D_Negative_new_1D), dim_2)
WE_points_Low1(indices(0,0),indices(0,1)) = WE_points_Low1(indices(0,0),indices(0,1)) + 1.0
delete(([/cor_4D_Negative_new,cor_4D_Negative_new_1D,indices/]))

end if
end do
end do




EE_points_High1 = new((/dim_0(0),dim_0(1)/), float)
EE_points_Low1 = new((/dim_0(0),dim_0(1)/), float)
EE_points_High1 = 0.0
EE_points_Low1 = 0.0

do ilat = 45-10,65-10-1   ;; eastern Europe
print(" ilat = "+ilat)
do jlon = 25,55-1
if(.not.ismissing(dim_avg(dim_avg(cor_4D_Positive(ilat,jlon,:,:)))))then

cor_4D_Positive_new = where(distance(ilat,jlon,:,:).le.1000.0, cor_4D_Positive(ilat,jlon,:,:), 0.0)
cor_4D_Positive_new_1D = ndtooned(cor_4D_Positive_new)
indices  = ind_resolve(maxind(cor_4D_Positive_new_1D), dim_2)
EE_points_High1(indices(0,0),indices(0,1)) = EE_points_High1(indices(0,0),indices(0,1)) + 1.0
delete(([/cor_4D_Positive_new,cor_4D_Positive_new_1D,indices/]))

cor_4D_Negative_new = where(distance(ilat,jlon,:,:).gt.1500.0 .and. distance(ilat,jlon,:,:).le.4000.0, cor_4D_Negative(ilat,jlon,:,:), 0.0)
cor_4D_Negative_new_1D = ndtooned(cor_4D_Negative_new)
indices  = ind_resolve(maxind(cor_4D_Negative_new_1D), dim_2)
EE_points_Low1(indices(0,0),indices(0,1)) = EE_points_Low1(indices(0,0),indices(0,1)) + 1.0
delete(([/cor_4D_Negative_new,cor_4D_Negative_new_1D,indices/]))

end if
end do
end do








;;-----------------------                                                  ----------------------;;
;;-----------------------   most-connected points (Concurrence networks)   ----------------------;;
f5 = addfile("data/Networks_Tmax90_H500Positive_lag0_land_Significant99.nc","r")
cur_4D_Positive = f5->networks0
dim_2 = dimsizes(cor_4D_Positive(0,0,:,:))
f5 = addfile("data/Networks_Tmax90_H500Negative_lag0_land_Significant99.nc","r")
cur_4D_Negative = f5->networks0

WE_points_High2 = new((/dim_0(0),dim_0(1)/), float)
WE_points_Low2 = new((/dim_0(0),dim_0(1)/), float)
WE_points_High2 = 0.0
WE_points_Low2 = 0.0

do ilat = 35-10,55-10-1   ;; western Europe
print(" ilat = "+ilat)
do jlon = 0,25-1
if(.not.ismissing(dim_avg(dim_avg(cur_4D_Positive(ilat,jlon,:,:)))))then

cur_4D_Positive_new = where(distance(ilat,jlon,:,:).le.1000.0, cur_4D_Positive(ilat,jlon,:,:), 0.0)
cur_4D_Positive_new_1D = ndtooned(cur_4D_Positive_new)
indices  = ind_resolve(maxind(cur_4D_Positive_new_1D), dim_2)
WE_points_High2(indices(0,0),indices(0,1)) = WE_points_High2(indices(0,0),indices(0,1)) + 1.0
delete(([/cur_4D_Positive_new,cur_4D_Positive_new_1D,indices/]))

cur_4D_Negative_new = where(distance(ilat,jlon,:,:).gt.1500.0 .and. distance(ilat,jlon,:,:).le.4000.0, cur_4D_Negative(ilat,jlon,:,:), 0.0)
cur_4D_Negative_new_1D = ndtooned(cur_4D_Negative_new)
indices  = ind_resolve(maxind(cur_4D_Negative_new_1D), dim_2)
WE_points_Low2(indices(0,0),indices(0,1)) = WE_points_Low2(indices(0,0),indices(0,1)) + 1.0
delete(([/cur_4D_Negative_new,cur_4D_Negative_new_1D,indices/]))

end if
end do
end do
do ilat = 35-10,55-10-1   ;; western Europe
print(" ilat = "+ilat)
do jlon = 350,360-1
if(.not.ismissing(dim_avg(dim_avg(cur_4D_Positive(ilat,jlon,:,:)))))then

cur_4D_Positive_new = where(distance(ilat,jlon,:,:).le.1000.0, cur_4D_Positive(ilat,jlon,:,:), 0.0)
cur_4D_Positive_new_1D = ndtooned(cur_4D_Positive_new)
indices  = ind_resolve(maxind(cur_4D_Positive_new_1D), dim_2)
WE_points_High2(indices(0,0),indices(0,1)) = WE_points_High2(indices(0,0),indices(0,1)) + 1.0
delete(([/cur_4D_Positive_new,cur_4D_Positive_new_1D,indices/]))

cur_4D_Negative_new = where(distance(ilat,jlon,:,:).gt.1500.0 .and. distance(ilat,jlon,:,:).le.4000.0, cur_4D_Negative(ilat,jlon,:,:), 0.0)
cur_4D_Negative_new_1D = ndtooned(cur_4D_Negative_new)
indices  = ind_resolve(maxind(cur_4D_Negative_new_1D), dim_2)
WE_points_Low2(indices(0,0),indices(0,1)) = WE_points_Low2(indices(0,0),indices(0,1)) + 1.0
delete(([/cur_4D_Negative_new,cur_4D_Negative_new_1D,indices/]))

end if
end do
end do



EE_points_High2 = new((/dim_0(0),dim_0(1)/), float)
EE_points_Low2 = new((/dim_0(0),dim_0(1)/), float)
EE_points_High2 = 0.0
EE_points_Low2 = 0.0

do ilat = 45-10,65-10-1   ;; eastern Europe
print(" ilat = "+ilat)
do jlon = 25,55-1
if(.not.ismissing(dim_avg(dim_avg(cur_4D_Positive(ilat,jlon,:,:)))))then

cur_4D_Positive_new = where(distance(ilat,jlon,:,:).le.1000.0, cur_4D_Positive(ilat,jlon,:,:), 0.0)
cur_4D_Positive_new_1D = ndtooned(cur_4D_Positive_new)
indices  = ind_resolve(maxind(cur_4D_Positive_new_1D), dim_2)
EE_points_High2(indices(0,0),indices(0,1)) = EE_points_High2(indices(0,0),indices(0,1)) + 1.0
delete(([/cur_4D_Positive_new,cur_4D_Positive_new_1D,indices/]))

cur_4D_Negative_new = where(distance(ilat,jlon,:,:).gt.1500.0 .and. distance(ilat,jlon,:,:).le.4000.0, cur_4D_Negative(ilat,jlon,:,:), 0.0)
cur_4D_Negative_new_1D = ndtooned(cur_4D_Negative_new)
indices  = ind_resolve(maxind(cur_4D_Negative_new_1D), dim_2)
EE_points_Low2(indices(0,0),indices(0,1)) = EE_points_Low2(indices(0,0),indices(0,1)) + 1.0
delete(([/cur_4D_Negative_new,cur_4D_Negative_new_1D,indices/]))

end if
end do
end do


copy_VarCoords(h500Negative_yearly(0,:,:), WE_points_High1)
copy_VarCoords(h500Negative_yearly(0,:,:), WE_points_High2)
copy_VarCoords(h500Negative_yearly(0,:,:), WE_points_Low1)
copy_VarCoords(h500Negative_yearly(0,:,:), WE_points_Low2)
copy_VarCoords(h500Negative_yearly(0,:,:), EE_points_High1)
copy_VarCoords(h500Negative_yearly(0,:,:), EE_points_High2)
copy_VarCoords(h500Negative_yearly(0,:,:), EE_points_Low1)
copy_VarCoords(h500Negative_yearly(0,:,:), EE_points_Low2)
print(" points ")

WE_points_High0 = WE_points_High1 + WE_points_High2
WE_points_Low0  = WE_points_Low1 + WE_points_Low2

EE_points_High0 = EE_points_High1 + EE_points_High2
EE_points_Low0 = EE_points_Low1 + EE_points_Low2

low_CA = dim_sum(dim_sum(EE_points_Low0(25-10:75-10, 58:85))) + dim_sum(dim_sum(EE_points_Low0(60-10:75-10, 86:102)))
low_NAf = dim_sum(dim_sum(WE_points_Low0(22-10:38-10, 5:25)))
low_NAt = dim_sum(dim_sum(WE_points_Low0(40-10:60-10, 315:359)))
low_all_EE = dim_sum(dim_sum(EE_points_Low0(:,:)))
low_all_WE = dim_sum(dim_sum(WE_points_Low0(:,:)))
print(" CA "+ low_CA/low_all_EE)
print(" NAf "+ low_NAf/low_all_WE)
print(" NAt "+ low_NAt/low_all_WE)



WE_points_High0 = WE_points_High0 + 0.0000001
WE_points_Low0 = WE_points_Low0 + 0.0000001
EE_points_High0 = EE_points_High0 + 0.0000001
EE_points_Low0 = EE_points_Low0 + 0.0000001



WE_points_HighA = new((/2,dim_0(0),dim_0(1)/), float)
WE_points_LowA = new((/2,dim_0(0),dim_0(1)/), float)
EE_points_HighA = new((/2,dim_0(0),dim_0(1)/), float)
EE_points_LowA = new((/2,dim_0(0),dim_0(1)/), float)


do i=1,dim_0(0)-2
WE_points_HighA(0,i,:) = dim_max_n_Wrap(WE_points_High0(i-1:i+1,:), 0) 
WE_points_LowA(0,i,:) = dim_max_n_Wrap(WE_points_Low0(i-1:i+1,:), 0)
EE_points_HighA(0,i,:) = dim_max_n_Wrap(EE_points_High0(i-1:i+1,:), 0) 
EE_points_LowA(0,i,:) = dim_max_n_Wrap(EE_points_Low0(i-1:i+1,:), 0) 
end do
do j=1,dim_0(1)-2
WE_points_HighA(1,:,j) = dim_max_n_Wrap(WE_points_High0(:,j-1:j+1), 1) 
WE_points_LowA(1,:,j) = dim_max_n_Wrap(WE_points_Low0(:,j-1:j+1), 1)
EE_points_HighA(1,:,j) = dim_max_n_Wrap(EE_points_High0(:,j-1:j+1), 1) 
EE_points_LowA(1,:,j) = dim_max_n_Wrap(EE_points_Low0(:,j-1:j+1), 1) 
end do

WE_points_High = dim_max_n_Wrap(WE_points_HighA, 0)
WE_points_Low = dim_max_n_Wrap(WE_points_LowA, 0)
EE_points_High = dim_max_n_Wrap(EE_points_HighA, 0)
EE_points_Low = dim_max_n_Wrap(EE_points_LowA, 0)

copy_VarCoords(h500Negative_yearly(0,:,:), WE_points_High)
copy_VarCoords(h500Negative_yearly(0,:,:), WE_points_Low)
copy_VarCoords(h500Negative_yearly(0,:,:), EE_points_High)
copy_VarCoords(h500Negative_yearly(0,:,:), EE_points_Low)



WE_points = WE_points_High - WE_points_Low
EE_points = EE_points_High - EE_points_Low
copy_VarCoords(h500Negative_yearly(0,:,:), WE_points)
copy_VarCoords(h500Negative_yearly(0,:,:), EE_points)










;;--------------                                                     -------------;;
;;--------------                    wave train (WE)                  -------------;;
;;--------------                                                     -------------;;
f00 = addfile("data/topo_adapt_era5_air2m.nc","r")
lon_WE = new(76,integer)
lon_WE(0:44) = ispan(315,359,1)
lon_WE(45:75) = ispan(0,30,1)
delete([/topo0,lat,lon/])
topo0 = f00->topo({20:70},{lon_WE})
lat = topo0&lat
lon = topo0&lon



f0 = addfile("data/vwind.2004-01.daily.nc","r")
v0 = short2flt(f0->v(0,{500},{20:70},{lon_WE}))
dim_0 = dimsizes(v0)
delete(v0)
v = new((/44,92,dim_0(0),dim_0(1)/), float)

do i = 0,43
print(" i = "+i)
f1 = addfile("data/vwind."+(i+1979)+"-06.daily.nc","r")
v(i,0:29,:,:)  = short2flt(f1->v(:,{500},{20:70},{lon_WE}))
f2 = addfile("data/vwind."+(i+1979)+"-07.daily.nc","r")
v(i,30:60,:,:) = short2flt(f2->v(:,{500},{20:70},{lon_WE}))
f3 = addfile("data/vwind."+(i+1979)+"-08.daily.nc","r")
v(i,61:91,:,:) = short2flt(f3->v(:,{500},{20:70},{lon_WE}))
end do
v!0 = "year"
v!1 = "day"

v_era_JJA = dim_avg_n_Wrap(v, 1)
v_era_JJA = dim_standardize_n_Wrap(v_era_JJA, 0, 0)
v_era_JJA!0 = "time"
delete(v)

;;----------------  EOF  ----------------;;
w = cos(0.01745329*v_era_JJA&latitude)
wp= v_era_JJA*conform(v_era_JJA, w, 1)
copy_VarCoords(v_era_JJA, wp)

x = wp(latitude|:,longitude|:,time|:)
eof0 = eofunc_Wrap(x, 1, False) 
eof_ts_yearly0 = eofunc_ts_Wrap(x, eof0, False)
pcWest_yearly = dim_standardize_n_Wrap(eof_ts_yearly0(0,:), 1, 0)
pcWest_dtrend = dtrend_n(pcWest_yearly, False, 0)
pcWest_trend = regCoef_n(fspan(1,44,44), pcWest_yearly, 0, 0)
pcWest_cor = escorc_n(fspan(1,44,44), pcWest_yearly, 0, 0)
pcWest_t = rtest(pcWest_cor, 44, 0)
delete([/w,wp,x,eof0,v_era_JJA/])
print("West "+pcWest_trend+" "+pcWest_cor+" "+pcWest_t)




;;----------------------------------  wave train (EE)  ----------------------------------------;;
f0 = addfile("data/vwind.2004-01.daily.nc","r")
v0 = short2flt(f0->v(0,{500},{20:70},{20:90}))
dim_0 = dimsizes(v0)
delete(v0)
v = new((/44,92,dim_0(0),dim_0(1)/), float)

do i = 0,43
print(" i = "+i)
f1 = addfile("data/vwind."+(i+1979)+"-06.daily.nc","r")
v(i,0:29,:,:)  = short2flt(f1->v(:,{500},{20:70},{20:90}))
f2 = addfile("data/vwind."+(i+1979)+"-07.daily.nc","r")
v(i,30:60,:,:) = short2flt(f2->v(:,{500},{20:70},{20:90}))
f3 = addfile("data/vwind."+(i+1979)+"-08.daily.nc","r")
v(i,61:91,:,:) = short2flt(f3->v(:,{500},{20:70},{20:90}))
end do
v!0 = "year"
v!1 = "day"

v_era_JJA = dim_avg_n_Wrap(v, 1)
v_era_JJA = dim_standardize_n_Wrap(v_era_JJA, 0, 0)
v_era_JJA!0 = "time"
delete(v)

;;----------------  EOF  ----------------;;
w = cos(0.01745329*v_era_JJA&latitude)
wp= v_era_JJA*conform(v_era_JJA, w, 1)
copy_VarCoords(v_era_JJA, wp)

x = wp(latitude|:,longitude|:,time|:)
eof0 = eofunc_Wrap(x, 1, False) 
eof_ts_yearly0 = eofunc_ts_Wrap(x, eof0, False)
pcEast_yearly = dim_standardize_n_Wrap(eof_ts_yearly0(0,:), 1, 0)
pcEast_dtrend = dtrend_n(pcEast_yearly, False, 0)
pcEast_trend = regCoef_n(fspan(1,44,44), pcEast_yearly, 0, 0)
pcEast_cor = escorc_n(fspan(1,44,44), pcEast_yearly, 0, 0)
pcEast_t = rtest(pcEast_cor, 44, 0)
print("West "+pcEast_trend+" "+pcEast_cor+" "+pcEast_t)






;;----------------------------------------------------------;;
;;--------------    read heatwave intensity    -------------;;
delete([/topo0, lat, lon/])
f00 = addfile("data/topo_adapt_era5_air2m.nc","r")
topo0 = f00->topo({10:85},{0:360})
lat = topo0&lat
lon = topo0&lon

f3 = addfile("data/era5_Tmax_heatwave_intensity_yearly.nc","r")
hwi_yearly = doubletofloat(f3->heatwave_intensity(:,10:85,:))
hwi_yearly@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(hwi_yearly))) then
      replace_ieeenan (hwi_yearly, value, 0)
      hwi_yearly@_FillValue = value
end if

hwi_yearly!1 = "lat"
hwi_yearly!2 = "lon"
hwi_yearly&lat = lat
hwi_yearly&lon = lon
hwi_yearly2 = hwi_yearly

;;-------------- remove ocean, remove trend -------------;;
topo2 = conform(hwi_yearly, topo0, (/1,2/))
hwi_yearly = where(topo2.ge.0.0, hwi_yearly, 1.0e20)
hwi_yearly@_FillValue = 1.0e20
copy_VarCoords(hwi_yearly2, hwi_yearly)
hwi_detrend = dtrend_msg_n(fspan(1, 44, 44), hwi_yearly, False, False, 0)
copy_VarCoords(hwi_yearly, hwi_detrend)
print(" --- hwi --- ")
printMinMax(hwi_yearly, 1)
printMinMax(hwi_detrend, 1)


;;---------------------------------------------------------------------;;
;;--------------------------    regressions    ------------------------;;
f0 = addfile("data/vwind.2004-01.daily.nc","r")
v0 = short2flt(f0->v(0,{500},{10:85},{0:360}))
dim_1 = dimsizes(v0)
delete(v0)
h = new((/44,92,dim_1(0),dim_1(1)/), float)
do i = 0,43
print(" i = "+i)
f1 = addfile("data/zg."+(i+1979)+"-06.daily.nc","r")
h(i,0:29,:,:)  = short2flt(f1->z(:,{200},{10:85},{0:360}))
f2 = addfile("data/zg."+(i+1979)+"-07.daily.nc","r")
h(i,30:60,:,:) = short2flt(f2->z(:,{200},{10:85},{0:360}))
f3 = addfile("data/zg."+(i+1979)+"-08.daily.nc","r")
h(i,61:91,:,:) = short2flt(f3->z(:,{200},{10:85},{0:360}))
end do
h_JJA = dim_avg_n_Wrap(h, 1)
delete(h)


u = new((/44,92,dim_1(0),dim_1(1)/), float)
do i = 0,43
print(" i = "+i)
f1 = addfile("data/uwind."+(i+1979)+"-06.daily.nc","r")
u(i,0:29,:,:)  = short2flt(f1->u(:,{200},{10:85},{0:360}))
f2 = addfile("data/uwind."+(i+1979)+"-07.daily.nc","r")
u(i,30:60,:,:) = short2flt(f2->u(:,{200},{10:85},{0:360}))
f3 = addfile("data/uwind."+(i+1979)+"-08.daily.nc","r")
u(i,61:91,:,:) = short2flt(f3->u(:,{200},{10:85},{0:360}))
end do
u_JJA = dim_avg_n_Wrap(u, 1)
delete(u)


v = new((/44,92,dim_1(0),dim_1(1)/), float)
do i = 0,43
print(" i = "+i)
f1 = addfile("data/vwind."+(i+1979)+"-06.daily.nc","r")
v(i,0:29,:,:)  = short2flt(f1->v(:,{200},{10:85},{0:360}))
f2 = addfile("data/vwind."+(i+1979)+"-07.daily.nc","r")
v(i,30:60,:,:) = short2flt(f2->v(:,{200},{10:85},{0:360}))
f3 = addfile("data/vwind."+(i+1979)+"-08.daily.nc","r")
v(i,61:91,:,:) = short2flt(f3->v(:,{200},{10:85},{0:360}))
end do
v_JJA = dim_avg_n_Wrap(v, 1)
delete(v)




h_cor1 = escorc_n(pcWest_dtrend, h_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), h_cor1)
u_cor1 = escorc_n(pcWest_dtrend, u_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), u_cor1)
v_cor1 = escorc_n(pcWest_dtrend, v_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), v_cor1)
hwi_cor1 = escorc_n(pcWest_dtrend, hwi_detrend, 0, 0)
copy_VarCoords(h_JJA(0,:,:), hwi_cor1)

h_reg1 = regCoef_n(pcWest_dtrend, h_JJA, 0, 0) / 9.8
copy_VarCoords(h_JJA(0,:,:), h_reg1)
u_reg1 = regCoef_n(pcWest_dtrend, u_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), u_reg1)
v_reg1 = regCoef_n(pcWest_dtrend, v_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), v_reg1)
hwi_reg1 = regCoef_n(pcWest_dtrend, hwi_detrend, 0, 0)
copy_VarCoords(h_JJA(0,:,:), hwi_reg1)

h_reg1 = mask(h_reg1, h_cor1.ge.0.2974.or.h_cor1.le.-0.2974, True)
copy_VarCoords(h_JJA(0,:,:), h_reg1)
u_reg1 = mask(u_reg1, u_cor1.ge.0.2974.or.u_cor1.le.-0.2974.or.v_cor1.ge.0.2974.or.v_cor1.le.-0.2974, True)
copy_VarCoords(h_JJA(0,:,:), u_reg1)
v_reg1 = mask(v_reg1, u_cor1.ge.0.2974.or.u_cor1.le.-0.2974.or.v_cor1.ge.0.2974.or.v_cor1.le.-0.2974, True)
copy_VarCoords(h_JJA(0,:,:), v_reg1)




h_cor2 = escorc_n(pcEast_dtrend, h_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), h_cor2)
u_cor2 = escorc_n(pcEast_dtrend, u_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), u_cor2)
v_cor2 = escorc_n(pcEast_dtrend, v_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), v_cor2)
hwi_cor2 = escorc_n(pcEast_dtrend, hwi_detrend, 0, 0)
copy_VarCoords(h_JJA(0,:,:), hwi_cor2)

h_reg2 = regCoef_n(pcEast_dtrend, h_JJA, 0, 0) / 9.8
copy_VarCoords(h_JJA(0,:,:), h_reg2)
u_reg2 = regCoef_n(pcEast_dtrend, u_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), u_reg2)
v_reg2 = regCoef_n(pcEast_dtrend, v_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), v_reg2)
hwi_reg2 = regCoef_n(pcEast_dtrend, hwi_detrend, 0, 0)
copy_VarCoords(h_JJA(0,:,:), hwi_reg2)

h_reg2 = mask(h_reg2, h_cor2.ge.0.2974.or.h_cor2.le.-0.2974, True)
copy_VarCoords(h_JJA(0,:,:), h_reg2)
u_reg2 = mask(u_reg2, u_cor2.ge.0.2974.or.u_cor2.le.-0.2974.or.v_cor2.ge.0.2974.or.v_cor2.le.-0.2974, True)
copy_VarCoords(h_JJA(0,:,:), u_reg2)
v_reg2 = mask(v_reg2, u_cor2.ge.0.2974.or.u_cor2.le.-0.2974.or.v_cor2.ge.0.2974.or.v_cor2.le.-0.2974, True)
copy_VarCoords(h_JJA(0,:,:), v_reg2)

















;************************************************** plots **********************************************
res = True
res@tmBorderThicknessF = 1.2
res@gsnDraw      =  False
res@gsnFrame     =  False
res@gsnAddCyclic =  True
res@gsnRightString       = ""
res@gsnLeftString        = ""
res@gsnLeftStringFontHeightF   = 0.04
res@gsnRightStringFontHeightF  = 0.04
res@gsnCenterStringFontHeightF  = 0.018
res@tmXTLabelsOn  = False
res@tmYRLabelsOn  = False
vcres = res
res@tmXTOn        = False
res@tmYROn        = False
res@tmXBOn        = True
res@tmYLOn        = True
res@tmXBLabelFontHeightF = 0.017
res@tmYLLabelFontHeightF = 0.017
res@tmXBTickSpacingF = 45.0
res@tmYLTickSpacingF = 20.0

res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/) + 165+3
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
res@tmXBMinorOn = True
res@tmXBMinorValues = fspan(0, 360, 25)+3
; res@tmYLMode = "Explicit"
; res@tmYLValues = (/20,30,40/)
; res@tmYLLabels = (/"20N","30N","40N"/)

res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.013
res@tmXBMajorLengthF = 0.013
res@tmYLMinorLengthF = 0.01
res@tmXBMinorLengthF = 0.01
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.013
res@tmYLMajorOutwardLengthF = 0.013
res@tmXBMinorOutwardLengthF = 0.01
res@tmYLMinorOutwardLengthF = 0.01




res@mpMinLatF = 15
res@mpMaxLatF = 85
res@mpMinLonF = 0
res@mpMaxLonF = 360
res@mpCenterLonF = 180-165-3
res@mpFillOn                    = False
res@mpOutlineOn                 = True
res@mpGeophysicalLineThicknessF = 1.0
res@mpGeophysicalLineColor      = "gray20";"gray20"
;res@mpGridAndLimbOn = True
;res@mpGridLatSpacingF = 90
;res@mpGridLonSpacingF = 360
;res@gsnMaximize = True
res@cnFillOn             = True
res@cnLinesOn            = False
;res@cnLineThicknessF     = 6.0
;res@cnLineColor          = "red"
res@cnLineLabelsOn       = False
res@lbLabelBarOn         = False
res@lbOrientation        = "Vertical"
res@pmLabelBarWidthF     = 0.045
res@pmLabelBarHeightF    = 0.14
res@pmLabelBarOrthogonalPosF = 0.012
res@lbLabelFontHeightF  = 0.015
;res@cnMissingValFillColor = "white"
;res@gsnYRefLine = 0.0
;res@gsnYRefLineColor = "black"
;res@gsnYRefLineThicknessF = 5.0
res@cnLevelSelectionMode = "ExplicitLevels"              
; res@cnLevels             = (/140.0,155.0,170.0,185.0,200.0,215.0,230.0,245.0/)
; res@cnFillColors         = (/16,15,13,12,8,6,5,3,-1/)
res@cnFillPalette        = "precip2_17lev" 



;vcres@mpFillDrawOrder         = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.255
;vcres@vcRefAnnoSide           = "TopRight"
vcres@vcGlyphStyle            = "LineArrow"
vcres@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
vcres@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref
vcres@vcMinDistanceF          = 0.025             ; thin out vectors
vcres@vcLineArrowColor        = "black"           ; change vector color
vcres@vcRefAnnoOn             = True
; vcres@vcRefAnnoString1        = "4"
; vcres@vcRefMagnitudeF         = 4
vcres@vcRefLengthF            = 0.025            ; ref vec length
vcres@vcRefAnnoFontHeightF    = 0.015   ;参考箭头字体大小
vcres@vcRefAnnoString1On      = True
vcres@vcRefAnnoString2On      = False
vcres@vcLineArrowThicknessF   = 1.5            ; make vectors larger
vcres@vcVectorDrawOrder = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -1.0 ;;正向下
vcres@vcRefAnnoParallelPosF   = 1.0  ;;正向右



sres = True
sres@cnLineLabelsOn   = False
sres@cnConstFLabelOn  = False
sres@cnInfoLabelOn    = False 
sres@gsnDraw          = False                   
sres@gsnFrame         = False
sres@gsnLeftString    = ""
sres@gsnRightString   = ""
sres@lbLabelBarOn = False
sres@cnFillOn  = False
sres@cnLinesOn = True
sres2 = sres
sres3 = sres
sres@cnLineColor = "green4"
sres@cnLevelSelectionMode = "ExplicitLevels"
sres@cnLevels         = (/0.5/)
sres@cnFillColors     = (/"black","gray20"/)
sres@cnMissingValFillColor= -1
sres@cnLineThicknessF = 2.0


sres2@cnLevelSelectionMode = "ExplicitLevels"
sres2@cnLevels         = (/2/)
sres2@cnLineColor = "blue"
sres2@cnLineDashPattern = 14
sres2@cnLineThicknessF = 4
; sres2@cnLineColors     = (/"red","red","red","red"/)
; sres2@cnLineDashPatterns = (/16,16,16,0,0,0,0/)
; sres2@cnLineThicknesses = (/4,4,4,4,4,4,4/)




lires = True
lires@tmBorderThicknessF = 1.2
lires@gsnDraw        = False
lires@gsnFrame       = False
lires@tiXAxisOn      = False
lires@tiYAxisOn      = False
lires@tmXTLabelsOn   = False
lires@tmXTOn         = False

lires@tmXBLabelFontHeightF = 0.0155
lires@tmYLLabelFontHeightF = 0.0155
lires@tmYRLabelFontHeightF = 0.0155
lires@gsnLeftStringFontHeightF   = 0.022
lires@gsnRightStringFontHeightF  = 0.022
lires@tmYLLabelsOn   = True
lires@tmYLOn         = True



lires@lgPerimOn = False
lires@lgLabelFontHeightF = 0.018
lires@vpWidthF  = 0.5
lires@vpHeightF = 0.15
lires@vpXF      = 0.06
lires@vpYF      = 0.8

lires@gsnLeftString  = ""
lires@gsnRightString = ""
lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)
lires@xyMarkers = (/1,1,1/)
lires@xyMarkerSizes = (/0.012,0.001,0.001/)
lires@xyMarkerColors= (/"red","black","black"/)
lires@xyMarkerThicknesses = (/4.0,2.0,1.0/)
lires@xyLineThicknesses  = (/3.5,8.0,8.0/)
; lires@xyLineColors       = (/"red","blue","orange"/)
lires@xyDashPatterns     = (/0,4,2/)

lires@tmYLMajorThicknessF = 1.6
lires@tmYRMajorThicknessF = 1.6
lires@tmXBMajorThicknessF = 1.5
lires@tmYLMinorThicknessF = 1.0
lires@tmYRMinorThicknessF = 1.0
lires@tmXBMinorThicknessF = 1.0
lires@tmYLMajorLengthF = 0.012
lires@tmYRMajorLengthF = 0.012
lires@tmXBMajorLengthF = 0.012
lires@tmYLMinorLengthF = 0.009
lires@tmYRMinorLengthF = 0.009
lires@tmXBMinorLengthF = 0.009
lires@tmXBLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.012
lires@tmYLMajorOutwardLengthF = 0.012
lires@tmYRMajorOutwardLengthF = 0.012
lires@tmXBMinorOutwardLengthF = 0.009
lires@tmYLMinorOutwardLengthF = 0.009
lires@tmYRMinorOutwardLengthF = 0.009




;;************************* dotting **************************
ores                 = True            
ores@gsnDraw         = False               
ores@gsnFrame        = False                 
ores@cnLineLabelsOn  = False               
ores@cnLinesOn       = False                 
ores@cnInfoLabelOn   = False                                                                                                        
ores@cnFillOn        = True                               
ores@lbLabelBarOn    = False                                                                                                            
ores@cnLevelSelectionMode = "ExplicitLevels"                                                                                                                                                                                
; ores@cnLevels        = (/14.007/)   ;;90%
; ores@cnLevels        = (/-0.312,0.312/)   ;;95%
ores@cnMonoFillPattern    = False            ; want multiple patterns                                                               
; ores@cnFillPatterns       = (/-1,11/)     ; the patterns                                                                         
ores@cnMonoFillScale      = False            ; want different densities                                                             
ores@cnFillScales         = (/0.6,0.6/)    ; change densities                                                                         
ores@cnMonoFillColor      =True                                                                                                     
ores@cnFillDotSizeF       = 0.003    
ores@cnFillDrawOrder ="postDraw"
ores@cnFillColor = "gray30"
; ores2 = ores
; ores2@cnFillColor = "gray20"
; ores2@cnFillPatterns       = (/6,-1,-1/) 




;**********************plot**************************
plot = new(14, "graphic")
vector = new(14, "graphic")
contour1 = new(14, "graphic")
contour2 = new(14, "graphic")
contour3 = new(14, "graphic")
contour4 = new(14, "graphic")
contour5 = new(14, "graphic")
contour6 = new(14, "graphic")
contour7 = new(14, "graphic")
contour8 = new(14, "graphic")
contour9 = new(14, "graphic")
topooo = new(14, "graphic") 











;;-----------------------------------------------------
pltType = "pdf"
pltName = "Figure4"
wks = gsn_open_wks(pltType, pltName)


res@cnFillPalette = "MPL_PiYG"
res@cnLevels      = (/-250.0,-120,-40,40,120,250.0/)
res@cnFillColors  = (/106,93,80,-1,45,32,19/)
res@tmXBLabelsOn = False
plot(0) = gsn_csm_contour_map_ce(wks, trend_Positive, res)

ores@cnLevels        = (/0.2512,0.2973/)   ;; 90%, 95%                                                           
ores@cnFillPatterns       = (/-1,11,-1/) 
ores@cnFillScales         = (/0.5,0.5/)    ; change densities   
ores@cnFillColor = "red4"

ores@cnLevels        = (/0.2973,1.0/)   ;; 90%, 95%                                                           
ores@cnFillPatterns       = (/-1,11,17/) 
contour4(0) = gsn_csm_contour(wks, cor_High_Posit, ores)    
overlay(plot(0), contour4(0))
contour6(0) = gsn_csm_contour(wks, cor_High_Posit, ores)    
overlay(plot(0), contour6(0))
contour8(0) = gsn_csm_contour(wks, cor_High_Posit, ores)    
overlay(plot(0), contour8(0))


ores@cnLevels        = (/-0.2973,-0.2512/)   ;; 90%, 95%                                                              
ores@cnFillPatterns       = (/-1,11,-1/) 
ores@cnFillColor = "darkgreen"

ores@cnLevels        = (/-1.0,-0.2973/)   ;; 90%, 95%                                                              
ores@cnFillPatterns       = (/17,11,-1/) 
contour5(0) = gsn_csm_contour(wks, cor_High_Negat, ores)    
overlay(plot(0), contour5(0))
contour7(0) = gsn_csm_contour(wks, cor_High_Negat, ores)    
overlay(plot(0), contour7(0))
contour9(0) = gsn_csm_contour(wks, cor_High_Negat, ores)    
overlay(plot(0), contour9(0))




res@cnFillPalette = "MPL_PiYG"
res@cnLevels      = (/-250.0,-120,-40,40,120,250.0/)
res@cnFillColors  = (/106,93,80,-1,45,32,19/)
res@tmXBLabelsOn = True
plot(1) = gsn_csm_contour_map_ce(wks, trend_Negative, res)


ores@cnFillColor = "red4"
ores@cnLevels        = (/0.2512,1.0/)   ;; 90%, 95%  
; ores@cnLevels        = (/0.2973,1.0/)   ;; 90%, 95%                                                           
ores@cnFillPatterns       = (/-1,11,11/) 
contour4(1) = gsn_csm_contour(wks, cor_Low_Posit, ores)    
overlay(plot(1), contour4(1))
contour6(1) = gsn_csm_contour(wks, cor_Low_Posit, ores)    
overlay(plot(1), contour6(1))
contour8(1) = gsn_csm_contour(wks, cor_Low_Posit, ores)    
overlay(plot(1), contour8(1))



ores@cnFillColor = "darkgreen"
ores@cnLevels        = (/-1.0,-0.2512/)   ;; 90%, 95%      
; ores@cnLevels        = (/-1.0,-0.2973/)   ;; 90%, 95%                                                              
ores@cnFillPatterns       = (/11,11,-1/) 
contour5(1) = gsn_csm_contour(wks, cor_Low_Negat, ores)    
overlay(plot(1), contour5(1))
contour7(1) = gsn_csm_contour(wks, cor_Low_Negat, ores)    
overlay(plot(1), contour7(1))
contour9(1) = gsn_csm_contour(wks, cor_Low_Negat, ores)    
overlay(plot(1), contour9(1))

















res@mpMinLatF = 15
res@mpMaxLatF = 88
res@mpMinLonF = -55
res@mpMaxLonF = 75
res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/)
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
delete(res@tmXBMinorValues)
res@tmXBMinorValues = fspan(-60, 120, 13)
res@tmXBLabelFontHeightF = 0.0385
res@tmYLLabelFontHeightF = 0.0385

res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.03
res@tmXBMajorLengthF = 0.03
res@tmYLMinorLengthF = 0.02
res@tmXBMinorLengthF = 0.02
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.03
res@tmYLMajorOutwardLengthF = 0.03
res@tmXBMinorOutwardLengthF = 0.02
res@tmYLMinorOutwardLengthF = 0.02

delete(res@tmXBValues)
delete(res@tmXBLabels)
delete(res@tmXBMinorValues)
res@tmXBMode = "Explicit"
res@tmXBValues = (/-60,-30,0,30,60,90,120,150/)
res@tmXBLabels = (/"60W","30W","0","30E","60E","90E","120E","150E"/)
res@tmXBMinorOn = True
res@tmXBMinorValues = fspan(-360, 360, 49)



res@cnFillPalette = "cmocean_curl" ;"MPL_PiYG"
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnLevels      = (/-12.0,-5.0,-1.0,1.0,5.0,12.0/)
res@cnFillColors  = (/10,50,90,-1,150,190,230/)
res@cnFillMode = "CellFill"

plot(2) = gsn_csm_contour_map_ce(wks, WE_points, res)
res@mpMinLonF = -25
res@mpMaxLonF = 105
res@tmYLLabelsOn = True
plot(3) = gsn_csm_contour_map_ce(wks, EE_points, res)


  txres               = True
  txres@txFont = "helvetica"
  txres@txFontColor = "blue"
  txres@txFontHeightF = 0.032
  eum1 = gsn_add_text(wks,plot(2), "29%", -40,33.0, txres)
  eum2 = gsn_add_text(wks,plot(2), "36%", -5.0,25, txres)
  eum3 = gsn_add_text(wks,plot(3), "77%", 95.0,40, txres)


;;------------------------ figs (a-b) boxes -----------------------
y1=25.0
y2=75
x1=58
x2=85
delete(x)
x:=(/x1,x2, 102,102,x2,    x2,x1,x1/)  ;; central Asia
y:=(/y2,y2, y2, 60, 60,    y1,y1,y2/)
   lnres   =    True
   lnres@gsLineThicknessF=2.8
   lnres@gsnFrame = False

   lnres@gsLineColor="blue"
   lnres@gsLineDashPattern=2
   dum1=gsn_add_polyline(wks,plot(3),x,y,lnres)
   dum2=gsn_add_polyline(wks,plot(1),x,y,lnres)

y1=22.0
y2=38
x1=5
x2=25
x:=(/x1,x2,x2,x1,x1/)  ;; northern Africa
y:=(/y2,y2,y1,y1,y2/)
   dum3=gsn_add_polyline(wks,plot(2),x,y,lnres)
   dum4=gsn_add_polyline(wks,plot(1),x,y,lnres)
y1=40.0
y2=60
x1=-45
x2=0
x=(/x1,x2,x2,x1,x1/)  ;; northeastern Atlantic
y=(/y2,y2,y1,y1,y2/)
   dum5=gsn_add_polyline(wks,plot(2),x,y,lnres)
   dum6=gsn_add_polyline(wks,plot(1),x,y,lnres)

y1=45.0
y2=65
x1=25
x2=55
x=(/x1,x2,x2,x1,x1/)  ;; eastern Europe
y=(/y2,y2,y1,y1,y2/)
   lnres   =    True
   lnres@gsLineThicknessF=2.8
   lnres@gsnFrame = False
   lnres@gsLineColor="red4"
   lnres@gsLineDashPattern=2
   dum7=gsn_add_polyline(wks,plot(0),x,y,lnres)

y1=35.0
y2=55
x1=-10
x2=25
x=(/x1,x2,x2,x1,x1/)  ;; western Europe
y=(/y2,y2,y1,y1,y2/)
   lnres@gsLineColor="red4"
   lnres@gsLineDashPattern=2
   dum8=gsn_add_polyline(wks,plot(0),x,y,lnres)







;;----------------------- EOF PC regressed hwi & UV500 -----------------------
res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.018
res@tmXBMajorLengthF = 0.018
res@tmYLMinorLengthF = 0.012
res@tmXBMinorLengthF = 0.012
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.018
res@tmYLMajorOutwardLengthF = 0.018
res@tmXBMinorOutwardLengthF = 0.012
res@tmYLMinorOutwardLengthF = 0.012


res@mpMinLatF = 15
res@mpMaxLatF = 80
res@mpMinLonF = -100
res@mpMaxLonF = 115
res@tmXBLabelFontHeightF = 0.0245
res@tmYLLabelFontHeightF = 0.0245
delete(res@tmXBValues)
delete(res@tmXBLabels)
delete(res@tmXBMinorValues)
res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/)
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
res@tmXBMinorOn = True
res@tmXBMinorValues = fspan(-360, 360, 49)
res@cnFillPalette = "BlueWhiteOrangeRed"
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnLevels      = (/-5.0,-1.0,1.0,5.0/)
res@cnFillColors  = (/75,100,-1,155,175/)
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 1.0
plot(4) = gsn_csm_contour_map_ce(wks, hwi_reg1, res)

ores@cnFillColor = "red4"
delete(ores@cnLevels)
delete(ores@cnFillPatterns)
ores@cnLevels        = (/0.2973/)   ;; 95%                                                             
ores@cnFillPatterns       = (/-1,11/) 
contour2(4) = gsn_csm_contour(wks, hwi_cor1, ores)    
overlay(plot(4), contour2(4))
contour3(4) = gsn_csm_contour(wks, hwi_cor1, ores)    
overlay(plot(4), contour3(4))
ores@cnFillColor = "blue"
ores@cnLevels        = (/-0.2973/)   ;; 95%                                                             
ores@cnFillPatterns       = (/11,-1/) 
contour4(4) = gsn_csm_contour(wks, hwi_cor1, ores)    
overlay(plot(4), contour4(4))
contour5(4) = gsn_csm_contour(wks, hwi_cor1, ores)    
overlay(plot(4), contour5(4))

vcres@vcRefAnnoString1        = "1"
vcres@vcRefMagnitudeF         = 1
vcres@vcLineArrowColor        = "darkgreen"
vcres@vcRefAnnoFontColor      = "darkgreen"
vcres@vcRefAnnoArrowFillColor = "darkgreen"
vcres@vcRefAnnoArrowEdgeColor = "darkgreen"
vcres@vcRefAnnoArrowUseVecColor = True
vcres@vcLineArrowColor        = "darkgreen"
vector(4)=gsn_csm_vector(wks, u_reg1, v_reg1, vcres)
overlay(plot(4), vector(4))



res@mpMinLonF = -20
res@mpMaxLonF = 195
delete(res@mpCenterLonF)
res@mpCenterLonF = 175/2.0
res@cnFillPalette = "BlueWhiteOrangeRed"
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnLevels      = (/-5.0,-1.0,1.0,5.0/)
res@cnFillColors  = (/75,100,-1,155,175/)
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 1.0
res@tmYLLabelsOn = False
plot(5) = gsn_csm_contour_map_ce(wks, hwi_reg2, res)

ores@cnFillColor = "red4"
ores@cnLevels        = (/0.2973/)   ;; 95%                                                             
ores@cnFillPatterns       = (/-1,11/) 
contour2(5) = gsn_csm_contour(wks, hwi_cor2, ores)    
overlay(plot(5), contour2(5))
contour3(5) = gsn_csm_contour(wks, hwi_cor2, ores)    
overlay(plot(5), contour3(5))
ores@cnFillColor = "blue"
ores@cnLevels        = (/-0.2973/)   ;; 95%                                                             
ores@cnFillPatterns       = (/11,-1/) 
contour4(5) = gsn_csm_contour(wks, hwi_cor2, ores)
overlay(plot(5), contour4(5))
contour5(5) = gsn_csm_contour(wks, hwi_cor2, ores)    
overlay(plot(5), contour5(5))

vcres@vcRefAnnoString1        = "1"
vcres@vcRefMagnitudeF         = 1
vcres@vcLineArrowColor        = "darkgreen"
vcres@vcRefAnnoFontColor      = "darkgreen"
vcres@vcRefAnnoArrowFillColor = "darkgreen"
vcres@vcRefAnnoArrowEdgeColor = "darkgreen"
vcres@vcRefAnnoArrowUseVecColor = True
vcres@vcLineArrowColor        = "darkgreen"
vector(5)=gsn_csm_vector(wks, u_reg2, v_reg2, vcres)
overlay(plot(5), vector(5))

print(" West  trend "+ pcWest_trend + " "+pcWest_t)
print(" East  trend "+ pcEast_trend + " "+pcEast_t)


;;------------------------ figs (e-f) EOF boxes -----------------------
y1=20.0
y2=70
x1=-45
x2=30
delete(x)
x=(/x1,x2,x2,x1,x1/)  ;; western Europe-extended
y=(/y2,y2,y1,y1,y2/)
   lnres   =    True
   lnres@gsLineThicknessF=2.8
   lnres@gsnFrame = False
   lnres@gsLineColor="red4"
   lnres@gsLineDashPattern=2
   dumeof1=gsn_add_polyline(wks,plot(4),x,y,lnres)
y1=20.0
y2=70
x1=20
x2=90
x=(/x1,x2,x2,x1,x1/)  ;; eastern Europe-extended
y=(/y2,y2,y1,y1,y2/)
   dumeof2=gsn_add_polyline(wks,plot(5),x,y,lnres)





lires@trYMinF = -3.0
lires@trYMaxF = 3.0
lires@trXMinF = 1979
lires@trXMaxF = 2022
lires@tmXBMode = "Explicit"
lires@tmXBValues = fspan(1980, 2020, 5)
lires@tmXBLabels = fspan(1980, 2020, 5)
lires@tmXBMinorValues = fspan(1980, 2020, 9)
lires@tmYLMode = "Explicit"
lires@tmYLValues = fspan(-3.0, 3.0, 7)
lires@tmYLLabels = fspan(-3.0, 3.0, 7)
lires@tmYLMinorValues = fspan(-3.0, 3.0, 13)
lires@xyMarkLineModes = (/"MarkLines","Lines","Lines"/)
lires@xyLineColors       = (/"orange","blue","orange"/)
lires@xyLineThicknesses  = (/2.0,8.0,8.0/)
lires@xyDashPatterns  = (/14,14,14/)
plot(6) = gsn_csm_xy(wks, fspan(1979,2022,44), pcWest_yearly, lires)
lires@xyLineColors       = (/"gray30","blue","orange"/)
lires@xyLineThicknesses  = (/1.0,8.0,8.0/)
lires@xyDashPatterns  = (/0,14,14/)
lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)
plot6 = gsn_csm_xy(wks, (/1978,2023/), (/0,0/), lires)
overlay(plot(6), plot6)

lires@xyLineColors       = (/"red4","blue","orange"/)
lires@xyLineThicknesses  = (/2.5,8.0,8.0/)
lires@xyDashPatterns  = (/0,14,14/)
plot8 = gsn_csm_xy(wks, (/2000.5-19,2000.5+19/), (/-19*pcWest_trend, 19*pcWest_trend/), lires)
overlay(plot(6), plot8)






lires@trYMinF = -3.0
lires@trYMaxF = 3.0
lires@trXMinF = 1979
lires@trXMaxF = 2022
lires@tmXBMode = "Explicit"
lires@tmXBValues = fspan(1980, 2020, 5)
lires@tmXBLabels = fspan(1980, 2020, 5)
lires@tmXBMinorValues = fspan(1980, 2020, 9)
lires@tmYLLabelsOn = False
lires@tmYRLabelsOn = True
lires@tmYRMode = "Explicit"
lires@tmYRValues = fspan(-3.0, 3.0, 7)
lires@tmYRLabels = fspan(-3.0, 3.0, 7)
lires@tmYRMinorValues = fspan(-3.0, 3.0, 13)
lires@xyMarkLineModes = (/"MarkLines","Lines","Lines"/)
lires@xyLineColors       = (/"orange","blue","orange"/)
lires@xyLineThicknesses  = (/2.0,8.0,8.0/)
lires@xyDashPatterns  = (/14,14,14/)
plot(7) = gsn_csm_xy(wks, fspan(1979,2022,44), pcEast_yearly, lires)
lires@xyLineColors       = (/"gray30","blue","orange"/)
lires@xyLineThicknesses  = (/1.0,8.0,8.0/)
lires@xyDashPatterns  = (/0,14,14/)
lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)
plot7 = gsn_csm_xy(wks, (/1978,2023/), (/0,0/), lires)
overlay(plot(7), plot7)

lires@xyLineColors       = (/"red4","blue","orange"/)
lires@xyLineThicknesses  = (/2.5,8.0,8.0/)
lires@xyDashPatterns  = (/0,14,14/)
plot9 = gsn_csm_xy(wks, (/2000.5-19,2000.5+19/), (/-19*pcEast_trend, 19*pcWest_trend/), lires)
overlay(plot(7), plot9)


;;------------------ denote high, low, slope ------------------;;
  txres               = True
  txres@txFont = "helvetica"
  txres@txFontColor = "black"
  txres@txFontHeightF = 0.016
  eur5 = gsn_add_text(wks,plot(0), "High", 162,24, txres)
  eur6 = gsn_add_text(wks,plot(1), "Low", 162,24, txres)
  txres@txFontHeightF = 0.014
  eur3 = gsn_add_text(wks,plot(6), "slope : 0.024  (p<0.05)", 1998,-2.3, txres)
  eur4 = gsn_add_text(wks,plot(7), "slope : 0.043  (p<0.001)", 2004,-2.3, txres)







  pres = True
  pres@gsnFrame = False

  pres@gsnPanelLabelBar = False
  pres@gsnPanelTop  = 0.865
  pres@gsnPanelBottom = 0.725
  pres@gsnPanelLeft  = 0.075
  pres@gsnPanelRight = 0.65
  pres@amJust  ="BottomRight"
  pres@gsnPanelFigureStringsFontHeightF = 0.012
  gsn_panel(wks,plot(0),(/1,1/),pres)
  
  pres@gsnPanelLabelBar = True
  pres@lbOrientation        = "Horizontal"
  pres@pmLabelBarWidthF     = 0.35
  pres@pmLabelBarHeightF    = 0.045
  pres@pmLabelBarOrthogonalPosF = -0.0
  ; pres@pmLabelBarParallelPosF = 0.07
  pres@lbLabelFontHeightF  = 0.011
  pres@gsnPanelTop  = 0.76-0.01
  pres@gsnPanelBottom = 0.55-0.01
  pres@gsnPanelLeft  = 0.075
  pres@gsnPanelRight = 0.65
  gsn_panel(wks,plot(1),(/1,1/),pres)


  pres@gsnPanelLabelBar = False
  pres@gsnPanelTop  = 0.86-0.002
  pres@gsnPanelBottom = 0.71-0.002
  pres@gsnPanelLeft  = 0.665
  pres@gsnPanelRight = 0.90
  pres@amJust  ="TopLeft"
  gsn_panel(wks,plot(2),(/1,1/),pres)
  pres@gsnPanelLabelBar = True
  pres@lbOrientation        = "Vertical"
  pres@pmLabelBarWidthF     = 0.04
  pres@pmLabelBarHeightF    = 0.18
  pres@pmLabelBarOrthogonalPosF = 0.003
  pres@pmLabelBarParallelPosF = 0.08
  pres@lbLabelFontHeightF  = 0.0113
  pres@gsnPanelTop  = 0.72-0.002
  pres@gsnPanelBottom = 0.57-0.002
  pres@gsnPanelLeft  = 0.665
  pres@gsnPanelRight = 0.95
  gsn_panel(wks,plot(3),(/1,1/),pres)



  pres@gsnPanelLabelBar = False
  pres@gsnPanelTop  = 0.533+0.005
  pres@gsnPanelBottom = 0.333+0.005
  pres@gsnPanelLeft  = 0.075
  pres@gsnPanelRight = 0.49
  pres@amJust  ="BottomRight"
  gsn_panel(wks,plot(4),(/1,1/),pres)
  pres@gsnPanelLabelBar = True
  pres@lbOrientation        = "Vertical"
  pres@pmLabelBarWidthF     = 0.04
  pres@pmLabelBarHeightF    = 0.12
  pres@pmLabelBarOrthogonalPosF = 0.002
  pres@pmLabelBarParallelPosF = 0.012
  pres@lbLabelFontHeightF  = 0.0115
  pres@gsnPanelTop  = 0.53+0.005
  pres@gsnPanelBottom = 0.33+0.005
  pres@gsnPanelLeft  = 0.51+0.006
  pres@gsnPanelRight = 0.977+0.006
  gsn_panel(wks,plot(5),(/1,1/),pres)


  pres@gsnPanelLabelBar = False
  pres@gsnPanelTop  = 0.375
  pres@gsnPanelBottom = 0.175
  pres@gsnPanelLeft  = 0.088
  pres@gsnPanelRight = 0.498
  pres@amJust  ="BottomRight"
  gsn_panel(wks,plot(6),(/1,1/),pres)
  pres@gsnPanelTop  = 0.375
  pres@gsnPanelBottom = 0.175
  pres@gsnPanelLeft  = 0.511+0.006
  pres@gsnPanelRight = 0.925+0.006
  gsn_panel(wks,plot(7),(/1,1/),pres)



  txres               = True
  txres@txFontHeightF = 0.0145
  txres@txFont = "helvetica"
  gsn_text_ndc(wks, "Zonally  asymmetric  trends  in  pressure  extremes",  0.380, 0.880, txres)
  gsn_text_ndc(wks, "Most-connected  points",  0.80, 0.880, txres)
  gsn_text_ndc(wks, "Western  Europe",  0.298, 0.525, txres)
  gsn_text_ndc(wks, "Eastern  Europe",  0.715, 0.525, txres)


;;------------------------- units ------------------------;; 
txres@txFontHeightF = 0.0105
txres@txFont = "helvetica"
txres@txFontColor = "black"
gsn_text_ndc(wks, "gpm/decade",  0.59, 0.563, txres)  
gsn_text_ndc(wks, "K/decade",  0.92, 0.377, txres)  


txres@txFontColor = "black"
txres@txFontHeightF = 0.0155
txres@txAngleF = 0.0
txres@txFont = "helvetica-bold"
gsn_text_ndc(wks, "a",  0.13, 0.86, txres)
gsn_text_ndc(wks, "b",  0.13, 0.732, txres)
gsn_text_ndc(wks, "c",  0.71, 0.86, txres)
gsn_text_ndc(wks, "d",  0.71, 0.721, txres)

gsn_text_ndc(wks, "e",  0.125, 0.52, txres)
gsn_text_ndc(wks, "f",  0.533, 0.52, txres)
gsn_text_ndc(wks, "g",  0.125, 0.355, txres)
gsn_text_ndc(wks, "h",  0.533, 0.355, txres)


frame(wks)
   delete(wks)
  ;  system("convert -geometry 1800x1800 -density 1200x1200 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
;   system("rm " + pltName + "." + pltType + " -f") 






end






