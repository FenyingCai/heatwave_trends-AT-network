load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin
;;********************************** 读取经纬度 ************************************
f00 = addfile("/public/home/fcai/data0/topography/topo_adapt_era5_air2m.nc","r")
topo0 = f00->topo({10:85},{0:360})
lat = topo0&lat
lon = topo0&lon
dim_0 = dimsizes(topo0)



;;----------------------------------------------------------------;;
;;--------------    读取热浪逐年强度、高压逐年强度    -------------;;
f1 = addfile("/public/home/fcai/extreme1_AT/NC2_1980_2010/1_threshold/era5_Tmax_heatwave_intensity_yearly.nc","r")
hwi_yearly = doubletofloat(f1->heatwave_intensity(:,10:85,:))
hwi_yearly@_FillValue = 1.0e20

f2 = addfile("/public/home/fcai/extreme1_AT/NC2_1980_2010/1_threshold/hgt_ano_44years_JJAyearly.nc","r")
h500Positive_yearly = f2->hgt_Positive_intensity(:,{500},{10:85},:)
h500Positive_yearly@_FillValue = 1.0e20

f3 = addfile("/public/home/fcai/extreme1_AT/NC2_1980_2010/1_threshold/hgt_ano_44years_JJAyearly.nc","r")
h500Negative_yearly = f3->hgt_Negative_intensity(:,{500},{10:85},:)
h500Negative_yearly@_FillValue = 1.0e20





value = 1.0e20
if (any(isnan_ieee(hwi_yearly))) then
      replace_ieeenan (hwi_yearly, value, 0)
      hwi_yearly@_FillValue = value
end if
hwi_yearly!1 = "lat"
hwi_yearly!2 = "lon"
hwi_yearly&lat = lat
hwi_yearly&lon = lon



;;-------------- 去掉海洋区域、不去趋势、标准化 -------------;;
topo2 = conform(hwi_yearly, topo0, (/1,2/))
hwi_yearly = where(topo2.ge.0.0, hwi_yearly, 1.0e20)
hwi_yearly@_FillValue = 1.0e20

; hwi  = dtrend_msg_n(fspan(1, 44, 44), hwi_yearly, False, False, 0)
hwi  = dim_standardize_n_Wrap(hwi_yearly , 1, 0)


;;----------------------- remove zonal mean  ----------------------;;
; h500Positive_zonal = dim_avg_n_Wrap(h500Positive_yearly, 2)
; h500Positive_zonal_trend = regCoef_n(fspan(1,44,44), h500Positive_zonal, 0,0)
; h500Positive_zonal_trend_44years = new((/44,dim_0(0),dim_0(1)/), float)
; do iyear =0,43
; h500Positive_zonal_trend_44years(iyear,:,:) = conform(h500Positive_yearly(0,:,:), (h500Positive_zonal_trend*iyear), 0)
; end do
; h500Positive_yearly = h500Positive_yearly - h500Positive_zonal_trend_44years
h500Positive_yearly = dim_rmvmean_n_Wrap(h500Positive_yearly, 2)


; h500Negative_zonal = dim_avg_n_Wrap(h500Negative_yearly, 2)
; h500Negative_zonal_trend = regCoef_n(fspan(1,44,44), h500Negative_zonal, 0,0)
; h500Negative_zonal_trend_44years = new((/44,dim_0(0),dim_0(1)/), float)
; do iyear =0,43
; h500Negative_zonal_trend_44years(iyear,:,:) = conform(h500Negative_yearly(0,:,:), (h500Negative_zonal_trend*iyear), 0)
; end do
; h500Negative_yearly = h500Negative_yearly - h500Negative_zonal_trend_44years
h500Negative_yearly = dim_rmvmean_n_Wrap(h500Negative_yearly, 2)




; h500Positive  = dtrend_msg_n(fspan(1, 44, 44), h500Positive_yearly, False, False, 0)  ;; 不去趋势
h500Positive  = dim_standardize_n_Wrap(h500Positive_yearly , 1, 0)

; h500Negative  = dtrend_msg_n(fspan(1, 44, 44), h500Negative_yearly, False, False, 0)
h500Negative  = dim_standardize_n_Wrap(h500Negative_yearly , 1, 0)





;;------------  读取距离数据  ------------;;
f0 = addfile("/public/home/fcai/extreme1_AT/NC2_1980_2010/1_threshold/distance_i_j.nc","r")
distance = f0->distance(10:85,:,10:85,:)


;;-------------- 读取相关系数网络(暖-高压) -------------;;
; f5 = addfile("/home/ys17-23/cai_fy/data4_heatwaves/yearly_correlation/H500/CorNetorks_Tmax90th_H500Positive.nc","r")
; cor_4D_Positive = f5->cor_4D

f5 = addfile("/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/part3_heathigh_daily/Networks_Tmax90_H500Positive_lag0_land_Significant99.nc","r")
cor_4D_Positive = f5->networks0
dim_2 = dimsizes(cor_4D_Positive(0,0,:,:))

;;-------------- 求高压中心 (<1.5kkm) 序列 -------------;;
Teleconnection_High = new((/44,dim_0(0),dim_0(1)/), float)
do ilat = 0,dim_0(0)-1
print(" ilat = "+ilat)
do jlon = 0,dim_0(1)-1
; print(" jlon = "+jlon)

if(.not.ismissing(dim_avg(hwi (:,ilat,jlon))))then

cor_4D_Positive_new = where(distance(ilat,jlon,:,:).le.1500.0, cor_4D_Positive(ilat,jlon,:,:), 0.0)
cor_4D_Positive_new_1D = ndtooned(cor_4D_Positive_new)
indices  = ind_resolve(maxind(cor_4D_Positive_new_1D), dim_2)
Teleconnection_High(:,ilat,jlon) = h500Positive (:,indices(0,0),indices(0,1))

; print("max cor "+cor_4D_Positive_new(indices(0,0),indices(0,1)))
; print("max cor (True) "+max(cor_4D_Positive_new_1D))
delete(([/cor_4D_Positive_new,cor_4D_Positive_new_1D,indices/]))
end if
end do
end do
delete(cor_4D_Positive)





;;-------------- 读取相关系数网络(暖-低压) -------------;;
; f5 = addfile("/home/ys17-23/cai_fy/data4_heatwaves/yearly_correlation/H500/CorNetorks_Tmax90th_H500Negative.nc","r")
; cor_4D_Negative = f5->cor_4D

f5 = addfile("/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/part4_heatlow_daily/Networks_Tmax90_H500Negative_lag0_land_Significant99.nc","r")
cor_4D_Negative = f5->networks0
dim_2 = dimsizes(cor_4D_Negative(0,0,:,:))

;;-------------- 求高压中心 ( 1.5-4 kkm ) 序列 -------------;;
Teleconnection_Low = new((/44,dim_0(0),dim_0(1)/), float)
do ilat = 0,dim_0(0)-1
print(" ilat = "+ilat)
do jlon = 0,dim_0(1)-1
; print(" jlon = "+jlon)

if(.not.ismissing(dim_avg(hwi (:,ilat,jlon))))then

cor_4D_Negative_new = where(distance(ilat,jlon,:,:).gt.1500.0 .and. distance(ilat,jlon,:,:).le.4000.0, cor_4D_Negative(ilat,jlon,:,:), 0.0)
cor_4D_Negative_new_1D = ndtooned(cor_4D_Negative_new)
indices  = ind_resolve(maxind(cor_4D_Negative_new_1D), dim_2)
Teleconnection_Low(:,ilat,jlon) = h500Negative (:,indices(0,0),indices(0,1))

; print("max cor "+cor_4D_Negative_new(indices(0,0),indices(0,1)))
; print("max cor (True) "+max(cor_4D_Negative_new_1D))
delete(([/cor_4D_Negative_new,cor_4D_Negative_new_1D,indices/]))
end if
end do
end do
delete(cor_4D_Negative)
delete(distance)



;;-------------- 遥相关序列 (44yrs × lat76 × lon360) -------------;;
Teleconnection = Teleconnection_High + Teleconnection_Low
copy_VarCoords(hwi_yearly, Teleconnection)
copy_VarCoords(hwi_yearly, Teleconnection_High)
copy_VarCoords(hwi_yearly, Teleconnection_Low)

printVarSummary(Teleconnection)





;;------------  书写 degree 文件  ------------;;
output="/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/"
name1="Teleconnection_Tmax90th_H500_removeZonalTrend_ConcurrentDays"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->Teleconnection=Teleconnection
ncdf1->Teleconnection_High=Teleconnection_High
ncdf1->Teleconnection_Low=Teleconnection_Low




end






