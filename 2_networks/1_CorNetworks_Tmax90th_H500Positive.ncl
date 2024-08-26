load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin
;;********************************** 读取经纬度 ************************************
f00 = addfile("/home/ys17-23/Extension/cai_fy/data0/topography/topo_adapt_era5_air2m.nc","r")
topo0 = f00->topo({10:86},{0:360})
lat = topo0&lat
lon = topo0&lon
; dim_0 = dimsizes(topo0)



;;----------------------------------------------------------------;;
;;--------------    读取热浪逐年强度和冷时间逐年强度    -------------;;
f1 = addfile("/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/1_threshold/era5_Tmax_heatwave_intensity_yearly.nc","r")
hwi_yearly = doubletofloat(f1->heatwave_intensity(:,10:86,:))
hwi_yearly@_FillValue = 1.0e20
printMinMax(hwi_yearly, 1)

f2 = addfile("/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/1_threshold/hgt_ano_44years_JJAyearly.nc","r")
h500Positive_yearly = f2->hgt_Positive_intensity(:,{500},{10:85},:)
h500Positive_yearly@_FillValue = 1.0e20
; dim_1 = dimsizes(h500Positive_yearly(0,:,:))


value = 1.0e20
if (any(isnan_ieee(hwi_yearly))) then
      replace_ieeenan (hwi_yearly, value, 0)
      hwi_yearly@_FillValue = value
end if
hwi_yearly!1 = "lat"
hwi_yearly!2 = "lon"
hwi_yearly&lat = lat
hwi_yearly&lon = lon






;;-------------- 去掉海洋区域、去趋势、标准化 -------------;;
topo2 = conform(hwi_yearly, topo0, (/1,2/))
hwi_yearly = where(topo2.ge.0.0, hwi_yearly, 1.0e20)
hwi_yearly@_FillValue = 1.0e20

hwi_detrend = dtrend_msg_n(fspan(1, 44, 44), hwi_yearly, False, False, 0)
hwi_detrend = dim_standardize_n_Wrap(hwi_detrend, 1, 0)

h500Positive_detrend = dtrend_msg_n(fspan(1, 44, 44), h500Positive_yearly, False, False, 0)
h500Positive_detrend = dim_standardize_n_Wrap(h500Positive_detrend, 1, 0)




;;-----------------------------------------------------------;;
;;----------------------- 求相关系数网络  --------------------;;
cor_4D = escorc_n(hwi_detrend, h500Positive_detrend, 0, 0)
printVarSummary(cor_4D)
printMinMax(cor_4D, 0)






;;------------  书写 degree 文件  ------------;;
output="/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/2_networks/"
name1="CorNetorks_Tmax90th_H500Positive"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->cor_4D=cor_4D(0:75,:,:,:)





end






