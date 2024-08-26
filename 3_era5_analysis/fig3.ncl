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

topo22 = f00->topo({0:85},{0:360})
lat2 = topo22&lat
lon2 = topo22&lon
dim_2 = dimsizes(topo22)



;;----------------------------------------------------------;;
;;--------------    read heatwave intensity    -------------;;
f4 = addfile("data/era5_Tmax_heatwave_intensity_yearly.nc","r")
hwi_yearly0 = doubletofloat(f4->heatwave_intensity(:,0:85,:))
hwi_yearly0@_FillValue = 1.0e20
printMinMax(hwi_yearly0,1)

value = 1.0e20
if (any(isnan_ieee(hwi_yearly0))) then
      replace_ieeenan (hwi_yearly0, value, 0)
      hwi_yearly0@_FillValue = value
end if

hwi_yearly0!1 = "lat"
hwi_yearly0!2 = "lon"
hwi_yearly0&lat = lat2
hwi_yearly0&lon = lon2
;;-------------- exclude ocean -------------;;
topo2 = conform(hwi_yearly0, topo22, (/1,2/))
hwi_yearly = where(topo2.ge.0.0, hwi_yearly0, 1.0e20)
hwi_yearly@_FillValue = 1.0e20
copy_VarCoords(hwi_yearly0, hwi_yearly)


time = fspan(1, 44, 44)
hwi_trend = regCoef_n(time, hwi_yearly, 0, 0) * 10.0
copy_VarCoords(hwi_yearly(0,:,:), hwi_trend)
hwi_trend_cor = escorc_n(time, hwi_yearly, 0, 0)
copy_VarCoords(hwi_yearly(0,:,:), hwi_trend_cor)
hwi_trend_t = rtest(hwi_trend_cor, 44, 0)
copy_VarCoords(hwi_yearly(0,:,:), hwi_trend_t)


;;-------------- long-term trend -------------;;
hwi_trend = where(topo22.ge.0.0, hwi_trend, 1.0e20)
hwi_trend_t = where(topo22.ge.0.0, hwi_trend_t, 1.0e20)
hwi_trend_t = where(hwi_trend_cor.ge.0.0, hwi_trend_t, 1.0e20)
hwi_trend@_FillValue = 1.0e20
hwi_trend_t@_FillValue = 1.0e20


;;-------------- remove zonal mean (31 degree moving) -------------;;
hwi_trend_avg = dim_avg_n_Wrap(hwi_trend, 1)
hwi_trend_31moving = new((/dim_2(0),dim_2(1)/), float)


do ilat = 0,85
  lat_min = max((/0,ilat-15/))
  lat_max = min((/85,ilat+15/))
  hwi_trend_31moving(ilat,:) = dim_avg(hwi_trend_avg(lat_min:lat_max))
end do

hwi_trend_ano = hwi_trend(10:85,:) - hwi_trend_31moving(10:85,:)
copy_VarCoords(hwi_yearly(0,10:85,:), hwi_trend_ano)








;;-----------------------------------------------------------;;
;;---------------------- read AT index  ---------------------;;
f1 = addfile("data/Teleconnection_Tmax90th_H500_removeZonalTrend.nc","r")
degree11 = f1->Teleconnection(:,0:,:)  ;; 10-85N
f2 = addfile("data/Teleconnection_Tmax90th_H500_removeZonalTrend_ConcurrentDays.nc","r")
degree12 = f2->Teleconnection(:,0:,:)  ;; 10-85N
degree1 = degree11 + degree12
delete([/degree11,degree12/])
degree1@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(degree1))) then
      replace_ieeenan (degree1, value, 0)
      degree1@_FillValue = value
end if


;;------------- exclude ocean -------------;;
topo3 = conform(degree1, topo0, (/1,2/))
degree1 = where(topo3.ge.0.0, degree1, 1.0e20)
degree1@_FillValue = 1.0e20


degree1!1 = "lat"
degree1!2 = "lon"
degree1&lat = lat
degree1&lon = lon
print(" ----------- degree ---------- ")
printMinMax(degree1, 0)




;;-------------- long-term trend -------------;;
time = fspan(1, 44, 44)
trend1 = regCoef_n(time, degree1, 0, 0) * 10.0
copy_VarCoords(degree1(0,:,:), trend1)
trend1@_FillValue = 1.0e20
print(" ----------- trend ---------- ")
printMinMax(trend1, 0)



;;--------------  significant trends -------------;;
cor1 = escorc_n(time, degree1, 0, 0)
copy_VarCoords(degree1(0,:,:), cor1)
cor1_t = rtest(cor1, 44, 0)
copy_VarCoords(degree1(0,:,:), cor1_t)


cor1_t = where(topo0.ge.0.0, cor1_t, 1.0e20)
copy_VarCoords(degree1(0,:,:), cor1_t)

cor1_t1 = mask(cor1_t, cor1.ge.0.0, True)
cor1_t2 = mask(cor1_t, cor1.le.0.0, True)
copy_VarCoords(degree1(0,:,:), cor1_t1)
copy_VarCoords(degree1(0,:,:), cor1_t2)






;;---------------------------------------------------;;
;;----------------- Estimated trends ----------------;;
time = fspan(1979, 2022, 44)
degree1_dtrend = dtrend_msg_n(time, degree1, False, False, 0)
hwi_dtrend = dtrend_msg_n(time, hwi_yearly(:,10:85,:), False, False, 0)
copy_VarCoords(degree1(:,:,:), degree1_dtrend)
copy_VarCoords(degree1(:,:,:), hwi_dtrend)

do ilat = 0,dim_0(0)-1
      hwi_dtrend(:, ilat, :) = hwi_dtrend(:, ilat, :) * cos((ilat+10) / 180.0 * 3.14159)
end do


cor_hwi_tele = new((/dim_0(0), dim_0(1)/), float)
reg_hwi_tele = new((/dim_0(0), dim_0(1)/), float)
do ilat = 0,dim_0(0)-1
do jlon = 0,dim_0(1)-1
      cor_hwi_tele(ilat, jlon) = escorc_n(hwi_dtrend(:,ilat,jlon), degree1_dtrend(:,ilat,jlon), 0, 0)
      reg_hwi_tele(ilat, jlon) = regCoef_n(degree1_dtrend(:,ilat,jlon), hwi_dtrend(:,ilat,jlon), 0, 0)
end do
end do
copy_VarCoords(degree1(0,:,:), cor_hwi_tele)
copy_VarCoords(degree1(0,:,:), reg_hwi_tele)




;;-------------- explained variances R*R -------------;;
cor_hwi_tele2 = cor_hwi_tele * cor_hwi_tele
cor_hwi_tele2 = mask(cor_hwi_tele2, cor_hwi_tele.ge.0.304, True)
copy_VarCoords(degree1(0,:,:), cor_hwi_tele2)





;;-------------- Estimated trend / Observed -------------;;
Estimated_hwi = reg_hwi_tele * trend1
copy_VarCoords(degree1(0,:,:), Estimated_hwi)

Contribution_hwi = Estimated_hwi / hwi_trend_ano
copy_VarCoords(degree1(0,:,:), Contribution_hwi)

print(" -----------  Estimated / Observed  ---------- ")
Contribution_hwi_all = where(abs(hwi_trend_ano).ge.0.0, 1.0, 0.0)
Contribution_hwi_sam_sign = where(abs(hwi_trend_ano).ge.0.0.and.Contribution_hwi.gt.0.0, 1.0, 0.0)
Contribution_hwi_opposite = where(abs(hwi_trend_ano).ge.0.0.and.Contribution_hwi.lt.0.0, 1.0, 0.0)
do ilat = 0,75
Contribution_hwi_all(ilat,:) = Contribution_hwi_all(ilat,:) * cos((ilat+10) /180.0 * 3.14159)
Contribution_hwi_sam_sign(ilat,:) = Contribution_hwi_sam_sign(ilat,:) * cos((ilat+10) /180.0 * 3.14159)
Contribution_hwi_opposite(ilat,:) = Contribution_hwi_opposite(ilat,:) * cos((ilat+10) /180.0 * 3.14159)
end do
print("all : "+ dim_sum(dim_sum(Contribution_hwi_all(40-10:,:))))
print("same sign : "+ dim_sum(dim_sum(Contribution_hwi_sam_sign(40-10:,:))))
print("opposite sign : "+ dim_sum(dim_sum(Contribution_hwi_opposite(40-10:,:))))
print("same sign / all : "+ dim_sum(dim_sum(Contribution_hwi_sam_sign(40-10:,:))) / dim_sum(dim_sum(Contribution_hwi_all(40-10:,:))))

ratio_lat = dim_sum(Contribution_hwi_sam_sign) / dim_sum(Contribution_hwi_all) * 100.0
print(ratio_lat)

ratio_lat_11moving = new((/76/), float)
do ilat = 0,75
  if (ilat.le.5) then
    ratio_lat_11moving(ilat) = dim_avg(ratio_lat(0:ilat+5))
  end if
  if (ilat.ge.70) then
    ratio_lat_11moving(ilat) = dim_avg(ratio_lat(ilat-5:75))
  end if
  if (ilat.gt.5 .and. ilat.lt.70) then
    ratio_lat_11moving(ilat) = dim_avg(ratio_lat(ilat-5:ilat+5))
  end if
end do
print(ratio_lat_11moving)






Contribution_hwi = where(abs(hwi_trend_ano).ge.0.01, Contribution_hwi, 0.0)  ;; remove margin
copy_VarCoords(degree1(0,:,:), Contribution_hwi)

Contribution_hwi_P = mask(Contribution_hwi, hwi_trend_ano.ge.0.0, True)
Contribution_hwi_N = mask(Contribution_hwi, hwi_trend_ano.le.0.0, True)
copy_VarCoords(degree1(0,:,:), Contribution_hwi_P)
copy_VarCoords(degree1(0,:,:), Contribution_hwi_N)




print(" -----------  Trend  MIN  MAX  ---------- ")
printMinMax(hwi_trend_ano, 0)
printMinMax(Estimated_hwi, 0)
print(" -----------                   ---------- ")







;;---------------------------------------------------------------------------------------------------------;;
;;--------------------------------------------- Seven regions ---------------------------------------------;;
hwi_trend_ano_weighted = new((/dim_0(0),dim_0(1)/), float)
Estimated_hwi_weighted = new((/dim_0(0),dim_0(1)/), float)
grids65 = new((/dim_0(0),dim_0(1)/), float)
printVarSummary(hwi_trend_ano)
printVarSummary(Estimated_hwi)
do ilat = 0,dim_0(0)-1
  hwi_trend_ano_weighted(ilat,:) = hwi_trend_ano(ilat,:) * cos((ilat+10) / 180.0 * 3.14159)
  Estimated_hwi_weighted(ilat,:) = Estimated_hwi(ilat,:) * cos((ilat+10) / 180.0 * 3.14159)
  grids65(ilat,:) = cos((ilat+10) / 180.0 * 3.14159)
end do
grids65 = where(topo0.ge.0.0, grids65, 1.0e20)
grids65@_FillValue = 1.0e20

WE_lon = new(36, integer)
WE_lon(0:9) = ispan(350,359,1)
WE_lon(10:35) = ispan(0,25,1)



hwi_EE = dim_sum(dim_sum(hwi_trend_ano_weighted(45-10:65-10, 25:55))) / dim_sum(dim_sum(grids65(45-10:65-10, 25:55)))
hwi_WE = dim_sum(dim_sum(hwi_trend_ano_weighted(35-10:55-10, WE_lon))) / dim_sum(dim_sum(grids65(35-10:55-10, WE_lon)))
hwi_EA = dim_sum(dim_sum(hwi_trend_ano_weighted(28-10:34-10, 100:120))) / dim_sum(dim_sum(grids65(28-10:34-10, 100:120)))
hwi_SCA = dim_sum(dim_sum(hwi_trend_ano_weighted(22-10:45-10, 65:80))) / dim_sum(dim_sum(grids65(22-10:45-10, 65:80)))
hwi_NWNA = dim_sum(dim_sum(hwi_trend_ano_weighted(60-10:70-10, 200:260))) / dim_sum(dim_sum(grids65(60-10:70-10, 200:260)))
hwi_ENA = dim_sum(dim_sum(hwi_trend_ano_weighted(35-10:50-10, 260:290))) / dim_sum(dim_sum(grids65(35-10:50-10, 260:290)))



estimate_EE = dim_sum(dim_sum(Estimated_hwi_weighted(45-10:65-10, 25:55))) / dim_sum(dim_sum(grids65(45-10:65-10, 25:55)))
estimate_WE = dim_sum(dim_sum(Estimated_hwi_weighted(35-10:55-10, WE_lon))) / dim_sum(dim_sum(grids65(35-10:55-10, WE_lon)))
estimate_EA = dim_sum(dim_sum(Estimated_hwi_weighted(28-10:34-10, 100:120))) / dim_sum(dim_sum(grids65(28-10:34-10, 100:120)))
estimate_SCA = dim_sum(dim_sum(Estimated_hwi_weighted(22-10:45-10, 65:80))) / dim_sum(dim_sum(grids65(22-10:45-10, 65:80)))
estimate_NWNA = dim_sum(dim_sum(Estimated_hwi_weighted(60-10:70-10, 200:260))) / dim_sum(dim_sum(grids65(60-10:70-10, 200:260)))
estimate_ENA = dim_sum(dim_sum(Estimated_hwi_weighted(35-10:50-10, 260:290))) / dim_sum(dim_sum(grids65(35-10:50-10, 260:290)))


hwi_Europe = (dim_sum(dim_sum(hwi_trend_ano_weighted(45-10:65-10, 25:55))) + dim_sum(dim_sum(hwi_trend_ano_weighted(35-10:55-10, WE_lon)))) / \
(dim_sum(dim_sum(grids65(45-10:65-10, 25:55)))+dim_sum(dim_sum(grids65(35-10:55-10, WE_lon))))
estimate_Europe = (dim_sum(dim_sum(Estimated_hwi_weighted(45-10:65-10, 25:55))) + dim_sum(dim_sum(Estimated_hwi_weighted(35-10:55-10, WE_lon)))) /\\
(dim_sum(dim_sum(grids65(45-10:65-10, 25:55)))+ dim_sum(dim_sum(grids65(35-10:55-10, WE_lon))))



hwi_trend_ano_weighted2 = hwi_trend_ano_weighted
Estimated_hwi_weighted2 = Estimated_hwi_weighted
grids2 = grids65
do ilat = 0,dim_0(0)-1
do jlon = 0,dim_0(1)-1
if ((jlon-125.lt.(ilat-(60-10)))) then
  hwi_trend_ano_weighted2(ilat,jlon) = 0.0
  Estimated_hwi_weighted2(ilat,jlon) = 0.0
  grids2(ilat,jlon) = 0.0
end if
end do
end do
hwi_NEA = dim_sum(dim_sum(hwi_trend_ano_weighted2(50-10:75-10, 115:160))) / dim_sum(dim_sum(grids2(50-10:75-10, 115:160)))
estimate_NEA = dim_sum(dim_sum(Estimated_hwi_weighted2(50-10:75-10, 115:160))) / dim_sum(dim_sum(grids2(50-10:75-10, 115:160)))





hwi_trend_ano_weighted3 = hwi_trend_ano_weighted
Estimated_hwi_weighted3 = Estimated_hwi_weighted
grids3 = grids65
do ilat = 0,dim_0(0)-1
do jlon = 0,dim_0(1)-1
if ((ilat-30+10).le.(0.75*(268-jlon))) then
  hwi_trend_ano_weighted3(ilat,jlon) = 0.0
  Estimated_hwi_weighted3(ilat,jlon) = 0.0
  grids3(ilat,jlon) = 0.0
end if
end do
end do
hwi_ENA2 = dim_sum(dim_sum(hwi_trend_ano_weighted3(30-10:50-10,243:290))) / dim_sum(dim_sum(grids3(30-10:50-10,243:290)))
estimate_ENA2 = dim_sum(dim_sum(Estimated_hwi_weighted3(30-10:50-10,243:290))) / dim_sum(dim_sum(grids3(30-10:50-10,243:290)))






print(dim_sum(dim_sum(hwi_trend_ano_weighted(45-10:65-10, 25:55))))
print(dim_sum(dim_sum(hwi_trend_ano_weighted(35-10:55-10, WE_lon))))
print(dim_sum(grids65(45-10:65-10, 25:55)))
print(dim_sum(dim_sum(grids65(35-10:55-10, WE_lon))))

print(" -----------  Regional Contribution  ---------- ")

region_trend = new((/2,8/), float)
region_trend(0,0) = hwi_NWNA * (-1.0)
region_trend(0,1) = hwi_ENA2 * (-1.0)
region_trend(0,2) = hwi_WE
region_trend(0,3) = hwi_Europe
region_trend(0,4) = hwi_EE
region_trend(0,5) = hwi_SCA * (-1.0)
region_trend(0,6) = hwi_EA
region_trend(0,7) = hwi_NEA

region_trend(1,0) = estimate_NWNA * (-1.0)
region_trend(1,1) = estimate_ENA2 * (-1.0)
region_trend(1,2) = estimate_WE
region_trend(1,3) = estimate_Europe
region_trend(1,4) = estimate_EE
region_trend(1,5) = estimate_SCA * (-1.0)
region_trend(1,6) = estimate_EA
region_trend(1,7) = estimate_NEA

region_contribution = region_trend(1,:)/region_trend(0,:)
print(""+region_trend)
print(""+region_contribution)




;;--------------------------------------------------------------------------------------;;
;;-------------------------- Seven regions: significant test ---------------------------;;
degree1_weighted = new((/44, dim_0(0),dim_0(1)/), float)
do ilat = 0,dim_0(0)-1
  degree1_weighted(:,ilat,:) = degree1(:,ilat,:) * cos((ilat+10) / 180.0 * 3.14159)
end do



teleconnection_EE = dim_sum(dim_sum(degree1_weighted(:,45-10:65-10, 25:55)))
teleconnection_WE = dim_sum(dim_sum(degree1_weighted(:,35-10:55-10, WE_lon)))
teleconnection_Europe = teleconnection_EE+teleconnection_WE
teleconnection_EA = dim_sum(dim_sum(degree1_weighted(:,28-10:34-10, 100:120)))
teleconnection_SCA = dim_sum(dim_sum(degree1_weighted(:,22-10:45-10, 65:80)))
teleconnection_NWNA = dim_sum(dim_sum(degree1_weighted(:,60-10:70-10, 200:260)))
teleconnection_ENA = dim_sum(dim_sum(degree1_weighted(:,35-10:50-10, 260:290)))


degree1_weighted2 = degree1_weighted
do ilat = 0,dim_0(0)-1
do jlon = 0,dim_0(1)-1
if ((jlon-125.lt.(ilat-(60-10)))) then
  degree1_weighted2(:,ilat,jlon) = 0.0
end if
end do
end do
teleconnection_NEA = dim_sum(dim_sum(degree1_weighted2(:,50-10:75-10, 115:160)))

degree1_weighted3 = degree1_weighted
do ilat = 0,dim_0(0)-1
do jlon = 0,dim_0(1)-1
if ((ilat-30+10).le.(0.75*(268-jlon))) then
  degree1_weighted3(:,ilat,jlon) = 0.0
end if
end do
end do
teleconnection_ENA2 = dim_sum(dim_sum(degree1_weighted3(:,30-10:50-10,245:290)))



teleconnection_Cor_EE = escorc_n(fspan(1, 44, 44), teleconnection_EE, 0, 0)
teleconnection_Cor_WE = escorc_n(fspan(1, 44, 44), teleconnection_WE, 0, 0)
teleconnection_Cor_Europe = escorc_n(fspan(1, 44, 44), teleconnection_Europe, 0, 0)
teleconnection_Cor_EA = escorc_n(fspan(1, 44, 44), teleconnection_EA, 0, 0)
teleconnection_Cor_NEA = escorc_n(fspan(1, 44, 44), teleconnection_NEA, 0, 0)
teleconnection_Cor_SCA = escorc_n(fspan(1, 44, 44), teleconnection_SCA, 0, 0)
teleconnection_Cor_NWNA = escorc_n(fspan(1, 44, 44), teleconnection_NWNA, 0, 0)
teleconnection_Cor_ENA = escorc_n(fspan(1, 44, 44), teleconnection_ENA, 0, 0)
teleconnection_Cor_ENA2 = escorc_n(fspan(1, 44, 44), teleconnection_ENA2, 0, 0)



teleconnection_rtest_EE = rtest(teleconnection_Cor_EE, 44, 0)
teleconnection_rtest_WE = rtest(teleconnection_Cor_WE, 44, 0)
teleconnection_rtest_Europe = rtest(teleconnection_Cor_Europe, 44, 0)
teleconnection_rtest_EA = rtest(teleconnection_Cor_EA, 44, 0)
teleconnection_rtest_NEA = rtest(teleconnection_Cor_NEA, 44, 0)
teleconnection_rtest_SCA = rtest(teleconnection_Cor_SCA, 44, 0)
teleconnection_rtest_NWNA = rtest(teleconnection_Cor_NWNA, 44, 0)
teleconnection_rtest_ENA = rtest(teleconnection_Cor_ENA, 44, 0)
teleconnection_rtest_ENA2 = rtest(teleconnection_Cor_ENA2, 44, 0)
print(" CC "+estimate_ENA2/hwi_ENA2+" "+teleconnection_rtest_ENA2)



print(" -----------  Regional Trend Test (regression)  ---------- ")
print(teleconnection_rtest_NWNA + " ENA2 " + teleconnection_rtest_ENA2 + " EUR " + teleconnection_rtest_WE + "  " + teleconnection_rtest_Europe +  "  "+ teleconnection_rtest_EE + \
" CR " + teleconnection_rtest_SCA + " EA " + teleconnection_rtest_EA + " NEA " + teleconnection_rtest_NEA )










;;----------------- Pattern correlation between fig (a) and (b) ---------------;;
Observed_1D = new((/dim_0(0)*dim_0(1)/), float)
Estimated_1D = new((/dim_0(0)*dim_0(1)/), float)
do i = 0, dim_0(0)-1
    Observed_1D(i*dim_0(1):i*dim_0(1)+dim_0(1)-1) = hwi_trend_ano(i,:) * cos((i+10)/180*3.14158)
    Estimated_1D(i*dim_0(1):i*dim_0(1)+dim_0(1)-1) = Estimated_hwi(i,:) * cos((i+10)/180*3.14158)  ;; 纬度加权
end do

Observed_number = equiv_sample_size(Observed_1D, 0.001, 0)
Estimated_number = equiv_sample_size(Estimated_1D, 0.001, 0)

cor = escorc_n(Observed_1D, Estimated_1D, 0, 0)
reg = regCoef_n(Observed_1D, Estimated_1D, 0, 0)

print("degree of free  ,  "+Observed_number+" , "+Estimated_number)
print("Correlation  ,  "+cor+", Regression  ,  "+reg)
print("avg = "+dim_avg(Observed_1D)+"  "+dim_avg(Estimated_1D))
print("b = "+(dim_avg(Observed_1D)-reg*dim_avg(Estimated_1D)))








;;------------------------------------------------------------------------------;;
;;---------------- HWI PDF in each summer (1979-2012, 2013-2022) ---------------;;
weight_matrix = new((/44,86,360/), float)
do ilat = 0,85
weight_matrix(:,ilat,:) = cos(ilat / 180.0 * 3.14159)
end do





print("----------  HWI  standardized  by  1979-2000  ---------")
hwi_avg34 = conform(hwi_yearly, dim_avg_n_Wrap(hwi_yearly(0:21,:,:), 0), (/1,2/))
hwi_std34 = conform(hwi_yearly, dim_stddev_n_Wrap(hwi_yearly(0:21,:,:), 0), (/1,2/))
hwi_std34 = where(hwi_std34.le.0.01, 0.01, hwi_std34)
hwi_stz = (hwi_yearly-hwi_avg34)/hwi_std34
copy_VarCoords(hwi_yearly, hwi_stz)
printMinMax(hwi_stz, 1)


strength = fspan(-1.75,9.75,24)
print(strength)
pdf_first22 = new(24, float)
pdf_last22 = new(24, float)
pdfP_first22 = new(24, float)
pdfP_last22 = new(24, float)
pdfN_first22 = new(24, float)
pdfN_last22 = new(24, float)
Estimated_hwi_3D = conform(hwi_stz(:,10:85,:), Estimated_hwi, (/1,2/))

do i = 0,23
hwi_new = where(hwi_stz.ge.strength(i)-0.25.and.hwi_stz.lt.strength(i)+0.25, 1.0, 0.0)
if(i.eq.23)then
  hwi_new = where(hwi_stz.ge.strength(i)-0.25, 1.0, 0.0)
end if
hwi_new = hwi_new * weight_matrix

hwi_newP = where(Estimated_hwi_3D.gt.2.0, hwi_new(:,10:85,:), 0.0)    ;; 正负遥相关变化
hwi_newN = where(Estimated_hwi_3D.lt.-2.0, hwi_new(:,10:85,:), 0.0)

pdf_first22(i) = dim_sum(dim_sum(dim_sum(hwi_new(0:21,:,:))))
pdf_last22(i) = dim_sum(dim_sum(dim_sum(hwi_new(22:43,:,:))))
pdfP_first22(i) = dim_sum(dim_sum(dim_sum(hwi_newP(0:21,:,:))))
pdfP_last22(i) = dim_sum(dim_sum(dim_sum(hwi_newP(22:43,:,:))))
pdfN_first22(i) = dim_sum(dim_sum(dim_sum(hwi_newN(0:21,:,:))))
pdfN_last22(i) = dim_sum(dim_sum(dim_sum(hwi_newN(22:43,:,:))))
end do
pdf_first22 = pdf_first22 / dim_sum(pdf_first22) * 100.0
pdf_last22 = pdf_last22 / dim_sum(pdf_last22) * 100.0
pdfP_first22 = pdfP_first22 / dim_sum(pdfP_first22) * 100.0
pdfP_last22 = pdfP_last22 / dim_sum(pdfP_last22) * 100.0
pdfN_first22 = pdfN_first22 / dim_sum(pdfN_first22) * 100.0
pdfN_last22 = pdfN_last22 / dim_sum(pdfN_last22) * 100.0
pdf_first22 = where(pdf_first22.le.0.0, 0.001, pdf_first22)
pdf_last22 = where(pdf_last22.le.0.0, 0.001, pdf_last22)
pdfP_first22 = where(pdfP_first22.le.0.0, 0.001, pdfP_first22)
pdfP_last22 = where(pdfP_last22.le.0.0, 0.001, pdfP_last22)
pdfN_first22 = where(pdfN_first22.le.0.0, 0.001, pdfN_first22)
pdfN_last22 = where(pdfN_last22.le.0.0, 0.001, pdfN_last22)

print(pdf_last22 - pdf_first22)
print(pdfP_last22 - pdfP_first22)
print(pdfN_last22 - pdfN_first22)



print("        ---------------------        ")




hwi_new = where(hwi_stz.ge.-2.0, 1.0, 0.0)
Sum = dim_sum(dim_sum(dim_sum(hwi_new(0:21,40:65,:) * weight_matrix(0:21,40:65,:))))
hwi_new = where(hwi_stz.ge.2.0, 1.0, 0.0)
hwi_new = hwi_new * weight_matrix
pdf_A = dim_sum(dim_sum(dim_sum(hwi_new(0:21,40:65,:)))) / Sum
pdf_B = dim_sum(dim_sum(dim_sum(hwi_new(22:43,40:65,:)))) / Sum
print(" P+N, 2001-2022 > 2.0std "+ pdf_A + " "+ pdf_B +" "+ pdf_B/pdf_A)


hwi_new = where(hwi_stz.ge.-2.0, 1.0, 0.0)
hwi_newP = where(Estimated_hwi_3D.gt.0.0, hwi_new(:,10:85,:), 0.0)
Sum = dim_sum(dim_sum(dim_sum(hwi_newP(0:21,10-10:85-10,:) * weight_matrix(0:21,10-10:85-10,:))))
hwi_new = where(hwi_stz.ge.2.0, 1.0, 0.0)
hwi_new = hwi_new * weight_matrix
hwi_newP = where(Estimated_hwi_3D.gt.0.0, hwi_new(:,10:85,:), 0.0)
pdf_A = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newP(0:21,10-10:85-10,:),2),1),0) / Sum
pdf_B = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newP(22:43,10-10:85-10,:),2),1),0) / Sum
print(" P, 2001-2022 > 2.0std "+ pdf_A + " "+ pdf_B +" "+ pdf_B/pdf_A)


hwi_new = where(hwi_stz.ge.-2.0, 1.0, 0.0)
hwi_newN = where(Estimated_hwi_3D.lt.0.0, hwi_new(:,10:85,:), 0.0)
Sum = dim_sum(dim_sum(dim_sum(hwi_newN(0:21,10-10:85-10,:) * weight_matrix(0:21,10-10:85-10,:))))
hwi_new = where(hwi_stz.ge.2.0, 1.0, 0.0)
hwi_new = hwi_new * weight_matrix
hwi_newN = where(Estimated_hwi_3D.lt.0.0, hwi_new(:,10:85,:), 0.0)
pdf_A = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newN(0:21,10-10:85-10,:),2),1),0) / Sum
pdf_B = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newN(22:43,10-10:85-10,:),2),1),0) / Sum
print(" N, 2001-2022 > 2.0std "+ pdf_A + " "+ pdf_B +" "+ pdf_B/pdf_A)





hwi_new = where(hwi_stz.ge.-2.0, 1.0, 0.0)
hwi_newP = where(Estimated_hwi_3D.gt.0.5, hwi_new(:,10:85,:), 0.0)
Sum = dim_sum(dim_sum(dim_sum(hwi_newP(0:21,10-10:85-10,:) * weight_matrix(0:21,10-10:85-10,:))))
hwi_new = where(hwi_stz.ge.2.0, 1.0, 0.0)
hwi_new = hwi_new * weight_matrix
hwi_newP = where(Estimated_hwi_3D.gt.0.5, hwi_new(:,10:85,:), 0.0)
pdf_A = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newP(0:21,10-10:85-10,:),2),1),0) / Sum
pdf_B = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newP(22:43,10-10:85-10,:),2),1),0) / Sum
print(" P>0.5, 2001-2022 > 2.0std "+ pdf_A + " "+ pdf_B +" "+ pdf_B/pdf_A)

hwi_new = where(hwi_stz.ge.-2.0, 1.0, 0.0)
hwi_newP = where(Estimated_hwi_3D.gt.1.0, hwi_new(:,10:85,:), 0.0)
Sum = dim_sum(dim_sum(dim_sum(hwi_newP(0:21,10-10:85-10,:) * weight_matrix(0:21,10-10:85-10,:))))
hwi_new = where(hwi_stz.ge.2.0, 1.0, 0.0)
hwi_new = hwi_new * weight_matrix
hwi_newP = where(Estimated_hwi_3D.gt.1.0, hwi_new(:,10:85,:), 0.0)
pdf_A = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newP(0:21,10-10:85-10,:),2),1),0) / Sum
pdf_B = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newP(22:43,10-10:85-10,:),2),1),0) / Sum
print(" P>1, 2001-2022 > 2.0std "+ pdf_A + " "+ pdf_B +" "+ pdf_B/pdf_A)

hwi_new = where(hwi_stz.ge.-2.0, 1.0, 0.0)
hwi_newP = where(Estimated_hwi_3D.gt.1.5, hwi_new(:,10:85,:), 0.0)
Sum = dim_sum(dim_sum(dim_sum(hwi_newP(0:21,10-10:85-10,:) * weight_matrix(0:21,10-10:85-10,:))))
hwi_new = where(hwi_stz.ge.2.0, 1.0, 0.0)
hwi_new = hwi_new * weight_matrix
hwi_newP = where(Estimated_hwi_3D.gt.1.5, hwi_new(:,10:85,:), 0.0)
pdf_A = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newP(0:21,10-10:85-10,:),2),1),0) / Sum
pdf_B = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newP(22:43,10-10:85-10,:),2),1),0) / Sum
print(" P>1.5, 2001-2022 > 2.0std "+ pdf_A + " "+ pdf_B +" "+ pdf_B/pdf_A)

hwi_new = where(hwi_stz.ge.-2.0, 1.0, 0.0)
hwi_newP = where(Estimated_hwi_3D.gt.2.0, hwi_new(:,10:85,:), 0.0)
Sum = dim_sum(dim_sum(dim_sum(hwi_newP(0:21,10-10:85-10,:) * weight_matrix(0:21,10-10:85-10,:))))
hwi_new = where(hwi_stz.ge.2.0, 1.0, 0.0)
hwi_new = hwi_new * weight_matrix
hwi_newP = where(Estimated_hwi_3D.gt.2.0, hwi_new(:,10:85,:), 0.0)
pdf_A = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newP(0:21,10-10:85-10,:),2),1),0) / Sum
pdf_B = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newP(22:43,10-10:85-10,:),2),1),0) / Sum
print(" P>2.0, 2001-2022 > 2.0std "+ pdf_A + " "+ pdf_B +" "+ pdf_B/pdf_A)




hwi_new = where(hwi_stz.ge.-2.0, 1.0, 0.0)
hwi_newN = where(Estimated_hwi_3D.lt.-0.5, hwi_new(:,10:85,:), 0.0)
Sum = dim_sum(dim_sum(dim_sum(hwi_newN(0:21,10-10:85-10,:) * weight_matrix(0:21,10-10:85-10,:))))
hwi_new = where(hwi_stz.ge.2.0, 1.0, 0.0)
hwi_new = hwi_new * weight_matrix
hwi_newN = where(Estimated_hwi_3D.lt.-0.5, hwi_new(:,10:85,:), 0.0)
pdf_A = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newN(0:21,10-10:85-10,:),2),1),0) / Sum
pdf_B = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newN(22:43,10-10:85-10,:),2),1),0) / Sum
print(" N<0.5, 2001-2022 > 2.0std "+ pdf_A + " "+ pdf_B +" "+ pdf_B/pdf_A)

hwi_new = where(hwi_stz.ge.-2.0, 1.0, 0.0)
hwi_newN = where(Estimated_hwi_3D.lt.-1.0, hwi_new(:,10:85,:), 0.0)
Sum = dim_sum(dim_sum(dim_sum(hwi_newN(0:21,10-10:85-10,:) * weight_matrix(0:21,10-10:85-10,:))))
hwi_new = where(hwi_stz.ge.2.0, 1.0, 0.0)
hwi_new = hwi_new * weight_matrix
hwi_newN = where(Estimated_hwi_3D.lt.-1.0, hwi_new(:,10:85,:), 0.0)
pdf_A = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newN(0:21,10-10:85-10,:),2),1),0) / Sum
pdf_B = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newN(22:43,10-10:85-10,:),2),1),0) / Sum
print(" N<1, 2001-2022 > 2.0std "+ pdf_A + " "+ pdf_B +" "+ pdf_B/pdf_A)

hwi_new = where(hwi_stz.ge.-2.0, 1.0, 0.0)
hwi_newN = where(Estimated_hwi_3D.lt.-1.5, hwi_new(:,10:85,:), 0.0)
Sum = dim_sum(dim_sum(dim_sum(hwi_newN(0:21,10-10:85-10,:) * weight_matrix(0:21,10-10:85-10,:))))
hwi_new = where(hwi_stz.ge.2.0, 1.0, 0.0)
hwi_new = hwi_new * weight_matrix
hwi_newN = where(Estimated_hwi_3D.lt.-1.5, hwi_new(:,10:85,:), 0.0)
pdf_A = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newN(0:21,10-10:85-10,:),2),1),0) / Sum
pdf_B = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newN(22:43,10-10:85-10,:),2),1),0) / Sum
print(" N<1.5, 2001-2022 > 2.0std "+ pdf_A + " "+ pdf_B +" "+ pdf_B/pdf_A)

hwi_new = where(hwi_stz.ge.-2.0, 1.0, 0.0)
hwi_newN = where(Estimated_hwi_3D.lt.-2.0, hwi_new(:,10:85,:), 0.0)
Sum = dim_sum(dim_sum(dim_sum(hwi_newN(0:21,10-10:85-10,:) * weight_matrix(0:21,10-10:85-10,:))))
hwi_new = where(hwi_stz.ge.2.0, 1.0, 0.0)
hwi_new = hwi_new * weight_matrix
hwi_newN = where(Estimated_hwi_3D.lt.-2.0, hwi_new(:,10:85,:), 0.0)
pdf_A = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newN(0:21,10-10:85-10,:),2),1),0) / Sum
pdf_B = dim_sum_n(dim_sum_n(dim_sum_n(hwi_newN(22:43,10-10:85-10,:),2),1),0) / Sum
print(" N<2, 2001-2022 > 2.0std "+ pdf_A + " "+ pdf_B +" "+ pdf_B/pdf_A)






;************************************************** plots **********************************************
res = True
res@tmBorderThicknessF = 1.2
res@gsnDraw      =  False
res@gsnFrame     =  False
res@gsnRightString       = ""
res@gsnLeftString        = ""
res@gsnLeftStringFontHeightF   = 0.04
res@gsnRightStringFontHeightF  = 0.04
res@gsnCenterStringFontHeightF  = 0.018
res@tmXTLabelsOn  = False
res@tmYRLabelsOn  = False
vcres = res
bres = res
cnres = res

res@gsnAddCyclic =  True
res@tmXTOn        = False
res@tmYROn        = False
res@tmXBOn        = True
res@tmYLOn        = True
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015
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




res@mpMinLatF = 10
res@mpMaxLatF = 85
res@mpMinLonF = 0
res@mpMaxLonF = 360
res@mpCenterLonF = 180-165-3
res@mpFillOn                    = False
res@mpOutlineOn                 = True
res@mpGeophysicalLineThicknessF = 0.6
res@mpGeophysicalLineColor      = "gray10";"gray20"
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
res@pmLabelBarWidthF     = 0.08
res@pmLabelBarHeightF    = 0.145
res@pmLabelBarOrthogonalPosF = 0.012
res@lbLabelFontHeightF  = 0.012
;res@cnMissingValFillColor = "white"
;res@gsnYRefLine = 0.0
;res@gsnYRefLineColor = "black"
;res@gsnYRefLineThicknessF = 5.0
res@cnLevelSelectionMode = "ExplicitLevels"              
; res@cnLevels             = (/140.0,155.0,170.0,185.0,200.0,215.0,230.0,245.0/)
; res@cnFillColors         = (/16,15,13,12,8,6,5,3,-1/)
res@cnFillPalette        = "precip2_17lev" 


res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.009
res@tmXBMajorLengthF = 0.009
res@tmYLMinorLengthF = 0.006
res@tmXBMinorLengthF = 0.006
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.01
res@tmYLMajorOutwardLengthF = 0.01
res@tmXBMinorOutwardLengthF = 0.007
res@tmYLMinorOutwardLengthF = 0.007





;vcres@mpFillDrawOrder         = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.255
;vcres@vcRefAnnoSide           = "TopRight"
vcres@vcGlyphStyle            = "CurlyVector"
vcres@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
vcres@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref
vcres@vcMinDistanceF          = 0.06             ; thin out vectors
vcres@vcLineArrowColor        = "black"           ; change vector color
vcres@vcRefAnnoOn             = True
vcres@vcRefAnnoString1        = "4"
vcres@vcRefMagnitudeF         = 4
vcres@vcRefLengthF            = 0.034            ; ref vec length
vcres@vcRefAnnoFontHeightF    = 0.022   ;参考箭头字体大小
vcres@vcRefAnnoString1On      = True
vcres@vcRefAnnoString2On      = False
vcres@vcLineArrowThicknessF   = 3.0            ; make vectors larger
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
lires@tmBorderThicknessF = 1.0
lires@gsnDraw        = False
lires@gsnFrame       = False
lires@lgPerimOn = False
lires@lgLabelFontHeightF = 0.018
lires@vpWidthF  = 0.473
lires@vpHeightF = 0.15
lires@vpXF      = 0.06
lires@vpYF      = 0.55


lires@tiXAxisOn      = False
lires@tiYAxisOn      = False
lires@tmXTLabelsOn   = False
lires@tmXTOn         = False


lires@gsnXYBarChart = True
lires@gsnXYBarChartBarWidth = 0.8
lires@gsnYRefLine = 0.0
; lires@gsnXYBarChartColors = (/"violetred2","royalblue3","violetred2","violetred2","royalblue3","royalblue3","violetred2"/)

lires@trYMinF = 0.0
lires@trYMaxF = 7.0
lires@trXMinF = -0.2 - 0.2
lires@trXMaxF = 24.2 + 0.2


lires@tmXBLabelFontHeightF = 0.012
lires@tmYLLabelFontHeightF = 0.0135
lires@tmYRLabelFontHeightF = 0.0125
lires@gsnLeftStringFontHeightF   = 0.022
lires@gsnRightStringFontHeightF  = 0.022
lires@tmYRLabelsOn   = False
lires@tmYROn         = True



lires@tmYLMode = "Explicit"
lires@tmYLValues = fspan(0.0, 6.0, 4)
lires@tmYLLabels = fspan(0.0, 6.0, 4)
lires@tmYLMinorOn = True
lires@tmYLMinorValues = fspan(0.0, 7.0, 15)

lires@tmYRMode = "Explicit"
lires@tmYRValues = fspan(0.0, 6.0, 4)
lires@tmYRLabels = fspan(0.0, 6.0, 4)
lires@tmYRMinorOn = True
lires@tmYRMinorValues = fspan(0.0, 7.0, 15)

lires@tmXBMode = "Explicit"
lires@tmXBValues = (/2.0,5.0,8.0,11.0,14.0,17.0,20.0,23.0/) - 0.5
lires@tmXBLabels = (/"NWNA","ENA","WE","WE+EE","EE","SCA","EA","NEA"/)
lires@tmXBMinorOn = False

lires@tmYLMajorThicknessF = 1.5
lires@tmYRMajorThicknessF = 1.5
lires@tmXBMajorThicknessF = 1.5
lires@tmYLMinorThicknessF = 1.0
lires@tmYRMinorThicknessF = 1.0
lires@tmXBMinorThicknessF = 1.0
lires@tmYLMajorLengthF = 0.009
lires@tmXBMajorLengthF = 0.009
lires@tmYLMinorLengthF = 0.006
lires@tmXBMinorLengthF = 0.006
lires@tmXBLabelDeltaF = -0.4
lires@tmYLLabelDeltaF = -0.4
lires@tmXBMajorOutwardLengthF = 0.009
lires@tmYLMajorOutwardLengthF = 0.009
lires@tmXBMinorOutwardLengthF = 0.006
lires@tmYLMinorOutwardLengthF = 0.006







leres = True
leres@tmBorderThicknessF = 1.2
leres@gsnDraw        = False
leres@gsnFrame       = False
leres@tiXAxisOn      = False
leres@tiYAxisOn      = False
leres@tmXTLabelsOn   = False
leres@tmXTOn         = False

leres@trYMinF        = 10
leres@trYMaxF        = 85

leres@tmXBLabelFontHeightF = 0.023
leres@tmYLLabelFontHeightF = 0.023
leres@tmYRLabelFontHeightF = 0.023
leres@gsnLeftStringFontHeightF   = 0.022
leres@gsnRightStringFontHeightF  = 0.022
leres@tmYLLabelsOn   = False
leres@tmYLOn         = False

leres@tmYROn         = True
leres@tmYRLabelsOn   = True
leres@tmYRMode = "Explicit"
leres@tmYRValues = fspan(0, 80, 5)
leres@tmYRLabels = (/"0","20N","40N","60N","80N"/)
leres@tmYRMinorOn = True
leres@tmYRMinorValues = fspan(0, 80, 9)

leres@tmYLValues = fspan(0, 80, 5)
leres@tmYLLabels = (/"0","20N","40N","60N","80N"/)
leres@tmYLMinorOn = True
leres@tmYLMinorValues = fspan(0, 80, 9)
;leres@tmYLMode = "Explicit"
;leres@tmYLValues = fspan(1, 71, 8)
;leres@tmYLLabels = (/"1-April","11-April","21-April","1-May","11-May","21-May","31-May","10-June"/)

leres@xyMarkLineModes = (/"Lines","Lines","Lines"/)
;leres@xyMarkers = (/4,1,1/)
;leres@xyMarkerSizes = (/0.015,0.001,0.001/)
;leres@xyMarkerColors= (/"black","black","black"/)
;leres@xyMarkerThicknesses = (/4.0,2.0,1.0/)

leres@lgPerimOn = False
leres@lgLabelFontHeightF = 0.018
leres@vpWidthF  = 0.34
leres@vpHeightF = 0.275
leres@vpXF      = 0.06
leres@vpYF      = 0.75

leres@gsnLeftString  = ""
leres@gsnRightString = ""
leres@xyMarkLineModes = (/"Lines","Lines","Lines"/)
; leres@xyMarkers = (/1,1,1/)
; leres@xyMarkerSizes = (/0.015,0.001,0.001/)
; leres@xyMarkerColors= (/"black","black","black"/)
; leres@xyMarkerThicknesses = (/4.0,2.0,1.0/)
leres@xyLineThicknesses  = (/5.0,8.0,8.0/)
leres@xyLineColors       = (/"red","blue","orange"/)
leres@xyDashPatterns     = (/0,4,2/)

leres@tmYLMajorThicknessF = 1.5
leres@tmYRMajorThicknessF = 1.5
leres@tmXBMajorThicknessF = 1.5
leres@tmYLMinorThicknessF = 1.0
leres@tmYRMinorThicknessF = 1.0
leres@tmXBMinorThicknessF = 1.0
leres@tmYLMajorLengthF = 0.016
leres@tmYRMajorLengthF = 0.016
leres@tmXBMajorLengthF = 0.016
leres@tmYLMinorLengthF = 0.012
leres@tmYRMinorLengthF = 0.012
leres@tmXBMinorLengthF = 0.012
leres@tmXBLabelDeltaF = -0.25
leres@tmYLLabelDeltaF = -0.25
leres@tmXBMajorOutwardLengthF = 0.016
leres@tmYLMajorOutwardLengthF = 0.016
leres@tmYRMajorOutwardLengthF = 0.016
leres@tmXBMinorOutwardLengthF = 0.012
leres@tmYLMinorOutwardLengthF = 0.012
leres@tmYRMinorOutwardLengthF = 0.012





cnres@tmBorderThicknessF = 1.2
cnres@gsnDraw        = False
cnres@gsnFrame       = False
cnres@lgPerimOn = False
cnres@lgLabelFontHeightF = 0.018
cnres@vpWidthF  = 0.33
cnres@vpHeightF = 0.22
cnres@vpXF      = 0.06
cnres@vpYF      = 0.5

cnres@trYMinF = -4.0
cnres@trYMaxF = 4.0
cnres@trXMinF = -7.0
cnres@trXMaxF = 7.0


cnres@tmXTOn        = False
cnres@tmYROn        = False
cnres@tmXBOn        = True
cnres@tmYLOn        = True
cnres@tmXBLabelFontHeightF = 0.015
cnres@tmYLLabelFontHeightF = 0.015
cnres@tmYRLabelFontHeightF = 0.015
cnres@gsnLeftStringFontHeightF   = 0.022
cnres@gsnRightStringFontHeightF  = 0.022
cnres@tmYRLabelsOn   = False
cnres@tmYROn         = True

cnres@tmYLMode = "Explicit"
cnres@tmYLValues = fspan(-4.0, 4.0, 5)
cnres@tmYLLabels = fspan(-4.0, 4.0, 5)
cnres@tmYLMinorOn = True
cnres@tmYLMinorValues = fspan(-3.5, 4.5, 17)

cnres@tmXBMode = "Explicit"
cnres@tmXBValues = fspan(-6.0, 6.0, 7)
cnres@tmXBLabels = fspan(-6.0, 6.0, 7)
cnres@tmXBMinorOn = True
cnres@tmXBMinorValues = fspan(-7.0, 7.0, 15)

cnres@cnFillOn       = True
cnres@cnLinesOn      = False
cnres@cnLineLabelsOn = False
cnres@lbLabelBarOn   = False
cnres@cnInfoLabelOn  = False
cnres@cnLevelSelectionMode = "ExplicitLevels"

cnres@tmYLMajorThicknessF = 1.5
cnres@tmXBMajorThicknessF = 1.5
cnres@tmYLMinorThicknessF = 1.0
cnres@tmXBMinorThicknessF = 1.0
cnres@tmYLMajorLengthF = 0.009
cnres@tmXBMajorLengthF = 0.009
cnres@tmYLMinorLengthF = 0.0065
cnres@tmXBMinorLengthF = 0.0065
cnres@tmXBLabelDeltaF = -0.4
cnres@tmYLLabelDeltaF = -0.4
cnres@tmXBMajorOutwardLengthF = 0.01
cnres@tmYLMajorOutwardLengthF = 0.01
cnres@tmXBMinorOutwardLengthF = 0.0065
cnres@tmYLMinorOutwardLengthF = 0.0065



;************************* dotting ***************************
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
contour10 = new(14, "graphic")
contour11 = new(14, "graphic")
topooo = new(14, "graphic") 



;;---------------------- color map ----------------------;;
colors1 = read_colormap_file("MPL_OrRd")
colors2 = read_colormap_file("WhiteBlue")
colors = colors1(0:20,:)
colors(1,:) = colors1(11,:)
colors(2,:) = colors1(31,:)
colors(3,:) = colors1(51,:)
colors(4,:) = colors1(56,:)

colors(11,:) = colors2(45,:)
colors(12,:) = colors2(75,:)
colors(13,:) = colors2(115,:)
colors(14,:) = colors2(130,:)
delete(res@cnFillPalette)
delete(res@cnLevels)
delete(res@cnFillColors)




colorsA = colors1(0:20,:)
colorsA(1,:) = colors1(2,:)
colorsA(2,:) = colors1(10,:)
colorsA(3,:) = colors1(23,:)
colorsA(4,:) = colors1(37,:)
colorsA(5,:) = colors1(50,:)
colorsA(6,:) = colors1(61,:)
colorsA(7,:) = colors1(72,:)
colorsA(8,:) = colors1(83,:)
colorsA(9,:) = colors1(94,:)
colorsA(11,:) = colors2(10,:)
colorsA(12,:) = colors2(30,:)
colorsA(13,:) = colors2(55,:)
colorsA(14,:) = colors2(85,:)
colorsA(15,:) = colors2(103+1,:)
colorsA(16,:) = colors2(121+2,:)
colorsA(17,:) = colors2(139+4,:)
colorsA(18,:) = colors2(157+5,:)
colorsA(19,:) = colors2(185,:)

; res@cnFillColors = (/11,1,2,3,4/)
colorsB = colors1(0:20,:)
colorsB(1,:) = colors1(2,:)
colorsB(2,:) = colors1(15,:)
colorsB(3,:) = colors1(35,:)
colorsB(4,:) = colors1(50,:)
colorsB(5,:) = colors1(60,:)
colorsB(6,:) = colors1(70,:)
colorsB(7,:) = colors1(82,:)
colorsB(8,:) = colors1(94,:)
colorsB(11,:) = colors2(20,:)
colorsB(12,:) = colors2(55,:)
colorsB(13,:) = colors2(85,:)
colorsB(14,:) = colors2(108,:)
colorsB(15,:) = colors2(112+4,:)
colorsB(16,:) = colors2(130+4,:)
colorsB(17,:) = colors2(148+4,:)
colorsB(18,:) = colors2(166,:)





;;-----------------------------------------------------
pltType = "pdf"
; pltType = "eps"
pltName = "Figure3"
wks = gsn_open_wks(pltType, pltName)
res@cnFillPalette = colorsA
res@cnLevels     = (/-5.0,-4.0,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0/)
res@cnFillColors = (/19,18,17,16,15,14,13,12,11,1,2,3,4,5,6,7,8,9/)
res@mpMinLatF = 10
res@mpMaxLatF = 90

res@tmXBLabelsOn = False
plot(0) = gsn_csm_contour_map_ce(wks, hwi_trend_ano, res)




;;------------------------ fig (a) boxes -----------------------
y1=45.0
y2=65
x1=25
x2=55
x=(/x1,x2,x2,x1,x1/)  ;; eastern Europe
y=(/y2,y2,y1,y1,y2/)
   lnres   =    True
   lnres@gsLineThicknessF=2.8
   lnres@gsnFrame = False

   lnres@gsLineColor="white"
   lnres@gsLineDashPattern=0
   dum11=gsn_add_polyline(wks,plot(0),x,y,lnres)
   lnres@gsLineColor="violetred3"
   lnres@gsLineDashPattern=14
   dum1=gsn_add_polyline(wks,plot(0),x,y,lnres)

y1=35.0
y2=55
x1=-10
x2=25
x=(/x1,x2,x2,x1,x1/)  ;; western Europe
y=(/y2,y2,y1,y1,y2/)
   lnres@gsLineColor="white"
   lnres@gsLineDashPattern=0
   dum22=gsn_add_polyline(wks,plot(0),x,y,lnres)
   lnres@gsLineColor="violetred3"
   lnres@gsLineDashPattern=14
   dum2=gsn_add_polyline(wks,plot(0),x,y,lnres)

; 28-10:34-10, 100:120
y1=27.0
y2=35
x1=100
x2=120
x=(/x1,x2,x2,x1,x1/)  ;; East Asia
y=(/y2,y2,y1,y1,y2/)
   lnres@gsLineColor="violetred3"
   dum3=gsn_add_polyline(wks,plot(0),x,y,lnres)

y1=22.0
y2=45
x1=65
x2=80
x=(/x1,x2,x2,x1,x1/)  ;; South Asia-Central Asia
y=(/y2,y2,y1,y1,y2/)
   lnres@gsLineColor="royalblue4"
   dum4=gsn_add_polyline(wks,plot(0),x,y,lnres)

; 50-10:75-10, 115:160
y1=50.0
y2=75
x1=125
x2=160
xt=(/140,x2,x2,x1, 115, x1,140/)  ;; Northeastern Asia
yt=(/y2,y2,y1,y1, 50, 60,y2/)
   lnres@gsLineColor="violetred3"
   dum5=gsn_add_polyline(wks,plot(0), xt,yt,lnres)



; 60-10:70-10, 200:260
y1=60.0
y2=70
x1=200
x2=260
x=(/x1,x2,x2,x1,x1/)  ;; Northwestern North America
y=(/y2,y2,y1,y1,y2/)
   lnres@gsLineColor="royalblue4"
   lnres@gsLineDashPattern=14
   dum6=gsn_add_polyline(wks,plot(0),x,y,lnres)


xc = (/290,250-2,270-2,290,290/)  ;; eastern North America
yc = (/45,45,30,30,45/)
   lnres@gsLineColor="royalblue4"
   lnres@gsLineDashPattern=14
   dum7=gsn_add_polyline(wks,plot(0),xc,yc,lnres)




;;------------------------ fig (a) boxes -----------------------
  txres               = True
  txres@txFontHeightF = 0.016
  txres@txFont = "helvetica-bold"
  txres@txFontHeightF = 0.013

  txres@txFontColor = "royalblue4"
  eur1 = gsn_add_text(wks,plot(0), "SCA", 63,16.0, txres)
  eur7 = gsn_add_text(wks,plot(0), "NWNA", 213,52.0, txres)
  eur8 = gsn_add_text(wks,plot(0), "ENA", 295,24.5, txres)

  txres@txFontColor = "violetred3"
  eur4 = gsn_add_text(wks,plot(0), "EE", 45,71.5, txres)
  eur5 = gsn_add_text(wks,plot(0), "WE", -20,42.0, txres)
  eur6 = gsn_add_text(wks,plot(0), "EA", 129,30.0, txres)
  eur7 = gsn_add_text(wks,plot(0), "NEA", 165,43.5, txres)







res@tmXBLabelsOn = True
delete(res@cnFillPalette)
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnFillPalette = colorsB
res@cnLevels     = (/-2.0,-1.0,-0.5,0.0,0.5,1.0,2.0/)
res@cnFillColors = (/14,13,12,11,1,2,3,4/)
res@mpMinLatF = 10
res@mpMaxLatF = 86
plot(1) = gsn_csm_contour_map_ce(wks, Estimated_hwi, res)

ores@cnLevels    = (/0.0/)
ores@cnFillColor = "green4"
ores@cnFillPatterns       = (/-1,11/) 


ores@cnLevels    = (/0.0/)
ores@cnFillColor = "violetred3"
ores@cnFillPatterns       = (/-1,11/) 
contour3(1) = gsn_csm_contour(wks, Contribution_hwi_P, ores)    
overlay(plot(1), contour3(1))
contour8(1) = gsn_csm_contour(wks, Contribution_hwi_P, ores)    
overlay(plot(1), contour8(1))
ores@cnFillColor = "royalblue4"
contour4(1) = gsn_csm_contour(wks, Contribution_hwi_N, ores)    
overlay(plot(1), contour4(1))
contour9(1) = gsn_csm_contour(wks, Contribution_hwi_N, ores)    
overlay(plot(1), contour9(1))








liresL = lires
liresR = lires
liresL@gsnXYBarChartColors = (/"royalblue3","royalblue3","violetred2","violetred2","violetred2","royalblue3","violetred2","violetred2"/)
liresR@gsnXYBarChartColors = (/"royalblue3","royalblue3","violetred2","violetred2","violetred2","royalblue3","violetred2","violetred2"/)
liresR@gsnXYBarChartPatterns = (/6,6,6,6,6,6,6/)
liresR@gsnXYBarChartFillScaleF = 0.5
x11 = fspan(1.0, 22.0, 8) + 0.1
x22 = fspan(2.0, 23.0, 8) - 0.1
liresL@tmBorderThicknessF = 1.0
liresR@tmBorderThicknessF = 1.0

plot(2) = gsn_csm_x2y2(wks, x11, x22, region_trend(0,:), region_trend(1,:), liresL, liresR)
contour2(2) = gsn_csm_x2y2(wks, (/-5.0,-4.0/), x22, (/100.0,100.0/), region_trend(1,:), liresL, liresR)  ;; 再画一遍条纹框
overlay(plot(2), contour2(2))


;;----------------- fig (c) bars -------------------
lbres = True
lbres@vpWidthF = 0.06
lbres@vpHeightF= 0.03
lbres@lbBoxMajorExtentF = 0.4
lbres@lbLabelFontHeightF = 0.0105
lbres@lbPerimColor = False
lbres@lbFillColors = "violetred2"
lbres@lbFillPatterns = 0
gsn_labelbar_ndc(wks, 1, "", 0.131, 0.428, lbres)
lbres@lbFillColors = "royalblue3"
gsn_labelbar_ndc(wks, 1, "  Observed", 0.176, 0.428, lbres)

lbres@lbFillColors = "violetred2"
lbres@lbFillPatterns = 6
lbres@lbFillScaleF = 0.4
lbres@lbBoxMajorExtentF = 0.39
gsn_labelbar_ndc(wks, 1, "", 0.131, 0.407, lbres)
gsn_labelbar_ndc(wks, 1, "", 0.131, 0.407, lbres)
lbres@lbFillColors = "royalblue3"
; lbres@lbFillScaleF = 5.0
gsn_labelbar_ndc(wks, 1, "           ", 0.176, 0.407, lbres)
gsn_labelbar_ndc(wks, 1, "  Estimated", 0.176, 0.407, lbres)






;;------------------------ figs (a, b) text -----------------------
  txres               = True
  txres@txFont = "helvetica"
  txres@txFontColor = "black"
  txres@txFontHeightF = 0.014
  eum1 = gsn_add_text(wks,plot(0), "Observed", 165,20.0, txres)
  eum3 = gsn_add_text(wks,plot(1), "Estimated", 165,20.0, txres)


  txres@txFontHeightF = 0.012
  eum2 = gsn_add_text(wks,plot(1), "R((a),(b)) = 0.46", -45,32.0, txres)
  eum5 = gsn_add_text(wks,plot(1), "(p<0.001)", -45,22.0, txres)
  eum4 = gsn_add_text(wks,plot(1), "hatched:", -147,47.0, txres)
  euq1 = gsn_add_text(wks,plot(1), "(b)/(a)", -152,37.0, txres)
  txres@txFontColor = "royalblue4"
  euq2 = gsn_add_text(wks,plot(1), "> 0", -135,38.0, txres)
  txres@txFontColor = "violetred3"
  euq3 = gsn_add_text(wks,plot(1), "> 0", -135,31.0, txres)




;;------------------------ fig (c) text -----------------------
  txres@txFontHeightF = 0.0105
  txres@txFontColor = "royalblue4"
  euc1 = gsn_add_text(wks,plot(2), "29%", 1.5, region_trend(0,0)+0.55, txres)
  euc2 = gsn_add_text(wks,plot(2), "42%", 4.5, region_trend(0,1)+0.55, txres) 
  txres@txFontColor = "violetred3"
  euc3 = gsn_add_text(wks,plot(2), "52%", 7.5, region_trend(0,2)+0.55, txres)
  euc4 = gsn_add_text(wks,plot(2), "64%", 10.5, region_trend(0,3)+0.55, txres)
  euc5 = gsn_add_text(wks,plot(2), "72%", 13.5, region_trend(0,4)+0.55, txres)
  txres@txFontColor = "royalblue4"
  euc6 = gsn_add_text(wks,plot(2), "29%", 16.5, region_trend(0,5)+0.55, txres)
  txres@txFontColor = "violetred3"
  euc7 = gsn_add_text(wks,plot(2), "118%", 19.5, region_trend(1,6)+1.1, txres)
  euc8 = gsn_add_text(wks,plot(2), "101%", 22.5, region_trend(0,7)+0.55, txres)


;;------------------------ fig (c) significant circles -----------------------
  pmres = True
  pmres@gsMarkerSizeF    =  0.01  ;-- set size of marker
  pmres@gsLineThicknessF =  3.    ;-- marker line thickness

  pmres@gsMarkerColor = "royalblue3"
  pmres@gsMarkerIndex = 14   ;; 16实心，14半空心，4空心
  embb11 = gsn_add_polymarker(wks,plot(2), 1.9, region_trend(1,0)+0.43, pmres)
  embb12 = gsn_add_polymarker(wks,plot(2), 1.9, region_trend(1,0)+0.43, pmres)
  embb13 = gsn_add_polymarker(wks,plot(2), 1.9, region_trend(1,0)+0.43, pmres)
  pmres@gsMarkerColor = "royalblue3"
  pmres@gsMarkerIndex = 14   ;; 16实心，4空心
  embb21 = gsn_add_polymarker(wks,plot(2), 4.9, region_trend(1,1)+0.43, pmres)
  embb22 = gsn_add_polymarker(wks,plot(2), 4.9, region_trend(1,1)+0.43, pmres)
  embb23 = gsn_add_polymarker(wks,plot(2), 4.9, region_trend(1,1)+0.43, pmres)
  pmres@gsMarkerColor = "violetred2"
  pmres@gsMarkerIndex = 16   ;; 16实心，4空心
  embb31 = gsn_add_polymarker(wks,plot(2), 7.9, region_trend(1,2)+0.43, pmres)
  ; embb32 = gsn_add_polymarker(wks,plot(2), 7.9, region_trend(1,2)+0.43, pmres)
  ; embb33 = gsn_add_polymarker(wks,plot(2), 7.9, region_trend(1,2)+0.43, pmres)
  pmres@gsMarkerColor = "violetred2"
  pmres@gsMarkerIndex = 16   ;; 16实心，4空心
  embb41 = gsn_add_polymarker(wks,plot(2), 10.9, region_trend(1,3)+0.43, pmres)
  pmres@gsMarkerColor = "violetred2"
  pmres@gsMarkerIndex = 16   ;; 16实心，14半实心，4空心
  embb51 = gsn_add_polymarker(wks,plot(2), 13.9, region_trend(1,4)+0.43, pmres)
  ; embb52 = gsn_add_polymarker(wks,plot(2), 13.9, region_trend(1,4)+0.33, pmres)
  ; embb53 = gsn_add_polymarker(wks,plot(2), 13.9, region_trend(1,4)+0.33, pmres)
  pmres@gsMarkerColor = "royalblue3"
  pmres@gsMarkerIndex = 4   ;; 16实心，4空心
  emmm61 = gsn_add_polymarker(wks,plot(2), 16.9, region_trend(1,5)+0.43, pmres)
  emmm62 = gsn_add_polymarker(wks,plot(2), 16.9, region_trend(1,5)+0.43, pmres)
  emmm63 = gsn_add_polymarker(wks,plot(2), 16.9, region_trend(1,5)+0.43, pmres)
  pmres@gsMarkerColor = "violetred2"
  pmres@gsMarkerIndex = 14   ;; 16实心，4空心
  emcc71 = gsn_add_polymarker(wks,plot(2), 19.9, region_trend(1,6)+0.43, pmres)
  emcc72 = gsn_add_polymarker(wks,plot(2), 19.9, region_trend(1,6)+0.43, pmres)
  emcc73 = gsn_add_polymarker(wks,plot(2), 19.9, region_trend(1,6)+0.43, pmres)
  ; pmres@gsMarkerColor = "violetred2"
  ; pmres@gsMarkerIndex = 4   ;; 16实心，4空心
  ; emcc81 = gsn_add_polymarker(wks,plot(2), 22.9, region_trend(1,7)+0.33, pmres)
  ; emcc82 = gsn_add_polymarker(wks,plot(2), 22.9, region_trend(1,7)+0.33, pmres)
  ; emcc83 = gsn_add_polymarker(wks,plot(2), 22.9, region_trend(1,7)+0.33, pmres)




;;----------------- fig (c) markers  -------------------
  pmres@gsMarkerColor = "violetred2"
  pmres@gsMarkerIndex = 16   ;; 16实心，4空心
  emcc11 = gsn_add_polymarker(wks,plot(2), 17.5, 4.7+1.35, pmres)
  pmres@gsMarkerColor = "royalblue3"
  pmres@gsMarkerIndex = 16   ;; 16实心，4空心
  emcc12 = gsn_add_polymarker(wks,plot(2), 18.5, 4.7+1.35, pmres)

  pmres@gsMarkerColor = "violetred2"
  pmres@gsMarkerIndex = 14   ;; 14半空心
  emcc21 = gsn_add_polymarker(wks,plot(2), 17.5, 3.95+1.2, pmres)
  emcc22 = gsn_add_polymarker(wks,plot(2), 17.5, 3.95+1.2, pmres)
  emcc23 = gsn_add_polymarker(wks,plot(2), 17.5, 3.95+1.2, pmres)
  pmres@gsMarkerColor = "royalblue3"
  pmres@gsMarkerIndex = 14   ;; 14半空心
  emcc24 = gsn_add_polymarker(wks,plot(2), 18.5, 3.95+1.2, pmres)
  emcc25 = gsn_add_polymarker(wks,plot(2), 18.5, 3.95+1.2, pmres)
  emcc26 = gsn_add_polymarker(wks,plot(2), 18.5, 3.95+1.2, pmres)

  pmres@gsMarkerIndex = 4   ;; 16实心，4空心
  pmres@gsMarkerColor = "violetred2"
  emcc31 = gsn_add_polymarker(wks,plot(2), 17.5, 3.2+1.05, pmres)
  emcc32 = gsn_add_polymarker(wks,plot(2), 17.5, 3.2+1.05, pmres)
  emcc33 = gsn_add_polymarker(wks,plot(2), 17.5, 3.2+1.05, pmres)
  pmres@gsMarkerColor = "royalblue3"
  emcc34 = gsn_add_polymarker(wks,plot(2), 18.5, 3.2+1.05, pmres)
  emcc35 = gsn_add_polymarker(wks,plot(2), 18.5, 3.2+1.05, pmres)
  emcc36 = gsn_add_polymarker(wks,plot(2), 18.5, 3.2+1.05, pmres)

  txres@txFontHeightF = 0.0105
  emcc41 = gsn_add_text(wks,plot(2), "p<0.01", 20.35, 4.7+1.35, txres)
  emcc42 = gsn_add_text(wks,plot(2), "p<0.05", 20.35, 3.95+1.2, txres)
  emcc43 = gsn_add_text(wks,plot(2), "p<0.10", 20.35, 3.2+1.05, txres)






leres@trXMinF = 35.0
leres@trXMaxF = 90.0
leres@tmXBMode = "Explicit"
leres@tmXBValues = fspan(40.0, 80.0, 3)
; leres@tmXBLabels = fspan(40.0, 80.0, 3)
leres@tmXBLabels = (/"40%","60%","80%"/)
leres@tmXBMinorValues = fspan(30, 90, 7)
leres@xyLineColors       = (/"gray40","royalblue3","orange"/)
leres@xyLineThicknesses  = (/0.6,8.0,8.0/)
plot(3) = gsn_csm_xy(wks, ratio_lat(0:65), lat(0:65), leres)
leres@xyLineColors       = (/"green4","blue","orange"/)
leres@xyLineThicknesses  = (/2.2,8.0,8.0/)
plot2 = gsn_csm_xy(wks, ratio_lat_11moving(0:60), lat(0:60), leres)
overlay(plot(3), plot2)

delete([/x,y/])
x=(/70.0, 70.0/)
y=(/10.0, 85.0/)
   lnres   =    True
   lnres@gsLineColor="red"
   lnres@gsLineThicknessF=2.0
   lnres@gsLineDashPattern=2
   lnres@gsnFrame = False
   cvm1=gsn_add_polyline(wks,plot(3),x,y,lnres)









;;-----------------------------  (e-g) PDF  ------------------------------;;
leres@vpWidthF  = 0.34
leres@vpHeightF = 0.285
leres@trXMinF = -2.0
leres@trXMaxF = 8.0
leres@tmXBLabelFontHeightF = 0.0225
leres@tmYLLabelFontHeightF = 0.0225
leres@tmYRLabelFontHeightF = 0.0225
delete(leres@trYMinF)
delete(leres@trYMaxF)
leres@trYMinF = 0.4
leres@trYMaxF = 100.0
leres@tmYLDataBottomF = 0.1
leres@tmYLDataTopF = 50.0
leres@tmXBMode = "Explicit"
delete(leres@tmXBValues)
delete(leres@tmXBLabels)
delete(leres@tmXBMinorValues)
leres@tmXBValues = fspan(-2.0, 10.0, 7)
leres@tmXBLabels = fspan(-2.0, 10.0, 7)
leres@tmXBMinorValues = fspan(-2.0, 10.0, 13)
leres@xyLineColors       = (/"orange","blue","orange"/)
leres@xyLineThicknesses  = (/2.3,8.0,8.0/)
leres@xyDashPatterns     = (/16,0,0/)
leres@tmXBLabelsOn = False
; leres@tmYUseLeft = False
; leres@tmYLStyle = "Log"
leres@tmYLLabelsOn = True
leres@tmYLOn = True
leres@tmYRLabelsOn = False
leres@xyYStyle          = "Log"
delete(leres@tmYLValues)
delete(leres@tmYLLabels)
leres@tmYLMode          = "Explicit"
leres@tmYLValues = (/0.1,0.4,1.0,4.0,10.0,40.0,100.0/)
leres@tmYLLabels = (/0.1,0.4,1.0,4.0,10.0,40.0,100.0/)
leres@tmYLMinorOn = True
delete(leres@tmYLMinorValues)
leres@tmYLMinorValues = (/0.6,2.0,6.0,20.0,60.0/)
; leres@tmYLMinorValues = (/0.4,0.6,0.8,1.0,2.0,4.0,6.0,8.0,10.0,20.0,40.0,60.0,80.0,100.0/)
; leres@tmYRtyle = "Log"

plot(4) = gsn_csm_xy(wks, strength, pdf_first22, leres)
leres@xyLineColors       = (/"red","blue","orange"/)
leres@xyDashPatterns     = (/0,0,0/)
contour2(4) = gsn_csm_xy(wks, strength, pdf_last22, leres)
overlay(plot(4), contour2(4))



leres@xyLineColors       = (/"orange","blue","orange"/)
leres@xyLineThicknesses  = (/2.3,8.0,8.0/)
leres@xyDashPatterns     = (/16,0,0/)
plot(5) = gsn_csm_xy(wks, strength, pdfP_first22, leres)
leres@xyLineColors       = (/"red","blue","orange"/)
leres@xyDashPatterns     = (/0,0,0/)
contour2(5) = gsn_csm_xy(wks, strength, pdfP_last22, leres)
overlay(plot(5), contour2(5))




leres@xyLineColors       = (/"orange","blue","orange"/)
leres@xyLineThicknesses  = (/2.3,8.0,8.0/)
leres@xyDashPatterns     = (/16,0,0/)
leres@tmXBLabelsOn = True
plot(6) = gsn_csm_xy(wks, strength, pdfN_first22, leres)
leres@xyLineColors       = (/"red","blue","orange"/)
leres@xyDashPatterns     = (/0,0,0/)
contour2(6) = gsn_csm_xy(wks, strength, pdfN_last22, leres)
overlay(plot(6), contour2(6))

x=(/2.0, 2.0/)
y=(/0.4, 5.0/)
   lnres   =    True
   lnres@gsLineColor="blue"
   lnres@gsLineThicknessF=2.0
   lnres@gsLineDashPattern=2
   lnres@gsnFrame = False
   cvm4=gsn_add_polyline(wks,plot(4),x,y,lnres)
y=(/0.4, 5.8/)
   cvm5=gsn_add_polyline(wks,plot(5),x,y,lnres)
y=(/0.4, 2.3/)
   cvm6=gsn_add_polyline(wks,plot(6),x,y,lnres)


;;------------------------ figs (e-h) text -----------------------
  txres               = True
  txres@txFont = "helvetica"
  txres@txFontHeightF = 0.019
  txres@txFontColor = "orange"
  eum1 = gsn_add_text(wks,plot(5), "1979-2000", 1.35,0.7, txres)
  txres@txFontColor = "red"
  eum3 = gsn_add_text(wks,plot(5), "2001-2022", 5.65,4.0, txres)

  txres@txFontColor = "black"
  txres@txFontHeightF = 0.021
  euma1 = gsn_add_text(wks,plot(4), "All", 2.8,55.0, txres)
  euma2 = gsn_add_text(wks,plot(5), "Estimated > 2", 2.8,55.0,txres)
  euma3 = gsn_add_text(wks,plot(6), "Estimated < -2", 2.8,55.0, txres)



  txres@txFontHeightF = 0.021
  txres@txFontColor = "orange"
  eumc1 = gsn_add_text(wks,plot(4), "5.6%", 3.0,10.0, txres)
  txres@txFontColor = "black"
  eumc2 = gsn_add_text(wks,plot(4), "to", 4.7,10.0,txres)
  txres@txFontColor = "red"
  eumc3 = gsn_add_text(wks,plot(4), "19.7%", 6.5,10.0, txres)

  txres@txFontHeightF = 0.021
  txres@txFontColor = "orange"
  eumc4 = gsn_add_text(wks,plot(5), "5.4%", 2.8,14.0, txres)
  txres@txFontColor = "black"
  eumc5 = gsn_add_text(wks,plot(5), "to", 4.5,14.0,txres)
  txres@txFontColor = "red"
  eumc6 = gsn_add_text(wks,plot(5), "24.4%", 6.3,14.0, txres)

  txres@txFontColor = "orange"
  eumc7 = gsn_add_text(wks,plot(6), "5.2%", 3.0,8.0, txres)
  txres@txFontColor = "black"
  eumc8 = gsn_add_text(wks,plot(6), "to", 4.7,8.0,txres)
  txres@txFontColor = "red"
  eumc9 = gsn_add_text(wks,plot(6), "5.7%", 6.3,8.0, txres)




  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar = True
  pres@lbOrientation        = "Vertical"
  pres@pmLabelBarWidthF     = 0.05
  pres@pmLabelBarHeightF    = 0.16
  pres@pmLabelBarOrthogonalPosF = -0.0
  ; pres@pmLabelBarParallelPosF = -0.1
  pres@lbLabelFontHeightF  = 0.0115

  pres@gsnPanelTop  = 0.81
  pres@gsnPanelBottom = 0.61
  pres@gsnPanelLeft  = 0.0
  pres@gsnPanelRight = 0.85
  pres@amJust  ="BottomLeft"

  pres@gsnPanelFigureStringsFontHeightF = 0.012
  ; pres@gsnPanelFigureStringsPerimOn = False
  pres@pmLabelBarOrthogonalPosF = 0.0025
  pres@pmLabelBarParallelPosF = 0.01
  gsn_panel(wks,plot(0),(/1,1/),pres)

  pres@pmLabelBarWidthF     = 0.05
  pres@pmLabelBarHeightF    = 0.145
  pres@pmLabelBarOrthogonalPosF = 0.0025
  pres@pmLabelBarParallelPosF = 0.01
  pres@lbLabelFontHeightF  = 0.0115
  ; pres@pmLabelBarHeightF    = 0.154
  pres@gsnPanelTop  = 0.635
  pres@gsnPanelBottom = 0.435
  pres@gsnPanelLeft  = 0.0
  pres@gsnPanelRight = 0.85
  gsn_panel(wks,plot(1),(/1,1/),pres)


  pres@gsnPanelLabelBar = False
  pres@gsnPanelTop  = 0.44
  pres@gsnPanelBottom = 0.26
  pres@gsnPanelLeft  = 0.025
  pres@gsnPanelRight = 0.5
  pres@amJust  ="TopLeft"
  pres@pmLabelBarOrthogonalPosF = -0.0
  gsn_panel(wks,plot(2),(/1,1/),pres)

  pres@gsnPanelTop  = 0.44
  pres@gsnPanelBottom = 0.26
  pres@gsnPanelLeft  = 0.52
  pres@gsnPanelRight = 0.74
  pres@amJust  ="TopRight"
  pres@pmLabelBarOrthogonalPosF = -0.0
  gsn_panel(wks,plot(3),(/1,1/),pres)




  pres@gsnPanelTop  = 0.795
  pres@gsnPanelBottom = 0.630
  pres@gsnPanelLeft  = 0.77
  pres@gsnPanelRight = 1.0
  pres@amJust  ="TopRight"
  gsn_panel(wks,plot(4),(/1,1/),pres)

  pres@gsnPanelTop  = 0.6225
  pres@gsnPanelBottom = 0.4575
  pres@gsnPanelLeft  = 0.77
  pres@gsnPanelRight = 1.0
  gsn_panel(wks,plot(5),(/1,1/),pres)

  pres@gsnPanelTop  = 0.45
  pres@gsnPanelBottom = 0.265
  pres@gsnPanelLeft  = 0.77
  pres@gsnPanelRight = 1.0
  gsn_panel(wks,plot(6),(/1,1/),pres)



  txres               = True
  txres@txFontHeightF = 0.0155
  txres@txFont = "helvetica"
  txres@txFontColor = "black"
  txres@txAngleF = 0.0
  gsn_text_ndc(wks, "Zonally   Asymmetric   Trends   in   HWI",  0.373, 0.81, txres)
  gsn_text_ndc(wks, "PDF   in   HWI",  0.9, 0.81, txres)



;;------------------------- units ------------------------;; 
txres@txFontHeightF = 0.0105
txres@txFont = "helvetica"
txres@txFontColor = "black"
gsn_text_ndc(wks, "K/decade",  0.723, 0.63, txres)
gsn_text_ndc(wks, "K/decade",  0.723, 0.463, txres)

txres@txFontColor = "black"
txres@txFontHeightF = 0.0155
txres@txAngleF = 0.0
txres@txFont = "helvetica-bold"
gsn_text_ndc(wks, "a",  0.055, 0.80, txres)
gsn_text_ndc(wks, "b",  0.055, 0.628, txres)
gsn_text_ndc(wks, "c",  0.055, 0.445, txres)
gsn_text_ndc(wks, "d",  0.525, 0.445, txres)

gsn_text_ndc(wks, "e",  0.815, 0.80, txres)
gsn_text_ndc(wks, "f",  0.815, 0.627, txres)
gsn_text_ndc(wks, "g",  0.815, 0.452, txres)



frame(wks)
   delete(wks)  
  ; system("convert -geometry 1800x1800 -density 1200x1200 -trim " + pltName + "."+pltType + " " + pltName + ".png")       
;   system("rm " + pltName + "." + pltType + " -f") 





end






