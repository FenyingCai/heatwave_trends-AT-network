load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin
f00 = addfile("data/topo_adapt_era5_air2m.nc","r")
topo0 = f00->topo({20:60},{0:150})
lat = topo0&lat
lon = topo0&lon


;;------------  read ERA5 V200, Silk Road patterm  ------------;;
f0 = addfile("data/vwind.2004-01.daily.nc","r")
v0 = short2flt(f0->v(0,{200},{20:60},{0:150}))
dim_0 = dimsizes(v0)
delete(v0)
v = new((/44,92,dim_0(0),dim_0(1)/), float)

do i = 0,43
print(" i = "+i)
f1 = addfile("data/vwind."+(i+1979)+"-06.daily.nc","r")
v(i,0:29,:,:)  = short2flt(f1->v(:,{200},{20:60},{0:150}))
f2 = addfile("data/vwind."+(i+1979)+"-07.daily.nc","r")
v(i,30:60,:,:) = short2flt(f2->v(:,{200},{20:60},{0:150}))
f3 = addfile("data/vwind."+(i+1979)+"-08.daily.nc","r")
v(i,61:91,:,:) = short2flt(f3->v(:,{200},{20:60},{0:150}))
end do
v!0 = "year"
v!1 = "day"

v_era_JJA = dim_avg_n_Wrap(v, 1)
v_era_JJA!0 = "time"
delete(v)

;;----------------  EOF  ----------------;;
w = cos(0.01745329*v_era_JJA&latitude)
wp= v_era_JJA*conform(v_era_JJA, w, 1)
copy_VarCoords(v_era_JJA, wp)

x = wp(latitude|:,longitude|:,time|:)
eof0 = eofunc_Wrap(x, 1, False) 
eof_ts_yearly0 = eofunc_ts_Wrap(x, eof0, False)
eof_ts_yearly = dim_standardize_n_Wrap(eof_ts_yearly0(0,:), 1, 0)
printVarSummary(eof_ts_yearly)



eof_ts_dtrend = dtrend_n(eof_ts_yearly, False, 0)
srp_trend = regCoef_n(fspan(1,44,44), eof_ts_yearly, 0, 0)
srp_cor = escorc_n(fspan(1,44,44), eof_ts_yearly, 0, 0)
srp_t = rtest(srp_cor, 44, 0)



;;----------------------------------------------;;
;;---------------    read HWI    ---------------;;
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


;;-------------- residual trends，minus 31-degree moving means -------------;;
hwi_trend = regCoef_n(fspan(1,44,44), hwi_yearly, 0, 0) * 10.0
copy_VarCoords(hwi_yearly(0,:,:), hwi_trend)
hwi_trend = where(topo0.ge.0.0, hwi_trend, 1.0e20)
hwi_trend@_FillValue = 1.0e20
hwi_trend_avg = dim_avg_n_Wrap(hwi_trend, 1)
hwi_trend_31moving = new((/76,360/), float)
do ilat = 0,75
  if (ilat.le.15) then
    hwi_trend_31moving(ilat,:) = dim_avg(hwi_trend_avg(0:ilat+15))
  end if
  if (ilat.ge.60) then
    hwi_trend_31moving(ilat,:) = dim_avg(hwi_trend_avg(ilat-15:75))
  end if
  if (ilat.gt.15 .and. ilat.lt.60) then
    hwi_trend_31moving(ilat,:) = dim_avg(hwi_trend_avg(ilat-15:ilat+15))
  end if
end do
hwi_trend_ano = hwi_trend - hwi_trend_31moving
copy_VarCoords(hwi_yearly(0,:,:), hwi_trend_ano)






;;---------------------------------------------------------------------;;
;;-----------------     H500, UV500, HWI regressions    ---------------;;
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



;;--------------  Cor, Reg  --------------;;
h_cor = escorc_n(eof_ts_dtrend, h_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), h_cor)
u_cor = escorc_n(eof_ts_dtrend, u_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), u_cor)
v_cor = escorc_n(eof_ts_dtrend, v_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), v_cor)
hwi_cor = escorc_n(eof_ts_dtrend, hwi_detrend, 0, 0)
copy_VarCoords(h_JJA(0,:,:), hwi_cor)

h_reg = regCoef_n(eof_ts_dtrend, h_JJA, 0, 0) / 9.8
copy_VarCoords(h_JJA(0,:,:), h_reg)
u_reg = regCoef_n(eof_ts_dtrend, u_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), u_reg)
v_reg = regCoef_n(eof_ts_dtrend, v_JJA, 0, 0)
copy_VarCoords(h_JJA(0,:,:), v_reg)
hwi_reg = regCoef_n(eof_ts_dtrend, hwi_detrend, 0, 0)
copy_VarCoords(h_JJA(0,:,:), hwi_reg)

h_reg = mask(h_reg, h_cor.ge.0.2974.or.h_cor.le.-0.2974, True)
copy_VarCoords(h_JJA(0,:,:), h_reg)
u_reg = mask(u_reg, u_cor.ge.0.2974.or.u_cor.le.-0.2974.or.v_cor.ge.0.2974.or.v_cor.le.-0.2974, True)
copy_VarCoords(h_JJA(0,:,:), u_reg)
v_reg = mask(v_reg, u_cor.ge.0.2974.or.u_cor.le.-0.2974.or.v_cor.ge.0.2974.or.v_cor.le.-0.2974, True)
copy_VarCoords(h_JJA(0,:,:), v_reg)


printMinMax(hwi_reg, 1)


;;--------------  estimate HWI   --------------;;
hwi_estimated = srp_trend * hwi_reg * 10.0
copy_VarCoords(h_JJA(0,:,:), hwi_estimated)
print("spr " + srp_trend + "  " + srp_cor + " " + srp_t)
printMinMax(hwi_estimated, 1)


Contribution_hwi = hwi_estimated / hwi_trend_ano
copy_VarCoords(h_JJA(0,:,:), Contribution_hwi)

Contribution_hwi = where(abs(hwi_trend_ano).ge.0.01, Contribution_hwi, 0.0)  ;; remove margin
copy_VarCoords(h_JJA(0,:,:), Contribution_hwi)

Contribution_hwi_P = mask(Contribution_hwi, hwi_trend_ano.ge.0.0, True)
Contribution_hwi_N = mask(Contribution_hwi, hwi_trend_ano.le.0.0, True)
copy_VarCoords(h_JJA(0,:,:), Contribution_hwi_P)
copy_VarCoords(h_JJA(0,:,:), Contribution_hwi_N)


;;--------------  same sign  --------------;;
Contribution_hwi_all = where(abs(hwi_trend_ano).ge.0.0, 1.0, 0.0)
Contribution_hwi_sam_sign = where(abs(hwi_trend_ano).ge.0.0.and.Contribution_hwi.gt.0.0, 1.0, 0.0)
do ilat = 0,75
Contribution_hwi_all(ilat,:) = Contribution_hwi_all(ilat,:) * cos((ilat+10) /180.0 * 3.14159)
Contribution_hwi_sam_sign(ilat,:) = Contribution_hwi_sam_sign(ilat,:) * cos((ilat+10) /180.0 * 3.14159)
end do
ratio_lat = dim_sum(Contribution_hwi_sam_sign) / dim_sum(Contribution_hwi_all) * 100.0
print(ratio_lat)

ratio_lat_31moving = new((/76/), float)
do ilat = 0,75
  if (ilat.le.15) then
    ratio_lat_31moving(ilat) = dim_avg(ratio_lat(0:ilat+15))
  end if
  if (ilat.ge.60) then
    ratio_lat_31moving(ilat) = dim_avg(ratio_lat(ilat-15:75))
  end if
  if (ilat.gt.15 .and. ilat.lt.60) then
    ratio_lat_31moving(ilat) = dim_avg(ratio_lat(ilat-15:ilat+15))
  end if
end do
print(ratio_lat_31moving)





;*************************************** plots **************************************
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
res@tmXBLabelFontHeightF = 0.014
res@tmYLLabelFontHeightF = 0.014
res@tmXBTickSpacingF = 45.0
res@tmYLTickSpacingF = 20.0




res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.007
res@tmXBMajorLengthF = 0.007
res@tmYLMinorLengthF = 0.005
res@tmXBMinorLengthF = 0.005
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.008
res@tmYLMajorOutwardLengthF = 0.008
res@tmXBMinorOutwardLengthF = 0.005
res@tmYLMinorOutwardLengthF = 0.005
vcres = res
sres  = res


res@mpFillOn                    = False
res@mpOutlineOn                 = True
res@mpGeophysicalLineThicknessF = 0.8
res@mpGeophysicalLineColor      = "gray35";"gray20"
res@mpMinLatF = 10
res@mpMaxLatF = 85
res@mpMinLonF = 0
res@mpMaxLonF = 360
res@mpCenterLonF = 180-165-3
res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/) + 165+3
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
res@tmXBMinorOn = True
delete(res@tmXBMinorValues)
res@tmXBMinorValues = fspan(0, 360, 25)+3

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
res@cnFillPalette        = "precip2_17lev"



lires = True
lires@tmBorderThicknessF = 1.2
lires@gsnDraw        = False
lires@gsnFrame       = False
lires@tiXAxisOn      = False
lires@tiYAxisOn      = False
lires@tmXTLabelsOn   = False
lires@tmXTOn         = False

lires@trYMinF        = 2.5
lires@trYMaxF        = -3.5

lires@tmXBLabelFontHeightF = 0.020
lires@tmYLLabelFontHeightF = 0.020
lires@tmYRLabelFontHeightF = 0.020
lires@gsnLeftStringFontHeightF   = 0.022
lires@gsnRightStringFontHeightF  = 0.022
lires@tmYLLabelsOn   = True
lires@tmYLOn         = True

lires@tmYLMode = "Explicit"
lires@tmYLValues = fspan(-3, 3, 7)
lires@tmYLLabels = fspan(-3, 3, 7)
lires@tmYLMinorOn = True
lires@tmYLMinorValues = fspan(-3.5, 3, 14)

lires@lgPerimOn = False
lires@lgLabelFontHeightF = 0.018
lires@vpWidthF  = 0.82
lires@vpHeightF = 0.235
lires@vpXF      = 0.06
lires@vpYF      = 0.75

lires@gsnLeftString  = ""
lires@gsnRightString = ""
lires@xyMarkLineModes = (/"MarkLines","Lines","Lines"/)
lires@xyMarkers = (/1,1,1/)
lires@xyMarkerSizes = (/0.015,0.001,0.001/)
lires@xyMarkerColors= (/"red","black","black"/)
lires@xyMarkerThicknesses = (/4.0,2.0,1.0/)
lires@xyLineThicknesses  = (/3.5,8.0,8.0/)
lires@xyLineColors       = (/"red","blue","orange"/)
lires@xyDashPatterns     = (/0,4,2/)

lires@tmYLMajorThicknessF = 1.5
lires@tmYRMajorThicknessF = 1.5
lires@tmXBMajorThicknessF = 1.5
lires@tmYLMinorThicknessF = 1.0
lires@tmYRMinorThicknessF = 1.0
lires@tmXBMinorThicknessF = 1.0
lires@tmYLMajorLengthF = 0.014
lires@tmYRMajorLengthF = 0.014
lires@tmXBMajorLengthF = 0.014
lires@tmYLMinorLengthF = 0.008
lires@tmYRMinorLengthF = 0.008
lires@tmXBMinorLengthF = 0.008
lires@tmXBLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.013
lires@tmYLMajorOutwardLengthF = 0.013
lires@tmYRMajorOutwardLengthF = 0.013
lires@tmXBMinorOutwardLengthF = 0.007
lires@tmYLMinorOutwardLengthF = 0.007
lires@tmYRMinorOutwardLengthF = 0.007




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
leres@tmYLMinorValues = fspan(0, 80, 9);leres@tmYLMode = "Explicit"


leres@xyMarkLineModes = (/"Lines","Lines","Lines"/)

leres@lgPerimOn = False
leres@lgLabelFontHeightF = 0.018
leres@vpWidthF  = 0.42
leres@vpHeightF = 0.265
leres@vpXF      = 0.06
leres@vpYF      = 0.75

leres@gsnLeftString  = ""
leres@gsnRightString = ""
leres@xyMarkLineModes = (/"Lines","Lines","Lines"/)
leres@xyLineThicknesses  = (/5.0,8.0,8.0/)
leres@xyLineColors       = (/"red","blue","orange"/)
leres@xyDashPatterns     = (/0,4,2/)

leres@tmYLMajorThicknessF = 1.5
leres@tmYRMajorThicknessF = 1.5
leres@tmXBMajorThicknessF = 1.5
leres@tmYLMinorThicknessF = 1.0
leres@tmYRMinorThicknessF = 1.0
leres@tmXBMinorThicknessF = 1.0
leres@tmYLMajorLengthF = 0.014
leres@tmYRMajorLengthF = 0.014
leres@tmXBMajorLengthF = 0.014
leres@tmYLMinorLengthF = 0.008
leres@tmYRMinorLengthF = 0.008
leres@tmXBMinorLengthF = 0.008
leres@tmXBLabelDeltaF = -0.25
leres@tmYLLabelDeltaF = -0.25
leres@tmXBMajorOutwardLengthF = 0.013
leres@tmYLMajorOutwardLengthF = 0.013
leres@tmYRMajorOutwardLengthF = 0.013
leres@tmXBMinorOutwardLengthF = 0.007
leres@tmYLMinorOutwardLengthF = 0.007
leres@tmYRMinorOutwardLengthF = 0.007




;vcres@mpFillDrawOrder         = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.255
;vcres@vcRefAnnoSide           = "TopRight"
vcres@vcGlyphStyle            = "LineArrow"
vcres@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
vcres@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref
vcres@vcMinDistanceF          = 0.015             ; thin out vectors
vcres@vcLineArrowColor        = "black"           ; change vector color
vcres@vcRefAnnoOn             = True
; vcres@vcRefAnnoString1        = "5"
; vcres@vcRefMagnitudeF         = 5
vcres@vcRefLengthF            = 0.015            ; ref vec length
vcres@vcRefAnnoFontHeightF    = 0.01   ;参考箭头字体大小
vcres@vcRefAnnoString1On      = True
vcres@vcRefAnnoString2On      = False
vcres@vcLineArrowThicknessF   = 1.6            ; make vectors larger
vcres@vcVectorDrawOrder = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -1.0 ;;正向下
vcres@vcRefAnnoParallelPosF   = 1.0  ;;正向右


sres@cnFillOn  = False
sres@cnLinesOn = True
sres@cnLineColor = "gray40"
sres@cnLevelSelectionMode = "ExplicitLevels"
sres@cnLevels         = (/0.01,0.05/)
sres@cnLineColor     = "white"
sres@cnLineThicknessF = 4.0
sres@cnLineDashPattern = 0
sres@cnLineLabelsOn   = True
sres@cnLineLabelPerimOn   = False
sres@cnLineLabelBackgroundColor = -1
sres@cnLineLabelFontHeightF = 0.01
sres@cnLineLabelFontColor = "white"
sres@cnConstFLabelOn  = False
sres@cnInfoLabelOn    = False 
sres@cnMissingValFillColor= -1
sres@cnLineLabelInterval  = 1
sres@cnLabelMasking       = True
sres@cnLineLabelDensityF  = 1.8
sres@lbLabelBarOn = False
sres@cnFillDrawOrder = "PostDraw"
; sres@cnLineDrawOrder = "PostDraw"



;*************************dottings***************************
ores                 = True            
ores@gsnDraw         = False               
ores@gsnFrame        = False                 
ores@cnLineLabelsOn  = False               
ores@cnLinesOn       = False                 
ores@cnInfoLabelOn   = False                                                                                                        
ores@cnFillOn        = True                               
ores@lbLabelBarOn    = False                                                                                                            
ores@cnLevelSelectionMode = "ExplicitLevels"                                                                                                                                                                                
ores@cnLevels        = (/0.05/)   ;;90%
ores@cnMonoFillPattern    = False            ; want multiple patterns                                                               
ores@cnFillPatterns       = (/17,-1/)     ; the patterns                                                                         
ores@cnMonoFillScale      = False            ; want different densities                                                                        
ores@cnMonoFillColor      =True                                                                                                     
ores@cnFillDotSizeF       = 0.003    
ores@cnFillDrawOrder ="postDraw"
ores@cnFillColor = "white"
ores@cnFillScales         = (/0.6,0.6/) 




;**********************plot**************************
plot = new(12, "graphic")
vector = new(12, "graphic")
vector2 = new(12, "graphic")
contour = new(12, "graphic")
contour1 = new(12, "graphic")
contour2 = new(12, "graphic")
contour3 = new(12, "graphic")
contour4 = new(12, "graphic")
contour5 = new(12, "graphic")
contour6 = new(12, "graphic")
contour7 = new(12, "graphic")
contour8 = new(12, "graphic")
contour9 = new(12, "graphic")
contour10 = new(12, "graphic")

pltType = "eps"
pltName = "FigureS14_Silk_road_pattern"
wks = gsn_open_wks(pltType, pltName)


;;----------------------- 回归 H500， UV500 -----------------------
res@cnFillPalette = "MPL_PiYG"
res@cnLevels      = (/-20,-2,2,20.0/)
res@cnFillColors  = (/95,85,-1,45,30/)
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 1.0
res@tmXBLabelsOn  = False
plot(0) = gsn_csm_contour_map_ce(wks, h_reg, res)


vcres@vcRefAnnoString1        = "1"
vcres@vcRefMagnitudeF         = 1
vcres@vcLineArrowColor        = "royalblue3"
vcres@vcRefAnnoFontColor      = "royalblue3"
vcres@vcRefAnnoArrowFillColor = "royalblue3"
vcres@vcRefAnnoArrowEdgeColor = "royalblue3"
vcres@vcRefAnnoArrowUseVecColor = True
vcres@vcLineArrowColor        = "royalblue3"
vector(0)=gsn_csm_vector(wks, u_reg, v_reg, vcres)
overlay(plot(0), vector(0))



;;----------------------- reg hwi -----------------------
res@cnFillPalette = "BlueWhiteOrangeRed"
res@cnLevels      = (/-5.0,-1.0,1.0,5.0/)
res@cnFillColors  = (/75,100,-1,155,175/)
res@tmXBLabelsOn  = False
res@mpGeophysicalLineThicknessF = 0.7
plot(1) = gsn_csm_contour_map_ce(wks, hwi_reg, res)

ores@cnFillColor = "red4"
ores@cnLevels        = (/0.2973/)   ;; 95%                                                             
ores@cnFillPatterns       = (/-1,11/) 
contour2(1) = gsn_csm_contour(wks, hwi_cor, ores)    
overlay(plot(1), contour2(1))
contour3(1) = gsn_csm_contour(wks, hwi_cor, ores)    
overlay(plot(1), contour3(1))
ores@cnFillColor = "blue"
ores@cnLevels        = (/-0.2973/)   ;; 95%                                                             
ores@cnFillPatterns       = (/11,-1/) 
contour4(1) = gsn_csm_contour(wks, hwi_cor, ores)    
overlay(plot(1), contour4(1))
contour5(1) = gsn_csm_contour(wks, hwi_cor, ores)    
overlay(plot(1), contour5(1))





res@cnFillPalette = "BlueWhiteOrangeRed"
res@cnLevels      = (/-3.0,-0.5,0.5,3.0/)
res@cnFillColors  = (/75,100,-1,155,175/)
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 0.7
res@tmXBLabelsOn  = True
plot(2) = gsn_csm_contour_map_ce(wks, hwi_estimated, res)




ores@cnLevels    = (/0.0/)
ores@cnFillColor = "violetred3"
ores@cnFillPatterns       = (/-1,11/) 
contour3(2) = gsn_csm_contour(wks, Contribution_hwi_P, ores)    
overlay(plot(2), contour3(2))
contour8(2) = gsn_csm_contour(wks, Contribution_hwi_P, ores)    
overlay(plot(2), contour8(2))
ores@cnFillColor = "royalblue4"
contour4(2) = gsn_csm_contour(wks, Contribution_hwi_N, ores)    
overlay(plot(2), contour4(2))
contour9(2) = gsn_csm_contour(wks, Contribution_hwi_N, ores)    
overlay(plot(2), contour9(2))



lires@trXMinF = 1979
lires@trXMaxF = 2022
lires@tmXBMode = "Explicit"
lires@tmXBValues = fspan(1980, 2020, 5)
lires@tmXBLabels = fspan(1980, 2020, 5)
lires@tmXBMinorValues = fspan(1980, 2020, 9)
lires@xyLineColors       = (/"orange","blue","orange"/)
lires@xyLineThicknesses  = (/2.0,8.0,8.0/)
lires@xyDashPatterns  = (/14,14,14/)
plot(3) = gsn_csm_xy(wks, fspan(1979,2022,44), eof_ts_yearly, lires)
lires@xyLineColors       = (/"gray30","blue","orange"/)
lires@xyLineThicknesses  = (/0.8,8.0,8.0/)
lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)
plot3 = gsn_csm_xy(wks, (/1978,2023/), (/0,0/), lires)
overlay(plot(3), plot3)

lires@xyLineColors       = (/"red4","blue","orange"/)
lires@xyLineThicknesses  = (/2.8,8.0,8.0/)
lires@xyDashPatterns  = (/0,14,14/)
plot5 = gsn_csm_xy(wks, (/2000.5-19,2000.5+19/), (/-19*srp_trend, 19*srp_trend/), lires)
overlay(plot(3), plot5)




leres@trXMinF = 25.0
leres@trXMaxF = 70.0
leres@tmXBMode = "Explicit"
leres@tmXBValues = fspan(20.0, 80.0, 7)
; leres@tmXBLabels = fspan(20.0, 80.0, 7)
leres@tmXBLabels = (/"20%","30%","40%","50%","60%","70%","80%"/)
leres@tmXBMinorValues = fspan(20, 80, 13)
leres@xyLineColors       = (/"gray40","blue","orange"/)
leres@xyLineThicknesses  = (/0.8,8.0,8.0/)
plot(4) = gsn_csm_xy(wks, ratio_lat(0:65), lat(0:65), leres)
leres@xyLineColors       = (/"green4","blue","orange"/)
leres@xyLineThicknesses  = (/2.2,8.0,8.0/)
plot2 = gsn_csm_xy(wks, ratio_lat_31moving(0:65), lat(0:65), leres)
overlay(plot(4), plot2)

x2=(/50.0,50.0/)
y2=(/10, 85.0/)
   lnres   =    True
   lnres@gsLineColor="red"
   lnres@gsLineThicknessF=1.8
   lnres@gsLineDashPattern=2
   lnres@gsnFrame = False
   cvm1=gsn_add_polyline(wks,plot(4),x2,y2,lnres)
  




  txres               = True
  txres@txFontHeightF = 0.016
  txres@txFont = "helvetica"
  txres@txFontHeightF = 0.012

  txres@txFontColor = "black"
  eur0 = gsn_add_text(wks,plot(0), "Regressed (H500)", 158,17.0, txres)
  eur1 = gsn_add_text(wks,plot(1), "Regressed (HWI)", 158,17.0, txres)
  eur2 = gsn_add_text(wks,plot(2), "Estimated (HWI)", 158,17.0, txres)
  txres@txFontHeightF = 0.019
  eur3 = gsn_add_text(wks,plot(3), "slope: 0.043", 2016,-1.6, txres)
  eur4 = gsn_add_text(wks,plot(3), "(p<0.001)", 2016,-2.6, txres)


  pres =True
  pres@gsnPanelLabelBar = True
  pres@gsnFrame         = False   
  pres@lbOrientation    = "Vertical"
  pres@lbLabelFontHeightF   = 0.0113
  pres@pmLabelBarWidthF     = 0.04
  pres@pmLabelBarHeightF    = 0.13
  ; pres@pmLabelBarOrthogonalPosF = -0.005
  pres@amJust  ="BottomLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.012

  pres@gsnPanelLeft = 0.05
  pres@gsnPanelRight = 0.95
  pres@gsnPanelTop = 0.85
  pres@gsnPanelBottom = 0.65
  pres@gsnPanelFigureStrings = (/"(a)","(c)"/)
  gsn_panel(wks,plot(0),(/1,1/),pres)
  pres@gsnPanelTop = 0.69
  pres@gsnPanelBottom = 0.49
  pres@gsnPanelFigureStrings = (/"(b)","(c)"/)
  gsn_panel(wks,plot(1),(/1,1/),pres)
  pres@gsnPanelTop = 0.52
  pres@gsnPanelBottom = 0.32
  pres@gsnPanelFigureStrings = (/"(c)","(c)"/)
  gsn_panel(wks,plot(2),(/1,1/),pres)

  pres@gsnPanelLeft = 0.065
  pres@gsnPanelRight = 0.565
  pres@gsnPanelTop = 0.35
  pres@gsnPanelBottom = 0.12
  pres@amJust  ="TopLeft"
  pres@gsnPanelFigureStrings = (/"(d)","(c)"/)
  gsn_panel(wks,plot(3),(/1,1/),pres)

  pres@gsnPanelLeft = 0.58
  pres@gsnPanelRight = 0.83
  pres@gsnPanelTop = 0.35
  pres@gsnPanelBottom = 0.12
  pres@amJust  ="TopRight"
  pres@gsnPanelFigureStrings = (/"(e)","(c)"/)
  gsn_panel(wks,plot(4),(/1,1/),pres)

  txres               = True
  txres@txFontHeightF = 0.016
  txres@txFont = "helvetica"
  gsn_text_ndc(wks, "Circum-global   teleconnection   (CGT)   pattern",  0.44, 0.85, txres)



frame(wks)
   delete(wks)  
   system("convert -geometry 1600x1600 -density 800x800 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
end




