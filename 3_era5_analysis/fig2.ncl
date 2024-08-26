load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin
;;**********************************  read lat, lon  ************************************;;
f00 = addfile("data/topo_adapt_era5_air2m.nc","r")
topo0 = f00->topo({10:85},{0:360})
lat = topo0&lat
lon = topo0&lon
dim_0 = dimsizes(topo0)




;;----------------------------------------------------------;;
;;--------------    read heatwave intensity    -------------;;
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

;;-------------- exclude ocean, remove trend -------------;;
topo2 = conform(hwi_yearly, topo0, (/1,2/))
hwi_yearly = where(topo2.ge.0.0, hwi_yearly, 1.0e20)
hwi_yearly@_FillValue = 1.0e20
copy_VarCoords(hwi_yearly2, hwi_yearly)

hwi_detrend = dtrend_msg_n(fspan(1, 44, 44), hwi_yearly, False, False, 0)
copy_VarCoords(hwi_yearly, hwi_detrend)







;;-----------------------------------------------------------;;
;;--------------------   read AT index   --------------------;;
f1 = addfile("data/Teleconnection_Tmax90th_H500_removeZonalTrend.nc","r")
degree11 = f1->Teleconnection(:,0:,:)  ;; 10-85N
f2 = addfile("data/Teleconnection_Tmax90th_H500_removeZonalTrend_ConcurrentDays.nc","r")
degree12 = f2->Teleconnection(:,0:,:)  ;; 10-85N
Teleconnection = degree11 + degree12
delete([/degree11,degree12/])

Teleconnection@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(Teleconnection))) then
      replace_ieeenan (Teleconnection, value, 0)
      Teleconnection@_FillValue = value
end if

;;------------- exclude ocean -------------;;
topo3 = conform(Teleconnection, topo0, (/1,2/))
Teleconnection = where(topo3.ge.0.0, Teleconnection, 1.0e20)
Teleconnection@_FillValue = 1.0e20

Teleconnection!1 = "lat"
Teleconnection!2 = "lon"
Teleconnection&lat = lat
Teleconnection&lon = lon


print(" ----------- Teleconnection ---------- ")
printMinMax(Teleconnection, 0)

Teleconnection_detrend = dtrend_msg_n(fspan(1, 44, 44), Teleconnection, False, False, 0)
copy_VarCoords(Teleconnection, Teleconnection_detrend)





;;----------------------------------------------------------------------;;
;;----------------------- explained variances R2  ----------------------;;
cor = escorc_n(hwi_yearly, Teleconnection, 0, 0)
R2 = cor * cor * 100.0
copy_VarCoords(Teleconnection(0,:,:), cor)
copy_VarCoords(Teleconnection(0,:,:), R2)

cor_detrend = escorc_n(hwi_detrend, Teleconnection_detrend, 0, 0)
R2_detrend = cor_detrend * cor_detrend * 100.0
copy_VarCoords(Teleconnection(0,:,:), cor_detrend)
copy_VarCoords(Teleconnection(0,:,:), R2_detrend)



R2_lat = dim_median_n(R2, 1)
cor_detrend_lat = dim_median_n(cor_detrend, 1)
R2_detrend_lat = dim_median_n(R2_detrend, 1)


fn = addfile("data/R2_local_high.nc","r")
R2_detrend_lat_High = fn->R2_detrend_lat_High  ;; 0-85N



;;---------------------------------------------------------------;;
;;----------------------- long-term trend  ----------------------;;
trend = regCoef_n(fspan(1,44,44), Teleconnection, 0,0) * 1.0e2
copy_VarCoords(Teleconnection(0,:,:), trend)
printMinMax(trend, 1)

trend_significantCor = escorc_n(fspan(1,44,44), Teleconnection, 0,0)
copy_VarCoords(Teleconnection(0,:,:), trend_significantCor)







;;********************************** read lat, lon ************************************
f00 = addfile("data/topo_adapt_era5_air2m.nc","r")
topoC = f00->topo({0:90},{0:360})
;;-----------------------------------------------------------;;
;;--------------------  eastern Europe  --------------------;;
f1 = addfile("data/era5_EE_composite_T2m_UVH500.nc","r")   ;;; ???
t_diff = doubletofloat(f1->t_diff)
h_diff = doubletofloat(f1->h_diff)
u_diff = doubletofloat(f1->u_diff)
v_diff = doubletofloat(f1->v_diff)
t_diff@_FillValue = 1.0e20
h_diff@_FillValue = 1.0e20
u_diff@_FillValue = 1.0e20
v_diff@_FillValue = 1.0e20
if (any(isnan_ieee(t_diff))) then
      replace_ieeenan (t_diff, value, 0)
      t_diff@_FillValue = value
end if
if (any(isnan_ieee(h_diff))) then
      replace_ieeenan (h_diff, value, 0)
      h_diff@_FillValue = value
end if
if (any(isnan_ieee(u_diff))) then
      replace_ieeenan (u_diff, value, 0)
      u_diff@_FillValue = value
end if
if (any(isnan_ieee(v_diff))) then
      replace_ieeenan (v_diff, value, 0)
      v_diff@_FillValue = value
end if

;;------------- exclude ocean -------------;;
t_diff = where(topoC.ge.0.0, t_diff, 1.0e20)
t_diff@_FillValue = 1.0e20
h_diff = h_diff/9.8
copy_VarCoords(topoC, t_diff)
copy_VarCoords(topoC, h_diff)
copy_VarCoords(topoC, u_diff)
copy_VarCoords(topoC, v_diff)
print("-------------- T,H,U,V --------------")
printMinMax(t_diff, 1)
printMinMax(h_diff, 1)


;;--------------------  significant (eastern Europe)  --------------------;;
f1 = addfile("data/ttest_EE_T2m_UVH500.nc","r")
t_t = doubletofloat(f1->t_t)
u_t = doubletofloat(f1->u_t)
v_t = doubletofloat(f1->v_t)
h_t = doubletofloat(f1->h_t)

u_diff = mask(u_diff, u_t.le.0.05.or.v_t.le.0.05, True)
v_diff = mask(v_diff, u_t.le.0.05.or.v_t.le.0.05, True)

t_t_Posit = mask(t_t, t_diff.ge.0.0, True)
t_t_Negat = mask(t_t, t_diff.le.0.0, True)
h_t_Posit = mask(h_t, h_diff.ge.0.0, True)
h_t_Negat = mask(h_t, h_diff.le.0.0, True)

;;--------------------  wave activity flux (eastern Europe)  --------------------;;
f1 = addfile("data/era5_EE_WAF.nc","r")
Fx = doubletofloat(f1->Fx({500},:,:))
Fy = doubletofloat(f1->Fy({500},:,:))
Fx@_FillValue = 1.0e20
if (any(isnan_ieee(Fx))) then
      replace_ieeenan (Fx, value, 0)
      Fx@_FillValue = value
end if
Fy@_FillValue = 1.0e20
if (any(isnan_ieee(Fy))) then
      replace_ieeenan (Fy, value, 0)
      Fy@_FillValue = value
end if
printMinMax(Fx, 1)
printMinMax(Fy, 1)

Fx = mask(Fx, u_t.le.0.05.or.v_t.le.0.05, True)
Fy = mask(Fy, u_t.le.0.05.or.v_t.le.0.05, True)

Fx = mask(Fx, Fx.ge.0.5.or.Fy.ge.0.5, True)
Fy = mask(Fy, Fx.ge.0.5.or.Fy.ge.0.5, True)
copy_VarCoords(topoC, Fx)
copy_VarCoords(topoC, Fy)



;;-----------------------------------------------------------;;
;;--------------------  western Europe  --------------------;;
f1 = addfile("data/era5_WE_composite_T2m_UVH500.nc","r")
t2_diff = doubletofloat(f1->t_diff)
h2_diff = doubletofloat(f1->h_diff)
u2_diff = doubletofloat(f1->u_diff)
v2_diff = doubletofloat(f1->v_diff)
t2_diff@_FillValue = 1.0e20
h2_diff@_FillValue = 1.0e20
u2_diff@_FillValue = 1.0e20
v2_diff@_FillValue = 1.0e20
if (any(isnan_ieee(t2_diff))) then
      replace_ieeenan (t2_diff, value, 0)
      t2_diff@_FillValue = value
end if
if (any(isnan_ieee(h2_diff))) then
      replace_ieeenan (h2_diff, value, 0)
      h2_diff@_FillValue = value
end if
if (any(isnan_ieee(u2_diff))) then
      replace_ieeenan (u2_diff, value, 0)
      u2_diff@_FillValue = value
end if
if (any(isnan_ieee(v2_diff))) then
      replace_ieeenan (v2_diff, value, 0)
      v2_diff@_FillValue = value
end if

;;------------- exclude ocean -------------;;
t2_diff = where(topoC.ge.0.0, t2_diff, 1.0e20)
t2_diff@_FillValue = 1.0e20
h2_diff = h2_diff/9.8
copy_VarCoords(topoC, t2_diff)
copy_VarCoords(topoC, h2_diff)
copy_VarCoords(topoC, u2_diff)
copy_VarCoords(topoC, v2_diff)
print("-------------- T,H,U,V --------------")
printMinMax(t2_diff, 1)
printMinMax(h2_diff, 1)


;;--------------------  significant (western Europe)  --------------------;;
f1 = addfile("data/ttest_WE_T2m_UVH500.nc","r")
t2_t = doubletofloat(f1->t_t)
u2_t = doubletofloat(f1->u_t)
v2_t = doubletofloat(f1->v_t)
h2_t = doubletofloat(f1->h_t)

u2_diff = mask(u2_diff, u2_t.le.0.05.or.v2_t.le.0.05, True)
v2_diff = mask(v2_diff, u2_t.le.0.05.or.v2_t.le.0.05, True)

t2_t_Posit = mask(t2_t, t2_diff.ge.0.0, True)
t2_t_Negat = mask(t2_t, t2_diff.le.0.0, True)
h2_t_Posit = mask(h2_t, h2_diff.ge.0.0, True)
h2_t_Negat = mask(h2_t, h2_diff.le.0.0, True)


;;--------------------  wave activity flux (western Europe)  --------------------;;
f1 = addfile("data/era5_WE_WAF.nc","r")
Fx2 = doubletofloat(f1->Fx({500},:,:))
Fy2 = doubletofloat(f1->Fy({500},:,:))
Fx2@_FillValue = 1.0e20
if (any(isnan_ieee(Fx2))) then
      replace_ieeenan (Fx2, value, 0)
      Fx2@_FillValue = value
end if
Fy2@_FillValue = 1.0e20
if (any(isnan_ieee(Fy2))) then
      replace_ieeenan (Fy2, value, 0)
      Fy2@_FillValue = value
end if
printMinMax(Fx2, 1)
printMinMax(Fy2, 1)

Fx2 = mask(Fx2, u2_t.le.0.05.or.v2_t.le.0.05, True)
Fy2 = mask(Fy2, u2_t.le.0.05.or.v2_t.le.0.05, True)

Fx2 = mask(Fx2, Fx2.ge.0.5.or.Fy2.ge.0.5, True)
Fy2 = mask(Fy2, Fx2.ge.0.5.or.Fy2.ge.0.5, True)
copy_VarCoords(topoC, Fx2)
copy_VarCoords(topoC, Fy2)




;;-------------- geometric center for eastern and western Europe -------------;;
grids65 = new((/dim_0(0),dim_0(1)/), float)
lat_2D = conform(grids65, lat, 0)
lon_2D = conform(grids65, lon, 1)
lon_2D = where(lon_2D.ge.180.0, (lon_2D-360.0), lon_2D)
do ilat = 0,dim_0(0)-1
  grids65(ilat,:) = cos((ilat+10) / 180.0 * 3.14159)
end do
grids65 = where(topo0.ge.0.0, grids65, 1.0e20)
grids65@_FillValue = 1.0e20
lat_2D = lat_2D * grids65
lon_2D = lon_2D * grids65


lat_WE_index = new(36,integer)
lat_WE_index(0:9) = ispan(350,359,1)
lat_WE_index(10:35) = ispan(0,25,1)
lat_WE_avg = dim_sum(dim_sum(lat_2D(35-10:55-10, lat_WE_index))) / dim_sum(dim_sum(grids65(35-10:55-10, lat_WE_index)))
lon_WE_avg = dim_sum(dim_sum(lon_2D(35-10:55-10, lat_WE_index))) / dim_sum(dim_sum(grids65(35-10:55-10, lat_WE_index)))
print("WE, center, "+ lat_WE_avg + " " + lon_WE_avg)


lat_EE_avg = (dim_sum(dim_sum(lat_2D(45-10:65-10, 25:55)))) / dim_sum(dim_sum(grids65(45-10:65-10, 25:55)))
lon_EE_avg = (dim_sum(dim_sum(lon_2D(45-10:65-10, 25:55)))) / dim_sum(dim_sum(grids65(45-10:65-10, 25:55)))
print("EE, center, "+ lat_EE_avg + " " + lon_EE_avg)




;-------------- distance to EE and WE -------------;;
f0 = addfile("data/distance_i_j_land.nc","r")
distance = f0->distance(55,40,0:90,:)
copy_VarCoords(topoC, distance)

f0 = addfile("data/distance_i_j.nc","r")
distance2 = f0->distance(46,9,0:90,:)
distance2 = where(topoC.le.1.0, distance2, 1.0e20)
distance2@_FillValue = 1.0e20
copy_VarCoords(topoC, distance2)










;;-----------------------------------------------------------;;
;;--------------------   read  degree   --------------------;;
f0 = addfile("data/Degree_High_EasternEurope_WesternEurope.nc","r")
degree_EE_High = (f0->degree_EE)/1000.0
degree_WE_High = (f0->degree_WE)/1000.0

f1 = addfile("data/Degree_Low_EasternEurope_WesternEurope.nc","r")
degree_EE_Low = (f1->degree_EE)/1000.0
degree_WE_Low = (f1->degree_WE)/1000.0


degree_EE_High!0 = "lat"
degree_EE_High!1 = "lon"
degree_EE_High&lat = lat
degree_EE_High&lon = lon
copy_VarCoords(degree_EE_High, degree_WE_High)
copy_VarCoords(degree_EE_High, degree_EE_Low)
copy_VarCoords(degree_EE_High, degree_WE_Low)


printMinMax(degree_EE_High, 1)
printMinMax(degree_WE_High, 1)
printMinMax(degree_EE_Low, 1)
printMinMax(degree_WE_Low, 1)







;************************************************** plots **********************************************
res = True
; res@gsnPaperWidth  = 8.27  ; A4 page
; res@gsnPaperHeight = 11.69 ; A4 page
; res@gsnMaximize = True

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
res@tmXBLabelFontHeightF = 0.016
res@tmYLLabelFontHeightF = 0.016
res@tmXBTickSpacingF = 45.0
res@tmYLTickSpacingF = 20.0




res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.0012
res@tmXBMajorLengthF = 0.0012
res@tmYLMinorLengthF = 0.0012
res@tmXBMinorLengthF = 0.009
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.012
res@tmYLMajorOutwardLengthF = 0.012
res@tmXBMinorOutwardLengthF = 0.009
res@tmYLMinorOutwardLengthF = 0.009





res@mpFillOn                    = False
res@mpOutlineOn                 = True
res@mpGeophysicalLineThicknessF = 0.8
res@mpGeophysicalLineColor      = "gray30";"gray20"
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
vcres@vcMinDistanceF          = 0.03             ; thin out vectors
vcres@vcLineArrowColor        = "black"           ; change vector color
vcres@vcRefAnnoOn             = True
vcres@vcRefLengthF            = 0.022            ; ref vec length
vcres@vcRefAnnoFontHeightF    = 0.017   ;参考箭头字体大小
vcres@vcRefAnnoString1On      = True
vcres@vcRefAnnoString2On      = False
vcres@vcLineArrowThicknessF   = 1.8            ; make vectors larger
vcres@vcVectorDrawOrder = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.295 ;;正向下
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
sres@cnLineColor = "red"
sres@cnLevelSelectionMode = "ExplicitLevels"
sres@cnLevels         = (/2700.0/)
sres@cnFillColors     = (/"red","gray20"/)
sres@cnMissingValFillColor= -1
sres@cnLineThicknessF = 1.8
sres@cnLineDashPattern = 14


sres2@cnLevelSelectionMode = "ExplicitLevels"
sres2@cnLevels         = (/2/)
sres2@cnLineColor = "blue"
sres2@cnLineDashPattern = 14
sres2@cnLineThicknessF = 2.5







lires = True
lires@tmBorderThicknessF = 1.2
lires@gsnDraw        = False
lires@gsnFrame       = False
lires@tiXAxisOn      = False
lires@tiYAxisOn      = False
lires@tmXTLabelsOn   = False
lires@tmXTOn         = False

lires@trYMinF        = 10
lires@trYMaxF        = 85

lires@tmXBLabelFontHeightF = 0.0225
lires@tmYLLabelFontHeightF = 0.0225
lires@tmYRLabelFontHeightF = 0.0225
lires@gsnLeftStringFontHeightF   = 0.022
lires@gsnRightStringFontHeightF  = 0.022
lires@tmYLLabelsOn   = False
lires@tmYLOn         = False

lires@tmYROn         = True
lires@tmYRLabelsOn   = True
lires@tmYRMode = "Explicit"
lires@tmYRValues = fspan(0, 80, 5)
lires@tmYRLabels = (/"0","20N","40N","60N","80N"/)
lires@tmYRMinorOn = True
lires@tmYRMinorValues = fspan(0, 80, 9)

lires@tmYLValues = fspan(0, 80, 5)
lires@tmYLLabels = (/"0","20N","40N","60N","80N"/)
lires@tmYLMinorOn = True
lires@tmYLMinorValues = fspan(0, 80, 9);lires@tmYLMode = "Explicit"


lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)
lires@lgPerimOn = False
lires@lgLabelFontHeightF = 0.018
lires@vpWidthF  = 0.3
lires@vpHeightF = 0.26
lires@vpXF      = 0.06
lires@vpYF      = 0.75

lires@gsnLeftString  = ""
lires@gsnRightString = ""
lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)
; lires@xyMarkers = (/1,1,1/)
; lires@xyMarkerSizes = (/0.015,0.001,0.001/)
; lires@xyMarkerColors= (/"black","black","black"/)
; lires@xyMarkerThicknesses = (/4.0,2.0,1.0/)
lires@xyLineThicknesses  = (/2.5,8.0,8.0/)
lires@xyLineColors       = (/"red","blue","orange"/)
lires@xyDashPatterns     = (/0,4,2/)

lires@tmYLMajorThicknessF = 1.5
lires@tmYRMajorThicknessF = 1.5
lires@tmXBMajorThicknessF = 1.5
lires@tmYLMinorThicknessF = 1.0
lires@tmYRMinorThicknessF = 1.0
lires@tmXBMinorThicknessF = 1.0
lires@tmYLMajorLengthF = 0.015
lires@tmYRMajorLengthF = 0.015
lires@tmXBMajorLengthF = 0.015
lires@tmYLMinorLengthF = 0.01
lires@tmYRMinorLengthF = 0.01
lires@tmXBMinorLengthF = 0.01
lires@tmXBLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.015
lires@tmYLMajorOutwardLengthF = 0.015
lires@tmYRMajorOutwardLengthF = 0.015
lires@tmXBMinorOutwardLengthF = 0.01
lires@tmYLMinorOutwardLengthF = 0.01
lires@tmYRMinorOutwardLengthF = 0.01






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
; ores@cnFillScales         = (/0.6,0.6/)    ; change densities                                                                         
ores@cnMonoFillColor      =True                                                                                                     
ores@cnFillDotSizeF       = 0.003    
ores@cnFillDrawOrder ="postDraw"
ores@cnFillColor = "gray30"

ores@cnLevels        = (/0.05/)   ;;90%
; ores@cnLevels        = (/-0.312,0.312/)   ;;95%
ores@cnMonoFillPattern    = False            ; want multiple patterns                                                               
ores@cnFillPatterns       = (/5,-1/)     ; the patterns                                                                         
ores@cnMonoFillScale      = False            ; want different densities                                                             
ores@cnFillScales         = (/0.8,0.8/)    ; change densities   




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
topooo = new(14, "graphic") 






R2_all = where(abs(R2_detrend).ge.0.0, 1.0, 0.0)
R2_significant = where(abs(R2_detrend).ge.0.1478, 1.0, 0.0)
do ilat = 0,75
R2_all(ilat,:) = R2_all(ilat,:) * cos((ilat+10) /180.0 * 3.14159)
R2_significant(ilat,:) = R2_significant(ilat,:) * cos((ilat+10) /180.0 * 3.14159)
end do
ratio = dim_sum(dim_sum(R2_significant)) / dim_sum(dim_sum(R2_all)) * 100.0
print(ratio)








;;-----------------------------------------------------
; pltType = "eps"
pltType = "pdf"
pltName = "Figure2"
wks = gsn_open_wks(pltType, pltName)

;;----------------------- (a-b) -----------------------
res@mpMinLatF = 20
res@mpMaxLatF = 88
res@mpMinLonF = -15
res@mpMaxLonF = 180
res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/)
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
res@tmXBMinorValues = fspan(-15, 180, 14)
res@tmXBLabelFontHeightF = 0.032
res@tmYLLabelFontHeightF = 0.032
res@tmYLMode = "Explicit"
res@tmYLValues = (/20,40,60,80/)
res@tmYLLabels = (/"20N","40N","60N","80N"/)
res@tmYLMinorValues = fspan(10, 90, 9)

res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.025
res@tmXBMajorLengthF = 0.025
res@tmYLMinorLengthF = 0.018
res@tmXBMinorLengthF = 0.018
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.025
res@tmYLMajorOutwardLengthF = 0.025
res@tmXBMinorOutwardLengthF = 0.018
res@tmYLMinorOutwardLengthF = 0.018

res@cnFillPalette = "BlueWhiteOrangeRed"
res@cnLevels      = (/-1.5,-0.5,0.5,1.5/)
res@cnFillColors  = (/75,100,-1,155,175/)
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 0.9
plot(0) = gsn_csm_contour_map_ce(wks, t_diff, res)
ores@cnFillColor = "red4"
contour2(0) = gsn_csm_contour(wks, t_t_Posit, ores)    
overlay(plot(0), contour2(0))
contour3(0) = gsn_csm_contour(wks, t_t_Posit, ores)    
overlay(plot(0), contour3(0))
ores@cnFillColor = "blue"
contour4(0) = gsn_csm_contour(wks, t_t_Negat, ores)    
overlay(plot(0), contour4(0))
contour5(0) = gsn_csm_contour(wks, t_t_Negat, ores)    
overlay(plot(0), contour5(0))

sres@cnLevels         = (/2700.0/)
plot2 = gsn_csm_contour(wks, distance, sres)
overlay(plot(0), plot2)
sres@cnLevels         = (/5200.0/)
plot3 = gsn_csm_contour(wks, distance, sres)
overlay(plot(0), plot3)


vcres@vcRefAnnoString1        = "2"
vcres@vcRefMagnitudeF         = 2
vcres@vcLineArrowColor        = "darkgreen"
vcres@vcRefAnnoFontColor      = "darkgreen"
vcres@vcRefAnnoArrowFillColor = "darkgreen"
vcres@vcRefAnnoArrowEdgeColor = "darkgreen"
vcres@vcRefAnnoArrowUseVecColor = True
vcres@vcLineArrowColor        = "darkgreen"
vcres@vcLineArrowThicknessF   = 1.4
vector(0)=gsn_csm_vector(wks, u_diff, v_diff, vcres)
overlay(plot(0), vector(0))


res@cnFillPalette = "MPL_PiYG"
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnLevels      = (/-30,-5,5,30.0/)
res@cnFillColors  = (/95,85,-1,45,30/)
res@tmXBLabelsOn  = True
res@tmYLLabelsOn = False
res@mpGeophysicalLineThicknessF = 1.0
res@mpGeophysicalLineColor      = "gray15"
plot(1) = gsn_csm_contour_map_ce(wks, h_diff, res)

ores@cnFillColor = "red4"
contour2(1) = gsn_csm_contour(wks, h_t_Posit, ores)    
overlay(plot(1), contour2(1))
contour3(1) = gsn_csm_contour(wks, h_t_Posit, ores)    
overlay(plot(1), contour3(1))
ores@cnFillColor = "darkgreen"
contour4(1) = gsn_csm_contour(wks, h_t_Negat, ores)    
overlay(plot(1), contour4(1))
contour5(1) = gsn_csm_contour(wks, h_t_Negat, ores)    
overlay(plot(1), contour5(1))

vcres@vcRefAnnoString1        = "2"
vcres@vcRefMagnitudeF         = 2
vcres@vcLineArrowColor        = "blue1"
vcres@vcRefAnnoFontColor      = "blue1"
vcres@vcRefAnnoArrowFillColor = "blue1"
vcres@vcRefAnnoArrowEdgeColor = "blue1"
vcres@vcRefAnnoArrowUseVecColor = True
vcres@vcLineArrowColor        = "blue1"
vcres@vcLineArrowThicknessF   = 1.3
vector(1)=gsn_csm_vector(wks, Fx, Fy, vcres)
overlay(plot(1), vector(1))





;;----------------------- (d-e) -----------------------
res@mpMinLatF = 20
res@mpMaxLatF = 88
res@mpMinLonF = -60
res@mpMaxLonF = 135
res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/)
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
res@tmXBMinorValues = fspan(-75, 120, 14)
res@tmXBLabelsOn  = True
res@tmYLLabelsOn = True
res@cnFillPalette = "BlueWhiteOrangeRed"
res@cnLevels      = (/-1.5,-0.5,0.5,1.5/)
res@cnFillColors  = (/75,100,-1,155,175/)
res@mpGeophysicalLineThicknessF = 0.9
res@mpGeophysicalLineColor      = "gray35"
plot(2) = gsn_csm_contour_map_ce(wks, t2_diff, res)
ores@cnFillColor = "red4"
contour2(2) = gsn_csm_contour(wks, t2_t_Posit, ores)    
overlay(plot(2), contour2(2))
contour3(2) = gsn_csm_contour(wks, t2_t_Posit, ores)    
overlay(plot(2), contour3(2))
ores@cnFillColor = "blue"
contour4(2) = gsn_csm_contour(wks, t2_t_Negat, ores)    
overlay(plot(2), contour4(2))
contour5(2) = gsn_csm_contour(wks, t2_t_Negat, ores)    
overlay(plot(2), contour5(2))

sres@cnLevels         = (/2700.0/)
plot2 = gsn_csm_contour(wks, distance2, sres)
overlay(plot(2), plot2)



vcres@vcRefAnnoString1        = "2"
vcres@vcRefMagnitudeF         = 2
vcres@vcLineArrowColor        = "darkgreen"
vcres@vcRefAnnoFontColor      = "darkgreen"
vcres@vcRefAnnoArrowFillColor = "darkgreen"
vcres@vcRefAnnoArrowEdgeColor = "darkgreen"
vcres@vcRefAnnoArrowUseVecColor = True
vcres@vcLineArrowColor        = "darkgreen"
vcres@vcLineArrowThicknessF   = 1.4
vector(2)=gsn_csm_vector(wks, u2_diff, v2_diff, vcres)
overlay(plot(2), vector(2))


res@cnFillPalette = "MPL_PiYG"
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnLevels      = (/-30,-5,5,30.0/)
res@cnFillColors  = (/95,85,-1,45,30/)
res@tmXBLabelsOn  = True
res@tmYLLabelsOn = False
res@mpGeophysicalLineThicknessF = 1.0
res@mpGeophysicalLineColor      = "gray15"
plot(3) = gsn_csm_contour_map_ce(wks, h2_diff, res)

ores@cnFillColor = "red4"
contour2(3) = gsn_csm_contour(wks, h2_t_Posit, ores)    
overlay(plot(3), contour2(3))
contour3(3) = gsn_csm_contour(wks, h2_t_Posit, ores)    
overlay(plot(3), contour3(3))
ores@cnFillColor = "darkgreen"
contour4(3) = gsn_csm_contour(wks, h2_t_Negat, ores)    
overlay(plot(3), contour4(3))
contour5(3) = gsn_csm_contour(wks, h2_t_Negat, ores)    
overlay(plot(3), contour5(3))

vcres@vcRefAnnoString1        = "2"
vcres@vcRefMagnitudeF         = 2
vcres@vcLineArrowColor        = "blue1"
vcres@vcRefAnnoFontColor      = "blue1"
vcres@vcRefAnnoArrowFillColor = "blue1"
vcres@vcRefAnnoArrowEdgeColor = "blue1"
vcres@vcRefAnnoArrowUseVecColor = True
vcres@vcLineArrowColor        = "blue1"
vcres@vcLineArrowThicknessF   = 1.3
vector(3)=gsn_csm_vector(wks, Fx2, Fy2, vcres)
overlay(plot(3), vector(3))



;;------------------------ figs (a,c)，text -----------------------
txres               = True
txres@txFont = "helvetica"
txres@txFontColor = "red"
txres@txFontHeightF = 0.024
eum1 = gsn_add_text(wks,plot(0), "2700 km", 80,77, txres)
eum2 = gsn_add_text(wks,plot(0), "5200 km", 140,76, txres)
eum3 = gsn_add_text(wks,plot(2), "2700 km", -33,25, txres)
; eum4 = gsn_add_text(wks,plot(2), "5200 km", 85,76, txres)







degree_EE_Low = smth9_Wrap(degree_EE_Low, 0.5, 0.25, True)
degree_WE_Low = smth9_Wrap(degree_WE_Low, 0.5, 0.25, True)
;;--------------------- figures (c,f) degree inside western and eastern Europe ----------------------
res@mpMinLatF = 20
res@mpMaxLatF = 88
res@mpMinLonF = -15
res@mpMaxLonF = 180
res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/)
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
res@tmXBMinorValues = fspan(-15, 180, 14)

res@cnFillPalette = "WhiteYellowOrangeRed"
res@cnLevels     := (/60.0,85,110/)
res@cnFillColors := (/-1,60,140,215/)
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 1.0
plot(4) = gsn_csm_contour_map_ce(wks, degree_EE_High, res)
sres@cnLineColor = "blue"
sres@cnLevels    := (/42.0,50.0/)
contour2(4) = gsn_csm_contour(wks, degree_EE_Low, sres)    
overlay(plot(4), contour2(4))





res@mpMinLatF = 20
res@mpMaxLatF = 88
res@mpMinLonF = -60
res@mpMaxLonF = 135
res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/)
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
res@tmXBMinorValues = fspan(-75, 120, 14)

res@cnFillPalette = "WhiteYellowOrangeRed"
res@cnLevels     := (/30.0,40,50/)
res@cnFillColors := (/-1,60,140,215/)
res@tmXBLabelsOn  = True

plot(5) = gsn_csm_contour_map_ce(wks, degree_WE_High, res)
sres@cnLineColor = "blue"
sres@cnLevels    := (/21.0,25.0/)
contour2(5) = gsn_csm_contour(wks, degree_WE_Low, sres)    
overlay(plot(5), contour2(5))










;;----------------------- (g-h) explained variances (R2) -----------------------
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
res@tmXBLabelsOn = True
res@tmXBLabelFontHeightF = 0.014
res@tmYLLabelsOn = True
res@tmYLLabelFontHeightF = 0.014

res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.01
res@tmXBMajorLengthF = 0.01
res@tmYLMinorLengthF = 0.007
res@tmXBMinorLengthF = 0.007
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.01
res@tmYLMajorOutwardLengthF = 0.01
res@tmXBMinorOutwardLengthF = 0.007
res@tmYLMinorOutwardLengthF = 0.007

res@cnFillPalette = "WhiteYellowOrangeRed"
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnLevels     = (/30,40.0,50.0,60,70/)
res@cnFillColors = (/-1,40,85,130,185,235/)
res@tmXBLabelsOn = True
res@mpGeophysicalLineThicknessF = 0.7
plot(6) = gsn_csm_contour_map_ce(wks, R2_detrend, res)





lires@trXMinF = 15.0
lires@trXMaxF = 55.0
lires@tmXBMode = "Explicit"
delete(lires@tmXBValues)
delete(lires@tmXBLabels)
lires@tmXBValues = fspan(20.0, 50.0, 4)
lires@tmXBLabels = (/"20","30","40","50%"/)
lires@xyLineColors       = (/"orange","blue","orange"/)
lires@xyLineThicknesses  = (/1.3,8.0,8.0/)
lires@xyDashPatterns     = (/14,4,2/)
plot(7) = gsn_csm_xy(wks, R2_detrend_lat_High(0:65), lat(0:65), lires)

lires@xyLineColors       = (/"red","blue","orange"/)
lires@xyLineThicknesses  = (/2.5,8.0,8.0/)
lires@xyDashPatterns     = (/0,4,2/)
contour1(7) = gsn_csm_xy(wks, R2_detrend_lat(0:65), lat(0:65), lires)  
overlay(plot(7), contour1(7))


;;----------------------------
x=(/45.0,45.0/)
y=(/10, 85.0/)
   lnres   =    True
   lnres@gsLineColor="blue"
   lnres@gsLineThicknessF=2.0
   lnres@gsLineDashPattern=2
   lnres@gsnFrame = False
   dum1=gsn_add_polyline(wks,plot(7),x,y,lnres)


;;----------------------------
   txres = True
   txres@txFont = "helvetica"
   txres@txFontHeightF = 0.021
   txres@txFontColor = "red"
   eur5 = gsn_add_text(wks,plot(7), "AT", 35, 20, txres)
   txres@txFontColor = "orange"
   eur6 = gsn_add_text(wks,plot(7), "local-high", 25, 48, txres)




  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar = False

  pres@gsnPanelTop  = 0.875
  pres@gsnPanelBottom = 0.645
  pres@gsnPanelLeft  = 0.01
  pres@gsnPanelRight = 0.33
  pres@amJust  ="TopLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.013
  gsn_panel(wks,plot(0),(/1,1/),pres)
  pres@gsnPanelLeft  = 0.33
  pres@gsnPanelRight = 0.625
  gsn_panel(wks,plot(1),(/1,1/),pres)
  pres@gsnPanelLeft  = 0.64
  pres@gsnPanelRight = 0.995
  pres@gsnPanelLabelBar = True
  pres@lbOrientation        = "Vertical"
  pres@pmLabelBarWidthF     = 0.03
  pres@pmLabelBarHeightF    = 0.09
  pres@pmLabelBarOrthogonalPosF = 0.0
  ; pres@pmLabelBarParallelPosF = -0.1
  pres@lbLabelFontHeightF  = 0.011
  gsn_panel(wks,plot(4),(/1,1/),pres)


  pres@lbOrientation        = "Horizontal"
  pres@pmLabelBarWidthF     = 0.25
  pres@pmLabelBarHeightF    = 0.04
  pres@pmLabelBarOrthogonalPosF = -0.005
  ; pres@pmLabelBarParallelPosF = -0.1
  pres@lbLabelFontHeightF  = 0.011
  pres@gsnPanelLabelBar = True
  pres@gsnPanelTop  = 0.73
  pres@gsnPanelBottom = 0.5
  pres@gsnPanelLeft  = 0.01
  pres@gsnPanelRight = 0.334
  ; pres@amJust  ="BottomLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.013
  gsn_panel(wks,plot(2),(/1,1/),pres)
  pres@gsnPanelTop  = 0.73-0.001
  pres@gsnPanelBottom = 0.5-0.001
  pres@gsnPanelLeft  = 0.33
  pres@gsnPanelRight = 0.63
  gsn_panel(wks,plot(3),(/1,1/),pres)
  pres@gsnPanelTop  = 0.74
  pres@gsnPanelBottom = 0.51
  pres@gsnPanelLeft  = 0.64
  pres@gsnPanelRight = 1.0
  pres@lbOrientation        = "Vertical"
  pres@pmLabelBarWidthF     = 0.03
  pres@pmLabelBarHeightF    = 0.09
  pres@pmLabelBarOrthogonalPosF = -0.002
  ; pres@pmLabelBarParallelPosF = -0.1
  pres@lbLabelFontHeightF  = 0.011
  gsn_panel(wks,plot(5),(/1,1/),pres)



  pres@lbOrientation        = "Vertical"
  pres@pmLabelBarWidthF     = 0.045
  pres@pmLabelBarHeightF    = 0.130
  pres@pmLabelBarOrthogonalPosF = -0.005
  pres@lbLabelFontHeightF  = 0.0115
  pres@gsnPanelTop  = 0.50
  pres@gsnPanelBottom = 0.24
  pres@gsnPanelLeft  = 0.005
  pres@gsnPanelRight = 0.865
  pres@amJust  ="BottomLeft"
  gsn_panel(wks,plot(6),(/1,1/),pres)
  pres@gsnPanelLeft  = 0.77
  pres@gsnPanelRight = 0.97
  pres@amJust  ="TopLeft"
  gsn_panel(wks,plot(7),(/1,1/),pres)




  txres               = True
  txres@txFontHeightF = 0.016
  txres@txFont = "helvetica"
  txres@txFontColor = "black"
  gsn_text_ndc(wks, "Composite   for   Eastern   and   Western   European   Heatwaves",  0.490, 0.87, txres)
  gsn_text_ndc(wks, "Interannual   Explained   Variances   between   HWI   and   AT    (detrended)",  0.495, 0.488, txres)


  txres@txFontHeightF = 0.012
  txres@txFont = "helvetica"
  txres@txFontColor = "black"
  gsn_text_ndc(wks, "T2m  &  UV500",  0.175, 0.835, txres)
  gsn_text_ndc(wks, "H500  &  WAF500",  0.47, 0.835, txres)
  gsn_text_ndc(wks, "Cross-degree  centrality",  0.783, 0.835, txres)

;;------------------------- units ------------------------;; 
  txres@txFontHeightF = 0.0105
  txres@txFont = "helvetica"
  txres@txFontColor = "black"
  gsn_text_ndc(wks, "K",  0.305, 0.532, txres)
  gsn_text_ndc(wks, "gpm",  0.605, 0.532, txres)
  gsn_text_ndc(wks, "x1000",  0.95, 0.585, txres)
  gsn_text_ndc(wks, "x1000",  0.95, 0.72, txres)
  gsn_text_ndc(wks, "%",  0.742, 0.32, txres)



txres@txFontColor = "black"
txres@txFontHeightF = 0.0155
txres@txAngleF = 0.0
txres@txFont = "helvetica-bold"
gsn_text_ndc(wks, "a",  0.055, 0.835, txres)
gsn_text_ndc(wks, "b",  0.345, 0.835, txres)
gsn_text_ndc(wks, "c",  0.655, 0.835, txres)
gsn_text_ndc(wks, "d",  0.055, 0.70, txres)
gsn_text_ndc(wks, "e",  0.345, 0.70, txres)
gsn_text_ndc(wks, "f",  0.655, 0.70, txres)
gsn_text_ndc(wks, "g",  0.055, 0.465, txres)
gsn_text_ndc(wks, "h",  0.775, 0.465, txres)



frame(wks)
   delete(wks)   
;    system("convert -geometry 1800x1800 -density 1200x1200 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
;   system("rm " + pltName + "." + pltType + " -f") 





end





