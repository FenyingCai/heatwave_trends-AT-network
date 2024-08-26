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




;;********************************** read lat, lon ************************************
f00 = addfile("data/topo_adapt_era5_air2m.nc","r")
topoC = f00->topo({0:90},{0:360})

;;--------------------  Composite analysis for EE  --------------------;;
f1 = addfile("data/era5_EEtop5_composite_T2m_UVH500.nc","r")
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

;;------------- remove ocean -------------;;
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



;;--------------------  significant test  --------------------;;
f1 = addfile("data/ttest_EEtop5_T2m_UVH500.nc","r")
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








;;---------------------------------------------------------------------;;
;;--------------------  Composite analysis for WE  --------------------;;
f1 = addfile("data/era5_WEtop5_composite_T2m_UVH500.nc","r")
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

;;------------- remove ocean -------------;;
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


;;--------------------  significant test  --------------------;;
f1 = addfile("data/ttest_WEtop5_T2m_UVH500.nc","r")
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







;-------------- geograhpical center for EE, WE -------------;;
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




;-------------- Distance to geographical centers -------------;;
f0 = addfile("data/distance_i_j_land.nc","r")
distance = f0->distance(55,40,0:90,:)
copy_VarCoords(topoC, distance)

f0 = addfile("data/distance_i_j.nc","r")
distance2 = f0->distance(46,9,0:90,:)
distance2 = where(topoC.le.1.0, distance2, 1.0e20)
distance2@_FillValue = 1.0e20
copy_VarCoords(topoC, distance2)
















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
res@tmXBLabelFontHeightF = 0.016
res@tmYLLabelFontHeightF = 0.016
res@tmXBTickSpacingF = 45.0
res@tmYLTickSpacingF = 20.0



res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.001
res@tmXBMajorLengthF = 0.001
res@tmYLMinorLengthF = 0.001
res@tmXBMinorLengthF = 0.007
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.01
res@tmYLMajorOutwardLengthF = 0.01
res@tmXBMinorOutwardLengthF = 0.007
res@tmYLMinorOutwardLengthF = 0.007





res@mpFillOn                    = False
res@mpOutlineOn                 = True
res@mpGeophysicalLineThicknessF = 0.9
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
res@cnLevelSelectionMode = "ExplicitLevels"
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
vcres@vcRefAnnoFontHeightF    = 0.014
vcres@vcRefAnnoString1On      = True
vcres@vcRefAnnoString2On      = False
vcres@vcLineArrowThicknessF   = 1.7            ; make vectors larger
vcres@vcVectorDrawOrder = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.26 ;; 正向下
vcres@vcRefAnnoParallelPosF   = 1.0   ;; 正向右



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
sres@cnLineThicknessF = 3.0
sres@cnLineDashPattern = 14


sres2@cnLevelSelectionMode = "ExplicitLevels"
sres2@cnLevels         = (/2/)
sres2@cnLineColor = "blue"
sres2@cnLineDashPattern = 14
sres2@cnLineThicknessF = 4





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
lires@xyLineThicknesses  = (/4.5,8.0,8.0/)
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
lires@tmYLMinorLengthF = 0.006
lires@tmYRMinorLengthF = 0.006
lires@tmXBMinorLengthF = 0.006
lires@tmXBLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.015
lires@tmYLMajorOutwardLengthF = 0.015
lires@tmYRMajorOutwardLengthF = 0.015
lires@tmXBMinorOutwardLengthF = 0.006
lires@tmYLMinorOutwardLengthF = 0.006
lires@tmYRMinorOutwardLengthF = 0.006






;************************* dottings ***************************
ores                 = True            
ores@gsnDraw         = False               
ores@gsnFrame        = False                 
ores@cnLineLabelsOn  = False               
ores@cnLinesOn       = False                 
ores@cnInfoLabelOn   = False                                                                                                        
ores@cnFillOn        = True                               
ores@lbLabelBarOn    = False                                                                                                            
ores@cnLevelSelectionMode = "ExplicitLevels"
ores@cnMonoFillPattern    = False            ; want multiple patterns                                                            
ores@cnMonoFillScale      = False            ; want different densities                                                                         
ores@cnMonoFillColor      =True                                                                                                     
ores@cnFillDotSizeF       = 0.003    
ores@cnFillDrawOrder ="postDraw"
ores@cnFillColor = "gray30"
ores@cnLevels        = (/0.05/)   ;; 90%
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








;;-----------------------------------------------------
pltType = "eps"
pltName = "FigureS16_EE_WE_composite"
wks = gsn_open_wks(pltType, pltName)


res@mpMinLatF = 20
res@mpMaxLatF = 88
res@mpMinLonF = -15
res@mpMaxLonF = 180
res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/)
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
res@tmXBMinorValues = fspan(-15, 180, 14)
res@tmXBLabelFontHeightF = 0.030
res@tmYLLabelFontHeightF = 0.030
res@tmYLMode = "Explicit"
res@tmYLValues = (/20,40,60,80/)
res@tmYLLabels = (/"20N","40N","60N","80N"/)
res@tmYLMinorValues = fspan(10, 90, 9)

res@tmYLMajorThicknessF = 1.6
res@tmXBMajorThicknessF = 1.6
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
res@mpGeophysicalLineThicknessF = 0.7
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
vcres@vcLineArrowThicknessF   = 1.6
vector(0)=gsn_csm_vector(wks, u_diff, v_diff, vcres)
overlay(plot(0), vector(0))


res@cnFillPalette = "MPL_PiYG"
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnLevels      = (/-30,-5,5,30.0/)
res@cnFillColors  = (/95,85,-1,45,30/)
res@tmXBLabelsOn  = True
res@tmYLLabelsOn = False
res@mpGeophysicalLineThicknessF = 0.9
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
res@mpGeophysicalLineThicknessF = 0.7
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


vcres@vcRefAnnoOrthogonalPosF = -0.275 ;; 正向下
vcres@vcRefAnnoParallelPosF   = 1.0  ;; 正向右
vcres@vcRefAnnoString1        = "2"
vcres@vcRefMagnitudeF         = 2
vcres@vcLineArrowColor        = "darkgreen"
vcres@vcRefAnnoFontColor      = "darkgreen"
vcres@vcRefAnnoArrowFillColor = "darkgreen"
vcres@vcRefAnnoArrowEdgeColor = "darkgreen"
vcres@vcRefAnnoArrowUseVecColor = True
vcres@vcLineArrowColor        = "darkgreen"
vcres@vcLineArrowThicknessF   = 1.6
vector(2)=gsn_csm_vector(wks, u2_diff, v2_diff, vcres)
overlay(plot(2), vector(2))


res@cnFillPalette = "MPL_PiYG"
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnLevels      = (/-30,-5,5,30.0/)
res@cnFillColors  = (/95,85,-1,45,30/)
res@tmXBLabelsOn  = True
res@tmYLLabelsOn = False
res@mpGeophysicalLineThicknessF = 0.9
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





;;--------------------- (a,c) text -----------------------
txres               = True
txres@txFont = "helvetica"
txres@txFontColor = "red"
txres@txFontHeightF = 0.024
eum1 = gsn_add_text(wks,plot(0), "2700 km", 80,77, txres)
eum2 = gsn_add_text(wks,plot(0), "5200 km", 140,76, txres)
eum3 = gsn_add_text(wks,plot(2), "2700 km", -33,25, txres)













  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar = True
  pres@lbOrientation        = "Vertical"
  pres@pmLabelBarWidthF     = 0.035
  pres@pmLabelBarHeightF    = 0.12
  pres@pmLabelBarOrthogonalPosF = -0.004
  ; pres@pmLabelBarParallelPosF = -0.1
  pres@lbLabelFontHeightF  = 0.0115

  pres@gsnPanelTop  = 0.75
  pres@gsnPanelBottom = 0.5
  pres@gsnPanelLeft  = 0.01
  pres@gsnPanelRight = 0.54
  pres@amJust  ="TopLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.013
  pres@gsnPanelFigureStrings = (/"(a)","(c)"/)
  gsn_panel(wks,plot(0),(/1,1/),pres)
  pres@gsnPanelLeft  = 0.5
  pres@gsnPanelRight = 0.99
  pres@gsnPanelFigureStrings = (/"(b)","(c)"/)
  gsn_panel(wks,plot(1),(/1,1/),pres)



  pres@gsnPanelTop  = 0.5
  pres@gsnPanelBottom = 0.25
  pres@gsnPanelLeft  = 0.01
  pres@gsnPanelRight = 0.54
  pres@gsnPanelFigureStringsFontHeightF = 0.013
  pres@gsnPanelFigureStrings = (/"(c)","(c)"/)
  gsn_panel(wks,plot(2),(/1,1/),pres)
  pres@gsnPanelLeft  = 0.5
  pres@gsnPanelRight = 0.99
  pres@gsnPanelFigureStrings = (/"(d)","(c)"/)
  gsn_panel(wks,plot(3),(/1,1/),pres)
 





  txres               = True
  txres@txFontHeightF = 0.016
  txres@txFont = "helvetica"
  txres@txFontColor = "black"
  gsn_text_ndc(wks, "Composite   for   Eastern   European   Heatwaves   (top 5%)",  0.49, 0.74, txres)
  gsn_text_ndc(wks, "Composite   for   Western   European   Heatwaves   (top 5%)",  0.49, 0.49, txres)


frame(wks)
   delete(wks)  
   system("convert -geometry 1600x1600 -density 800x800 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
;   system("rm " + pltName + "." + pltType + " -f") 





end






