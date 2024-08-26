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




;;----------------------------------------------;;
;;--------------    read HWI    -------------;;
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

time = fspan(1, 44, 44)
hwi_trend = regCoef_n(time, hwi_yearly, 0, 0) * 10.0
copy_VarCoords(hwi_yearly(0,:,:), hwi_trend)
hwi_trend_cor = escorc_n(time, hwi_yearly, 0, 0)
copy_VarCoords(hwi_yearly(0,:,:), hwi_trend_cor)
hwi_trend_t = rtest(hwi_trend_cor, 44, 0)
copy_VarCoords(hwi_yearly(0,:,:), hwi_trend_t)

hwi_trend = where(topo0.ge.0.0, hwi_trend, 1.0e20)
hwi_trend_t = where(topo0.ge.0.0, hwi_trend_t, 1.0e20)
hwi_trend_t = where(hwi_trend_cor.ge.0.0, hwi_trend_t, 1.0e20)
hwi_trend@_FillValue = 1.0e20
hwi_trend_t@_FillValue = 1.0e20

hwi_trend_lat = dim_avg_n_Wrap(hwi_trend, 1)


;;-------------- minus 31-degree moving mean -------------;;
hwi_trend_31moving = new((/76/), float)
do ilat = 0,75
  if (ilat.le.15) then
    hwi_trend_31moving(ilat) = dim_avg(hwi_trend_lat(0:ilat+15))
  end if
  if (ilat.ge.60) then
    hwi_trend_31moving(ilat) = dim_avg(hwi_trend_lat(ilat-15:75))
  end if
  if (ilat.gt.15 .and. ilat.lt.60) then
    hwi_trend_31moving(ilat) = dim_avg(hwi_trend_lat(ilat-15:ilat+15))
  end if
end do



;;-------------- X times of zonal mean -------------;;
hwi_trend_times = new((/76,360/), float)
do ilat = 0,75
    hwi_trend_times(ilat,:) = hwi_trend(ilat,:) / hwi_trend_31moving(ilat)
end do
copy_VarCoords(hwi_trend, hwi_trend_times)





;;---------------------------------------------------------------;;
;;--------------    read high-pressure intensity    -------------;;
f2 = addfile("data/hgt_ano_44years_JJAyearly.nc","r")
high_yearly = doubletofloat(f2->hgt_Positive_intensity(:,{500},{10:85},:))
high_yearly@_FillValue = 1.0e20

high_yearly@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(high_yearly))) then
      replace_ieeenan (high_yearly, value, 0)
      high_yearly@_FillValue = value
end if

high_yearly!1 = "lat"
high_yearly!2 = "lon"
high_yearly&lat = lat
high_yearly&lon = lon
high_yearly2 = high_yearly



time = fspan(1, 44, 44)
high_trend = regCoef_n(time, high_yearly, 0, 0) * 10.0 / 9.8
copy_VarCoords(high_yearly(0,:,:), high_trend)
high_trend_cor = escorc_n(time, high_yearly, 0, 0)
copy_VarCoords(high_yearly(0,:,:), high_trend_cor)
high_trend_t = rtest(high_trend_cor, 44, 0)
copy_VarCoords(high_yearly(0,:,:), high_trend_t)

high_trend_t = where(high_trend_cor.ge.0.0, high_trend_t, 1.0e20)
high_trend@_FillValue = 1.0e20
high_trend_t@_FillValue = 1.0e20

high_trend_lat = dim_avg_n_Wrap(high_trend, 1)






;;--------------------------------------------------------------;;
;;--------------    read low-pressure intensity    -------------;;
f2 = addfile("data/hgt_ano_44years_JJAyearly.nc","r")
low_yearly = doubletofloat(f2->hgt_Negative_intensity(:,{500},{10:85},:))
low_yearly@_FillValue = 1.0e20

low_yearly@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(low_yearly))) then
      replace_ieeenan (low_yearly, value, 0)
      low_yearly@_FillValue = value
end if

low_yearly!1 = "lat"
low_yearly!2 = "lon"
low_yearly&lat = lat
low_yearly&lon = lon
low_yearly2 = low_yearly


time = fspan(1, 44, 44)
low_trend = regCoef_n(time, low_yearly, 0, 0) * 10.0 / 9.8
copy_VarCoords(low_yearly(0,:,:), low_trend)
low_trend_cor = escorc_n(time, low_yearly, 0, 0)
copy_VarCoords(low_yearly(0,:,:), low_trend_cor)
low_trend_t = rtest(low_trend_cor, 44, 0)
copy_VarCoords(low_yearly(0,:,:), low_trend_t)


low_trend_t = where(low_trend_cor.ge.0.0, low_trend_t, 1.0e20)
low_trend@_FillValue = 1.0e20
low_trend_t@_FillValue = 1.0e20

low_trend_lat = dim_avg_n_Wrap(low_trend, 1)








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
res@tmYLMajorLengthF = 0.012
res@tmXBMajorLengthF = 0.012
res@tmYLMinorLengthF = 0.009
res@tmXBMinorLengthF = 0.009
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.013
res@tmYLMajorOutwardLengthF = 0.013
res@tmXBMinorOutwardLengthF = 0.01
res@tmYLMinorOutwardLengthF = 0.01





res@mpFillOn                    = False
res@mpOutlineOn                 = True
res@mpGeophysicalLineThicknessF = 0.8
res@mpGeophysicalLineColor      = "gray35";"gray20"
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
vcres@vcRefAnnoFontHeightF    = 0.014   ;参考箭头字体大小
vcres@vcRefAnnoString1On      = True
vcres@vcRefAnnoString2On      = False
vcres@vcLineArrowThicknessF   = 3.0            ; make vectors larger
vcres@vcVectorDrawOrder = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.26 ;;正向下
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
sres@cnLineThicknessF = 6.0
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

lires@tmXBLabelFontHeightF = 0.0255
lires@tmYLLabelFontHeightF = 0.0255
lires@tmYRLabelFontHeightF = 0.0255
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
lires@vpHeightF = 0.29
lires@vpXF      = 0.06
lires@vpYF      = 0.75

lires@gsnLeftString  = ""
lires@gsnRightString = ""
lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)

lires@xyLineThicknesses  = (/2.8,8.0,8.0/)
lires@xyLineColors       = (/"red","blue","orange"/)
lires@xyDashPatterns     = (/0,4,2/)

lires@tmYLMajorThicknessF = 1.5
lires@tmYRMajorThicknessF = 1.5
lires@tmXBMajorThicknessF = 1.5
lires@tmYLMinorThicknessF = 1.0
lires@tmYRMinorThicknessF = 1.0
lires@tmXBMinorThicknessF = 1.0
lires@tmYLMajorLengthF = 0.018
lires@tmYRMajorLengthF = 0.018
lires@tmXBMajorLengthF = 0.018
lires@tmYLMinorLengthF = 0.013
lires@tmYRMinorLengthF = 0.013
lires@tmXBMinorLengthF = 0.013
lires@tmXBLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.019
lires@tmYLMajorOutwardLengthF = 0.019
lires@tmYRMajorOutwardLengthF = 0.019
lires@tmXBMinorOutwardLengthF = 0.014
lires@tmYLMinorOutwardLengthF = 0.014
lires@tmYRMinorOutwardLengthF = 0.014






;*************************打点***************************
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

; ores@cnLevels        = (/0.05/)   ;;90%
; ores@cnLevels        = (/-0.312,0.312/)   ;;95%
ores@cnMonoFillPattern    = False            ; want multiple patterns                                                               
; ores@cnFillPatterns       = (/5,-1/)     ; the patterns                                                                         
ores@cnMonoFillScale      = False            ; want different densities                                                             
; ores@cnFillScales         = (/0.8,0.8/)    ; change densities   




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
topooo = new(14, "graphic") 










;;-----------------------------------------------------
pltType = "eps"
pltName = "FigureS1__Original_trends"
wks = gsn_open_wks(pltType, pltName)


;;----------------------- (a) heatwave trends -----------------------
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
res@tmXBLabelFontHeightF = 0.016
res@tmYLLabelsOn = True
res@tmYLLabelFontHeightF = 0.016

res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.010
res@tmXBMajorLengthF = 0.010
res@tmYLMinorLengthF = 0.007
res@tmXBMinorLengthF = 0.007
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.011
res@tmYLMajorOutwardLengthF = 0.011
res@tmXBMinorOutwardLengthF = 0.008
res@tmYLMinorOutwardLengthF = 0.008

res@cnFillPalette = "WhiteYellowOrangeRed"
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnLevels     = (/0,2.0,4.0,6,8.0/)
res@cnFillColors = (/-1,40,85,130,185,235/)
res@tmXBLabelsOn = True
plot(0) = gsn_csm_contour_map_ce(wks, hwi_trend, res)

ores@cnLevels    = (/-0.297,0.297/)
ores@cnFillColor = "green4"
ores@cnFillPatterns       = (/11,-1,11/) 
ores@cnFillScales         = (/0.8,0.8,0.8/)    ; change densities  
contour3(0) = gsn_csm_contour(wks, hwi_trend_cor, ores)    
overlay(plot(0), contour3(0))
contour8(0) = gsn_csm_contour(wks, hwi_trend_cor, ores)    
overlay(plot(0), contour8(0))

ores@cnLevels    = (/0.5,2.0/)
ores@cnFillColor = "red4"
ores@cnFillPatterns       = (/-1,-1,17/) 
ores@cnFillScales         = (/0.5,0.5,0.5/)    ; change densities  
contour4(0) = gsn_csm_contour(wks, hwi_trend_times, ores)    
overlay(plot(0), contour4(0))
ores@cnLevels    = (/0.0,2.0/)
ores@cnFillColor = "blue"
ores@cnFillPatterns       = (/17,-1,-1/) 
ores@cnFillScales         = (/0.5,0.5,0.5/)    ; change densities  
contour5(0) = gsn_csm_contour(wks, hwi_trend_times, ores)    
overlay(plot(0), contour5(0))




lires@trXMinF = 0.0
lires@trXMaxF = 8.0
lires@tmXBMode = "Explicit"
delete(lires@tmXBValues)
delete(lires@tmXBLabels)
lires@tmXBValues = fspan(0.0, 8.0, 5)
lires@tmXBLabels = fspan(0.0, 8.0, 5)
lires@tmXBMinorValues = fspan(0.0, 8.0, 9)
lires@xyLineThicknesses  = (/1.0,8.0,8.0/)
lires@xyLineColors       = (/"orange","blue","orange"/)
lires@xyDashPatterns     = (/0,4,2/)
res@mpGeophysicalLineThicknessF = 1.0
res@mpGeophysicalLineColor      = "gray5"
plot(1) = gsn_csm_xy(wks, hwi_trend_lat(0:75), lat(0:75), lires)
lires2 = lires
lires2@xyLineColors       = (/"red","blue","orange"/)
lires2@xyLineThicknesses  = (/2.8,8.0,8.0/)
lires2@xyDashPatterns     = (/0,4,2/)
plot2 = gsn_csm_xy(wks, hwi_trend_31moving(0:75), lat(0:75), lires2)
overlay(plot(1), plot2)
plot3 = gsn_csm_xy(wks, hwi_trend_lat(0:75), lat(0:75), lires)
overlay(plot(1), plot3)




res@cnFillPalette = "BlueWhiteOrangeRed"
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnLevels     = (/-300,-100,0.0,100,300/)
res@cnFillColors = (/20,50,90,165,205,235/)
plot(2) = gsn_csm_contour_map_ce(wks, high_trend, res)
ores@cnLevels    = (/-0.297,0.297/)
ores@cnFillColor = "white"
ores@cnFillPatterns       = (/11,-1,11/) 
ores@cnFillScales         = (/0.8,0.8,0.8/)    ; change densities  
contour3(2) = gsn_csm_contour(wks, high_trend_cor, ores)    
overlay(plot(2), contour3(2))
contour8(2) = gsn_csm_contour(wks, high_trend_cor, ores)    
overlay(plot(2), contour8(2))


lires@trXMinF = 100.0
lires@trXMaxF = 400.0
lires@tmXBMode = "Explicit"
delete(lires@tmXBValues)
delete(lires@tmXBLabels)
delete(lires@tmXBMinorValues)
lires@tmXBValues = fspan(100.0, 400.0, 4)
lires@tmXBLabels = fspan(100.0, 400.0, 4)
lires@tmXBMinorValues = fspan(100.0, 400.0, 7)
lires@xyLineThicknesses  = (/2.8,8.0,8.0/)
lires@xyLineColors       = (/"red","blue","orange"/)
plot(3) = gsn_csm_xy(wks, high_trend_lat(0:75), lat(0:75), lires)





res@cnFillPalette = "BlueWhiteOrangeRed"
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnLevels     = (/-300,-100,0.0,100,300/)
res@cnFillColors = (/10,50,90,165,205,245/)
plot(4) = gsn_csm_contour_map_ce(wks, low_trend, res)
contour3(4) = gsn_csm_contour(wks, low_trend_cor, ores)    
overlay(plot(4), contour3(4))
contour8(4) = gsn_csm_contour(wks, low_trend_cor, ores)    
overlay(plot(4), contour8(4))

lires@trXMinF = -300.0
lires@trXMaxF = -0.0
lires@tmXBMode = "Explicit"
delete(lires@tmXBValues)
delete(lires@tmXBLabels)
delete(lires@tmXBMinorValues)
lires@tmXBValues = fspan(-300.0, 0.0, 4)
lires@tmXBLabels = fspan(-300.0, 0.0, 4)
lires@tmXBMinorValues = fspan(-350.0, -50.0, 7)
lires@xyLineColors       = (/"red","blue","orange"/)
plot(5) = gsn_csm_xy(wks, low_trend_lat(0:75), lat(0:75), lires)




  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar = True
  pres@lbOrientation        = "Vertical"
  pres@pmLabelBarWidthF     = 0.045
  pres@pmLabelBarHeightF    = 0.132
  pres@pmLabelBarOrthogonalPosF = -0.005
  pres@lbLabelFontHeightF  = 0.0115

  pres@gsnPanelTop  = 0.84
  pres@gsnPanelBottom = 0.6
  pres@gsnPanelLeft  = 0.03
  pres@gsnPanelRight = 0.89
  pres@amJust  ="BottomLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.013
  pres@gsnPanelFigureStrings = (/"(a)","(c)"/)
  gsn_panel(wks,plot(0),(/1,1/),pres)

  pres@gsnPanelLabelBar = False
  pres@gsnPanelLeft  = 0.795
  pres@gsnPanelRight = 0.995
  pres@amJust  ="TopLeft"
  pres@gsnPanelFigureStrings = (/"(b)","(c)"/)
  gsn_panel(wks,plot(1),(/1,1/),pres)



  pres@gsnPanelLabelBar = True
  pres@gsnPanelTop  = 0.61
  pres@gsnPanelBottom = 0.37
  pres@gsnPanelLeft  = 0.03
  pres@gsnPanelRight = 0.89
  pres@amJust  ="BottomLeft"
  pres@gsnPanelFigureStrings = (/"(c)","(c)"/)
  gsn_panel(wks,plot(2),(/1,1/),pres)

  pres@gsnPanelLabelBar = False
  pres@gsnPanelLeft  = 0.795
  pres@gsnPanelRight = 0.995
  pres@amJust  ="TopLeft"
  pres@gsnPanelFigureStrings = (/"(d)","(c)"/)
  gsn_panel(wks,plot(3),(/1,1/),pres)




  pres@gsnPanelLabelBar = True
  pres@gsnPanelTop  = 0.38
  pres@gsnPanelBottom = 0.12
  pres@gsnPanelLeft  = 0.03
  pres@gsnPanelRight = 0.89
  pres@amJust  ="BottomLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.013
  pres@gsnPanelFigureStrings = (/"(e)","(c)"/)
  gsn_panel(wks,plot(4),(/1,1/),pres)

  pres@gsnPanelLabelBar = False
  pres@gsnPanelLeft  = 0.795
  pres@gsnPanelRight = 0.995
  pres@amJust  ="BottomLeft"
  pres@gsnPanelFigureStrings = (/"(f)","(c)"/)
  gsn_panel(wks,plot(5),(/1,1/),pres)




  txres               = True
  txres@txFontHeightF = 0.019
  txres@txFont = "helvetica"
  txres@txFontColor = "black"
  gsn_text_ndc(wks, "Trends   in   HWI     (ERA5,   1979-2022)",  0.505, 0.825, txres)
  gsn_text_ndc(wks, "Trends   in   high-pressure   intensity",  0.505, 0.595, txres)
  gsn_text_ndc(wks, "Trends   in   low-pressure   intensity",  0.505, 0.355, txres)


frame(wks)
   delete(wks)  
   system("convert -geometry 1600x1600 -density 800x800 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
;   system("rm " + pltName + "." + pltType + " -f") 





end






