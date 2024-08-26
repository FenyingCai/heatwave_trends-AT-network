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
hwi_yearly = doubletofloat(f3->heatwave_intensity(0:35,10:85,:))
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

;;-------------- remove ocean, remove trends -------------;;
topo2 = conform(hwi_yearly, topo0, (/1,2/))
hwi_yearly = where(topo2.ge.0.0, hwi_yearly, 1.0e20)
hwi_yearly@_FillValue = 1.0e20
copy_VarCoords(hwi_yearly2, hwi_yearly)

time = fspan(1, 36, 36)
hwi_trend = regCoef_n(time, hwi_yearly, 0, 0) * 10.0
copy_VarCoords(hwi_yearly(0,:,:), hwi_trend)
hwi_trend_cor = escorc_n(time, hwi_yearly, 0, 0)
copy_VarCoords(hwi_yearly(0,:,:), hwi_trend_cor)
hwi_trend_t = rtest(hwi_trend_cor, 36, 0)
copy_VarCoords(hwi_yearly(0,:,:), hwi_trend_t)

hwi_trend = where(topo0.ge.0.0, hwi_trend, 1.0e20)
hwi_trend_t = where(topo0.ge.0.0, hwi_trend_t, 1.0e20)
hwi_trend_t = where(hwi_trend_cor.ge.0.0, hwi_trend_t, 1.0e20)
hwi_trend@_FillValue = 1.0e20
hwi_trend_t@_FillValue = 1.0e20

hwi_trend_lat = dim_avg_n_Wrap(hwi_trend, 1)


;;-------------- remove 31-degree moving means -------------;;
hwi_trend_31moving = new((/76,360/), float)
do ilat = 0,75
  if (ilat.le.15) then
    hwi_trend_31moving(ilat,:) = dim_avg(hwi_trend_lat(0:ilat+15))
  end if
  if (ilat.ge.60) then
    hwi_trend_31moving(ilat,:) = dim_avg(hwi_trend_lat(ilat-15:75))
  end if
  if (ilat.gt.15 .and. ilat.lt.60) then
    hwi_trend_31moving(ilat,:) = dim_avg(hwi_trend_lat(ilat-15:ilat+15))
  end if
end do


hwi_trend_ano = hwi_trend - hwi_trend_31moving
copy_VarCoords(hwi_trend, hwi_trend_ano)





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
res@tmXBLabelFontHeightF = 0.014
res@tmYLLabelFontHeightF = 0.014
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
res@mpGeophysicalLineThicknessF = 0.6
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
lires@tmBorderThicknessF = 2.2
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
lires@xyLineThicknesses  = (/5.0,8.0,8.0/)
lires@xyLineColors       = (/"red","blue","orange"/)
lires@xyDashPatterns     = (/0,4,2/)

lires@tmYLMajorThicknessF = 2.5
lires@tmYRMajorThicknessF = 2.5
lires@tmXBMajorThicknessF = 2.5
lires@tmYLMinorThicknessF = 1.5
lires@tmYRMinorThicknessF = 1.5
lires@tmXBMinorThicknessF = 1.5
lires@tmYLMajorLengthF = 0.019
lires@tmYRMajorLengthF = 0.019
lires@tmXBMajorLengthF = 0.019
lires@tmYLMinorLengthF = 0.015
lires@tmYRMinorLengthF = 0.015
lires@tmXBMinorLengthF = 0.015
lires@tmXBLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.019
lires@tmYLMajorOutwardLengthF = 0.019
lires@tmYRMajorOutwardLengthF = 0.019
lires@tmXBMinorOutwardLengthF = 0.015
lires@tmYLMinorOutwardLengthF = 0.015
lires@tmYRMinorOutwardLengthF = 0.015






;*************************dotting***************************
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
pltType = "eps"
pltName = "FigureS8__trends_2014"
wks = gsn_open_wks(pltType, pltName)


;;----------------------- 热浪趋势 -----------------------
res@mpMinLatF = 10
res@mpMaxLatF = 88
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
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelsOn = True
res@tmYLLabelFontHeightF = 0.015

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






delete(res@cnFillPalette)
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnFillPalette = colorsA
res@cnLevels     = (/-5.0,-4.0,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0/)
res@cnFillColors = (/19,18,17,16,15,14,13,12,11,1,2,3,4,5,6,7,8,9/)
res@tmXBLabelsOn = True
plot(1) = gsn_csm_contour_map_ce(wks, hwi_trend_ano, res)






  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar = True
  pres@lbOrientation        = "Vertical"
  pres@pmLabelBarWidthF     = 0.048
  pres@pmLabelBarHeightF    = 0.145
  pres@pmLabelBarOrthogonalPosF = -0.005
  pres@lbLabelFontHeightF  = 0.0115

  pres@gsnPanelTop  = 0.75
  pres@gsnPanelBottom = 0.5
  pres@gsnPanelLeft  = 0.03
  pres@gsnPanelRight = 0.97
  pres@amJust  ="BottomLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.013
  pres@gsnPanelFigureStrings = (/"(a)","(c)"/)
  gsn_panel(wks,plot(0),(/1,1/),pres)

  pres@pmLabelBarOrthogonalPosF = 0.003
  pres@gsnPanelTop  = 0.5
  pres@gsnPanelBottom = 0.25
  pres@gsnPanelLeft  = 0.03
  pres@gsnPanelRight = 0.97
  ; pres@amJust  ="TopLeft"
  pres@gsnPanelFigureStrings = (/"(b)","(c)"/)
  gsn_panel(wks,plot(1),(/1,1/),pres)







  txres               = True
  txres@txFontHeightF = 0.016
  txres@txFont = "helvetica"
  txres@txFontColor = "black"
  gsn_text_ndc(wks, "Trends   in   HWI     (ERA5,   1979-2014)",  0.43, 0.740, txres)
  gsn_text_ndc(wks, "Zonally   asymmetric   trends",  0.43, 0.490, txres)



frame(wks)
   delete(wks)  
   system("convert -geometry 1200x1200 -density 800x800 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
;   system("rm " + pltName + "." + pltType + " -f") 





end






