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

;;--------------  remove ocean, dtrended  -------------;;
topo2 = conform(hwi_yearly, topo0, (/1,2/))
hwi_yearly = where(topo2.ge.0.0, hwi_yearly, 1.0e20)
hwi_yearly@_FillValue = 1.0e20
copy_VarCoords(hwi_yearly2, hwi_yearly)

hwi_detrend = dtrend_msg_n(fspan(1, 44, 44), hwi_yearly, False, False, 0)
copy_VarCoords(hwi_yearly, hwi_detrend)







;;-----------------------------------------------------------;;
;;-----------------------   read AT   -----------------------;;
f1 = addfile("data/Teleconnection_Tmax90th_H500_1_5std_removeZonalTrend.nc","r")
degree21 = f1->Teleconnection(:,0:,:) 
f2 = addfile("data/Teleconnection_Tmax90th_H500_1_5std_removeZonalTrend_ConcurrentDays.nc","r")
degree22 = f2->Teleconnection(:,0:,:)
Teleconnection_1_5std = degree21 + degree22
delete([/degree21,degree22/])

f1 = addfile("data/Teleconnection_Tmax90th_H500_90th_removeZonalTrend.nc","r")
degree21 = f1->Teleconnection(:,0:,:) 
f2 = addfile("data/Teleconnection_Tmax90th_H500_90th_removeZonalTrend_ConcurrentDays.nc","r")
degree22 = f2->Teleconnection(:,0:,:)
Teleconnection_90th = degree21 + degree22
delete([/degree21,degree22/])


Teleconnection_1_5std@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(Teleconnection_1_5std))) then
      replace_ieeenan (Teleconnection_1_5std, value, 0)
      Teleconnection_1_5std@_FillValue = value
end if
Teleconnection_90th@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(Teleconnection_90th))) then
      replace_ieeenan (Teleconnection_90th, value, 0)
      Teleconnection_90th@_FillValue = value
end if




;;-------------  remove ocean  -------------;;
topo3 = conform(Teleconnection_1_5std, topo0, (/1,2/))
Teleconnection_1_5std = where(topo3.ge.0.0, Teleconnection_1_5std, 1.0e20)
Teleconnection_1_5std@_FillValue = 1.0e20
Teleconnection_1_5std!1 = "lat"
Teleconnection_1_5std!2 = "lon"
Teleconnection_1_5std&lat = lat
Teleconnection_1_5std&lon = lon
Teleconnection_90th = where(topo3.ge.0.0, Teleconnection_90th, 1.0e20)
Teleconnection_90th@_FillValue = 1.0e20
copy_VarCoords(Teleconnection_1_5std, Teleconnection_90th)


Teleconnection_1_5std_detrend = dtrend_msg_n(fspan(1, 44, 44), Teleconnection_1_5std, False, False, 0)
copy_VarCoords(Teleconnection_1_5std, Teleconnection_1_5std_detrend)
Teleconnection_90th_detrend = dtrend_msg_n(fspan(1, 44, 44), Teleconnection_90th, False, False, 0)
copy_VarCoords(Teleconnection_90th, Teleconnection_90th_detrend)





;;---------------------------------------------------------------------;;
;;----------------------- explained variance R2  ----------------------;;
cor_1_5std = escorc_n(hwi_yearly, Teleconnection_1_5std, 0, 0)
R2_1_5std = cor_1_5std * cor_1_5std * 100.0
copy_VarCoords(Teleconnection_1_5std(0,:,:), cor_1_5std)
copy_VarCoords(Teleconnection_1_5std(0,:,:), R2_1_5std)

cor_1_5std_detrend = escorc_n(hwi_detrend, Teleconnection_1_5std_detrend, 0, 0)
R2_1_5std_detrend = cor_1_5std_detrend * cor_1_5std_detrend * 100.0
copy_VarCoords(Teleconnection_1_5std(0,:,:), cor_1_5std_detrend)
copy_VarCoords(Teleconnection_1_5std(0,:,:), R2_1_5std_detrend)

R2_1_5std_detrend_lat = dim_median_n(R2_1_5std_detrend, 1)




cor_90th = escorc_n(hwi_yearly, Teleconnection_90th, 0, 0)
R2_90th = cor_90th * cor_90th * 100.0
copy_VarCoords(Teleconnection_1_5std(0,:,:), cor_90th)
copy_VarCoords(Teleconnection_1_5std(0,:,:), R2_90th)

cor_90th_detrend = escorc_n(hwi_detrend, Teleconnection_90th_detrend, 0, 0)
R2_90th_detrend = cor_90th_detrend * cor_90th_detrend * 100.0
copy_VarCoords(Teleconnection_1_5std(0,:,:), cor_90th_detrend)
copy_VarCoords(Teleconnection_1_5std(0,:,:), R2_90th_detrend)

R2_90th_detrend_lat = dim_median_n(R2_90th_detrend, 1)












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
res@mpGeophysicalLineThicknessF = 0.7
res@mpGeophysicalLineColor      = "gray30";"gray20"
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
lires@xyLineThicknesses  = (/2.3,8.0,8.0/)
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
pltName = "FigureS21_R2_1_5std"
wks = gsn_open_wks(pltType, pltName)




;;----------------------- explained variances R2 -----------------------
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
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelsOn = True
res@tmYLLabelFontHeightF = 0.015

res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.012
res@tmXBMajorLengthF = 0.012
res@tmYLMinorLengthF = 0.008
res@tmXBMinorLengthF = 0.008
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.012
res@tmYLMajorOutwardLengthF = 0.012
res@tmXBMinorOutwardLengthF = 0.008
res@tmYLMinorOutwardLengthF = 0.008

res@cnFillPalette = "WhiteYellowOrangeRed"
res@cnLevels     := (/30,40.0,50.0,60,70/)
res@cnFillColors := (/-1,40,85,130,185,235/)
res@tmXBLabelsOn = True
plot(0) = gsn_csm_contour_map_ce(wks, R2_1_5std_detrend, res)



lires@trXMinF = 12.0
lires@trXMaxF = 58.0
lires@tmXBMode = "Explicit"
lires@tmXBValues := fspan(20.0, 50.0, 4)
lires@tmXBLabels := (/"20","30","40","50%"/)
lires@xyLineColors       = (/"red","blue","orange"/)
plot(1) = gsn_csm_xy(wks, R2_1_5std_detrend_lat(0:65), lat(0:65), lires)


x=(/45.0,45.0/)
y=(/10, 85.0/)
   lnres   =    True
   lnres@gsLineColor="blue"
   lnres@gsLineThicknessF=1.8
   lnres@gsLineDashPattern=2
   lnres@gsnFrame = False
   dum1=gsn_add_polyline(wks,plot(1),x,y,lnres)








res@cnFillPalette = "WhiteYellowOrangeRed"
res@cnLevels     := (/30,40.0,50.0,60,70/)
res@cnFillColors := (/-1,40,85,130,185,235/)
res@tmXBLabelsOn = True
plot(2) = gsn_csm_contour_map_ce(wks, R2_90th_detrend, res)



lires@trXMinF = 12.0
lires@trXMaxF = 58.0
lires@tmXBMode = "Explicit"
lires@tmXBValues := fspan(20.0, 50.0, 4)
lires@tmXBLabels := (/"20","30","40","50%"/)
lires@xyLineColors       = (/"red","blue","orange"/)
plot(3) = gsn_csm_xy(wks, R2_90th_detrend_lat(0:65), lat(0:65), lires)


x=(/45.0,45.0/)
y=(/10, 85.0/)
   lnres   =    True
   lnres@gsLineColor="blue"
   lnres@gsLineThicknessF=1.8
   lnres@gsLineDashPattern=2
   lnres@gsnFrame = False
   dum3=gsn_add_polyline(wks,plot(3),x,y,lnres)






  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar = True
  pres@lbOrientation        = "Horizontal"
  pres@pmLabelBarWidthF     = 0.3
  pres@pmLabelBarHeightF    = 0.045
  pres@pmLabelBarOrthogonalPosF = -0.005
  ; pres@pmLabelBarParallelPosF = -0.1
  pres@lbLabelFontHeightF  = 0.0125




  pres@lbOrientation        = "Vertical"
  pres@pmLabelBarWidthF     = 0.045
  pres@pmLabelBarHeightF    = 0.130
  pres@pmLabelBarOrthogonalPosF = 0.0
  pres@lbLabelFontHeightF  = 0.0115
  pres@gsnPanelTop  = 0.76
  pres@gsnPanelBottom = 0.5
  pres@gsnPanelLeft  = 0.02
  pres@gsnPanelRight = 0.875
  pres@amJust  ="BottomLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.016
  pres@gsnPanelFigureStrings = (/"(a)","(c)"/)
  gsn_panel(wks,plot(0),(/1,1/),pres)
  pres@gsnPanelLeft  = 0.79
  pres@gsnPanelRight = 0.995
  pres@amJust  ="TopLeft"
  pres@gsnPanelFigureStrings = (/"(b)","(c)"/)
  gsn_panel(wks,plot(1),(/1,1/),pres)


  pres@gsnPanelTop  = 0.50
  pres@gsnPanelBottom = 0.24
  pres@gsnPanelLeft  = 0.02
  pres@gsnPanelRight = 0.875
  pres@amJust  ="BottomLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.016
  pres@gsnPanelFigureStrings = (/"(c)","(c)"/)
  gsn_panel(wks,plot(2),(/1,1/),pres)
  pres@gsnPanelLeft  = 0.79
  pres@gsnPanelRight = 0.995
  pres@amJust  ="TopLeft"
  pres@gsnPanelFigureStrings = (/"(d)","(c)"/)
  gsn_panel(wks,plot(3),(/1,1/),pres)




  txres               = True
  txres@txFontHeightF = 0.016
  txres@txFont = "helvetica"
  txres@txFontColor = "black"

  gsn_text_ndc(wks, "R~S~2~N~  (HWI,  AT  based  on  1.5  std)",  0.513, 0.74, txres)
  gsn_text_ndc(wks, "R~S~2~N~  (HWI,  AT  based  on  90th  or  10th  threshold)",  0.513, 0.48, txres)

  

frame(wks)
   delete(wks)  
   system("convert -geometry 1800x1800 -density 1200x1200 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
;   system("rm " + pltName + "." + pltType + " -f") 





end






