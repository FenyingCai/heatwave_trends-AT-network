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




;;----------------------------------------------------------------------------;;
;;----------------------  network_max_distance (1-4kkm)  ---------------------;;
f1 = addfile("data/spatial_heatlow_strongest_distance.nc","r")
network_max_distance = f1->network_max_distance
printVarSummary(network_max_distance)
printMinMax(network_max_distance, 1)




network_max_distance_1_4kkm = where(network_max_distance.ge.1500.0.and.network_max_distance.le.4000, 1, 0.0)
lat_1_4kkm = dim_sum_n_Wrap(network_max_distance_1_4kkm, 1)
network_max_distance_all = where(network_max_distance.ge.0.0, 1, 0.0)
lat_all = dim_sum_n_Wrap(network_max_distance_all, 1)
lat_all = where(lat_all.gt.0.0, lat_all, 1.0e20)
lat_all@_FillValue = 1.0e20
percentage_1_4kkm = (lat_1_4kkm / lat_all) * 100.0

print(percentage_1_4kkm)







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
res@mpGeophysicalLineThicknessF = 0.7
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

lires@tmXBLabelFontHeightF = 0.02
lires@tmYLLabelFontHeightF = 0.02
lires@tmYRLabelFontHeightF = 0.02
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
lires@vpHeightF = 0.22
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








;;-----------------------------------------------------
pltType = "eps"
pltName = "FigureS7_network_max_distance"
wks = gsn_open_wks(pltType, pltName)

;;----------------------- (a) -----------------------
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

res@tmYLMajorThicknessF = 1.6
res@tmXBMajorThicknessF = 1.6
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.015
res@tmXBMajorLengthF = 0.015
res@tmYLMinorLengthF = 0.01
res@tmXBMinorLengthF = 0.01
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.015
res@tmYLMajorOutwardLengthF = 0.015
res@tmXBMinorOutwardLengthF = 0.01
res@tmYLMinorOutwardLengthF = 0.01

res@cnFillPalette = "WhiteBlueGreenYellowRed"
res@cnLevels      = (/500,1000,1500,2000,2500,3000,3500,4000,4500,5000.0/)
res@cnFillColors  = (/255,252,225,180,160,140,120,100,60,30,-1/)
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 0.75
plot(0) = gsn_csm_contour_map_ce(wks, network_max_distance, res)



lires@trXMinF = 0.0
lires@trXMaxF = 100.0
lires@tmXBMode = "Explicit"
lires@tmXBValues := fspan(20.0, 100.0, 5)
lires@tmXBLabels := (/"20","40","60","80","100"/)
lires@xyLineColors       = (/"springgreen4","blue","orange"/)
plot(1) = gsn_csm_xy(wks, percentage_1_4kkm(0:65), lat(0:65), lires)

x=(/85.0,85.0/)
y=(/10, 85.0/)
   lnres   =    True
   lnres@gsLineColor="blue"
   lnres@gsLineThicknessF=1.8
   lnres@gsLineDashPattern=2
   lnres@gsnFrame = False
   dum1=gsn_add_polyline(wks,plot(1),x,y,lnres)










  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar = True
  pres@lbOrientation        = "Horizontal"
  pres@pmLabelBarWidthF     = 0.60
  pres@pmLabelBarHeightF    = 0.045
  pres@pmLabelBarOrthogonalPosF = -0.002
  ; pres@pmLabelBarParallelPosF = -0.1
  pres@lbLabelFontHeightF  = 0.0115

  pres@gsnPanelTop  = 0.65
  pres@gsnPanelBottom = 0.3
  pres@gsnPanelLeft  = 0.01
  pres@gsnPanelRight = 0.74
  pres@amJust  ="BottomLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.013
  pres@gsnPanelFigureStrings = (/"(a)","(c)"/)
  gsn_panel(wks,plot(0),(/1,1/),pres)

  pres@gsnPanelTop  = 0.64
  pres@gsnPanelBottom = 0.35
  pres@gsnPanelLeft  = 0.75
  pres@gsnPanelRight = 0.99
  pres@amJust  ="TopLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.013
  pres@gsnPanelFigureStrings = (/"(b)","(c)"/)
  gsn_panel(wks,plot(1),(/1,1/),pres)



  txres               = True
  txres@txFontHeightF = 0.016
  txres@txFont = "helvetica"
  txres@txFontColor = "black"
  gsn_text_ndc(wks, "Distance   of   the   strongest   links   (km)",  0.39, 0.61, txres)
  gsn_text_ndc(wks, "Percentage  (1.5-4kkm)",  0.865, 0.61, txres)


frame(wks)
   delete(wks)  
   system("convert -geometry 2400x2400 -density 1600x1600 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
;   system("rm " + pltName + "." + pltType + " -f") 





end






