load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin
;;********************************** read lat, lon ************************************;;
f00 = addfile("data/topo_adapt_era5_air2m.nc","r")
topo0 = f00->topo({30:85},{0:360})
lat = topo0&lat
lon = topo0&lon
dim_0 = dimsizes(topo0)




;;-----------------------------------------------;;
;;----------------    read pdf    ---------------;;
f1 = addfile("data/PDF_heatrain_daily.nc","r")
pdf = f1->pdf




;;----------------------------------------------;;
;;--------------    read degree    -------------;;
f3 = addfile("data/Degree_HighRain_EasternEurope_WesternEurope.nc","r")
degree_EE1 = f3->degree_EE1
degree_EE2 = f3->degree_EE2
degree_WE1 = f3->degree_WE1
degree_WE2 = f3->degree_WE2

printVarSummary(degree_EE1)
printMinMax(degree_EE1, 1)
printMinMax(degree_EE2, 1)
printMinMax(degree_WE1, 1)
printMinMax(degree_WE2, 1)



;-------------- distance to EE and WE -------------;;
f0 = addfile("data/distance_i_j_land.nc","r")
distance1 = f0->distance(55,40,30:85,:)
copy_VarCoords(topo0, distance1)

f0 = addfile("data/distance_i_j_land.nc","r")
distance2 = f0->distance(46,9,30:85,:)
copy_VarCoords(topo0, distance2)











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
res@tmYLMajorLengthF = 0.00115
res@tmXBMajorLengthF = 0.00115
res@tmYLMinorLengthF = 0.00115
res@tmXBMinorLengthF = 0.009
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.0115
res@tmYLMajorOutwardLengthF = 0.0115
res@tmXBMinorOutwardLengthF = 0.009
res@tmYLMinorOutwardLengthF = 0.009





res@mpFillOn                    = False
res@mpOutlineOn                 = True
res@mpGeophysicalLineThicknessF = 0.7
res@mpGeophysicalLineColor      = "gray30" ;"gray20"
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
vcres@vcLineArrowThicknessF   = 3.0            ; make vectors larger
vcres@vcVectorDrawOrder = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.26  ;; 正向下
vcres@vcRefAnnoParallelPosF   = 1.0    ;; 正向右



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
sres@cnLevels         = (/2200.0/)
sres@cnFillColors     = (/"green4","gray20"/)
sres@cnMissingValFillColor= -1
sres@cnLineThicknessF = 2.3
sres@cnLineDashPattern = 2


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


lires@tmXBLabelFontHeightF = 0.01
lires@tmYLLabelFontHeightF = 0.01
lires@tmYRLabelFontHeightF = 0.01
lires@gsnLeftStringFontHeightF   = 0.022
lires@gsnRightStringFontHeightF  = 0.022
lires@tmYLLabelsOn   = True
lires@tmYLOn         = True

lires@tmYROn         = False
lires@tmYRLabelsOn   = False
lires@tmYRMode = "Explicit"
lires@tmYRValues = fspan(0, 80, 5)
lires@tmYRLabels = (/"0","20N","40N","60N","80N"/)
lires@tmYRMinorOn = True
lires@tmYRMinorValues = fspan(0, 80, 9)

lires@tmYLMode = "Explicit"
lires@tmYLValues = fspan(0, 0.10, 6)
lires@tmYLLabels = fspan(0, 0.10, 6)
lires@tmYLMinorOn = True
lires@tmYLMinorValues = fspan(0, 0.10, 11)

lires@xyMarkLineModes = (/"MarkLines","Lines","Lines"/)
lires@xyMarkers = (/4,1,1/)
lires@xyMarkerSizes = (/0.01,0.001,0.001/)
lires@xyMarkerColors= (/"blue","black","black"/)
lires@xyMarkerThicknesses = (/1.6,2.0,1.0/)

lires@lgPerimOn = False
lires@lgLabelFontHeightF = 0.018
lires@vpWidthF  = 0.65
lires@vpHeightF = 0.10
lires@vpXF      = 0.02
lires@vpYF      = 0.75

lires@gsnLeftString  = ""
lires@gsnRightString = ""
lires@xyLineThicknesses  = (/1.2,8.0,8.0/)
lires@xyLineColors       = (/"blue","blue","orange"/)
lires@xyDashPatterns     = (/11,4,2/)

lires@tmYLMajorThicknessF = 1.5
lires@tmYRMajorThicknessF = 1.5
lires@tmXBMajorThicknessF = 1.5
lires@tmYLMinorThicknessF = 1.0
lires@tmYRMinorThicknessF = 1.0
lires@tmXBMinorThicknessF = 1.0
lires@tmYLMajorLengthF = 0.0095
lires@tmYRMajorLengthF = 0.0095
lires@tmXBMajorLengthF = 0.0095
lires@tmYLMinorLengthF = 0.0065
lires@tmYRMinorLengthF = 0.0065
lires@tmXBMinorLengthF = 0.0065
lires@tmXBLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.0095
lires@tmYLMajorOutwardLengthF = 0.0095
lires@tmYRMajorOutwardLengthF = 0.0095
lires@tmXBMinorOutwardLengthF = 0.0065
lires@tmYLMinorOutwardLengthF = 0.0065
lires@tmYRMinorOutwardLengthF = 0.0065


lires@lgLegendOn  = False
lires@pmLegendDisplayMode  = "Always"
lires@xyExplicitLegendLabels = (/"  Heatwaves - Extreme Rainfall  (p<0.01)",""/)
lires@pmLegendHeightF = 0.018
lires@pmLegendWidthF  = 0.05
lires@lgLabelFontHeightF = 0.012
lires@pmLegendOrthogonalPosF = -1.13
lires@pmLegendParallelPosF = 0.65





;************************* dottings ***************************;
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
ores@cnMonoFillPattern    = False            ; want multiple patterns                                                               
ores@cnFillPatterns       = (/5,-1/)     ; the patterns                                                                         
ores@cnMonoFillScale      = False            ; want different densities                                                             
ores@cnFillScales         = (/0.8,0.8/)    ; change densities   




;**********************plot**************************;
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
pltName = "FigureS15_ExRain_pdf_Degree"
wks = gsn_open_wks(pltType, pltName)


;;----------------------- (a) -----------------------;;
lires@trYMinF = 0.0
lires@trYMaxF = 0.06
lires@trXMinF = 0.0
lires@trXMaxF = 10.0
lires@tmXBMode = "Explicit"
lires@tmXBValues := fspan(0.0, 10.0, 6)
lires@tmXBLabels := (/"0","2","4","6","8","10"/)
lires@tmXBMinorValues := fspan(0.0, 10.0, 11)
lires@xyLineColors       = (/"blue","blue","orange"/)
plot(0) = gsn_csm_xy(wks, fspan(0.0,10.2,52), pdf(:), lires)





;;----------------------- (b-e) -----------------------;;
res@mpMinLatF = 30
res@mpMaxLatF = 90
res@mpMinLonF = 0
res@mpMaxLonF = 360
res@mpCenterLonF = 180-165-3
res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/) + 165+3
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
res@tmXBMinorOn = True
res@tmXBMinorValues := fspan(0, 360, 25)+3
res@tmYLMinorOn = True
res@tmYLMinorValues := fspan(0, 90, 10)
res@tmXBLabelsOn = True
res@tmXBLabelFontHeightF = 0.012
res@tmYLLabelsOn = True
res@tmYLLabelFontHeightF = 0.012

res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.011
res@tmXBMajorLengthF = 0.011
res@tmYLMinorLengthF = 0.007
res@tmXBMinorLengthF = 0.007
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.011
res@tmYLMajorOutwardLengthF = 0.011
res@tmXBMinorOutwardLengthF = 0.007
res@tmYLMinorOutwardLengthF = 0.007

res@cnFillPalette = "WhiteYellowOrangeRed"
res@cnLevels     := (/50.0,100.0,200,400,600/)
res@cnFillColors := (/-1,40,85,130,185,235/)
res@tmXBLabelsOn = True

plot(1) = gsn_csm_contour_map_ce(wks, degree_EE2, res)
sres@cnLevels         = (/2400.0/)
contour1(1) = gsn_csm_contour(wks, distance1, sres)
overlay(plot(1), contour1(1))


res@cnLevels     := (/100.0,150.0,200,250,300/)
plot(2) = gsn_csm_contour_map_ce(wks, degree_WE2, res)
sres@cnLevels         = (/1800.0/)
contour1(2) = gsn_csm_contour(wks, distance2, sres)
overlay(plot(2), contour1(2))




;;--------------------- figs (a) text, line --------------------;;
txres               = True
txres@txFont = "helvetica"
txres@txFontHeightF = 0.011
txres@txFontColor = "blue"
eur1 = gsn_add_text(wks,plot(0), "2100 km", 2.1,0.051, txres)
x=(/2.4, 2.4/)
y=(/0, 0.037/)
   lnres   =    True
   lnres@gsLineThicknessF=2.2
   lnres@gsnFrame = False
   lnres@gsLineColor="blue"
   lnres@gsLineDashPattern=2
   dum10=gsn_add_polyline(wks,plot(0),x,y,lnres)
x=(/1.8, 1.8/)
y=(/0, 0.037/)
   dun10=gsn_add_polyline(wks,plot(0),x,y,lnres)




;;--------------------- figs (b-c) boxes --------------------;;
y1=45.0
y2=65
x1=25
x2=55
x:=(/x1,x2,x2,x1,x1/)  ;; eastern Europe
y:=(/y2,y2,y1,y1,y2/)
   lnres   =    True
   lnres@gsLineThicknessF=2.4
   lnres@gsnFrame = False
   lnres@gsLineColor="blue"
   lnres@gsLineDashPattern=14
   dum11=gsn_add_polyline(wks,plot(1),x,y,lnres)

y1=35.0
y2=55
x1=-10
x2=25
x=(/x1,x2,x2,x1,x1/)  ;; western Europe
y=(/y2,y2,y1,y1,y2/)
   lnres@gsLineColor="blue"
   dum22=gsn_add_polyline(wks,plot(2),x,y,lnres)



;;------------------------ figs (b-c)，text -------------------------;;
txres               = True
txres@txFont = "helvetica"
txres@txFontColor = "green4"
txres@txFontHeightF = 0.013
eum1 = gsn_add_text(wks,plot(1), "2400 km", 72,78, txres)
eum2 = gsn_add_text(wks,plot(2), "1800 km", -1,69, txres)







  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar = False
  pres@lbOrientation        = "Vertical"
  pres@pmLabelBarWidthF     = 0.04
  pres@pmLabelBarHeightF    = 0.13
  pres@pmLabelBarOrthogonalPosF = 0.0
  pres@lbLabelFontHeightF  = 0.0115

  pres@gsnPanelLabelBar = False
  pres@gsnPanelTop  = 0.8
  pres@gsnPanelBottom = 0.6
  pres@gsnPanelLeft  = 0.05
  pres@gsnPanelRight = 0.83
  pres@amJust  ="TopLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.014
  pres@gsnPanelFigureStrings = (/"(a)","(c)"/)
  gsn_panel(wks,plot(0),(/1,1/),pres)

  pres@gsnPanelLabelBar = True
  pres@gsnPanelTop  = 0.6
  pres@gsnPanelBottom = 0.4
  pres@gsnPanelLeft  = 0.05
  pres@gsnPanelRight = 0.99
  pres@amJust  ="BottomLeft"
  pres@gsnPanelFigureStrings = (/"(b)","(c)"/)
  gsn_panel(wks,plot(1),(/1,1/),pres)
  pres@gsnPanelTop  = 0.4
  pres@gsnPanelBottom = 0.2
  pres@gsnPanelFigureStrings = (/"(c)","(c)"/)
  gsn_panel(wks,plot(2),(/1,1/),pres)




  txres               = True
  txres@txFontHeightF = 0.016
  txres@txFont = "helvetica"
  txres@txFontColor = "black"
  gsn_text_ndc(wks, "Distance   distribution   of   significant   links",  0.46, 0.788, txres)
  gsn_text_ndc(wks, "Cross-degree   inside   eastern   European   extreme   rainfall",  0.46, 0.59, txres)
  gsn_text_ndc(wks, "Cross-degree   inside   western   European   extreme   rainfall",  0.46, 0.39, txres)


frame(wks)
   delete(wks)  
   system("convert -geometry 1800x1800 -density 1200x1200 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
;   system("rm " + pltName + "." + pltType + " -f") 





end






