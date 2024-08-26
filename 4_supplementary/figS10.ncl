load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin
;;********************************** read lat, lon ************************************;;
f00 = addfile("data/topo_adapt_era5_air2m.nc","r")
topo0 = f00->topo({10:85},{0:360})
lat = topo0&lat
lon = topo0&lon
dim_0 = dimsizes(topo0)

topo22 = f00->topo({0:85},{0:360})
lat2 = topo22&lat
lon2 = topo22&lon
dim_2 = dimsizes(topo22)



;;--------------------------------------------------------------------------------------------------------------------;;
;;----------------------------------------    CMIP6 pressure intensity    --------------------------------------------;;
;;--------------------------------------------------------------------------------------------------------------------;;
f2 = addfile("data/historical_HighLow_intensity_yearly_29models.nc","r")
h500Positive_yearlyC = doubletofloat(f2->zg_Positive_intensity(:,:,{10:85},:))
h500Positive_yearlyC@_FillValue = 1.0e20
h500Negative_yearlyC = doubletofloat(f2->zg_Negative_intensity(:,:,{10:85},:))
h500Negative_yearlyC@_FillValue = 1.0e20

value = 1.0e20
if (any(isnan_ieee(h500Positive_yearlyC))) then
      replace_ieeenan (h500Positive_yearlyC, value, 0)
      h500Positive_yearlyC@_FillValue = value
end if
h500Positive_yearlyC@_FillValue = 1.0e20
if (any(isnan_ieee(h500Negative_yearlyC))) then
      replace_ieeenan (h500Negative_yearlyC, value, 0)
      h500Negative_yearlyC@_FillValue = value
end if
h500Negative_yearlyC@_FillValue = 1.0e20

h500Positive_yearlyC!2 = "lat"
h500Positive_yearlyC!3 = "lon"
h500Positive_yearlyC&lat = lat
h500Positive_yearlyC&lon = lon
copy_VarCoords(h500Positive_yearlyC, h500Negative_yearlyC)


;;-----------------------  H500 trend  ----------------------;;
trend_HighC = regCoef_n(fspan(1,36,36), h500Positive_yearlyC, 0,1)
trend_HighC = dim_rmvmean_n_Wrap(trend_HighC, 2)
copy_VarCoords(h500Positive_yearlyC(:,0,:,:), trend_HighC)
printVarSummary(trend_HighC)
printMinMax(trend_HighC, 1)

trend_LowC = regCoef_n(fspan(1,36,36), h500Negative_yearlyC, 0,1)
trend_LowC = dim_rmvmean_n_Wrap(trend_LowC, 2)
copy_VarCoords(h500Negative_yearlyC(:,0,:,:), trend_LowC)







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
res@tmXBLabelFontHeightF = 0.024
res@tmYLLabelFontHeightF = 0.024
res@tmXBTickSpacingF = 45.0
res@tmYLTickSpacingF = 20.0

res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/) + 165+3
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
res@tmXBMinorOn = True
res@tmXBMinorValues = fspan(0, 360, 25)+3





res@mpMinLatF = 10
res@mpMaxLatF = 83
res@mpMinLonF = 0
res@mpMaxLonF = 360
res@mpCenterLonF = 180-165-3
res@mpFillOn                    = False
res@mpOutlineOn                 = True
res@mpGeophysicalLineThicknessF = 0.9
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
res@cnLevelSelectionMode = "ExplicitLevels"
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





lires = True
lires@tmBorderThicknessF = 1.5
lires@gsnDraw        = False
lires@gsnFrame       = False
lires@lgPerimOn = False
lires@lgLabelFontHeightF = 0.018
lires@vpWidthF  = 0.473
lires@vpHeightF = 0.15
lires@vpXF      = 0.06
lires@vpYF      = 0.5


lires@tiXAxisOn      = False
lires@tiYAxisOn      = False
lires@tmXTLabelsOn   = False
lires@tmXTOn         = False


lires@gsnXYBarChart = True
lires@gsnXYBarChartBarWidth = 0.8
lires@gsnYRefLine = 0.0

lires@trYMinF = 0.0
lires@trYMaxF = 5.5
lires@trXMinF = -0.2 - 0.2
lires@trXMaxF = 24.2 + 0.2


lires@tmXBLabelFontHeightF = 0.012
lires@tmYLLabelFontHeightF = 0.0125
lires@tmYRLabelFontHeightF = 0.0125
lires@gsnLeftStringFontHeightF   = 0.022
lires@gsnRightStringFontHeightF  = 0.022
lires@tmYRLabelsOn   = False
lires@tmYROn         = True



lires@tmYLMode = "Explicit"
lires@tmYLValues = fspan(0.0, 5.0, 6)
lires@tmYLLabels = fspan(0.0, 5.0, 6)
lires@tmYLMinorOn = True
lires@tmYLMinorValues = fspan(0.0, 5.5, 12)

lires@tmYRMode = "Explicit"
lires@tmYRValues = fspan(0.0, 5.0, 6)
lires@tmYRLabels = fspan(0.0, 5.0, 6)
lires@tmYRMinorOn = True
lires@tmYRMinorValues = fspan(0.0, 5.5, 12)

lires@tmXBMode = "Explicit"
lires@tmXBValues = (/2.0,5.0,8.0,11.0,14.0,17.0,20.0,23.0/) - 0.5
lires@tmXBLabels = (/"NWNA","ENA","WE","WE+EE","EE","CC","EA","NEA"/)
lires@tmXBMinorOn = False

lires@tmYLMajorThicknessF = 2.5
lires@tmXBMajorThicknessF = 2.5
lires@tmYLMinorThicknessF = 1.5
lires@tmXBMinorThicknessF = 1.5
lires@tmYLMajorLengthF = 0.009
lires@tmXBMajorLengthF = 0.009
lires@tmYLMinorLengthF = 0.005
lires@tmXBMinorLengthF = 0.005
lires@tmXBLabelDeltaF = -0.4
lires@tmYLLabelDeltaF = -0.4
lires@tmXBMajorOutwardLengthF = 0.009
lires@tmYLMajorOutwardLengthF = 0.009
lires@tmXBMinorOutwardLengthF = 0.005
lires@tmYLMinorOutwardLengthF = 0.005







leres = True
leres@tmBorderThicknessF = 2.2
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


leres@xyMarkLineModes = (/"Lines","Lines","Lines"/)

leres@lgPerimOn = False
leres@lgLabelFontHeightF = 0.018
leres@vpWidthF  = 0.34
leres@vpHeightF = 0.275
leres@vpXF      = 0.06
leres@vpYF      = 0.75

leres@gsnLeftString  = ""
leres@gsnRightString = ""
leres@xyMarkLineModes = (/"Lines","Lines","Lines"/)
leres@xyLineThicknesses  = (/5.0,8.0,8.0/)
leres@xyLineColors       = (/"red","blue","orange"/)
leres@xyDashPatterns     = (/0,4,2/)

leres@tmYLMajorThicknessF = 2.5
leres@tmYRMajorThicknessF = 2.5
leres@tmXBMajorThicknessF = 2.5
leres@tmYLMinorThicknessF = 1.5
leres@tmYRMinorThicknessF = 1.5
leres@tmXBMinorThicknessF = 1.5
leres@tmYLMajorLengthF = 0.016
leres@tmYRMajorLengthF = 0.016
leres@tmXBMajorLengthF = 0.016
leres@tmYLMinorLengthF = 0.0115
leres@tmYRMinorLengthF = 0.0115
leres@tmXBMinorLengthF = 0.0115
leres@tmXBLabelDeltaF = -0.25
leres@tmYLLabelDeltaF = -0.25
leres@tmXBMajorOutwardLengthF = 0.016
leres@tmYLMajorOutwardLengthF = 0.016
leres@tmYRMajorOutwardLengthF = 0.016
leres@tmXBMinorOutwardLengthF = 0.0115
leres@tmYLMinorOutwardLengthF = 0.0115
leres@tmYRMinorOutwardLengthF = 0.0115





cnres@tmBorderThicknessF = 2.2
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

cnres@tmYLMajorThicknessF = 2.5
cnres@tmXBMajorThicknessF = 2.5
cnres@tmYLMinorThicknessF = 1.5
cnres@tmXBMinorThicknessF = 1.5
cnres@tmYLMajorLengthF = 0.009
cnres@tmXBMajorLengthF = 0.009
cnres@tmYLMinorLengthF = 0.006
cnres@tmXBMinorLengthF = 0.006
cnres@tmXBLabelDeltaF = -0.4
cnres@tmYLLabelDeltaF = -0.4
cnres@tmXBMajorOutwardLengthF = 0.01
cnres@tmYLMajorOutwardLengthF = 0.01
cnres@tmXBMinorOutwardLengthF = 0.006
cnres@tmYLMinorOutwardLengthF = 0.006



;************************* dotting ***************************;
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




 




;**********************plot**************************
plot = new(30, "graphic")
vector = new(30, "graphic")
contour1 = new(30, "graphic")
contour2 = new(30, "graphic")
contour3 = new(30, "graphic")
contour4 = new(30, "graphic")
contour5 = new(18, "graphic")
contour6 = new(18, "graphic")
contour7 = new(18, "graphic")
contour8 = new(18, "graphic")
contour9 = new(18, "graphic")
contour10 = new(18, "graphic")
contour11 = new(18, "graphic")
topooo = new(30, "graphic") 




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
pltName = "FigureS10_High_Rtrends_29models"
wks = gsn_open_wks(pltType, pltName)
res@cnFillPalette = colorsA
res@cnLevels     = (/-40,-30,-20,-15,-10,-5,-3.0,-1.0,0.0,1.0,3.0,5,10,15,20,30,40/)
res@cnFillColors = (/19,18,17,16,15,14,13,12,11,1,2,3,4,5,6,7,8,9/)
res@mpMinLatF = 10
res@mpMaxLatF = 85


res@tmXBLabelsOn = True
res@gsnCenterStringFontHeightF = 0.028
res@gsnLeftStringFontHeightF = 0.028
CMIP_name = (/"ACCESS-CM2", "AWI-ESM-1-1-LR","BCC-ESM1","CNRM-CM6-1", "CNRM-CM6-1-HR", \
              "CNRM-ESM2-1", "CanESM5", "EC-Earth3-Veg", "EC-Earth3", "FGOALS-f3-L", \
              "FGOALS-g3", "GFDL-CM4", "GISS-E2-1-G", "HadGEM3-GC31-LL","HadGEM3-GC31-MM", \
              "IITM-ESM", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", \
              "MIROC6", "MPI-ESM1-2-HAM", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", \
              "NorESM2-LM", "NorESM2-MM", "TaiESM1", "UKESM1-0-LL"/)
; left_label = (/"(a)","(b)","(c)","(d)","(e)","(f)",\
;               "(g)","(h)","(i)","(j)","(k)","(l)",\
;               "(m)","(n)","(o)","(p)","(q)","(r)"/)
do i =0,28
;; Bottom two plots.
    if(any(i.eq.(/26,27,28/))) then
      res@tmXBOn       = True
      res@tmXBLabelsOn = True
    else
      res@tmXBOn       = True
      res@tmXBLabelsOn = False
    end if
;; Leftmost four plots.
    if(any(i.eq.(/0,3,6,9,12,15,18,21,24,27/))) then
      res@tmYLOn       = True
      res@tmYLLabelsOn = True
    else
      res@tmYLOn       = True
      res@tmYLLabelsOn = False
    end if

res@gsnCenterString = CMIP_name(i)
; res@gsnLeftString = " "+left_label(i)
plot(i) = gsn_csm_contour_map_ce(wks, trend_HighC(i,:,:), res)
end do


;;------------------------ boxes -----------------------
y1=45.0
y2=65
x1=25
x2=55
x=(/x1,x2,x2,x1,x1/)
y=(/y2,y2,y1,y1,y2/)
   lnres   =    True
   lnres@gsLineThicknessF=2.6
   lnres@gsnFrame = False
dumm = new(30, "graphic")

   lnres@gsLineColor="green"
   lnres@gsLineDashPattern=0
do i=0,28

   dumm(i) =gsn_add_polyline(wks,plot(i),x,y,lnres)
end do





  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar = False
  pres@lbOrientation        = "Vertical"
  pres@pmLabelBarWidthF     = 0.045
  pres@pmLabelBarHeightF    = 0.58
  pres@pmLabelBarOrthogonalPosF = -0.0
  ; pres@pmLabelBarParallelPosF = -0.1
  pres@lbLabelFontHeightF  = 0.012


  pres@gsnPanelTop  = 0.98
  pres@gsnPanelBottom = 0.0
  pres@gsnPanelLeft  = 0.034
  pres@gsnPanelRight = 0.35
  pres@amJust  ="BottomLeft"

  pres@pmLabelBarOrthogonalPosF = 0.005
  gsn_panel(wks,plot(0:28:3),(/10,1/),pres)

  pres@gsnPanelLeft  = 0.35
  pres@gsnPanelRight = 0.64
  pres@pmLabelBarOrthogonalPosF = 0.005
  gsn_panel(wks,plot(1:28:3),(/10,1/),pres)

  pres@gsnPanelLabelBar = True
  pres@gsnPanelTop  = 0.981
  pres@gsnPanelBottom = 0.0
  pres@gsnPanelLeft  = 0.64
  pres@gsnPanelRight = 0.993
  pres@pmLabelBarOrthogonalPosF = 0.005
  gsn_panel(wks,plot(2:28:3),(/10,1/),pres)




  txres               = True
  txres@txFontHeightF = 0.016
  txres@txFont = "helvetica"
  txres@txFontColor = "red"
  txres@txAngleF = 0.0
  gsn_text_ndc(wks, "Zonally   Asymmetric   Trends   in   High-pressure   Intensity",  0.498, 0.90, txres)





frame(wks)
   delete(wks)  
   system("convert -geometry 1800x1800 -density 1200x1200 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
  ; system("rm " + pltName + "." + pltType + " -f") 





end






