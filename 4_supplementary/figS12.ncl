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





;;--------------------------------------------------------------------------------------------------------------------;;
;;------------------------------------------------    CMIP6 HWI    ---------------------------------------------------;;
;;--------------------------------------------------------------------------------------------------------------------;;
f4 = addfile("data/historical_HWI_yearly_29models.nc","r")
hwiC_yearly0 = doubletofloat(f4->heatwave_intensity(:,:,0:85,:))
hwiC_yearly0@_FillValue = 1.0e20
printMinMax(hwiC_yearly0,1)

value = 1.0e20
if (any(isnan_ieee(hwiC_yearly0))) then
      replace_ieeenan (hwiC_yearly0, value, 0)
      hwiC_yearly0@_FillValue = value
end if

hwiC_yearly0!2 = "lat"
hwiC_yearly0!3 = "lon"
hwiC_yearly0&lat = lat2
hwiC_yearly0&lon = lon2
;;-------------- remove ocean -------------;;
topo2 = conform(hwiC_yearly0, topo22, (/2,3/))
hwiC_yearly = where(topo2.ge.0.0, hwiC_yearly0, 1.0e20)
hwiC_yearly@_FillValue = 1.0e20
copy_VarCoords(hwiC_yearly0, hwiC_yearly)


time = fspan(1, 36, 36)
hwiC_trend = regCoef_n(time, hwiC_yearly, 0, 1) * 10.0
copy_VarCoords(hwiC_yearly(:,0,:,:), hwiC_trend)
hwiC_trend_cor = escorc_n(time, hwiC_yearly, 0, 1)
copy_VarCoords(hwiC_yearly(:,0,:,:), hwiC_trend_cor)
hwiC_trend_t = rtest(hwiC_trend_cor, 36, 0)
copy_VarCoords(hwiC_yearly(:,0,:,:), hwiC_trend_t)


;;-------------- trends -------------;;
topo3 = conform(hwiC_trend, topo22, (/1,2/))
hwiC_trend = where(topo3.ge.0.0, hwiC_trend, 1.0e20)
hwiC_trend_t = where(topo3.ge.0.0, hwiC_trend_t, 1.0e20)
hwiC_trend_t = where(hwiC_trend_cor.ge.0.0, hwiC_trend_t, 1.0e20)
hwiC_trend@_FillValue = 1.0e20
hwiC_trend_t@_FillValue = 1.0e20


;;-------------- minus 31-degree moving mean -------------;;
hwiC_trend_avg = dim_avg_n_Wrap(hwiC_trend, 1)
hwiC_trend_31moving = new((/29,dim_2(0),dim_2(1)/), float)

do kmodel = 0,28
do ilat = 0,84
  if (ilat.le.15) then
    hwiC_trend_31moving(kmodel,ilat,:) = dim_avg(hwiC_trend_avg(kmodel,0:ilat+15))
  end if
  if (ilat.ge.70) then
    hwiC_trend_31moving(kmodel,ilat,:) = dim_avg(hwiC_trend_avg(kmodel,ilat-15:85))
  end if
  if (ilat.gt.15 .and. ilat.lt.70) then
    hwiC_trend_31moving(kmodel,ilat,:) = dim_avg(hwiC_trend_avg(kmodel,ilat-15:ilat+15))
  end if
end do
end do

hwiC_trend_ano = hwiC_trend(:,10:85,:) - hwiC_trend_31moving(:,10:85,:)
copy_VarCoords(hwiC_yearly(:,0,10:85,:), hwiC_trend_ano)




;;--------------------------------------------------------------------------------------------------------------------;;
;;------------------------------------------------    CMIP6 pressure    ---------------------------------------------------;;
;;--------------------------------------------------------------------------------------------------------------------;;
f2 = addfile("/public/home/fcai/CMIP6/historical_zg_day/zg_for1AT/historical_HighLow_intensity_yearly_29models.nc","r")
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


;;-----------------------  H500 trends  ----------------------;;
trend_HighC = regCoef_n(fspan(1,36,36), h500Positive_yearlyC, 0,1)
trend_HighC = dim_rmvmean_n_Wrap(trend_HighC, 2)
copy_VarCoords(h500Positive_yearlyC(:,0,:,:), trend_HighC)

trend_LowC = regCoef_n(fspan(1,36,36), h500Negative_yearlyC, 0,1)
trend_LowC = dim_rmvmean_n_Wrap(trend_LowC, 2)
copy_VarCoords(h500Negative_yearlyC(:,0,:,:), trend_LowC)

trend_HighLow = trend_HighC - trend_LowC
copy_VarCoords(trend_LowC, trend_HighLow)



;;--------------------------------------------------------------------------------------------------------------------;;
;;----------------------------------------------    Eastern Europe    ------------------------------------------------;;
;;--------------------------------------------------------------------------------------------------------------------;;
hwiC_trend_ano_weighted = new((/29,dim_0(0),dim_0(1)/), float)
trend_HighC_weighted = new((/29,dim_0(0),dim_0(1)/), float)
trend_LowC_weighted = new((/29,dim_0(0),dim_0(1)/), float)
grids65 = new((/dim_0(0),dim_0(1)/), float)
do ilat = 0,dim_0(0)-1
  hwiC_trend_ano_weighted(:,ilat,:) = hwiC_trend_ano(:,ilat,:) * cos((ilat+10) / 180.0 * 3.14159)
  trend_HighC_weighted(:,ilat,:) = trend_HighC(:,ilat,:) * cos((ilat+10) / 180.0 * 3.14159)
  trend_LowC_weighted(:,ilat,:) = trend_LowC(:,ilat,:) * cos((ilat+10) / 180.0 * 3.14159)
  grids65(ilat,:) = cos((ilat+10) / 180.0 * 3.14159)
end do


hwiC_EE = dim_sum(dim_sum(hwiC_trend_ano_weighted(:,45-10:65-10, 25:55))) / dim_sum(dim_sum(grids65(45-10:65-10, 25:55)))
HighC_EE = dim_sum(dim_sum(trend_HighC_weighted(:,45-10:65-10, 25:55))) / dim_sum(dim_sum(grids65(45-10:65-10, 25:55)))
WE_lon = new(26, integer)
WE_lon(0:14) = ispan(345,359,1)
WE_lon(15:25) = ispan(0,10,1)
LowC_EE1 = dim_sum(dim_sum(trend_LowC_weighted(:,45-10:60-10, 330:359)))
LowC_EE2 = dim_sum(dim_sum(trend_LowC_weighted(:,25-10:40-10, 55:75)))
LowC_EE = (LowC_EE1+LowC_EE2) / (dim_sum(dim_sum(grids65(45-10:60-10, 330:359)))+dim_sum(dim_sum(grids65(25-10:40-10, 55:75))))
AT_EE = (HighC_EE+LowC_EE) / 2.0
printMinMax(hwiC_EE, 1)
printMinMax(HighC_EE, 1)
printMinMax(AT_EE, 1)
print("heat-high "+escorc(hwiC_EE, HighC_EE)+"   "+rtest(escorc(hwiC_EE, HighC_EE), 29, 0))
print("heat-low "+escorc(hwiC_EE, LowC_EE)+"   "+rtest(escorc(hwiC_EE, LowC_EE), 29, 0))
print("high-low "+escorc(HighC_EE, LowC_EE)+"   "+rtest(escorc(HighC_EE, LowC_EE), 29, 0))
print("heat-AT "+escorc(hwiC_EE, AT_EE)+"   "+rtest(escorc(hwiC_EE, AT_EE), 29, 0))

hwiC_EE_avg = dim_avg_n_Wrap(hwiC_EE, 0)
HighC_EE_avg = dim_avg_n_Wrap(HighC_EE, 0)
AT_EE_avg = dim_avg_n_Wrap(AT_EE, 0)
reg_hwiC_HighC = regCoef_n(hwiC_EE, HighC_EE, 0, 0)
reg_hwiC_AT = regCoef_n(hwiC_EE, AT_EE, 0, 0)

print("   Min-Max  ")
printMinMax(hwiC_EE, 1)
printMinMax(HighC_EE, 1)
printMinMax(AT_EE, 1)








;;--------------- EE ranking --------------;;
print(hwiC_EE)
hwiC_mintomax = hwiC_EE
qsort(hwiC_mintomax)
print(hwiC_mintomax)
do i=0,28
print("i = "+(ind(hwiC_EE.eq.hwiC_mintomax(i))+1))
end do



;;--------------------------------------------------------------------------------------------------------------------;;
;;------------------------------------------------    four models    -------------------------------------------------;;
;;--------------------------------------------------------------------------------------------------------------------;;
Positive0 = (/28,21,11,15,19,    12,26,22,1,9/)-1
Negative0 = (/6,2,24,5,29,       25,3,17,10,13/)-1
tsize = 8
Positive = Positive0(0:tsize-1)
Negative = Negative0(0:tsize-1)


thres = 5.5
hwiC_Positive = dim_avg_n_Wrap(hwiC_trend_ano(Positive,:,:), 0)
hwiC_Negative = dim_avg_n_Wrap(hwiC_trend_ano(Negative,:,:), 0)
copy_VarCoords(hwiC_trend_ano(0,:,:), hwiC_Positive)
copy_VarCoords(hwiC_trend_ano(0,:,:), hwiC_Negative)

HighLow_trend = trend_HighC - trend_LowC
HighLow_Positive = dim_avg_n_Wrap(HighLow_trend(Positive,:,:), 0)
HighLow_Negative = dim_avg_n_Wrap(HighLow_trend(Negative,:,:), 0)
copy_VarCoords(hwiC_trend_ano(0,:,:), HighLow_Positive)
copy_VarCoords(hwiC_trend_ano(0,:,:), HighLow_Negative)







;;-----------------------------------------    inter-model agreements    -------------------------------------------;;
hwiC_Positive_P = where(hwiC_trend_ano(Positive,:,:).gt.0, 1, 0)
hwiC_Positive_PP = dim_sum_n_Wrap(hwiC_Positive_P, 0)
hwiC_Positive_N = where(hwiC_trend_ano(Positive,:,:).lt.0, 1, 0)
hwiC_Positive_NN = dim_sum_n_Wrap(hwiC_Positive_N, 0)

HighLow_Positive_P = where(HighLow_trend(Positive,:,:).gt.0, 1, 0)
HighLow_Positive_PP = dim_sum_n_Wrap(HighLow_Positive_P, 0)
HighLow_Positive_N = where(HighLow_trend(Positive,:,:).lt.0, 1, 0)
HighLow_Positive_NN = dim_sum_n_Wrap(HighLow_Positive_N, 0)

copy_VarCoords(hwiC_trend_ano(0,:,:), hwiC_Positive_PP)
copy_VarCoords(hwiC_trend_ano(0,:,:), hwiC_Positive_NN)
copy_VarCoords(hwiC_trend_ano(0,:,:), HighLow_Positive_PP)
copy_VarCoords(hwiC_trend_ano(0,:,:), HighLow_Positive_NN)




hwiC_Negative_P = where(hwiC_trend_ano(Negative,:,:).gt.0, 1, 0)
hwiC_Negative_PP = dim_sum_n_Wrap(hwiC_Negative_P, 0)
hwiC_Negative_N = where(hwiC_trend_ano(Negative,:,:).lt.0, 1, 0)
hwiC_Negative_NN = dim_sum_n_Wrap(hwiC_Negative_N, 0)

HighLow_Negative_P = where(HighLow_trend(Negative,:,:).gt.0, 1, 0)
HighLow_Negative_PP = dim_sum_n_Wrap(HighLow_Negative_P, 0)
HighLow_Negative_N = where(HighLow_trend(Negative,:,:).lt.0, 1, 0)
HighLow_Negative_NN = dim_sum_n_Wrap(HighLow_Negative_N, 0)

copy_VarCoords(hwiC_trend_ano(0,:,:), hwiC_Negative_PP)
copy_VarCoords(hwiC_trend_ano(0,:,:), hwiC_Negative_NN)
copy_VarCoords(hwiC_trend_ano(0,:,:), HighLow_Negative_PP)
copy_VarCoords(hwiC_trend_ano(0,:,:), HighLow_Negative_NN)




printMinMax(HighLow_Positive_PP, 1)










;************************************************** plots **********************************************
res = True
res@tmBorderThicknessF = 1.2
res@gsnDraw      =  False
res@gsnFrame     =  False
res@gsnAddCyclic =  True
res@gsnRightString       = ""
res@gsnLeftString        = ""
res@gsnLeftStringFontHeightF   = 0.0255
res@gsnRightStringFontHeightF  = 0.0255
res@gsnCenterStringFontHeightF  = 0.018
res@tmXTLabelsOn  = False
res@tmYRLabelsOn  = False
vcres = res
res@tmXTOn        = False
res@tmYROn        = False
res@tmXBOn        = True
res@tmYLOn        = True
res@tmXBLabelFontHeightF = 0.018
res@tmYLLabelFontHeightF = 0.018
res@tmXBTickSpacingF = 45.0
res@tmYLTickSpacingF = 20.0

res@tmXBMode = "Explicit"
res@tmXBMinorOn = True
res@tmXBMinorValues = fspan(0, 360, 25)+3


res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.014
res@tmXBMajorLengthF = 0.014
res@tmYLMinorLengthF = 0.010
res@tmXBMinorLengthF = 0.010
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.014
res@tmYLMajorOutwardLengthF = 0.014
res@tmXBMinorOutwardLengthF = 0.010
res@tmYLMinorOutwardLengthF = 0.010




res@mpMinLatF = 15
res@mpMaxLatF = 85
res@mpMinLonF = 0
res@mpMaxLonF = 360
res@mpCenterLonF = 180-165-3
res@mpFillOn                    = False
res@mpOutlineOn                 = True
res@mpGeophysicalLineThicknessF = 0.8
res@mpGeophysicalLineColor      = "gray20";"gray20"
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
res@pmLabelBarWidthF     = 0.06
res@pmLabelBarHeightF    = 0.205
res@pmLabelBarOrthogonalPosF = 0.012
res@lbLabelFontHeightF  = 0.018
;res@cnMissingValFillColor = "white"
;res@gsnYRefLine = 0.0
;res@gsnYRefLineColor = "black"
;res@gsnYRefLineThicknessF = 5.0
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnFillPalette        = "precip2_17lev" 



;vcres@mpFillDrawOrder         = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.255
;vcres@vcRefAnnoSide           = "TopRight"
vcres@vcGlyphStyle            = "LineArrow"
vcres@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
vcres@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref
vcres@vcMinDistanceF          = 0.025             ; thin out vectors
vcres@vcLineArrowColor        = "black"           ; change vector color
vcres@vcRefAnnoOn             = True
; vcres@vcRefAnnoString1        = "4"
; vcres@vcRefMagnitudeF         = 4
vcres@vcRefLengthF            = 0.025            ; ref vec length
vcres@vcRefAnnoFontHeightF    = 0.015   ;参考箭头字体大小
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
lires@tmBorderThicknessF = 2.2
lires@gsnDraw        = False
lires@gsnFrame       = False
lires@tiXAxisOn      = False
lires@tiYAxisOn      = False
lires@tmXTLabelsOn   = False
lires@tmXTOn         = False

lires@tmXBLabelFontHeightF = 0.04
lires@tmYLLabelFontHeightF = 0.04
lires@tmYRLabelFontHeightF = 0.04
lires@gsnLeftStringFontHeightF   = 0.022
lires@gsnRightStringFontHeightF  = 0.022
lires@tmYLLabelsOn   = True
lires@tmYLOn         = True



lires@lgPerimOn = False
lires@lgLabelFontHeightF = 0.018
lires@vpWidthF  = 0.5
lires@vpHeightF = 0.5
lires@vpXF      = 0.06
lires@vpYF      = 0.75

lires@gsnLeftString  = ""
lires@gsnRightString = ""
lires@xyMarkLineModes = (/"Markers","Lines","Lines"/)
lires@xyMarkers = (/1,1,1/)
lires@xyMarkerSizes = (/0.012,0.001,0.001/)
lires@xyMarkerColors= (/"red","black","black"/)
lires@xyMarkerThicknesses = (/7.0,2.0,1.0/)
lires@xyLineThicknesses  = (/6.0,8.0,8.0/)
; lires@xyLineColors       = (/"red","blue","orange"/)
lires@xyDashPatterns     = (/0,4,2/)

lires@tmYLMajorThicknessF = 2.5
lires@tmYRMajorThicknessF = 2.5
lires@tmXBMajorThicknessF = 2.5
lires@tmYLMinorThicknessF = 1.5
lires@tmYRMinorThicknessF = 1.5
lires@tmXBMinorThicknessF = 1.5
lires@tmYLMajorLengthF = 0.025
lires@tmYRMajorLengthF = 0.025
lires@tmXBMajorLengthF = 0.025
lires@tmYLMinorLengthF = 0.018
lires@tmYRMinorLengthF = 0.018
lires@tmXBMinorLengthF = 0.018
lires@tmXBLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.025
lires@tmYLMajorOutwardLengthF = 0.025
lires@tmYRMajorOutwardLengthF = 0.025
lires@tmXBMinorOutwardLengthF = 0.018
lires@tmYLMinorOutwardLengthF = 0.018
lires@tmYRMinorOutwardLengthF = 0.018




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
plot = new(14, "graphic")
plot2 = new(30, "graphic")
plot3 = new(30, "graphic")
vector = new(30, "graphic")
contour1 = new(30, "graphic")
contour2 = new(30, "graphic")
contour3 = new(30, "graphic")
contour4 = new(30, "graphic")
contour5 = new(30, "graphic")
contour6 = new(30, "graphic")
contour7 = new(30, "graphic")
contour8 = new(30, "graphic")
contour9 = new(30, "graphic")
topooo = new(30, "graphic") 



;---------------------- color map ----------------------;;
colors1 = read_colormap_file("MPL_OrRd")
colors2 = read_colormap_file("WhiteBlue")
colorsA = colors1(0:20,:)
colorsA(1,:) = colors1(15,:)
colorsA(2,:) = colors1(50,:)
colorsA(3,:) = colors1(70,:)
colorsA(4,:) = colors1(56,:)

colorsA(11,:) = colors2(50,:)
colorsA(12,:) = colors2(100,:)
colorsA(13,:) = colors2(145,:)
colorsA(14,:) = colors2(130,:)







;;-----------------------------------------------------
pltType = "eps"
pltName = "FigureS12_EE_Positive_Negative"
wks = gsn_open_wks(pltType, pltName)




;;------------------------ (a-b) -----------------------;;
res@mpMinLatF = 15
res@mpMaxLatF = 80
res@mpMinLonF = -50
res@mpMaxLonF = 192
res@tmXBLabelFontHeightF = 0.021
res@tmYLLabelFontHeightF = 0.021
delete(res@tmXBValues)
delete(res@tmXBLabels)
delete(res@tmXBMinorValues)
res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/)
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
res@tmXBMinorOn = True
res@tmXBMinorValues = fspan(-360, 360, 49)
delete(res@cnFillPalette)
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnFillPalette = colorsA
res@cnLevels      = (/-5.0,-3.0,-1.0,1.0,3.0,5.0/)
res@cnFillColors  = (/13,12,11,-1,1,2,3/)
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 0.5
res@gsnLeftString = " (a)   Trend  in  heatwave  intensity  (HWI)"
plot(0) = gsn_csm_contour_map_ce(wks, hwiC_Positive, res)

ores@cnFillColor = "violetred3"
ores@cnLevels    = (/thres/)
ores@cnFillPatterns       = (/-1,11/) 
contour3(0) = gsn_csm_contour(wks, hwiC_Positive_PP, ores)    
overlay(plot(0), contour3(0))
contour8(0) = gsn_csm_contour(wks, hwiC_Positive_PP, ores)    
overlay(plot(0), contour8(0))
ores@cnFillColor = "royalblue4"
contour4(0) = gsn_csm_contour(wks, hwiC_Positive_NN, ores)
overlay(plot(0), contour4(0))
contour9(0) = gsn_csm_contour(wks, hwiC_Positive_NN, ores)    
overlay(plot(0), contour9(0))




delete(res@cnFillPalette)
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnFillPalette = "MPL_PiYG"
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 0.9
res@mpMinLatF = 15
res@mpMaxLatF = 82
res@gsnLeftString = " (b)   Trend  in  intensity  (high  -  low)"
res@cnLevels      = (/-15,-10,-4,4,10,15.0/)
res@cnFillColors  = (/104,92,78,-1,50,35,20/)
plot(1) = gsn_csm_contour_map_ce(wks, HighLow_Positive, res)

ores@cnFillColor = "violetred4"
ores@cnLevels    = (/thres/)
ores@cnFillPatterns       = (/-1,11/) 
contour3(1) = gsn_csm_contour(wks, HighLow_Positive_PP, ores)    
overlay(plot(1), contour3(1))
contour8(1) = gsn_csm_contour(wks, HighLow_Positive_PP, ores)    
overlay(plot(1), contour8(1))
ores@cnFillColor = "green4"
contour4(1) = gsn_csm_contour(wks, HighLow_Positive_NN, ores)    
overlay(plot(1), contour4(1))
contour9(1) = gsn_csm_contour(wks, HighLow_Positive_NN, ores)    
overlay(plot(1), contour9(1))





;;----------------------------------------------------------------
delete(res@cnFillPalette)
res@cnFillPalette = colorsA
res@cnLevels      = (/-5.0,-3.0,-1.0,1.0,3.0,5.0/)
res@cnFillColors  = (/13,12,11,-1,1,2,3/)
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 0.5
res@gsnLeftString = " (c)   Trend  in  heatwave  intensity  (HWI)"
plot(2) = gsn_csm_contour_map_ce(wks, hwiC_Negative, res)

ores@cnFillColor = "violetred3"
contour3(2) = gsn_csm_contour(wks, hwiC_Negative_PP, ores)    
overlay(plot(2), contour3(2))
contour8(2) = gsn_csm_contour(wks, hwiC_Negative_PP, ores)    
overlay(plot(2), contour8(2))
ores@cnFillColor = "royalblue4"
contour4(2) = gsn_csm_contour(wks, hwiC_Negative_NN, ores)
overlay(plot(2), contour4(2))
contour9(2) = gsn_csm_contour(wks, hwiC_Negative_NN, ores)    
overlay(plot(2), contour9(2))




delete(res@cnFillPalette)
delete(res@cnLevels)
delete(res@cnFillColors)
res@cnFillPalette = "MPL_PiYG"
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 0.9
res@mpMinLatF = 15
res@mpMaxLatF = 82
res@gsnLeftString = " (d)   Trend  in  intensity  (high  -  low)"
res@cnLevels      = (/-15,-10,-4,4,10,15.0/)
res@cnFillColors  = (/104,92,78,-1,50,35,20/)
plot(3) = gsn_csm_contour_map_ce(wks, HighLow_Negative, res)

ores@cnFillColor = "violetred4"
contour3(3) = gsn_csm_contour(wks, HighLow_Negative_PP, ores)    
overlay(plot(3), contour3(3))
contour8(3) = gsn_csm_contour(wks, HighLow_Negative_PP, ores)    
overlay(plot(3), contour8(3))
ores@cnFillColor = "green4"
contour4(3) = gsn_csm_contour(wks, HighLow_Negative_NN, ores)    
overlay(plot(3), contour4(3))
contour9(3) = gsn_csm_contour(wks, HighLow_Negative_NN, ores)    
overlay(plot(3), contour9(3))








;;------------------------ boxes -----------------------
y1=45.0
y2=65
x1=25
x2=55
xx1=(/x1,x2,x2,x1,x1/)  ;; EE
yy1=(/y2,y2,y1,y1,y2/)
y1=60.0
y2=45
x1=-30
x2=0
xx3=(/x1,x2,x2,x1,x1/)  ;; WE
yy3=(/y2,y2,y1,y1,y2/)
y1=30
y2=45
x1=55
x2=75
xx5=(/x1,x2,x2,x1,x1/)  ;; CA
yy5=(/y2,y2,y1,y1,y2/)
   lnres   =    True
   lnres@gsLineThicknessF=2.8
   lnres@gsnFrame = False
   lnres@gsLineColor="red"
   lnres@gsLineDashPattern=14
   dum1 = gsn_add_polyline(wks,plot(0),xx1,yy1,lnres)

   dum3 = gsn_add_polyline(wks,plot(2),xx1,yy1,lnres)
   lnres@gsLineColor="red"
   dum2 = gsn_add_polyline(wks,plot(1),xx1,yy1,lnres)
   dum4 = gsn_add_polyline(wks,plot(3),xx1,yy1,lnres)

   lnres@gsLineColor="blue"
   dum5 = gsn_add_polyline(wks,plot(1),xx3,yy3,lnres)
   dum6 = gsn_add_polyline(wks,plot(1),xx5,yy5,lnres)
   dum7 = gsn_add_polyline(wks,plot(3),xx3,yy3,lnres)
   dum8 = gsn_add_polyline(wks,plot(3),xx5,yy5,lnres)






   










;;------------------------ FIG (a) TEXT -----------------------
  txres               = True
  txres@txFont = "helvetica"
  txres@txFontHeightF = 0.018
  txres@txFontColor = "black"

  dutext2 = gsn_add_text(wks,plot(0), "hatched:", 160,33.0, txres)
  dutext3 = gsn_add_text(wks,plot(0), ">=6 models (75%)", 160,23.0, txres)













  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar = True
  pres@lbOrientation        = "Horizontal"
  pres@pmLabelBarWidthF     = 0.35
  pres@pmLabelBarHeightF    = 0.044
  pres@pmLabelBarOrthogonalPosF = -0.001
  ; pres@pmLabelBarParallelPosF = 0.07
  pres@lbLabelFontHeightF  = 0.0125


  pres@gsnPanelTop  = 0.8
  pres@gsnPanelBottom = 0.57
  pres@gsnPanelLeft  = 0.0
  pres@gsnPanelRight = 0.49
  pres@amJust  ="BottomLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.0105
  pres@gsnPanelFigureStrings = (/"","(b)"/)
  gsn_panel(wks,plot(0),(/1,1/),pres)

  pres@gsnPanelTop  = 0.555
  pres@gsnPanelBottom = 0.325
  gsn_panel(wks,plot(1),(/1,1/),pres)


  pres@gsnPanelTop  = 0.8
  pres@gsnPanelBottom = 0.57
  pres@gsnPanelLeft  = 0.50
  pres@gsnPanelRight = 0.99
  gsn_panel(wks,plot(2),(/1,1/),pres)

  pres@gsnPanelTop  = 0.555
  pres@gsnPanelBottom = 0.325
  gsn_panel(wks,plot(3),(/1,1/),pres)




  







  txres               = True
  txres@txFontHeightF = 0.016
  txres@txFont = "helvetica"
  txres@txFontColor = "red"
  gsn_text_ndc(wks, "High-EE-trend  group  (8)",  0.27, 0.83, txres)
  gsn_text_ndc(wks, "Low-EE-trend  group  (8)",  0.76, 0.83, txres)

  

frame(wks)
   delete(wks)  
   system("convert -geometry 1600x1600 -density 800x800 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
;   system("rm " + pltName + "." + pltType + " -f") 






end




