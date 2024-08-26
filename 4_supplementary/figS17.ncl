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




;;----------------------------------------------;;
;;----------------    read HWI    --------------;;
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

;;-------------- remove ocean, remove trends -------------;;
topo2 = conform(hwi_yearly, topo0, (/1,2/))
hwi_yearly = where(topo2.ge.0.0, hwi_yearly, 1.0e20)
hwi_yearly@_FillValue = 1.0e20
copy_VarCoords(hwi_yearly2, hwi_yearly)

hwi_detrend = dtrend_msg_n(fspan(1, 44, 44), hwi_yearly, False, False, 0)
copy_VarCoords(hwi_yearly, hwi_detrend)







;;-----------------------------------------------------------;;
;;--------------------- read yearly AT  ---------------------;;
f1 = addfile("data/Teleconnection_Tmax90th_H500_removeZonalTrend.nc","r")
degree11 = f1->Teleconnection_High(:,0:,:)  ;; 10-85N
degree21 = f1->Teleconnection(:,0:,:) 
f2 = addfile("data/Teleconnection_Tmax90th_H500_removeZonalTrend_ConcurrentDays.nc","r")
degree12 = f2->Teleconnection_High(:,0:,:)  ;; 10-85N
degree22 = f2->Teleconnection(:,0:,:)
Teleconnection_High = degree11 + degree12
Teleconnection = degree21 + degree22
delete([/degree11,degree12,degree21,degree22/])

Teleconnection_High@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(Teleconnection_High))) then
      replace_ieeenan (Teleconnection_High, value, 0)
      Teleconnection_High@_FillValue = value
end if
Teleconnection@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(Teleconnection))) then
      replace_ieeenan (Teleconnection, value, 0)
      Teleconnection@_FillValue = value
end if


;;------------- remove ocean -------------;;
topo3 = conform(Teleconnection_High, topo0, (/1,2/))
Teleconnection_High = where(topo3.ge.0.0, Teleconnection_High, 1.0e20)
Teleconnection = where(topo3.ge.0.0, Teleconnection, 1.0e20)
Teleconnection_High@_FillValue = 1.0e20
Teleconnection@_FillValue = 1.0e20
Teleconnection_High!1 = "lat"
Teleconnection_High!2 = "lon"
Teleconnection_High&lat = lat
Teleconnection_High&lon = lon
copy_VarCoords(Teleconnection_High, Teleconnection)

print(" ----------- Teleconnection_High ---------- ")
printMinMax(Teleconnection_High, 0)

Teleconnection_High_detrend = dtrend_msg_n(fspan(1, 44, 44), Teleconnection_High, False, False, 0)
copy_VarCoords(Teleconnection_High, Teleconnection_High_detrend)
Teleconnection_detrend = dtrend_msg_n(fspan(1, 44, 44), Teleconnection, False, False, 0)
copy_VarCoords(Teleconnection, Teleconnection_detrend)





;;------------------------------------------------------------;;
;;----------------------- calculate R2  ----------------------;;
cor_High = escorc_n(hwi_yearly, Teleconnection_High, 0, 0)
R2_High = cor_High * cor_High * 100.0
cor = escorc_n(hwi_yearly, Teleconnection, 0, 0)
R2 = cor * cor * 100.0
copy_VarCoords(Teleconnection_High(0,:,:), cor_High)
copy_VarCoords(Teleconnection_High(0,:,:), R2_High)
copy_VarCoords(Teleconnection_High(0,:,:), cor)
copy_VarCoords(Teleconnection_High(0,:,:), R2)


cor_detrend_High = escorc_n(hwi_detrend, Teleconnection_High_detrend, 0, 0)
R2_detrend_High = cor_detrend_High * cor_detrend_High * 100.0
cor_detrend = escorc_n(hwi_detrend, Teleconnection_detrend, 0, 0)
R2_detrend = cor_detrend * cor_detrend * 100.0
copy_VarCoords(Teleconnection_High(0,:,:), cor_detrend_High)
copy_VarCoords(Teleconnection_High(0,:,:), R2_detrend_High)
copy_VarCoords(Teleconnection_High(0,:,:), cor_detrend)
copy_VarCoords(Teleconnection_High(0,:,:), R2_detrend)



R2_lat_High = dim_median_n(R2_High, 1)
cor_detrend_lat_High = dim_median_n(cor_detrend_High, 1)
R2_detrend_lat_High = dim_median_n(R2_detrend_High, 1)
R2_lat = dim_median_n(R2, 1)
cor_detrend_lat = dim_median_n(cor_detrend, 1)
R2_detrend_lat = dim_median_n(R2_detrend, 1)


R2_detrend_High = where(R2_detrend_High.gt.0.0, R2_detrend_High, 1.0e20)
R2_detrend_High@_FillValue = 1.0e20
R2_detrend_diff = R2_detrend / R2_detrend_High
copy_VarCoords(Teleconnection_High(0,:,:), R2_detrend_diff)
R2_detrend_lat_diff = R2_detrend_lat / R2_detrend_lat_High








;;------------------------------------           ------------------------------------;;
;;------------------------------------ pdf in R2 ------------------------------------;;
lat_matrix = new((/dim_0(0),dim_0(1)/), float)
lat_matrix = 0.0
do ilat = 0,dim_0(0)-1
lat_matrix(ilat,:) = cos((ilat+10) / 180.0 * 3.14159)
end do

R2_interval = fspan(0.0, 100, 21)
print(R2_interval)
pdf = new(21, float)
pdf_High = new(21, float)

do i = 0,20
print(" i = "+i)

minT = R2_interval(i)-2.5
maxT = R2_interval(i)+2.5
if (i.eq.0) then
minT = 0.0
end if
if (i.eq.20) then
maxT = 100.0
end if

R2_detrend_new = where(R2_detrend.ge.minT.and.R2_detrend.lt.maxT, 1.0, 0.0)
pdf(i) = dim_sum_n_Wrap(dim_sum_n_Wrap(R2_detrend_new*lat_matrix, 1), 0)

R2_detrend_High_new = where(R2_detrend_High.ge.minT.and.R2_detrend_High.lt.maxT, 1.0, 0.0)
pdf_High(i) = dim_sum_n_Wrap(dim_sum_n_Wrap(R2_detrend_High_new*lat_matrix, 1), 0)

end do

land_area = where(R2_detrend.ge.0.0, 1.0, 0.0)


pdf = pdf / dim_sum(pdf)
pdf_High = pdf_High / dim_sum(pdf_High)
print(pdf_High - pdf)












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
res@mpGeophysicalLineThicknessF = 0.6
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
vcres@vcLineArrowThicknessF   = 3.0            ; make vectors larger
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
lires@xyLineThicknesses  = (/2.1,8.0,8.0/)
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
lires@tmYLMinorLengthF = 0.009
lires@tmYRMinorLengthF = 0.009
lires@tmXBMinorLengthF = 0.009
lires@tmXBLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.015
lires@tmYLMajorOutwardLengthF = 0.015
lires@tmYRMajorOutwardLengthF = 0.015
lires@tmXBMinorOutwardLengthF = 0.009
lires@tmYLMinorOutwardLengthF = 0.009
lires@tmYRMinorOutwardLengthF = 0.009









lcres = True
lcres@tmBorderThicknessF = 1.0
lcres@gsnDraw        = False
lcres@gsnFrame       = False
lcres@lgPerimOn = False
lcres@lgLabelFontHeightF = 0.018
lcres@vpWidthF  = 0.5
lcres@vpHeightF = 0.3
lcres@vpXF      = 0.06
lcres@vpYF      = 0.55


lcres@tiXAxisOn      = False
lcres@tiYAxisOn      = False
lcres@tmXTLabelsOn   = False
lcres@tmXTOn         = False


lcres@gsnXYBarChart = True
lcres@gsnXYBarChartBarWidth = 0.8
lcres@gsnYRefLine = 0.0

lcres@trYMinF = 0.65
lcres@trYMaxF = 1.0
lcres@trXMinF = -1.0
lcres@trXMaxF = 6.3


lcres@tmXBLabelFontHeightF = 0.023
lcres@tmYLLabelFontHeightF = 0.023
lcres@tmYRLabelFontHeightF = 0.023
lcres@gsnLeftStringFontHeightF   = 0.022
lcres@gsnRightStringFontHeightF  = 0.022
lcres@tmYRLabelsOn   = False
lcres@tmYROn         = True



lcres@tmYLMode = "Explicit"
lcres@tmYLValues = fspan(0.6, 1.0, 5)
lcres@tmYLLabels = (/"60%","70%","80%","90%","100%"/)
lcres@tmYLMinorOn = True
lcres@tmYLMinorValues = fspan(0.5, 1.0, 11)

lcres@tmYRMode = "Explicit"
lcres@tmYRValues = fspan(0.0, 6.0, 4)
lcres@tmYRLabels = fspan(0.0, 6.0, 4)
lcres@tmYRMinorOn = True
lcres@tmYRMinorValues = fspan(0.0, 7.0, 15)

lcres@tmXBMode = "Explicit"
lcres@tmXBValues = (/1.5, 4.5/)
lcres@tmXBLabels = (/"p<0.01","p<0.001"/)
lcres@tmXBMinorOn = False

lcres@tmYLMajorThicknessF = 1.5
lcres@tmXBMajorThicknessF = 1.5
lcres@tmYLMinorThicknessF = 1.0
lcres@tmXBMinorThicknessF = 1.0
lcres@tmYLMajorLengthF = 0.016
lcres@tmXBMajorLengthF = 0.016
lcres@tmYLMinorLengthF = 0.01
lcres@tmXBMinorLengthF = 0.01
lcres@tmXBLabelDeltaF = -0.4
lcres@tmYLLabelDeltaF = -0.4
lcres@tmXBMajorOutwardLengthF = 0.016
lcres@tmYLMajorOutwardLengthF = 0.016
lcres@tmXBMinorOutwardLengthF = 0.01
lcres@tmYLMinorOutwardLengthF = 0.01











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



R2_all = where(abs(R2_detrend).ge.0.0, 1.0, 0.0)
R2_significant = where(abs(R2_detrend).ge.0.1478, 1.0, 0.0)
do ilat = 0,75
R2_all(ilat,:) = R2_all(ilat,:) * cos((ilat+10) /180.0 * 3.14159)
R2_significant(ilat,:) = R2_significant(ilat,:) * cos((ilat+10) /180.0 * 3.14159)
end do
ratio = dim_sum(dim_sum(R2_significant)) / dim_sum(dim_sum(R2_all)) * 100.0
print(ratio)





;;-----------------------------------------------------
pltType = "eps"
pltName = "data/FigureS17_explain_variance_local_High2"
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
plot(0) = gsn_csm_contour_map_ce(wks, R2_detrend_High, res)






lires@trXMinF = 12.0
lires@trXMaxF = 58.0
lires@tmXBMode = "Explicit"
lires@tmXBValues := fspan(20.0, 50.0, 4)
lires@tmXBLabels := (/"20","30","40","50%"/)
lires@tmXBMinorValues = fspan(0, 100, 21)
lires@xyLineColors       = (/"orange","blue","orange"/)
plot(1) = gsn_csm_xy(wks, R2_detrend_lat_High(0:65), lat(0:65), lires)
lires@xyLineColors       = (/"red","blue","orange"/)
contour1(1) = gsn_csm_xy(wks, R2_detrend_lat(0:65), lat(0:65), lires)
overlay(plot(1), contour1(1))

x=(/45.0,45.0/)
y=(/10, 85.0/)
   lnres   =    True
   lnres@gsLineColor="blue"
   lnres@gsLineThicknessF=2.0
   lnres@gsLineDashPattern=2
   lnres@gsnFrame = False
   dum1=gsn_add_polyline(wks,plot(1),x,y,lnres)








lires@vpWidthF  = 0.7
lires@vpHeightF = 0.18

lires@lgLegendOn  = False
lires@pmLegendDisplayMode  = "Always"
lires@xyExplicitLegendLabels = (/"   dipolar  AT",""/)
lires@pmLegendHeightF = 0.022
lires@pmLegendWidthF  = 0.06
lires@lgLabelFontHeightF = 0.014
lires@pmLegendOrthogonalPosF = -1.05  ;; 负为向上
lires@pmLegendParallelPosF = 0.8
lires@tmXBLabelFontHeightF = 0.015
lires@tmYLLabelFontHeightF = 0.015
lires@tmYROn         = False
lires@tmYRLabelsOn   = False
lires@tmYLOn         = True
lires@tmYLLabelsOn   = True

lires@tmYLMajorLengthF = 0.013
lires@tmYRMajorLengthF = 0.013
lires@tmXBMajorLengthF = 0.013
lires@tmYLMinorLengthF = 0.008
lires@tmYRMinorLengthF = 0.008
lires@tmXBMinorLengthF = 0.008
lires@tmXBLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.013
lires@tmYLMajorOutwardLengthF = 0.013
lires@tmYRMajorOutwardLengthF = 0.013
lires@tmXBMinorOutwardLengthF = 0.008
lires@tmYLMinorOutwardLengthF = 0.008
lires@tmYRMinorOutwardLengthF = 0.008


lires@xyMarkLineModes = (/"MarkLines","Lines","Lines"/)
lires@xyMarkers = (/4,1,1/)
lires@xyMarkerSizes = (/0.013,0.013,0.013/)
lires@xyMarkerColors= (/"blue","black","black"/)
lires@xyMarkerThicknesses = (/1.8,1.0,1.0/)
lires@xyMarkerThicknessF = 1.3
lires@xyLineThicknesses = (/1.1,2.0,2.0/)

lires@trYMinF := 0.0
lires@trYMaxF := 0.15
lires@trXMinF := 0.0
lires@trXMaxF := 20.0
lires@tmXBMode = "Explicit"
lires@tmXBValues := fspan(0.0, 20.0, 6)
lires@tmXBLabels := (/"0","20%","40%","60%","80%","100%"/)
lires@tmXBMinorValues = fspan(0, 20, 21)
lires@xyLineColors       = (/"red","blue","orange"/)
lires@xyMarkerColors= (/"red","black","black"/)
lires@xyMarkers = (/7,1,1/)
plot(2) = gsn_csm_xy(wks, fspan(0.0,20,21), pdf(:), lires)
lires@xyLineColors       = (/"orange","blue","orange"/)
lires@xyMarkerColors= (/"orange","black","black"/)
lires@xyMarkers = (/4,1,1/)
lires@xyExplicitLegendLabels = (/"   local-high",""/)
lires@pmLegendOrthogonalPosF = -0.85  ;; 负为向上
contour1(2)  = gsn_csm_xy (wks, fspan(0, 20, 21), pdf_High, lires)
overlay(plot(2), contour1(2))
   

;;----------------- fig (c) vertical line -------------------
x:=(/35.53888/5, 35.53888/5/)
y:=(/0, 0.12/)
   lnres@gsLineColor="red"
   lnres@gsLineThicknessF=2.2
   lnres@gsLineDashPattern=2
   dum2=gsn_add_polyline(wks,plot(2),x,y,lnres)
x:=(/33.1583/5, 33.1583/5/)
y:=(/0, 0.093/)
   lnres@gsLineColor="orange"
   dum3=gsn_add_polyline(wks,plot(2),x,y,lnres)






lcres@gsnXYBarChartColors = (/"royalblue3","royalblue3","violetred2","violetred2","violetred2","royalblue3","violetred2","violetred2"/)
lcres@gsnXYBarChartColors = (/"royalblue3","royalblue3","violetred2","violetred2","violetred2","royalblue3","violetred2","violetred2"/)
lcres@gsnXYBarChartPatterns = (/0,0,6,6,6,6,6/)
lcres@gsnXYBarChartFillScaleF = 1.0
liresL = lcres
liresR = lcres
liresL@gsnXYBarChartColors := (/"red","red"/)
liresR@gsnXYBarChartColors := (/"orange","orange"/)
plot(3) = gsn_csm_x2y2(wks, (/1.0,4.0/), (/2.0,5.0/), (/0.8961,0.77499/), (/0.81624,0.691132/), liresL, liresR)
;;----------------- fig (d) bars -------------------
lbres = True
lbres@vpWidthF = 0.06
lbres@vpHeightF= 0.03
lbres@lbBoxMajorExtentF = 0.4
lbres@lbLabelFontHeightF = 0.0105
lbres@lbPerimColor = False
lbres@lbFillColors = "red"
lbres@lbFillPatterns = 0
gsn_labelbar_ndc(wks, 1, "  dipolar AT", 0.865, 0.455, lbres)

lbres@lbFillColors = "orange"
lbres@lbFillPatterns = 0
gsn_labelbar_ndc(wks, 1, "  local-high", 0.865, 0.435, lbres)


;;------------------------ fig (b,d) texts -----------------------
txres               = True
txres@txFont = "helvetica"
txres@txFontHeightF = 0.021

txres@txFontColor = "red"
eur1 = gsn_add_text(wks,plot(3), "90%", 1.0,0.93, txres)
eur2 = gsn_add_text(wks,plot(3), "77%", 4.0,0.80, txres)

txres@txFontColor = "orange"
eur3 = gsn_add_text(wks,plot(3), "82%", 2.0,0.85, txres)
eur4 = gsn_add_text(wks,plot(3), "69%", 5.0,0.72, txres)



txres@txFontColor = "red"
eur5 = gsn_add_text(wks,plot(1), "dipolar AT", 42, 19, txres)
txres@txFontColor = "orange"
eur6 = gsn_add_text(wks,plot(1), "local-high", 25, 47, txres)





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


  pres@gsnPanelTop  = 0.51
  pres@gsnPanelBottom = 0.25
  pres@gsnPanelLeft  = 0.02
  pres@gsnPanelRight = 0.65
  pres@amJust  ="TopLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.016
  pres@gsnPanelFigureStrings = (/"(c)","(c)"/)
  gsn_panel(wks,plot(2),(/1,1/),pres)
  pres@gsnPanelTop  = 0.51
  pres@gsnPanelBottom = 0.252
  pres@gsnPanelLeft  = 0.685
  pres@gsnPanelRight = 0.99
  pres@amJust  ="TopLeft"
  pres@gsnPanelFigureStrings = (/"(d)","(c)"/)
  gsn_panel(wks,plot(3),(/1,1/),pres)




  txres               = True
  txres@txFontHeightF = 0.016
  txres@txFont = "helvetica"
  txres@txFontColor = "black"
  gsn_text_ndc(wks, "Interannual   Explained   Variances   between   HWI   and   local-high     (detrended)",  0.519, 0.74, txres)
  gsn_text_ndc(wks, "PDF   in   R~S~2~N~",  0.35, 0.495, txres)
  gsn_text_ndc(wks, "Significant  Percentage",  0.855, 0.495, txres)


frame(wks)
   delete(wks)  
   system("convert -geometry 1800x1800 -density 1200x1200 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
;   system("rm " + pltName + "." + pltType + " -f") 





end






