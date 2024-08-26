load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin
;;-----------------  read lat, lon  ----------------;;
f00 = addfile("data/topo_adapt_era5_air2m.nc","r")
topo0 = f00->topo({10:85},{0:360})
lat = topo0&lat
lon = topo0&lon
dim_0 = dimsizes(topo0)



;;------------  degree  ------------;;
degree_high = new((/10,dim_0(0),dim_0(1)/), float)
degree_low = new((/10,dim_0(0),dim_0(1)/), float)
f1 = addfile("data/Networks_Tmax90_H500Positive_lag0_land_Significant99.nc","r")
networks1 = (f1->networks0)
printVarSummary(networks1)


degree_high(0,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(50-10:60-10,70:80,:,:), 1), 0)
degree_high(2,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(50-10:60-10,100:110,:,:), 1), 0)
degree_high(4,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(50-10:60-10,130:140,:,:), 1), 0)
degree_high(6,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(50-10:60-10,240:250,:,:), 1), 0)
degree_high(8,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(50-10:60-10,285:295,:,:), 1), 0)

degree_high(1,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(30-10:40-10,60:70,:,:), 1), 0)
degree_high(3,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(30-10:40-10,90:100,:,:), 1), 0)
degree_high(5,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(30-10:40-10,110:120,:,:), 1), 0)
degree_high(7,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(30-10:40-10,240:250,:,:), 1), 0)
degree_high(9,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(30-10:40-10,270:280,:,:), 1), 0)
delete(networks1)



f2 = addfile("data/Networks_Tmax90_H500Negative_lag0_land_Significant99.nc","r")
networks2 = (f2->networks0)

degree_low(0,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks2(50-10:60-10,70:80,:,:), 1), 0)
degree_low(2,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks2(50-10:60-10,100:110,:,:), 1), 0)
degree_low(4,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks2(50-10:60-10,130:140,:,:), 1), 0)
degree_low(6,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks2(50-10:60-10,240:250,:,:), 1), 0)
degree_low(8,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks2(50-10:60-10,285:295,:,:), 1), 0)

degree_low(1,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks2(30-10:40-10,60:70,:,:), 1), 0)
degree_low(3,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks2(30-10:40-10,90:100,:,:), 1), 0)
degree_low(5,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks2(30-10:40-10,110:120,:,:), 1), 0)
degree_low(7,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks2(30-10:40-10,240:250,:,:), 1), 0)
degree_low(9,:,:) = dim_sum_n_Wrap(dim_sum_n_Wrap(networks2(30-10:40-10,270:280,:,:), 1), 0)
delete(networks2)


degree_high!1 = "lat"
degree_high!2 = "lon"
degree_high&lat = lat
degree_high&lon = lon
copy_VarCoords(degree_high, degree_low)











;;------------  to nc file  ------------;;
output="data/"
name1="Degree8"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->degree_high=degree_high
ncdf1->degree_low=degree_low


f3 = addfile("data/Degree8.nc","r")
degree_high = (f3->degree_high)/1000.0
degree_low = (f3->degree_low)/1000.0



degree_high!1 = "lat"
degree_high!2 = "lon"
degree_high&lat = lat
degree_high&lon = lon
copy_VarCoords(degree_high, degree_low)
printMinMax(degree_high, 1)






;-------------- distance -------------;;
distance = new((/10,dim_0(0),dim_0(1)/), double)
f0 = addfile("data/distance_i_j.nc","r")
printVarSummary(f0->distance(55,75,10:85,:))
printVarSummary(distance(0,:,:))
distance(0,:,:) = f0->distance(55,75,10:85,:)
distance(2,:,:) = f0->distance(55,105,10:85,:)
distance(4,:,:) = f0->distance(55,135,10:85,:)
distance(6,:,:) = f0->distance(55,245,10:85,:)
distance(8,:,:) = f0->distance(55,290,10:85,:)

distance(1,:,:) = f0->distance(35,65,10:85,:)
distance(3,:,:) = f0->distance(35,95,10:85,:)
distance(5,:,:) = f0->distance(35,115,10:85,:)
distance(7,:,:) = f0->distance(35,245,10:85,:)
distance(9,:,:) = f0->distance(35,275,10:85,:)


copy_VarCoords(degree_high, distance)






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
res@tmXBLabelFontHeightF = 0.018
res@tmYLLabelFontHeightF = 0.018
res@tmXBTickSpacingF = 45.0
res@tmYLTickSpacingF = 20.0




res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmXBMajorLengthF = 0.0012
res@tmYLMinorLengthF = 0.0012
res@tmXBMinorLengthF = 0.009
res@tmYLMinorLengthF = 0.009
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.012
res@tmYLMajorOutwardLengthF = 0.012
res@tmXBMinorOutwardLengthF = 0.009
res@tmYLMinorOutwardLengthF = 0.009





res@mpMinLatF = 10
res@mpMaxLatF = 85
res@mpMinLonF = 0
res@mpMaxLonF = 360
res@mpCenterLonF = 180-165-3

res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/) + 165+3
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
res@tmXBMinorOn = True
res@tmXBMinorValues = fspan(0, 360, 25)+3

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
lires@xyLineThicknesses  = (/4.5,8.0,8.0/)
lires@xyLineColors       = (/"red","blue","orange"/)
lires@xyDashPatterns     = (/0,4,2/)

lires@tmYLMajorThicknessF = 2.5
lires@tmYRMajorThicknessF = 2.5
lires@tmXBMajorThicknessF = 2.5
lires@tmYLMinorThicknessF = 1.5
lires@tmYRMinorThicknessF = 1.5
lires@tmXBMinorThicknessF = 1.5
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











;;-----------------------------------------------------
pltType = "eps"
pltName = "data/FigureS4_high_degree8"
wks = gsn_open_wks(pltType, pltName)




;;--------------------- figures (c,f) degree inside western and eastern Europe ----------------------
res@cnFillPalette = "WhiteYellowOrangeRed" ;"WhiteYellowOrangeRed"
res@cnLevels     := (/7.0,12,17/)
res@cnFillColors := (/-1,60,140,215/)
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 0.8

do i=0,9
plot(i) = gsn_csm_contour_map_ce(wks, degree_high(i,:,:), res)

contour2(i) = gsn_csm_contour(wks, distance(i,:,:), sres)    
overlay(plot(i), contour2(i))
end do





;;------------------------ boxes -----------------------
x11 = (/70,60,100,90,130,110,240,240,285,270/)
y11 = (/50,30,50,30,50,30,50,30,50,30/)
dum11 = new(14, "graphic") 
do i=0,9
y1=y11(i)
y2=y1+10
x1=x11(i)
x2=x1+10
x=(/x1,x2,x2,x1,x1/)  ;; eastern Europe
y=(/y2,y2,y1,y1,y2/)
   lnres   =    True
   lnres@gsLineThicknessF=1.2
   lnres@gsnFrame = False

   lnres@gsLineColor="green4"
   lnres@gsLineDashPattern=0
   dum11(i)=gsn_add_polyline(wks,plot(i),x,y,lnres)
end do





pres = True
pres@gsnFrame = False
pres@gsnPanelLabelBar = True
pres@lbOrientation        = "Vertical"
pres@pmLabelBarWidthF     = 0.04
pres@pmLabelBarHeightF    = 0.3
pres@pmLabelBarParallelPosF = 0.11
pres@pmLabelBarOrthogonalPosF = 0.0105
pres@lbLabelFontHeightF  = 0.0115

pres@gsnPanelTop  = 0.95
pres@gsnPanelBottom = 0.05
pres@gsnPanelLeft  = 0.01
pres@gsnPanelRight = 0.99
pres@amJust  ="BottomRight"
pres@gsnPanelFigureStringsFontHeightF = 0.0115
pres@gsnPanelFigureStrings = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)"/)
gsn_panel(wks,plot(0:9),(/7,2/),pres)



txres               = True
txres@txFontHeightF = 0.016
txres@txFont = "helvetica"
txres@txFontColor = "black"
gsn_text_ndc(wks, "High-pressure   cross-degree   centrality   inside   regional   heatwaves",  0.485, 0.905, txres)

frame(wks)
   delete(wks)  
   system("convert -geometry 1800x1800 -density 1200x1200 -trim " + pltName + "."+pltType + " " + pltName + ".png")


end




