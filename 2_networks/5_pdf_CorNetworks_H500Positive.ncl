load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin
;;------------------------------------------------------------;;
;;-------------------    read networks    --------------------;;
f1 = addfile("/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/CorNetorks_Tmax90th_H500Positive.nc","r")
cor_4D = f1->cor_4D
cor_4D@_FillValue = 1.0e20

cor_4D = where(cor_4D.ge.0.2973, 1.0, 0.0)    ;; 0.384 (p=0.01)  ;; 0.350 (p=0.02)   ;; 0.2973 (p=0.05)
printMinMax(cor_4D, 0)

;;------------  remove ocean  ------------;;
f00 = addfile("/public/home/fcai/data0/topography/topo_adapt_era5_air2m.nc","r")
topo = f00->topo({10:85},:)
topo_4D = conform(cor_4D, topo, (/0,1/))
cor_4D = where(topo_4D.ge.0.0, cor_4D, 1.0e20) 
cor_4D@_FillValue = 1.0e20
delete(topo_4D)



;;------------  lattiude weighted  ------------;;
lat_matrix = new((/76,360,76,360/), float)
do i = 0,75
  do j = 0,75
    lat_matrix(i,:,j,:) = cos((i+10+j+10) /2.0 /180.0 * 3.14159)
  end do
end do
cor_all = dim_sum(dim_sum(dim_sum(dim_sum(cor_4D*lat_matrix))))
print(" R all = "+ cor_all)




;;------------  read distance  ------------;;
f0 = addfile("/public/home/fcai/extreme1_AT/NC2_1980_2010/1_threshold/distance_i_j.nc","r")
distance = f0->distance(10:85,:,10:85,:)

dist = fspan(0.0,10200.0,52)

cor_pdf = new(52, float)

do i =0,51

if (dist(i).eq.10200.0) then
cor_4D_new = where((distance.ge.dist(i)-100), cor_4D, 0.0)
cor_pdf(i) = dim_sum(dim_sum(dim_sum(dim_sum(cor_4D_new*lat_matrix)))) / cor_all
end if
if (dist(i).lt.10200.0) then
cor_4D_new = where((distance.ge.dist(i)-100 .and. distance.lt.dist(i)+100), cor_4D, 0.0)
cor_pdf(i) = dim_sum(dim_sum(dim_sum(dim_sum(cor_4D_new*lat_matrix)))) / cor_all
end if

delete(cor_4D_new)
print(" i="+i+", dist= "+dist(i)+", pdf = "+cor_pdf(i))
end do

delete([/lat_matrix,distance,cor_4D/])


print(dist)
print("---- Correlation  heat - h500 Positive > 99% ----")
print(cor_pdf)




;;------------  to nc file  ------------;;
output="/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/"
name1="PDF_heathigh_interannualR"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->cor_pdf=cor_pdf



end
