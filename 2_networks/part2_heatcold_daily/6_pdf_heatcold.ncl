load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin
;;------------------------------------------------------------;;
;;--------------    读取热浪-冷事件相关系数网络    -------------;;
f1 = addfile("/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/part2_heatcold_daily/Networks_Tmax90_Tave10_lag0_part2_heatcold_land_Significant99.nc","r")
networks0 = f1->networks0
networks0 = where(networks0.gt.0.0, 1.0, 1.0e20)
networks0@_FillValue = 1.0e20
dim_0 = dimsizes(networks0)



;;------------  纬度加权矩阵  ------------;;
lat_matrix = new((/dim_0(0),dim_0(1), dim_0(0),dim_0(1)/), float)
do i = 0,dim_0(0)-1
  do j = 0,dim_0(0)-1
    lat_matrix(i,:,j,:) = cos((i+10+j+10) /2.0 /180.0 * 3.14159)
  end do
end do




;;------------  读取距离数据  ------------;;
f0 = addfile("/public/home/fcai/extreme1_AT/NC2_1980_2010/1_threshold/distance_i_j_land.nc","r")
distance = f0->distance(10:85,:,10:85,:)
dist = fspan(0.0,10200.0,52)


pdf = new(52, float)
do i =0,51

if (dist(i).eq.10200.0) then
networks0_new = where((distance.ge.dist(i)-100), networks0, 0.0)
pdf(i) = dim_sum(dim_sum(dim_sum(dim_sum(networks0_new*lat_matrix))))
end if
if (dist(i).lt.10200.0) then
networks0_new = where((distance.ge.dist(i)-100 .and. distance.lt.dist(i)+100), networks0, 0.0)
pdf(i) = dim_sum(dim_sum(dim_sum(dim_sum(networks0_new*lat_matrix))))
end if

delete(networks0_new)
print(" i="+i+", dist= "+dist(i)+", pdf = "+pdf(i))
end do

delete([/lat_matrix,distance,networks0/])

pdf = pdf/dim_sum(pdf)
print(dist)
print("----  heat-cold,  concurrent days > 99% significant  ----")
print(pdf)


output="/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/part2_heatcold_daily/"
name1="PDF_heatcold_daily"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->pdf=pdf


end






