load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin
;;------------  读取数据  ------------;;
f0 = addfile("/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/2_networks/part1_heatheat_daily/Networks_Tmax90_lag0_part1_heatheat_land.nc","r")
networks0 = f0->networks0(:,:,:,:)
print(dim_avg(dim_avg(dim_avg(dim_avg(networks0)))))
dim_0 = dimsizes(networks0)



;;------------  mask 掉不显著 links  ------------;;
f1 = addfile("/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/2_networks/part1_heatheat_daily/shuffle99_Tmax90_lag0_part1_heatheat.nc","r")
matrix4D_99 = f1->matrix4D_99  ;;??
; printVarSummary(matrix4D_99)

networks0 = where(networks0.ge.matrix4D_99, networks0, 1.0e20)
networks0@_FillValue = 1.0e20
print(dim_avg(dim_avg(dim_avg(dim_avg(networks0)))))


;;------------  转换成01矩阵  ------------;;
networks01 = new((/dim_0(0),dim_0(1),dim_0(2),dim_0(3)/), integer)

networks01 = where(networks0.ge.0.0, 1, 0)
; networks01@_FillValue = 1.0e20
printVarSummary(networks01)
printMinMax(networks01, 0)
delete(networks0)


output="/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/2_networks/part1_heatheat_daily/"
name1="Networks_Tmax90_lag0_part1_heatheat_land_Significant99_01"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->networks01=networks01


end

