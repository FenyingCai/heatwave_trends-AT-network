load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

dim_0 = (/76,360/)


network_44yrs = new((/dim_0(0),dim_0(1),dim_0(0),dim_0(1)/), float)
network_44yrs = 0.0

do iyear=0,43
print(" iyear = " + (iyear+1979))
;;------------  读取每年的 concurrence 矩阵  ------------;;
f1 = addfile("/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/2_networks/part2_heatcold_daily/Networks_Tmax90_Tave10_lag0_part2_heatcold_"+int2flt(iyear+1979)+".nc","r")
networks0 = doubletofloat(f1->networks0)

networks0@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(networks0))) then
      replace_ieeenan (networks0, value, 0)
      networks0@_FillValue = value
end if


network_44yrs = network_44yrs + networks0
delete(networks0)
printMinMax(network_44yrs, 0)
end do



output="/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/2_networks/part2_heatcold_daily/"
name1="Networks_Tmax90_Tave10_lag0_part2_heatcold_44yrs"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->networks0=network_44yrs


end

