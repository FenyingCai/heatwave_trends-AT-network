load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin
;;------------  读取数据  ------------;;
f00 = addfile("/home/ys17-23/Extension/cai_fy/data0/topography/topo_adapt_era5_air2m.nc","r")
topo = f00->topo({10:85},:)

f0 = addfile("/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/2_networks/part1_heatheat_daily/Networks_Tmax90_lag0_part1_heatheat_44yrs.nc","r")
networks0 = f0->hw_networks0(:,:,:,:)
printMinMax(networks0, 0)


;;------------  mask掉海洋  ------------;;
topo_networks1 = conform(networks0, topo, (/0,1/))
topo_networks2 = conform(networks0, topo, (/2,3/))

networks0 = where(topo_networks1.ge.0.0, networks0, 1.0e20)
networks0 = where(topo_networks2.ge.0.0, networks0, 1.0e20)
networks0@_FillValue = 1.0e20
printMinMax(networks0, 0)



output="/home/ys17-23/cai_fy/data4_heatwaves/NCC2_1980_2010/2_networks/part1_heatheat_daily/"
name1="Networks_Tmax90_lag0_part1_heatheat_land"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->networks0=networks0
end

