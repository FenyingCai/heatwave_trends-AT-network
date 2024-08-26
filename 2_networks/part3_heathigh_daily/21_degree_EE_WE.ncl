load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin
;;-----------------  读取经纬度  ----------------;;
f00 = addfile("/public/home/fcai/data0/topography/topo_adapt_era5_air2m.nc","r")
topo0 = f00->topo({10:85},{0:360})
lat = topo0&lat
lon = topo0&lon
dim_0 = dimsizes(topo0)



;;------------  mask 掉不显著 links  ------------;;
f1 = addfile("/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/part3_heathigh_daily/Networks_Tmax90_H500Positive_lag0_land_Significant99.nc","r")
networks0 = f1->networks0
topo_4D = conform(networks0, topo0, (/0,1/))
networks0 = where(topo_4D.gt.0.0, networks0, 0.0)


degree_EE = dim_sum_n_Wrap(dim_sum_n_Wrap(networks0(45-10:65-10,25:65,:,:), 1), 0)   ;; eastern Europe
degree_EE!0 = "lat"
degree_EE!1 = "lon"
degree_EE&lat = lat
degree_EE&lon = lon
printVarSummary(degree_EE)
printMinMax(degree_EE, 1)


WE_lon = new(36, integer)
WE_lon(0:9) = ispan(350,359,1)
WE_lon(10:35) = ispan(0,25,1)
degree_WE = dim_sum_n_Wrap(dim_sum_n_Wrap(networks0(35-10:55-10, WE_lon,:,:), 1), 0)   ;; western Europe
copy_VarCoords(degree_EE, degree_WE)
printMinMax(degree_WE, 1)
delete(networks0)



;;------------  to nc file  ------------;;
output="/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/part3_heathigh_daily/"
name1="Degree_High_EasternEurope_WesternEurope"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->degree_EE=degree_EE
ncdf1->degree_WE=degree_WE


end




