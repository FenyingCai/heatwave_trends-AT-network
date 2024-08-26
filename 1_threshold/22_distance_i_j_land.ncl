load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

;-------------- read togpography -------------;;
f00 = addfile("/public/home/fcai/data0/topography/topo_adapt_era5_air2m.nc","r")
topo0 = f00->topo({0:90},{0:360})


;-------------- read distance -------------;;
f0 = addfile("/public/home/fcai/extreme1_AT/NC2_1980_2010/5_response/11_distance_test/distance_i_j.nc","r")
distance = f0->distance(:,:,:,:)
printMinMax(distance, 1)


topo1 = conform(distance, topo0, (/0,1/))
distance = where(topo1.gt.0.0, distance, 1.0e20)
distance@_FillValue = 1.0e20
delete(topo1)
topo2 = conform(distance, topo0, (/2,3/))
distance = where(topo2.gt.0.0, distance, 1.0e20)
distance@_FillValue = 1.0e20

printMinMax(distance, 1)




;-------------- to nc file -------------;;
output="/public/home/fcai/extreme1_AT/NC2_1980_2010/1_threshold/"
system("/bin/rm -f "+output+"distance_i_j_land.nc")
ncdf=addfile(output+"distance_i_j_land.nc","c")
ncdf->distance=distance


end

