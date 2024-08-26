load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

hgt = new((/44,92,4,91,360/), float)

;;-------------- read data --------------;;
do i = 0,43
print(" year = " + (i+1979))

f6 = addfile("data/zg."+flt2string(i+1979)+"-06.daily.nc", "r")
hgt(i,0:29,:,:,:) = short2flt(f6->z(:,{(/850,500,300,200/)},{0:90},:))
f7 = addfile("data/zg."+flt2string(i+1979)+"-07.daily.nc", "r")
hgt(i,30:60,:,:,:) = short2flt(f7->z(:,{(/850,500,300,200/)},{0:90},:))
f8 = addfile("data/zg."+flt2string(i+1979)+"-08.daily.nc", "r")
hgt(i,61:91,:,:,:) = short2flt(f8->z(:,{(/850,500,300,200/)},{0:90},:))
end do
lev = hgt&level


;;-------------- anomaly relative to 1981-2010 --------------;;
hgt_clm = dim_avg_n_Wrap(hgt(2:31,:,:,:,:), 0)   ;; 1981-2010
hgt_clm_5D = conform(hgt, hgt_clm, (/1,2,3,4/))
delete(hgt_clm)
hgt_ano = hgt - hgt_clm_5D
copy_VarCoords(hgt, hgt_ano)
delete(hgt_clm_5D)


;;-------------- 1981-2010 std --------------;;
hgt_std = dim_stddev_n_Wrap(hgt(2:31,:,:,:,:), 0)   ;; 1981-2010
hgt_std_5D = conform(hgt, hgt_std, (/1,2,3,4/))
delete(hgt_std)
copy_VarCoords(hgt, hgt_std_5D)
delete(hgt)


;;-------------- >1std --------------;;
hgt_Positive = where(hgt_ano.ge.hgt_std_5D, hgt_ano, 0.0)
copy_VarCoords(hgt_ano, hgt_Positive)
hgt_Positive@_FillValue=1e+20

hgt_Negative = where(hgt_ano.le.hgt_std_5D*(-1.0), hgt_ano, 0.0)
hgt_Negative@_FillValue=1e+20
copy_VarCoords(hgt_ano, hgt_Negative)
delete([/hgt_ano, hgt_std_5D/])

hgt_Positive_intensity = dim_sum_n_Wrap(hgt_Positive, 1)
hgt_Negative_intensity = dim_sum_n_Wrap(hgt_Negative, 1) * (-1.0)
copy_VarCoords(hgt_Positive_intensity, hgt_Negative_intensity)
printMinMax(hgt_Positive_intensity, 1)
printMinMax(hgt_Negative_intensity, 1)


hgt_Positive = where(hgt_Positive.gt.0.0, 1.0, 0.0)
hgt_Negative = where(hgt_Negative.lt.0.0, 1.0, 0.0)

hgt_Positive_day= dim_sum_n_Wrap(hgt_Positive, 1)
hgt_Negative_day = dim_sum_n_Wrap(hgt_Negative, 1)
copy_VarCoords(hgt_Positive_intensity, hgt_Positive_day)
copy_VarCoords(hgt_Positive_intensity, hgt_Negative_day)
delete([/hgt_Positive, hgt_Negative/])
printMinMax(hgt_Positive_day, 1)
printMinMax(hgt_Negative_day, 1)



output="data/"
system("/bin/rm -f "+output+"hgt_ano_44years_JJAyearly.nc")
ncdf=addfile(output+"hgt_ano_44years_JJAyearly.nc","c")
ncdf->hgt_Positive_intensity=hgt_Positive_intensity
ncdf->hgt_Negative_intensity=hgt_Negative_intensity

ncdf->hgt_Positive_day=hgt_Positive_day
ncdf->hgt_Negative_day=hgt_Negative_day

end

