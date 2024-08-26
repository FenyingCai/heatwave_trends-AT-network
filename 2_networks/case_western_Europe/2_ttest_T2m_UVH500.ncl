load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin
f00 = addfile("/public/home/fcai/data0/topography/topo_adapt_era5_air2m.nc","r")
topo0 = f00->topo({0:90},:)
lat = topo0&lat
lon = topo0&lon


;;------------  read regional heatwaves: T2m, Hgt, U, V  ------------;;
f1 = addfile("/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/case_western_Europe/era5_WE_composite_T2m_UVH500.nc", "r")
t_diff = f1->t_diff({0:90},:)
u_diff = f1->u_diff({0:90},:)
v_diff = f1->v_diff({0:90},:)
h_diff = f1->h_diff({0:90},:)
u_diff&lat = lat
u_diff&lon = lon
copy_VarCoords(u_diff, t_diff)
copy_VarCoords(u_diff, v_diff)
copy_VarCoords(u_diff, h_diff)


t_WE = f1->t_WE(:,{0:90},:)
u_WE = f1->u_WE(:,{0:90},:)
v_WE = f1->v_WE(:,{0:90},:)
h_WE = f1->h_WE(:,{0:90},:)
u_WE&lat = lat
u_WE&lon = lon
copy_VarCoords(u_WE, t_WE)
copy_VarCoords(u_WE, v_WE)
copy_VarCoords(u_WE, h_WE)

t_WE_clm = dim_avg_n_Wrap(t_WE, 0)
u_WE_clm = dim_avg_n_Wrap(u_WE, 0)
v_WE_clm = dim_avg_n_Wrap(v_WE, 0)
h_WE_clm = dim_avg_n_Wrap(h_WE, 0)
t_WE_size = equiv_sample_size(t_WE(lat|:,lon|:,day|:), 0.05, 0)
u_WE_size = equiv_sample_size(u_WE(lat|:,lon|:,day|:), 0.05, 0)
v_WE_size = equiv_sample_size(v_WE(lat|:,lon|:,day|:), 0.05, 0)
h_WE_size = equiv_sample_size(h_WE(lat|:,lon|:,day|:), 0.05, 0)
copy_VarCoords(u_diff, t_WE_size)
copy_VarCoords(u_diff, u_WE_size)
copy_VarCoords(u_diff, v_WE_size)
copy_VarCoords(u_diff, h_WE_size)
t_WE_var = dim_variance_n_Wrap(t_WE, 0)
u_WE_var = dim_variance_n_Wrap(u_WE, 0)
v_WE_var = dim_variance_n_Wrap(v_WE, 0)
h_WE_var = dim_variance_n_Wrap(h_WE, 0)
delete([/t_WE,u_WE,v_WE,h_WE/])




t = f1->t(:,{0:90},:)
u = f1->u(:,{0:90},:)
v = f1->v(:,{0:90},:)
h = f1->h(:,{0:90},:)
u&lat = lat
u&lon = lon
copy_VarCoords(u, t)
copy_VarCoords(u, v)
copy_VarCoords(u, h)
t_clm = dim_avg_n_Wrap(t, 0)
u_clm = dim_avg_n_Wrap(u, 0)
v_clm = dim_avg_n_Wrap(v, 0)
h_clm = dim_avg_n_Wrap(h, 0)
t_size = equiv_sample_size(t(lat|:,lon|:,day2|:), 0.05, 0)
u_size = equiv_sample_size(u(lat|:,lon|:,day2|:), 0.05, 0)
v_size = equiv_sample_size(v(lat|:,lon|:,day2|:), 0.05, 0)
h_size = equiv_sample_size(h(lat|:,lon|:,day2|:), 0.05, 0)
copy_VarCoords(u_diff, t_size)
copy_VarCoords(u_diff, u_size)
copy_VarCoords(u_diff, v_size)
copy_VarCoords(u_diff, h_size)
t_var = dim_variance_n_Wrap(t, 0)
u_var = dim_variance_n_Wrap(u, 0)
v_var = dim_variance_n_Wrap(v, 0)
h_var = dim_variance_n_Wrap(h, 0)
delete([/t,u,v,h/])










;;------------  t test   ------------;;
t_t = ttest(t_WE_clm, t_WE_var, t_WE_size, t_clm, t_var, t_size, False, False)
u_t = ttest(u_WE_clm, u_WE_var, u_WE_size, u_clm, u_var, u_size, False, False)
v_t = ttest(v_WE_clm, v_WE_var, v_WE_size, v_clm, v_var, v_size, False, False)
h_t = ttest(h_WE_clm, h_WE_var, h_WE_size, h_clm, h_var, h_size, False, False)

copy_VarCoords(u_diff, t_t)
copy_VarCoords(u_diff, u_t)
copy_VarCoords(u_diff, v_t)
copy_VarCoords(u_diff, h_t)



output="/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/case_western_Europe/"
name1="ttest_WE_T2m_UVH500"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->t_t=t_t
ncdf1->u_t=u_t
ncdf1->v_t=v_t
ncdf1->h_t=h_t


end




