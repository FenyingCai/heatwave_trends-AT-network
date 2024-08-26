
import numpy as np
import xarray as xr

##--##--##--##--  read data   --##--##--##--##
ds2 = xr.open_dataset("/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/case_eastern_Europe/era5_UVTH_clm_5levels.nc")
h_clm = ds2.h_clm
u_clm = ds2.u_clm
v_clm = ds2.v_clm
t_clm = ds2.t_clm
level = h_clm.level;  lat = h_clm.lat;  lon = h_clm.lon

ds3 = xr.open_dataset("/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/case_eastern_Europe/era5_EE_composite_UVH_5levels.nc")
h_diff = ds3.h_diff

print(" h_clm ", np.nanmax(h_clm), np.nanmin(h_clm))
print(" h_diff ", np.nanmax(h_diff), np.nanmin(h_diff))






##--##--##--##--  TN wave activvity flux   --##--##--##--##
###### constant ######
g     = 9.806650   # gravitational acceleration
gc    = 290.0      # gas constant
radius= 6371009.0  # Radius of Earth
sclhgt= 8000.0     # atmospheric scale height
omega = 7.292e-5   # angular velocity


###### 3D (level,lat,lon), DataArray
def WAF_flux(Tc,Uc,Vc,GEOc,GEOa):

    data_shape =Tc.shape
    data_coords=Tc.coords
    UVc=np.sqrt(Uc**2+Vc**2)


    Uc  =xr.where(abs(Uc  ['lat'])<=20,np.nan,Uc  ).transpose('level','lat','lon')
    Vc  =xr.where(abs(Vc  ['lat'])<=20,np.nan,Vc  ).transpose('level','lat','lon')
    Tc  =xr.where(abs(Tc  ['lat'])<=20,np.nan,Tc  ).transpose('level','lat','lon')
    UVc =xr.where(abs(UVc ['lat'])<=20,np.nan,UVc ).transpose('level','lat','lon')
    GEOc=xr.where(abs(GEOc['lat'])<=20,np.nan,GEOc).transpose('level','lat','lon')
    GEOa=xr.where(abs(GEOa['lat'])<=20,np.nan,GEOa).transpose('level','lat','lon')

    lon=np.array(Tc['lon'  ])[np.newaxis,np.newaxis,:         ]
    lat=np.array(Tc['lat'  ])[np.newaxis,:         ,np.newaxis]
    pp =np.array(Tc['level'])[:         ,np.newaxis,np.newaxis]


    Uc  =np.array(Uc  )
    Vc  =np.array(Vc  )
    Tc  =np.array(Tc  )
    UVc =np.array(UVc )
    GEOc=np.array(GEOc)
    GEOa=np.array(GEOa)


    Uc  =np.where(Uc>=0,Uc  ,np.nan)
    Vc  =np.where(Uc>=0,Vc  ,np.nan)
    Tc  =np.where(Uc>=0,Tc  ,np.nan)
    UVc =np.where(Uc>=0,UVc ,np.nan)
    GEOc=np.where(Uc>=0,GEOc,np.nan)
    GEOa=np.where(Uc>=0,GEOa,np.nan)

    ##### Finite difference
    dlon  =np.deg2rad(np.gradient(lon,axis=2))
    dlat  =np.deg2rad(np.gradient(lat,axis=1))
    coslat=np.array(np.cos(np.deg2rad(lat)))
    sinlat=np.array(np.sin(np.deg2rad(lat)))
    if not data_shape[0]==1:
        dlev=np.gradient(-sclhgt*np.log(pp/1000.0),axis=0)

    f=2*omega*sinlat

    ###### N^2
    if not data_shape[0]==1:
        N2=np.array(gc*(pp/1000.0)**0.286)/sclhgt*np.gradient(Tc*(1000.0/pp)**0.286,axis=0)/\
            (np.gradient(-sclhgt*np.log(pp/1000.0),axis=0))

    ## PSI
    PSIa=GEOa/f

    ##### Finite difference
    dzdlon=np.gradient(PSIa,axis=2)/dlon
    dzdlat=np.gradient(PSIa,axis=1)/dlat

    ddzdlonlon=np.gradient(dzdlon,axis=2)/dlon
    ddzdlonlat=np.gradient(dzdlon,axis=1)/dlat
    ddzdlatlat=np.gradient(dzdlat,axis=1)/dlat

    if not data_shape[0]==1:
        dzdlev    =np.gradient(PSIa  ,axis=0)/dlev
        ddzdlonlev=np.gradient(dzdlon,axis=0)/dlev
        ddzdlatlev=np.gradient(dzdlat,axis=0)/dlev

    xuterm=dzdlon*dzdlon-PSIa*ddzdlonlon
    xvterm=dzdlon*dzdlat-PSIa*ddzdlonlat

    yuterm=xvterm
    yvterm=dzdlat*dzdlat-PSIa*ddzdlatlat

    if not data_shape[0]==1:
        zuterm=dzdlon*dzdlev-PSIa*ddzdlonlev
        zvterm=dzdlat*dzdlev-PSIa*ddzdlatlev

    coef=pp*coslat/1000.0/2.0/UVc
    Fx=coef*(xuterm*Uc/(radius*coslat)**2+xvterm*Vc/(radius**2*coslat))
    Fy=coef*(yuterm*Uc/(radius**2*coslat)+yvterm*Vc/radius**2)
    if not data_shape[0]==1:
        Fz=coef*(f**2/N2*(zuterm*Uc/(radius*coslat)+zvterm*Vc/radius))

    ###### to dataarray
    Fx=xr.DataArray(
        Fx,
        dims  =('level','lat','lon'), coords=data_coords
    )
    Fy=xr.DataArray(
        Fy,
        dims  =('level','lat','lon'), coords=data_coords
    )
    if not data_shape[0]==1:
        Fz=xr.DataArray(
            Fz,
            dims  =('level','lat','lon'), coords=data_coords
        )

    ### return Fx, Fy
    if not data_shape[0]==1:
        return Fx, Fy
    else:
        return Fx, Fy





F0 = WAF_flux(t_clm, u_clm, v_clm, h_clm, h_diff)
Fx0=F0[0];  Fy0=F0[1]


##--##--##--   to nc file   --##--##--##
ds8 = xr.Dataset(data_vars=dict(Fx=Fx0, Fy=Fy0))
ds8.to_netcdf("/public/home/fcai/extreme1_AT/NC2_1980_2010/2_networks/case_eastern_Europe/era5_EE_WAF.nc")
ds8.close()




