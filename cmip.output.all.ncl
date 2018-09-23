load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
; Nino 1+2 (0-10S, 90W-80W), Nino 3 (5N-5S, 150W-90W), 
; Nino 3.4/ONI (5N-5S, 170W-120W) and Nino 4 (5N-5S, 160E-150W)
;(0) 100000 (1) 92500 (2) 85000 (3) 70000 (4) 60000 (5) 50000 (6) 40000 (7) 30000 (8) 25000 (9) 20000 (10)  15000 (11)  10000 (12)  7000 (13)  5000 (14)  3000 (15)  2000 (16)  1000 
begin
;diri = "/nuist/p/work/feil/xc/"
diri = "/data2/home/liufeiclim2014/xing/cmip_vol/"
diri2 = "/data2/data/cmip5/historical_regrid/atmos/mon/"
diri3 = "/data2/data/cmip5/historical/atmos/mon/"
diri0 = "/data2/home/liufeiclim2014/xing/"
;name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\
;       "CESM1-BGC", "CESM1-FASTCHEM","CESM1-WACCM","CNRM-CM5-2",  \; "CESM1-CAM5-1-FV2",va
;      "CNRM-CM5",  \;,"CSIRO-Mk3-6-0" va
;        "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
;        "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES","MIROC5","MIROC-ESM-CHEM",\
;         "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\;"bcc-csm1-1-m","bcc-csm1-1",
name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\;
            "CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM",\;
            "CNRM-CM5-2",  "CNRM-CM5","CSIRO-Mk3-6-0",  \
             "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
             "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES",\
             "MIROC5","MIROC-ESM-CHEM",\
              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);;all
;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\
;            "CNRM-CM5","CNRM-CM5-2","CSIRO-Mk3-6-0", \;  
;             "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC",\;
;             "HadGEM2-ES",\;"HadGEM2-AO","GISS-E2-H","HadCM3","HadGEM2-CC","MIROC-ESM","MRI-CGCM3",
;              "MIROC5","MIROC-ESM-CHEM",\
;              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P",\
;              "MRI-ESM1","NorESM1-M","NorESM1-ME"/); rsdscs
NLAT = 50;90;70;5;5;
SLAT = 30;-90;-70;-5;-5;
WLON = 0;0; 190 ;
ELON = 360;360;240 ;
  ymstrt= 187001  
  ymlast= 200412
  yrstrt= 1870  
  yrlast= 2004
nm = dimsizes(name_model)
nyr = yrlast-yrstrt+1

volyear = (/1883,1902,1963,1982,1991/);

LEV = 3
LEVEL = (/1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10/)

sign = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)" ,"(l)"/)


MO = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
SEASON = (/"DJF","MAM", "JJA","SON" /)

f0 = addfile(diri0+"sst.mnmean.v4.1854-2016.nc","r")
TIME0=f0->time
yyyymm0=cd_calendar(TIME0, -1)
ymstrt0 =ind(yyyymm0.eq.ymstrt)
ymlast0 = ind(yyyymm0.eq.ymlast)
ss = f0->sst(ymstrt0:ymlast0,{SLAT:NLAT},{WLON:ELON})
ntime = dimsizes(ss&time)


ny1 = 2 ; year before year0
ny2 = 2 ; year after year0
nyr = ny1+ny2
nmon = ny1*12+ny2*12
;nyear = (/"year(-5)","year(-4)","year(-3)","year(-2)","year(-1)","year(0)","year(+1)",\
;         "year(+2)","year(+3)","year(+4)","year(+5)","year(+6)","year(+7)"/);"year(-6)",
files = systemfunc(" ls -1 "+diri2+"tas/"+name_model(0)+".*.nc")
;files = systemfunc(" ls -1 "+diri3+"zos_Omon_"+name_model(0)+"_historical*.nc")
f   = addfiles(files,"r") 
TIME=f[:]->time
yyyymm=cd_calendar(TIME, -1)
ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)
print(ymstrt1)
print(ymlast1)
x0  = f[:]->tas(ymstrt1(0):ymlast1(0),{SLAT:NLAT},{WLON:ELON});:,:)
nlat = dimsizes(x0&lat)
nlon = dimsizes(x0&lon)
ndim = dimsizes(x0)

;xa =  new((/nm,ndim(0),ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
;xavol =   new((/nm,n,  nmon,17,ndim(1)/), typeof(x0), x0@_FillValue) 
;xseavol = new((/nm,n, 4,nyr,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
xa_p =new((/nm,ndim(0),17,ndim(2)/), typeof(x0), x0@_FillValue) 
;--------------------------
delete(files)
do ik =0,nm-1
;files = systemfunc(" ls -1 "+diri2+"hgt/"+name_model(ik)+".hgt.historical.mon.nc") ;zg
;files = systemfunc(" ls -1 "+diri3+"ta/"+name_model(ik)+"/r1i1p1/*.nc")
files = systemfunc(" ls -1 "+diri2+"ua/"+name_model(ik)+".ua.historical.nc") 

;(0) 100000(1) 92500(2) 85000(3) 70000(4) 60000(5) 50000(6) 40000(7) 30000(8) 25000(9) 20000(10)  15000(11)  10000(12)  7000(13)  5000(14)  3000(15)  2000(16)  1000
;files = systemfunc(" ls -1 "+diri3+"zos_Omon_"+name_model(ik)+"_historical*.nc")


f   = addfiles(files,"r") 
TIME=f[:]->time
yyyymm=cd_calendar(TIME, -1)
ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)
print(ymstrt1)
print(ymlast1)
x  = f[:]->ua(ymstrt1(0):ymstrt1(0)+1619,0:16,{SLAT:NLAT},{WLON:ELON});:,:);ymlast1(0)

;--------------resolution--------------
lat = x&lat
lon = x&lon 
klon = ispan(0,nlon-1,1)*2.5+WLON  
klat = ispan(0,nlat-1,1)*2.5+SLAT
xx= linint2_Wrap(lon, lat, x, True, klon, klat, 0) 
    xx&LAT@units="degrees_north"
    xx&LON@units="degrees_east"

;xx  = x 
;xx = xx*24*3600 ;prate-->precip
;xx = xx/100. ;slp Pa-->hPa
delete([/f,files,TIME,x,yyyymm,ymstrt1,ymlast1,lat,lon/]);

xx= dtrend_msg_n(xx&time,xx,False,False,0)

;sc1 = clmMonTLL(sst1)
;sc2 = clmMonTLL(sst2)
xc  = clmMonTLLL(xx)
xa  = calcMonAnomTLLL(xx, xc)
xa_p(ik,:,:,:) = dim_avg_n_Wrap(xa, 2)
;delete([/f1,f2,files1,files2,TIME,yyyymm,ymstrt1,ymlast1/]);,x1,lat,lon

print(ik)
end do 
printVarSummary(xa_p)
delete(ndim)
ndim = dimsizes(xa_p)
;-------------output---------------------------
fout_name = "cmip.ua.alllev.30-50N.nc"
   system("/bin/rm -f "+diri+fout_name)   ; remove any pre-existing file
   fout = addfile(diri+fout_name ,"c")  ; open output netCDF file

   fAtt               = True            ; assign file attributes
       fAtt@title         = "ua anomaly in CMIP5 ensemble"
;       fAtt@source_file   =  "original-file.nc"
       fAtt@Conventions   = "None"
       fAtt@creation_date = systemfunc ("date")
       fileattdef( fout, fAtt )            ; copy file attributes
  ;    filedimdef(fout,"time",-1,True) 

    fout->xa = xa_p

  delete(fout)
  end 