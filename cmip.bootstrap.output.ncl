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
name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\;
            "CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM",\;
            "CNRM-CM5-2",  "CNRM-CM5","CSIRO-Mk3-6-0",  \
             "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
             "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES",\
             "MIROC5","MIROC-ESM-CHEM",\
              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);;all
NLAT = 90;40;70;5;5;
SLAT = -90;-40;-70;-5;-5;
WLON = 0;100; 190 ;
ELON = 360;300;240 ;
  ymstrt= 187001  
  ymlast= 200412
  yrstrt= 1870  
  yrlast= 2004
nm = dimsizes(name_model)
nyr = yrlast-yrstrt+1

LEV = 0
LEVEL = (/1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10/)

sign = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)" ,"(l)"/)

volyear_tp = (/1883,1902,1963,1982,1991/);1896(sh) (1968),1886,1968,

volmon_tp = (/188301,190201,196301,198201,199101/);1896(sh) (1968),196801,188601,
volyear = volyear_tp;ELyear;
volmon = volmon_tp;EL;+100
nv = dimsizes(volyear)

n = nv

MO = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
SEASON = (/"DJF","MAM", "JJA","SON" /)

f0 = addfile(diri0+"sst.mnmean.v4.1854-2016.nc","r")
TIME0=f0->time
yyyymm0=cd_calendar(TIME0, -1)
ymstrt0 =ind(yyyymm0.eq.ymstrt)
ymlast0 = ind(yyyymm0.eq.ymlast)
ss = f0->sst(ymstrt0:ymlast0,{SLAT:NLAT},{WLON:ELON})
ntime = dimsizes(ss&time)


ny1 = 5 ; year before year0
ny2 = 5 ; year after year0
nyr = ny1+ny2
nmon = ny1*12+ny2*12
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
delete([/f,files,TIME,yyyymm,ymstrt1,ymlast1/]);

nynino = 2004-1870+1
nino3 =  asciiread(diri+"CMIP5.nino34.unnorm.1870-2004.txt", (/nm*nynino*12/), "float")
nino3_a = onedtond(nino3, (/nm,nynino*12/))
nino3_a!0 = "number"
nino3_a!1 = "time"

xsea =   new((/nm, ndim(0)/12, 17,ndim(1)/), typeof(x0), x0@_FillValue)
;xseaj =   new((/nm, ndim(0)/12, ndim(1),ndim(2)/), typeof(x0), x0@_FillValue)
;--------------------------
do ik =0,nm-1
;  if (ik.eq.13)then 
;    files = systemfunc(" ls -1 "+diri3+"psl/"+name_model(ik)+"/r1i1p1/*.nc")
;  else
;  if (ik.eq.12)then 
;    files = systemfunc(" ls -1 "+diri3+"hus/"+name_model(ik)+"/r1i1p1/*.nc")
;  else
 ; if (ik.eq.11)then 
 ;   files = systemfunc(" ls -1 "+diri3+"va/"+name_model(ik)+"/r1i1p1/*.nc")
 ; else
;files = systemfunc(" ls -1 "+diri+"rsds/rsds_Amon_"+name_model(ik)+"*.nc")
;files = systemfunc(" ls -1 "+diri2+"pr/"+name_model(ik)+".pr.historical.nc");.mon
;files = systemfunc(" ls -1 "+diri2+"hgt/"+name_model(ik)+".hgt.historical.mon.nc") ;zg
files = systemfunc(" ls -1 "+diri3+"ta/"+name_model(ik)+"/r1i1p1/*.nc")

;files = systemfunc(" ls -1 "+diri+"hus/hus_*"+name_model(ik)+".nc")
;files = systemfunc(" ls -1 "+diri2+"slp/"+name_model(ik)+".slp.historical.mon.nc");psl
;files = systemfunc(" ls -1 "+diri2+"tas/"+name_model(ik)+".tas.historical.mon.nc")
;files = systemfunc(" ls -1 "+diri2+"tos/"+name_model(ik)+".tos.historical.mon.nc")
;files = systemfunc(" ls -1 "+diri2+"ua/"+name_model(ik)+".ua.historical.nc")
;files = systemfunc(" ls -1 "+diri2+"wap/"+name_model(ik)+"*.nc")

;print(files)
;files = systemfunc(" ls -1 "+diri+"va/va_regrid"+name_model(ik)+".nc")
;end if 
;files = systemfunc(" ls -1 "+diri3+"zos_Omon_"+name_model(ik)+"_historical*.nc")

;(0) 100000(1) 92500(2) 85000(3) 70000(4) 60000(5) 50000(6) 40000(7) 30000(8) 25000(9) 20000(10)  15000(11)  10000(12)  7000(13)  5000(14)  3000(15)  2000(16)  1000
LEV = 13
f   = addfiles(files,"r")
TIME=f[:]->time
yyyymm=cd_calendar(TIME, -1)
ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)
print(ymstrt1)
print(ymlast1)
x  = f[:]->ta(ymstrt1(0):ymstrt1(0)+1619,0:16,{SLAT:NLAT},{WLON:ELON});LEV,:,:);ymstrt1(0)+1619;ymlast1(0)
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
;xx = xx/100. ;slp Pa-->hPa  (name:psl)
delete([/f,files,TIME,x,yyyymm,ymstrt1,ymlast1,lat,lon/]);

xx= dtrend_msg_n(xx&time,xx,False,False,0)
xx1 = dim_avg_n_Wrap(xx, 3)

;ndim = dimsizes(xx)
;sc1 = clmMonTLL(sst1)
;sc2 = clmMonTLL(sst2)
xc  = clmMonTLL(xx1)
xa  = calcMonAnomTLL(xx1, xc);(ik,:,:,:,:)

;----------------filter----------------------PAY ATTENTION! This is month data, not years!
xf = xa;(ik,:,:,:,:);rem_sa;
;uf = ua(ik,:,:,:)
;vf = va(ik,:,:,:)
ihp   = 1       ; low-pass filter: ihp = 0; high-pass ihp = 1; band-pass ihp = 2.
sigma = 1.0                           ; Lanczos sigma
nWgt  = 51 ;必须是奇数 表示权重系数的总个数，越大效果越好，但是两端损失越多
fca   = 1./(9*12) ;                        
fcb   = -999;1./(2*12);1./8. *12  ;8——50年滤波
wgt  = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
xf = wgt_runave_n (xa, wgt, 0, 0 ) ;xa(ik,:,:,:)最后一个是指sst1第0维
xsea(ik,:,:,:) = month_to_season(xf, "DJF")
;xseaj(ik,:,:,:) = month_to_season(xf, "JJA")
end do 
printVarSummary(xsea)
;----------------original bootstrap----------------------------
delete(ndim)
ndim = dimsizes(xsea)
  nBoot   = 10000                               ; This could be time consuming
                                               ; Test with small nBoot (eg: nBoot=10)
  xsea1 = reshape(xsea, (/ndim(0)*ndim(1),ndim(2),ndim(3)/))
  xBt   = new ((/3,nBoot,  ndim(2),ndim(3)/), typeof(x0), x0@_FillValue)         ; array to hold nBoot results 

nsort = (/67,39,33/) ;type A; type B; type C
do ia = 0,2
  do ns=0,nBoot-1                              ; generate bootstrap results
           iw = random_uniform(0, ndim(0)*ndim(1)-1, nsort(ia)) ;第一年和最后三年废掉了，因为要提取前一年以及后三年
           jj = round(iw, 3)
           xBt(ia,ns,:,:)= dim_avg_n_Wrap(xsea1(jj,:,:),0)
  end do
  delete(iw)
  delete(jj)
end do

xBootLow =   new((/3,  ndim(2),ndim(3)/), typeof(x0), x0@_FillValue)
xBootHi =   new((/3,  ndim(2),ndim(3)/), typeof(x0), x0@_FillValue)
do ia = 0,2
do ix = 0,ndim(2)-1
  do iy = 0,ndim(3)-1

  xBoot = xBt(ia,:,ix,iy);, ndtooned((/nBoot*nyra/))
                                               ; mean of all bootstrap samples at each grid point
;at =  ind(.not.ismissing(xBoot))
nis =  num(ismissing(xBoot))

  idx  = dim_pqsort_n(xBoot, 2, 0) ;(at) sort bootstrap means into ascending order at each grid point
  i025 = toint(0.05*(nBoot-nis));*nyra
  i975 = toint(0.95*(nBoot-nis));*nyra
  xBootLow(ia,ix,iy)   = xBoot(i025+nis-1)           ;  2.5% level; (nlat,mlon)
  xBootHi(ia,ix,iy)   = xBoot(i975+nis-1)           ; 97.5% level; (nlat,mlon)
end do 
end do 
end do 

  units = xsea@units
  xBootLow@long_bame = "Bootstrap  2.5% Mean"
  xBootLow@units     = units
  xBootHi@long_bame = "Bootstrap 97.5% Mean"
  xBootHi@units     = units
;  xBoot_Range@long_name = "95% anomaly range"
;  xBoot_Range@units     = units

  copy_VarCoords(xsea(0:2,0,:,:), xBootLow(:,:,:))  
  copy_VarCoords(xsea(0:2,0,:,:), xBootHi(:,:,:))

  xBootLow!0 = "sort"
  xBootHi!0 = "sort"

printVarSummary(xBootLow)
;-------------output---------------------------
fout_name = "cmip.ta.allev.bootstrap.nc"
   system("/bin/rm -f "+diri+fout_name)   ; remove any pre-existing file
   fout = addfile(diri+fout_name ,"c")  ; open output netCDF file

   fAtt               = True            ; assign file attributes
       fAtt@title         = "bootstrap 90% confidence level in CMIP5"
;       fAtt@source_file   =  "original-file.nc"
       fAtt@Conventions   = "None"
       fAtt@creation_date = systemfunc ("date")
       fileattdef( fout, fAtt )            ; copy file attributes
  ;    filedimdef(fout,"time",-1,True) 

    fout->xBootLow = xBootLow
    fout->xBootHi = xBootHi

  delete(fout)


end