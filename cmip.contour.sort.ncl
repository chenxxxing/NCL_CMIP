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
diri = "/nuist/p/work/feil/xc/"
;diri = "/data2/home/liufeiclim2014/xing/cmip_vol/"
;diri2 = "/data2/data/cmip5/historical_regrid/atmos/mon/"
;diri3 = "/data2/home/liufeiclim2014/xing/cmip_vol/zos/"
;diri0 = "/data2/home/liufeiclim2014/xing/"
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
;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\3,4,7,8,9,10,20,
;            "CNRM-CM5-2","CNRM-CM5","CSIRO-Mk3-6-0", \;  
;             "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC",\;
;             "HadGEM2-ES",\;"HadGEM2-AO","GISS-E2-H","HadCM3","HadGEM2-CC","MIROC-ESM","MRI-CGCM3",
;              "MIROC5","MIROC-ESM-CHEM",\
;              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P",\
;              "MRI-ESM1","NorESM1-M","NorESM1-ME"/); rsdscs
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

volyear_nh_gao = (/1912, 1925,1943,1976,1980/);NH 1875, 1907,1933,1956,
;volyear_sh_gao = (/1886 /);SH 把sh放入tp？！,1932,1991, 2011
;vol_sh_mag_gao = (/1.93/);gao
volyear_tp_gao = (/1883,1886,1902,1963,1982,1991/);1913,模式中貌似没有这个的火山强迫……

volyear_nh_sato = (/1891,1912/);NH (1891or92)
volyear_tp_sato = (/1883,1902,1963,1982,1991/);1896(sh) (1968),1886,1968,

volmon_tp_sato = (/188301,190201,196301,198201,199101/);1896(sh) (1968),196801,188601,
volyear = volyear_tp_sato;ELyear;
volmon = volmon_tp_sato;EL;+100
nv = dimsizes(volyear)

n = nv

MO = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
SEASON = (/"DJF","MAM", "JJA","SON" /)

f0 = addfile(diri+"sst.mnmean.v4.1854-2016.nc","r")
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
;nyear = (/"year(-5)","year(-4)","year(-3)","year(-2)","year(-1)","year(0)","year(+1)",\
;         "year(+2)","year(+3)","year(+4)","year(+5)","year(+6)","year(+7)"/);"year(-6)",
files = systemfunc(" ls -1 "+diri+"tas/"+name_model(0)+".*.nc")
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

xa =  new((/nm,ndim(0),ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
xavol =   new((/nm,n,  nmon,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
xseavol = new((/nm,n, 4,nyr,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
xa_p =new((/nm,ndim(0),ndim(2)/), typeof(x0), x0@_FillValue) 
;ua =  new((/nm,ndim(0),ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
;uavol =   new((/nm,n,  nmon,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
;useavol = new((/nm,n, 4,nyr,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
;ua_p =new((/nm,ndim(0),ndim(2)/), typeof(x0), x0@_FillValue) 
;va =  new((/nm,ndim(0),ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
;vavol =   new((/nm,n,  nmon,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
;vseavol = new((/nm,n, 4,nyr,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
;va_p =new((/nm,ndim(0),ndim(2)/), typeof(x0), x0@_FillValue) 
;--------------------------
do ik =0,nm-1
;files = systemfunc(" ls -1 "+diri+"rsds/rsds_Amon_"+name_model(ik)+"*.nc")
;files = systemfunc(" ls -1 "+diri+"pr/"+name_model(ik)+".pr.historical.nc");.mon
;files = systemfunc(" ls -1 "+diri+"tos/tos_Amon_"+name_model(0)+"_his*.nc")
;files = systemfunc(" ls -1 "+diri+"tas/"+name_model(ik)+".tas.historical.mon.nc")
files = systemfunc(" ls -1 "+diri+"ua/"+name_model(ik)+".ua.historical.nc")
;files2 = systemfunc(" ls -1 "+diri+"va/va_regrid"+name_model(ik)+".nc")

;(0) 100000(1) 92500(2) 85000(3) 70000(4) 60000(5) 50000(6) 40000(7) 30000(8) 25000(9) 20000(10)  15000(11)  10000(12)  7000(13)  5000(14)  3000(15)  2000(16)  1000
;files = systemfunc(" ls -1 "+diri3+"zos_Omon_"+name_model(ik)+"_historical*.nc")

LEV = 5

f   = addfiles(files,"r") 
TIME=f[:]->time
yyyymm=cd_calendar(TIME, -1)
ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)
print(ymstrt1)
print(ymlast1)
x  = f[:]->ua(ymstrt1(0):ymlast1(0),LEV,{SLAT:NLAT},{WLON:ELON});:,:);

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

ndim = dimsizes(xx)
;sc1 = clmMonTLL(sst1)
;sc2 = clmMonTLL(sst2)
xc  = clmMonTLL(xx)
xa(ik,:,:,:)  = calcMonAnomTLL(xx, xc)
;f1   = addfiles(files1,"r") 
;f2   = addfiles(files2,"r") 
;TIME=f1[:]->time
;yyyymm=cd_calendar(TIME, -1)
;ymstrt1 =ind(yyyymm.eq.ymstrt)
;ymlast1 = ind(yyyymm.eq.ymlast)
;;0.1000E+06  0.9250E+05  0.8500E+05  0.7000E+05  0.6000E+05  0.5000E+05
;;  0.4000E+05  0.3000E+05  0.2500E+05  0.2000E+05  0.1500E+05  0.1000E+05
;;   7000.       5000.       3000.       2000.       1000.    
;u  = f1[:]->ua(ymstrt1(0):ymlast1(0),2,{SLAT:NLAT},{WLON:ELON});:,:);
;v  = f2[:]->va(ymstrt1(0):ymlast1(0),2,{SLAT:NLAT},{WLON:ELON});:,:);

;delete([/f1,f2,files1,files2,TIME,yyyymm,ymstrt1,ymlast1/]);,x1,lat,lon

;u= dtrend_msg_n(u&time,u,False,False,0)
;uc  = clmMonTLL(u)
;ua(ik,:,:,:)  = calcMonAnomTLL(u, uc)
;v= dtrend_msg_n(v&time,v,False,False,0)
;;-------------;改变精度---------------------
;lat = v&lat
;lon = v&lon 
;klon = ispan(0,nlon-1,1)*2.5+WLON  
;klat = ispan(0,nlat-1,1)*2.5+SLAT
;vv= linint2_Wrap(lon, lat, v, True, klon, klat, 0) 
;    vv&LAT@units="degrees_north"
;    vv&LON@units="degrees_east"
;vc  = clmMonTLL(vv)
;va(ik,:,:,:)  = calcMonAnomTLL(vv, vc)
;delete(v)
;----------------filter----------------------PAY ATTENTION! This is month data, not years!
xf = xa(ik,:,:,:)
;uf = ua(ik,:,:,:)
;vf = va(ik,:,:,:)
ihp   = 1       ; low-pass filter: ihp = 0; high-pass ihp = 1; band-pass ihp = 2.
sigma = 1.0                           ; Lanczos sigma
nWgt  = 51 ;必须是奇数 表示权重系数的总个数，越大效果越好，但是两端损失越多
fca   = 1./(9*12) ;                        
fcb   = -999;1./(2*12);1./8. *12  ;8——50年滤波
wgt  = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
xf = wgt_runave_n (xa(ik,:,:,:), wgt, 0, 0 ) ;最后一个是指sst1第0维
;uf = wgt_runave_n (ua(ik,:,:,:), wgt, 0, 0 )
;vf = wgt_runave_n (va(ik,:,:,:), wgt, 0, 0 )

;------------normalize------------------------------
;xf = dim_standardize_n_Wrap(xf, 0, 0)
;---------------------------------------------
TIMEs2=xa&time
yyyymms2=cd_calendar(TIMEs2, -1)
ymvol2 = new(nv, integer)
do iv = 0,n-1
ymvol2(iv) = ind(yyyymms2.eq.volmon(iv));EL(iv));
xavol(ik,iv,:,:,:) = xf(ymvol2(iv)-ny1*12:ymvol2(iv)+ny2*12-1,:,:);ik,
;uavol(ik,iv,:,:,:) = uf(ymvol2(iv)-ny1*12:ymvol2(iv)+ny2*12-1,:,:);ik,
;vavol(ik,iv,:,:,:) = vf(ymvol2(iv)-ny1*12:ymvol2(iv)+ny2*12-1,:,:);ik,
;xseavol(ik,iv,:,:,:,:) = month_to_seasonN(xavol(ik,iv,:,:,:), SEASON)
end do 
xavol!0 = "model"
xavol!1 = "vol"
;uavol!0 = "model"
;uavol!1 = "vol"
;vavol!0 = "model"
;vavol!1 = "vol"
print(ymvol2)
;xa_p(ik,:,:) = dim_avg_n_Wrap(xa(:,{-5:5},:), 1)
print(ik)
end do 
printVarSummary(xavol)
;-------------remove nearby annual---------------------
xavolc = new((/nm,nv,12,ndim(1),ndim(2)/), typeof(xa), xa@_FillValue)
xavol1 = new((/nm,nv,nmon,ndim(1),ndim(2)/), typeof(xa), xa@_FillValue)
;uavolc = new((/nm,nv,12,ndim(1),ndim(2)/), typeof(xa), xa@_FillValue)
;uavol1 = new((/nm,nv,nmon,ndim(1),ndim(2)/), typeof(xa), xa@_FillValue)
;vavolc = new((/nm,nv,12,ndim(1),ndim(2)/), typeof(xa), xa@_FillValue)
;vavol1 = new((/nm,nv,nmon,ndim(1),ndim(2)/), typeof(xa), xa@_FillValue)

do i=0,11
xavolc(:,:,i,:,:) = dim_avg_n_Wrap(xavol(:,:,i::12,:,:), 2)
;uavolc(:,:,i,:,:) = dim_avg_n_Wrap(uavol(:,:,i::12,:,:), 2)
;vavolc(:,:,i,:,:) = dim_avg_n_Wrap(vavol(:,:,i::12,:,:), 2)
end do 
xavolc!0 = "sample"
xavol1 = xavol
;uavolc!0 = "sample"
;uavol1 = uavol
;vavolc!0 = "sample"
;vavol1 = vavol
do i = 0,11
do j = 0,nyr-1
xavol1(:,:,j*12+i,:,:)  = xavol(:,:,i+j*12,:,:)-xavolc(:,:,i,:,:)
;uavol1(:,:,j*12+i,:,:)  = uavol(:,:,i+j*12,:,:)-uavolc(:,:,i,:,:)
;vavol1(:,:,j*12+i,:,:)  = vavol(:,:,i+j*12,:,:)-vavolc(:,:,i,:,:)
end do 
end do 
printVarSummary(xavol1)

;------------------EL----------
;cem = new((/nm,nmon,ndim(1),ndim(2)/),float)
;EL = asciiread(diri+"", dimensions, datatype)
;do im=0,nm-1
;ne = num(EL(im,:).not.1e+20)
;ALLe1 = new((/ne,nmon,ndim(1),ndim(2)/), float, xa@_FillValue)
;ALLe1_clm = new((/ne,12,dim(1),ndim(2)/),float)
;ALLea1  = new((/ne,nmon,ndim(1),ndim(2)/),float)
;ymvol3 = new(ne, integer)
;do ie =0,ne-1
;    ymvol3(ie) = ind(yyyymms3.eq.EL(im,ie));EL(iv));
;    ALLe1(ie,:,:,:) = xa(im,ymvol2(ie)-ny1*12:ymvol2(ie)+ny2*12-1,:,:)
;    do j=0,11
;    ALLe1_clm(ie,j,:,:) = dim_avg_n_Wrap(ALLe1(ie,j::12,:,:), 1)
;    end do 
;
;    do i = 0,11
;    do j = 0,nyr-1
;    ALLea1(ie,j*12+i,:,:)  = ALLe1(ie,i+j*12,:,:)-ALLe1_clm(ie,i,:,:)
;    end do 
;    end do 
;
;end do 
;cem(:,im,:,:) = dim_avg_n_Wrap(ALLea1, 0)
;delete([/ALLe1,ALLe1_clm,ALLea1,ymvol3/])
;end do 
;printVarSummary(cem)
;
;ce_p = new((/nmon,ndim(2)/),float)
;ce_p = dim_avg_n_Wrap(cem1, (/1,2/))
;printVarSummary(ce_p)
;;--------------sort----------------------
a1 = readAsciiTable(diri+"cmip.nino.sort.cla.filter.txt", 5, "float", 0)
a2 = readAsciiTable(diri+"cmip.nino.sort.cel.filter.txt", 5, "float", 0)
a3 = readAsciiTable(diri+"cmip.nino.sort.wel.filter.txt", 5, "float", 0)
a1!0 = "model"
a1!1 = "vol"
a2!0 = "model"
a2!1 = "vol"
a3!0 = "model"
a3!1 = "vol"
printVarSummary(a2)

A = new((/3,nm,5/), "float", a1@_FillValue)
A(0,:,:) = a1
A(1,:,:) = a2
A(2,:,:) = a3
AC = (/"A","B","C"/)
do ia = 0,2
xxx = A(ia,:,:)
;uuu = ua1
;vvv = va1

xaa = xavol1
;xa2 = xavol1
;xa3 = xavol1
;ua1 = uavol1
;ua2 = uavol1
;ua3 = uavol1
;va1 = vavol1
;va2 = vavol1
;va3 = vavol1
do imon = 0,nmon-1
  do ilat = 0,ndim(1)-1
    do ilon = 0,ndim(2)-1
      xaa(:,:,imon,ilat,ilon) = xxx
   ;   ua1(:,:,imon,ilat,ilon) = cla
   ;   ua2(:,:,imon,ilat,ilon) = cel      
   ;   ua3(:,:,imon,ilat,ilon) = wel      
   ;   va1(:,:,imon,ilat,ilon) = cla
   ;   va2(:,:,imon,ilat,ilon) = cel
   ;   va3(:,:,imon,ilat,ilon) = wel
    end do 
end do 
end do 
      xaa = where(xaa.ne.0.0, xavol1, xa@_FillValue)
    ;  ua1 = where(ua1.ne.0.0, uavol1, ua@_FillValue)
    ;  ua2 = where(ua2.ne.0.0, uavol1, ua@_FillValue)
    ;  ua3 = where(ua3.ne.0.0, uavol1, ua@_FillValue)
    ;  va1 = where(va1.ne.0.0, vavol1, va@_FillValue)
    ;  va2 = where(va2.ne.0.0, vavol1, va@_FillValue)
    ;  va3 = where(va3.ne.0.0, vavol1, va@_FillValue)
printVarSummary(xaa)



;---------------composite---------------------
;xxx = xa1
;uuu = ua1
;vvv = va1
xvolavg  = dim_avg_n_Wrap(xaa, 1)     ;一个模式所有火山的合成
xmodavg = dim_avg_n_Wrap(xaa, 0)      ;一个火山所有模式的合成
;uvolavg  = dim_avg_n_Wrap(uaa, 1)     
;umodavg = dim_avg_n_Wrap(uaa, 0)      
;vvolavg  = dim_avg_n_Wrap(vaa, 1)     
;vmodavg = dim_avg_n_Wrap(vaa, 0)      
;svolavg1  = dim_avg_n_Wrap(sseav_p, 0)     
;uvolavg1  = dim_avg_n_Wrap(useavol, 0)     
;vvolavg1  = dim_avg_n_Wrap(vseavol, 0) 
printVarSummary(xvolavg)
printVarSummary(xmodavg)

;-------------t test ------------注意!!!样本是异常场!!总体也得是!!!!
;avesumxv = new((/nm,nmon,ndim(1),ndim(2)/), float, xa@_FillValue)
;varsumxv = new((/nm,nmon,ndim(1),ndim(2)/), float, xa@_FillValue)
;avesumxm = new((/nv,nmon,ndim(1),ndim(2)/), float, xa@_FillValue)
;varsumxm = new((/nv,nmon,ndim(1),ndim(2)/), float, xa@_FillValue)
;        
;avesumxm = 0 ;既然是异常场的平均值,那就应该是零啊!!
;avesumxv = 0
;
;varxm = dim_variance_n_Wrap(xxx, 0)
;varxv = dim_variance_n_Wrap(xxx, 1)
;;varu = dim_variance_n_Wrap(uavol, 0)
;;varv = dim_variance_n_Wrap(vavol, 0)
;siglvl=0.05
;do iv=0,nv-1
;do ik=0,nm-1
;do j =0,(ny1+ny2)-1
;do i =0,11
;varsumxv(ik,i+j*12,:,:) = dim_variance_n_Wrap(xa(ik,i::12,:,:), 0)
;varsumxm(iv,i+j*12,:,:) = dim_variance_n_Wrap(xa(ik,i::12,:,:), 0)
;
;;varsumu(i+j*12,:,:) = dim_variance_n_Wrap(ua(i::12,:,:), 0)
;;varsumv(i+j*12,:,:) = dim_variance_n_Wrap(va(i::12,:,:), 0)
;
;end do 
;end do
;end do 
;end do 
;printVarSummary(varsumxv)
;printVarSummary(varsumxm)
;
;                                                    ;所有sst不同月的方差
;
;;sxep = ssumsf
;
;ssum = 2004-1870
;snm = nm
;snv = nv
;probv = xvolavg
;probm = xmodavg
;
;;probu = uvolavg
;;probv = vvolavg
;
;    iflag= True              ; population variance similar 
;
;  probv=(1- ttest(xvolavg,varxv,snv, avesumxv,varsumxv,ssum, iflag, False))*100
;  probm=(1- ttest(xmodavg,varxm,snm, avesumxm,varsumxm,ssum, iflag, False))*100
; ; probu=(1- ttest(uvolavg,varu,snp, avesumu,varsumu,ssum, iflag, False))*100
; ; probv=(1- ttest(vvolavg,varv,snp, avesumv,varsumv,ssum, iflag, False))*100
;
;  copy_VarMeta(xvolavg,probv)
;  copy_VarMeta(xmodavg,probm)
;
;; copy_VarMeta(uvolavg,probu)
;; copy_VarMeta(vvolavg,probv)  
;
;printVarSummary(probm)
;printVarSummary(probv)
;
;-----------------season-------------
xmsea = new((/nv,4,nyr,ndim(1),ndim(2)/), float, xa@_FillValue)
xvsea = new((/nm,4,nyr,ndim(1),ndim(2)/), float, xa@_FillValue)
;umsea = new((/nv,4,nyr,ndim(1),ndim(2)/), float, xa@_FillValue)
;uvsea = new((/nm,4,nyr,ndim(1),ndim(2)/), float, xa@_FillValue)
;vmsea = new((/nv,4,nyr,ndim(1),ndim(2)/), float, xa@_FillValue)
;vvsea = new((/nm,4,nyr,ndim(1),ndim(2)/), float, xa@_FillValue)
;probm_sea = new((/nv,4,nyr,ndim(1),ndim(2)/), float, xa@_FillValue)
;probv_sea = new((/nm,4,nyr,ndim(1),ndim(2)/), float, xa@_FillValue)

do ik = 0,nm-1
do iv = 0,n-1
xmsea(iv,:,:,:,:) = month_to_seasonN(xmodavg(iv,:,:,:), SEASON)
xvsea(ik,:,:,:,:) = month_to_seasonN(xvolavg(ik,:,:,:), SEASON)
;umsea(iv,:,:,:,:) = month_to_seasonN(umodavg(iv,:,:,:), SEASON)
;uvsea(ik,:,:,:,:) = month_to_seasonN(uvolavg(ik,:,:,:), SEASON)
;vmsea(iv,:,:,:,:) = month_to_seasonN(vmodavg(iv,:,:,:), SEASON)
;vvsea(ik,:,:,:,:) = month_to_seasonN(vvolavg(ik,:,:,:), SEASON)
;probv_sea(ik,:,:,:,:) = month_to_seasonN(probv(ik,:,:,:), SEASON)
;probm_sea(iv,:,:,:,:) = month_to_seasonN(probm(iv,:,:,:), SEASON)
end do 
end do

printVarSummary(xvsea)
printVarSummary(xmsea)
;printVarSummary(probv_sea)
;printVarSummary(probm_sea)

xsea = dim_avg_n_Wrap(xvsea, 0)
;usea = dim_avg_n_Wrap(uvsea, 0)
;vsea = dim_avg_n_Wrap(vvsea, 0)

;psea = dim_avg_n_Wrap(probv_sea, 0)
printVarSummary(xsea)
;print(psea)
;--------------plot-------------------
;wks = gsn_open_wks("ps", diri+"obs1.composite.90%.TP.volcano.u10.lon-time."+ny1+ny2);(eruptmon)
;do ik=0,nm-1
;do iv=0,nv-1
;wks = gsn_open_wks("ps", diri+"cmip.composite.TP.volcano.tas.volavg"+name_model(ik));."+volyear(iv)
;wks = gsn_open_wks("ps", diri+"cmip.composite.TP.volcano.sst.modavg"+volyear(iv));"+name_model(ik));.
wks = gsn_open_wks("ps", diri+"cmip.composite.volcano.1000u.sort."+AC(ia));"+name_model(ik));.

gsn_define_colormap(wks,"MPL_PRGn");"BlueDarkRed18""BrownBlue12""MPL_BrBG" "BlueWhiteOrangeRed" 

;gsn_define_colormap(wks, "BlueWhiteOrangeRed")
do iyr = 4,8
plot = new(4, graphic)
res2= True
res2@gsnDraw = False
res2@gsnFrame = False
res2@mpCenterLonF = 180.
res2@gsnLeftString = " "
sdres = res2
vcres = res2  


res2@cnLinesOn     = False
res2@gsnAddCyclic         =  False ;True;        ; data not cyclic
  res2@mpMaxLatF                   = NLAT         ; choose subregion           
  res2@mpMinLatF                   = SLAT
  res2@mpMaxLonF                   = ELON
  res2@mpMinLonF                   = WLON
res2@tmXMajorGrid                = True
res2@tmYMajorGrid                = True
res2@tmXMajorGridThicknessF      = 1.0
res2@tmYMajorGridThicknessF      = 1.0
res2@tmXMajorGridLineDashPattern = 2
res2@tmYMajorGridLineDashPattern = 2

  opt = res2

res2@cnFillOn             = True         ; turn on color fill
sdres@cnFillOn = False
sdres@cnLinesOn = False
sdres@cnLineLabelsOn        = False; True   
res2@cnInfoLabelOn            = False
  res2@cnLineLabelsOn           = False ;True
 
        sdres@cnLevelSelectionMode = "ManualLevels"  
        sdres@cnMinLevelValF        = 90
      sdres@cnLevelSpacingF      =  5         ;间隔5个画
        sdres@cnLevels = (/80,90,95,99.9/)
sdres@cnSmoothingOn = True
sdres@cnSmoothingTensionF = 1
sdres@cnInfoLabelOn = False
sdres@gsnLeftString = " "

res2@tiMainString = " "

res2@mpProjection         = "Robinson" ;"Mollweide"
  res2@mpPerimOn            = False     ;Turns off the box perimeter.
  res2@mpGridLatSpacingF    =  90                ; change latitude  line spacing
  res2@mpGridLonSpacingF    = 360.               ; change longitude line spacing
 ;res2@mpGridLineColor      = "transparent"      ; trick ncl into drawing perimeter
  res2@mpGridAndLimbOn      = True               ; turn on lat/lon lines
  res2@mpFillOn          = False;True ;

  res2@lbLabelBarOn        = False;True ;            ; turn off individual cb's
res2@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
res2@cnMinLevelValF        = -2.0;280-273.16
res2@cnMaxLevelValF        =  2.0 ; 310-273.16
res2@cnLevelSpacingF       = 0.1
;res2@cnLevels = fspan(-2.0, 2.0, 20)
; res2@cnFillDrawOrder      = "Predraw"  ; draw contours first 可以用地图把陆地盖住
;plot(i+12*j) = gsn_csm_contour_map_ce(wks2,sepavg(j*12+i,:,:),res2)
;plot(2*i+j) = gsn_csm_contour_map_ce(wks2,probep(j*12+i,:,:),res2)



  opt@gsnShadeFillType = "pattern"      ; pattern fill
  opt@gsnShadeHigh = 17          ; 高于95的use pattern #17
  opt@gsnshadelow = -1       ;低于95为空白  
  opt@cnLineLabelsOn        = False   
  opt@gsnLeftString  = " "
  opt@cnInfoLabelOn  = False
;  ptest0 = gsn_csm_contour(wks2, prob0, sdres)
;ptest0 = gsn_contour_shade(ptest0,-999,95,opt) ;Pattern shade contours above 95.

ressf = True
ressf@gsnFrame = False
ressf@gsnDraw = False
ressf@cnLinesOn = True
ressf@cnFillOn  = False
ressf@cnLineLabelsOn = False
ressf@cnInfoLabelOn       = False
    ressf@gsnContourNegLineDashPattern = 1
 ;ressf@cnLineThicknessF = 1.5

vcres@vcMinDistanceF           = 0.017            ; thin vectors
;res2@cnFillPalette           = "gui_default"     ; set color map
vcres@vcRefAnnoOrthogonalPosF = -.535    ; move ref vector into plot
  vcres@vcRefAnnoFontHeightF    = 0.005
vcres@vcRefMagnitudeF          = 3.0             ; define vector ref mag
vcres@vcRefLengthF             = 0.1            ; define length of vec ref
vcres@vcRefAnnoOrthogonalPosF  = .1               ; move ref vector down
;vcres@vcRefLengthF             = 0.045            ; define length of vec ref
vcres@gsnAddCyclic         = False        ; data not cyclic 
vcres@gsnLeftString = ""
vcres@gsnRightString = ""
vcres@lbLabelBarOn = False
vcres@vcGlyphStyle = "CurlyVector"

res2@gsnLeftStringFontHeightF = 0.02
res2@gsnRightStringFontHeightF = 0.02
;res2@tmXBLabelFontHeightF = 0.014
;res2@tmYBLabelFontHeightF = 0.014
resk = True
resk@txFontHeightF = 0.02
resk@gsnPanelLabelBar   = True                 ; add common colorbar
resk@gsnMaximize = True
resk@gsnPaperOrientation = "portrait"   ; force portrait
;resk@lbBoxEndCapStyle = "TriangleBothEnds"
do isea = 1,3

res2@gsnLeftString = SEASON(isea);name_model(ik) ;"A";MO(j); LABEL(isea)+
iyr1 = iyr-5
res2@gsnRightString = "year("+iyr1+")"
;resk@txString = "composite surface pr "+name_model(ik)+" year("+iyr1+")";
;plot(isea)  = gsn_csm_contour_map_ce(wks, xvsea(ik,isea,iyr,:,:), res2) 
;ptest = gsn_csm_hov(wks, probv_sea(ik,isea,iyr,:,:), sdres)
resk@txString = "composite 1000 ua year("+iyr1+") in "+AC(ia);
kk = isea-1
plot(kk)  = gsn_csm_contour_map(wks, xsea(isea,iyr,:,:), res2) 
;vcplot = gsn_csm_vector(wks, usea(isea,iyr,:,:),vsea(isea,iyr,:,:), vcres)

;ptest = gsn_csm_hov(wks, psea(isea,iyr,:,:), sdres)
;ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
;overlay(plot(isea), ptest)
;overlay(plot(kk), vcplot)
end do 
res2@gsnLeftString = "DJF";name_model(ik) ;"A";MO(j); LABEL(isea)+
plot(3)  = gsn_csm_contour_map(wks, xsea(0,iyr+1,:,:), res2) 
;vcplot = gsn_csm_vector(wks, usea(0,iyr+1,:,:),vsea(0,iyr+1,:,:), vcres)
;overlay(plot(3), vcplot)



gsn_panel(wks, plot, (/4,1/), resk)
delete(res2)
end do 
;end do 
;delete([/res2,x,xx,xc,xa,TIME,yyyymm,ymstrt1,ymlast1,files,f,xavol,xseavol,xavolc,xavol1,xseavol1/])
;print(ik)
end do 
end