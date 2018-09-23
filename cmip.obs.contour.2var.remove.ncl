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
diri = "/data2/home/liufeiclim2014/xing/cmip_vol/";"/nuist/p/work/feil/xc/";
diri2 = "/data2/data/cmip5/historical_regrid/atmos/mon/"
diri3 = "/data2/data/cmip5/historical/atmos/mon/"
diri0 = "/data2/home/liufeiclim2014/xing/";"/nuist/p/work/feil/xc/";
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

ny1 = 5 ; year before year0
ny2 = 5 ; year after year0
nyr = ny1+ny2
nmon = ny1*12+ny2*12

;f3= addfile(diri0+"HadISST_1870-2016.nc", "r");sst.mnmean.v4.1854-2016.nc", "r");            
f4= addfile(diri+"prate.mon.mean.vc2.1851-2014.nc", "r")
;f5 = addfile(diri0+"air.sfc.mon.mean.vc2.1851-2014.nc", "r")
;f6 = addfile(diri+"uswrf.sfc.mon.mean.vc2.1851-2014.nc", "r")
;f7 = addfile(diri+"uswrf.ntat.mon.mean.vc2.1851-2014.nc", "r")
;f8 = addfile(diri+"dswrf.sfc.mon.mean.vc2.1851-2014.nc", "r")
;f9 = addfile(diri+"i20.anom.1871-2008.nc", "r")
;f10 = addfile(diri2+"ssh.soda.1871-2008.TP.nc", "r")
;f11 = addfile("/data/slp(prmsl).mon.mean.1851-2014.nc" , "r" )
;f12 = addfile(diri0+"hgt.mon.mean.vc2.1851-2014.nc" , "r" )
;f13 = addfile(diri+"shum.2m.mon.mean.v2.1851-2014.nc" , "r" )
f14 = addfile(diri0+"air.mon.mean.vc2.1851-2014.nc", "r")
MO = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
SEASON = (/"DJF","MAM", "JJA","SON" /)
;---------sst--------------
TIME1=f4->time
yyyymm=cd_calendar(TIME1, -1)
printMinMax(yyyymm, True)

ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)


sss = f4->prate(ymstrt1:ymlast1,:,:);LEV,
 sss2 = sss
;sss2 = sss*1000 ; shum kg/kg->g/kg
;sss2 = sss/100. ;slp Pa-->hPa
sss2 = sss*24*3600 ;prate-->precip
;sss1 = lonPivot(sss,0.5) ;HadISST center change
; printVarSummary(sss1)
;-------去掉神奇的-1000（不是missing value）-----
;dimo = dimsizes(sss1)
;sstone = ndtooned(sss1)
;Y = ind(sstone.lt.-10);.eq.-1000)
;sstone(Y) = sss1@_FillValue
;sss2 = reshape(sstone, dimo)
;copy_VarMeta(sss1, sss2)
;printVarSummary(sss2)
;sss2!0 = "time"
;sss2!1 = "lat"
;sss2!2 = "lon"
;printVarSummary(sss2)

sst = sss2(:,{SLAT:NLAT},{WLON:ELON})
sst!0 = "time"
sst!1 = "lat"
sst!2 = "lon"
sst= dtrend_msg_n(sst&time,sst,False,False,0)
;  copy_VarMeta(sst,sstd)
LAT=sst&lat
LON=sst&lon
sndim = dimsizes(sst)
;sc1 = clmMonTLL(sst1)
;sc2 = clmMonTLL(sst2)
sc  = clmMonTLL(sst)
sa  = calcMonAnomTLL(sst, sc)
;ssea = month_to_seasonN(sa, SEASON)
printVarSummary(sa)

;---------air-----------------------------
TIME2=f14->time
yyyymm2=cd_calendar(TIME2, -1)

ymstrt2 =ind(yyyymm2.eq.ymstrt)
ymlast2 = ind(yyyymm2.eq.ymlast)

LEV = 20 ;500hPa:10 200hPa:16 850hPa:3 50hPa:20 300hPa:14
u = f14->air(ymstrt2:ymlast2,LEV,:,:);16,{SLAT:NLAT},{WLON:ELON}
;v = f2->vwnd(ymstrt2:ymlast2,LEV,{SLAT:NLAT},{WLON:ELON})
u= dtrend_msg_n(u&time,u,False,False,0)
;v= dtrend_msg_n(v&time,v,False,False,0)
uLAT=u&lat
uLON=u&lon
undim = dimsizes(u)
uc  = clmMonTLL(u)
;vc  = clmMonTLL(v)
ua  = calcMonAnomTLL(u, uc)
;va  = calcMonAnomTLL(v, vc)
usea = month_to_seasonN(ua, SEASON)
;vsea = month_to_seasonN(va, SEASON)
printVarSummary(u)



;-------------reg--------------------先remove还是先滤波？

nynino = 2004-1870+1
obs_nino34 =  asciiread(diri+"obs.ersst.nino34.anomaly.unnorm.1870-2016.txt", (/nynino*12/), "float")
reg_sa = regCoef(obs_nino34, sa(lat|:,lon|:,time|:))
reg_ua = regCoef(obs_nino34, ua(lat|:,lon|:,time|:))
printVarSummary(reg_sa)
rem_sa = sa
rem_ua = ua
do it = 0,nynino*12-1
rem_sa(it,:,:) = sa(it,:,:)-reg_sa*obs_nino34(it)
rem_ua(it,:,:) = ua(it,:,:)-reg_ua*obs_nino34(it)
end do 
copy_VarMeta(sa, rem_sa)
copy_VarMeta(ua, rem_ua)

printVarSummary(rem_sa)

;----------------filter----------------------PAY ATTENTION! This is month data, not years!
xf = rem_sa
uf = rem_ua
;vf = va(ik,:,:,:)
ihp   = 1       ; low-pass filter: ihp = 0; high-pass ihp = 1; band-pass ihp = 2.
sigma = 1.0                           ; Lanczos sigma
nWgt  = 71 ;必须是奇数 表示权重系数的总个数，越大效果越好，但是两端损失越多
fca   = 1./(10*12) ;                        
fcb   = -999;1./(2*12);1./8. *12  ;8——50年滤波
wgt  = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
xf = wgt_runave_n (rem_sa, wgt, 0, 0 ) ;最后一个是指sst1第0维
uf = wgt_runave_n (rem_ua, wgt, 0, 0 ) ;最后一个是指sst1第0维

n = nv

ny1 = 5 ; year before year0
ny2 = 5 ; year after year0
nyr = ny1+ny2
nmon = ny1*12+ny2*12
nyear = (/"year(-5)","year(-4)","year(-3)","year(-2)","year(-1)","year(0)","year(+1)",\
         "year(+2)","year(+3)","year(+4)","year(+5)","year(+6)","year(+7)"/);"year(-6)",

savol =   new((/n,  nmon,sndim(1),sndim(2)/), typeof(sa), sa@_FillValue) 
sseavol = new((/n, 4,nyr,sndim(1),sndim(2)/), typeof(sa), sa@_FillValue) 
uavol =   new((/n,  nmon,undim(1),undim(2)/), typeof(ua), ua@_FillValue) 
;vavol =   new((/n,  12*(ny1+ny2),undim(1),undim(2)/), typeof(va), va@_FillValue) 
useavol = new((/n, 4,nyr,undim(1),undim(2)/), typeof(ua), ua@_FillValue) 
;vseavol = new((/n, 4,ny1+ny2,undim(1),undim(2)/), typeof(va), va@_FillValue) 

;TIMEu2=ua&time
;yyyymmu2=cd_calendar(TIMEu2, -1)
;;yyyyu2=cd_calendar(TIMEu2, 4)
;do iv = 0,n-1
;ymvol = ind(yyyymmu2.eq.volmon(iv));EL(iv));
;uavol(iv,:,:,:) = ua(ymvol:ymvol+ny*12-1,:,:)
;vavol(iv,:,:,:) = va(ymvol:ymvol+ny*12-1,:,:)
;useavol(iv,:,:,:,:) = month_to_seasonN(uavol(iv,:,:,:), SEASON)
;vseavol(iv,:,:,:,:) = month_to_seasonN(vavol(iv,:,:,:), SEASON)
;
;end do 

xxx = xf

TIMEs2=sa&time
yyyymms2=cd_calendar(TIMEs2, -1)
;yyyys2=cd_calendar(TIMEs2, 4)
do iv = 0,n-1
ymvol2 = ind(yyyymms2.eq.volmon(iv));EL(iv));
savol(iv,:,:,:) = xxx(ymvol2-ny1*12:ymvol2+ny2*12-1,:,:)
sseavol(iv,:,:,:,:) = month_to_seasonN(savol(iv,:,:,:), SEASON)
uavol(iv,:,:,:) = uf(ymvol2-ny1*12:ymvol2+ny2*12-1,:,:)
useavol(iv,:,:,:,:) = month_to_seasonN(uavol(iv,:,:,:), SEASON)
end do 
savol!0 = "sample"
sseavol!0 = "sample"
uavol!0 = "sample"
useavol!0 = "sample"
;printVarSummary(u)
;printVarSummary(va)
;printVarSummary(uavol)
;-------------remove nearby annual---------------------
savolc = new((/n,12,sndim(1),sndim(2)/), typeof(sa), sa@_FillValue)
savol1 = new((/n,nmon,sndim(1),sndim(2)/), typeof(sa), sa@_FillValue)
sseavol1 = new((/n, 4,nyr,sndim(1),sndim(2)/), typeof(sa), sa@_FillValue) 
uavolc = new((/n,12,undim(1),undim(2)/), typeof(ua), ua@_FillValue)
uavol1 = new((/n,nmon,undim(1),undim(2)/), typeof(ua), ua@_FillValue)
useavol1 = new((/n, 4,nyr,undim(1),undim(2)/), typeof(ua), ua@_FillValue) 

do i=0,11
  savolc(:,i,:,:) = dim_avg_n_Wrap(savol(:,i::12,:,:), 1)
  uavolc(:,i,:,:) = dim_avg_n_Wrap(uavol(:,i::12,:,:), 1)
end do 
savolc!0 = "sample"
savol1 = savol
uavolc!0 = "sample"
uavol1 = uavol
do i = 0,11
do j = 0,nyr-1
savol1(:,j*12+i,:,:)  = savol(:,i+j*12,:,:)-savolc(:,i,:,:)
uavol1(:,j*12+i,:,:)  = uavol(:,i+j*12,:,:)-uavolc(:,i,:,:)
end do 
end do 
printVarSummary(savol1)
do iv = 0,n-1
sseavol1(iv,:,:,:,:) = month_to_seasonN(savol1(iv,:,:,:), SEASON)
useavol1(iv,:,:,:,:) = month_to_seasonN(uavol1(iv,:,:,:), SEASON)
end do 

;---------------composite---------------------
svolavg  = dim_avg_n_Wrap(savol1, 0)     
uvolavg  = dim_avg_n_Wrap(uavol1, 0)     
;uvolavg  = dim_avg_n_Wrap(uavol, 0)     
;vvolavg  = dim_avg_n_Wrap(vavol, 0) 
svolavg1  = dim_avg_n_Wrap(sseavol1, 0)     
uvolavg1  = dim_avg_n_Wrap(useavol, 0)     
;vvolavg1  = dim_avg_n_Wrap(vseavol, 0) 
printVarSummary(svolavg)
;-------------t test ------------注意!!!样本是异常场!!总体也得是!!!!
avesumsst = new((/12*nyr,sndim(1),sndim(2)/), float, sst@_FillValue)
varsumsst = new((/12*nyr,sndim(1),sndim(2)/), float, sst@_FillValue)
avesumu =   new((/12*nyr,undim(1),undim(2)/), float, ua@_FillValue)
varsumu =   new((/12*nyr,undim(1),undim(2)/), float, ua@_FillValue)
;avesumv =   new((/12*ny,undim(1),undim(2)/), float, va@_FillValue)
;varsumv =   new((/12*ny,undim(1),undim(2)/), float, va@_FillValue)

        
avesumsst = 0 ;既然是异常场的平均值,那就应该是零啊!!
avesumu = 0
avesumv = 0

varsst = dim_variance_n_Wrap(savol1, 0)
varu = dim_variance_n_Wrap(uavol1, 0)
;varv = dim_variance_n_Wrap(vavol, 0)

siglvl=0.05

do j =0,(ny1+ny2)-1
do i =0,11
varsumsst(i+j*12,:,:) = dim_variance_n_Wrap(sa(i::12,:,:), 0)
varsumu(i+j*12,:,:) = dim_variance_n_Wrap(ua(i::12,:,:), 0)
;varsumv(i+j*12,:,:) = dim_variance_n_Wrap(va(i::12,:,:), 0)

end do 
end do
printVarSummary(varsumsst)   
                                                    ;所有sst不同月的方差

;sxep = ssumsf

ssum = 2016-1870
snp = n
probsst = svolavg
probu = uvolavg
;probv = vvolavg

    iflag= True              ; population variance similar 

  probsst=(1- ttest(svolavg,varsst,snp, avesumsst,varsumsst,ssum, iflag, False))*100
  probu=(1- ttest(uvolavg,varu,snp, avesumu,varsumu,ssum, iflag, False))*100
 ; probv=(1- ttest(vvolavg,varv,snp, avesumv,varsumv,ssum, iflag, False))*100

  copy_VarMeta(svolavg,probsst)
 copy_VarMeta(uvolavg,probu)
; copy_VarMeta(vvolavg,probv)  

printVarSummary(probsst)
probsst1 = month_to_seasonN(probsst, SEASON)
probu1 = month_to_seasonN(probu, SEASON)
;ssst = svolavg1
;do i1 = 0,
;ssst = where(.not.all(ismissing(probsst1).and.probsst1.lt.80), svolavg1, svolavg@_FillValue)
;printVarSummary(ssst)

;============CMIP5=======================================
;--------------------------
f1 = addfile(diri+"cmip.vol.pr.remove.nc", "r")
xavol = f1->xavol(:,:,:,{SLAT:NLAT},{WLON:ELON})
;xavol = xavol*1000 ; hus kg/kg->g/kg
printVarSummary(xavol)
ndim = dimsizes(xavol)

f2 = addfile(diri+"cmip.vol.50ta.remove.nc", "r")
xavol2 = f2->xavol(:,:,:,{SLAT:NLAT},{WLON:ELON})
;xavol = xavol*1000 ; hus kg/kg->g/kg
printVarSummary(xavol2)
ndim2 = dimsizes(xavol2)
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
xaa = new((/3,nm,nv,nmon,ndim(3),ndim(4)/),"float", xavol@_FillValue)
xaa2 = new((/3,nm,nv,nmon,ndim(3),ndim(4)/),"float", xavol@_FillValue)
delete(xxx)
do ia = 0,2
xxx = A(ia,:,:)
;uuu = ua1
;vvv = va1


;xa3 = xavol
;ua1 = uavol
;ua2 = uavol
;ua3 = uavol
;va1 = vavol
;va2 = vavol
;va3 = vavol
do imon = 0,nmon-1
  do ilat = 0,ndim(3)-1
    do ilon = 0,ndim(4)-1
      xaa(ia,:,:,imon,ilat,ilon) = xxx
      xaa2(ia,:,:,imon,ilat,ilon) = xxx
   ;   ua1(:,:,imon,ilat,ilon) = cla
   ;   ua2(:,:,imon,ilat,ilon) = cel      
   ;   ua3(:,:,imon,ilat,ilon) = wel      
   ;   va1(:,:,imon,ilat,ilon) = cla
   ;   va2(:,:,imon,ilat,ilon) = cel
   ;   va3(:,:,imon,ilat,ilon) = wel
    end do 
end do 
end do 
      xaa(ia,:,:,:,:,:) = where(xaa(ia,:,:,:,:,:).ne.0.0, xavol, xavol@_FillValue)
      xaa2(ia,:,:,:,:,:) = where(xaa2(ia,:,:,:,:,:).ne.0.0, xavol2, xavol2@_FillValue)
    ;  ua1 = where(ua1.ne.0.0, uavol, ua@_FillValue)
    ;  ua2 = where(ua2.ne.0.0, uavol, ua@_FillValue)
    ;  ua3 = where(ua3.ne.0.0, uavol, ua@_FillValue)
    ;  va1 = where(va1.ne.0.0, vavol, va@_FillValue)
    ;  va2 = where(va2.ne.0.0, vavol, va@_FillValue)
    ;  va3 = where(va3.ne.0.0, vavol, va@_FillValue)
copy_VarMeta(xavol, xaa2(ia,:,:,:,:,:))
copy_VarMeta(xavol, xaa(ia,:,:,:,:,:))
end do 
xaa!0 = "sort"
xaa2!0 = "sort"
printVarSummary(xaa)



;-------------t test ------------注意!!!样本是异常场!!总体也得是!!!!
xva_1 = reshape(xavol, (/nm*nv,nmon,ndim(3),ndim(4)/))
xaa_1 = reshape(xaa, (/3,nm*nv,nmon,ndim(3),ndim(4)/))
xva_2 = reshape(xavol2, (/nm*nv,nmon,ndim2(3),ndim2(4)/))
xaa_2 = reshape(xaa2, (/3,nm*nv,nmon,ndim2(3),ndim2(4)/))

avesum = new((/3,nmon,ndim(3),ndim(4)/), float, xavol@_FillValue)
varsum = new((/3,nmon,ndim(3),ndim(4)/), float, xavol@_FillValue)
        
avesum = 0 ;既然是异常场的平均值,那就应该是零啊!!
avesum2 = avesum
varsum2 = varsum
prob = xaa_1(:,0,:,:,:)
prob2 = xaa_1(:,0,:,:,:)
xaa_11 = dim_avg_n_Wrap(xaa_1, 1)
xaa_22 = dim_avg_n_Wrap(xaa_2, 1)

do ia = 0,2
varx = dim_variance_n_Wrap(xaa_1(ia,:,:,:,:), 0)
varx2 = dim_variance_n_Wrap(xaa_2(ia,:,:,:,:), 0)

siglvl=0.05
;do ik=0,nm*nv-1
;do j =0,(ny1+ny2)-1
;do i =0,11
varsum(ia,:,:,:) = dim_variance_n_Wrap(xva_1, 0) ;?!(i+j*12,:,:)(:,i::12,:,:)
varsum2(ia,:,:,:) = dim_variance_n_Wrap(xva_2, 0) ;?!(i+j*12,:,:)(:,i::12,:,:)

;varsumu(i+j*12,:,:) = dim_variance_n_Wrap(ua(i::12,:,:), 0)
;varsumv(i+j*12,:,:) = dim_variance_n_Wrap(va(i::12,:,:), 0)

;end do
;end do 
;end do 
                                                    ;所有sst不同月的方差

;sxep = ssumsf

;ssum = 2004-1870
ssum = nm*nv
sn   = num(A(ia,:,:).ne.0.0)
print(sn)

    iflag= True              ; population variance similar 

  prob(ia,:,:,:)=(1- ttest(xaa_11(ia,:,:,:),varx,sn, avesum(ia,:,:,:),varsum(ia,:,:,:),ssum, iflag, False))*100

     iflag= True              ; population variance similar 

  prob2(ia,:,:,:)=(1- ttest(xaa_22(ia,:,:,:),varx2,sn, avesum2(ia,:,:,:),varsum2(ia,:,:,:),ssum, iflag, False))*100
delete(sn)
end do 
  copy_VarMeta(xaa(:,0,0,:,:,:),prob2)
  copy_VarMeta(xaa(:,0,0,:,:,:),prob)

printVarSummary(prob)
 

;asciiwrite(diri+"cmip.prob.air"+AC(ia)+".txt", prob)
;---------------composite & season---------------------
xavg  = dim_avg_n_Wrap(xaa, (/1,2/))     ;合成
xavg2  = dim_avg_n_Wrap(xaa2, (/1,2/))     ;合成
xsea = new((/3,4,nmon/12,ndim(3),ndim(4)/), float, xavol@_FillValue)
psea = new((/3,4,nmon/12,ndim(3),ndim(4)/), float, xavol@_FillValue)
xsea2 = new((/3,4,nmon/12,ndim(3),ndim(4)/), float, xavol@_FillValue)
psea2 = new((/3,4,nmon/12,ndim(3),ndim(4)/), float, xavol@_FillValue)
do ia = 0,2
xsea(ia,:,:,:,:) = month_to_seasonN(xavg(ia,:,:,:), SEASON);(time|:,lev|:,LON|:)
xsea2(ia,:,:,:,:) = month_to_seasonN(xavg2(ia,:,:,:), SEASON);(time|:,lev|:,LON|:)
psea(ia,:,:,:,:) = month_to_seasonN(prob(ia,:,:,:), SEASON)
psea2(ia,:,:,:,:) = month_to_seasonN(prob2(ia,:,:,:), SEASON)
end do 
printVarSummary(psea)
printVarSummary(xsea2)

;xxxt = xsea
;xxxt = where(.not.all(ismissing(psea).and.psea.lt.80), xsea, xsea@_FillValue)
;;xxxt = xsea
;;xxxt = where((.not.all(ismissing(psea).and.psea.lt.80), xsea, xsea@_FillValue)
;printVarSummary(xxxt)

;p1 = ndtooned(psea)
;miss = ind(p1.gt.110.)
;p1(miss) = psea@_FillValue
;delete(miss)
;psea1 = reshape(p1, (/3,4,nmon/12,ndim(3),ndim(4)/))
;copy_VarMeta(psea, psea1)
;printVarSummary(psea1)
;
;
;p2 = ndtooned(psea2)
;miss2 = ind(p2.gt.110.)
;p2(miss2) = psea2@_FillValue
;delete(miss2)
;psea2 = reshape(p2, (/3,4,nmon/12,ndim(3),ndim(4)/))
;copy_VarMeta(psea2, psea2)
;printVarSummary(psea2)
;asciiwrite(diri+"cmip.prob.sea.air"+AC(ia)+".txt", psea1)
;--------------plot-------------------
;wks = gsn_open_wks("ps", diri+"obs1.composite.90%.TP.volcano.u10.lon-time."+ny1+ny2);(eruptmon)
;do ik=0,nm-1
;do iv=0,nv-1
;wks = gsn_open_wks("ps", diri+"cmip.composite.TP.volcano.tas.volavg"+name_model(ik));."+volyear(iv)
;wks = gsn_open_wks("ps", diri+"cmip.composite.TP.volcano.sst.modavg"+volyear(iv));"+name_model(ik));.
wks = gsn_open_wks("ps", diri+"cmip&obs.composite.volcano.50ta.sort.removeENSO");"+name_model(ik));.

gsn_define_colormap(wks,"BlueDarkRed18");"BlueDarkRed18""precip_diff_12lev""MPL_BrBG" "BlueWhiteOrangeRed" 

;gsn_define_colormap(wks, "BlueWhiteOrangeRed")
plot = new(4, graphic)
do isea = 0,3
  if (isea.eq.0) then
    iyr = 6
  else 
    iyr = 5
  end if
res2= True
res2@gsnDraw = False
res2@gsnFrame = False
res2@gsnLeftString = ""
res2@gsnRightString = ""
;res2@mpProjection         = "Robinson" ;"Mollweide"
sdres = res2
vcres = res2  
  opt = res2
res3  = res2

res2@cnLinesOn     = False
res2@gsnAddCyclic         =  True;False ;        ; data not cyclic
  res2@mpMaxLatF                   = NLAT         ; choose subregion           
  res2@mpMinLatF                   = SLAT;30 ;
  res2@mpMaxLonF                   = ELON
  res2@mpMinLonF                   = WLON
res2@mpCenterLonF = 180.;-90;
;res2@tmXMajorGrid                = True
;res2@tmYMajorGrid                = True
;res2@tmXMajorGridThicknessF      = 1.0
;res2@tmYMajorGridThicknessF      = 1.0
;res2@tmXMajorGridLineDashPattern = 2
;res2@tmYMajorGridLineDashPattern = 2


sdres@cnFillOn = False
sdres@cnLinesOn = False
sdres@cnLineLabelsOn        = False; True   
res2@cnFillOn             = True         ; turn on color fill
res2@cnInfoLabelOn            = False
  res2@cnLineLabelsOn           = False ;True
 
   ;     sdres@cnLevelSelectionMode = "ManualLevels"  
   ;     sdres@cnMinLevelValF        = 90
   ;   sdres@cnLevelSpacingF      =  5         ;间隔5个画
   ;     sdres@cnLevels = (/80,90,95,99.9/)
sdres@cnSmoothingOn = True
sdres@cnSmoothingTensionF = 1
sdres@cnInfoLabelOn = False
sdres@gsnLeftString = " "

res2@tiMainString = " "

 ; res2@mpPerimOn            = False     ;Turns off the box perimeter.
  res2@mpGridLatSpacingF    =  180.                ; change latitude  line spacing
  res2@mpGridLonSpacingF    = 360.               ; change longitude line spacing
 ;res2@mpGridLineColor      = "transparent"      ; trick ncl into drawing perimeter
 ; res2@mpGridAndLimbOn      = True               ; turn on lat/lon lines
  res2@mpFillOn          = False;True ;

  res2@lbLabelBarOn        = False;True ;            ; turn off individual cb's
res2@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
res2@cnMinLevelValF        = -2.0;280-273.16
res2@cnMaxLevelValF        =  2.0 ; 310-273.16
res2@cnLevelSpacingF       =  0.2
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

res3@cnLinesOn = True
res3@cnFillOn  = False
res3@cnLineLabelsOn = False;True ;
res3@cnLineLabelBackgroundColor = -1
res3@cnInfoLabelOn            = False
;res3@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
;res3@cnMinLevelValF        = -1.0;280-273.16
;res3@cnMaxLevelValF        =  1.0 ; 310-273.16
res3@cnLevelSpacingF       =  0.5
;res3@cnLineLabelInterval   =  1

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
;do isea = 1,3

;res2@gsnLeftString = SEASON(isea);name_model(ik) ;"A";MO(j); LABEL(isea)+
;iyr1 = iyr-5
;res2@gsnRightString = "year("+iyr1+")"
;resk@txString = "composite surface pr "+name_model(ik)+" year("+iyr1+")";
;plot(isea)  = gsn_csm_contour_map_ce(wks, xvsea(ik,isea,iyr,:,:), res2) 
;ptest = gsn_csm_hov(wks, probv_sea(ik,isea,iyr,:,:), sdres)
resk@txString = "composite 50hPa ta in year(0)"+SEASON(isea)+" of OBS and CMIP5 ENSO removed";hus(g/kg)

res2@gsnLeftString = "(a) observation"
res2@gsnRightString = ""
;res2@gsnPolar   = "NH"                          ; specify the hemisphere
;plot(0) = gsn_csm_contour_map_polar(wks, uvolavg1(0,6,:,:), res2)    ; {0:90}create the plot
plot(0)  = gsn_csm_contour_map(wks, uvolavg1(isea,iyr,:,:), res2) ;0,iyr+1
res3@gsnLeftString = ""
res2@gsnLeftString = ""
;plot0 = gsn_csm_contour(wks, svolavg1(0,6,:,:), res3);res3svolavg1
;plot0 = ColorNegDashZeroPosContour(plot0,"blue","transparent","red")
;ptest1 = gsn_csm_contour(wks, probsst1(0,6,:,:), sdres)
;ptest1 = gsn_contour_shade(ptest1,-999,90,opt) ;Pattern shade contours above 95.
ptest2 = gsn_csm_contour(wks, probu1(isea,iyr,:,:), sdres)
ptest2 = gsn_contour_shade(ptest2,-999,90,opt) ;Pattern shade contours above 95.
;overlay(plot(0), ptest1)
;overlay(plot(0), plot0)
overlay(plot(0), ptest2)
;************************************************
; use landsea_mask function to mask out land and then ocean data
;************************************************

     b = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")

     lsm_sst= landsea_mask(b->LSMASK,xsea&LAT,xsea&LON)

     xsea_land = mask(xsea,lsm_sst.eq.0, False)      ; mask out ocean points
     xsea_ocean = mask(xsea,lsm_sst.eq.1, False)   ; mask out ocean points

     copy_VarMeta(xsea,xsea_land)        ;建议把这两句加上，这样新的land 和 ocean都有
     copy_VarMeta(xsea,xsea_ocean)     ;了和原始数据 T一样的meta信息，并且海洋(或陆地)设为了却测值-999

res2@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
res2@cnMinLevelValF        = -2.0;280-273.16
res2@cnMaxLevelValF        =  2.0 ; 310-273.16
res2@cnLevelSpacingF       =  0.2

do ia = 0,2
	j = ia+1
res2@gsnLeftString = sign(j)+" "+AC(ia)
plot(j)  = gsn_csm_contour_map(wks, xsea2(ia,isea,iyr,:,:), res2) ;0,iyr+1
;res2@gsnPolar   = "NH"                          ; specify the hemisphere
;plot(j) = gsn_csm_contour_map_polar(wks, xsea2(ia,0,6,:,:), res2)    ;{0:90} create the plot

res3@gsnLeftString = ""
res2@gsnLeftString = ""
;plot0 = gsn_csm_contour(wks, xsea_ocean(ia,0,6,:,:) ,res3);xsea,res2
;plot0 = ColorNegDashZeroPosContour(plot0,"blue","transparent","red")

;ptest1 = gsn_csm_contour(wks, psea(ia,0,6,:,:), sdres)
;ptest1 = gsn_contour_shade(ptest1,-999,90,opt) ;Pattern shade contours above 95.
ptest2 = gsn_csm_contour(wks, psea2(ia,isea,iyr,:,:), sdres)
ptest2 = gsn_contour_shade(ptest2,-999,90,opt) ;Pattern shade contours above 95.
;overlay(plot(j), ptest1)
;overlay(plot(j), plot0)
overlay(plot(j), ptest2)
end do 


gsn_panel(wks, plot, (/2,2/), resk)
delete(res2)
end do 
;end do 
;delete([/res2,x,xx,xc,xa,TIME,yyyymm,ymstrt1,ymlast1,files,f,xavol,xseavol,xavolc,xavol1,xseavol1/])
;print(ik)
end