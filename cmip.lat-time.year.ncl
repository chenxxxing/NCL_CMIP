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

sign = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)" ,"(l)"/)

;volyear_nh_gao = (/1912, 1925,1943,1976,1980/);NH 1875, 1907,1933,1956,
;;volyear_sh_gao = (/1886 /);SH 把sh放入tp？！,1932,1991, 2011
;;vol_sh_mag_gao = (/1.93/);gao
;volyear_tp_gao = (/1883,1886,1902,1963,1982,1991/);1913,模式中貌似没有这个的火山强迫……
;
;volyear_nh_sato = (/1891,1912/);NH (1891or92)
;volyear_tp_sato = (/1883,1902,1963,1982,1991/);1896(sh) (1968),1886,1968,
;
;volmon_tp_sato = (/188301,190201,196301,198201,199101/);1896(sh) (1968),196801,188601,
volyear = (/1891,1912, 1925,1943,1976,1980/);volyear_tp_sato;ELyear;
volmon = (/189101,191201, 192501,194301,197601,198001/);volmon_tp_sato;EL;+100
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
delete([/f,files,TIME,yyyymm,ymstrt1,ymlast1/]);

;xa =  new((/nm,ndim(0),17,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
xavol =   new((/nm,n,  nmon,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) ;数组占内存太大就会读取不了文件啊！！！

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
x  = f[:]->ta(ymstrt1(0):ymstrt1(0)+1619,LEV,{SLAT:NLAT},{WLON:ELON});:,:);ymstrt1(0)+1619;ymlast1(0)

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

;ndim = dimsizes(xx)
;sc1 = clmMonTLL(sst1)
;sc2 = clmMonTLL(sst2)
xc  = clmMonTLL(xx)
xa  = calcMonAnomTLL(xx, xc);(ik,:,:,:,:)
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
;uf = wgt_runave_n (ua(ik,:,:,:), wgt, 0, 0 )
;vf = wgt_runave_n (va(ik,:,:,:), wgt, 0, 0 )

;---------------------------------------------
TIMEs2=x0&time ;xa用x0？？！！提取出来之后time的ind应该通用吧！！！
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
delete(xa)
end do 
printVarSummary(xavol)

;-------------remove nearby annual---------------------
xavolc = new((/nm,nv,12,ndim(1),ndim(2)/), typeof(xavol), xavol@_FillValue)
xavol1 = new((/nm,nv,nmon,ndim(1),ndim(2)/), typeof(xavol), xavol@_FillValue)
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

;-----------------------------------------
NLAT = 90
SLAT = -90

xavol2_p = dim_avg_n_Wrap(xavol1(:,:,:,{SLAT:NLAT},:), 4)
printVarSummary(xavol2_p)

;--------------plot-------------------
wks = gsn_open_wks("ps", diri+"cmip.composite.vol.50ta.lat-time.year(NH)")
;wks = gsn_open_wks("ps", diri+"cmip.TP.volcano.ssta1.lon-time"+name_model(ik));."+volyear(iv)
gsn_define_colormap(wks,"BlueDarkRed18");"BlueDarkRed18""BlueWhiteOrangeRed"
 
;gsn_define_colormap(wks, "BlueWhiteOrangeRed")
do iv=0,nv-1

plot = new(nm, graphic)
res                  = True                ; plot mods desired
res@gsnDraw           = False
res@gsnFrame          = False
res@gsnLeftString = ""
sdres = res
vcres = res  
  opt = res
res@tmXMajorGrid                = True
res@tmYMajorGrid                = True
res@tmXMajorGridThicknessF      = 1.0
res@tmYMajorGridThicknessF      = 1.0
res@tmXMajorGridLineDashPattern = 2
res@tmYMajorGridLineDashPattern = 2
 res@cnFillOn         = True                ; color on
 res@lbLabelStride    = 4                   ; every other label
 ;res@lbOrientation    = "Vertical"          ; vertical label bar
 ;res@pmLabelBarOrthogonalPosF = -0.04       ; move label bar closer to axis
 res@cnLinesOn        = False  ;True;              ; turn off contour lines
; res@cnFillPalette    = "BlueDarkRed18"    ; set color map
res@gsnLeftString     = " "
res@gsnRightString    = " "
res@tmYROn            = False   ;右侧坐标轴不显示？
 res@tiMainString     = "  "  ; title
res@lbLabelBarOn      = False
res@cnInfoLabelOn     = False
res@cnLineLabelsOn    = False
; res@vpXF             = 0.8               ; default is 0.2 (aspect ratio)
; res@vpYF             = 0.8                 ; default is 0.8
res@vpHeightF        = 0.6                 ; default is 0.6
res@vpWidthF         = 0.8; 0.7;               ; default is 0.6

res@gsnAddCyclic         = False        ; data not cyclic 
 ;res@vcGlyphStyle    = "CurlyVector"      ; turn on curly vectors
 res@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
res@cnMinLevelValF        = -1.5
res@cnMaxLevelValF        =  1.5
res@cnLevelSpacingF       =  0.15
   

 res1 = True

 res1@gsnDraw = False
 res1@gsnFrame = False

res1@gsnLeftString = " "
res1@cnLinesOn = True
res1@cnFillOn  = False
res1@cnLineLabelsOn = False
res1@cnInfoLabelOn       = False
    res1@gsnContourNegLineDashPattern = 1
res1@cnLineColor  = "red"
;level@units = "hPa"
sdres@cnFillOn = False
sdres@cnLinesOn = False
sdres@cnLineLabelsOn        = False; True   
 
   ;     sdres@cnLevelSelectionMode = "ManualLevels"  
   ;     sdres@cnMinLevelValF        = 90
   ;   sdres@cnLevelSpacingF      =  5         ;间隔5个画
   ;     sdres@cnLevels = (/80,90,95,99.9/)

sdres@cnSmoothingOn = True
sdres@cnSmoothingTensionF = 1
sdres@cnInfoLabelOn = False
sdres@gsnLeftString = ""
sdres@gsnRightString = ""
  opt@gsnShadeFillType = "pattern"      ; pattern fill
  opt@gsnShadeHigh = 17          ; 高于95的use pattern #17
  opt@gsnshadelow = -1       ;低于95为空白  
  opt@cnLineLabelsOn        = False   
  opt@gsnLeftString  = " "
  opt@cnInfoLabelOn  = False
res@tiMainString = "";


res@tmXBMode   = "Explicit"                                                
;aa = ymvol2(0);-ny1*12
;bb = ymvol2(0)+2*12-1;ny2*12-1 
ys = (ny1-0)*12-1
yl = ys+2*12-6
tt = xavol&time
res@tmXBValues = tt(ys:yl:3)
res@tmXBLabels = (/"Jan(0)","Apr(0)","Jul(0)","Oct(0)", "Jan(+1)","Apr(+1)","Jul(+1)" /) ;,"Oct(+1)","Jan(+2)","Jul(+2)""Jan(-2)","Jul(-2)",

;--------lat-time--------------------------------

;res@gsnLeftString = "(a) observation"
;res@gsnRightString = ""


volyear_nh_gao = (/1912, 1925,1943,1976,1980/);NH 1875, 1907,1933,1956,

volyear_nh_sato = (/1891,1912/);NH (1891or92)

volyear_c = (/"1891","1912", "1925","1943","1976","1980"/);(/"1883", "1902", "1963", "1982","1991"/);"1913",
do im = 0,nm-1
res@gsnLeftString = name_model(im) ;sign(iv)+volyear_c(iv)
plot(im) = gsn_csm_lat_time(wks, xavol2_p(model|im,vol|iv,LAT|:,time|ys:yl), res) ;xvolavg(ik,:,:)
end do 



resk   = True
resk@txString = volyear_c(iv);"composite erupted events";ENSO removed   hus(g/kg)

resk@txFontHeightF = 0.03
;resk@gsnPanelLabelBar   = True                 ; add common colorbar

resk@gsnMaximize = True
resk@gsnPaperOrientation = "portrait"   ; force portrait
resk@lbLabelFontHeightF = 0.008
;resk@lbBoxEndCapStyle = "TriangleBothEnds"
gsn_panel(wks, plot, (/6,5/), resk)
;delete([/res/])
;print(iv)
end do 
end