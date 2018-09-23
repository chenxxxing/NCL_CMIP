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
diri = "/data2/home/liufeiclim2014/xing/cmip_vol/"
diri2 = "/data2/data/cmip5/historical_regrid/atmos/mon/"
diri3 = "/data2/home/liufeiclim2014/xing/cmip_vol/zos/"
diri0 = "/data2/home/liufeiclim2014/xing/"
name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\;
            \
            "CNRM-CM5-2","CSIRO-Mk3-6-0", \;  "CNRM-CM5",
            "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\; 
             "GISS-E2-R","GISS-E2-R-CC","HadGEM2-ES", "MIROC-ESM-CHEM","MRI-ESM1","MIROC5",\;
              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "NorESM1-ME"/);"NorESM1-M",

;name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\
;       "CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM","CNRM-CM5-2",  \
;      "CNRM-CM5","CSIRO-Mk3-6-0",  \
;        "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
;        "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES","MIROC5","MIROC-ESM-CHEM",\
;         "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\;"bcc-csm1-1-m","bcc-csm1-1",

NLAT = 5;40;70;5;
SLAT = -5;-40;-70;-5;
WLON =  120;0;190 ;
ELON = 285;360;240 ;
  ymstrt= 187001  
  ymlast= 200412
  yrstrt= 1870  
  yrlast= 2004
nm = dimsizes(name_model)
nyr = yrlast-yrstrt+1

LEV = 3
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

f0 = addfile(diri0+"sst.mnmean.v4.1854-2016.nc","r")
TIME0=f0->time
yyyymm0=cd_calendar(TIME0, -1)
ymstrt0 =ind(yyyymm0.eq.ymstrt)
ymlast0 = ind(yyyymm0.eq.ymlast)
ss = f0->sst(ymstrt0:ymlast0,{SLAT:NLAT},{WLON:ELON})
ntime = dimsizes(ss&time)
nlat = dimsizes(ss&lat)
nlon = dimsizes(ss&lon)

ny1 = 5 ;1 ; year before eruption
ny2 = 5 ;4 ; year after eruption
nyr = ny1+ny2
nmon = ny1*12+ny2*12
;nyear = (/"year(-5)","year(-4)","year(-3)","year(-2)","year(-1)","year(0)","year(+1)",\
;         "year(+2)","year(+3)","year(+4)","year(+5)","year(+6)","year(+7)"/);"year(-6)",

files = systemfunc(" ls -1 "+diri2+"tos/"+name_model(0)+"*.nc")
;files = systemfunc(" ls -1 "+diri3+"zos_Omon_"+name_model(0)+"_historical*.nc")
f   = addfiles(files,"r") 
TIME=f[:]->time
yyyymm=cd_calendar(TIME, -1)
ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)
print(ymstrt1)
print(ymlast1)
x0  = f[:]->tos(ymstrt1(0):ymlast1(0),{SLAT:NLAT},{WLON:ELON});:,:)
;-------------;改变精度---------------------
lat = x0&lat
lon = x0&lon 
klon = ispan(0,nlon-1,1)*2.+WLON  
 klat = ispan(0,nlat-1,1)*2.+SLAT
xx0= linint2_Wrap(lon, lat, x0, True, klon, klat, 0) 
    xx0&LAT@units="degrees_north"
    xx0&LON@units="degrees_east"
printVarSummary(xx0)
ndim = dimsizes(xx0)
delete([/lat,lon/])

xa =  new((/nm,ndim(0),ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
xavol =   new((/nm,n,  nmon,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
xseavol = new((/nm,n, 4,nyr,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
xa_p =new((/nm,ndim(0),ndim(2)/), typeof(x0), x0@_FillValue) 
;--------------------------
do ik =0,nm-1
;files = systemfunc(" ls -1 "+diri+"rsds/rsds_Amon_"+name_model(ik)+"*.nc")
;files = systemfunc(" ls -1 "+diri2+"pr/"+name_model(ik)+".pr.historical.nc");.mon
files = systemfunc(" ls -1 "+diri2+"tos/"+name_model(ik)+".tos.historical.mon.nc")
;(0) 100000(1) 92500(2) 85000(3) 70000(4) 60000(5) 50000(6) 40000(7) 30000(8) 25000(9) 20000(10)  15000(11)  10000(12)  7000(13)  5000(14)  3000(15)  2000(16)  1000
;files = systemfunc(" ls -1 "+diri3+"zos_Omon_"+name_model(ik)+"_historical*.nc")


f   = addfiles(files,"r") 
TIME=f[:]->time
yyyymm=cd_calendar(TIME, -1)
ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)
print(ymstrt1)
print(ymlast1)
x  = f[:]->tos(ymstrt1(0):ymlast1(0),{SLAT:NLAT},{WLON:ELON});:,:);
;xx  = x 
;xx = xx*24*3600 ;prate-->precip
;if(x!1.eq."lat") then
;x1 = x(:,{SLAT:NLAT},{WLON:ELON});:,:);
;else 
;x@lat2d =  f[:]->lat
;x@lon2d =  f[:]->lon
;printVarSummary(x)
;    Opt                  = True
;   Opt@ForceOverwrite      = True
;   Opt@SrcInputFileName = f[:]
;   Opt@SrcTitle          = "NCEP Grid"
;   Opt@WgtFileName      = diri+"a.nc"   ; destination file name
;   Opt@InterpMethod      =    "patch"
;
;   Opt@DstGridType      = "1x1"              ; Destination grid
;   Opt@DstTitle         = "World Grid 1-degree Resolution"
;   Opt@DstLLCorner      = (/-89.75d,   0.00d /)
;   Opt@DstURCorner      = (/ 89.75d, 359.75d /) 
;    ;;Opt@PrintTimings   = True
;    ;;Opt@Debug          = True
;      zz = ESMF_regrid(x,Opt) 
;      printVarSummary(zz)
;x1  = zz(:,{SLAT:NLAT},{WLON:ELON})
;    delete(Opt)
;delete(zz)
;end if
;    system(" rm "+diri+"source_grid_file.nc" )
;    system(" rm "+diri+"destination_grid_file.nc" )
;    system(" rm "+diri+"a.nc" )
;printVarSummary(x1)
;-------------;改变精度---------------------
lat = x&lat
lon = x&lon 
klon = ispan(0,nlon-1,1)*2.+WLON  
 klat = ispan(0,nlat-1,1)*2.+SLAT
xx= linint2_Wrap(lon, lat, x, True, klon, klat, 0) 
    xx&LAT@units="degrees_north"
    xx&LON@units="degrees_east"
printVarSummary(xx)
;xx = xx*24*3600 ;prate-->precip
delete([/f,files,TIME,x,yyyymm,ymstrt1,ymlast1,lat,lon/]);x1,

xx= dtrend_msg_n(xx&time,xx,False,False,0)

ndim = dimsizes(xx)
;sc1 = clmMonTLL(sst1)
;sc2 = clmMonTLL(sst2)
xc  = clmMonTLL(xx)
xa(ik,:,:,:)  = calcMonAnomTLL(xx, xc)

;--------------------------------------------
TIMEs2=xa&time
yyyymms2=cd_calendar(TIMEs2, -1)
ymvol2 = new(nv, integer)
do iv = 0,n-1
ymvol2(iv) = ind(yyyymms2.eq.volmon(iv));EL(iv));
xavol(ik,iv,:,:,:) = xa(ik,ymvol2(iv)-ny1*12:ymvol2(iv)+ny2*12-1,:,:)
;xseavol(ik,iv,:,:,:,:) = month_to_seasonN(xavol(ik,iv,:,:,:), SEASON)
end do 
xavol!0 = "model"
xavol!1 = "vol"
;xseavol!0 = "model"
;xseavol!1 = "vol"
;printVarSummary(u)
;printVarSummary(va)
;printVarSummary(uavol)
print(ymvol2)
;xa_p(ik,:,:) = dim_avg_n_Wrap(xa(:,{-5:5},:), 1)
print(ik)
end do 
printVarSummary(xavol)
;-------------remove nearby annual---------------------
xavolc = new((/nm,nv,12,ndim(1),ndim(2)/), typeof(xa), xa@_FillValue)
xavol1 = new((/nm,nv,nmon,ndim(1),ndim(2)/), typeof(xa), xa@_FillValue)

do i=0,11
xavolc(:,:,i,:,:) = dim_avg_n_Wrap(xavol(:,:,i::12,:,:), 2)
end do 
xavolc!0 = "sample"
xavol1 = xavol
do i = 0,11
do j = 0,nyr-1
xavol1(:,:,j*12+i,:,:)  = xavol(:,:,i+j*12,:,:)-xavolc(:,:,i,:,:)
end do 
end do 
printVarSummary(xavol1)

xav_p      = dim_avg_n_Wrap(xavol1(:,:,:,:,:), 3)
;xseav_p = dim_avg_n_Wrap(xseavol1(:,:,:,:,:,:), 4)
printVarSummary(xav_p)

;---------------composite---------------------
xvolavg  = dim_avg_n_Wrap(xav_p, 1)     ;一个模式所有火山的合成
xmodavg = dim_avg_n_Wrap(xav_p, 0)      ;一个火山所有模式的合成
printVarSummary(xvolavg)
printVarSummary(xmodavg)

;------------------EL----------
cem = new((/nm,nmon,ndim(1),ndim(2)/),float)

do im=0,nm-1
ne = num(EL(im,:).not.1e+20)
ALLe1 = new((/ne,nmon,ndim(1),ndim(2)/), float, xa@_FillValue)
ALLe1_clm = new((/ne,12,dim(1),ndim(2)/),float)
ALLea1  = new((/ne,nmon,ndim(1),ndim(2)/),float)
ymvol3 = new(ne, integer)
do ie =0,ne-1
    ymvol3(ie) = ind(yyyymms3.eq.EL(im,ie));EL(iv));
    ALLe1(ie,:,:,:) = xa(im,ymvol2(ie)-ny1*12:ymvol2(ie)+ny2*12-1,:,:)
    do j=0,11
    ALLe1_clm(ie,j,:,:) = dim_avg_n_Wrap(ALLe1(ie,j::12,:,:), 1)
    end do 

    do i = 0,11
    do j = 0,nyr-1
    ALLea1(ie,j*12+i,:,:)  = ALLe1(ie,i+j*12,:,:)-ALLe1_clm(ie,i,:,:)
    end do 
    end do 

end do 
cem(:,im,:,:) = dim_avg_n_Wrap(ALLea1, 0)
delete([/ALLe1,ALLe1_clm,ALLea1,ymvol3/])
end do 
printVarSummary(cem)

ce_p = new((/nmon,ndim(2)/),float)
ce_p = dim_avg_n_Wrap(cem1, (/1,2/))
printVarSummary(ce_p)
;-------------t test ------------注意!!!样本是异常场!!总体也得是!!!!
avesumxv = new((/nm,nmon,ndim(2)/), float, xa@_FillValue)
varsumxv = new((/nm,nmon,ndim(2)/), float, xa@_FillValue)
avesumxm = new((/nv,nmon,ndim(2)/), float, xa@_FillValue)
varsumxm = new((/nv,nmon,ndim(2)/), float, xa@_FillValue)
        
avesumxm = 0 ;既然是异常场的平均值,那就应该是零啊!!
avesumxv = 0

varxm = dim_variance_n_Wrap(xav_p, 0)
varxv = dim_variance_n_Wrap(xav_p, 1)
;varu = dim_variance_n_Wrap(uavol, 0)
;varv = dim_variance_n_Wrap(vavol, 0)
siglvl=0.05
xa_p = dim_avg_n_Wrap(xa, 2)
do ik=0,nm-1
do iv=0,nv-1
do j =0,(ny1+ny2)-1
do i =0,11
varsumxv(ik,i+j*12,:) = dim_variance_n_Wrap(xa_p(ik,i::12,:), 0)
varsumxm(iv,i+j*12,:) = dim_variance_n_Wrap(xa_p(ik,i::12,:), 0)

;varsumu(i+j*12,:,:) = dim_variance_n_Wrap(ua(i::12,:,:), 0)
;varsumv(i+j*12,:,:) = dim_variance_n_Wrap(va(i::12,:,:), 0)

end do 
end do
end do 
end do 
printVarSummary(varsumxv)
printVarSummary(varsumxm)

                                                    ;所有sst不同月的方差

;sxep = ssumsf

ssum = 2004-1870
snm = nm
snv = nv
probv = xvolavg
probm = xmodavg

;probu = uvolavg
;probv = vvolavg

    iflag= True              ; population variance similar 

  probv=(1- ttest(xvolavg,varxv,snv, avesumxv,varsumxv,ssum, iflag, False))*100
  probm=(1- ttest(xmodavg,varxm,snm, avesumxm,varsumxm,ssum, iflag, False))*100
 ; probu=(1- ttest(uvolavg,varu,snp, avesumu,varsumu,ssum, iflag, False))*100
 ; probv=(1- ttest(vvolavg,varv,snp, avesumv,varsumv,ssum, iflag, False))*100

  copy_VarMeta(xvolavg,probv)
  copy_VarMeta(xmodavg,probm)

; copy_VarMeta(uvolavg,probu)
; copy_VarMeta(vvolavg,probv)  

printVarSummary(probm)
printVarSummary(probv)

;--------------plot-------------------
wks = gsn_open_wks("ps", diri+"cmip.composite.EL.sst.lon-time");+volyear(iv)););(eruptmon)
;wks = gsn_open_wks("ps", diri+"cmip.TP.volcano.ssta1.lon-time"+name_model(ik));."+volyear(iv)
gsn_define_colormap(wks,"BlueDarkRed18");"MPL_bwr""MPL_BrBG" "BrownBlue12""BlueDarkRed18""BlueDarkOrange18"
do iv=0,nv-1
;do ik=0,nm-1
 
;gsn_define_colormap(wks, "BlueWhiteOrangeRed")

plot = new(1, graphic)
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
 res@pmLabelBarOrthogonalPosF = -0.04       ; move label bar closer to axis
 res@cnLinesOn        = False               ; turn off contour lines
; res@cnFillPalette    = "BlueDarkRed18"    ; set color map
res@gsnLeftString     = " "
res@gsnRightString    = " "
res@tmYROn            = False   ;右侧坐标轴不显示？
 res@tiMainString     = "  "  ; title
res@lbLabelBarOn      = False
res@cnInfoLabelOn     = False
res@cnLineLabelsOn    = False
 res@vpXF             = 0.12                ; default is 0.2 (aspect ratio)
 res@vpYF             = 0.8                 ; default is 0.8
; res@vpHeightF        = 0.4                 ; default is 0.6
; res@vpWidthF         = 0.75                ; default is 0.6

res@gsnAddCyclic         = False        ; data not cyclic 
 ;res@vcGlyphStyle    = "CurlyVector"      ; turn on curly vectors
 res@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
res@cnMinLevelValF        = -0.5
res@cnMaxLevelValF        = 0.5
res@cnLevels = fspan(-0.5, 0.5, 20)

res@tmYLMode   = "Explicit"                                                
aa = ymvol2(nv-1);-ny1*12
bb = ymvol2(nv-1)+4*12-1;ny2*12-1 
tt = xa&time
res@tmYLValues = tt(aa:bb:6)
res@tmYLLabels = (/"Jan(-1)","Jul(-1)","Jan(0)","Jul(0)","Jan(+1)","Jul(+1)","Jan(+2)","Jul(+2)" /) ;"Jan(-2)","Jul(-2)",
   

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
 
        sdres@cnLevelSelectionMode = "ManualLevels"  
        sdres@cnMinLevelValF        = 90
      sdres@cnLevelSpacingF      =  5         ;间隔5个画
        sdres@cnLevels = (/80,90,95,99.9/)

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
;--------lon-time--------------------------------
res@tiMainString = "";
resk   = True
ys = (ny1-0)*12
yl = ys+4*12
;resk@txString = "SSTA(2S-2N) vol "+volyear(iv)+" "+name_model(ik) ;composite";"
;plot(0) = gsn_csm_hov(wks, xav_p(iv,:,:), res)
resk@txString = "composite sst "+volyear(iv);
plot(0) = gsn_csm_hov(wks, xmodavg(iv,ys:yl,:), res) ;xvolavg(ik,:,:)
;resk@txString = "composite precip "+name_model(ik);
;plot(0) = gsn_csm_hov(wks, xvolavg(ik,:,:), res) ;xvolavg(ik,:,:)
ptest = gsn_csm_hov(wks, probm(iv,ys:yl,:), sdres)
ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
overlay(plot(0), ptest)

resk@txFontHeightF = 0.03
resk@gsnPanelLabelBar   = True                 ; add common colorbar
resk@gsnMaximize = True
resk@gsnPaperOrientation = "portrait"   ; force portrait
;resk@lbBoxEndCapStyle = "TriangleBothEnds"
gsn_panel(wks, plot, (/1,1/), resk)
delete([/res/])
print(ik)
end do 
;end do 
end