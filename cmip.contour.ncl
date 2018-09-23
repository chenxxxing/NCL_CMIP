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
;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\;
;            \
;            "CNRM-CM5-2","CSIRO-Mk3-6-0", \;  "CNRM-CM5",
;            "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\; 
;             "GISS-E2-R","GISS-E2-R-CC","HadGEM2-ES", "MIROC-ESM-CHEM","MRI-ESM1","MIROC5",\;
;              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "NorESM1-ME"/);"NorESM1-M",
;
name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\
       "CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM","CNRM-CM5-2",  \
      "CNRM-CM5","CSIRO-Mk3-6-0",  \
        "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
        "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES","MIROC5","MIROC-ESM-CHEM",\
         "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\;"bcc-csm1-1-m","bcc-csm1-1",

NLAT = 40;90;70;5;5;
SLAT = -40;-90;-70;-5;-5;
WLON = 100;0; 190 ;
ELON = 300;360;240 ;
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
;-------------;改变精度---------------------
;lat = x0&lat
;lon = x0&lon 
;klon = ispan(0,nlon-1,1)*2.+WLON  
; klat = ispan(0,nlat-1,1)*2.+SLAT
;xx0= linint2_Wrap(lon, lat, x0, True, klon, klat, 0) 
;    xx0&LAT@units="degrees_north"
;    xx0&LON@units="degrees_east"
;printVarSummary(xx0)
;ndim = dimsizes(xx0)
;delete([/lat,lon/])
ndim = dimsizes(x0)

xa =  new((/nm,ndim(0),ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
xavol =   new((/nm,n,  nmon,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
xseavol = new((/nm,n, 4,nyr,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
xa_p =new((/nm,ndim(0),ndim(2)/), typeof(x0), x0@_FillValue) 
;--------------------------
do ik =0,nm-1
;files = systemfunc(" ls -1 "+diri+"rsds/rsds_Amon_"+name_model(ik)+"*.nc")
files = systemfunc(" ls -1 "+diri2+"pr/"+name_model(ik)+".pr.historical.nc");.mon
;files = systemfunc(" ls -1 "+diri2+"tos/"+name_model(ik)+".tos.historical.mon.nc")
;(0) 100000(1) 92500(2) 85000(3) 70000(4) 60000(5) 50000(6) 40000(7) 30000(8) 25000(9) 20000(10)  15000(11)  10000(12)  7000(13)  5000(14)  3000(15)  2000(16)  1000
;files = systemfunc(" ls -1 "+diri3+"zos_Omon_"+name_model(ik)+"_historical*.nc")


f   = addfiles(files,"r") 
TIME=f[:]->time
yyyymm=cd_calendar(TIME, -1)
ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)
print(ymstrt1)
print(ymlast1)
x  = f[:]->pr(ymstrt1(0):ymlast1(0),{SLAT:NLAT},{WLON:ELON});:,:);
xx  = x 
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
;;-------------;改变精度---------------------
;lat = x1&lat
;lon = x1&lon 
;klon = ispan(0,nlon-1,1)*2.+WLON  ;机智！
; klat = ispan(0,nlat-1,1)*2.+SLAT
;xx= linint2_Wrap(lon, lat, x1, True, klon, klat, 0) 
;    xx&LAT@units="degrees_north"
;    xx&LON@units="degrees_east"
;printVarSummary(xx)
;xx = xx*24*3600 ;prate-->precip
delete([/f,files,TIME,x,yyyymm,ymstrt1,ymlast1/]);,x1,lat,lon

xx= dtrend_msg_n(xx&time,xx,False,False,0)

ndim = dimsizes(xx)
;sc1 = clmMonTLL(sst1)
;sc2 = clmMonTLL(sst2)
xc  = clmMonTLL(xx)
xa(ik,:,:,:)  = calcMonAnomTLL(xx, xc)

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

;---------------composite---------------------
xvolavg  = dim_avg_n_Wrap(xavol1, 1)     ;一个模式所有火山的合成
xmodavg = dim_avg_n_Wrap(xavol1, 0)      ;一个火山所有模式的合成
;uvolavg  = dim_avg_n_Wrap(uavol, 0)     
;vvolavg  = dim_avg_n_Wrap(vavol, 0) 
;svolavg1  = dim_avg_n_Wrap(sseav_p, 0)     
;uvolavg1  = dim_avg_n_Wrap(useavol, 0)     
;vvolavg1  = dim_avg_n_Wrap(vseavol, 0) 
printVarSummary(xvolavg)
printVarSummary(xmodavg)

;-------------t test ------------注意!!!样本是异常场!!总体也得是!!!!
avesumxv = new((/nm,nmon,ndim(1),ndim(2)/), float, xa@_FillValue)
varsumxv = new((/nm,nmon,ndim(1),ndim(2)/), float, xa@_FillValue)
avesumxm = new((/nv,nmon,ndim(1),ndim(2)/), float, xa@_FillValue)
varsumxm = new((/nv,nmon,ndim(1),ndim(2)/), float, xa@_FillValue)
        
avesumxm = 0 ;既然是异常场的平均值,那就应该是零啊!!
avesumxv = 0

varxm = dim_variance_n_Wrap(xavol1, 0)
varxv = dim_variance_n_Wrap(xavol1, 1)
;varu = dim_variance_n_Wrap(uavol, 0)
;varv = dim_variance_n_Wrap(vavol, 0)
siglvl=0.05
do ik=0,nm-1
do iv=0,nv-1
do j =0,(ny1+ny2)-1
do i =0,11
varsumxv(ik,i+j*12,:,:) = dim_variance_n_Wrap(xa(ik,i::12,:,:), 0)
varsumxm(iv,i+j*12,:,:) = dim_variance_n_Wrap(xa(ik,i::12,:,:), 0)

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

;-----------------season-------------
xmsea = new((/nv,4,nyr,ndim(1),ndim(2)/), float, xa@_FillValue)
xvsea = new((/nm,4,nyr,ndim(1),ndim(2)/), float, xa@_FillValue)
probm_sea = new((/nv,4,nyr,ndim(1),ndim(2)/), float, xa@_FillValue)
probv_sea = new((/nm,4,nyr,ndim(1),ndim(2)/), float, xa@_FillValue)

do ik = 0,nm-1
do iv = 0,n-1
xmsea(iv,:,:,:,:) = month_to_seasonN(xmodavg(iv,:,:,:), SEASON)
xvsea(ik,:,:,:,:) = month_to_seasonN(xvolavg(ik,:,:,:), SEASON)
probv_sea(ik,:,:,:,:) = month_to_seasonN(probv(ik,:,:,:), SEASON)
probm_sea(iv,:,:,:,:) = month_to_seasonN(probm(iv,:,:,:), SEASON)
end do 
end do

printVarSummary(xvsea)
printVarSummary(xmsea)
printVarSummary(probv_sea)
printVarSummary(probm_sea)

xsea = dim_avg_n_Wrap(xvsea, 0)
psea = dim_avg_n_Wrap(probv_sea, 0)
printVarSummary(xsea)
;print(psea)
;--------------plot-------------------
;wks = gsn_open_wks("ps", diri+"obs1.composite.90%.TP.volcano.u10.lon-time."+ny1+ny2);(eruptmon)
;do ik=0,nm-1
;do iv=0,nv-1
;wks = gsn_open_wks("ps", diri+"cmip.composite.TP.volcano.tas.volavg"+name_model(ik));."+volyear(iv)
;wks = gsn_open_wks("ps", diri+"cmip.composite.TP.volcano.sst.modavg"+volyear(iv));"+name_model(ik));.
wks = gsn_open_wks("ps", diri+"cmip.composite.TP.volcano.pr");"+name_model(ik));.

gsn_define_colormap(wks,"BrownBlue12");"BlueDarkRed18""BrownBlue12""MPL_BrBG" "temp_19lev" 

;gsn_define_colormap(wks, "BlueWhiteOrangeRed")
do iyr = 5,8
plot = new(4, graphic)
res2= True
res2@gsnDraw = False
res2@gsnFrame = False
res2@mpCenterLonF = 180.
res2@gsnLeftString = " "
sdres = res2


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

vcres = res2  
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

;res2@mpProjection         = "Robinson" ;"Mollweide"
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
res2@cnLevels = fspan(-2.0, 2.0, 20)
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

vcres@vcMinDistanceF           = 0.017            ; thin vectors
;res2@cnFillPalette           = "gui_default"     ; set color map
vcres@vcRefAnnoOrthogonalPosF = -.535    ; move ref vector into plot
  vcres@vcRefAnnoFontHeightF    = 0.005
vcres@vcRefMagnitudeF          = 7.0             ; define vector ref mag
vcres@vcRefLengthF             = 0.1            ; define length of vec ref
vcres@vcRefAnnoOrthogonalPosF  = .1               ; move ref vector down
;vcres@vcRefLengthF             = 0.045            ; define length of vec ref
vcres@gsnAddCyclic         = False        ; data not cyclic 
vcres@gsnLeftString   = " "

ressf = True
ressf@gsnFrame = False
ressf@gsnDraw = False
ressf@cnLinesOn = True
ressf@cnFillOn  = False
ressf@cnLineLabelsOn = False
ressf@cnInfoLabelOn       = False
    ressf@gsnContourNegLineDashPattern = 1
 ;ressf@cnLineThicknessF = 1.5

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
do isea = 0,3

res2@gsnLeftString = SEASON(isea);name_model(ik) ;"A";MO(j); LABEL(isea)+
iyr1 = iyr-5
res2@gsnRightString = "year("+iyr1+")"
;resk@txString = "composite surface air "+name_model(ik)+" year("+iyr1+")";
;plot(isea)  = gsn_csm_contour_map_ce(wks, xvsea(ik,isea,iyr,:,:), res2) 
;ptest = gsn_csm_hov(wks, probv_sea(ik,isea,iyr,:,:), sdres)
resk@txString = "composite pr year("+iyr1+")";
plot(isea)  = gsn_csm_contour_map(wks, xsea(isea,iyr,:,:), res2) 
ptest = gsn_csm_hov(wks, psea(isea,iyr,:,:), sdres)
ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
overlay(plot(isea), ptest)
end do 



gsn_panel(wks, plot, (/4,1/), resk)
delete(res2)
end do 
;end do 
;delete([/res2,x,xx,xc,xa,TIME,yyyymm,ymstrt1,ymlast1,files,f,xavol,xseavol,xavolc,xavol1,xseavol1/])
;print(ik)
;end do 
end