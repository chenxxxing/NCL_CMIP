load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"

; Nino 1+2 (0-10S, 90W-80W), Nino 3 (5N-5S, 150W-90W), 
; Nino 3.4/ONI (5N-5S, 170W-120W) and Nino 4 (5N-5S, 160E-150W)
;(0) 100000 (1) 92500 (2) 85000 (3) 70000 (4) 60000 (5) 50000 (6) 40000 (7) 30000 (8) 25000 (9) 20000 (10)  15000 (11)  10000 (12)  7000 (13)  5000 (14)  3000 (15)  2000 (16)  1000 
begin
NLAT = 40;90;10;5;
SLAT = -40;-90;-10;-5;
WLON = 120;0; 190 ;
ELON = 300;360;240 ;
  ymstrt= 187001  
  ymlast= 200812
; yrstrt = 190001
; yrlast = 190012
LEV = 3
LEVEL = (/1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10/)

sign = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)" ,"(l)"/)
EL1=(/1876,1881,1888,1895,1896,1899,1902,1911,1918,1923,1925,1930,1939\
,1940,1941,1951,1953,1957,1968,1969,1972,1976,1977,1986,1987,1994,1997\
,2002/);1884,1904,1905,1914,1963,1965,1982,1991,,2009,2004,2006
EL=(/187601,188101,188801,189501,189601,189901,190201,191101,191801,192301,192501,193001,193901\
,194001,194101,195101,195301,195701,196801,196901,197201,197601,197701,198601,198701,199401,199701\
,200201/);1884,1904,1905,1914,1963,1965,1982,1991,,2009,200401,200601
ne = dimsizes(EL)
;volmon_tp = (/188308,190210,191301,196303,198203,199106/);
;volmon_nh = (/187503, 190703,191206, 193301,195603,198005/);NH
volmon_nh = (/187501, 190701,191201, 193301,195601,198001/);NH
;volmon_sh = (/188606,193204,199108,201106/);SH
volmon_sh = (/188601,193201,199101,201101/);SH
volmon_tp = (/188301,190201,196301,198201,199101/);189001,191301,
volyear_tp = (/1883,1902,1963,1982,1991/);1890,1913,

volyear = volyear_tp
volmon = volmon_tp;EL;+100
nv = dimsizes(volmon)

diri  = "/data2/home/liufeiclim2014/xing/" ;"/data/"; "/nuist/p/work/feil/xc/"
diri2 = "/data2/home/liufeiclim2014/xing/cmip_vol/";"/work/cmip_vol/";
;diri="/data/"
;f1=addfile("/data/uwnd.mon.mean.1948-2016.nc", "r")               ;from 1948 to 2016 dec 
;f2=addfile("/data/vwnd.mon.mean.1948-2016.nc", "r")  
;f1=addfile(diri+"uwnd.mon.mean.vc2.1851-2014.nc", "r")  
;f2=addfile(diri+"vwnd.mon.mean.vc2.1851-2014.nc", "r")   
f3= addfile(diri+"HadISST_1870-2016.nc", "r");sst.mnmean.v4.1854-2016.nc", "r");            
;f4= addfile(diri+"prate.mon.mean.vc2.1851-2014.nc", "r")
;f5 = addfile(diri+"air.sfc.mon.mean.vc2.1851-2014.nc", "r")
;f6 = addfile(diri+"uswrf.sfc.mon.mean.vc2.1851-2014.nc", "r")
;f7 = addfile(diri+"uswrf.ntat.mon.mean.vc2.1851-2014.nc", "r")
;f8 = addfile(diri+"dswrf.sfc.mon.mean.vc2.1851-2014.nc", "r")
;f9 = addfile(diri+"i20.anom.1871-2008.nc", "r")
;f10 = addfile(diri+"ssh.soda.1871-2008.TP.nc", "r")
;f11 = addfile("/data/slp(prmsl).mon.mean.1851-2014.nc" , "r" )
;f12 = addfile(diri+"hgt.mon.mean.vc2.1851-2014.nc" , "r" )
;f13 = addfile(diri+"shum.2m.mon.mean.v2.1851-2014.nc" , "r" )
MO = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
SEASON = (/"DJF","MAM", "JJA","SON" /)
;---------sst--------------
TIME1=f3->time
yyyymm=cd_calendar(TIME1, -1)
printMinMax(yyyymm, True)

ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)


LEV = 14 ;500hPa:10 200hPa:16 850hPa:3 50hPa:20 300hPa:14
sss = f3->sst(ymstrt1:ymlast1,:,:); ;LEV,
 sss2 = sss
;sss2 = sss*1000 ; shum kg/kg->g/kg
;sss2 = sss/100. ;slp Pa-->hPa
;sss2 = sss*24*3600 ;prate-->precip
sss1 = lonPivot(sss,0.5) ;HadISST center change
; printVarSummary(sss1)
;-------去掉神奇的-1000（不是missing value）-----
dimo = dimsizes(sss1)
sstone = ndtooned(sss1)
Y = ind(sstone.lt.-10);.eq.-1000)
sstone(Y) = sss1@_FillValue
sss2 = reshape(sstone, dimo)
copy_VarMeta(sss1, sss2)
printVarSummary(sss2)
sss2!0 = "time"
sss2!1 = "lat"
sss2!2 = "lon"
printVarSummary(sss2)

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

;;---------uv-----------------------------
;TIME2=f1->time
;yyyymm2=cd_calendar(TIME2, -1)
;
;ymstrt2 =ind(yyyymm2.eq.ymstrt)
;ymlast2 = ind(yyyymm2.eq.ymlast)
;
;u = f1->uwnd(ymstrt2:ymlast2,LEV,{SLAT:NLAT},{WLON:ELON})
;v = f2->vwnd(ymstrt2:ymlast2,LEV,{SLAT:NLAT},{WLON:ELON})
;u= dtrend_msg_n(u&time,u,False,False,0)
;v= dtrend_msg_n(v&time,v,False,False,0)
;uLAT=u&lat
;uLON=u&lon
;undim = dimsizes(u)
;uc  = clmMonTLL(u)
;vc  = clmMonTLL(v)
;ua  = calcMonAnomTLL(u, uc)
;va  = calcMonAnomTLL(v, vc)
;usea = month_to_seasonN(ua, SEASON)
;vsea = month_to_seasonN(va, SEASON)
;printVarSummary(u)
;----------------filter----------------------PAY ATTENTION! This is month data, not years!
xf = sa
;uf = ua(ik,:,:,:)
;vf = va(ik,:,:,:)
ihp   = 1       ; low-pass filter: ihp = 0; high-pass ihp = 1; band-pass ihp = 2.
sigma = 1.0                           ; Lanczos sigma
nWgt  = 71 ;必须是奇数 表示权重系数的总个数，越大效果越好，但是两端损失越多
fca   = 1./(10*12) ;                        
fcb   = -999;1./(2*12);1./8. *12  ;8——50年滤波
wgt  = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
xf = wgt_runave_n (sa, wgt, 0, 0 ) ;最后一个是指sst1第0维

n = nv

ny1 = 5 ; year before year0
ny2 = 5 ; year after year0
nyr = ny1+ny2
nmon = ny1*12+ny2*12
nyear = (/"year(-5)","year(-4)","year(-3)","year(-2)","year(-1)","year(0)","year(+1)",\
         "year(+2)","year(+3)","year(+4)","year(+5)","year(+6)","year(+7)"/);"year(-6)",
;nyear = (/"year(-3)","year(-2)","year(-1)","year(0)","year(+1)",\
;         "year(+2)","year(+3)"/);"year(-6)",

savol =   new((/n,  nmon,sndim(1),sndim(2)/), typeof(sa), sa@_FillValue) 
sseavol = new((/n, 4,nyr,sndim(1),sndim(2)/), typeof(sa), sa@_FillValue) 
;uavol =   new((/n,  12*(ny1+ny2),undim(1),undim(2)/), typeof(va), va@_FillValue) 
;vavol =   new((/n,  12*(ny1+ny2),undim(1),undim(2)/), typeof(va), va@_FillValue) 
;useavol = new((/n, 4,ny1+ny2,undim(1),undim(2)/), typeof(va), va@_FillValue) 
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
end do 
savol!0 = "sample"
sseavol!0 = "sample"
;printVarSummary(u)
;printVarSummary(va)
;printVarSummary(uavol)
;-------------remove nearby annual---------------------
savolc = new((/n,12,sndim(1),sndim(2)/), typeof(sa), sa@_FillValue)
savol1 = new((/n,nmon,sndim(1),sndim(2)/), typeof(sa), sa@_FillValue)
sseavol1 = new((/n, 4,nyr,sndim(1),sndim(2)/), typeof(sa), sa@_FillValue) 

do i=0,11
  savolc(:,i,:,:) = dim_avg_n_Wrap(savol(:,i::12,:,:), 1)
end do 
savolc!0 = "sample"
savol1 = savol
do i = 0,11
do j = 0,nyr-1
savol1(:,j*12+i,:,:)  = savol(:,i+j*12,:,:)-savolc(:,i,:,:)
end do 
end do 
printVarSummary(savol1)
do iv = 0,n-1
sseavol1(iv,:,:,:,:) = month_to_seasonN(savol1(iv,:,:,:), SEASON)
end do 
;-------------标准化？-------------------------
;ny3 = 2016-1870
;nino3_index =  asciiread(diri+"1obs.HadISST.nino34.anomaly.unnorm.1870-2016.txt", (/ny3*12/), "float")
;nino3_index!0 = "time"
;ndjf = month_to_season(nino3_index, "DJF")
;
;do iv=0, n-1
;  ;  do j=0,nmon-1
;
;    kk = (volyear(iv)-1870+1);(EL1(iv)-1870+1);
;        if(volyear(iv).eq.1883)then
;    savol1(iv,:,:,:) = savol1(iv,:,:,:)/ndjf(kk+1)
;    sseavol1(iv,:,:,:,:) = sseavol1(iv,:,:,:,:)/ndjf(kk+1) 
;    print(i)   
;  else
;    savol1(iv,:,:,:) = savol1(iv,:,:,:)/ndjf(kk);;VEI(i) ;vol_mag_sato(i) ;
;    sseavol1(iv,:,:,:,:) = sseavol1(iv,:,:,:,:)/ndjf(kk);;VEI(i) ;vol_mag_sato(i) ;
;  end if
;    print(ndjf(kk))
;end do 
;
;printVarSummary(sseavol1)
;---------------composite---------------------
svolavg  = dim_avg_n_Wrap(savol1, 0)     
;uvolavg  = dim_avg_n_Wrap(uavol, 0)     
;vvolavg  = dim_avg_n_Wrap(vavol, 0) 
svolavg1  = dim_avg_n_Wrap(sseavol1, 0)     
;uvolavg1  = dim_avg_n_Wrap(useavol, 0)     
;vvolavg1  = dim_avg_n_Wrap(vseavol, 0) 
printVarSummary(svolavg)
;-------------t test ------------注意!!!样本是异常场!!总体也得是!!!!
avesumsst = new((/12*(ny1+ny2),sndim(1),sndim(2)/), float, sst@_FillValue)
varsumsst = new((/12*(ny1+ny2),sndim(1),sndim(2)/), float, sst@_FillValue)
;avesumu =   new((/12*ny,undim(1),undim(2)/), float, ua@_FillValue)
;varsumu =   new((/12*ny,undim(1),undim(2)/), float, ua@_FillValue)
;avesumv =   new((/12*ny,undim(1),undim(2)/), float, va@_FillValue)
;varsumv =   new((/12*ny,undim(1),undim(2)/), float, va@_FillValue)

        
avesumsst = 0 ;既然是异常场的平均值,那就应该是零啊!!
avesumu = 0
avesumv = 0

varsst = dim_variance_n_Wrap(savol1, 0)
;varu = dim_variance_n_Wrap(uavol, 0)
;varv = dim_variance_n_Wrap(vavol, 0)

siglvl=0.05

do j =0,(ny1+ny2)-1
do i =0,11
varsumsst(i+j*12,:,:) = dim_variance_n_Wrap(sa(i::12,:,:), 0)
;varsumu(i+j*12,:,:) = dim_variance_n_Wrap(ua(i::12,:,:), 0)
;varsumv(i+j*12,:,:) = dim_variance_n_Wrap(va(i::12,:,:), 0)

end do 
end do
printVarSummary(varsumsst)   
                                                    ;所有sst不同月的方差

;sxep = ssumsf

ssum = 2009-1870
snp = n
probsst = svolavg
;probu = uvolavg
;probv = vvolavg

    iflag= True              ; population variance similar 

  probsst=(1- ttest(svolavg,varsst,snp, avesumsst,varsumsst,ssum, iflag, False))*100
 ; probu=(1- ttest(uvolavg,varu,snp, avesumu,varsumu,ssum, iflag, False))*100
 ; probv=(1- ttest(vvolavg,varv,snp, avesumv,varsumv,ssum, iflag, False))*100

  copy_VarMeta(svolavg,probsst)
; copy_VarMeta(uvolavg,probu)
; copy_VarMeta(vvolavg,probv)  

printVarSummary(probsst)
probsst1 = month_to_seasonN(probsst, SEASON)
;probu1 = month_to_seasonN(probu, SEASON)
;probv1 = month_to_seasonN(probv, SEASON)
;u2 = uvolavg1
;v2 = vvolavg1
;do i=0,3
;  do ii=0,ny-1
;  do j=0,undim(1)-1
;    do k=0,undim(2)-1
;if (.not.all(ismissing(probu1(i,ii,j,k))).and.probu1(i,ii,j,k).lt.95.and..not.all(ismissing(probv1(i,ii,j,k))).and.probv1(i,ii,j,k).lt.95)then 
;  u2(i,ii,j,k) = 0
;  v2(i,ii,j,k) = 0  
;end if 
;end do
;end do
;end do 
;end do 
;printVarSummary(v2)
;
;s2 = svolavg1
;do i=0,3
;  do ii=0,ny-1
;  do j=0,sndim(1)-1
;    do k=0,sndim(2)-1
;if (.not.all(ismissing(probsst1(i,ii,j,k))).and.probsst1(i,ii,j,k).lt.95)then 
;  s2(i,ii,j,k) = 0
;end if 
;end do
;end do
;end do 
;end do 
;----------plot------------------
;do iv=0,n-1
wks = gsn_open_wks("ps", diri2+"obs.vol.sst.mon.com(filter)");(eruptmon)diri+
;wks = gsn_open_wks("ps","/work/obs.volcano.sst.sea"+volyear(iv));(eruptmon) diri+

gsn_define_colormap(wks,"BlueDarkRed18");"BlueWhiteOrangeRed") ;"BlueDarkRed18");"MPL_bwr");"MPL_BrBG");"WhiteBlue"); precip_diff_12lev;BlueWhiteOrangeRed
;gsn_define_colormap(wks, "BrownBlue12")
;plot = new((/nm,12/),graphic)
plot = new((/3/),graphic)
;=========contour============================================
res2= True
res2@gsnDraw = False
res2@gsnFrame = False
res2@gsnLeftString = ""
res2@gsnRightString = ""

sdres = res2

vcres = res2  
  opt = res2
res2@vpWidthF    = 0.7        ; Make plots wider
res2@vpHeightF   = 0.3
res2@gsnLeftStringFontHeightF = 0.03
res2@gsnRightStringFontHeightF = 0.03
res2@cnLinesOn     = False
res2@cnFillOn             = True         ; turn on color fill
res2@gsnAddCyclic         = False        ; data not cyclic
res2@mpCenterLonF = 180.
res2@tmXMajorGrid                = True
res2@tmYMajorGrid                = True
res2@tmXMajorGridThicknessF      = 1.0
res2@tmYMajorGridThicknessF      = 1.0
res2@tmXMajorGridLineDashPattern = 2
res2@tmYMajorGridLineDashPattern = 2

  res2@mpMaxLatF                   = NLAT         ; choose subregion           
 res2@mpMinLatF                   = SLAT;30 ;
  res2@mpMaxLonF                   = ELON
  res2@mpMinLonF                   = WLON
;res2@mpProjection         = "Robinson" ;"Mollweide"
;  res2@mpPerimOn            = False     ;Turns off the box perimeter.
  res2@mpGridLatSpacingF    =  90                ; change latitude  line spacing
  res2@mpGridLonSpacingF    = 360.               ; change longitude line spacing
 ;res2@mpGridLineColor      = "transparent"      ; trick ncl into drawing perimeter
  res2@mpGridAndLimbOn      = True               ; turn on lat/lon lines
  res2@mpFillOn          = False;True ;

sdres@cnFillOn = False
sdres@cnLinesOn = False
sdres@cnLineLabelsOn        = False; True   
res2@cnInfoLabelOn            = False
  res2@cnLineLabelsOn           = False ;True
 
     ;   sdres@cnLevelSelectionMode = "ManualLevels"  
     ;   sdres@cnMinLevelValF        = 90
     ; sdres@cnLevelSpacingF      =  5         ;间隔5个画
     ;   sdres@cnLevels = (/80,90,95,99.9/)

sdres@cnSmoothingOn = True
sdres@cnSmoothingTensionF = 1
sdres@cnInfoLabelOn = False

res2@tiMainString = " "
  res2@lbLabelBarOn        = True;False;             ; turn off individual cb's
  res2@lbOrientation = "Vertical"
res2@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
res2@cnMinLevelValF        = -0.7;1-sst 4-pr 3-air 10-upsw.sfc 40-upsw.ntat 0-2downsw.sfc 6-850u 2-10u
res2@cnMaxLevelValF        =  0.7
res2@cnLevelSpacingF       = 0.07
;res2@cnLevels = fspan(-2.0, 2.0, 20)
res2@cnFillDrawOrder      = "Predraw"  ; draw contours first 可以用地图把陆地盖住
;res2@cnRasterSmoothingOn = True
;plot(i+12*j) = gsn_csm_contour_map_ce(wks2,sepavg(j*12+i,:,:),res2)
;plot(2*i+j) = gsn_csm_contour_map_ce(wks2,probep(j*12+i,:,:),res2)
  res2@tmXBLabelFontHeightF = 0.014  ;tick mark 坐标轴刻度！
  res2@tmYLLabelFontHeightF = 0.014

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
vcres@vcRefMagnitudeF          = 8.0             ; define vector ref mag
vcres@vcRefLengthF             = 0.1            ; define length of vec ref
vcres@vcRefAnnoOrthogonalPosF  = .1               ; move ref vector down
;vcres@vcRefLengthF             = 0.045            ; define length of vec ref
vcres@gsnAddCyclic         = False        ; data not cyclic 
vcres@gsnLeftString   = " "
vcres@gsnRightString   = " "

ressf = True
ressf@gsnFrame = False
ressf@gsnDraw = False
ressf@cnLinesOn = True
ressf@cnFillOn  = False
ressf@cnLineLabelsOn = False
ressf@cnInfoLabelOn       = False
    ressf@gsnContourNegLineDashPattern = 1
 ;ressf@cnLineThicknessF = 1.5
ressf@gsnLeftString   = " "
ressf@gsnRightString   = " "

res = True
;res@gsnPanelLabelBar   = True                 ; add common colorbar
res@txFontHeightF = 0.03
res@gsnMaximize = True
res@gsnPaperOrientation = "portrait"   ; force portrait
;res@lbBoxEndCapStyle = "TriangleBothEnds"

;  svolavg1  = where(abs(svolavg1).lt.0.3, svolavg1@_FillValue, svolavg1) ;设置白色（不想画）的值！！！
;=============SEASON====================================================
;do i=4,8
res2@gsnLeftString = "DJF";MO(j);"+"+lastmon ;name_model(i) ; 
res2@gsnRightString = "year(-1)" ;"year("+i+")"
plot(0)  = gsn_csm_contour_map(wks, svolavg1(0,5,:,:), res2) 
ptest = gsn_csm_contour(wks, probsst1(0,5,:,:), sdres)
ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
;;sfplot  = gsn_csm_contour(wks,sf(i,0,:,:),ressf)
;overlay(plot(j),vcplot)
overlay(plot(0), ptest)

res2@gsnLeftString = "JJA";MO(j);"+"+lastmon ;name_model(i) ; 
res2@gsnRightString = "year(0)" ;"year("+i+")"
plot(1)  = gsn_csm_contour_map(wks, svolavg1(2,5,:,:), res2) 
ptest = gsn_csm_contour(wks, probsst1(2,5,:,:), sdres)
ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
;;sfplot  = gsn_csm_contour(wks,sf(i,0,:,:),ressf)
;overlay(plot(j),vcplot)
overlay(plot(1), ptest)

;i=5
;do j=1,3
;kk = j
;;i1 = i-5
;res2@gsnRightString = nyear(i) ;"year("+i+")"
;;res2@gsnLeftString = " "
;lastmon = j+i*12
;res2@gsnLeftString = SEASON(j);sign(j)+MO(j);"+"+lastmon ;name_model(i) ; 
;;res@txString = "annual cycle index regressed onto sst in "+name_model(i);difference with observation"
;
;;plot(j)  = gsn_csm_contour_map_ce(wks, sseavol1(iv,j,i,:,:), res2) 
;plot(kk)  = gsn_csm_contour_map(wks, svolavg1(j,i,:,:), res2) 
;;res2@gsnPolar   = "NH"                          ; specify the hemisphere
;;  plot(kk) = gsn_csm_contour_map_polar(wks, svolavg1(j,i,{30:90},:), res2)    ; create the plot
;
;;vcplot = gsn_csm_vector(wks, u2(j,i,:,:),v2(j,i,:,:), vcres)
;ptest = gsn_csm_contour(wks, probsst1(j,i,:,:), sdres)
;ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
;;;sfplot  = gsn_csm_contour(wks,sf(i,0,:,:),ressf)
;;overlay(plot(j),vcplot)
;overlay(plot(kk), ptest)
;;plot(0)  = gsn_csm_contour_map_ce(wks, p0_jja, res2)
;;overlay(plot(0), ptest0)
;end do 
res2@gsnLeftString = "DJF";name_model(ik) ;"A";MO(j); LABEL(isea)+
res2@gsnRightString = "year(0)" ;"year("+i+")"
plot(2)  = gsn_csm_contour_map(wks, svolavg1(0,6,:,:), res2) 
;  plot(0) = gsn_csm_contour_map_polar(wks, svolavg1(0,6,{30:90},:), res2)    ;0,i+1 create the plot

ptest = gsn_csm_contour(wks, probsst1(0,6,:,:), sdres)
ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
overlay(plot(2), ptest)

;res@txString = "ssh in "+volyear(iv)+" event "+nyear(i);
res@txString ="(a) observation";+nyear(i);"surface down sw in "+volyear(iv)+" volcano event in "+nyear(i); "";
;(eruption mon) 
gsn_panel(wks,plot,(/3,1/),res);(without removing nearby annual cycle)

;;================MONTH=====================================================
;do i=2,3; vol:4,5 el(+-3):2,3
;do j=0,11
;kk = j
;;i1 = i-5
;res2@gsnRightString = nyear(i) ;"year("+i+")"
;;res2@gsnLeftString = " "
;lastmon = j+i*12
;res2@gsnLeftString = sign(j)+MO(j);SEASON(j);"+"+lastmon ;name_model(i) ; 
;;res@txString = "annual cycle index regressed onto sst in "+name_model(i);difference with observation"
;
;;plot(j)  = gsn_csm_contour_map_ce(wks, sseavol1(iv,j,i,:,:), res2) 
;plot(j)  = gsn_csm_contour_map(wks, svolavg(j+12*i,:,:), res2) 
;;res2@gsnPolar   = "NH"                          ; specify the hemisphere
;;  plot(kk) = gsn_csm_contour_map_polar(wks, svolavg1(j,i,{30:90},:), res2)    ; create the plot
;
;;vcplot = gsn_csm_vector(wks, u2(j,i,:,:),v2(j,i,:,:), vcres)
;ptest = gsn_csm_contour(wks, probsst(j+12*i,:,:), sdres)
;ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
;;;sfplot  = gsn_csm_contour(wks,sf(i,0,:,:),ressf)
;;overlay(plot(j),vcplot)
;overlay(plot(j), ptest)
;;plot(0)  = gsn_csm_contour_map_ce(wks, p0_jja, res2)
;;overlay(plot(0), ptest0)
;end do 
;
;;res@txString = "ssh in "+volyear(iv)+" event "+nyear(i);
;res@txString ="(a) observation";+nyear(i);"surface down sw in "+volyear(iv)+" volcano event in "+nyear(i); "";
;;(eruption mon) 
;gsn_panel(wks,plot,(/6,2/),res);(without removing nearby annual cycle)
;end do 

;end do 
;end do 

end