load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
; Nino 1+2 (0-10S, 90W-80W), Nino 3 (5N-5S, 150W-90W), 
; Nino 3.4/ONI (5N-5S, 170W-120W) and Nino 4 (5N-5S, 160E-150W)
begin
MO = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
SEASON = (/"DJF","MAM", "JJA","SON" /)

name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\;
            "CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM",\;
            "CNRM-CM5-2",  "CNRM-CM5","CSIRO-Mk3-6-0",  \
             "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
             "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES",\
             "MIROC5","MIROC-ESM-CHEM",\
              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);;"HadGEM2-CC"(没有tas资料）,"FIO-ESM","FGOALS-g2",  "MRI-CGCM3",,"GISS-E2-H"

;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\;"bcc-csm1-1-m","bcc-csm1-1",
;       "CNRM-CM5-2",  \;"CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM",
;      "CNRM-CM5","CSIRO-Mk3-6-0",  \
;        "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
;        "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES","MIROC5","MIROC-ESM-CHEM",\
;         "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/); thetao
NLAT = 5;2;
SLAT = -5;-2;
WLON = 120;150; 
ELON = 285;155;;
dd   = 300;30 ; depth

  ymstrt= 187001  
  ymlast= 200412
  yrstrt= 1870  
  yrlast= 2004
nm = dimsizes(name_model)
nyr = yrlast-yrstrt+1

sign = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)" ,"(l)"/)
nyear = (/"year(-5)","year(-4)","year(-3)","year(-2)","year(-1)","year(0)","year(+1)",\
         "year(+2)","year(+3)","year(+4)","year(+5)","year(+6)","year(+7)"/);"year(-6)",

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

diri = "/data2/home/liufeiclim2014/xing/cmip_vol/";"/nuist/p/work/feil/xc/" ;"/work/cmip_vol/" ;
diri2 = "/data2/data/cmip5/historical_regrid/atmos/mon/"
diri3 = "/data2/home/liufeiclim2014/xing/thetao/";"/nuist/p/work/feil/xc/thetao/" ;"/work/cmip_vol/thetao/" ;
diri0 = "/data2/home/liufeiclim2014/xing/";"/data/" ; 
f0 = addfile(diri0+"sst.mnmean.v4.1854-2016.nc","r")
TIME0=f0->time
yyyymm0=cd_calendar(TIME0, -1)
ymstrt0 =ind(yyyymm0.eq.ymstrt)
ymlast0 = ind(yyyymm0.eq.ymlast)
ss = f0->sst(ymstrt0:ymlast0,{SLAT:NLAT},{WLON:ELON})
ntime = dimsizes(ss&time)
nlat = dimsizes(ss&lat)
nlon = dimsizes(ss&lon)

;elm = asciiread(diri+"cmip.enso.year.less.txt",(/nm*35/),"float")
;EL = onedtond(elm, (/nm,35/))

ny1 = 5 ;1 ; year before eruption
ny2 = 5 ;4 ; year after eruption
nyr = ny1+ny2
nmon = ny1*12+ny2*12
;nyear = (/"year(-5)","year(-4)","year(-3)","year(-2)","year(-1)","year(0)","year(+1)",\
;         "year(+2)","year(+3)","year(+4)","year(+5)","year(+6)","year(+7)"/);"year(-6)",
;
files = systemfunc(" ls -1 "+diri3+"thetao_10level_TP_*.nc")
print(files)
nfiles = dimsizes(files)

f00 = systemfunc(" ls -1 "+diri3+"thetao_10level_TP_ACCESS1-0.nc")
;files = systemfunc(" ls -1 "+diri3+"zos_Omon_"+name_model(0)+"_historical*.nc")
f1   = addfiles(f00,"r") 
TIME=f1[:]->time
yyyymm=cd_calendar(TIME, -1)
ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)
print(ymstrt1)
print(ymlast1)
x0  = f1[:]->thetao(ymstrt1(0):ymlast1(0),{0:dd},{SLAT:NLAT},{WLON:ELON});:,:);{0:300} ;
printVarSummary(x0)
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
printVarSummary(xx0)
delete([/f1,x0/]);lat,lon,
;xavol =   new((/nm,n,  nmon,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
;xseavol = new((/nm,n, 4,nyr,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
;xa_p =new((/nm,ndim(0),ndim(2)/), typeof(x0), x0@_FillValue) 
;--------------------------


;t =  new((/nfiles,ndim(0),ndim(1),ndim(2),ndim(3)/), typeof(xx0), xx0@_FillValue) 
tc =  new((/nfiles,12,ndim(1),ndim(2),ndim(3)/), typeof(xx0), xx0@_FillValue)
ta =  new((/nfiles,ndim(0),ndim(1),ndim(2),ndim(3)/), typeof(xx0), xx0@_FillValue) 

do i = 0,nfiles-1
f   = addfiles(files(i),"r") 
TIME=f[:]->time
yyyymm=cd_calendar(TIME, -1)
ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)
print(ymstrt1)
print(ymlast1)
if (name_model(i).eq."bcc-csm1-1-m".or.name_model(i).eq."bcc-csm1-1".or.name_model(i).eq."CESM1-BGC"\
   .or.name_model(i).eq."CESM1-CAM5-1-FV2".or.name_model(i).eq."CESM1-FASTCHEM".or.name_model(i).eq."CESM1-WACCM") then
  x1 = f[:]->thetao(240:1859,{0:dd}, {SLAT:NLAT},{WLON:ELON})
else

x1  = f[:]->thetao(ymstrt1(0):ymlast1(0),{0:dd},{SLAT:NLAT},{WLON:ELON});;
end if 
printVarSummary(x1)
;-------------;改变精度---------------------
lat = x1&lat
lon = x1&lon 
klon = ispan(0,nlon-1,1)*2.+WLON  
 klat = ispan(0,nlat-1,1)*2.+SLAT
xx= linint2_Wrap(lon, lat, x1, True, klon, klat, 0) 
    xx&LAT@units="degrees_north"
    xx&LON@units="degrees_east"
printVarSummary(xx)

delete([/f,TIME,x1,yyyymm,ymstrt1,ymlast1,lat,lon/]);

xx= dtrend_msg_n(xx&time,xx,False,False,0)

t = xx
tc(i,:,:,:,:) = clmMonTLLL(t)
ta(i,:,:,:,:) = calcMonAnomTLLL(t, tc(i,:,:,:,:))
print(i)
end do 
ta!0 = "model"
printVarSummary(ta)

t_p = dim_avg_n_Wrap(ta, 3)
ndimp = dimsizes(t_p)
printVarSummary(t_p)
;---------------composite----------------------
ny1 = 5
ny2 = 5
nyr = ny1+ny2
nmon = ny1*12+ny2*12

tv = new((/nfiles,nv,ny1*12+ny2*12,ndimp(2),ndimp(3)/), float,ta@_FillValue)
;te = new((/nfiles,ne,ny1*12+ny2*12,ndimp(1),ndimp(2)/), float,ta@_FillValue)

TIME=ta&time
yyyymm=cd_calendar(TIME, -1)
ymvol2 = new(nv, integer)
;ymel2  = new(ne, integer)

do iv =0,nv-1
  ;  kk = (EL(i)-ny1-1871)*12-1
  ;  ll = (EL(i)+ny2-1871)*12-2
    ymvol2(iv) = ind(yyyymm.eq.volmon(iv));EL(iv));
    tv(:,iv,:,:,:) = t_p(:,ymvol2(iv)-ny1*12:ymvol2(iv)+ny2*12-1,:,:);(kk:ll,:,:)

end do 
;do iv=0, ne-1
;   ; kk = (volyear(i)-ny1-1871)*12-1
;   ; ll = (volyear(i)+ny2-1871)*12-2
;    ymel2(iv) = ind(yyyymm.eq.ELmon(iv));EL(iv));
;    te(:,iv,:,:,:) = t_p(:,ymel2(iv)-ny1*12:ymel2(iv)+ny2*12-1,:,:);(kk:ll,:,:)
;end do 
printVarSummary(tv)
tv!1 = "sample"
tv!2 = "time"
;te!1 = "sample"
;te!2 = "time"
;----------------------------------------
tv_clm = new((/nfiles,nv,12,ndimp(2),ndimp(3)/),float,ta@_FillValue)
;te_clm = new((/nfiles,ne,12,ndimp(1),ndimp(2)/),float,ta@_FillValue)
do i=0,11
  tv_clm(:,:,i,:,:) = dim_avg_n_Wrap(tv(:,:,i::12,:,:), 2)
 ; te_clm(:,:,i,:,:) = dim_avg_n_Wrap(te(:,:,i::12,:,:), 2)

end do 
printVarSummary(tv_clm)

tva  = new((/nfiles,nv,nmon,ndimp(2),ndimp(3)/),float,ta@_FillValue)
;tea  = new((/nfiles,ne,nmon,ndimp(1),ndimp(2)/),float,ta@_FillValue)

do i = 0,11
do j = 0,nyr-1
tva(:,:,j*12+i,:,:)  = tv(:,:,i+j*12,:,:)-tv_clm(:,:,i,:,:)
;tea(:,:,j*12+i,:,:)  = te(:,:,i+j*12,:,:)-te_clm(:,:,i,:,:)

end do 
end do 
copy_VarMeta(tv, tva)
;copy_VarMeta(te, tea)
printVarSummary(tva)


;-----------------------------------------------------
;ce = dim_avg_n_Wrap(tea, 0)
cv = dim_avg_n_Wrap(tva, (/0,1/))
cv_sea = month_to_seasonN(cv(time|:,lev|:,LON|:), SEASON)
printVarSummary(cv)
printVarSummary(cv_sea)


cm_sea = new((/nv,4,nmon/12,ndimp(2),ndimp(3)/), float, ta@_FillValue)
do i=0,nv-1
  cm = dim_avg_n_Wrap(tva(:,i,:,:,:), 0)
  cm_sea(i,:,:,:,:) = month_to_seasonN(cm, SEASON)
end do 
printVarSummary(cm_sea)
;tva_sea = month_to_seasonN(tva(time|:,model|:,sample|:,lev|:,LON|:), SEASON)
;;tea_sea = month_to_seasonN(tea(time|:,sample|:,depth|:,lon|:), SEASON)
;printVarSummary(tva_sea)
;;ce_sea = dim_avg_n_Wrap(tea_sea, 2)
;cv_sea = dim_avg_n_Wrap(tva_sea, (/2,3/))

;------------------EL----------
;cem = new((/nm,nmon,ndimp(2),ndimp(3)/),float,ta@_FillValue)
;
;do im=0,nm-1
;ne = num(EL(im,:).gt.1e+20)
;ALLe1 = new((/ne,nmon,ndimp(2),ndimp(3)/), float, t_p@_FillValue)
;ALLe1_clm = new((/ne,12,ndimp(2),ndimp(3)/),float, t_p@_FillValue)
;ALLea1  = new((/ne,nmon,ndimp(2),ndimp(3)/),float, t_p@_FillValue)
;ymvol3 = new(ne, integer)
;do ie =0,ne-1
;    ymvol3(ie) = ind(yyyymms3.eq.EL(im,ie));EL(iv));
;    ALLe1(ie,:,:,:) = t_p(im,ymvol2(ie)-ny1*12:ymvol2(ie)+ny2*12-1,:,:)
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
;ce_p = new((/nmon,ndimp(2),ndimp(3)/),float)
;ce_p = dim_avg_n_Wrap(cem, (/0/))
;printVarSummary(ce_p)
;  ce_sea = month_to_seasonN(ce_p, SEASON)
;printVarSummary(ce_sea)
;-------------output---------------------------
;fout_name = "thetao_vol.nc"
;   system("/bin/rm -f "+diri+fout_name)   ; remove any pre-existing file
;   fout = addfile(diri+fout_name ,"c")  ; open output netCDF file
;
;   fAtt               = True            ; assign file attributes
;       fAtt@title         = "vertical thetao in CMIP5 ensemble"
;;       fAtt@source_file   =  "original-file.nc"
;       fAtt@Conventions   = "None"
;       fAtt@creation_date = systemfunc ("date")
;       fileattdef( fout, fAtt )            ; copy file attributes
;  ;    filedimdef(fout,"time",-1,True) 
;
;    fout->tva = tva
;
;;---Close file. Not necessary, but a good idea.
;  delete(fout)
;
;--------------plot-------------------
;do iv=0,nv-1
;wks = gsn_open_wks("ps", diri+"obs.ocean-temp.profile");+ny1+ny2);(eruptmon)+volyear(iv)
wks = gsn_open_wks("ps", diri+"cmip.composite.vol.ocean-temp");+volyear(iv)););."
gsn_define_colormap(wks,"BlueDarkRed18");"MPL_bwr""MPL_BrBG" "BlueDarkRed18"
;gsn_define_colormap(wks, "BlueWhiteOrangeRed")

plot = new(3, graphic)
res                  = True                ; plot mods desired
res@gsnDraw           = False
res@gsnFrame          = False
res@gsnLeftString = ""
sdres = res
vcres = res  
  opt = res

res@vpWidthF    = 0.7        ; Make plots wider
res@vpHeightF   = 0.3
res@gsnLeftStringFontHeightF = 0.03
res@gsnRightStringFontHeightF = 0.03
;res@tmXMajorGrid                = True
;res@tmYMajorGrid                = True
;res@tmXMajorGridThicknessF      = 1.0
;res@tmYMajorGridThicknessF      = 1.0
;res@tmXMajorGridLineDashPattern = 2
;res@tmYMajorGridLineDashPattern = 2
 res@cnFillOn         = True                ; color on
 res@lbLabelStride    = 4                   ; every other label
 ;res@lbOrientation    = "Vertical"          ; vertical label bar
 ;res@pmLabelBarOrthogonalPosF = -0.04       ; move label bar closer to axis
 res@cnLinesOn        = False ; True;              ; turn off contour lines
; res@cnFillPalette    = "BlueDarkRed18"    ; set color map
res@gsnLeftString     = ""
res@gsnRightString    = ""
res@tmYROn            = False   ;右侧坐标轴不显示？
 res@tiMainString     = ""  ; title
res@lbLabelBarOn      = False;True ;
;res@lbOrientation = "vertical"

res@cnInfoLabelOn     = False
res@cnLineLabelsOn    = False
 res@vpXF             = 0.12                ; default is 0.2 (aspect ratio)
 res@vpYF             = 0.8                 ; default is 0.8
; res@vpHeightF        = 0.4                 ; default is 0.6
; res@vpWidthF         = 0.75                ; default is 0.6

res@gsnAddCyclic         = False        ; data not cyclic 
 ;res@vcGlyphStyle    = "CurlyVector"      ; turn on curly vectors
res@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
res@cnMinLevelValF        = -0.3
res@cnMaxLevelValF        =  0.3
res@cnLevelSpacingF       =  0.03
;res@cnLevels = fspan(-0.5, 0.5, 20)

;res@tmYLMode   = "Explicit"                                                
;aa = ymvol2(nv-1)-ny1*12
;bb = ymvol2(nv-1)+ny2*12-1 
;tt = sa&time
;res@tmYLValues = tt(aa:bb:6)
;res@tmYLLabels = (/"Jan(-1)","Jul(-1)","Jan(0)","Jul(0)","Jan(+1)","Jul(+1)","Jan(+2)","Jul(+2)" /) 
   

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
 
     ;   sdres@cnLevelSelectionMode = "ManualLevels"  
     ;   sdres@cnMinLevelValF        = 90
     ; sdres@cnLevelSpacingF      =  5         ;间隔5个画
     ;   sdres@cnLevels = (/80,90,95,99.9/)

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
;--------lon-depth--------------------------------
res@tiMainString = "";
resk   = True
 res@trYReverse                      = True                  ; reverse y-axis
 res@tiYAxisString = "depth"
res@gsnLeftString = "DJF";MO(j);"+"+lastmon ;name_model(i) ; 
res@gsnRightString = "year(-1)" ;"year("+i+")"
plot(0)  = gsn_csm_contour(wks, cv_sea(0,5,:,:), res) 

res@gsnLeftString = "JJA";MO(j);"+"+lastmon ;name_model(i) ; 
res@gsnRightString = "year(0)" ;"year("+i+")"
plot(1)  = gsn_csm_contour(wks, cv_sea(2,5,:,:), res) 
;i=5
;do j=1,3
;kk = j
;;i1 = i-5
;res@gsnRightString = nyear(i) ;"year("+i+")"
;;res@gsnLeftString = " "
;lastmon = j+i*12
;res@gsnLeftString = SEASON(j);MO(j);"+"+lastmon ;name_model(i) ; 
;;res@txString = "annual cycle index regressed onto sst in "+name_model(i);difference with observation"
;
;;plot(j)  = gsn_csm_contour_map_ce(wks, sseavol1(iv,j,i,:,:), res) 
;plot(kk)  = gsn_csm_contour(wks, cv_sea(j,i,:,:), res) 
;;res@gsnPolar   = "NH"                          ; specify the hemisphere
;;  plot(kk) = gsn_csm_contour_map_polar(wks, svolavg1(j,i,{30:90},:), res)    ; create the plot
;
;;vcplot = gsn_csm_vector(wks, u2(j,i,:,:),v2(j,i,:,:), vcres)
;ptest = gsn_csm_contour(wks, prob_sea(j,i,:,:), sdres)
;ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
;;;sfplot  = gsn_csm_contour(wks,sf(i,0,:,:),ressf)
;;overlay(plot(j),vcplot)
;overlay(plot(kk), ptest)
;;plot(0)  = gsn_csm_contour_map_ce(wks, p0_jja, res)
;;overlay(plot(0), ptest0)
;end do 
res@gsnLeftString = "DJF";name_model(ik) ;"A";MO(j); LABEL(isea)+
res@gsnRightString = "year(0)" ;"year("+i+")"
plot(2)  = gsn_csm_contour(wks, cv_sea(0,6,:,:), res) 
;  plot(0) = gsn_csm_contour_map_polar(wks, svolavg1(0,6,{30:90},:), res)    ;0,i+1 create the plot


;res@txString = "ssh in "+volyear(iv)+" event "+nyear(i);
;(eruption mon) 
;resk@gsnFrame = False
resk@gsnMaximize = True
resk@gsnPanelLabelBar   = True                 ; add common colorbar
;resk@lbJustification = "BottomRight"

;resk@lbOrientation = "vertical"
resk@gsnPaperOrientation = "portrait"   ; force portrait
resk@txString = "CMIP5 (all)";+volyear(iv);
;resk@txFontHeightF
gsn_panel(wks,plot,(/3,1/),resk);(without removing nearby annual cycle)
end
