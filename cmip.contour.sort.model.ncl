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

WEAK = (/"ACCESS1-0","ACCESS1-3","CanESM2","CSIRO-Mk3-6-0","GFDL-CM3","GISS-E2-R-CC", "HadGEM2-AO",\
 "MIROC5","MIROC-ESM-CHEM","MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P","MRI-ESM1"/)
;LA ="HadGEM2-AO",
;EL ="ACCESS1-3","CanESM2","CSIRO-Mk3-6-0","GISS-E2-R-CC","MPI-ESM-LR","MPI-ESM-P","MRI-ESM1",
;N = "ACCESS1-0","GFDL-CM3","MIROC5","MIROC-ESM-CHEM", "MPI-ESM-MR"
W_EL =(/"ACCESS1-3","CanESM2","CSIRO-Mk3-6-0","GISS-E2-R-CC","MPI-ESM-LR","MPI-ESM-P","MRI-ESM1"/)
W_LA = (/"ACCESS1-0","GFDL-CM3","HadGEM2-AO","MIROC5","MIROC-ESM-CHEM", "MPI-ESM-MR"/)


STRONG = (/"bcc-csm1-1-m","bcc-csm1-1","CCSM4","CESM1-BGC", "CESM1-CAM5-1-FV2","CESM1-FASTCHEM",\
          "CESM1-WACCM","CNRM-CM5-2","CNRM-CM5","GFDL-ESM2G", "GFDL-ESM2M","GISS-E2-H-CC", \
          "GISS-E2-R", "HadGEM2-ES","NorESM1-M","NorESM1-ME"/)

S_LA = (/"bcc-csm1-1-m","CCSM4","CESM1-BGC", "CESM1-CAM5-1-FV2","CESM1-WACCM", "GFDL-ESM2G",\
 "GISS-E2-R", "HadGEM2-ES","NorESM1-M","NorESM1-ME"/)
S_EL = (/"bcc-csm1-1","CESM1-FASTCHEM", "CNRM-CM5-2","CNRM-CM5", "GFDL-ESM2M","GISS-E2-H-CC"/)

A_EL = (/"ACCESS1-3","CanESM2","CSIRO-Mk3-6-0","GISS-E2-R-CC","MPI-ESM-LR","MPI-ESM-P",\
        "MRI-ESM1","bcc-csm1-1","CESM1-FASTCHEM", "CNRM-CM5-2","CNRM-CM5", "GFDL-ESM2M","GISS-E2-H-CC"/)
A_LA = (/"ACCESS1-0","GFDL-CM3","HadGEM2-AO","MIROC5","MIROC-ESM-CHEM", "MPI-ESM-MR","bcc-csm1-1-m",\
        "CCSM4","CESM1-BGC", "CESM1-CAM5-1-FV2","CESM1-WACCM", "GFDL-ESM2G",\
        "GISS-E2-R", "HadGEM2-ES","NorESM1-M","NorESM1-ME"/)
;E = (/1,2,8,10,11,15,16,23/)
;L = (/3,5,6,7,9,14,17,19,20,27/)
;N = (/0,4,12,13,18,21,22,24,25,26,28/)

A = W_EL
B = W_LA
C = WEAK

nm = dimsizes(name_model)
nme = dimsizes(A)
nml = dimsizes(B)
nmn = dimsizes(C)

NLAT = 90;40;70;5;5;
SLAT = -90;-40;-70;-5;-5;
WLON = 100;120; 190 ;
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

MO = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
SEASON = (/"DJF","MAM", "JJA","SON" /)

ny1 = 5 ; year before year0
ny2 = 5 ; year after year0
nyr = ny1+ny2
nmon = ny1*12+ny2*12

;--------------------------
f1 = addfile(diri+"cmip.vol.sst.nc", "r")
xavol = f1->xavol(:,:,:,{SLAT:NLAT},{WLON:ELON})
;xavol = xavol*1000 ; hus kg/kg->g/kg
printVarSummary(xavol)
ndim = dimsizes(xavol)

;;--------------sort----------------------


A = new((/3,nm,5/), "float", a1@_FillValue)
A(0,:,:) = a1
A(1,:,:) = a2
A(2,:,:) = a3
AC = (/"A","B","C"/)
do ia = 0,2
xxx = A(ia,:,:)
;uuu = ua1
;vvv = va1

xaa = xavol
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
      xaa(:,:,imon,ilat,ilon) = xxx
   ;   xaa2(:,:,imon,ilat,ilon) = xxx
   ;   ua1(:,:,imon,ilat,ilon) = cla
   ;   ua2(:,:,imon,ilat,ilon) = cel      
   ;   ua3(:,:,imon,ilat,ilon) = wel      
   ;   va1(:,:,imon,ilat,ilon) = cla
   ;   va2(:,:,imon,ilat,ilon) = cel
   ;   va3(:,:,imon,ilat,ilon) = wel
    end do 
end do 
end do 
      xaa = where(xaa.ne.0.0, xavol, xavol@_FillValue)
    ;  xaa2 = where(xaa2.ne.0.0, xavol2, xavol2@_FillValue)
    ;  ua1 = where(ua1.ne.0.0, uavol, ua@_FillValue)
    ;  ua2 = where(ua2.ne.0.0, uavol, ua@_FillValue)
    ;  ua3 = where(ua3.ne.0.0, uavol, ua@_FillValue)
    ;  va1 = where(va1.ne.0.0, vavol, va@_FillValue)
    ;  va2 = where(va2.ne.0.0, vavol, va@_FillValue)
    ;  va3 = where(va3.ne.0.0, vavol, va@_FillValue)
printVarSummary(xaa)



;-------------t test ------------注意!!!样本是异常场!!总体也得是!!!!
xva_1 = reshape(xavol, (/nm*nv,nmon,ndim(3),ndim(4)/))
xaa_1 = reshape(xaa, (/nm*nv,nmon,ndim(3),ndim(4)/))
;xva_2 = reshape(xavol2, (/nm*nv,nmon,ndim2(3),ndim2(4)/))
;xaa_2 = reshape(xaa2, (/nm*nv,nmon,ndim2(3),ndim2(4)/))

avesum = new((/nmon,ndim(3),ndim(4)/), float, xavol@_FillValue)
varsum = new((/nmon,ndim(3),ndim(4)/), float, xavol@_FillValue)
        
avesum = 0 ;既然是异常场的平均值,那就应该是零啊!!
;avesum2 = avesum
;varsum2 = varsum

varx = dim_variance_n_Wrap(xaa_1, 0)
;varx2 = dim_variance_n_Wrap(xaa_2, 0)

siglvl=0.05
;do ik=0,nm*nv-1
;do j =0,(ny1+ny2)-1
;do i =0,11
varsum = dim_variance_n_Wrap(xva_1, 0) ;?!(i+j*12,:,:)(:,i::12,:,:)
;varsum2 = dim_variance_n_Wrap(xva_2, 0) ;?!(i+j*12,:,:)(:,i::12,:,:)

;varsumu(i+j*12,:,:) = dim_variance_n_Wrap(ua(i::12,:,:), 0)
;varsumv(i+j*12,:,:) = dim_variance_n_Wrap(va(i::12,:,:), 0)

;end do
;end do 
;end do 
printVarSummary(varsum)

                                                    ;所有sst不同月的方差

;sxep = ssumsf

;ssum = 2004-1870
ssum = nm*nv
sn   = num(A(ia,:,:).ne.0.0)
print(sn)
xaa_11 = dim_avg_n_Wrap(xaa_1, 0)
prob = xaa_11

    iflag= True              ; population variance similar 

  prob=(1- ttest(xaa_11,varx,sn, avesum,varsum,ssum, iflag, False))*100

  copy_VarMeta(xavol(0,0,:,:,:),prob)
 
printVarSummary(prob)

;asciiwrite(diri+"cmip.prob.air"+AC(ia)+".txt", prob)
;---------------composite & season---------------------
xavg  = dim_avg_n_Wrap(xaa, (/0,1/))     ;合成
xsea = month_to_seasonN(xavg, SEASON);(time|:,lev|:,LON|:)
printVarSummary(xsea)
psea = month_to_seasonN(prob, SEASON)
printVarSummary(psea)
p1 = ndtooned(psea)
miss = ind(p1.gt.110.)
p1(miss) = psea@_FillValue
delete(miss)
psea1 = reshape(p1, (/4,nmon/12,ndim(3),ndim(4)/))
copy_VarMeta(psea, psea1)
printVarSummary(psea1)

;asciiwrite(diri+"cmip.prob.sea.air"+AC(ia)+".txt", psea1)
;--------------plot-------------------
;wks = gsn_open_wks("ps", diri+"obs1.composite.90%.TP.volcano.u10.lon-time."+ny1+ny2);(eruptmon)
;do ik=0,nm-1
;do iv=0,nv-1
;wks = gsn_open_wks("ps", diri+"cmip.composite.TP.volcano.tas.volavg"+name_model(ik));."+volyear(iv)
;wks = gsn_open_wks("ps", diri+"cmip.composite.TP.volcano.sst.modavg"+volyear(iv));"+name_model(ik));.
wks = gsn_open_wks("ps", diri+"cmip.composite.volcano.sst.year0.sort"+AC(ia));"+name_model(ik));.

gsn_define_colormap(wks,"BlueDarkRed18");"BlueDarkRed18""precip_diff_12lev""MPL_BrBG" "BlueWhiteOrangeRed" 

;gsn_define_colormap(wks, "BlueWhiteOrangeRed")
;do iyr = 4,8
plot = new(5, graphic)
res2= True
res2@gsnDraw = False
res2@gsnFrame = False
res2@gsnLeftString = ""
res2@gsnRightString = ""
sdres = res2
vcres = res2  
  opt = res2


res2@cnLinesOn     = False
res2@gsnAddCyclic         =  False ;True;        ; data not cyclic
  res2@mpMaxLatF                   = NLAT         ; choose subregion           
 res2@mpMinLatF                   = SLAT;30 ;
  res2@mpMaxLonF                   = ELON
  res2@mpMinLonF                   = WLON
res2@mpCenterLonF = 180.
res2@tmXMajorGrid                = True
res2@tmYMajorGrid                = True
res2@tmXMajorGridThicknessF      = 1.0
res2@tmYMajorGridThicknessF      = 1.0
res2@tmXMajorGridLineDashPattern = 2
res2@tmYMajorGridLineDashPattern = 2


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

;res2@mpProjection         = "Robinson" ;"Mollweide"
;  res2@mpPerimOn            = False     ;Turns off the box perimeter.
  res2@mpGridLatSpacingF    =  90                ; change latitude  line spacing
  res2@mpGridLonSpacingF    = 360.               ; change longitude line spacing
 ;res2@mpGridLineColor      = "transparent"      ; trick ncl into drawing perimeter
  res2@mpGridAndLimbOn      = True               ; turn on lat/lon lines
  res2@mpFillOn          = False;True ;

  res2@lbLabelBarOn        = False;True ;            ; turn off individual cb's
res2@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
res2@cnMinLevelValF        = -1.0;280-273.16
res2@cnMaxLevelValF        =  1.0 ; 310-273.16
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
res2@gsnLeftString = "(a)DJF";MO(j);"+"+lastmon ;name_model(i) ; 
res2@gsnRightString = "year(-1)" ;"year("+i+")"
plot(0)  = gsn_csm_contour_map_ce(wks, xsea(0,5,:,:), res2) 
ptest = gsn_csm_contour(wks, psea1(0,5,:,:), sdres)
ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
;;sfplot  = gsn_csm_contour(wks,sf(i,0,:,:),ressf)
;overlay(plot(j),vcplot)
overlay(plot(0), ptest)
i=5
do j=1,3
kk = j
;i1 = i-5
res2@gsnRightString = nyear(i) ;"year("+i+")"
;res2@gsnLeftString = " "
lastmon = j+i*12
res2@gsnLeftString = sign(j)+SEASON(j);MO(j);"+"+lastmon ;name_model(i) ; 
;res@txString = "annual cycle index regressed onto sst in "+name_model(i);difference with observation"

;plot(j)  = gsn_csm_contour_map_ce(wks, sseavol1(iv,j,i,:,:), res2) 
plot(kk)  = gsn_csm_contour_map_ce(wks, xsea(j,i,:,:), res2) 
;res2@gsnPolar   = "NH"                          ; specify the hemisphere
;  plot(kk) = gsn_csm_contour_map_polar(wks, svolavg1(j,i,{30:90},:), res2)    ; create the plot

;vcplot = gsn_csm_vector(wks, u2(j,i,:,:),v2(j,i,:,:), vcres)
ptest = gsn_csm_contour(wks, psea1(j,i,:,:), sdres)
ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
;;sfplot  = gsn_csm_contour(wks,sf(i,0,:,:),ressf)
;overlay(plot(j),vcplot)
overlay(plot(kk), ptest)
;plot(0)  = gsn_csm_contour_map_ce(wks, p0_jja, res2)
;overlay(plot(0), ptest0)
end do 
res2@gsnLeftString = "(e)DJF";name_model(ik) ;"A";MO(j); LABEL(isea)+
res2@gsnRightString = "year(0)" ;"year("+i+")"
plot(4)  = gsn_csm_contour_map_ce(wks, xsea(0,6,:,:), res2) 
;  plot(0) = gsn_csm_contour_map_polar(wks, svolavg1(0,6,{30:90},:), res2)    ;0,i+1 create the plot

ptest = gsn_csm_contour(wks, psea1(0,6,:,:), sdres)
ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
overlay(plot(4), ptest)

;res@txString = "ssh in "+volyear(iv)+" event "+nyear(i);
;(eruption mon) 
resk@txString = "composite sst year(0) in "+AC(ia);hus(g/kg)
gsn_panel(wks,plot,(/5,1/),resk);(without removing nearby annual cycle)

;do isea = 1,3

;res2@gsnLeftString = SEASON(isea);name_model(ik) ;"A";MO(j); LABEL(isea)+
;iyr1 = iyr-5
;res2@gsnRightString = "year("+iyr1+")"
;resk@txString = "composite surface pr "+name_model(ik)+" year("+iyr1+")";
;plot(isea)  = gsn_csm_contour_map_ce(wks, xvsea(ik,isea,iyr,:,:), res2) 
;ptest = gsn_csm_hov(wks, probv_sea(ik,isea,iyr,:,:), sdres)
;kk = isea-1
;;plot(kk)  = gsn_csm_contour_map(wks, xsea(isea,iyr,:,:), res2) 
; res2@gsnPolar   = "NH"                          ; specify the hemisphere
;  plot(kk) = gsn_csm_contour_map_polar(wks, xsea(isea,iyr,{30:90},:), res2)    ; create the plot
;
;;vcplot = gsn_csm_vector(wks, usea(isea,iyr,:,:),vsea(isea,iyr,:,:), vcres)
;
;ptest = gsn_csm_hov(wks, psea1(isea,iyr,:,:), sdres)
;ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
;overlay(plot(kk), ptest)
;overlay(plot(kk), vcplot)
;end do 
;res2@gsnLeftString = "DJF";name_model(ik) ;"A";MO(j); LABEL(isea)+
;plot(0)  = gsn_csm_contour_map(wks, xsea(0,6,:,:), res2) ;0,iyr+1
;;plot(0) = gsn_csm_contour_map_polar(wks, xsea(0,6,{30:90},:), res2)
;;vcplot = gsn_csm_vector(wks, usea(0,iyr+1,:,:),vsea(0,iyr+1,:,:), vcres)
;;overlay(plot(3), vcplot)
;ptest = gsn_csm_hov(wks, psea1(0,6,:,:), sdres)
;ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
;overlay(plot(0), ptest)
;
;
;
;gsn_panel(wks, plot, (/5,1/), resk)
;end do 
;end do 
;delete([/res2,x,xx,xc,xa,TIME,yyyymm,ymstrt1,ymlast1,files,f,xavol,xseavol,xavolc,xavol1,xseavol1/])
;print(ik)
delete(res2)
end do 
end