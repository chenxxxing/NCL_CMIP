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
name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\
       "CESM1-BGC", "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM","CNRM-CM5-2",  \
      "CNRM-CM5","CSIRO-Mk3-6-0",  \
        "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
        "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES","MIROC5","MIROC-ESM-CHEM",\
         "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);
E_c = (/"ACCESS1-3","bcc-csm1-1-m","CESM1-FASTCHEM","CNRM-CM5-2","CNRM-CM5",\
  "GFDL-ESM2M","GISS-E2-H-CC", "MPI-ESM-LR"/)
L_c = (/"bcc-csm1-1","CCSM4","CESM1-BGC", "CESM1-CAM5-1-FV2","CESM1-WACCM",\
  "GFDL-ESM2G","GISS-E2-R","HadGEM2-AO","HadGEM2-ES","NorESM1-M"/)
N_c = (/"ACCESS1-0","CanESM2","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2-R-CC","MIROC5",\
  "MIROC-ESM-CHEM","MPI-ESM-MR","MPI-ESM-P","MRI-ESM1","NorESM1-ME"/)

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

NLAT =  5;
SLAT = -5;
WLON = 40;120; 190 ;
ELON = 280;285;240 ;
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

ny1 = 5 ; year before year0
ny2 = 5 ; year after year0
nyr = ny1+ny2
nmon = ny1*12+ny2*12
;nyear = (/"year(-5)","year(-4)","year(-3)","year(-2)","year(-1)","year(0)","year(+1)",\
;         "year(+2)","year(+3)","year(+4)","year(+5)","year(+6)","year(+7)"/);"year(-6)",
;--------------------------
f5 = addfile(diri+"cmip.vol.1000u.nc", "r")
xavol = f5->uavol(:,:,:,{SLAT:NLAT},{WLON:ELON})
;xavol = xavol*1000 ; hus kg/kg->g/kg
printVarSummary(xavol)
ndim = dimsizes(xavol)

f4 = addfile(diri+"cmip.vol.sst.nc", "r")
sst = f4->xavol(:,:,:,{SLAT:NLAT},{WLON:ELON})
;xavol = xavol*1000 ; hus kg/kg->g/kg
printVarSummary(sst)
;ndim = dimsizes(sst)
 
;time = f5->time
;yyyymm = cd_calendar(time, -1)
;print(yyyymm)
;printVarSummary(data)
;;--------------sort----------------------
xE = new((/nme,nv,nmon,ndim(3),ndim(4)/), float,xavol@_FillValue)
xL = new((/nml,nv,nmon,ndim(3),ndim(4)/), float,xavol@_FillValue)
xN = new((/nmn,nv,nmon,ndim(3),ndim(4)/), float,xavol@_FillValue)
sE = new((/nme,nv,nmon,ndim(3),ndim(4)/), float,sst@_FillValue)
sL = new((/nml,nv,nmon,ndim(3),ndim(4)/), float,sst@_FillValue)
sN = new((/nmn,nv,nmon,ndim(3),ndim(4)/), float,sst@_FillValue)

do i =0,nme-1
  ei = ind(name_model.eq.A(i))
xE(i,:,:,:,:) = xavol(ei,:,:,:,:)
sE(i,:,:,:,:) = sst(ei,:,:,:,:)
end do

do i =0,nml-1
  li = ind(name_model.eq.B(i))
xL(i,:,:,:,:) = xavol(li,:,:,:,:)
sL(i,:,:,:,:) = sst(li,:,:,:,:)

end do

do i =0,nmn-1
  ni = ind(name_model.eq.C(i))
xN(i,:,:,:,:) = xavol(ni,:,:,:,:)
sN(i,:,:,:,:) = sst(ni,:,:,:,:)
end do
printVarSummary(xE)

  xsort = xN
  ssort = sN
  nmsort = nmn
  k = "WEAK"

xsort1 = new((/nmsort,nv,nmon,ndim(3),ndim(4)/),float,xavol@_FillValue)
xsort1 = xavol@_FillValue
;if (mk.eq."EL") then 
nmsort1 = 0
do im=0,nmsort-1
  do iv=0,nv-1
  aa = month_to_season(ssort(im,iv,:,:,:), "DJF")
  nino34 = wgt_areaave_Wrap(aa(:,{-5:5},{190:240}),1.0, 1.0, 1)
    if (nino34(5).le.0.5)then
  xsort1(im,iv,:,:,:) = xsort(im,iv,:,:,:)
  nmsort1 = nmsort1+1
  end if 
delete([/aa,nino34/])
end do 
end do
print(nmsort1)
;else
;  xsort1 = xsort
;end if 

;-------------t test ------------注意!!!样本是异常场!!总体也得是!!!!
xva_1 = reshape(xavol, (/nm*nv,nmon,ndim(3),ndim(4)/))
xaa_1 = reshape(xsort1, (/nmsort*nv,nmon,ndim(3),ndim(4)/))

avesum = new((/nmon,ndim(3),ndim(4)/), float, xavol@_FillValue)
varsum = new((/nmon,ndim(3),ndim(4)/), float, xavol@_FillValue)
        
avesum = 0 ;既然是异常场的平均值,那就应该是零啊!!

varx = dim_variance_n_Wrap(xaa_1, 0)

siglvl=0.05
;do ik=0,nm*nv-1
;do j =0,(ny1+ny2)-1
;do i =0,11
varsum = dim_variance_n_Wrap(xva_1, 0) ;?!(i+j*12,:,:)(:,i::12,:,:)

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
sn   = nmsort1
print(sn)
xaa_2 = dim_avg_n_Wrap(xaa_1, 0)
prob = xaa_2

    iflag= True              ; population variance similar 

  prob=(1- ttest(xaa_2,varx,sn, avesum,varsum,ssum, iflag, False))*100

  copy_VarMeta(xavol(0,0,:,:,:),prob)
 
printVarSummary(prob)
;asciiwrite(diri+"cmip.prob.air"+AC(ia)+".txt", prob)

;---------------composite---------------------
xavg  = dim_avg_n_Wrap(xsort1, (/0,1/))     ;合成

xav_p      = dim_avg_n_Wrap(xavg(:,:,:), 1)
printVarSummary(xav_p)

prob_p      = dim_avg_n_Wrap(prob, 1)
prob_p1 = prob_p
;p1 = ndtooned(prob_p)
;miss = ind(p1.gt.110.)
;print(miss)
;p1(miss) = prob_p@_FillValue
;delete(miss)
;prob_p1 = reshape(p1, (/nmon,ndim(4)/))
;copy_VarMeta(prob_p, prob_p1)
printVarSummary(prob_p1)
;--------------plot-------------------
wks = gsn_open_wks("ps", diri+"cmip.composite.vol.1000u.lon-time.sort."+k)
;wks = gsn_open_wks("ps", diri+"cmip.TP.volcano.ssta1.lon-time"+name_model(ik));."+volyear(iv)
gsn_define_colormap(wks,"BlueWhiteOrangeRed");"BlueDarkRed18""precip_diff_12lev""MPL_BrBG" "BlueWhiteOrangeRed" 
;do iv=0,nv-1
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
res@cnMinLevelValF        = -1.0
res@cnMaxLevelValF        =  1.0
res@cnLevelSpacingF       = 0.1

res@tmYLMode   = "Explicit"                                                
aa = (ny1-1)*12;;
bb = (ny1+2)*12-1;ny2*12-1 
tt = xavol&time
res@tmYLValues = tt(aa:bb:6)
;es@tmYLValues = ispan(4*12, (4+4)*12-1, 6)
res@tmYLLabels = (/"Jan(-1)","Jul(-1)","Jan(0)","Jul(0)","Jan(+1)","Jul(+1)" /) ;,"Jan(+2)","Jul(+2)""Jan(-2)","Jul(-2)",
   

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
;--------lon-time--------------------------------
res@tiMainString = "";
resk   = True
ys = (ny1-1)*12
yl = ys+3*12
;resk@txString = "SSTA(2S-2N) vol "+volyear(iv)+" "+name_model(ik) ;composite";"
;plot(0) = gsn_csm_hov(wks, xav_p(iv,:,:), res)
resk@txString = "composite 1000hPa u in "+k
plot(0) = gsn_csm_hov(wks, xav_p(ys:yl,:), res) ;xvolavg(ik,:,:)
;resk@txString = "composite precip "+name_model(ik);
;plot(0) = gsn_csm_hov(wks, xvolavg(ik,:,:), res) ;xvolavg(ik,:,:)
ptest = gsn_csm_hov(wks, prob_p1(ys:yl,:), sdres)
ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
overlay(plot(0), ptest)

resk@txFontHeightF = 0.03
resk@gsnPanelLabelBar   = True                 ; add common colorbar
resk@gsnMaximize = True
resk@gsnPaperOrientation = "portrait"   ; force portrait
;resk@lbBoxEndCapStyle = "TriangleBothEnds"
gsn_panel(wks, plot, (/1,1/), resk)
delete([/res,ssort,xsort,xsort1/])
;print(ik)
;end do 
end