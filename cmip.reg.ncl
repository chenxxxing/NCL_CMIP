load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"

begin
SEASON = (/"DJF","MAM","JJA","SON"/)
MO = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

diri = "/data2/home/liufeiclim2014/xing/cmip_vol/"
diri2 = "/data2/data/cmip5/historical_regrid/atmos/mon/"
;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\;
;            \
;            "CNRM-CM5-2","CSIRO-Mk3-6-0", \;  "CNRM-CM5",
;            "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\; 
;             "GISS-E2-R","GISS-E2-R-CC","HadGEM2-ES", "MIROC-ESM-CHEM","MRI-ESM1","MIROC5",\;
;              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "NorESM1-ME"/);"NorESM1-M",

name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\
       "CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM","CNRM-CM5-2",  \
      "CNRM-CM5","CSIRO-Mk3-6-0",  \
        "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
        "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES","MIROC5","MIROC-ESM-CHEM",\
         "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);

NLAT = 70;40;5;5;
SLAT = -70;-40;-5;-5;
WLON =  0;100;190 ;
ELON = 360;300;240 ;
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

ym1 = 1870
ym2 = 2004
ny2 = ym2-ym1+1
print(ny2)
nino3 =  asciiread(diri+"CMIP5.nino34.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")
tt_a = asciiread(diri+"CMIP5.nino34.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");wait for shortwave!
;air_a_c = asciiread(diri+"CMIP5.part.top.sw.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");wait for shortwave!
;air_a_c =  asciiread(diri+"CMIP5.part.SSHAI.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");借用
;gmst_a =  asciiread(diri+"CMIP5.part.RDI.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");借用

gmst_obs = asciiread(diri+"CMIP5.GMSST.org.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")
gmst_a = asciiread(diri+"CMIP5.GMSST.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")

air_a_c =  asciiread(diri+"CMIP5.GMST.anomaly.unnorm.1870-2004.txt", (/29,ny2*12/), "float");借用
air =  asciiread(diri+"CMIP5.GMST.org.unnorm.1870-2004.txt", (/29,ny2*12/), "float")
;;------------------是否标准化--------------------------
;nino3_a_obs= dim_standardize_n_Wrap(nino3_a_obs, 0, 0)
;nino3_a_obs1= dim_standardize_n_Wrap(nino3_a_obs1, 0, 0)
;tt_a_obs = dim_standardize_n_Wrap(tt_a_obs, 0, 0)
;tt_a_obs1= dim_standardize_n_Wrap(tt_a_obs1, 0, 0)
;gmst_a_obs = dim_standardize_n_Wrap(gmst_a_obs, 0, 0)
;gmst_a_obs1 = dim_standardize_n_Wrap(gmst_a_obs1, 0, 0)
;;air_a = dim_standardize_n_Wrap(air_a, 0, 0)

;---------是否去趋势---------------------------
air_dt   =dtrend(air_a_c, True)
printVarSummary(air_dt)

rcs = new((/nm,56,144/),float)
probs = new((/nm,56,144/),float)
do i=0,nm-1
;files = systemfunc(" ls -1 "+diri2+"tos/"+name_model(i)+"*.nc")
;f   = addfiles(files,"r") 
f = addfile(diri2+"tos/"+name_model(i)+".tos.historical.mon.nc", "r")
TIME=f->time
yyyymm=cd_calendar(TIME, -1)
ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)
print(ymstrt1)
print(ymlast1)
delete([/TIME,yyyymm/])
xx  = f->tos(ymstrt1(0):ymlast1(0),{SLAT:NLAT},{WLON:ELON});:,:)

xx= dtrend_msg_n(xx&time,xx,False,False,0)

printVarSummary(xx)
ndimsst = dimsizes(xx)
xc  = clmMonTLL(xx)
xa  = calcMonAnomTLL(xx, xc)

;------gmst index regression----------
X  =  air_dt(i,:) ;air1;NAT3;
printVarSummary(X)
;X2  =  AIRI2 ;NAT;

 tval =new((/ndimsst(1),ndimsst(2)/),typeof(xx)) ;12,
 df = new((/ndimsst(1),ndimsst(2)/),"integer") ;12,
tval = xx(0,:,:)

;do mon = 0,11
s_reg = xx(lat|:,lon|:,time|:)

SEA = 2   ;JJA
printVarSummary(X)
printVarSummary(s_reg)
rcs(i,:,:) = regcoef(X, s_reg, tval, df)

;end do
copy_VarCoords(xx(0,:,:), rcs(i,:,:))

df = df-2
b = tval               ; b must be same size as tval (and df) 
b = 0.5 

probs(i,:,:) = (1-(betainc(df/(df+tval^2),df/2.0,b)))*100 
delete([/f,xx,xc,xa,s_reg,tval,df,ymstrt1,ymlast1/]);files,
print(i)
end do 
rcs!0 = "model"
  ;
  copy_VarMeta(rcs,probs)

   printVarSummary(probs)
printVarSummary(rcs)

cr = dim_avg_n_Wrap(rcs, 0)
cp = dim_avg_n_Wrap(probs, 0)

;------plot---------
wks1 = gsn_open_wks("ps", diri+"reg.cmip.gmst.sst.em");+name_model(i))
plot = new(1,graphic)                ; create graphic array
gsn_define_colormap(wks1, "BlueDarkRed18");"BrownBlue12");"BlueGreen14");"BlueYellowRed");"sunshine_9lev")"BrownBlue12";
;do i=0,nm-1
res2 = True
res2@gsnDraw = False        ;dont draw
res2@gsnFrame = False        ;dont advance frame
res2@gsnLeftString = " "
res2@gsnRightString = " "
res = res2
sdres         =res2
opt  = res2
res3 = res2
res3@cnFillOn =  True;False
res3@cnLinesOn =False ;True ;
res3@lbLabelBarOn = False

res2@mpCenterLonF = 180.
res2@mpFillOn = False
res2@mpFillDrawOrder = "PostDraw"
  res2@mpMaxLatF                   = NLAT         ; choose subregion           
  res2@mpMinLatF                   = SLAT
  res2@mpMaxLonF                   = ELON
  res2@mpMinLonF                   = WLON
res2@lbLabelBarOn = False
res2@gsnAddCyclic = False         ;not 360 需要设置？
res=res2
;res2@cnFillPalette = 
res2@cnFillOn = True;False ;
res2@cnLinesOn = False;True ;
res2@cnLineLabelsOn = False
;res2@tiMainString = "differ "+LEVEL(LEV)+"hPa wind regressed in nino-A index 80-10"

;res2@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
;res2@cnMinLevelValF        = -0.4
;res2@cnMaxLevelValF        = 0.4
res2@gsnLeftStringFontHeightF = 0.02
;res2@cnFillMode = "RasterFill"
;res2@cnRasterSmoothingOn = True
;res2@vpHeightF = 0.7
;res2@vpWidthF = 0.7
res2@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
res2@cnMinLevelValF        = -0.8;280-273.16
res2@cnMaxLevelValF        =  2. ; 310-273.16
res2@cnLevels = fspan(-0.8, 2., 20)

sdres@cnFillOn = False
sdres@cnLinesOn = False
sdres@cnLineLabelsOn        = False; True   
        sdres@cnLevelSelectionMode = "ManualLevels"  
        sdres@cnMinLevelValF        = 90
      sdres@cnLevelSpacingF      =  5         ;间隔5个画
   ;     sdres@cnLevels = (/80,90,95,99.9/)
sdres@cnSmoothingOn = True
sdres@cnSmoothingTensionF = 1
sdres@cnInfoLabelOn = False
sdres@gsnLeftString = " "
res@vcMinDistanceF           = 0.017            ; thin vectors
;res@cnFillPalette           = "gui_default"     ; set color map
  res@vcRefMagnitudeF          = 0.3             ; define vector ref mag
  res@vcRefLengthF             = 0.045            ; define length of vec ref
    res@vcRefAnnoOrthogonalPosF  = .1               ; move ref vector down
;  res@gsnMaximize       = True                     ; maximize plot
ressf = True
ressf@gsnFrame = False
ressf@gsnDraw = False
ressf@cnLinesOn = True
ressf@cnFillOn  = False
ressf@cnLineLabelsOn = False
ressf@cnInfoLabelOn       = False
    ressf@gsnContourNegLineDashPattern = 1
ressf@gsnContourLineThicknessesScale = 0.3
ressf@gsnLeftString = " "
ressf@gsnRightString = "　　"

  opt@gsnShadeFillType = "pattern"      ; pattern fill
  opt@gsnShadeHigh = 17          ; 高于95的use pattern #17
  opt@gsnshadelow = -1       ;低于95为空白  
  opt@cnLineLabelsOn        = False   
  opt@gsnLeftString  = " "
  opt@cnInfoLabelOn  = False
;do i = 0,11
;res2@gsnLeftString = MO(i)
plot(0) = gsn_csm_contour_map_ce(wks1, cr, res2)
 ptest2 = gsn_csm_contour(wks1, cp, sdres)
ptest2 = gsn_contour_shade(ptest2,-999,95,opt) ;Pattern shade contours above 95.
;aplot(i) = gsn_csm_contour(wks1, rca(i,:,:), res2)
;vcplot = gsn_csm_vector(wks1, rcu(i,:,:), rcv(i,:,:), res)
;sfplot  = gsn_csm_contour(wks1, m(i,:,:), ressf)

;overlay(plot(0),ptest1)
;overlay(plot(i), vcplot)
overlay(plot(0), ptest2)

;end do

resk = True
resk@gsnMaximize = True
resk@gsnPanelLabelBar   = True                 ; add common colorbar
resk@gsnPaperOrientation = "portrait"   ; force portrait
resk@lbBoxEndCapStyle = "TriangleBothEnds"

resk@txString ="sst regressed onto GMST in CMIP5(95%) ";+name_model(i); ;NAT(MAY-MAR)"; 1950-1989
gsn_panel(wks1,plot,(/1,1/),resk)    
;delete([/files,f,xx,xc,xa,rcs,probs,s_reg,tval,df,ymstrt1,ymlast1/])
;print(i)
;end do 


end