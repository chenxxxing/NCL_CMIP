load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
; Nino 1+2 (0-10S, 90W-80W), Nino 3 (5N-5S, 150W-90W), 
; Nino 3.4/ONI (5N-5S, 170W-120W) and Nino 4 (5N-5S, 160E-150W)
begin

name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC",\
  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM","CMCC-CESM","CMCC-CM","CMCC-CMS","CNRM-CM5-2",\
  "CNRM-CM5","CSIRO-Mk3-6-0", "EC-EARTH",  "FGOALS-g2","FIO-ESM","GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC",\
  "inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MIROC-ESM-CHEM", \
  "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "NorESM1-M","NorESM1-ME"/);  "MRI-CGCM3",,"GISS-E2-H"
;17 1752 18 1740 19 1740,"GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M" 23 1752 24 1752 25 1753 "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES",1860"MRI-ESM1",
nm = dimsizes(name_model)
;print("name_model="+nm)

volmon = (/188308,190210,191301,196303,198203,199106/);
volyear = (/1883,1902,1913,1963,1982,1991/)
volmo = (/1,1,1,1,1,1/);(/5,10,1,3,3,6/);
volmo = volmo-1
volyear = volyear
nv = dimsizes(volyear)


n3 =  asciiread("/work/CMIP5.nino34.mon.anomly.1854-2005.standardized.txt", (/33*152*12/), "float")
;n3 =  asciiread("/work/CMIP5.nino3.norm.1854-2005.txt", (/33*152*12/), "float")
nino3_obs =  asciiread("/work/obs.nino34.mon.anomly.1854-2005.standardized.txt", (/152*12/), "float")

tt =  asciiread("/work/CMIP5.TNI.norm.1854-2005.txt", (/33*152*12/), "float")
tt_obs = asciiread("/work/obs.TNI.norm.1854-2005.txt", (/152*12/), "float")
;tt_djf=asciiread("\work\cmip_TNI_JJA_1854_2005.txt", (/33*152/), "float");
;n3_djf=asciiread("\work\cmip_nino3_DJF_1854_2005.txt", (/33*152/), "float");
;tt_all=asciiread("\work\cmip_TNI_allyear_1854_2005.txt", (/33*152/), "float");
;ne_all=asciiread("\work\cmip_nino3_allyear_1854_2005.txt", (/33*152/), "float");

nino3 = onedtond(n3, (/nm,152*12/))
nino3!0 = "number"
nino3!1 = "time"
printVarSummary(nino3)
TNI = onedtond(tt, (/nm,152*12/))
TNI!0 = "number"
TNI!1 = "time"

n3_djf = new((/nm,152/), float)
n3_allyr = new((/nm,152/), float)
tt_jja = new((/nm,152/), float)
tt_allyr = new((/nm,152/), float)

do j=1,151
n3_djf(:,j) = (nino3(:,11+(j-1)*12)+nino3(:,0+j*12)+nino3(:,1+j*12))/3
end do 
n3_djf(:,0) = (nino3(:,0)+nino3(:,1))/2
;tt_djf(:,0) = (TNI(:,0)+TNI(:,1))/2
do j=0,151
tt_jja(:,j) = (TNI(:,5+j*12)+TNI(:,6+j*12)+TNI(:,7+j*12))/3
n3_allyr(:,j) = dim_avg_n_Wrap(nino3(:,j*12:j*12+11),1)
tt_allyr(:,j) = dim_avg_n_Wrap(TNI(:,j*12:j*12+11),1)
end do 

printVarSummary(tt_jja)
printVarSummary(n3_allyr)

;-------------trend--------------------
dt = 1.
dndjfdt = center_finite_diff_n(n3_djf, dt, False, 0, 1)
dtjjadt = center_finite_diff_n(tt_jja, dt, False, 0, 1)
dnalldt = center_finite_diff_n(n3_allyr, dt, False, 0, 1)
dtalldt = center_finite_diff_n(tt_allyr, dt, False, 0, 1)
;dnodt = center_finite_diff_n(nino3_obs, dt, False, 0, 0)
;dtodt = center_finite_diff_n(tt_obs, dt, False, 0, 0)
printVarSummary(dnalldt)

;---------volcano year-----------------
x = dndjfdt
y = dtjjadt
xa = dnalldt
ya = dtalldt
nino3_v = new((/nm,nv,15/),float)
tt_v =    new((/nm,nv,15/),float)
nino3a_v = new((/nm,nv,15/),float)
tta_v =    new((/nm,nv,15/),float)
nino3_ov =   new((/nv,15/),float)
tt_ov =      new((/nv,15/),float)
do i = 0,nv-1
  kk = volyear(i)-1854
  nino3_v(:,i,:) = x(:,kk-7+volmo(i):kk+7+volmo(i))
  tt_v(:,i,:) =    y(:,kk-7+volmo(i):kk+7+volmo(i))
  nino3a_v(:,i,:) = xa(:,kk-7+volmo(i):kk+7+volmo(i))
  tta_v(:,i,:) =    ya(:,kk-7+volmo(i):kk+7+volmo(i))
;  nino3_ov(i,:) =  n3(kk*12-60+volmo(i):kk*12+119+volmo(i))
;  tt_ov(i,:) =     dtodt(kk*12-60+volmo(i):kk*12+119+volmo(i))
end do 
printVarSummary(nino3_v)
;do i=0,nv-1
;asciiwrite("/work/cmip_nino3_DJF_trend_"+volyear(i)+"_yr1.txt", nino3_v(:,i,9))
;asciiwrite("/work/cmip_nino3_allyear_trend_"+volyear(i)+"_yr1.txt", nino3a_v(:,i,9))
;asciiwrite("/work/cmip_TNI_JJA_trend_"+volyear(i)+"_yr1.txt", tt_v(:,i,9))
;asciiwrite("/work/cmip_TNI_allyear_trend_"+volyear(i)+"_yr1.txt", tta_v(:,i,9))
;end do
opt=True
opt@fout="/work/cmip_vol_trend/cmip_nino3_DJF_trend_1.txt"
write_matrix(nino3_v(:,:,8), "33f7.2",opt) ; 2D array to be printed.
opt@fout="/work/cmip_vol_trend/cmip_nino3_allyear_trend_1.txt"
write_matrix(nino3a_v(:,:,8), "33f7.2",opt)
opt@fout="/work/cmip_vol_trend/cmip_TNI_JJA_trend_1.txt"
write_matrix(tt_v(:,:,8), "33f7.2",opt)
opt@fout="/work/cmip_vol_trend/cmip_TNI_allyear_trend_1.txt"
write_matrix(tta_v(:,:,8), "33f7.2",opt)
;-----------moving average---------------------
;nrun = 3
; endopt    = 1 ;起始结束不损失滑动的月份
;
;  nino3_v = runave_n_Wrap (nino3_v, nrun, endopt, 0)
;  nino3_ov = runave_n_Wrap (nino3_ov, nrun, endopt, 0)
;  tt_v = runave_n_Wrap (tt_v, nrun, endopt, 0)
;  tt_ov = runave_n_Wrap (tt_ov, nrun, endopt, 0)
;
;printVarSummary(tt_v)
;
;nino3_v_avg = dim_avg_n(nino3_v, 0)
;tt_v_avg    = dim_avg_n(tt_v, 0)
;--------plot-------------------------
do j=0,nv-1
wks1 = gsn_open_wks("pdf", "/work/cmip.trendyr0.nino34&TNI.xy"+volyear(j));

;wks2 = gsn_open_wks("pdf", "/work/obs.trend.vol-5-+10.nino34&TNI.mon.xy")
x_axis = ispan(1, 33, 1)
plot = new(2,graphic)
;do i=0,nm-1
res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnXYBarChart = True

res@vpWidthF    = 0.8        ; Make plots wider
res@vpHeightF   = 0.4
res@tmYROn                 = False    ; right off
res@gsnRightString = " "
;res@tiXAxisString   = " "
;res@trXMinF            = -mxlag+1
;res@trXMaxF            = mxlag-1
res@trYMinF            = -3
res@trYMaxF            = 3
;res@gsnYRefLine       = (/-0.13,0.13/)
;res@gsnXRefLine       = 0
res@gsnYRefLine       = 0
res@gsnAboveYRefLineColor = "coral"
res@gsnBelowYRefLineColor = "lightskyblue1"
res@xyLineColors = (/"red","blue"/)
;res@xyLineThicknessF = 2
 ; res@xyMonoDashPattern      = True     ; force all solid lines
colors = (/"coral",   "lightblue1", "blue",  "green","magenta2","navajowhite2"/)

resk = True
resk@gsnMaximize = True
resk@gsnPaperOrientation = "portrait"   ; force portrait
resk@lbBoxEndCapStyle = "TriangleBothEnds"
resk@txString = "nino34(DJF)&TNI(JJA) trend in "+volyear(j)+" year(0)" ;+name_model(i)" 

;res@lbFillColors = colors(j)
res@gsnLeftString = "nino3.4(DJF)"
  plot(0) = gsn_csm_xy(wks1,x_axis,nino3_v(:,j,7), res) ; ,nino3_v_avg(i,:),res2
res@gsnLeftString = "TNI(JJA)"
  plot(1) = gsn_csm_xy(wks1,x_axis,tt_v(:,j,7), res)  ;,tt_v_avg(i,:),res2
gsn_panel(wks1, plot, (/2,1/), resk)

;end do 

;do j = 0,nv-1
;res@gsnLeftString = "nino3.4"
;  plot(0) = gsn_csm_xy(wks2,x_axis,nino3_ov(:,:),res)  
;res@gsnLeftString = "TNI"
;  plot(1) = gsn_csm_xy(wks2,x_axis,tt_ov(:,:),res)  
; resk@txString = "nino34&TNI(trend) in volcano(observation)"
; 
;gsn_panel(wks2, plot, (/2,1/), resk)
end do 
end