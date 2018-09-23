load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
; Nino 1+2 (0-10S, 90W-80W), Nino 3 (5N-5S, 150W-90W), 
; Nino 3.4/ONI (5N-5S, 170W-120W) and Nino 4 (5N-5S, 160E-150W)
begin
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
;         "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);thetao
;A_n = (/"CanESM2","GFDL-CM3","GFDL-ESM2M","GISS-E2-R-CC",\
; "MIROC5", "MPI-ESM-LR", "MPI-ESM-P",/) ;正相关 主要限定在二、四象限"NorESM1-M",,"HadGEM2-AO"
;B_n = (/"CNRM-CM5-2",\
;  "GFDL-ESM2G","GISS-E2-H-CC",  "NorESM1-ME" , "GISS-E2-R",\
; "MPI-ESM-MR"/) ;负相关 主要限定在一、三象限"bcc-csm1-1-m","bcc-csm1-1", "CESM1-CAM5-1-FV2",
;C_n = (/"ACCESS1-0","ACCESS1-3","CCSM4","CNRM-CM5","CSIRO-Mk3-6-0""HadGEM2-ES",\
;  "MIROC-ESM-CHEM","MRI-ESM1" /) ;杂乱"CESM1-BGC"6,"CESM1-WACCM"9,"CESM1-FASTCHEM",
nm = dimsizes(name_model)
EL=(/1876,1880,1884,1888,1895,1896,1899,1902,1904,1905,1911,1914,1918,1923,1925,1930,1939\
,1940,1941,1951,1953,1957,1963,1965,1968,1969,1972,1976,1977,1982,1986,1987,1991,1994,1997\
/);  没除掉vol年,2002,2004,2006,2009
;EL = (/1876,1880,1884,1895,1896,1918,1982,1997/);EP
;EL = (/1904,1914,1940,1941,1957,1963,1965,1968,1969,1977,1987,1991,1994,2002,2004,2009/);CP
;EL = (/1888,1899,1902,1905,1911,1923,1925,1930,1939,1951,1953,1972,1976,1986,2006/);MIX 
;EL=(/1881,1895,1896,1899,1904,1905,1911,1918,1923,1925,1930,1939\
;,1940,1941,1951,1953,1965,1968,1969,1972,1976,1977,1986,1987,1994,1997\
;,2002,2004,2006,2009/);  除去tp vol年
ne = dimsizes(EL)
diri = "/work/cmip_vol/"

;volmon = (/188308,190210,191301,196303,198203,199106/);
volyear = (/1883,1902,1963,1982,1991/);1886,1968, 1896(sh) (1968),
;vol_tp_mag_sato = (/0.1860,0.0931,0.1406,0.1338,0.1627/);0.0692,0.0491, sato 0.0557 (0.0491)

volyear_c = (/"1883", "1902","1963", "1982","1991"/);

nv = dimsizes(volyear)

volmo = volyear ;(/1,1,1,1,1,1/);(/5,10,1,3,3,6/);
volmo = 1
volmo = volmo-1

ymstrt = 1870
ymlast2 = 2004
ny2 = ymlast2-ymstrt+1
print(ny2)
nino3 =  asciiread(diri+"CMIP5.nino34.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")
;nino3 =  asciiread(diri+"CMIP5.nino34.anomaly.unnorm.1870-2004.19.txt", (/nm,ny2*12/), "float")
tt_a = asciiread(diri+"CMIP5.nino34.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");wait for shortwave!
;tt_a = asciiread(diri+"CMIP5.thetao100-200.nino34.anomaly.unnorm.1870-2004.19.txt", (/nm,ny2*12/), "float");wait for shortwave!
;air_a_c = asciiread(diri+"CMIP5.part.top.sw.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");wait for shortwave!

;gmst_obs = asciiread(diri+"CMIP5.GMSST.org.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")
gmst_a = asciiread(diri+"CMIP5.GMSST.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")

air_a_c =  asciiread(diri+"CMIP5.GMST.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")
;air_a_c =  asciiread(diri+"CMIP5.GMST.anomaly.unnorm.1870-2004.19.txt", (/nm,ny2*12/), "float")
;air =  asciiread(diri+"CMIP5.GMST.org.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")
nino3_a_obs1 =  asciiread("/work/1obs.HadISST.nino34.anomaly.unnorm.1870-2016.txt", (/(ny2+12)*12/), "float")
air_a_obs =  asciiread("/work/obs.air.anomaly.unnorm.1870-2014.txt", (/(ny2+10)*12/), "float")
t_a_obs =  asciiread("/work/obs.temp100-200.nino34.anom.unnorm.187102-2008.txt", (/(ny2+3)*12/), "float")
nobs = nino3_a_obs1(:ny2*12-1)
nobs!0 = "time"
aobs = air_a_obs(:ny2*12-1)
aobs!0 = "time"
tobs = aobs
tobs(0:11) = 0
tobs(12:) = air_a_obs(:(ny2-1)*12-1)
tobs!0 = "time"
printVarSummary(aobs)
printVarSummary(tobs)
;nobs = dim_standardize_n_Wrap(nobs, 0, 0)
;aobs = dim_standardize_n_Wrap(aobs, 0, 0)
;tobs = dim_standardize_n_Wrap(tobs, 0, 0)
;air_a_c = dim_standardize_n_Wrap(air_a_c, 0, 0)
;nino3 = dim_standardize_n_Wrap(nino3, 0, 0)
;tt_a = dim_standardize_n_Wrap(tt_a, 0, 0)
;gmst_a = dim_standardize_n_Wrap(gmst_a, 0, 0)
;-----------------TOA net SW averaged in tropic (30S-30N)----------------------
;rsdscs =  asciiread(diri+"CMIP5.rsdscs.anomaly.unnorm.1870-2004.22.txt", (/22,ny2*12/), "float")
;rsuscs =  asciiread(diri+"CMIP5.rsuscs.anomaly.unnorm.1870-2004.22.txt", (/22,ny2*12/), "float")
;raa = rsdscs 
;raa = -1*(rsdscs-rsuscs)
;printVarSummary(raa)
;;---------volcano year-----------------
;nyear1 = 5
;nyear2 = 5
;nyr = nyear1+nyear2
;nmon = nyear1*12+nyear2*12
;aod1 =    new((/nm,nv,nmon/),float)
;do i = 0,nv-1
;  kk = volyear(i)-yrstrt
;  aa = kk*12-nyear1*12+volmo(i);
;  bb = kk*12+nyear2*12-1+volmo(i);
;  aod1(:,i,:) =     raa(:,aa:bb)  
;end do 
;printVarSummary(aod1)
;aod2 = dim_avg_n_Wrap(aod1, 1)

;-------------------------------------------------------------------
nino3_a_obs = onedtond(nino3, (/nm,ny2*12/))
nino3_a_obs!0 = "number"
nino3_a_obs!1 = "time"
gmst_a_obs = onedtond(gmst_a, (/nm,ny2*12/))
gmst_a_obs!0 = "number"
gmst_a_obs!1 = "time"
air_a = onedtond(air_a_c, (/nm,ny2*12/))
air_a!0 = "number"
air_a!1 = "time"
tt_a_obs = onedtond(tt_a, (/nm,ny2*12/))
tt_a_obs!0 = "number"
tt_a_obs!1 = "time"

cor = escorc(air_a, gmst_a_obs) 
print(cor)

;---------------------------------------------
nyear1 = 5
nyear2 = 5
nyr = nyear1+nyear2
nmon = nyear1*12+nyear2*12
n3_all = new((/nm,nv,nyr/), float)
n_period = new((/nm,nv/), float)
a_all = new((/nm,nv,nyr/), float)
t_all = new((/nm,nv,nyr/), float)
ad_all = new((/nm,nv,2/), float)
do nn=0,nm-1;1;
nino3_om =nino3_a_obs(nn,:)
tt_om = tt_a_obs(nn,:)
g_om =gmst_a_obs(nn,:)
a_om =air_a(nn,:)  
nino3_a1 = nino3_a_obs(nn,:)
air_a1 = air_a(nn,:)
tt_a1 = tt_a_obs(nn,:)
;----------------filter----------------------PAY ATTENTION! This is month data, not years!
ihp   = 1       ; low-pass filter: ihp = 0; high-pass ihp = 1; band-pass ihp = 2.
sigma = 1.0                           ; Lanczos sigma
nWgt  = 51 ;必须是奇数 表示权重系数的总个数，越大效果越好，但是两端损失越多
fca   = 1./(9*12) ;                        
fcb   = -999;
wgt  = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
nino3_a1 = wgt_runave_n ( nino3_a_obs(nn,:), wgt, 0, 0 ) ;最后一个是指sst1第0维
tt_a1 = wgt_runave_n ( tt_a_obs(nn,:), wgt, 0, 0 ) ;最后一个是指sst1第0维
air_a1 = wgt_runave_n ( air_a(nn,:), wgt, 0, 0 ) ;最后一个是指sst1第0维
;-----------moving average---------------------
nrun = 3
 endopt    = 1 ;起始结束不损失滑动的月份

  nino3_om = runave_n_Wrap (nino3_a1, nrun, endopt, 0)
  t_om =    runave_n_Wrap (tt_a1, nrun, endopt, 0)
  g_om =     runave_n_Wrap (gmst_a_obs(nn,:), nrun, endopt, 0)
  a_om =runave_n_Wrap (air_a1, nrun, endopt, 0)
;printVarSummary(tt_om)

ndjf = month_to_season(nino3_om, "DJF")
;printVarSummary(ndjf)
;---------volcano year-----------------


n3_ov =    new((/nv,nmon/),float)
t_ov =       new((/nv,nmon/),float)
g_ov =        new((/nv,nmon/),float)
a_ov =       new((/nv,nmon/),float)
aobs_ov =   new((/nv,nmon/),float)
nobs_ov =   new((/nv,nmon/),float)
tobs_ov =   new((/nv,nmon/),float)
do i = 0,nv-1
  kk = volyear(i)-ymstrt
  aa = kk*12-nyear1*12+volmo(i)
  bb = kk*12+nyear2*12-1+volmo(i)

  t_ov(i,:)    =  t_om(aa:bb)
  n3_ov(i,:) =  nino3_om(aa:bb)
  
  g_ov(i,:) =     g_om(aa:bb)  

  a_ov(i,:) =     a_om(aa:bb)  
  aobs_ov(i,:) =     aobs(aa:bb)  

  nobs_ov(i,:) =     nobs(aa:bb) 
  tobs_ov(i,:) =     tobs(aa:bb) 
end do 
;nino3_ov_rm = dim_rmvmean_n_Wrap(nino3_ov, 1)
;nino3_ov_rm1 = dim_rmvmean_n_Wrap(nino3_ov1, 1)
;g_ov_rm = dim_rmvmean_n_Wrap(g_ov, 1)
;g_ov_rm1 = dim_rmvmean_n_Wrap(g_ov1, 1)
;;---------elnino year-----------------
;n3_e =   new((/ne,nmon/),float)
;t_e =      new((/ne,nmon/),float)
;g_e =       new((/ne,nmon/),float)
;a_e =       new((/ne,nmon/),float)
;
;do i = 0,ne-1
;  kk = EL(i)-ymstrt
;  aa = kk*12-nyear1*12
;  bb = kk*12+nyear2*12-1
;  n3_e(i,:) =  nino3_om(aa:bb)
;  t_e(i,:) =     t_om(aa:bb)
;  g_e(i,:) =       g_om(aa:bb)  
;  a_e(i,:) =       a_om(aa:bb)      
;end do 
;----------------------------------------
n3_clm = new((/nv,12/),float)
t_clm = new((/nv,12/),float)
g_clm = new((/nv,12/),float)
a_clm = new((/nv,12/),float)
nobs_clm = new((/nv,12/),float)
aobs_clm = new((/nv,12/),float)
tobs_clm = new((/nv,12/),float)
e_n3_clm = new((/ne,12/),float)
e_t_clm = new((/ne,12/),float)
e_g_clm = new((/ne,12/),float)
e_a_clm = new((/ne,12/),float)
do i=0,11
  n3_clm(:,i)  = dim_avg_n_Wrap(n3_ov(:,i::12), 1)
  t_clm(:,i)      = dim_avg_n_Wrap(t_ov(:,i::12), 1)
  g_clm(:,i)      = dim_avg_n_Wrap(g_ov(:,i::12), 1)
  a_clm(:,i)      = dim_avg_n_Wrap(a_ov(:,i::12), 1)
  nobs_clm(:,i)      = dim_avg_n_Wrap(nobs_ov(:,i::12), 1)
  aobs_clm(:,i)      = dim_avg_n_Wrap(aobs_ov(:,i::12), 1)  
  tobs_clm(:,i)      = dim_avg_n_Wrap(tobs_ov(:,i::12), 1)  
;  e_n3_clm(:,i)  = dim_avg_n_Wrap(n3_e(:,i::12), 1)
;  e_t_clm(:,i)      = dim_avg_n_Wrap(t_e(:,i::12), 1)
;  e_g_clm(:,i)      = dim_avg_n_Wrap(g_e(:,i::12), 1)
;  e_a_clm(:,i)      = dim_avg_n_Wrap(a_e(:,i::12), 1)
end do 
;printVarSummary(n3_clm)

n3_v_a  = new((/nv,nmon/),float)
t_v_a   = new((/nv,nmon/),float)
g_v_a   = new((/nv,nmon/),float)
a_v_a   = new((/nv,nmon/),float)
nobs_v_a   = new((/nv,nmon/),float)
aobs_v_a   = new((/nv,nmon/),float)
tobs_v_a   = new((/nv,nmon/),float)
n3_e_a  = new((/ne,nmon/),float)
t_e_a   = new((/ne,nmon/),float)
g_e_a   = new((/ne,nmon/),float)
a_e_a   = new((/ne,nmon/),float)
do i = 0,11
do j = 0,nyr-1
n3_v_a(:,j*12+i)  = n3_ov(:,i+j*12)-n3_clm(:,i)
t_v_a(:,j*12+i)   = t_ov(:,i+j*12)  -   t_clm(:,i)
g_v_a(:,j*12+i)   = g_ov(:,i+j*12)  -   g_clm(:,i)
a_v_a(:,j*12+i)   = a_ov(:,i+j*12)  -   a_clm(:,i)
nobs_v_a(:,j*12+i)   = nobs_ov(:,i+j*12)  -   nobs_clm(:,i)
aobs_v_a(:,j*12+i)   = aobs_ov(:,i+j*12)  -   aobs_clm(:,i)
tobs_v_a(:,j*12+i)   = tobs_ov(:,i+j*12)  -   tobs_clm(:,i)
;n3_e_a(:,j*12+i)  = n3_e(:,i+j*12)- e_n3_clm(:,i)
;t_e_a(:,j*12+i)   = t_e(:,i+j*12)  -   e_t_clm(:,i)
;g_e_a(:,j*12+i)   = g_e(:,i+j*12)  -   e_g_clm(:,i)
;a_e_a(:,j*12+i)   = a_e(:,i+j*12)  -   e_a_clm(:,i)
end do 
end do 
;printVarSummary(n3_v_a)

nino3_ov_avg =  dim_avg_n(n3_v_a, 0)
t_ov_avg =      dim_avg_n(t_v_a, 0)
g_ov_avg =      dim_avg_n(g_v_a, 0)
a_ov_avg =      dim_avg_n(a_v_a, 0)
;nino3_e_avg  =  dim_avg_n(n3_e_a, 0)
;t_e_avg  =      dim_avg_n(t_e_a, 0)
;g_e_avg  =      dim_avg_n(g_e_a, 0)
;a_e_avg  =      dim_avg_n(a_e_a, 0)


n3_v_a!0 = "vol"
n3_v_a!1 = "time"
n3_djf = new((/nv,nyr/),float)
a_v_a!0 = "vol"
a_v_a!1 = "time"
a_djf = new((/nv,nyr/),float)
t_v_a!0 = "vol"
t_v_a!1 = "time"
t_djf = new((/nv,nyr/),float)
nobs_v_a!0 = "vol"
nobs_v_a!1 = "time"
nobs_djf = new((/nv,nyr/),float)
aobs_v_a!0 = "vol"
aobs_v_a!1 = "time"
aobs_djf = new((/nv,nyr/),float)
tobs_v_a!0 = "vol"
tobs_v_a!1 = "time"
tobs_djf = new((/nv,nyr/),float)
do kk=0,nv-1
n3_djf(kk,:) = month_to_season(n3_v_a(kk,:), "DJF")
a_djf(kk,:) = month_to_season(a_v_a(kk,:), "DJF")
t_djf(kk,:) = month_to_season(t_v_a(kk,:), "DJF")
nobs_djf(kk,:) = month_to_season(nobs_v_a(kk,:), "DJF")
aobs_djf(kk,:) = month_to_season(aobs_v_a(kk,:), "DJF")
tobs_djf(kk,:) = month_to_season(tobs_v_a(kk,:), "DJF")
end do 
a_all(nn,:,:)  = a_djf
n3_all(nn,:,:) = n3_djf
t_all(nn,:,:) = t_djf
ad_all(nn,:,0) = a_v_a(:,6*12)-a_v_a(:,5*12)
ad_all(nn,:,1) = a_v_a(:,7*12)-a_v_a(:,6*12)
end do 
a_all!1 = "vol"
a_all!0 = "model"
n3_all!1 = "vol"
n3_all!0 = "model"
ad_all!1 = "vol"
ad_all!0 = "model"
printVarSummary(n3_all)

a_all_v = dim_avg_n_Wrap(a_all, 1)
n3_all_v = dim_avg_n_Wrap(n3_all, 1)
ad_all_v = dim_avg_n_Wrap(ad_all, 1)
printVarSummary(a_all_v)
printVarSummary(n3_all_v)
;--------plot-------------------------
wks = gsn_open_wks("ps", diri+"cmip.all.nino34.box.filter")
plot = new(2,graphic)
res = True
;res@gsnFrame = False
;res@gsnDraw  = False
res@vpWidthF    = 0.7        ; Make plots wider
res@vpHeightF   = 0.3
;res@vpXF           = 0.15            ; start plot at x ndc coord 

;res@xyMarkLineMode     = "Markers"  ; choose which have markers
;; res@xyDashPatterns      = 1                      ; solid line 
;  res@xyMarker         =  16                      ; choose type of marker  
;  res@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
;   res@xyLineThicknesses   = (/1,2/)                ; set second line to 2
;res@xyMarkerOpacityF = 0.7;0.5;
res@trYMinF            = -3.0
res@trYMaxF            =  3.0
;res@trXMinF = 1        
;res@trXMaxF = nm  

;res@gsnYRefLine       = 0.  ;Xrefline和Yrefline不能共存？！
;res@gsnYRefLineColor = "black"
;res@gsnYRefLineDashPattern = 0;0
;res@xyLineThicknessF = 2.
;res@xyMonoMarkerColor  = False ;True ;          ; Allow different colors for markers
res@tiXAxisString = "model"
res@tiYAxisString = "nino3.4(DJF) year(0)"
res@tiXAxisFontHeightF = 0.015
res@tiYAxisFontHeightF = 0.015

colors = (/"lightsalmon1","LightBlue","lightgoldenrod3","limegreen","orchid3", "lightslategray",\
 "brown","blue2","red2","orangered","goldenrod1","hotpink","greenyellow","lemonchiffon"/);colors
;gsn_define_colormap(wks,colors)
;colors = read_colormap_file("cb_9step")
;res@xyMarkerColors = colors

x_axis = ispan(0 , nm-1, 1)
res@tmXBMode   = "Explicit"                                                
res@tmXBValues = x_axis-1
res@tmXBLabels = name_model
res@tmXBLabelAngleF = 45.
res@tmXBLabelFontHeightF = 0.007

Y1 = n3_all(vol|:,model|:,time|6)
Y2 = n3_all(vol|:,model|:,time|7)
;printVarSummary(x_axis)
;printVarSummary(Y1)
;  res@xyMarkerColor  =  "Background"       ;也就是透明，作为背景？！！              
;res@xyMarkerOpacityF = 0.0;0.5;
;    do iv=0,nv-1
;
;plot(0) = gsn_csm_xy(wks, (/x_axis/), (/Y1/), res)
;res@tiYAxisString = "nino3.4(DJF) year(2)"
;plot(1) = gsn_csm_xy(wks, (/x_axis/), (/Y2/), res)

;-------------add box ------------------------------------
yval1 = new((/nm,5/),float,-999.)
yval2 = new((/nm,5/),float,-999.)

  opt = True
  opt@PrintStat = False
do im=0,nm-1
  statb1 = stat_dispersion(Y1(:,im), opt );Computes a number of robust statistics.什么统计量都有！平均中位标准差等
  statb2 = stat_dispersion(Y2(:,im), opt )
yval1(im,0) = statb1(8)-statb1(1) ;median - std?
yval1(im,1) = statb1(8)-statb1(1)/2;statb1(6)
yval1(im,2) = statb1(0) ;mean
yval1(im,3) = statb1(8)+statb1(1)/2;statb1(10)
yval1(im,4) = statb1(8)+statb1(1)

yval2(im,0) = statb2(8)-statb2(1) ;median - std?
yval2(im,1) = statb2(8)-statb2(1)/2;statb2(6)
yval2(im,2) = statb2(0) ;mean
yval2(im,3) = statb2(8)+statb2(1)/2;statb2(10)
yval2(im,4) = statb2(8)+statb2(1)
end do 
print(statb1)
;**********************************************
; resources for polylines that draws the boxes
;**********************************************  
  res2 = True
  res2@gsnFrame = False
  res2@gsnDraw  = False
  llres                   = True      
  llres@gsLineThicknessF  = 2.5                 ; line thickness 
;**********************************************
; resources that control color and width of boxes
;**********************************************  
  opti          = True      
  opti@boxWidth = .25       ; Width of box (x units)
  opti@boxColors = (/"grey51"/)   ; Color of box(es)
;***********************************************
printVarSummary(x_axis)
printVarSummary(yval1)
  plot(0) = boxplot(wks,x_axis,yval1,opti,res,llres) ; All 3 options used...
res@tiYAxisString = "nino3.4(DJF) year(1)"
  plot(1) = boxplot(wks,x_axis,yval2,opti,res,llres)
;overlay(plot(0), plot00)
;overlay(plot(1), plot11)
;------add polymarkers---------------------------------
  gsres     = True
  gsres@xyMarkerOpacityF = 1.0;0.5;

  do im=0,nm-1 
    do iv=0,nv-1
    xmarker = x_axis(im)
    ymarker1 = Y1(iv,im)
    ymarker2 = Y2(iv,im)
    
    gsres@gsMarkerColor = colors(iv)
    gsres@gsMarkerIndex = 16
   ;gsres@gsMarkerSizeF = 15.0
                        ; add (attach) polymarkers to existing plot object 
    plot@$unique_string("dum")$ = gsn_add_polymarker(wks,plot(0),xmarker,ymarker1,gsres)
    plot@$unique_string("dum")$ = gsn_add_polymarker(wks,plot(1),xmarker,ymarker2,gsres)
  end do
end do 

;==============================================
;  attach plots
;=============================================
;res = True
;res@txString = "BJ index heat budget 1982 (unit:0.01K/mon)"
;res@gsnPanelLabelBar   = True                 ; add common colorbar
;res@gsnMaximize = True
;res@gsnAttachBorderOn = True ; False ;图之间用线间隔
;
;res@gsnAttachPlotsXAxis = True ;False ; 是否以x轴为连接线
;
; newplot = gsn_attach_plots(plot(0),(/plot(1),plot(2)/),res,res)

; Set up resources for a customized legend.
;
  lgres                    = True
  lgres@lgMarkerColors       = colors
  lgres@lgBoxLineDashPattern = 1;(/1,0,0,1,1/)
  lgres@lgMarkerSizeF  = 0.01
  lgres@lgItemType         = "Markers";"Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .06            ; legend label font thickness
  lgres@vpWidthF           = 0.05           ; width of legend (NDC)
  lgres@vpHeightF          = 0.08           ; height of legend (NDC)
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter
  lgres@lgMonoMarkerIndex = True
  lgres@lgMarkerIndex = 16

  labels = volyear_c;name_model(sort_name) ;
;TP
  nitems = dimsizes(labels)                 ; number of legend items

; Create legend
  lbid = gsn_create_legend(wks,nitems,labels,lgres)
  ;
; Use gsn_add_annotation to attach this legend to our existing plot.
; This way, if we resize the plot, the legend will stay with the
; plot and be resized automatically.
;
  amres                  = True

;
; Point (0,0) is the dead center of the plot. Point (0,.5) is center,
; flush bottom. Point (0.5,0.5) is flush bottom, flush right.
;
 amres@amOrthogonalPosF = -0.28   ; This is above the top edge of the plot.
  amres@amJust           = "BottomRight"    ; Use bottom right corner of box
                                            ; for determining its location.
  amres@amParallelPosF   = 0.5              ; Move legend to right
  ;amres@amOrthogonalPosF = 0.5              ; Move legend down.

  annoid = gsn_add_annotation(plot(0),lbid,amres)  ; add legend to plot

 resk = True
resk@gsnMaximize = True
resk@gsnPaperOrientation = "portrait"   ; force portrait
resk@txString = "volcano nino3.4(DJF) filtered ";+volyear(iv) ;+name_model(ik);A"; with Elnino(grey)
gsn_panel(wks, plot, (/2,1/), resk)
;end do 
end