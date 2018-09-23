load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
; Nino 1+2 (0-10S, 90W-80W), Nino 3 (5N-5S, 150W-90W), 
; Nino 3.4/ONI (5N-5S, 170W-120W) and Nino 4 (5N-5S, 160E-150W)
;----------------------------------------------------------------------
; local procedures to draw 'error bars' in the x and y direction
;----------------------------------------------------------------------
undef("add_x_errorbar")
procedure add_x_errorbar(wks, plot, x[*]:numeric,y[*]:numeric, xerror[*]:numeric, res)
local dimz, xpoints, ypoints, gsres, i
begin
  gsres = res       
  dimz  = dimsizes(y)

  xpoints = (/x-xerror, x+xerror/)
  ypoints = (/ y, y/)

  do i=0, dimz(0)-1
     plot@$tostring("x_errorbar"+i)$ = gsn_add_polyline(wks, plot, xpoints(:,i), ypoints(:,i), gsres )
  end do
end
;---
undef("add_y_errorbar")
procedure add_y_errorbar(wks, plot, x[*]:numeric,y[*]:numeric, yerror[*]:numeric, res)
local dimz, xpoints, ypoints, gsres, i
begin
  gsres = res
  dimz  = dimsizes(y)

  xpoints = (/x, x/)
  ypoints = (/ y-yerror, y+yerror/)
  do i=0, dimz(0)-1
     plot@$tostring("y_errorbar"+i)$ = gsn_add_polyline(wks, plot, xpoints(:,i), ypoints(:,i), gsres )
  end do
end
;============================================
;             main 
;==============================================
begin
name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\;
            "CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM",\;
            "CNRM-CM5-2",  "CNRM-CM5","CSIRO-Mk3-6-0",  \
             "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
             "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES",\
             "MIROC5","MIROC-ESM-CHEM",\
              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);;"HadGEM2-CC"(没有tas资料）,"FIO-ESM","FGOALS-g2",  "MRI-CGCM3",,"GISS-E2-H"

Ammann = (/"bcc-csm1-1-m","bcc-csm1-1","CCSM4","CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM",\
          "CESM1-WACCM","CNRM-CM5-2",  "CNRM-CM5","NorESM1-M","NorESM1-ME"/)
Sato = (/"ACCESS1-0","ACCESS1-3","CanESM2","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H-CC", "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES",\
             "MIROC5","MIROC-ESM-CHEM" /)
others = (/ "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1"/) 

nm = dimsizes(name_model)
nam = dimsizes(Ammann)
nsa = dimsizes(Sato)
not = dimsizes(others)
EL=(/1876,1880,1884,1888,1895,1896,1899,1902,1904,1905,1911,1914,1918,1923,1925,1930,1939\
,1940,1941,1951,1953,1957,1963,1965,1968,1969,1972,1976,1977,1982,1986,1987,1991,1994,1997\
/);  没除掉vol年,2002,2004,2006,2009

ne = dimsizes(EL)
diri = "/work/cmip_vol/"

;volmon = (/188308,190210,191301,196303,198203,199106/);
volyear = (/1883,1902,1963,1982,1991/);1886,1968, 1896(sh) (1968),1913,
;vol_tp_mag_sato = (/0.1860,0.0931,0.1406,0.1338,0.1627/);0.0692,0.0491, sato 0.0557 (0.0491)

volyear_c = (/"1883", "1902", "1963", "1982","1991"/);"1913",

nv = dimsizes(volyear)

volmo = volyear ;(/1,1,1,1,1,1/);(/5,10,1,3,3,6/);
volmo = 1
volmo = volmo-1

ymstrt = 1870
ymlast2 = 2004
ny2 = ymlast2-ymstrt+1
print(ny2)
nino3 =  asciiread(diri+"CMIP5.nino34.unnorm.1870-2004.txt", (/nm*ny2*12/), "float")
;nino3 =  asciiread(diri+"CMIP5.nino3.anomaly.unnorm.1870-2004.txt", (/nm*ny2*12/), "float")
;nino3 =  asciiread(diri+"CMIP5.nino34.anomaly.unnorm.1870-2004.thetao.txt", (/nm,ny2*12/), "float")
;tt_a = asciiread(diri+"CMIP5.thetao100-150.nino34.anomaly.unnorm.1870-2004.23.txt", (/nm,ny2*12/), "float");wait for shortwave!
;tt_a = asciiread(diri+"CMIP5.i20.nino34.anomaly.unnorm.1870-2004.23.txt", (/nm,ny2*12/), "float");wait for shortwave!
;tt_a = asciiread(diri+"CMIP5.nino34.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");wait for shortwave!
;air_a_c = asciiread(diri+"CMIP5.part.top.sw.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");wait for shortwave!

;gmst_obs = asciiread(diri+"CMIP5.GMSST.org.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")
;gmst_a = asciiread(diri+"CMIP5.GMSST.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")

;air_a_c =  asciiread(diri+"CMIP5.GMST.anomaly.unnorm.1870-2004.thetao.txt", (/nm,ny2*12/), "float");借用
air_a_c =  asciiread(diri+"CMIP5.GMST.anomaly.unnorm.1870-2004.txt", (/nm*ny2*12/), "float");借用

gmst_a = asciiread(diri+"CMIP5.TP(-10-10.160-270).sst.anomaly.unnorm.1870-2004.txt", (/nm*ny2*12/), "float")
tt_a =  asciiread(diri+"CMIP5.air.MC.anomaly.unnorm.1870-2004.txt", (/nm*ny2*12/), "float");借用

;gmst_a = asciiread(diri+"CMIP5.air.EU.anomaly.unnorm.1870-2004.txt", (/nm*ny2*12/), "float")
;tt_a =  asciiread(diri+"CMIP5.air.NA.anomaly.unnorm.1870-2004.txt", (/nm*ny2*12/), "float");借用

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
g_all = new((/nm,nv,nyr/), float)
ad_all = new((/nm,nv,2/), float)
do nn=0,nm-1;1;
nino3_om =nino3_a_obs(nn,:)
tt_om = tt_a_obs(nn,:)
g_om =gmst_a_obs(nn,:)
a_om =air_a(nn,:)  
nino3_a1 = nino3_a_obs(nn,:)
air_a1 = air_a(nn,:)
tt_a1 = tt_a_obs(nn,:)
g_a1 = gmst_a_obs(nn,:)
;----------------filter----------------------PAY ATTENTION! This is month data, not years!
ihp   = 1       ; low-pass filter: ihp = 0; high-pass ihp = 1; band-pass ihp = 2.
sigma = 1.0                           ; Lanczos sigma
nWgt  = 51 ;必须是奇数 表示权重系数的总个数，越大效果越好，但是两端损失越多
fca   = 1./(10*12) ;                        
fcb   = -999;
wgt  = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
nino3_a1 = wgt_runave_n ( nino3_a_obs(nn,:), wgt, 0, 0 ) ;最后一个是指sst1第0维
tt_a1 = wgt_runave_n ( tt_a_obs(nn,:), wgt, 0, 0 ) ;最后一个是指sst1第0维
air_a1 = wgt_runave_n ( air_a(nn,:), wgt, 0, 0 ) ;最后一个是指sst1第0维
g_a1 = wgt_runave_n ( gmst_a_obs(nn,:), wgt, 0, 0 )
;-----------moving average---------------------
nrun = 3
 endopt    = 1 ;起始结束不损失滑动的月份

  nino3_om = runave_n_Wrap (nino3_a1, nrun, endopt, 0)
  t_om =    runave_n_Wrap (tt_a1, nrun, endopt, 0)
  g_om =     runave_n_Wrap (g_a1, nrun, endopt, 0)
  a_om =runave_n_Wrap (air_a1, nrun, endopt, 0)
;printVarSummary(tt_om)


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
g_v_a!0 = "vol"
g_v_a!1 = "time"
g_djf = new((/nv,nyr/),float)
nobs_v_a!0 = "vol"
nobs_v_a!1 = "time"
nobs_djf = new((/nv,nyr/),float)
aobs_v_a!0 = "vol"
aobs_v_a!1 = "time"
aobs_djf = new((/nv,nyr/),float)
tobs_v_a!0 = "vol"
tobs_v_a!1 = "time"
tobs_djf = new((/nv,nyr/),float)
;gobs_v_a!0 = "vol"
;gobs_v_a!1 = "time"
;gobs_djf = new((/nv,nyr/),float)
do kk=0,nv-1
n3_djf(kk,:) = month_to_season(n3_v_a(kk,:), "DJF")
a_djf(kk,:) = month_to_season(a_v_a(kk,:), "DJF")
t_djf(kk,:) = month_to_season(t_v_a(kk,:), "DJF")
g_djf(kk,:) = month_to_season(g_v_a(kk,:), "DJF")
nobs_djf(kk,:) = month_to_season(nobs_v_a(kk,:), "DJF")
aobs_djf(kk,:) = month_to_season(aobs_v_a(kk,:), "DJF")
tobs_djf(kk,:) = month_to_season(tobs_v_a(kk,:), "DJF")
end do 
a_all(nn,:,:)  = a_djf
n3_all(nn,:,:) = n3_djf
t_all(nn,:,:) = t_djf
g_all(nn,:,:) = g_djf
ad_all(nn,:,0) = a_v_a(:,6*12)-a_v_a(:,5*12)
ad_all(nn,:,1) = a_v_a(:,7*12)-a_v_a(:,6*12)
end do 
a_all!1 = "vol"
a_all!0 = "model"
n3_all!1 = "vol"
n3_all!0 = "model"
t_all!1 = "vol"
t_all!0 = "model"
g_all!1 = "vol"
g_all!0 = "model"
ad_all!1 = "vol"
ad_all!0 = "model"
printVarSummary(n3_all)

a_all_v = dim_avg_n_Wrap(a_all, 1)
n3_all_v = dim_avg_n_Wrap(n3_all, 1)
a_all_std = dim_stddev_n_Wrap(a_all, 1)
n3_all_std = dim_stddev_n_Wrap(n3_all, 1)
a_all_max = dim_max_n_Wrap(a_all, 1)
n3_all_max = dim_max_n_Wrap(n3_all, 1)

printVarSummary(a_all_v)
printVarSummary(n3_all_std)

;---------forcing sort--------------
n_am = new((/nam,nv,nyr/),float)
n_sa = new((/nsa,nv,nyr/),float)
n_ot = new((/not,nv,nyr/),float)
a_am = new((/nam,nv,nyr/),float)
a_sa = new((/nsa,nv,nyr/),float)
a_ot = new((/not,nv,nyr/),float)
;ia = 0
;is = 0
;io = 0
do i = 0,nm-1
  do ia = 0,nam-1
  if (name_model(i).eq.Ammann(ia)) then
  n_am(ia,:,:) = n3_all(i,:,:)
  a_am(ia,:,:) = a_all(i,:,:)
  ;ia = ia+1
  end if
end do 
do is = 0,nsa-1
  if (name_model(i).eq.Sato(is)) then
  n_sa(is,:,:) = n3_all(i,:,:)
  a_sa(is,:,:) = a_all(i,:,:)
  ;is = is+1
  end if
end do 
do io = 0,not-1
  if (name_model(i).eq.others(io)) then
  n_ot(io,:,:) = n3_all(i,:,:)
  a_ot(io,:,:) = a_all(i,:,:)
  ;io = io+1
  end if
end do 
end do 
printVarSummary(n_ot)
;--------plot-------------------------
wks = gsn_open_wks("ps", diri+"cmip.vol.nino34&gmst.box.forcing");+name_model(nn));.A2");&el
plot = new(1,graphic)

res = True
res@gsnFrame = False
res@gsnDraw  = False
 
res@xyMarkLineMode     = "Markers"  ; choose which have markers
; res@xyDashPatterns      = 1                      ; solid line 
  res@xyMarker         =  1                      ; choose type of marker  
 ; res@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
;   res@xyLineThicknesses   = (/1,2/)                ; set second line to 2

res@xyMarkerOpacityF = 0.5;1.0;0.5;
  res@tmYROn                 = False    ; right off
res2 = res
res3 = res
res4 = res
res@vpWidthF    = 0.5        ; Make plots wider
res@vpHeightF   = 0.5
res@tmXTOn              = False   ; Turn off top tickmarks
res@tmYLMinorOn         = False   ; Turn off left minor tickmarks
res@tmXMajorGrid                = True
res@tmYMajorGrid                = True
res@tmXMajorGridThicknessF      = 1.0
res@tmYMajorGridThicknessF      = 1.0
res@tmXMajorGridLineDashPattern = 2
res@tmYMajorGridLineDashPattern = 2

res@gsnRightString = ""
res@gsnLeftString = ""
res@gsnLeftStringFontHeightF = 0.03
;res@tiXAxisString   = " "
res@trXMinF            = -3.
res@trXMaxF            =  3.
res@trYMinF            = -0.5;-2.0;
res@trYMaxF            =  0.5;2.0;
res@gsnXRefLine       = 0;volyear;(/1980,1988/);1900;
res@gsnXRefLineThicknessF = 2.
res@gsnXRefLineColor = "black";(/"grey","grey","grey","grey","grey","grey"/);
res@gsnXRefLineDashPattern = 0;0

res@gsnYRefLine       = 0.  ;Xrefline和Yrefline不能共存？！
res@gsnYRefLineColor = "black"
res@gsnYRefLineDashPattern = 0;0
res@xyLineThicknessF = 2.

res2 = res
res3 = res
;colors = (/"lightsalmon1","LightBlue","lightgoldenrod3","limegreen","orchid3","lightslategray"/)
;gsn_define_colormap(wks2,colors)
;res@xyMarkerColors = (/"lightsalmon1","LightBlue","lightgoldenrod3","limegreen","orchid3",\
;                      "lightslategray","brown","red2","skyblue","black"/);colors
;res@xyMarkerSizes = vol_mag*0.003 
;colors = (/"lightsalmon1","LightBlue","lightgoldenrod3","limegreen","orchid3", "lightslategray",\
; "brown","blue2","red2","orangered","goldenrod1","hotpink","greenyellow","lemonchiffon"/);colors
colors = read_colormap_file("t2m_29lev") ;("cosam");("BkBlAqGrYeOrReViWh200");("default");("cosam");"GMT_paired"

res@gsnLeftString = "";"(a)";"

 resk = True
resk@gsnMaximize = True
resk@gsnPaperOrientation = "portrait"   ; force portrait

res@xyMarkerColor = "black" ;"grey70";"blue" ;colors(::3,:);(0);(is)
;res@xyMarkerColors = colors;"grey70";(::2,:);
res@xyMarker         = 1
res@xyMarkerSizeF    = 0.03
;resk@txFontHeightF = 0.08
resk@txString = "different forcing in CMIP5";"Ammann:red Sato:green Others:blue";+name_model(ik);A"; with Elnino(grey)
res@tiXAxisString = "nino 3.4 in year(0)"
res@tiYAxisString = "GMST in year(0)"
res@trXMinF             = -3.0 ;-0.5;
res@trXMaxF             =  3.0 ;0.5;
res@trYMinF             = -0.5;-0.3;
res@trYMaxF             =  0.5;0.3;
;res3@trXMinF            = -3.0;-0.5;
;res3@trXMaxF            =  3.0; 0.5;
;res3@trYMinF            = -3.0;-0.5;
;res3@trYMaxF            =  3.0;0.5; 
x1 = ndtooned(n_am(:,:,6));n3_all_v(:,5)
y1 = ndtooned(a_am(:,:,6));n3_all_v(:,6)
ux1 = dim_stddev_n_Wrap(x1, 0)
uy1 = dim_stddev_n_Wrap(y1, 0)
x2 = ndtooned(n_sa(:,:,6));n3_all_v(:,5)
y2 = ndtooned(a_sa(:,:,6));n3_all_v(:,6)
ux2 = dim_stddev_n_Wrap(x2, 0)
uy2 = dim_stddev_n_Wrap(y2, 0)
x3 = ndtooned(n_ot(:,:,6));n3_all_v(:,5)
y3 = ndtooned(a_ot(:,:,6));n3_all_v(:,6)
ux3 = dim_stddev_n_Wrap(x3, 0)
uy3 = dim_stddev_n_Wrap(y3, 0)
;u1 = n3_all_std(:,5);n3_all_max(:,6)
;u2 = n3_all_std(:,6);a_all_max(:,6)
;u1 = n3_all_max(:,6)-x
;u2 = a_all_max(:,6)-y
x0 = ndtooned(n3_all(:,:,6))
y0 = ndtooned(a_all(:,:,6))

res@xyMonoMarkerColor  = True ; False ;         ; Allow different colors for markers
plot(0) = gsn_csm_xy(wks, transpose((/x0/)),transpose((/y0/)),res) 
res3@xyMonoMarkerColor  = False ;True ; 
res3@xyMarkerOpacityF = 1.0;
res3@xyMarkerSizeF    = 0.01
res3@xyMarker         =  16
res3@xyMarkerColor = "red" ;colors
plot1 = gsn_csm_xy(wks, transpose((/x1/)),transpose((/y1/)),res3) 
res3@xyMarkerColor = "green" ;colors
plot2 = gsn_csm_xy(wks, transpose((/x2/)),transpose((/y2/)),res3) 
res3@xyMarkerColor = "blue" ;colors
plot3 = gsn_csm_xy(wks, transpose((/x3/)),transpose((/y3/)),res3) 
;res3@xyMarkerColor = "blue" 
;plot1   = gsn_csm_xy(wks, transpose((/none01/)),transpose((/none11/)),res3) 
overlay(plot(0), plot1)
overlay(plot(0), plot2)
overlay(plot(0), plot3)

;-----------error bar----------------------------------
;   gsresx = True                                          ; add x error bars
;   gsresx@gsLineColor      = "lightskyblue" ;"grey70";"";"red"
; ;  gsresx@gsLineColor      = colors
;   gsresx@gsLineThicknessF = 2.0
;   add_x_errorbar(wks, plot(0), x1, y1, ux1, gsresx)
;   add_x_errorbar(wks, plot(0), x2, y2, ux2, gsresx)
;   add_x_errorbar(wks, plot(0), x3, y3, ux3, gsresx)
;    
;   gsresy = True
;   gsresy@gsLineColor      = "lightskyblue" ;"grey70";                      ; add y error bars
;   ;gsresx@gsLineColor      = colors
;   gsresy@gsLineThicknessF = 2.0
;   add_y_errorbar(wks, plot(0), x1, y1, uy1, gsresy)
;   add_y_errorbar(wks, plot(0), x2, y2, uy2, gsresy)
;   add_y_errorbar(wks, plot(0), x3, y3, uy3, gsresy)

;-------------------------------------------------------
; Set up resources for a customized legend.
;
  lgres                    = True
  lgres@lgMarkerColors       = (/"red","green","blue"/) ;colors;(::2,:);(/"red","green" /) ; "blue"
  lgres@lgBoxLineDashPattern = 1;(/1,0,0,1,1/)
  lgres@lgMarkerSizeF  = 0.01
  lgres@lgItemType         = "Markers";"Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .06            ; legend label font thickness
  lgres@vpWidthF           = 0.1           ; width of legend (NDC)
  lgres@vpHeightF          = 0.2           ; height of legend (NDC)
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter
  lgres@lgMonoMarkerIndex = True
  lgres@lgMarkerIndex = 16

  labels = (/"Ammann","Sato","others" /) ;name_model ;volyear_c;(/"A","B"/) ;name_model(sort_name) ;,"C"
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
  amres@amJust           = "BottomRight"    ; Use bottom right corner of box
                                            ; for determining its location.
  amres@amParallelPosF   = 0.5              ; Move legend to right
  amres@amOrthogonalPosF = 0.5              ; Move legend down.

  annoid = gsn_add_annotation(plot(0),lbid,amres)  ; add legend to plot

delete(res)
delete(res2)
delete(res3)
delete(res4)
;print(nn)

gsn_panel(wks, plot, (/1,1/), resk)
end