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
;         "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/); thetao
;
;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\
;            "CNRM-CM5","CNRM-CM5-2","CSIRO-Mk3-6-0", \;  
;             "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC",\;
;             "HadGEM2-ES",\;"HadGEM2-AO","GISS-E2-H","HadCM3","HadGEM2-CC","MIROC-ESM","MRI-CGCM3",
;              "MIROC5","MIROC-ESM-CHEM",\
;              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P",\
;              "MRI-ESM1","NorESM1-M","NorESM1-ME"/); rsdscs

nm = dimsizes(name_model)
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

;gmst_a = asciiread(diri+"CMIP5.TP(-10-10.160-270).sst.anomaly.unnorm.1870-2004.txt", (/nm*ny2*12/), "float")
;tt_a =  asciiread(diri+"CMIP5.air.MC.anomaly.unnorm.1870-2004.txt", (/nm*ny2*12/), "float");借用

;nino3 =  asciiread(diri+"CMIP5.i20.nino34.anomaly.unnorm.1870-2004.txt", (/nm*ny2*12/), "float")

gmst_a = asciiread(diri+"CMIP5.thetao100-150.160-270.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")

tt_a =  asciiread(diri+"CMIP5.thetao0-300.120-280.anomaly.unnorm.1870-2004.txt", (/nm*ny2*12/), "float");借用

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
nyear1 = 3
nyear2 = 5
nyr = nyear1+nyear2
nmon = nyear1*12+nyear2*12
n3_all = new((/nm,nv,nyr/), float)
a_all = new((/nm,nv,nyr/), float)
t_all = new((/nm,nv,nyr/), float)
g_all = new((/nm,nv,nyr/), float)
n3_all_mon = new((/nm,nv,nmon/), float)
a_all_mon = new((/nm,nv,nmon/), float)
t_all_mon = new((/nm,nv,nmon/), float)
g_all_mon = new((/nm,nv,nmon/), float)
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
nWgt  = 71 ;必须是奇数 表示权重系数的总个数，越大效果越好，但是两端损失越多
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
;ad_all(nn,:,0) = a_v_a(:,(nyear1+1)*12)-a_v_a(:,nyear1*12)
;ad_all(nn,:,1) = a_v_a(:,(nyear1+2)*12)-a_v_a(:,(nyear1+1)*12)
a_all_mon(nn,:,:)  = a_v_a
n3_all_mon(nn,:,:) = n3_v_a
t_all_mon(nn,:,:) = t_v_a
g_all_mon(nn,:,:) = g_v_a
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
a_all_mon!1 = "vol"
a_all_mon!0 = "model"
n3_all_mon!1 = "vol"
n3_all_mon!0 = "model"
t_all_mon!1 = "vol"
t_all_mon!0 = "model"
g_all_mon!1 = "vol"
g_all_mon!0 = "model"
printVarSummary(n3_all)

;-------center finite difference -------------------
r = 3.
r_n = center_finite_diff_n(n3_all_mon, r, False, 0, 2)
r_a = center_finite_diff_n(a_all_mon, r, False, 0, 2)
r_g = center_finite_diff_n(g_all_mon, r, False, 0, 2)
r_t = center_finite_diff_n(t_all_mon, r, False, 0, 2)
printVarSummary(r_n)
;----二阶求导？
r_n = center_finite_diff_n(r_n, r, False, 0, 2)
r_a = center_finite_diff_n(r_a, r, False, 0, 2)
r_g = center_finite_diff_n(r_g, r, False, 0, 2)
r_t = center_finite_diff_n(r_t, r, False, 0, 2)
;r_n = n3_all_mon
;r_a = a_all_mon
;r_g = g_all_mon
;r_t = t_all_mon
;
;-------sort----------------------------
;nan = n3_all(:,:,6)
;naa = a_all(:,:,6)
;k = 0
;l = 0
;j = 0
;nlam = new(68,float);69  thetao:54
;nlav = new(68,float)
;;nelm = new(64,float)
;;nelv = new(64,float)
;nwelm = new(36,float);35 thetao:28
;nwelv = new(36,float)
;nselm = new(32,float);29 thetao:23
;nselv = new(32,float)
;nla = nan
;nla = 0
;;nel = nan
;;nel = 0
;nwel = nan
;nwel = 0
;nsel = nan
;nsel = 0
;do im = 0,nm-1
;  do iv=0,nv-1
;    if (nan(im,iv).lt.0.0.and.naa(im,iv).le.0.0) then
;      nla(im,iv) = nan(im,iv)
;      nlam(k) = im
;      nlav(k) = iv
;      k = k+1
;    else if(nan(im,iv).ge.0.0.and.naa(im,iv).le.0.0)then 
;      nwel(im,iv) = nan(im,iv)
;      nwelm(l) = im
;      nwelv(l) = iv 
;      l = l+1
;    else if(nan(im,iv).ge.0.0.and.naa(im,iv).gt.0.0)then
;      nsel(im,iv) = nan(im,iv)
;      nselm(j) = im
;      nselv(j) = iv   
;      j = j+1
;    end if
;  end if 
;  ;else 
;  ; nel(im,iv) = nan(im,iv)
;  ; nelm(l) = im
;  ; nelv(l) = iv 
;  ; l = l+1  
;end if 
;  end do 
;end do 
;print(k) ;81-7 = 74
;print(l) ;30
;print(j) ;34pooooooooooo
;;printVarSummary(data)
;;print(nsel)
;printVarSummary(nla)
;  N    = nm    ; number of rows
;  M    = nv     ; number of columns (N,M)
;fWidth = 10    ; specify the format width
;  fDec   = 4     ; specify the number to the right of decimal point
;  fmtx   = M + "f" + fWidth + "." + fDec              ; fmtx="50f10.3"
;  opt = True
;  opt@fout = diri+"cmip.nino.sort.cla.filter.txt"
;  write_matrix (nla, fmtx, opt)
; ; opt@fout = diri+"cmip.nino.sort.el.filter.23.txt"
; ; write_matrix (nel, fmtx, opt)
;  opt@fout = diri+"cmip.nino.sort.cel.filter.txt"
;  write_matrix (nwel, fmtx, opt)
;  opt@fout = diri+"cmip.nino.sort.wel.filter.txt"
;  write_matrix (nsel, fmtx, opt)
;;printVarSummary(data)
;--------------select----------------------------------
;;--------------sort----------------------
nla = readAsciiTable(diri+"cmip.nino.sort.cla.filter.txt", 5, "float", 0)
nwel = readAsciiTable(diri+"cmip.nino.sort.cel.filter.txt", 5, "float", 0)
nsel = readAsciiTable(diri+"cmip.nino.sort.wel.filter.txt", 5, "float", 0)
nla!0 = "model"
nla!1 = "vol"
nwel!0 = "model"
nwel!1 = "vol"
nsel!0 = "model"
nsel!1 = "vol"
printVarSummary(nwel)

t1 = r_t
t2 = r_t
t3 = r_t
a1 = r_a
a2 = r_a
a3 = r_a
g1 = r_g
g2 = r_g
g3 = r_g
n1 = r_n
n2 = r_n
n3 = r_n
do it = 0,nmon-1
  t1(:,:,it) = nla
  t2(:,:,it) = nwel
  t3(:,:,it) = nsel
  g1(:,:,it) = nla
  g2(:,:,it) = nwel
  g3(:,:,it) = nsel
  a1(:,:,it) = nla
  a2(:,:,it) = nwel
  a3(:,:,it) = nsel
  n1(:,:,it) = nla
  n2(:,:,it) = nwel
  n3(:,:,it) = nsel
end do 
  t1 = where(t1.ne.0.0 , r_t, t_all@_FillValue)
  t2 = where(t2.ne.0.0 , r_t, t_all@_FillValue)
  t3 = where(t3.ne.0.0 , r_t, t_all@_FillValue)

  g1 = where(g1.ne.0.0 , r_g, t_all@_FillValue)
  g2 = where(g2.ne.0.0 , r_g, t_all@_FillValue)
  g3 = where(g3.ne.0.0 , r_g, t_all@_FillValue)

  a1 = where(a1.ne.0.0 , r_a, t_all@_FillValue)
  a2 = where(a2.ne.0.0 , r_a, t_all@_FillValue)
  a3 = where(a3.ne.0.0 , r_a, t_all@_FillValue)

  n1 = where(n1.ne.0.0 , r_n, t_all@_FillValue)
  n2 = where(n2.ne.0.0 , r_n, t_all@_FillValue)
  n3 = where(n3.ne.0.0 , r_n, t_all@_FillValue)

printVarSummary(a2)
a1_avg = dim_avg_n_Wrap(a1, (/0,1/))
a2_avg = dim_avg_n_Wrap(a2, (/0,1/))
a3_avg = dim_avg_n_Wrap(a3, (/0,1/))

n1_avg = dim_avg_n_Wrap(n1, (/0,1/))
n2_avg = dim_avg_n_Wrap(n2, (/0,1/))
n3_avg = dim_avg_n_Wrap(n3, (/0,1/))

g1_avg = dim_avg_n_Wrap(g1, (/0,1/))
g2_avg = dim_avg_n_Wrap(g2, (/0,1/))
g3_avg = dim_avg_n_Wrap(g3, (/0,1/))

t1_avg = dim_avg_n_Wrap(t1, (/0,1/))
t2_avg = dim_avg_n_Wrap(t2, (/0,1/))
t3_avg = dim_avg_n_Wrap(t3, (/0,1/))

;--------plot-------------------------
wks2 = gsn_open_wks("ps", diri+"cmip.gmst.sort.xy.trend2.C")
;wks2 = gsn_open_wks("ps", "/work/obs.eachvol.nino34&air.ERSST&HadISST")
nys1 = 1
nys2 = 2
x_axis = ispan(-nys1*12, nys2*12-1, 1)
x_axis!0 = "time"
;x_axis = fspan(ymstrt, ymlast2, ny2*12)
;do j=0,nv-1
plot = new(3,graphic)
;do i=0,nm-1

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@xyLineOpacityF = 0.5
  res@tmYROn                 = False    ; right off
res2 = res
res3 = res
res4 = res
res@vpWidthF    = 1.0        ; Make plots wider
res@vpHeightF   = 0.3
res@tmXTOn              = False   ; Turn off top tickmarks
res@tmYLMinorOn         = False   ; Turn off left minor tickmarks
res@tmYLLabelFontHeightF = 0.018
res@tmYLMajorLengthF = 0.01
;res@tmXMajorGrid                = True
;res@tmYMajorGrid                = True
;res@tmXMajorGridThicknessF      = 1.0
;res@tmYMajorGridThicknessF      = 1.0
;res@tmXMajorGridLineDashPattern = 2
;res@tmYMajorGridLineDashPattern = 2

res@gsnRightString = ""
res@gsnLeftString = ""
res@gsnLeftStringFontHeightF = 0.03
;res@tiXAxisString   = " "
res@trXMinF            = -nys1*12;ymstrt;1870;
res@trXMaxF            = nys2*12-1;ymlast2;2016;
res@tmXBMode   = "Explicit"                                                
res@tmXBValues = (/-36,-24,-12,0,11,23,35,47,59/) 
res@tmXBLabels = (/"-3","-2","year(-1)","year(0)","year(+1)","year(+2)","year(+3)","year(+4)","year(+5)"/)      
res@tiXAxisString = ""
res@tiXAxisFontHeightF =0.03
res@gsnXRefLine       = 0;volyear;(/1980,1988/);1900;
res@gsnXRefLineThicknessF = 2.
res@gsnXRefLineColor = "black";(/"grey","grey","grey","grey","grey","grey"/);
res@gsnXRefLineDashPattern = 16;0

res@gsnYRefLine       = 0.  ;Xrefline和Yrefline不能共存？！
res@gsnYRefLineColor = "black"
res@gsnYRefLineDashPattern = 16;0
res@gsnYRefLineThicknessF = 2.

res@xyMonoDashPattern      = True     ; force all solid lines

resk = True
resk@gsnMaximize = True
resk@gsnPaperOrientation = "portrait"   ; force portrait
resk@lbBoxEndCapStyle = "TriangleBothEnds"

res@xyLineThicknessF = 1.
res@xyLineOpacityF = 0.6
res@xyLineColor = "grey44"
;res@xyDashPattern = 14

res2@xyLineThicknessF = 5.;5.;
res2@xyLineOpacityF = 1.0
res2@xyLineColor = "black"
;res2@xyDashPattern = 14.
;res@xyLineThicknessF = 3.
;res@xyLineColors = (/"lightsalmon1","LightBlue","lightgoldenrod3","limegreen","orchid3","lightslategray"/)
;res@xyLineOpacityF = 0.7

res@gsnLeftStringOrthogonalPosF = -0.14


res3@xyLineThicknessF = 10.
res3@xyLineOpacityF = 1.0;0.5
res3@xyLineColor = "black"
aa1 = reshape(a1,(/nm*nv,nmon/))
aa2 = reshape(a2,(/nm*nv,nmon/))
aa3 = reshape(a3,(/nm*nv,nmon/))
nn1 = reshape(n1,(/nm*nv,nmon/))
nn2 = reshape(n2,(/nm*nv,nmon/))
nn3 = reshape(n3,(/nm*nv,nmon/))
gg1 = reshape(g1,(/nm*nv,nmon/))
gg2 = reshape(g2,(/nm*nv,nmon/))
gg3 = reshape(g3,(/nm*nv,nmon/))
tt1 = reshape(t1,(/nm*nv,nmon/))
tt2 = reshape(t2,(/nm*nv,nmon/))
tt3 = reshape(t3,(/nm*nv,nmon/))
printVarSummary(aa1)
printVarSummary(x_axis)
aa = 12*(nyear1-nys1)
bb = 12*(nyear1+nys2)-1
res@gsnLeftString = "(a) GMST [C]";(air temperature all level averaged)";
res@trYMinF             = -0.008;-0.031;-0.41;
res@trYMaxF             =  0.008; 0.031;0.41; 
;res3@trYMinF            = -2.;-0.5;
;res3@trYMaxF            =  2.;0.5;
plotxy0 = gsn_csm_xy(wks2,x_axis,aa3(:,aa:bb),res); a_e_a,,res3
plot1  = gsn_csm_xy(wks2,x_axis,a3_avg(aa:bb),res2)
;plot2  = gsn_csm_xy(wks2,x_axis,a_e_avg,res4)
;plot11  = gsn_csm_xy(wks2,x_axis,a90,res3)

res@gsnLeftString = "(b) Nino 3.4 index [C]";(air)"
res@trYMinF             =  -0.02;-0.21;-2.5;
res@trYMaxF             =   0.02; 0.21;2.5; 
;res3@trYMinF            = -1.;-0.5;
;res3@trYMaxF            =  1.;0.5;
plotxy1 = gsn_csm_xy(wks2,x_axis,nn3(:,aa:bb),res); ,g_e_a,res3
plot3  = gsn_csm_xy(wks2,x_axis,n3_avg(aa:bb),res2)
;plot4  = gsn_csm_xy(wks2,x_axis,t_e_avg,res4)
;plot33  = gsn_csm_xy(wks2,x_axis,t90,res3)

res@gsnLeftString = "(c) heat content (0-300m) [C]";(air)"
res@trYMinF            =  -0.01;-0.051;-0.5;
res@trYMaxF            =   0.01; 0.051;0.5;
;res3@trYMinF            = -1.;-0.5;
;res3@trYMaxF            =  1.;0.5;
plotxy2 = gsn_csm_xy(wks2,x_axis,tt3(:,aa:bb),res); ,g_e_a,res3
plot5  = gsn_csm_xy(wks2,x_axis,t3_avg(aa:bb),res2)
;plot6  = gsn_csm_xy(wks2,x_axis,g_e_avg,res4)
;plot55  = gsn_csm_xy(wks2,x_axis,g90,res3)


;res@gsnLeftString = "TNI(ERSST)"
;res@trYMinF            = -3.;-0.5;
;res@trYMaxF            = 3.;0.5;
;res3@trYMinF            = -3.;-0.5;
;res3@trYMaxF            = 3.;0.5;
;  plot(3) = gsn_csm_xy(wks2,x_axis,t_v_a,res) ; ,res3t_e_a,
;plot7  = gsn_csm_xy(wks2,x_axis,t_ov_avg,res2)  
;  plot8  = gsn_csm_xy(wks2,x_axis,t_e_avg,res4) 
 resk@txString = "";with Elnino(grey)

 overlay(plotxy0, plot1) 
; overlay(plot(0), plot2) 
 overlay(plotxy1, plot3)
; overlay(plot(1), plot4) 
 overlay(plotxy2, plot5)
; overlay(plot(2), plot6) 
;==============================================
;  attach plots
;=============================================
;res = True
;res@txString = "BJ index heat budget 1982 (unit:0.01K/mon)"
;res@gsnPanelLabelBar   = True                 ; add common colorbar
res1 = True
res2 = True
res2@gsnMaximize = True
res2@gsnPaperOrientation = "portrait"   ; force portrait
res2@gsnAttachBorderOn = True ; False ;图之间用线间隔
res2@gsnAttachPlotsXAxis = True ;False ; 是否以x轴为连接线
 newplot = gsn_attach_plots(plotxy0,(/plotxy1,plotxy2/),res2,res1)
;gsn_panel(wks,plot,(/6,1/),res)
draw(plotxy0)
frame(wks2)


delete(res)
delete(res2)
delete(res3)
delete(res4)
;print(nn)

;gsn_panel(wks2, plot, (/3,1/), resk)
;end do 
end
