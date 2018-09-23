load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
; Nino 1+2 (0-10S, 90W-80W), Nino 3 (5N-5S, 150W-90W), 
; Nino 3.4/ONI (5N-5S, 170W-120W) and Nino 4 (5N-5S, 160E-150W)
begin
name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\
       "CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM","CNRM-CM5-2",  \
      "CNRM-CM5","CSIRO-Mk3-6-0",  \
        "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
        "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES","MIROC5","MIROC-ESM-CHEM",\
         "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\;"bcc-csm1-1-m","bcc-csm1-1",
;"CMCC-CESM", "CMCC-CM","CMCC-CMS"   CMCC & FGOALS没有加入火山强迫！！！"CMCC-CESM", "CMCC-CM","CMCC-CMS",
;A = (/4,13,14,19,20,23,26,28,30,32/)
;B = (/2,3,7,10,15,17,18,25,29,33/)
;C = (/0,1,5,6,8,9,11,12,21,22,24,27,31/)
;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\
;            "CNRM-CM5-2","CSIRO-Mk3-6-0", \;  "CNRM-CM5",
;             "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\;
;             "GISS-E2-R","GISS-E2-R-CC","HadGEM2-ES", "MIROC5","MIROC-ESM-CHEM",\
;              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-ME"/);"NorESM1-M",
;A_n = (/"CanESM2","GFDL-CM3","GFDL-ESM2M","GISS-E2-R-CC",\
; "MIROC5", "MPI-ESM-LR", "MPI-ESM-P",/) ;正相关 主要限定在二、四象限"NorESM1-M",,"HadGEM2-AO"
;B_n = (/"CNRM-CM5-2",\
;  "GFDL-ESM2G","GISS-E2-H-CC",  "NorESM1-ME" , "GISS-E2-R",\
; "MPI-ESM-MR"/) ;负相关 主要限定在一、三象限"bcc-csm1-1-m","bcc-csm1-1", "CESM1-CAM5-1-FV2",
;C_n = (/"ACCESS1-0","ACCESS1-3","CCSM4","CNRM-CM5","CSIRO-Mk3-6-0""HadGEM2-ES",\
;  "MIROC-ESM-CHEM","MRI-ESM1" /) ;杂乱"CESM1-BGC"6,"CESM1-WACCM"9,"CESM1-FASTCHEM",
Y1 = (/1,2,3,7,10,11,12,14,16,17,18,22,23,25,26,28/)
 Y2 = (/0,4,5,6,8,9,13,15,19,20,21,24,27/)
A = (/4,19,22,21,23,25,27/);13,18,
B = (/2,3,10,14,16,28/);7,24,
C = (/0,1,5,8,12,20,21,26/);6,9,11,
na = dimsizes(A)
nb = dimsizes(B)
nc = dimsizes(C)

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
wks2 = gsn_open_wks("ps", diri+"cmip.vol(TP)sato.nino&period.scatter.sort.A");+name_model(nn));.A2");&el
plot = new(3,graphic)
sort_name = A
ns = na
;volmon = (/188308,190210,191301,196303,198203,199106/);
;volmon_NH = (/187503, 190703,191206, 193301,195603,198005/)
;volmon_SH = (/188606,193204,199108,201106/)
volyear_nh_gao = (/1912, 1925,1943,1976,1980/);NH 1875, 1907,1933,1956,
vol_nh_mag_gao = (/11.04,11.15,6.61,4.72/);gao
;volyear_sh_gao = (/1886 /);SH 把sh放入tp？！,1932,1991, 2011
;vol_sh_mag_gao = (/1.93/);gao
volyear_tp_gao = (/1883,1886,1902,1963,1982,1991/);1913,模式中貌似没有这个的火山强迫……
vol_tp_mag_gao = (/21.87,1.93,3.77,17.00,14.00,30.10/);gao

volyear_nh_sato = (/1891,1912/);NH (1891or92)
vol_nh_mag_sato = (/ 0.0495,0.0800/);sato 0.0495
volyear_tp_sato = (/1883,1886,1902,1963,1968,1982,1991/);1896(sh) (1968),
vol_tp_mag_sato = (/0.1860,0.0692,0.0931,0.1406,0.0491,0.1338,0.1627/);sato 0.0557 (0.0491)
ratio_tp = vol_tp_mag_gao(0)/vol_tp_mag_sato(0)
ratio_nh = vol_nh_mag_gao(0)/vol_nh_mag_sato(1)

vol_nh_sato = vol_nh_mag_sato*ratio_nh
vol_tp_sato = vol_tp_mag_sato*ratio_tp
vol_tp_gao = vol_tp_mag_gao
vol_nh_gao = vol_nh_mag_gao
print(vol_tp_gao)
print(vol_nh_gao)
print(vol_tp_sato)
print(vol_nh_sato)
volyear_c_nh_gao =  (/"1912","1925", "1943","1976","1980"/);
volyear_c_tp_gao =  (/"1883","1886", "1902","1963","1982","1991"/);
volyear_c_nh_sato = (/"1891","1912"/);
volyear_c_tp_sato = (/"1883","1886", "1902","1963","1968", "1982","1991"/);

volyear     = volyear_tp_sato
volyear_c = volyear_c_tp_sato
vol_mag =         vol_tp_sato
nv = dimsizes(volyear)

volmo = volyear ;(/1,1,1,1,1,1/);(/5,10,1,3,3,6/);
volmo = 1
volmo = volmo-1

ymstrt = 1870
ymlast2 = 2004
ny2 = ymlast2-ymstrt+1
print(ny2)
nino3 =  asciiread(diri+"CMIP5.nino34.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")
tt_a = asciiread(diri+"CMIP5.nino34.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");wait for shortwave!
;air_a_c = asciiread(diri+"CMIP5.part.top.sw.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");wait for shortwave!
;air_a_c =  asciiread(diri+"CMIP5.part.SSHAI.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");借用
;gmst_a =  asciiread(diri+"CMIP5.part.RDI.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");借用

gmst_obs = asciiread(diri+"CMIP5.GMSST.org.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")
gmst_a = asciiread(diri+"CMIP5.GMSST.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")

air_a_c =  asciiread(diri+"CMIP5.GMST.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");借用
air =  asciiread(diri+"CMIP5.GMST.org.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")
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

nino3_a_obs= dim_standardize_n_Wrap(nino3_a_obs, 0, 1)
tt_a_obs = dim_standardize_n_Wrap(tt_a_obs, 0, 1)
gmst_a_obs = dim_standardize_n_Wrap(gmst_a_obs, 0, 1)
air_a = dim_standardize_n_Wrap(air_a, 0, 1)

tt_obs_dt = dtrend(tt_a_obs, True) ;rightmost dimension
nino3_obs_dt = dtrend(nino3_a_obs, True)
g_obs_dt =     dtrend(gmst_a_obs, True)   ;air,True)   
air_dt   =dtrend(air_a, True)

cor = escorc(air_a, gmst_a_obs) 
print(cor)

;---------------------------------------------
nyear1 = 1
nyear2 = 5
nyr = nyear1+nyear2
nmon = nyear1*12+nyear2*12
n3_all = new((/nm,nv,nyr/), float)
n_period = new((/nm,nv/), float)

do nn=0,nm-1;1;
nino3_om =nino3_obs_dt(nn,:)
tt_om = tt_obs_dt(nn,:)
g_om =g_obs_dt(nn,:)
a_om =air_dt(nn,:)  
;-----------moving average---------------------
nrun = 3
 endopt    = 1 ;起始结束不损失滑动的月份

  nino3_om = runave_n_Wrap (nino3_obs_dt(nn,:), nrun, endopt, 0)
  t_om =    runave_n_Wrap (tt_obs_dt(nn,:), nrun, endopt, 0)
  g_om =     runave_n_Wrap (g_obs_dt(nn,:), nrun, endopt, 0)
  a_om =runave_n_Wrap (air_dt(nn,:), nrun, endopt, 0)
printVarSummary(tt_om)


;---------volcano year-----------------


n3_ov =    new((/nv,nmon/),float)
t_ov =       new((/nv,nmon/),float)
g_ov =        new((/nv,nmon/),float)
a_ov =       new((/nv,nmon/),float)

do i = 0,nv-1
  kk = volyear(i)-ymstrt
  aa = kk*12-nyear1*12+volmo(i)
  bb = kk*12+nyear2*12-1+volmo(i)

  t_ov(i,:)    =  t_om(aa:bb)
  n3_ov(i,:) =  nino3_om(aa:bb)
  
  g_ov(i,:) =     g_om(aa:bb)  

  a_ov(i,:) =     a_om(aa:bb)  

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
e_n3_clm = new((/ne,12/),float)
e_t_clm = new((/ne,12/),float)
e_g_clm = new((/ne,12/),float)
e_a_clm = new((/ne,12/),float)
do i=0,11
  n3_clm(:,i)  = dim_avg_n_Wrap(n3_ov(:,i::12), 1)
  t_clm(:,i)      = dim_avg_n_Wrap(t_ov(:,i::12), 1)
  g_clm(:,i)      = dim_avg_n_Wrap(g_ov(:,i::12), 1)
  a_clm(:,i)      = dim_avg_n_Wrap(a_ov(:,i::12), 1)
;  e_n3_clm(:,i)  = dim_avg_n_Wrap(n3_e(:,i::12), 1)
;  e_t_clm(:,i)      = dim_avg_n_Wrap(t_e(:,i::12), 1)
;  e_g_clm(:,i)      = dim_avg_n_Wrap(g_e(:,i::12), 1)
;  e_a_clm(:,i)      = dim_avg_n_Wrap(a_e(:,i::12), 1)
end do 
printVarSummary(n3_clm)

n3_v_a  = new((/nv,nmon/),float)
t_v_a   = new((/nv,nmon/),float)
g_v_a   = new((/nv,nmon/),float)
a_v_a   = new((/nv,nmon/),float)
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

;n3_e_a(:,j*12+i)  = n3_e(:,i+j*12)- e_n3_clm(:,i)
;t_e_a(:,j*12+i)   = t_e(:,i+j*12)  -   e_t_clm(:,i)
;g_e_a(:,j*12+i)   = g_e(:,i+j*12)  -   e_g_clm(:,i)
;a_e_a(:,j*12+i)   = a_e(:,i+j*12)  -   e_a_clm(:,i)
end do 
end do 
printVarSummary(n3_v_a)

nino3_ov_avg =  dim_avg_n(n3_v_a, 0)
t_ov_avg =      dim_avg_n(t_v_a, 0)
g_ov_avg =      dim_avg_n(g_v_a, 0)
a_ov_avg =      dim_avg_n(a_v_a, 0)
;nino3_e_avg  =  dim_avg_n(n3_e_a, 0)
;t_e_avg  =      dim_avg_n(t_e_a, 0)
;g_e_avg  =      dim_avg_n(g_e_a, 0)
;a_e_avg  =      dim_avg_n(a_e_a, 0)
;--------------------period--------------------------
 ;  n3s = new((//))
do j=0,nv-1
n3_acr = esacr(n3_v_a(j,:),36)

n_period(nn,j) = ind(n3_acr.eq.min(n3_acr))
  ntt = esacr(n3_v_a(j,:),28)
  ptt = ind(ntt.eq.min(ntt))
if(n_period(nn,j).ge.30.and.ptt.lt.24) then 
    n_period(nn,j) = ptt
end if 
end do 


n3_v_a!0 = "vol"
n3_v_a!1 = "time"
n3_djf = new((/nv,nyr/),float)
do kk=0,nv-1
n3_djf(kk,:) = month_to_season(n3_v_a(kk,:), "DJF")
end do 

n3_all(nn,:,:) = n3_djf
end do 
printVarSummary(n3_all)
printVarSummary(n_period)
;-------------sort----------------------
nvA = new((/na,nv,nyr/), float)
nvB = new((/nb,nv,nyr/), float)
nvC = new((/nc,nv,nyr/), float)
pA = new((/na,nv/), float)
pB = new((/nb,nv/), float)
pC = new((/nc,nv/), float)
do i =0,na-1
nvA(i,:,:) = n3_all(A(i),:,:)
pA(i,:) = n_period(A(i),:)

end do

do i =0,nb-1
nvB(i,:,:) = n3_all(B(i),:,:)
pB(i,:) = n_period(B(i),:)
end do

do i =0,nc-1
nvC(i,:,:) = n3_all(C(i),:,:)
pC(i,:) = n_period(C(i),:)
end do
printVarSummary(nvA)
;--------plot-------------------------
;wks2 = gsn_open_wks("png", diri+"cmip.vol(NH).nino.scatter."+name_model(nn));&el
;wks2 = gsn_open_wks("ps", "/work/obs.eachvol.nino34&air.ERSST&HadISST")
;x_axis = ispan(-nyear1*12, nyear2*12-1, 1)
;x_axis = fspan(ymstrt, ymlast2, ny2*12)
;do j=0,nv-1
;plot = new(3,graphic)
;do i=0,nm-1

res = True
res@gsnFrame = False
res@gsnDraw  = False
 
res@xyMarkLineMode     = "Markers"  ; choose which have markers
; res@xyDashPatterns      = 1                      ; solid line 
  res@xyMarker         =  16                      ; choose type of marker  
  res@xyMarkerSizeF     = 0.03                     ; Marker size (default 0.01)
;   res@xyLineThicknesses   = (/1,2/)                ; set second line to 2

res@xyMarkerOpacityF = 0.7;0.5;
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
res@trYMinF            = -3.;-0.5;
res@trYMaxF            =  3.;0.5;
res@gsnXRefLine       = 0;volyear;(/1980,1988/);1900;
res@gsnXRefLineThicknessF = 2.
res@gsnXRefLineColor = "black";(/"grey","grey","grey","grey","grey","grey"/);
res@gsnXRefLineDashPattern = 0;0

res@gsnYRefLine       = 0.  ;Xrefline和Yrefline不能共存？！
res@gsnYRefLineColor = "black"
res@gsnYRefLineDashPattern = 0;0
res@xyLineThicknessF = 2.




res@xyMonoMarkerColor  =  True ;False          ; Allow different colors for markers
;colors = (/"lightsalmon1","LightBlue","lightgoldenrod3","limegreen","orchid3","lightslategray"/)
;gsn_define_colormap(wks2,colors)
;res@xyMarkerColors = (/"lightsalmon1","LightBlue","lightgoldenrod3","limegreen","orchid3",\
;                      "lightslategray","brown","red2","skyblue","black"/);colors
res@xyMarkerSizes = vol_mag*0.002 
colors = (/"lightsalmon1","LightBlue","lightgoldenrod3","limegreen","orchid3","blue2",\
  "lightslategray","brown","red2","orangered","goldenrod1","hotpink","greenyellow","lemonchiffon",\
  "seagreen1","plum3"/);colors
res@gsnLeftString = "";"(a)";"
res@tiXAxisString = "nino34 in year(0)"
res@tiYAxisString = "nino34 in year(1)"

 resk = True
resk@gsnMaximize = True
resk@gsnPaperOrientation = "portrait"   ; force portrait
 
      resk@txString = "volcano(TP)(gao) in A";+name_model(nn); with Elnino(grey)
          nsort = nvA
          psort =  pA

X = nsort(0,:,1)
Y = nsort(0,:,2)
;res@xyMarkerColor = colors(0);(is)
res@xyMarkerColors = colors
res@xyMonoMarkerColor  =  False ;True ;         ; Allow different colors for markers
plot(0) = gsn_csm_xy(wks2,transpose((/X/)),transpose((/Y/)),res);(wks2,X,Y,res); ,g_e_a,res3
res@tiYAxisString = "nino34 in year(2)"
Y = nsort(0,:,3)
  plot(1) = gsn_csm_xy(wks2,transpose((/X/)),transpose((/Y/)),res); ,g_e_a,res3
res@tiYAxisString = "nino34 period"
res@trYMinF            =  5.;-0.5;
res@trYMaxF            =  40.;0.5;
Y = psort(0,:)
  plot(2) = gsn_csm_xy(wks2,transpose((/X/)),transpose((/Y/)),res) ; n3_e_a,,res3
do is = 1,ns-1
X = nsort(is,:,1)
Y = nsort(is,:,2)
;res@xyMarkerColor = colors(is)
res@xyMarkerColors = colors
res@trYMinF            = -3.;-0.5;
res@trYMaxF            =  3.;0.5;
plot1 = gsn_csm_xy(wks2,transpose((/X/)),transpose((/Y/)),res);(wks2,X,Y,res) ,g_e_a,res3
res@tiYAxisString = "nino34 in year(2)"
Y = nsort(is,:,3)
  plot3 = gsn_csm_xy(wks2,transpose((/X/)),transpose((/Y/)),res); ,g_e_a,res3
res@tiYAxisString = "nino34 period"
Y = psort(is,:)
res@trYMinF            =  5.;-0.5;
res@trYMaxF            =  40.;0.5;
  plot5 = gsn_csm_xy(wks2,transpose((/X/)),transpose((/Y/)),res) ; n3_e_a,,res3
overlay(plot(0), plot1) 
overlay(plot(1), plot3)
overlay(plot(2), plot5)
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
  lgres@lgMarkerSizeF  = 0.004
  lgres@lgItemType         = "Markers";"Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .06            ; legend label font thickness
  lgres@vpWidthF           = 0.11           ; width of legend (NDC)
  lgres@vpHeightF          = 0.13           ; height of legend (NDC)
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter
  lgres@lgMonoMarkerIndex = True
  lgres@lgMarkerIndex = 16

  labels = volyear_c;name_model(sort_name) ;
;TP
  nitems = dimsizes(labels)                 ; number of legend items

; Create legend
  lbid = gsn_create_legend(wks2,nitems,labels,lgres)
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
 ; annoid = gsn_add_annotation(plot(1),lbid,amres)  ; add legend to plot
 ; annoid = gsn_add_annotation(plot(2),lbid,amres)  ; add legend to plot

delete(res)
delete(res2)
delete(res3)
delete(res4)
;print(nn)

;end do 
gsn_panel(wks2, plot, (/3,1/), resk)
end