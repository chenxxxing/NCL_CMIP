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

volyear_c = (/"1883", "1902","1913", "1963", "1982","1991"/);

nv = dimsizes(volyear)

volmo = volyear ;(/1,1,1,1,1,1/);(/5,10,1,3,3,6/);
volmo = 1
volmo = volmo-1

ymstrt = 1870
ymlast2 = 2004
ny2 = ymlast2-ymstrt+1
print(ny2)
;nino3 =  asciiread(diri+"CMIP5.nino34.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")
nino3 =  asciiread(diri+"CMIP5.nino3.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")
;nino3 =  asciiread(diri+"CMIP5.nino34.anomaly.unnorm.1870-2004.thetao.txt", (/nm,ny2*12/), "float")
;tt_a = asciiread(diri+"CMIP5.thetao100-150.nino34.anomaly.unnorm.1870-2004.23.txt", (/nm,ny2*12/), "float");wait for shortwave!
;tt_a = asciiread(diri+"CMIP5.i20.nino34.anomaly.unnorm.1870-2004.23.txt", (/nm,ny2*12/), "float");wait for shortwave!
;tt_a = asciiread(diri+"CMIP5.nino34.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");wait for shortwave!
;air_a_c = asciiread(diri+"CMIP5.part.top.sw.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");wait for shortwave!

;gmst_obs = asciiread(diri+"CMIP5.GMSST.org.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")
;gmst_a = asciiread(diri+"CMIP5.GMSST.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")

;air_a_c =  asciiread(diri+"CMIP5.GMST.anomaly.unnorm.1870-2004.thetao.txt", (/nm,ny2*12/), "float");借用
air_a_c =  asciiread(diri+"CMIP5.GMST.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");借用

gmst_a = asciiread(diri+"CMIP5.air.EU.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float")
tt_a =  asciiread(diri+"CMIP5.air.NA.anomaly.unnorm.1870-2004.txt", (/nm,ny2*12/), "float");借用

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
n1  = new((/nm,ny2*12/), float)
a1  = new((/nm,ny2*12/), float)
t1  = new((/nm,ny2*12/), float)
g1  = new((/nm,ny2*12/), float)
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
fca   = 1./(9*12) ;                        
fcb   = -999;
wgt  = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
nino3_a1 = wgt_runave_n ( nino3_a_obs(nn,:), wgt, 0, 0 ) ;最后一个是指sst1第0维
tt_a1 = wgt_runave_n ( tt_a_obs(nn,:), wgt, 0, 0 ) ;最后一个是指sst1第0维
air_a1 = wgt_runave_n ( air_a(nn,:), wgt, 0, 0 ) ;最后一个是指sst1第0维
g_a1 = wgt_runave_n ( gmst_a_obs(nn,:), wgt, 0, 0 )
;-----------moving average---------------------
nrun = 3
 endopt    = 1 ;起始结束不损失滑动的月份

  n1(nn,:) = runave_n_Wrap (nino3_a1, nrun, endopt, 0)
  t1(nn,:) =    runave_n_Wrap (tt_a1, nrun, endopt, 0)
  g1(nn,:) =     runave_n_Wrap (g_a1, nrun, endopt, 0)
  a1(nn,:) =runave_n_Wrap (air_a1, nrun, endopt, 0)
;printVarSummary(tt_om)

end do

;--------plot-------------------------
wks2 = gsn_open_wks("ps", diri+"cmip.scatter.filter");+name_model(nn));.A2");&el
plot = new(2,graphic)
;wks2 = gsn_open_wks("png", diri+"cmip.vol(NH).nino.scatter."+name_model(nn));&el
;wks2 = gsn_open_wks("ps", "/work/obs.eachvol.nino34&air.ERSST&HadISST")
;x_axis = ispan(-nyear1*12, nyear2*12-1, 1)
;x_axis = fspan(ymstrt, ymlast2, ny2*12)
;do j=0,nv-1
;plot = new(3,graphic)
do ik=0,nm-1

res = True
res@gsnFrame = False
res@gsnDraw  = False
 
res@xyMarkLineMode     = "Markers"  ; choose which have markers
; res@xyDashPatterns      = 1                      ; solid line 
  res@xyMarker         =  16                      ; choose type of marker  
  res@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
;   res@xyLineThicknesses   = (/1,2/)                ; set second line to 2

res@xyMarkerOpacityF = 0.5;0.5;
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
res@trXMinF            = -3.0
res@trXMaxF            =  3.0
res@trYMinF            = -0.3;-2.0;
res@trYMaxF            =  0.3;2.0;
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
res@xyMonoMarkerColor  = True ; False ;         ; Allow different colors for markers
;colors = (/"lightsalmon1","LightBlue","lightgoldenrod3","limegreen","orchid3","lightslategray"/)
;gsn_define_colormap(wks2,colors)
;res@xyMarkerColors = (/"lightsalmon1","LightBlue","lightgoldenrod3","limegreen","orchid3",\
;                      "lightslategray","brown","red2","skyblue","black"/);colors
;res@xyMarkerSizes = vol_mag*0.003 
colors = (/"lightsalmon1","LightBlue","lightgoldenrod3","limegreen","orchid3", "lightslategray",\
 "brown","blue2","red2","orangered","goldenrod1","hotpink","greenyellow","lemonchiffon"/);colors
;colors = read_colormap_file("cb_9step") ;("BkBlAqGrYeOrReViWh200");("default");("cosam");"GMT_paired"

res@gsnLeftString = "";"(a)";"

 resk = True
resk@gsnMaximize = True
resk@gsnPaperOrientation = "portrait"   ; force portrait
res2@xyMonoMarkerColor  = False ;True ;          ; Allow different colors for markers
res2@xyMarkerColors = colors;"grey70";
res2@xyMarkerSizeF     = 0.02  
res2@xyMarkerOpacityF = 1.0;0.5;

res@xyMarkerColor = "black" ;"grey70";"blue" ;colors(::3,:);(0);(is)

resk@txString = "box air vs GMST(filter) "+name_model(ik);;A"; with Elnino(grey)
res@tiXAxisString = "EU air"
res@tiYAxisString = "GMST"
res@trXMinF             = -5.
res@trXMaxF             =  5.
res@trYMinF             = -0.8;-3.0;
res@trYMaxF             =  0.8;3.0;
;res3@trXMinF            = -3.
;res3@trXMaxF            =  3.
;res3@trYMinF            = -3.0;-0.5;
;res3@trYMaxF            =  3.0;0.5;
plot(0) = gsn_csm_xy(wks2, transpose((/g1(ik,:)/)),transpose((/a1(ik,:)/)),res) 

;plot(0) = gsn_csm_xy(wks2, transpose((/gone1/)),transpose((/none1/)),res) 
res@tiXAxisString = "NA air"
res@trXMinF            = -5.0
res@trXMaxF            =  5.0
plot(1) = gsn_csm_xy(wks2, transpose((/t1(ik,:)/)),transpose((/a1(ik,:)/)),res)
;plot(1) = gsn_csm_xy(wks2, transpose((/tone1/)),transpose((/none1/)),res) 

;res@tiXAxisString = "nino3 in year(0)"
;res@tiYAxisString = "nino3 in year(1)"
;res@trXMinF            = -3.
;res@trXMaxF            =  3.
;res@trYMinF            = -3.;-0.5;
;res@trYMaxF            =  3.;0.5;
;res3@trXMinF           = -3.
;res3@trXMaxF           =  3.
;res3@trYMinF           = -3.;-0.5;
;res3@trYMaxF           =  3.;0.5;
;plot(1) = gsn_csm_xy(wks2, transpose((/none11/)),transpose((/none21/)),res) 
;res3@xyMarkerColor = "red" ;
;plot3 = gsn_csm_xy(wks2, transpose((/none12/)),transpose((/none22/)),res3) 
;res3@xyMarkerColor = "blue" ;
;plot4 = gsn_csm_xy(wks2, transpose((/none03/)),transpose((/tone13/)),res3) 
;plot0 = gsn_csm_xy(wks2,transpose(t5),transpose(t6),res2)
;plot1 = gsn_csm_xy(wks2,transpose(n5),transpose(n6),res2)
;overlay(plot(0), plot1)
;overlay(plot(0), plot2)
;overlay(plot(1), plot3)
;overlay(plot(1), plot4)
;resk@txString = "volcano nino3.4&mixed layer(DJF)";+name_model(ik);A"; with Elnino(grey)
;res@tiYAxisString = "mixed layer in year(-1)"
;res@tiXAxisString = "nino3.4 in year(-1)"
;plot(0) = gsn_csm_xy(wks2, transpose((/none01/)),transpose((/tone01/)),res) ;transpose((/n3_all_v(:,6)/)),transpose((/a_all_v(:,6)/)),res)
;res@tiYAxisString = "mixed layer in year(0)"
;res@tiXAxisString = "nino3.4 in year(0)"
;plot(1) = gsn_csm_xy(wks2, transpose((/none11/)),transpose((/tone11/)),res) ;transpose((/n3_all_v(:,6)/)),transpose((/a_all_v(:,6)/)),res)
;res@tiYAxisString = "mixed layer in year(1)"
;res@tiXAxisString = "nino3.4 in year(1)"
;;res@xyMarkerColors = "grey30";colors;(0);(is)
;plot(2) = gsn_csm_xy(wks2, transpose((/none21/)),transpose((/tone21/)),res) ;transpose((/n3_all_v(:,7)/)),transpose((/a_all_v(:,7)/)),res)

;res3@xyMarkerColor = "green" ;
;plot02 = gsn_csm_xy(wks2, transpose((/none02/)),transpose((/tone02/)),res3) ;transpose((/n3_all_v(:,6)/)),transpose((/a_all_v(:,6)/)),res)
;plot12 = gsn_csm_xy(wks2, transpose((/none12/)),transpose((/tone12/)),res3) ;transpose((/n3_all_v(:,6)/)),transpose((/a_all_v(:,6)/)),res)
;plot22 = gsn_csm_xy(wks2, transpose((/none22/)),transpose((/tone22/)),res3) 
;overlay(plot(0), plot02)
;overlay(plot(1), plot12)
;overlay(plot(2), plot22)
;res3@xyMarkerColor = "blue" ;
;plot03 = gsn_csm_xy(wks2, transpose((/none03/)),transpose((/tone03/)),res3) ;transpose((/n3_all_v(:,6)/)),transpose((/a_all_v(:,6)/)),res)
;plot13 = gsn_csm_xy(wks2, transpose((/none13/)),transpose((/tone13/)),res3) ;transpose((/n3_all_v(:,6)/)),transpose((/a_all_v(:,6)/)),res)
;plot23 = gsn_csm_xy(wks2, transpose((/none23/)),transpose((/tone23/)),res3) 
;overlay(plot(0), plot03)
;overlay(plot(1), plot13)
;overlay(plot(2), plot23)
;ploto0 = gsn_csm_xy(wks2,transpose(n6),transpose(a6),res2) ;transpose((/nobs_djf(:,6)/)), transpose((/aobs_djf(:,6)/)), res2)
;ploto1 = gsn_csm_xy(wks2,transpose(n6),transpose(n7),res2) ;transpose((/nobs_djf(:,6)/)), transpose((/aobs_djf(:,6)/)), res2)
;;plot2 = gsn_csm_xy(wks2,transpose(n7),transpose(t7),res2) ;transpose((/nobs_djf(:,7)/)), transpose((/aobs_djf(:,7)/)), res2)
;overlay(plot(0), ploto0)
;overlay(plot(1), ploto1)
;overlay(plot(2), plot2)
;------------------EL-----------------------------
;elm = asciiread(diri+"cmip.all.enso.year.txt",(/nm*70/),"float")
;EL = onedtond(elm, (/nm,70/))
;
;do im=0,nm-1
;ne = num(EL(im,:).gt.1e+20)
;nel = new((/ne,nmon/), float)
;nel_clm = new((/ne,12/),float)
;nela  = new((/ne,nmon/),float)
;;ymvol3 = new(ne, integer)
;do ie =0,ne-1
;    ;ymvol3(ie) = ind(yyyymms3.eq.EL(im,ie));EL(iv));
;      kk = EL(im,ie)-ymstrt
;      aa = kk*12-nyear1*12
;      bb = kk*12+nyear2*12-1
;    nel(ie,:) = nino3_om(im,aa:bb)
;    do j=0,11
;    nel_clm(ie,j) = dim_avg_n_Wrap(nel(ie,j::12), 1)
;    end do 
;
;    do i = 0,11
;    do j = 0,nyr-1
;    nela(ie,j*12+i)  = nel(ie,i+j*12)-nel_clm(ie,i)
;    end do 
;    end do 
;
;end do 
;delete([/nel,nel_clm,nela/])
;end do 

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
  lgres@lgMarkerColors       = colors;(/"red","green" /) ; "blue"
  lgres@lgBoxLineDashPattern = 1;(/1,0,0,1,1/)
  lgres@lgMarkerSizeF  = 0.01
  lgres@lgItemType         = "Markers";"Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .06            ; legend label font thickness
  lgres@vpWidthF           = 0.11           ; width of legend (NDC)
  lgres@vpHeightF          = 0.13           ; height of legend (NDC)
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter
  lgres@lgMonoMarkerIndex = True
  lgres@lgMarkerIndex = 16

  labels =  volyear_c;(/"A","B"/) ;name_model(sort_name) ;,"C"
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

 ; annoid = gsn_add_annotation(plot(0),lbid,amres)  ; add legend to plot
 ; annoid = gsn_add_annotation(plot(1),lbid,amres)  ; add legend to plot
 ; annoid = gsn_add_annotation(plot(2),lbid,amres)  ; add legend to plot

delete(res)
delete(res2)
delete(res3)
delete(res4)
;print(nn)

gsn_panel(wks2, plot, (/2,1/), resk)
end do 
end
