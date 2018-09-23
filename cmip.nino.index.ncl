load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; Nino 1+2 (0-10S, 90W-80W), Nino 3 (5N-5S, 150W-90W), 
; Nino 3.4/ONI (5N-5S, 170W-120W) and Nino 4 (5N-5S, 160E-150W)
; TNI is Nino 1+2 - Nino 4
begin
;  clatS  = -5.0;85*1.6 ; or 85-95?!
;  clatN  = 5.0; 95*1.6
;  clonL  = 210.0;190.0;120.0; 
;  clonR  = 270.0;240.0;250.0;
;
;  elatS  = -5.0
;  elatN  = 5.0
;  elonL  = 120 ;RDI
;  elonR  = 270
NLAT =  10
SLAT =  -10
WLON =   190
ELON = 240

  ymstrt= 187001  
  ymlast= 200412
  yrstrt= 1870  
  yrlast= 2004

;"bcc-csm1-1-m","bcc-csm1-1","CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM","EC-EARTH", "GISS-E2-H-CC",
;"HadGEM2-AO","IPSL-CM5A-LR",
diri = "/data2/data/cmip5/historical_regrid/atmos/mon/tos/"
diri2 = "/data2/data/cmip5/historical_regrid/atmos/mon/tas/"
diri3 = "/data2/home/liufeiclim2014/xing/cmip_vol/"
diri4 = "/data2/data/cmip5/historical/atmos/mon/"

;diri = "/nuist/p/work/feil/xc/"
;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\
;            "CNRM-CM5","CNRM-CM5-2","CSIRO-Mk3-6-0", \;  
;             "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC",\;
;             "HadGEM2-ES",\;"HadGEM2-AO","GISS-E2-H","HadCM3","HadGEM2-CC","MIROC-ESM","MRI-CGCM3",
;              "MIROC5","MIROC-ESM-CHEM",\
;              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P",\
;              "MRI-ESM1","NorESM1-M","NorESM1-ME"/); rsdscs

;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\;"bcc-csm1-1-m","bcc-csm1-1",
;       "CNRM-CM5-2",  \;"CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM",
;      "CNRM-CM5","CSIRO-Mk3-6-0",  \
;        "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
;        "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES","MIROC5","MIROC-ESM-CHEM",\
;         "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);thetao
;
name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\;
            "CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM",\;
            "CNRM-CM5-2",  "CNRM-CM5","CSIRO-Mk3-6-0",  \
             "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
             "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES",\
             "MIROC5","MIROC-ESM-CHEM",\
              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);;"HadGEM2-CC"(没有tas资料）,"FIO-ESM","FGOALS-g2",  "MRI-CGCM3",,"GISS-E2-H"
;"HadGEM2-CC"(没有tas资料）,"FIO-ESM","FGOALS-g2",  "MRI-CGCM3",,"GISS-E2-H"
; "CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM","CNRM-CM5-2",  "CNRM-CM5","EC-EARTH","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","HadGEM2-AO","MIROC5","MRI-ESM1",
nm = dimsizes(name_model)

nlev = 17

nyr = yrlast-yrstrt+1
SSHAI = new((/nm,nyr*12/), float)
RDI = new((/nm,nyr*12/), float)
gg = new((/nm,nyr*12/), float)
ga = new((/nm,nyr*12/), float)
ga_l = new((/nm,nyr*12,nlev/), float)
nino4 = new((/nm,nyr*12/), float)
eof_ts= new((/nm,nyr*12/), float)
eof= new((/nm,nyr*12,24,77/), float)
t2 = new((/nm,nyr*12/), float)
;sst = new((/36,1872,72,144/),  typeof(ss),ss@_FillValue)
;files = systemfunc(" ls -1 "+diri3+"thetao/thetao_10level_TP_*.nc")
;print(files)
;nfiles = dimsizes(files)
;print(nfiles)
;print(nm)
do i =0,nm-1
;files = systemfunc(" ls -1 "+diri3+"rsdt/rsdt_Amon_"+name_model(i)+"_*.nc");
;files = systemfunc(" ls -1 "+diri4+"ta/"+name_model(i)+"/r1i1p1/*.nc")
;files = systemfunc(" ls -1 /data2/data/cmip5/historical_regrid/atmos/mon/hgt/"+name_model(i)+"*.nc")
;print(files)
;files = systemfunc(" ls -1 "+diri2+name_model(i)+".tas.*.nc")
;f = addfile(diri2+name_model(i)+".tas.historical.mon.nc","r")

;f   = addfiles(files,"r") ;(i)
;f   = addfile(diri3+"thetao/thetao_10level_TP_"+name_model(i)+".nc", "r")
f   = addfile(diri+name_model(i)+".tos.historical.mon.nc", "r")
TIME=f->time;
yyyymm=cd_calendar(TIME, -1)
ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)
print(ymstrt1)
print(ymlast1)
;if (name_model(i).eq."GISS-E2-H-CC".or.name_model(i).eq."GISS-E2-R".or.name_model(i).eq."GISS-E2-R-CC") then
; x  = f[:]->rsdt(ymstrt1(0):1859,:,:)
;else 
;x  = f[:]->rsdscs(ymstrt1(0):ymstrt1(0)+1859-240,{SLAT:NLAT},{WLON:ELON});f[:]ymlast1(0)
;x  = f[:]->thetao(ymstrt1(0):ymlast1(0),{100:200},{SLAT:NLAT},{WLON:ELON});:,:)
;end if

;xx = x 
;if(x!1.eq."lat") then
;x1 = x(:,:,{SLAT:NLAT},{WLON:ELON});:,:);
;else 
;x@lat2d =  f[:]->lat
;x@lon2d =  f[:]->lon
;printVarSummary(x)
;    Opt                  = True
;   Opt@ForceOverwrite      = True
;   Opt@SrcInputFileName = f[:]
;   Opt@SrcTitle          = "NCEP Grid"
;   Opt@WgtFileName      = diri3+"a.nc"   ; destination file name
;   Opt@InterpMethod      =    "patch"
;
;   Opt@DstGridType      = "2x2"              ; Destination grid
;   Opt@DstTitle         = "World Grid 2-degree Resolution"
;   Opt@DstLLCorner      = (/-89.75d,   0.00d /)
;   Opt@DstURCorner      = (/ 89.75d, 359.75d /) 
;    ;;Opt@PrintTimings   = True
;   ;;Opt@Debug          = True
;      zz = ESMF_regrid(x,Opt) 
;      printVarSummary(zz)
;x1  = zz(:,:,{SLAT:NLAT},{WLON:ELON})
;    delete(Opt)
;delete(zz)
;print(i)
;end if
;printVarSummary(xx)
;    system(" rm "+diri3+"source_grid_file.nc" )
;    system(" rm "+diri3+"destination_grid_file.nc" )
;    system(" rm "+diri3+"a.nc" )
;copy_VarCoords(sst, x)

;x1 = new((/dimsizes(x&time),dimsizes(x&j),dimsizes(x&i)/),float,x@_FillValue)
;      nlat = dimsizes(x&j)
;      lat  = latGau(nlat, "lat", "latitude", "degrees_north")
;      mlon = dimsizes(x&i)
;      lon  = lonGlobeF(mlon, "lon", "longitude", "degrees_east")
;x1!1 = "lat"
;x1&lat = lat
;x1!2 = "lon"
;x1&lon = lon
;x1 = x 
;1&lat@units=sst&lat@units
;1&lat = fspan(-90, 90, 300)

;1&lon@units=sst&lon@units
;1&lon = fspan(0, 360, 360)
;f = addfile(diri+name_model(i)+".tos.historical.mon.nc", "r")
;TIME1=f->time
;yyyymm=cd_calendar(TIME1, -1)
;
;ymstrt1 =ind(yyyymm.eq.ymstrt)
;ymlast1 = ind(yyyymm.eq.ymlast)
;if(name_model(i).eq."FGOALS-g2")then
;ymstrt1 = (yrstrt-1850)*12
;ymlast1 = (yrlast-1850+1)*12-1
;end if 
print(ymstrt1)
print(ymlast1)
;;--------thetao--------------------------------------------
;if (name_model(i).eq."bcc-csm1-1-m".or.name_model(i).eq."bcc-csm1-1".or.name_model(i).eq."CESM1-BGC"\
;   .or.name_model(i).eq."CESM1-CAM5-1-FV2".or.name_model(i).eq."CESM1-FASTCHEM".or.name_model(i).eq."CESM1-WACCM") then
;  x = f->thetao(240:1859,{0:300}, {SLAT:NLAT},{WLON:ELON})
;else
;x = f->thetao(ymstrt1:ymlast1,{0:300}, {SLAT:NLAT},{WLON:ELON});:,:);{clatS:clatN},{clonL:clonR});136:152,190:240)

;(0) 100000(1) 92500(2) 85000(3) 70000(4) 60000(5) 50000(6) 40000(7) 30000(8) 25000(9) 20000(10)  15000(11)  10000(12)  7000(13)  5000(14)  3000(15)  2000(16)  1000
;卧槽！不能用这个！有的模式层数怎么不一样？！希望只是高度不一样，而不是间隔…… LEV = 13
;貌似只是low top 和 high top的差异……
x = f->tos(ymstrt1(0):ymstrt1(0)+1619,{SLAT:NLAT},{WLON:ELON});ymlast1;ymstrt1(0)+1619
;x1 = x
;x1 = x-273.15
; x1 = lonPivot(x,0.5)
;end if
;printVarSummary(x)
xx = x
xc=clmMonTLL(xx)
xa = calcMonAnomTLL(xx, xc)
printVarSummary(xa)

;t1 = wgt_areaave_Wrap(xa, 1.0, 1.0, 1)
;printVarSummary(t1)
;t2(i,:) = dim_avg_n_Wrap(t1, 1)
;delete(t1)
;  ddx = 6378388.*cos(0.0174533*0)*0.0174533
;  x_tend = center_finite_diff_n(xa(:,:,:),ddx,False,0,2)
;  copy_VarMeta(xa, x_tend)
;  printVarSummary(x_tend)

;  SSHAI(i,:)   = wgt_areaave_Wrap(xa(:,{clatS:clatN},{clonL:clonR}),1.0, 1.0, 1)
;RDI(i,:) = wgt_areaave_Wrap(x_tend(:,{clatS:clatN},{clonL:clonR}),1.0, 1.0, 1)
lat = x&lat
rad    = 4.0*atan(1.0)/180.0
clat   = cos(lat*rad)
ga(i,:) = wgt_areaave_Wrap(xa,clat, 1.0, 1)
;ga_l(i,:,:) = wgt_areaave_Wrap(xa,clat, 1.0, 1)

;gg(i,:) = wgt_areaave_Wrap(x,clat, 1.0, 1) ;top radiation不应该去掉年循环！！！
nino4(i,:) = wgt_areaave_Wrap(xa,1.0, 1.0, 1);(:,{clatS:clatN},{clonL:clonR})
;nino4(i,:) = dim_standardize_n_Wrap(nino4(i,:), 0, 0)
;nino12(i,:) = dim_standardize_n_Wrap(nino12(i,:), 0, 0)
;--------------EOF-------------------
 ; optEof = True
 ; neof = 1
 ; eof(i,0,:,:)    = eofunc_Wrap( xa(lat|:,lon|:,time|:), neof, optEof)    ;eofunc_n_Wrap 6.4.0
 ; eof_ts(i,:) = eofunc_ts_Wrap( xa(lat|:,lon|:,time|:), eof, False)   ;if eof_ts(x) rightmost must be time!!!!
 ; printVarSummary(eof)
 ; if(avg(eof(0,{0},{200:240})).lt.0)then
 ;   eof_ts(i,:) = (-1)*eof_ts(i,:)
 ;   print(i)
 ; end if

print(i)
delete([/x,xx,xc,xa,TIME,yyyymm,ymstrt1,ymlast1,f,lat,clat/]);,t1files,x_tend,,t1
;delete([/x,xx,xc,xa,TIME,yyyymm,ymstrt1,ymlast1,f,lat,clat/]);,eof

end do
;printVarSummary(t2)
;air = new((/nm,nyr*12/), float)
;aa = new((/nm,nyr*12/), float)
;
;do i =0,nm-1
;
;f2 = addfile(diri2+name_model(i)+".tas.historical.mon.nc", "r")
;TIME1=f2->time
;yyyymm=cd_calendar(TIME1, -1)
;
;ymstrt1 =ind(yyyymm.eq.ymstrt)
;ymlast1 = ind(yyyymm.eq.ymlast)
;;if(name_model(i).eq."FGOALS-g2")then
;;ymstrt1 = (yrstrt-1850)*12
;;ymlast1 = (yrlast-1850+1)*12-1
;;end if 
;print(ymstrt1)
;print(ymlast1)
;printMinMax(TIME1, 0)
;printMinMax(yyyymm, 0)
;x = f2->tas(ymstrt1:ymlast1,:,:);{clatS:clatN},{clonL:clonR});136:152,190:240)
;xc=clmMonTLL(x)
;xa = calcMonAnomTLL(x, xc)
;
;lat = x&lat
;rad    = 4.0*atan(1.0)/180.0
;clat   = cos(lat*rad)
;aa(i,:) = wgt_areaave_Wrap(xa,clat, 1.0, 1)
;air(i,:) = wgt_areaave_Wrap(x,clat, 1.0, 1)
;;nino4(i,:) = dim_standardize_n_Wrap(nino4(i,:), 0, 0)
;;nino12(i,:) = dim_standardize_n_Wrap(nino12(i,:), 0, 0)
;print(i)
;delete(x)
;delete(TIME1)
;delete(yyyymm)
;delete(ymstrt1)
;delete(ymlast1)
;end do
;gg1 = dim_standardize_n_Wrap(gg, 0, 1)

;printVarSummary(gg1)
;asciiwrite("CMIP5.part.top.sw.org.1870-2004.txt", gg)
;asciiwrite("CMIP5.pc1.1870-2004.txt", eof_ts)

asciiwrite("CMIP5.nino34(10).anomaly.unnorm.1870-2004.txt", nino4)
;asciiwrite("CMIP5.GMSST.anomaly.unnorm.1870-2004.txt", ga)
;asciiwrite("CMIP5.air.NA.anomaly.unnorm.1870-2004.txt", ga)
;asciiwrite("CMIP5.part.SSHAI.anomaly.unnorm.1870-2004.txt", SSHAI)
;asciiwrite("CMIP5.part.RDI.anomaly.unnorm.1870-2004.txt", RDI)
;asciiwrite("obs.nino3.anomly.filter.1yr.1854-2005.std.txt", var_obs) ;1940-2016.
;asciiwrite("CMIP5.rsdscs.anomaly.unnorm.1870-2004.22.txt", ga)
;asciiwrite("CMIP5.thetao100-150.160-270.anomaly.unnorm.1870-2004.txt", t2)
;asciiwrite("CMIP5.thetao0-300.120-280.anomaly.unnorm.1870-2004.txt", t2)
;asciiwrite("CMIP5.nino3.anomaly.unnorm.1870-2004.txt", nino4)
;asciiwrite("CMIP5.TP.50ta.anomaly.unnorm.1870-2004.txt", nino4)
;asciiwrite("CMIP5.polarvortex.50hgt.anomaly.unnorm.1870-2004.txt", ga)
;asciiwrite("CMIP5.polarvortex.50hgt.alllev.anomaly.unnorm.1870-2004.txt", ga_l)
;asciiwrite("CMIP5.air.SH30-90.anomaly.unnorm.1870-2004.txt", ga)

;asciiwrite("CMIP5.air.MC.anomaly.unnorm.1870-2004.txt", ga)

;  N    = 5    ; number of rows
;  M    = (nyear1+nyear2)*12     ; number of columns (N,M)
;  fWidth = 10    ; specify the format width
;  fDec   = 4     ; specify the number to the right of decimal point
;  fmtx   = M + "f" + fWidth + "." + fDec              ; fmtx="50f10.3"
;  opt = True
;  opt@fout = "CMIP5.polarvortex.50hgt.alllev.anomaly.unnorm.1870-2004.txt"
;  write_matrix (ga_l, fmtx, opt)
end