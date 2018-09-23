load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
;load "/public/software/ncl_640/lib/ncarg/nclscripts/csm/contributed.ncl"
; Nino 1+2 (0-10S, 90W-80W), Nino 3 (5N-5S, 150W-90W), 
; Nino 3.4/ONI (5N-5S, 170W-120W) and Nino 4 (5N-5S, 160E-150W)
; TNI is Nino 1+2 - Nino 4
begin
  clatS  = -90;-5.0;85*1.6 ; or 85-95?!
  clatN  = -30;5.0; 95*1.6
  clonL  =  0;160.0;190.0;
  clonR  = 360;280.0;240.0;

  elatS  = -10
  elatN  = 0
  elonL  = 270
  elonR  = 280

SLAT   = -10;65;-30;
NLAT   =  10;90; 30;
WLON   = 190;  0;120;
ELON   = 240;360;280;  

  ymstrt= 187001;190001;187102  
  ymlast= 200412;201412;200812

;ycstrt = 187001 ;climate
;yclast = 201212 ;climate
;yms  = 194001
;yml  = 201612
;yms1 = 194001
;yml1 = 197612
;yms2 = 197701
;yml2 = 201612

;f= new(36,"file" )
;printVarSummary(f)
;f0 = addfile("air.mon.mean.vc2.1851-2014.nc", "r");data2/home/liufeiclim2014/xing
;f0 = addfile("/data/dswrf.sfc.mon.mean.vc2.1851-2014.nc", "r")

;f0 = addfile("/data/i20.anom.1871-2008.nc", "r")
;1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10.
;f0 = addfile("/data/air.sfc.mon.mean.vc2.1851-2014.nc", "r")
f0 = addfile("/data2/home/liufeiclim2014/xing/HadISST_1870-2016.nc","r");HadSST2_1x1.nc", "r");sst.ersst.1854-2018.v5.nc","r");sst.mnmean.v4.1854-2016.nc", "r");/data2/home/liufeiclim2014/xing/
;f0 = addfile("/data/dswrf.sfc.mon.mean.vc2.1851-2014.nc", "r")
;f0 = addfile("/data/air.mon.mean.vc2.1851-2014.nc", "r")
;f0 = addfile("/data/hgt.mon.mean.vc2.1851-2014.nc", "r")
;diri = "/data2/home/liufeiclim2014/xing/";"/nuist/p/work/feil/xc/" ;"/work/";
;f = addfile(diri+"HadSST2_1x1.nc", "r");temp.soda.1871-2008.TP.nc
;sst_obs = f->sst(:,:,:);{SLAT:NLAT},{WLON:ELON}
;printVarSummary(sst_obs)

TIME0=f0->time
yyyymm0=cd_calendar(TIME0, -1)
ymstrt0 =ind(yyyymm0.eq.ymstrt)
ymlast0 = ind(yyyymm0.eq.ymlast)
print(yyyymm0(0:100))
LEV = 20 ;500hPa:10 200hPa:16 850hPa:3 50hPa:20 300hPa:14
sst_obs = f0->sst(ymstrt0:ymlast0,:,:);{SLAT:NLAT},{WLON:ELON});{clatS:clatN},{clonL:clonR}{SLAT:NLAT},:);
printVarSummary(sst_obs)
sst1_obs = lonPivot(sst_obs,0.5)
printVarSummary(sst1_obs)
sst2_obs = sst1_obs(:,{SLAT:NLAT},{WLON:ELON})
;yyyymm1=cd_calendar(TIME0, -1)
;ymstrt1 =ind(yyyymm1.eq.ycstrt) ;climate
;ymlast1 = ind(yyyymm1.eq.yclast)

;     plev = f0->level                                      ; q&plev                
;    
;     ptop = min(plev)    ; 300==>psfc at each grid point
;     psfc = sst_obs(:,0,:,:)
;     psfc = 1008.0                   
;                         ; layer thicknesses
;     dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0) ;dp也需要是随经纬度变化的全球变量？！
;  vopt  = 0     ; vertically weighted average
;sst2_obs = wgt_vertical_n(sst_obs,dp,vopt,1)  ;整层
;sst2_obs = sst_obs
;;-------去掉神奇的-1000（不是missing value）-----
;dimo = dimsizes(sst1_obs)
;sstone = ndtooned(sst1_obs)
;Y = ind(sstone.lt.-10);.eq.-1000)
;sstone(Y) = sst1_obs@_FillValue
;sst2_obs = reshape(sstone, dimo)
;copy_VarMeta(sst1_obs, sst2_obs)
;printVarSummary(sst2_obs)
;sst2_obs!0 = "time"
;sst2_obs!1 = "lat"
;sst2_obs!2 = "lon"
;printVarSummary(sst2_obs)
;************************************************
; use landsea_mask function to mask out land and then ocean data
;************************************************
;
;     b = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
;
;     lsm_sst= landsea_mask(b->LSMASK,sst_obs&lat,sst_obs&lon)
;
;     sst_land = mask(sst_obs,lsm_sst.eq.0, False)      ; mask out ocean points
;     sst_ocean = mask(sst_obs,lsm_sst.eq.1, False)   ; mask out ocean points
;
;     copy_VarMeta(sst_obs,sst_land)        ;建议把这两句加上，这样新的land 和 ocean都有
;     copy_VarMeta(sst_obs,sst_ocean)     ;了和原始数据 T一样的meta信息，并且海洋(或陆地)设为了却测值-999
;
;
;printVarSummary(sst_land)
;printVarSummary(sst_ocean)
;
;scl=clmMonTLL(sst_land);(ymstrt1:ymlast1,:,:)
;sstla = calcMonAnomTLL(sst_land, scl)
;printVarSummary(sstla)
;sst2_obs = sst1_obs;(:,{SLAT:NLAT},{WLON:ELON})
sc=clmMonTLL(sst2_obs);(ymstrt1:ymlast1,:,:)
sa = calcMonAnomTLL(sst2_obs, sc)
;print(sa)
printVarSummary(sa)
;level = f->depth
;t= dtrend_msg_n(t&time,t,False,False,0)

;************************************************
; Interpolate to specific [constant] TEMP levels                找特定温度！
; The default returned vertical coordinate is Z_T but change to 'tlev'
;************************************************
;  tlev    = (/ 20./)    ; same units [here, C] as TEMP
;  tlev@units = t@units
;  tlev!0     = "tlev"
;
;  depth1   = conform(t, level, 1) ;Expands an array or scalar so that it conforms to the shape of the given variable.
;  copy_VarCoords(t, depth1)
;  printVarSummary(depth1)
;  
;  isot    = int2p_n_Wrap(t,depth1, tlev, 0, 1)
;
; printVarSummary(isot) ; [time | 1] x [tlev | 2] x [nlat | 384] x [nlon | 320]  
; printMinMax(isot, True) 
;
; i20 = isot(:,0,:,:)
; i20@lat = f->lat         ; used in graphics
; i20@lon = f->lon
; printVarSummary(i20)
; ;print(i20(0:5,:,:))
; ;printVarSummary(data)
;
;tc = clmMonTLL(i20)
;ta = calcMonAnomTLL(i20, tc)
;
;
;t1 = wgt_areaave_Wrap(ta, 1.0, 1.0, 1)
;printVarSummary(t1)
;t2 = dim_avg_n_Wrap(t1, 1)

;printVarSummary(sa)
;nino4= wgt_areaave_Wrap(sa,1.0, 1.0, 1);(:,{SLAT:NLAT},{WLON:ELON})(:,{clatS:clatN},{clonL:clonR})
;nino4o= wgt_areaave_Wrap(sst_obs(:,{clatS:clatN},{clonL:clonR}),1.0, 1.0, 1)
;
;nino12= wgt_areaave_Wrap(soa(:,{elatS:elatN},{elonL:elonR}),1.0, 1.0, 1)
;nino12o= wgt_areaave_Wrap(sst2_obs(:,{elatS:elatN},{elonL:elonR}),1.0, 1.0, 1)

lat = sa&latitude
rad    = 4.0*atan(1.0)/180.0
clat   = cos(lat*rad)
;clat   = cos(lat*pi/180)
;GLA = wgt_areaave_Wrap(sstla,clat, 1.0, 1)
;GL = wgt_areaave_Wrap(sst_land,clat, 1.0, 1)
GOA = wgt_areaave_Wrap(sa,clat, 1.0, 1)
;GO = wgt_areaave_Wrap(sst_obs,clat, 1.0, 1)
;;printVarSummary(sst_obs)
;printVarSummary(sst1_obs)
;printVarSummary(GMST1)
;TNIo = nino12o
;TNIo = nino12o-nino4o
;TNI = nino12
;TNI = nino12-nino4
;TNIo = dim_standardize_n_Wrap(TNIo, 0, 0)
;GMST = dim_standardize_n_Wrap(GMST1, 0, 0)
;asciiwrite("obs.air.verticalwgt.anom.unnorm.1870-2014.txt", GMST1)
;asciiwrite("obs.air.verticalwgt.org.unnorm.1870-2014.txt", GMST)
;asciiwrite("/work/obs.air.landair.anom.unnorm.1870-2014.txt", GLA)
;asciiwrite("/work/obs.air.landair.org.unnorm.1870-2014.txt", GL)
;asciiwrite("/work/obs.dsw.sfc.trop.anom.unnorm.1870-2014.txt", GOA)
;asciiwrite("/work/obs.dsw.sfc.trop.org.unnorm.1870-2014.txt", GO)
;asciiwrite("/work/obs.i20.160-280.anom.unnorm.1872-2008.txt", nino4o)
;asciiwrite("obs.temp100-200.nino34.anom.unnorm.187102-2008.txt", nino4)
;asciiwrite("obs.temp0-300.120-280.anom.unnorm.187102-2008.txt", t2)
;asciiwrite("obs.i20.nino34.anom.unnorm.187102-2008.txt", t1)
;asciiwrite("/work/cmip_vol/obs.50ta.TP.anom.unnorm.1870-2014.txt", nino4)
;asciiwrite("/work/cmip_vol/obs.50hgt.polarvortex.anom.unnorm.1870-2014.txt", GOA)
;asciiwrite("/work/cmip_vol/obs.air.sfc.SH.anom.unnorm.1870-2014.txt", GOA)
;asciiwrite("/data2/home/liufeiclim2014/xing/cmip_monsoon/obs.HadSST2.nino34.anomaly.1900-2004.txt", GOA)
asciiwrite("/data2/home/liufeiclim2014/xing/cmip_vol/obs.HadSST2.nino34(10).anomaly.1870-2004.txt", GOA)
;asciiwrite("/work/cmip_vol/obs.HadISST1.TP.anomaly.unnorm.1870-2016.txt", nino4)

;asciiwrite("/work/obs.ersst.nino34.anomaly.unnorm.1870-2016.txt", nino4)
;asciiwrite("/work/obs.ersst.nino34.org.unnorm.1870-2016.txt", nino4o)
;asciiwrite("/work/obs.HadISST.TNI.anomaly.unnorm.1870-2016.txt", TNI)
;asciiwrite("/work/obs.HadISST.TNI.org.unnorm.1870-2016.txt", TNIo)
end