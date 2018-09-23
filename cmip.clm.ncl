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
;diri = "/nuist/p/work/feil/xc/"
diri = "/data2/home/liufeiclim2014/xing/cmip_vol/"
diri2 = "/data2/data/cmip5/historical_regrid/atmos/mon/"
diri3 = "/data2/data/cmip5/historical/atmos/mon/"
diri0 = "/data2/home/liufeiclim2014/xing/"
;name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\
;       "CESM1-BGC","CESM1-CAM5-1-FV2", "CESM1-FASTCHEM","CESM1-WACCM","CNRM-CM5-2",  \; va
;      "CNRM-CM5","CSIRO-Mk3-6-0",  \; va
;        "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
;        "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES","MIROC5","MIROC-ESM-CHEM",\
;         "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\;"bcc-csm1-1-m","bcc-csm1-1",
name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\;
            "CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM",\;
            "CNRM-CM5-2",  "CNRM-CM5","CSIRO-Mk3-6-0",  \
             "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
             "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES",\
             "MIROC5","MIROC-ESM-CHEM",\
              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);;all
;name_model = (/"ACCESS1-0","ACCESS1-3","CanESM2","CCSM4",\3,4,7,8,9,10,20,
;            "CNRM-CM5-2","CNRM-CM5","CSIRO-Mk3-6-0", \;  
;             "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC",\;
;             "HadGEM2-ES",\;"HadGEM2-AO","GISS-E2-H","HadCM3","HadGEM2-CC","MIROC-ESM","MRI-CGCM3",
;              "MIROC5","MIROC-ESM-CHEM",\
;              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P",\
;              "MRI-ESM1","NorESM1-M","NorESM1-ME"/); rsdscs
NLAT =  30;40;70;5;5;
SLAT = -30;-40;-70;-5;-5;
WLON = 0;100; 190 ;
ELON = 360;300;240 ;
  ymstrt= 187001  
  ymlast= 200412
  yrstrt= 1870  
  yrlast= 2004
nm = dimsizes(name_model)
nyr = yrlast-yrstrt+1

LEV = 0
LEVEL = (/1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10/)

sign = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)" ,"(l)"/)

volyear_tp = (/1883,1902,1963,1982,1991/);1896(sh) (1968),1886,1968,

volmon_tp = (/188301,190201,196301,198201,199101/);1896(sh) (1968),196801,188601,
volyear = volyear_tp;ELyear;
volmon = volmon_tp;EL;+100
nv = dimsizes(volyear)

n = nv

MO = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
SEASON = (/"DJF","MAM", "JJA","SON" /)


files = systemfunc(" ls -1 "+diri2+"tas/"+name_model(0)+".*.nc")
;files = systemfunc(" ls -1 "+diri3+"zos_Omon_"+name_model(0)+"_historical*.nc")
f   = addfiles(files,"r") 
TIME=f[:]->time
yyyymm=cd_calendar(TIME, -1)
ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)
print(ymstrt1)
print(ymlast1)
x0  = f[:]->tas(ymstrt1(0):ymlast1(0),{SLAT:NLAT},{WLON:ELON});:,:)
nlat = dimsizes(x0&lat)
nlon = dimsizes(x0&lon)
ndim = dimsizes(x0)
delete([/f,files,TIME,yyyymm,ymstrt1,ymlast1/]);

f0 = addfile(diri0+"prate.mon.mean.vc2.1851-2014.nc","r")
TIME0=f0->time
yyyymm0=cd_calendar(TIME0, -1)
ymstrt0 =ind(yyyymm0.eq.ymstrt)
ymlast0 = ind(yyyymm0.eq.ymlast)
pr_obs = f0->prate(ymstrt0:ymlast0,{SLAT:NLAT},{WLON:ELON})
pr_obs = pr_obs*24*3600 ;prate-->precip
pr_obs= dtrend_msg_n(pr_obs&time,pr_obs,False,False,0)
;--------------resolution--------------
lat0 = pr_obs&lat
lon0 = pr_obs&lon 
klon0 = ispan(0,nlon-1,1)*2.5+WLON  
klat0 = ispan(0,nlat-1,1)*2.5+SLAT
pp= linint2_Wrap(lon0, lat0, pr_obs, True, klon0, klat0, 0) 
    pp&LAT@units="degrees_north"
    pp&LON@units="degrees_east"
    
prc_obs = clmMonTLL(pp)
printVarSummary(prc_obs)
prsea_obs = month_to_seasonN(prc_obs, SEASON)
printVarSummary(prsea_obs)


;xa =  new((/nm,ndim(0),17,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
xc =   new((/nm,  12,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) ;数组占内存太大就会读取不了文件啊！！！
xc_sea =   new((/nm,4,  1,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
dif_sea =   new((/nm,4,  1,ndim(1),ndim(2)/), typeof(x0), x0@_FillValue) 
;--------------------------
do ik =0,nm-1
;  if (ik.eq.13)then 
;    files = systemfunc(" ls -1 "+diri3+"psl/"+name_model(ik)+"/r1i1p1/*.nc")
;  else
;  if (ik.eq.12)then 
;    files = systemfunc(" ls -1 "+diri3+"hus/"+name_model(ik)+"/r1i1p1/*.nc")
;  else
 ; if (ik.eq.11)then 
 ;   files = systemfunc(" ls -1 "+diri3+"va/"+name_model(ik)+"/r1i1p1/*.nc")
 ; else
;files = systemfunc(" ls -1 "+diri+"rsds/rsds_Amon_"+name_model(ik)+"*.nc")
files = systemfunc(" ls -1 "+diri2+"pr/"+name_model(ik)+".pr.historical.nc");.mon
;files = systemfunc(" ls -1 "+diri2+"hgt/"+name_model(ik)+".hgt.historical.mon.nc") ;zg
;files = systemfunc(" ls -1 "+diri3+"ta/"+name_model(ik)+"/r1i1p1/*.nc")

;files = systemfunc(" ls -1 "+diri+"hus/hus_*"+name_model(ik)+".nc")
;files = systemfunc(" ls -1 "+diri2+"slp/"+name_model(ik)+".slp.historical.mon.nc");psl
;files = systemfunc(" ls -1 "+diri2+"tas/"+name_model(ik)+".tas.historical.mon.nc")
;files = systemfunc(" ls -1 "+diri2+"tos/"+name_model(ik)+".tos.historical.mon.nc")
;files = systemfunc(" ls -1 "+diri2+"ua/"+name_model(ik)+".ua.historical.nc")
;files = systemfunc(" ls -1 "+diri2+"wap/"+name_model(ik)+"*.nc")

;print(files)
;files = systemfunc(" ls -1 "+diri+"va/va_regrid"+name_model(ik)+".nc")
;end if 
;files = systemfunc(" ls -1 "+diri3+"zos_Omon_"+name_model(ik)+"_historical*.nc")

;(0) 100000(1) 92500(2) 85000(3) 70000(4) 60000(5) 50000(6) 40000(7) 30000(8) 25000(9) 20000(10)  15000(11)  10000(12)  7000(13)  5000(14)  3000(15)  2000(16)  1000
LEV = 13
f   = addfiles(files,"r")
TIME=f[:]->time
yyyymm=cd_calendar(TIME, -1)
ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)
print(ymstrt1)
print(ymlast1)
x  = f[:]->pr(ymstrt1(0):ymlast1(0),{SLAT:NLAT},{WLON:ELON});LEV,:,:);ymstrt1(0)+1619;ymlast1(0)

;--------------resolution--------------
lat = x&lat
lon = x&lon 
klon = ispan(0,nlon-1,1)*2.5+WLON  
klat = ispan(0,nlat-1,1)*2.5+SLAT
xx= linint2_Wrap(lon, lat, x, True, klon, klat, 0) 
    xx&LAT@units="degrees_north"
    xx&LON@units="degrees_east"

;xx  = x 
xx = xx*24*3600 ;prate-->precip
;xx = xx/100. ;slp Pa-->hPa  (name:psl)
delete([/f,files,TIME,x,yyyymm,ymstrt1,ymlast1,lat,lon/]);

xx= dtrend_msg_n(xx&time,xx,False,False,0)

;ndim = dimsizes(xx)
;sc1 = clmMonTLL(sst1)
;sc2 = clmMonTLL(sst2)
xc(ik,:,:,:)  = clmMonTLL(xx)
xc_sea(ik,:,:,:,:) = month_to_seasonN(xc(ik,:,:,:), SEASON)
dif_sea(ik,:,:,:,:) = xc_sea(ik,:,:,:,:) 
dif_sea(ik,:,:,:,:) = xc_sea(ik,:,:,:,:) - prsea_obs
;xa  = calcMonAnomTLL(xx, xc);(ik,:,:,:,:)
end do 
xc!0 = "model"
xc_sea!0 = "model"
printVarSummary(xc)
printVarSummary(xc_sea)
printVarSummary(dif_sea)

i1 = new((/nm/), float)
i2 = new((/nm/), float)
i3 = new((/nm/), float)
do im = 0,nm-1
i1(im) = wgt_areaave_Wrap(dif_sea(im,2,0,{-5:5},{240:280}),1.0, 1.0, 1);-5S-5N 120W-80W
i2(im) = wgt_areaave_Wrap(dif_sea(im,2,0,{-10:0},{120:180}),1.0, 1.0, 1);-10S-0 120E-180E
i3(im) = wgt_areaave_Wrap(dif_sea(im,2,0,{0:10},{240:280}),1.0, 1.0, 1);0-10N 120E-180E
end do 
printVarSummary(i3)
asciiwrite("CMIP5.pr.clm.jja.dif.ct.txt", i1)
asciiwrite("CMIP5.pr.clm.jja.dif.sw.txt", i2)
asciiwrite("CMIP5.pr.clm.jja.dif.nw.txt", i3)

;--------------plot-------------------
;wks = gsn_open_wks("ps", diri+"cmip.model.clm.pr.dif")
;;wks = gsn_open_wks("ps", diri+"cmip.TP.volcano.ssta1.lon-time"+name_model(ik));."+volyear(iv)
;gsn_define_colormap(wks,"precip_11lev");"BlueDarkRed18""BlueWhiteOrangeRed"
;
;res2= True
;res2@gsnDraw = False
;res2@gsnFrame = False
;res2@gsnLeftString = ""
;res2@gsnRightString = ""
;;res2@mpProjection         = "Robinson" ;"Mollweide"
;sdres = res2
;vcres = res2  
;  opt = res2
;res3  = res2
;
;res2@cnLinesOn     = False
;res2@gsnAddCyclic         =  True;False ;        ; data not cyclic
;  res2@mpMaxLatF                   = NLAT         ; choose subregion           
;  res2@mpMinLatF                   = SLAT;30 ;
;  res2@mpMaxLonF                   = ELON
;  res2@mpMinLonF                   = WLON
;res2@mpCenterLonF = 180.;-90;
;;res2@tmXMajorGrid                = True
;;res2@tmYMajorGrid                = True
;;res2@tmXMajorGridThicknessF      = 1.0
;;res2@tmYMajorGridThicknessF      = 1.0
;;res2@tmXMajorGridLineDashPattern = 2
;;res2@tmYMajorGridLineDashPattern = 2
;
;
;sdres@cnFillOn = False
;sdres@cnLinesOn = False
;sdres@cnLineLabelsOn        = False; True   
;res2@cnFillOn             = True         ; turn on color fill
;res2@cnInfoLabelOn            = False
;  res2@cnLineLabelsOn           = False ;True
; 
;   ;     sdres@cnLevelSelectionMode = "ManualLevels"  
;   ;     sdres@cnMinLevelValF        = 90
;   ;   sdres@cnLevelSpacingF      =  5         ;间隔5个画
;   ;     sdres@cnLevels = (/80,90,95,99.9/)
;sdres@cnSmoothingOn = True
;sdres@cnSmoothingTensionF = 1
;sdres@cnInfoLabelOn = False
;sdres@gsnLeftString = " "
;
;res2@tiMainString = " "
;
; ; res2@mpPerimOn            = False     ;Turns off the box perimeter.
;  res2@mpGridLatSpacingF    =  180.                ; change latitude  line spacing
;  res2@mpGridLonSpacingF    = 360.               ; change longitude line spacing
; ;res2@mpGridLineColor      = "transparent"      ; trick ncl into drawing perimeter
; ; res2@mpGridAndLimbOn      = True               ; turn on lat/lon lines
;  res2@mpFillOn          = False;True ;
;
;  res2@lbLabelBarOn        = False;True ;            ; turn off individual cb's
;res2@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
;res2@cnMinLevelValF        = 5
;res2@cnMaxLevelValF        = 15
;res2@cnLevelSpacingF       = 0.5
;
;res3@cnLinesOn = True
;res3@cnFillOn  = False
;res3@cnLineLabelsOn = False;True ;
;res3@cnLineLabelBackgroundColor = -1
;res3@cnInfoLabelOn            = False
;;res3@cnLineColor = "navajowhite"
;;res3@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
;;res3@cnMinLevelValF        = -1.0;280-273.16
;;res3@cnMaxLevelValF        =  1.0 ; 310-273.16
;res3@cnLevelSpacingF       =  2.5
;;res3@cnLineLabelInterval   =  1
;
;res2@gsnLeftStringFontHeightF = 0.03
;res2@gsnRightStringFontHeightF = 0.02
;;res2@tmXBLabelFontHeightF = 0.014
;;res2@tmYBLabelFontHeightF = 0.014
;do isea = 0,3
;plot = new(nm+1, graphic)
;res2@gsnLeftString = "OBS"
;plot(0)  = gsn_csm_contour_map(wks, prsea_obs(isea,0,:,:), res2) ;0,iyr+1
;
;do im = 0,nm-1
;res2@gsnLeftString = name_model(im)
;plot(im+1)  = gsn_csm_contour_map(wks, xc_sea(im,isea,0,:,:), res2) ;0,iyr+1
;plot0 = gsn_csm_contour(wks, dif_sea(im,isea,0,:,:), res3);svolavg1
;plot0 = ColorNegDashZeroPosContour(plot0,"lightskyblue1","transparent","pink1")
;overlay(plot(im+1), plot0)
;end do 
;
;
;resk   = True
;resk@txString = "climate precipitation in "+SEASON(isea);"composite erupted events";ENSO removed   hus(g/kg)
;
;resk@txFontHeightF = 0.02
;resk@gsnPanelLabelBar   = True                 ; add common colorbar
;
;resk@gsnMaximize = True
;resk@gsnPaperOrientation = "portrait"   ; force portrait
;;resk@lbLabelFontHeightF = 0.008
;;resk@lbBoxEndCapStyle = "TriangleBothEnds"
;gsn_panel(wks, plot, (/6,5/), resk)
;;delete([/res/])
;;print(iv)
;end do 

end