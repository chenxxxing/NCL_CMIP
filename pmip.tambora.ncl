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
;Tambora 1815.4
begin
diri = "/data2/home/liufeiclim2014/xing/cmip_vol/"
diri2 = "/data/CMIP5-LME/cmip5-data/"
diri0 = "/data2/home/liufeiclim2014/xing/"
name_model = (/"bcc-csm1-1","CCSM4", "CSIRO-Mk3L-1-2"  ,"FGOALS-gl","GISS-E2-R","HadCM3",\
 "IPSL-CM5A-LR","MIROC-ESM","MPI-ESM-P"/);
NLAT = 90;
SLAT = -90;
WLON = 0;
ELON = 360;
  ymstrt= 180901  ;clmate state select which section?
  ymlast= 182012
  yrstrt= 1809  
  yrlast= 1820
  ycstrt= 180101
  yclast= 180812
nm = dimsizes(name_model)
nyr = yrlast-yrstrt+1
nyr1 = 1815-yrstrt
nym1 = nyr1*12+3 ;erupted month
nym = nyr*12

LEV = 3
LEVEL = (/1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10/)

sign = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)" ,"(l)"/)

MO = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
SEASON = (/"DJF","MAM", "JJA","SON" /)


gmst = new((/nm,nyr*12/),float)
wks = gsn_open_wks("ps", diri+"pmip.tambora.pr.eachmodel.season");"+name_model(ik));.

gsn_define_colormap(wks,"precip_diff_12lev");"BlueDarkRed18""BrownBlue12""precip_diff_12lev" "temp_19lev" 

;gsn_define_colormap(wks, "BlueWhiteOrangeRed")
plot = new(8, graphic)

;--------------------------
do ik =0,nm-1
if (name_model(ik).eq."IPSL-CM5A-LR") then
  files =  systemfunc(" ls -1 "+diri2+"pr_Amon_IPSL-CM5A-LR_past1000_r1i1p1_165001-184912.nc")
else if(name_model(ik).eq."GISS-E2-R") then
  files =  systemfunc(" ls -1 "+diri2+"pr_Amon_GISS-E2-R_past1000_r1i1p121_180101-185012.nc")
else 
;files = systemfunc(" ls -1 "+diri+"rsds/rsds_Amon_"+name_model(ik)+"*.nc")
files = systemfunc(" ls -1 "+diri2+"pr_Amon_"+name_model(ik)+"_past1000*");.mon
;files = systemfunc(" ls -1 "+diri2+"tos/"+name_model(ik)+".tos.historical.mon.nc")
;(0) 100000(1) 92500(2) 85000(3) 70000(4) 60000(5) 50000(6) 40000(7) 30000(8) 25000(9) 20000(10)  15000(11)  10000(12)  7000(13)  5000(14)  3000(15)  2000(16)  1000
;files = systemfunc(" ls -1 "+diri3+"zos_Omon_"+name_model(ik)+"_historical*.nc")
end if 
end if
print(files)

f   = addfiles(files,"r") 
TIME=f[:]->time
yyyymm=cd_calendar(TIME, -1)
printMinMax(yyyymm, False)
ymstrt1 =ind(yyyymm.eq.ymstrt)
ymlast1 = ind(yyyymm.eq.ymlast)
print(ymstrt1)
print(ymlast1)
ycstrt1 =ind(yyyymm.eq.ycstrt)
yclast1 = ind(yyyymm.eq.yclast)
print(ycstrt1)
print(yclast1)
;if(name_model(ik).eq."bcc-csm1-1".or.name_model(ik).eq."GISS-E2-R") then
;x  = f[:]->pr(11508:11651,{SLAT:NLAT},{WLON:ELON});1000hPa
;xc  = f[:]->pr(11400:11507,{SLAT:NLAT},{WLON:ELON});1000hPa
;else
x  = f[:]->pr(ymstrt1(0):ymlast1(0),{SLAT:NLAT},{WLON:ELON});1000hPa
xc  = f[:]->pr(ycstrt1(0):yclast1(0),{SLAT:NLAT},{WLON:ELON});1000hPa
;end if
delete([/f,files,TIME,yyyymm,ymstrt1,ymlast1/]);,x1,lat,lon
x = x*24*3600 ;prate-->precip
xc = xc*24*3600 ;prate-->precip

;x= dtrend_msg_n(x&time,x,False,False,0)

ndim = dimsizes(x)

xcc  = clmMonTLL(xc) 
xa  = calcMonAnomTLL(x, xcc)
printVarSummary(xa)
xsea = month_to_seasonN(xa, SEASON)
printVarSummary(xsea)
;-------------output---------------------------
fout_name = "pmip.pr.1809-1820."+name_model(ik)+".nc"
   system("/bin/rm -f "+diri+fout_name)   ; remove any pre-existing file
   fout = addfile(diri+fout_name ,"c")  ; open output netCDF file

   fAtt               = True            ; assign file attributes
       fAtt@title         = "precipitation anomaly in surface"
;       fAtt@source_file   =  "original-file.nc"
       fAtt@Conventions   = "None"
       fAtt@creation_date = systemfunc ("date")
       fileattdef( fout, fAtt )            ; copy file attributes
  ;    filedimdef(fout,"time",-1,True) 

    fout->pr = xa

  delete(fout)

;--------------plot-------------------
res2= True
res2@gsnDraw = False
res2@gsnFrame = False
res2@mpCenterLonF = 180.
res2@gsnLeftString = " "
sdres = res2


res2@cnLinesOn     = False
res2@gsnAddCyclic         =  True;False ;        ; data not cyclic
  res2@mpMaxLatF                   = NLAT         ; choose subregion           
  res2@mpMinLatF                   = SLAT
  res2@mpMaxLonF                   = ELON
  res2@mpMinLonF                   = WLON
res2@tmXMajorGrid                = True
res2@tmYMajorGrid                = True
res2@tmXMajorGridThicknessF      = 1.0
res2@tmYMajorGridThicknessF      = 1.0
res2@tmXMajorGridLineDashPattern = 2
res2@tmYMajorGridLineDashPattern = 2

vcres = res2  
  opt = res2

res2@cnFillOn             = True         ; turn on color fill
sdres@cnFillOn = False
sdres@cnLinesOn = False
sdres@cnLineLabelsOn        = False; True   
res2@cnInfoLabelOn            = False
  res2@cnLineLabelsOn           = False ;True
 
        sdres@cnLevelSelectionMode = "ManualLevels"  
        sdres@cnMinLevelValF        = 90
      sdres@cnLevelSpacingF      =  5         ;间隔5个画
        sdres@cnLevels = (/80,90,95,99.9/)
sdres@cnSmoothingOn = True
sdres@cnSmoothingTensionF = 1
sdres@cnInfoLabelOn = False
sdres@gsnLeftString = " "

res2@tiMainString = " "

;res2@mpProjection         = "Robinson" ;"Mollweide"
  res2@mpPerimOn            = False     ;Turns off the box perimeter.
  res2@mpGridLatSpacingF    =  90                ; change latitude  line spacing
  res2@mpGridLonSpacingF    = 360.               ; change longitude line spacing
 ;res2@mpGridLineColor      = "transparent"      ; trick ncl into drawing perimeter
  res2@mpGridAndLimbOn      = True               ; turn on lat/lon lines
  res2@mpFillOn          = False;True ;

  res2@lbLabelBarOn        = False;True ;            ; turn off individual cb's
res2@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
res2@cnMinLevelValF        = -4.0;280-273.16
res2@cnMaxLevelValF        =  4.0 ; 310-273.16
res2@cnLevelSpacingF       =  0.4
;res2@cnLevels = fspan(-2.0, 2.0, 20)
; res2@cnFillDrawOrder      = "Predraw"  ; draw contours first 可以用地图把陆地盖住
;plot(i+12*j) = gsn_csm_contour_map_ce(wks2,sepavg(j*12+i,:,:),res2)
;plot(2*i+j) = gsn_csm_contour_map_ce(wks2,probep(j*12+i,:,:),res2)



  opt@gsnShadeFillType = "pattern"      ; pattern fill
  opt@gsnShadeHigh = 17          ; 高于95的use pattern #17
  opt@gsnshadelow = -1       ;低于95为空白  
  opt@cnLineLabelsOn        = False   
  opt@gsnLeftString  = " "
  opt@cnInfoLabelOn  = False
;  ptest0 = gsn_csm_contour(wks2, prob0, sdres)
;ptest0 = gsn_contour_shade(ptest0,-999,95,opt) ;Pattern shade contours above 95.

vcres@vcMinDistanceF           = 0.017            ; thin vectors
;res2@cnFillPalette           = "gui_default"     ; set color map
vcres@vcRefAnnoOrthogonalPosF = -.535    ; move ref vector into plot
  vcres@vcRefAnnoFontHeightF    = 0.005
vcres@vcRefMagnitudeF          = 7.0             ; define vector ref mag
vcres@vcRefLengthF             = 0.1            ; define length of vec ref
vcres@vcRefAnnoOrthogonalPosF  = .1               ; move ref vector down
;vcres@vcRefLengthF             = 0.045            ; define length of vec ref
vcres@gsnAddCyclic         = False        ; data not cyclic 
vcres@gsnLeftString   = " "

ressf = True
ressf@gsnFrame = False
ressf@gsnDraw = False
ressf@cnLinesOn = True
ressf@cnFillOn  = False
ressf@cnLineLabelsOn = False
ressf@cnInfoLabelOn       = False
    ressf@gsnContourNegLineDashPattern = 1
 ;ressf@cnLineThicknessF = 1.5

res2@gsnLeftStringFontHeightF = 0.02
res2@gsnRightStringFontHeightF = 0.02
;res2@tmXBLabelFontHeightF = 0.014
;res2@tmYBLabelFontHeightF = 0.014
do j = 5,8 
do i = 0,1
  isea = i*2
  k = (j-5)*2+i
iyr = j-6
res2@gsnLeftString = sign(k)+SEASON(isea)+iyr;name_model(ik) ;"A";MO(j); LABEL(isea)+
res2@gsnRightString = "";"year("+iyr+")"
;resk@txString = "composite surface air "+name_model(ik)+" year("+iyr1+")";
;plot(isea)  = gsn_csm_contour_map_ce(wks, xvsea(ik,isea,iyr,:,:), res2) 
;ptest = gsn_csm_hov(wks, probv_sea(ik,isea,iyr,:,:), sdres)
plot(k)  = gsn_csm_contour_map(wks,xsea(isea,j,:,:) , res2) 
;ptest = gsn_csm_hov(wks, psea(isea,iyr,:,:), sdres)
;ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
;overlay(plot(isea), ptest)
end do 
end do 

resk = True
resk@txFontHeightF = 0.02
resk@gsnPanelLabelBar   = True                 ; add common colorbar
resk@gsnMaximize = True
resk@gsnPaperOrientation = "portrait"   ; force portrait
;resk@lbBoxEndCapStyle = "TriangleBothEnds"
resk@txString = "Tembora eruption in "+name_model(ik);

gsn_panel(wks, plot, (/4,2/), resk)
lat = xa&lat
rad    = 4.0*atan(1.0)/180.0
clat   = cos(lat*rad)
gmst(ik,:) =  wgt_areaave_Wrap(xa,clat, 1.0, 1)

delete(res2)
delete([/x,xa,xc,xcc,xsea,lat,clat/])
end do 
printVarSummary(gmst)
;---------------xy---------------------------
wks2 = gsn_open_wks("ps", diri+"pmip.tambora.pr.xy.eachmodel")
plot2 = new(1, graphic)
resxy = True
resxy@gsnFrame = False
resxy@gsnDraw = False
resxy@gsnYRefLine = 0
resxy@vpHeightF        = 0.6                 ; default is 0.6
resxy@vpWidthF         = 0.8; 0.7;               ; default is 0.6
resxy@tmXTOn = False
resxy@tmYROn = False
resxy@tiYAxisString = "GMST"
resxy@trYMinF             = -0.2
resxy@trYMaxF             =  0.1
resxy@gsnLeftString = "";(a) observation
resxy@gsnLeftStringOrthogonalPosF = -0.12
resxy@gsnRightString = ""
x_axis = ispan(0, nym-1, 1)
x_axis!0 = "time"
resxy@tmXBMode   = "Explicit"                                                
resxy@tmXBValues = (/0,35,71,107,143,179,215,251/) 
resxy@tmXBLabels = (/"year(-6)","year(-3)","year(0)","year(+3)","year(+6)","year(+9)","year(+12)","year(+15)"/)      
resxy@tmXBLabelFontHeightF = 0.01
colors = (/"lightsalmon1","LightBlue","lightgoldenrod3","limegreen","orchid3", "lightslategray",\
 "brown","blue2","red2","orangered","goldenrod1","hotpink","greenyellow","lemonchiffon"/)
resxy@xyLineColors = colors
resxy@xyDashPattern = 0
resxy@gsnXRefLine = nym1
resxy@gsnXRefLineDashPattern = 1
 plot2(0)  = gsn_csm_xy(wks2, x_axis, gmst, resxy)
;==============================================
;  attach plots
;=============================================
  lgres                    = True
  lgres@lgLineColors       = colors;(/"red","green" /) ; "blue"
  lgres@lgBoxLineDashPattern = 0
  lgres@lgMarkerSizeF  = 0.01
  lgres@lgItemType         = "Lines" ;"Markers";       ; show lines only (default)
  lgres@lgLabelFontHeightF = .09            ; legend label font thickness
  lgres@vpWidthF           = 0.16           ; width of legend (NDC)
  lgres@vpHeightF          = 0.20           ; height of legend (NDC)
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter
  lgres@lgMonoDashIndex    = True
  lgres@lgMarkerIndex = 16
 ; lgres@lgPerimFill = 0
 ; lgres@lgPerimFillColor = 0

  labels =  name_model
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

  annoid = gsn_add_annotation(plot2(0),lbid,amres)  ; add legend to plot
resxyk = True
resxyk@txString = "GMST Tambora"
resxyk@gsnMaximize = True
resxyk@gsnPaperOrientation = "portrait"   ; force portrait

gsn_panel(wks2, plot2, (/1,1/), resxyk)

;end do 
;end do 
;delete([/res2,x,xx,xc,xa,TIME,yyyymm,ymstrt1,ymlast1,files,f,xavol,xseavol,xavolc,xavol1,xseavol1/])
;print(ik)
;end do 
end