load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
; Nino 1+2 (0-10S, 90W-80W), Nino 3 (5N-5S, 150W-90W), 
; Nino 3.4/ONI (5N-5S, 170W-120W) and Nino 4 (5N-5S, 160E-150W)
begin
NLAT = 5
SLAT = -5
WLON = 120 ;210;should be fixed baced the regressed pattern
ELON = 280
MO = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
SEASON = (/"DJF","MAM", "JJA","SON" /)
sign = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)" ,"(l)"/)
nyear = (/"year(-5)","year(-4)","year(-3)","year(-2)","year(-1)","year(0)","year(+1)",\
         "year(+2)","year(+3)","year(+4)","year(+5)","year(+6)","year(+7)"/);"year(-6)",

;volmon = (/188308,190210,191301,196303,198203,199106/);
volmon = (/188301,190201,191301,196301,198201,199101/);
volyear = (/1883,1902,1963,1982,1991/);1913,
volyearF = (/1883.67,1902.83,1913.08,1963.25,1982.25,1991.5/)
volyear_c = (/"1883","1902","1913","1963","1982","1991"/);
;EL=(/1880,1884,1888,1895,1896,1899,1902,1904,1905,1911,1914,1918,1923,1925,1930,1939\
;    ,1940,1941,1951,1953, 1957,1963,1965,1968,1969,1972,1976,1977,1982,1986,1987,1991,\
;    1994,1997,2002/); 1876,,2004,2006;,2009

ELmon=(/188001,188801,189501,189601,189901,190201,191101,191801,192301,192501,193001,193901\
,194001,194101,195101,195301,195701,196801,196901,197201,197601,197701,198601,198701,199401,199701\
,200201/);187601,1884,1904,1905,1914,1963,1965,1982,1991,,2009,200601,200401


;EL = (/1876,1880,1884,1895,1896,1918,1982,1997/);EP
;EL = (/1904,1914,1940,1941,1957,1963,1965,1968,1969,1977,1987,1991,1994,2002,2004/);CP,2009
;EL = (/1888,1899,1902,1905,1911,1923,1925,1930,1939,1951,1953,1972,1976,1986,2006/);MIX 
;EL = (/1876,1880,1895,1896,1918,1997/);EP 1884,1982,
;EL = (/1940,1941,1957,1965,1968,1969,1977,1987,1994,2002,2004/);CP1904,1914,1963,1991,,2009
;EL = (/1888,1899,1902,1905,1911,1923,1925,1930,1939,1951,1953,1972,1976,1986,2006/);MIX 
ne = dimsizes(ELmon)
nv = dimsizes(volyear)

diri = "/data2/home/liufeiclim2014/xing/";"/nuist/p/work/feil/xc/" ;"/work/";
f = addfile(diri+"temp.soda.1871-2008.TP.nc", "r")
t = f->temp(:,:,{SLAT:NLAT},{WLON:ELON})
printVarSummary(t)
;----------------filter----------------------PAY ATTENTION! This is month data, not years!
tf = t
;uf = ua(ik,:,:,:)
;vf = va(ik,:,:,:)
ihp   = 1       ; low-pass filter: ihp = 0; high-pass ihp = 1; band-pass ihp = 2.
sigma = 1.0                           ; Lanczos sigma
nWgt  = 71 ;必须是奇数 表示权重系数的总个数，越大效果越好，但是两端损失越多
fca   = 1./(10*12) ;                        
fcb   = -999;1./(2*12);1./8. *12  ;8——50年滤波
wgt  = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
tf = wgt_runave_n (t, wgt, 0, 0 ) ;最后一个是指sst1第0维

tc = clmMonTLLL(tf)
ta = calcMonAnomTLLL(tf, tc)
t_p = dim_avg_n_Wrap(ta, 2)
ndim = dimsizes(t_p)
printVarSummary(t_p)
;---------------composite----------------------
ny1 = 5
ny2 = 5
nyr = ny1+ny2
nmon = ny1*12+ny2*12

tv = new((/nv,ny1*12+ny2*12,ndim(1),ndim(2)/), float,ta@_FillValue)
te = new((/ne,ny1*12+ny2*12,ndim(1),ndim(2)/), float,ta@_FillValue)

TIME=ta&time
yyyymm=cd_calendar(TIME, -1)
ymvol2 = new(nv, integer)
ymel2  = new(ne, integer)

do iv =0,nv-1
  ;  kk = (EL(i)-ny1-1871)*12-1
  ;  ll = (EL(i)+ny2-1871)*12-2
    ymvol2(iv) = ind(yyyymm.eq.volmon(iv));EL(iv));
    tv(iv,:,:,:) = t_p(ymvol2(iv)-ny1*12:ymvol2(iv)+ny2*12-1,:,:);(kk:ll,:,:)

end do 
do iv=0, ne-1
   ; kk = (volyear(i)-ny1-1871)*12-1
   ; ll = (volyear(i)+ny2-1871)*12-2
    ymel2(iv) = ind(yyyymm.eq.ELmon(iv));EL(iv));
    te(iv,:,:,:) = t_p(ymel2(iv)-ny1*12:ymel2(iv)+ny2*12-1,:,:);(kk:ll,:,:)
end do 
printVarSummary(tv)
tv!0 = "sample"
tv!1 = "time"
te!0 = "sample"
te!1 = "time"

tea = te 
tva = tv

;----------------------------------------
tv_clm = new((/nv,12,ndim(1),ndim(2)/),float,ta@_FillValue)
te_clm = new((/ne,12,ndim(1),ndim(2)/),float,ta@_FillValue)
do i=0,11
  tv_clm(:,i,:,:) = dim_avg_n_Wrap(tv(:,i::12,:,:), 1)
  te_clm(:,i,:,:) = dim_avg_n_Wrap(te(:,i::12,:,:), 1)

end do 
printVarSummary(te_clm)

tva  = new((/nv,nmon,ndim(1),ndim(2)/),float,ta@_FillValue)
tea  = new((/ne,nmon,ndim(1),ndim(2)/),float,ta@_FillValue)

do i = 0,11
do j = 0,nyr-1
tva(:,j*12+i,:,:)  = tv(:,i+j*12,:,:)-tv_clm(:,i,:,:)
tea(:,j*12+i,:,:)  = te(:,i+j*12,:,:)-te_clm(:,i,:,:)

end do 
end do 
copy_VarMeta(tv, tva)
copy_VarMeta(te, tea)

printVarSummary(tva)


ce = dim_avg_n_Wrap(tea, 0)
cv = dim_avg_n_Wrap(tva, 0)
printVarSummary(cv)

tva_sea = month_to_seasonN(tva(time|:,sample|:,depth|:,lon|:), SEASON)
tea_sea = month_to_seasonN(tea(time|:,sample|:,depth|:,lon|:), SEASON)
printVarSummary(tva_sea)
ce_sea = dim_avg_n_Wrap(tea_sea, 2)
cv_sea = dim_avg_n_Wrap(tva_sea, 2)
dif_sea = ce_sea
dif_sea = cv_sea-ce_sea
printVarSummary(dif_sea)

;-------------t test ------------注意!!!样本是异常场!!总体也得是!!!!
avesumsst = new((/nmon,ndim(1),ndim(2)/), float, ta@_FillValue)
varsumsst = new((/nmon,ndim(1),ndim(2)/), float, ta@_FillValue)
;avesumu =   new((/12*ny,undim(1),undim(2)/), float, ua@_FillValue)
;varsumu =   new((/12*ny,undim(1),undim(2)/), float, ua@_FillValue)
;avesumv =   new((/12*ny,undim(1),undim(2)/), float, va@_FillValue)
;varsumv =   new((/12*ny,undim(1),undim(2)/), float, va@_FillValue)

        
avesumsst = 0 ;既然是异常场的平均值,那就应该是零啊!!
avesumu = 0
avesumv = 0

varsst = dim_variance_n_Wrap(tva, 0)
;varu = dim_variance_n_Wrap(uavol, 0)
;varv = dim_variance_n_Wrap(vavol, 0)

siglvl=0.05

do j =0,(ny1+ny2)-1
do i =0,11
varsumsst(i+j*12,:,:) = dim_variance_n_Wrap(t_p(i::12,:,:), 0)
;varsumu(i+j*12,:,:) = dim_variance_n_Wrap(ua(i::12,:,:), 0)
;varsumv(i+j*12,:,:) = dim_variance_n_Wrap(va(i::12,:,:), 0)

end do 
end do
printVarSummary(varsumsst)   
                                                    ;所有sst不同月的方差

;sxep = ssumsf

ssum = 2016-1870
snp = nv
prob = cv
;probu = uvolavg
;probv = vvolavg

    iflag= True              ; population variance similar 

  prob=(1- ttest(cv,varsst,snp, avesumsst,varsumsst,ssum, iflag, False))*100
 ; probu=(1- ttest(uvolavg,varu,snp, avesumu,varsumu,ssum, iflag, False))*100
 ; probv=(1- ttest(vvolavg,varv,snp, avesumv,varsumv,ssum, iflag, False))*100

  copy_VarMeta(cv,prob)
; copy_VarMeta(uvolavg,probu)
; copy_VarMeta(vvolavg,probv)  

printVarSummary(prob)
prob1 = month_to_seasonN(prob, SEASON)

;--------------plot-------------------
;wks = gsn_open_wks("ps", diri+"obs.ocean-temp.profile");+ny1+ny2);(eruptmon)+volyear(iv)
wks = gsn_open_wks("ps", diri+"cmip_vol/obs.composite.ocean-temp");+volyear(iv));."
gsn_define_colormap(wks,"BlueDarkRed18");"MPL_bwr""MPL_BrBG" "BlueDarkRed18"
;gsn_define_colormap(wks, "BlueWhiteOrangeRed")
;do iv=0,nv-1

plot = new(3, graphic)
res                  = True                ; plot mods desired
res@gsnDraw           = False
res@gsnFrame          = False
res@gsnLeftString = ""
sdres = res
vcres = res  
  opt = res

res@vpWidthF    = 0.7        ; Make plots wider
res@vpHeightF   = 0.3
res@gsnLeftStringFontHeightF = 0.03
res@gsnRightStringFontHeightF = 0.03
;res@tmXMajorGrid                = True
;res@tmYMajorGrid                = True
;res@tmXMajorGridThicknessF      = 1.0
;res@tmYMajorGridThicknessF      = 1.0
;res@tmXMajorGridLineDashPattern = 2
;res@tmYMajorGridLineDashPattern = 2
 res@cnFillOn         = True                ; color on
 res@lbLabelStride    = 4                   ; every other label
 res@lbOrientation    = "Vertical"          ; vertical label bar
 ;res@pmLabelBarOrthogonalPosF = -0.04       ; move label bar closer to axis
 res@cnLinesOn        = False ; True;              ; turn off contour lines
; res@cnFillPalette    = "BlueDarkRed18"    ; set color map
res@gsnLeftString     = " "
res@gsnRightString    = " "
res@tmYROn            = False   ;右侧坐标轴不显示？
 res@tiMainString     = "  "  ; title
res@cnInfoLabelOn     = False
res@cnLineLabelsOn    = False
 res@vpXF             = 0.12                ; default is 0.2 (aspect ratio)
 res@vpYF             = 0.8                 ; default is 0.8
; res@vpHeightF        = 0.4                 ; default is 0.6
; res@vpWidthF         = 0.75                ; default is 0.6

res@gsnAddCyclic         = False        ; data not cyclic 
 ;res@vcGlyphStyle    = "CurlyVector"      ; turn on curly vectors
res@lbLabelBarOn      = True ;False
 res@cnLevelSelectionMode  = "ManualLevels"    ;修改范围需要设定selection mode！！
res@cnMinLevelValF        = -1.0;10;
res@cnMaxLevelValF        =  1.0; 35;
res@cnLevelSpacingF       =  0.10; 2;

;res@tmYLMode   = "Explicit"                                                
;aa = ymvol2(nv-1)-ny1*12
;bb = ymvol2(nv-1)+ny2*12-1 
;tt = sa&time
;res@tmYLValues = tt(aa:bb:6)
;res@tmYLLabels = (/"Jan(-1)","Jul(-1)","Jan(0)","Jul(0)","Jan(+1)","Jul(+1)","Jan(+2)","Jul(+2)" /) 
   

 res1 = True

 res1@gsnDraw = False
 res1@gsnFrame = False

res1@gsnLeftString = " "
res1@cnLinesOn = True
res1@cnFillOn  = False
res1@cnLineLabelsOn = False
res1@cnInfoLabelOn       = False
    res1@gsnContourNegLineDashPattern = 1
res1@cnLineColor  = "red"
;level@units = "hPa"
sdres@cnFillOn = False
sdres@cnLinesOn = False
sdres@cnLineLabelsOn        = False; True   
 
     ;   sdres@cnLevelSelectionMode = "ManualLevels"  
     ;   sdres@cnMinLevelValF        = 90
     ; sdres@cnLevelSpacingF      =  5         ;间隔5个画
     ;   sdres@cnLevels = (/80,90,95,99.9/)

sdres@cnSmoothingOn = True
sdres@cnSmoothingTensionF = 1
sdres@cnInfoLabelOn = False
sdres@gsnLeftString = ""
sdres@gsnRightString = ""
  opt@gsnShadeFillType = "pattern"      ; pattern fill
  opt@gsnShadeHigh = 17          ; 高于95的use pattern #17
  opt@gsnshadelow = -1       ;低于95为空白  
  opt@cnLineLabelsOn        = False   
  opt@gsnLeftString  = " "
  opt@cnInfoLabelOn  = False
;--------lon-depth--------------------------------
res@tiMainString = "";
resk   = True
;resk@gsnPanelLabelBar   = True                 ; add common colorbar
resk@gsnMaximize = True
resk@gsnPaperOrientation = "portrait"   ; force portrait

 res@trYReverse                      = True                  ; reverse y-axis
 res@tiYAxisString = "depth"
res@gsnLeftString = "DJF";MO(j);"+"+lastmon ;name_model(i) ; 
res@gsnRightString = "year(-1)" ;"year("+i+")"
plot(0)  = gsn_csm_contour(wks, cv_sea(0,5,:,:), res)
ptest = gsn_csm_hov(wks, prob1(0,5,:,:), sdres)
ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
;;sfplot  = gsn_csm_contour(wks,sf(i,0,:,:),ressf)
;overlay(plot(j),vcplot)
overlay(plot(0), ptest)

res@gsnLeftString = "JJA";MO(j);"+"+lastmon ;name_model(i) ; 
res@gsnRightString = "year(0)" ;"year("+i+")"
plot(1)  = gsn_csm_contour(wks, cv_sea(2,5,:,:), res)
ptest = gsn_csm_hov(wks, prob1(2,5,:,:), sdres)
ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
;;sfplot  = gsn_csm_contour(wks,sf(i,0,:,:),ressf)
;overlay(plot(j),vcplot)
overlay(plot(1), ptest)

;i=5
;do j=1,3
;kk = j
;;i1 = i-5
;res@gsnRightString = nyear(i) ;"year("+i+")"
;;res@gsnLeftString = " "
;lastmon = j+i*12
;res@gsnLeftString = SEASON(j);sign(j)+MO(j);"+"+lastmon ;name_model(i) ; 
;;res@txString = "annual cycle index regressed onto sst in "+name_model(i);difference with observation"
;
;;plot(j)  = gsn_csm_contour_map_ce(wks, sseavol1(iv,j,i,:,:), res) 
;plot(kk)  = gsn_csm_contour(wks, cv_sea(j,i,:,:), res) 
;ptest = gsn_csm_hov(wks, prob1(j,i,:,:), sdres)
;ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
;;;sfplot  = gsn_csm_contour(wks,sf(i,0,:,:),ressf)
;;overlay(plot(j),vcplot)
;overlay(plot(kk), ptest)
;;plot(0)  = gsn_csm_contour_ce(wks, p0_jja, res)
;;overlay(plot(0), ptest0)
;end do 
res@gsnLeftString = "DJF";name_model(ik) ;"A";MO(j); LABEL(isea)+
res@gsnRightString = "year(0)" ;"year("+i+")"
plot(2)  = gsn_csm_contour(wks, cv_sea(0,6,:,:), res) 
;  plot(0) = gsn_csm_contour_map_polar(wks, cv_sea(0,6,{30:90},:), res)    ;0,i+1 create the plot

ptest = gsn_csm_hov(wks, prob1(0,6,:,:), sdres)
ptest = gsn_contour_shade(ptest,-999,90,opt) ;Pattern shade contours above 95.
overlay(plot(2), ptest)

;res@txString = "ssh in "+volyear(iv)+" event "+nyear(i);
resk@txString ="(a) observation";+nyear(i);"surface down sw in "+volyear(iv)+" volcano event in "+nyear(i); "";
;(eruption mon) 
gsn_panel(wks,plot,(/3,1/),resk);(without removing nearby annual cycle)
delete(res)
;end do 

end
