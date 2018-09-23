load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
begin

name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\;
            "CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM",\;
            "CNRM-CM5-2",  "CNRM-CM5","CSIRO-Mk3-6-0",  \
             "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
             "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES",\
             "MIROC5","MIROC-ESM-CHEM",\
              "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);
nm = dimsizes(name_model)

diri = "/work/cmip_vol/"
ymstrt = 1870
ymlast2 = 2004
ny2 = ymlast2-ymstrt+1

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

tt_dt = dtrend(tt_a_obs, True) ;rightmost dimension
nino3_dt = dtrend(nino3_a_obs, True)
g_dt =     dtrend(gmst_a_obs, True)   ;air,True)   
air_dt   =dtrend(air_a, True)
;--------------------------------------------------

mxlag = 30
x = nino3_dt;g_dt ;
y = air_dt;g_dt ; 
x_Lead_y = new((/nm,mxlag+1/),float)
y_Lead_x = new((/nm,mxlag+1/),float)
do im = 0,nm-1
x_Lead_y(im,:) = esccr(x(im,:),y(im,:),mxlag)
y_Lead_x(im,:) = esccr(y(im,:),x(im,:),mxlag) 
end do
printVarSummary(x_Lead_y)
ccr = new ((/nm, 2*mxlag+1/), float)    
     ccr(:,0:mxlag-1) = y_Lead_x(:,1:mxlag:-1)  ; "negative lag", -1 reverses order
     ccr(:,mxlag:)    = x_Lead_y(:,0:mxlag)     ; "positive lag"
x_axis   = ispan(-mxlag,mxlag,1)                     ; define x axis
printVarSummary(ccr)
printVarSummary(x_axis)
cr = dim_avg_n_Wrap(ccr, 0)
printVarSummary(cr)
;--------plot-------------------------
wks = gsn_open_wks("ps", "/work/cmip_vol/cmip.leadlag.gmst&nino34.avg")

;do im=0,nm-1
res = True
;res@gsnFrame = False
;res@gsnDraw  = False
res@vpWidthF    = 1.0        ; Make plots wider
res@vpHeightF   = 0.3
res@gsnMaximize = True
res@gsnPaperOrientation = "portrait"   ; force portrait

res@gsnLeftString = "GMST(air) lead"
res@gsnRightString = "nino34 lead ";"GMSST(ocean) lead";
res@tiXAxisString   = "LAG(month)"
res@trXMinF            = -mxlag+1
res@trXMaxF            = mxlag-1
res@trYMinF            = -0.4
res@trYMaxF            = 1.0
;res@gsnYRefLine       = (/-0.13,0.13/)
res@gsnXRefLine       = 0
res@gsnYRefLine       = 0
res@xyLineThicknessF = 2
res@tiMainString = "CMIP5 ensemble mean" ;name_model(im)
  plot = gsn_csm_xy(wks,x_axis,cr,res)  
;end do 


end 