load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
; Nino 1+2 (0-10S, 90W-80W), Nino 3 (5N-5S, 150W-90W), 
; Nino 3.4/ONI (5N-5S, 170W-120W) and Nino 4 (5N-5S, 160E-150W)
begin
diri = "/data2/home/liufeiclim2014/xing/cmip_vol/" ;"/data2/home/liufeiclim2014/xing/thetao/"
diri2 = "/data2/home/liufeiclim2014/xing/thetao/"
diri3 = "/data2/data/cmip5/historical/atmos/mon/hus/"
name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\
       "CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM","CNRM-CM5-2",  \
      "CNRM-CM5","CSIRO-Mk3-6-0",  \
        "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;
        "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES","MIROC5","MIROC-ESM-CHEM",\
         "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P", "MRI-ESM1","NorESM1-M","NorESM1-ME"/);
nm = dimsizes(name_model)
do i=0,nm-1
;system("cdo genbil,n32 "+diri2+"thetao_*"+name_model(i)+"_*.TP.nc remapweights.nc")
;system("cdo remap,n32,remapweights.nc "+diri2+"thetao_*"+name_model(i)+"_*.TP.nc "+diri2+"thetao_a.nc")
;system("cdo intlevel,5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,155,165,175,185,195,205,215,225,235,245,255,265,275,285,295,305"+diri2+"thetao_a.nc "+diri2+"thetao_10level_"+name_model(i)+".nc")
;system("cdo mermean,5,-5 "+diri2+"thetao_10level_"+name_model(i)+".nc "+diri2+"thetao_10level_mermean_"+name_model(i)+".nc")
;system("cdo fldmean,190,300,5,-5 "+diri2+"thetao_10level_"+name_model(i)+".nc "+diri2+"thetao_1.nc")
;system("cdo vertmean,105,115,125,135,145,155,165,175,185,195 "+diri2+"thetao_1.nc "+diri2+"thetao_"+name_model(i)+".txt")
system("cdo cat "+diri3+name_model(i)+"/r1i1p1/hus_*.nc "+diri+"hus_aa.nc")
;files = systemfunc("ls -1 "+diri3+name_model(i)+"/r1i1p1/va_*.nc")
;print(files)
;nfile = dimsizes(files)
;do j =0,nfile-1
system("cdo genbil,n32 "+diri+"hus_aa.nc "+diri+"remapweights.nc")
system("cdo remap,n32,"+diri+"remapweights.nc "+diri+"hus_aa.nc "+diri+"hus/hus_regrid_"+name_model(i)+".nc")


;system(" rm "+diri2+"thetao_a.nc" )
system(" rm "+diri+"remapweights.nc")
system(" rm "+diri+"hus_aa.nc" )
print(name_model(i))
end do 
end