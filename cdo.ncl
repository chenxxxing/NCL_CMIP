load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
diri = "/disk/ipcc/cmip5_data/historical/ocean/mon/" 
diri2 = "/home/binliu/temp/xingchen/"
name_model = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1-m","bcc-csm1-1","CanESM2","CCSM4",\;
       "CESM1-BGC",  "CESM1-CAM5-1-FV2","CESM1-FASTCHEM","CESM1-WACCM","CNRM-CM5-2",  \;
      "CSIRO-Mk3-6-0",   "GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H-CC",\;"CNRM-CM5",
        "GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-ES","MIROC5","MIROC-ESM-CHEM",\;
         "MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P","MRI-ESM1", "NorESM1-M","NorESM1-ME"/);

nm = dimsizes(name_model)
NLAT = 5;2;
SLAT = -5;-2;
WLON = 120;150; 
ELON = 285;155;;
  ymstrt= 187001  
  ymlast= 200412
  yrstrt= 1870  
  yrlast= 2004


do i=0,nm-1
system("cdo cat "+diri+"thetao/"+name_model(i)+"/r1i1p1/thetao_Omon_"+name_model(i)+"_*.nc "+diri2+"thetao_pp.nc")
system("cdo genbil,n32 "+diri2+"thetao_pp.nc "+diri2+"remapweights.nc")
system("cdo remap,n32,"+diri2+"remapweights.nc "+diri2+"thetao_pp.nc "+diri2+"thetao_a.nc")
system("cdo intlevel,5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,155,165,175,185,195,205,215,225,235,245,255,265,275,285,295,305 "+diri2+"thetao_a.nc "+diri2+"thetao_10level_1.nc")
;system("cdo mermean,5,-5 "+diri2+"thetao_10level_1.nc "+diri2+"thetao_10level_mermean_"+name_model(i)+".nc")
system("cdo sellonlatbox,120,280,5,-5 "+diri2+"thetao_10level_1.nc "+diri2+"thetao_10level_TP_"+name_model(i)+".nc")

system(" rm "+diri2+"thetao_pp.nc" )
system(" rm "+diri2+"thetao_a.nc" )
system(" rm "+diri2+"thetao_10level_1.nc" )
system(" rm "+diri2+"remapweights.nc" )

system("cdo cat "+diri+"uo/"+name_model(i)+"/r1i1p1/uo_Omon_"+name_model(i)+"_*.nc "+diri2+"uo_pp.nc")
system("cdo genbil,n32 "+diri2+"uo_pp.nc "+diri2+"remapweights.nc")
system("cdo remap,n32,"+diri2+"remapweights.nc "+diri2+"uo_pp.nc "+diri2+"uo_a.nc")
system("cdo intlevel,5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,155,165,175,185,195,205,215,225,235,245,255,265,275,285,295,305 "+diri2+"uo_a.nc "+diri2+"uo_10level_1.nc")
;system("cdo mermean,5,-5 "+diri2+"uo_10level_1.nc "+diri2+"uo_10level_mermean_"+name_model(i)+".nc")
system("cdo sellonlatbox,120,280,5,-5 "+diri2+"uo_10level_1.nc "+diri2+"uo_10level_TP_"+name_model(i)+".nc")
system(" rm "+diri2+"uo_pp.nc" )
system(" rm "+diri2+"uo_a.nc" )
system(" rm "+diri2+"uo_10level_1.nc" )
system(" rm "+diri2+"remapweights.nc" )

system("cdo cat "+diri+"vo/"+name_model(i)+"/r1i1p1/vo_Omon_"+name_model(i)+"_*.nc "+diri2+"vo_pp.nc")
system("cdo genbil,n32 "+diri2+"vo_pp.nc "+diri2+"remapweights.nc")
system("cdo remap,n32,"+diri2+"remapweights.nc "+diri2+"vo_pp.nc "+diri2+"vo_a.nc")
system("cdo intlevel,5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,155,165,175,185,195,205,215,225,235,245,255,265,275,285,295,305 "+diri2+"vo_a.nc "+diri2+"vo_10level_1.nc")
;system("cdo mermean,5,-5 "+diri2+"vo_10level_1.nc "+diri2+"vo_10level_mermean_"+name_model(i)+".nc")
system("cdo sellonlatbox,120,280,5,-5 "+diri2+"vo_10level_1.nc "+diri2+"vo_10level_TP_"+name_model(i)+".nc")
system(" rm "+diri2+"vo_pp.nc" )
system(" rm "+diri2+"vo_a.nc" )
system(" rm "+diri2+"vo_10level_1.nc" )
system(" rm "+diri2+"remapweights.nc" )

system("cdo cat "+diri+"wo/"+name_model(i)+"/r1i1p1/wo_Omon_"+name_model(i)+"_*.nc "+diri2+"wo_pp.nc")
system("cdo genbil,n32 "+diri2+"wo_pp.nc "+diri2+"remapweights.nc")
system("cdo remap,n32,"+diri2+"remapweights.nc "+diri2+"wo_pp.nc "+diri2+"wo_a.nc")
system("cdo intlevel,5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,155,165,175,185,195,205,215,225,235,245,255,265,275,285,295,305 "+diri2+"wo_a.nc "+diri2+"wo_10level_1.nc")
;system("cdo mermean,5,-5 "+diri2+"wo_10level_1.nc "+diri2+"wo_10level_mermean_"+name_model(i)+".nc")
system("cdo sellonlatbox,120,280,5,-5 "+diri2+"wo_10level_1.nc "+diri2+"wo_10level_TP_"+name_model(i)+".nc")
system(" rm "+diri2+"wo_pp.nc" )
system(" rm "+diri2+"wo_a.nc" )
system(" rm "+diri2+"wo_10level_1.nc" )
system(" rm "+diri2+"remapweights.nc" )

end do 

end