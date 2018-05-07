args=(commandArgs(TRUE))


for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

library(mi)
source("021_mi_utils.R")


if (task == 1) {
  tmp = load_RD("data/RDallnew_v9_2001_2009.Rdata")
  RD = tmp[[1]]
  vars = tmp[[2]]
  
  RD$mFeel_Useless_Q1_refl = abs(RD$mFeel_Useless_Q1-5)
  RD$mFeel_NotProud_Q1_refl = abs(RD$mFeel_NotProud_Q1-5)

  for (i in 1:3) {
    v_name = paste("mDEP",i,sep = "")
  }
  
  ipdata = make_ipdata_scales(RD,vars)
  ipdata = ipdata[,-which(names(ipdata)=="valid_cADHD")]
  ipdata = ipdata[,-which(names(ipdata)=="fADHDdiag")]
  rm(tmp,RD)
  mdf = make_mdf(ipdata)
  save(mdf,file = "../IPW/data/new_all_mdf_pADHD_v9.Rdata")
} else if (task == 2){
  load("../IPW/data/new_all_mdf_pADHD_v9.Rdata")

  seeds = c(1017, 8664, 7458, 9720)
  fn = paste("../IPW/data/new_all_imputed_pADHD_",iter,"_v9s.Rdata",sep = "")
  mdf@workpath = paste(paste(getwd(),"/mi_tmp",iter,"_",sep = ""),gsub(":","",gsub(" ","_",Sys.time())),sep = "")
  dir.create(mdf@workpath)
  options(error = quote(dump.frames(paste("imp_error_dump",iter,sep = ""), TRUE)))
  IMP <- mi(mdf, n.chains = 1, max.minutes = 5500, n.iter = 50,parallel = F,verbose = T, seed = seeds[iter])
  save(IMP,file = fn)
  warnings()

}  else if (task == 3){
  fn = paste("data/new_all_imputed_",iter,".Rdata",sep = "")
  load(fn)
  IMP@total_iters = 27L
  IMP <- mi(IMP, max.minutes = 2100, n.iter = 13,parallel = F,verbose = T)
  fn = paste("data/new_all_imputed_",iter,"b.Rdata",sep = "")
  save(IMP,file = fn)
}
