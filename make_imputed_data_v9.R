
source("../imputation/make_imputed_data_utils.R")

nimps = 31
make_ADHD_data = function(nimps){
  load("data/new_all_imputed_pADHD_v9s.Rdata")
  analysis_vars = c("cADHD", "cADHDz", "Parity", "cSEX", "cBMI","weeks","preterm",
                    "small_fga","mAge", "mEDU", "mINC", "mSMOK.aftPreg.DailAvg.ord", "mSMOK.EVER", "mALC.EVER",
                    "mDRUGS","mALC.FREQ.durPreg.num","mALC.EFFUNITS", "mLTH", "mSCL", "mBMI", "fAge", "fEDU", "fINC",
                    "fSMOK.EVER", "fALC.EVER", "fSMOK.aftPreg.DailAvg.ord","fALC.FREQ.aftPreg.num","fALC.TYPUNITS",
                    "fDRUGS","HELSEREGION","CivilStatus","imp","uidx", "mADHD", "fADHD")
  
  
  impdata.df = get_impdata(IMP,nimps,analysis_vars)
  
  
  mips = length(unique(impdata.df$imp))
  vars = names(impdata.df)
  
  impdata.df = add_drugsuse(impdata.df)
  # take only those who participated in Q6
  impdata.df = impdata.df[impdata.df$group > 1,] 
  impdata.df$uidx = rep(1:table(impdata.df$imp)[1],nimps)
  
  impdata.df = add_vars(impdata.df,analysis_vars)
  
  imputed_data = data.table(impdata.df[,c(analysis_vars,"group","idx","has_missing")]);
  rm(impdata.df)
  imputed_data[,mALC.EVERx := mALC.EVER]
  imputed_data[,fALC.EVERx := fALC.EVER]
  
  imputed_data[mALC.EFFUNITS == 0, mALC.EVERx := 1]
  imputed_data[fALC.FREQ.aftPreg.num > 0, fALC.EVERx := 2]
  imputed_data[fALC.TYPUNITS > 1, fALC.EVERx := 2]
  imputed_data[fALC.TYPUNITS == 1, fALC.EVERx := 1]
  
  imputed_data[,fALC.EVERx := fALC.EVERx-1]
  imputed_data[,mALC.EVERx := mALC.EVERx-1]
  imputed_data[,fALC.EVER := fALC.EVER-1]
  imputed_data[,mALC.EVER := mALC.EVER-1]
  
  par(mfrow=c(3,3),mar=c(2,2,2,2))
  ridx = imputed_data$imp == 1
  for (v in setdiff(analysis_vars,c("imp","uidx"))){
    if (class(imputed_data[,get(v)]) == "numeric") {
      hist(imputed_data[ridx,get(v)],main = v,ylab = "freq",xlab = "")
    } else {
      plot(table(imputed_data[ridx,get(v)]),main = v,ylab = "freq",xlab = "")
    }
  }
  
  save(imputed_data,file = "data/imputed_preprocessed_ADHDdata_pADHD_v9s_21-50.Rdata")
  
  #########################################
  ###### save info about missingness ######
  #########################################
  d = IMP@data[[1]]
  missings = data.frame(matrix(0,nrow = dim(d)[1], ncol = dim(d)[2]))
  names(missings) = names(d)
  for (v in names(d)) {
    miss = as.numeric(d@variables[[v]]@which_miss)
    if (length(miss) > 0) {
      missings[miss,v] = 1
    }
  }
  missings[,"idx"] = d@variables[["idx"]]@raw_data
  missings[,"group"] = d@variables[["group"]]@raw_data
  
 
  save(missings,file = "data/imputed_preprocessed_ADHDdata_pADHD++missings_v9s_21-50.Rdata")

}

####################################################################################
####################################################################################
################## load and prepare imputed data (CONTINUATION) ####################
####################################################################################
####################################################################################
make_CP_data = function(nimps){
  load("data/new_all_imputed.Rdata")
  impdata.df = get_impdata(IMP,nimps)
  mips = length(unique(impdata.df$imp))
  vars = names(impdata.df)
  
  impdata.df = add_drugsuse(impdata.df)
  impdata.df$uidx = rep(1:table(impdata.df$imp)[1],nimps)
  
  
  ########################## complete and save data #########################
  
  analysis_vars = c("Parity", "cSEX", "cBMI","weeks","preterm",
                    "small_fga","mAge", "mEDU", "mINC", "mSMOK.aftPreg.DailAvg.ord",
                    "mDRUGS","mALC.FREQ.durPreg.num","mALC.EFFUNITS", "mLTH", "mSCL", "mBMI", "fAge", "fEDU", "fINC",
                    "fSMOK.aftPreg.DailAvg.ord","fALC.FREQ.aftPreg.num","fALC.TYPUNITS","fDRUGS","HELSEREGION","CivilStatus","imp","uidx","has_missing")
  
  impdata.df = add_vars(impdata.df,analysis_vars)
  
  imputed_data = impdata.df[,c(analysis_vars,"group","idx")];
  imputed_data$continuation = factor(impdata.df$group>1,labels = c("No","Yes"))
  
  save(imputed_data,file = "imputed_preprocessed_CPdata++.Rdata")
  save(imputed_data,file = "../stan_hlms/data/imputed_preprocessed_CPdata++.Rdata")
}
make_ADHD_data(nimps)