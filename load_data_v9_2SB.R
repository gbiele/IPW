library(rstan)
library(plyr)
library(data.table)

load_data = function() {
  predictors = c("Parity.cat", "cSEX", "small_fga", 
                 "mSMOKE", "mSMOKE.5Cs", 
                 "mDRUGS","mDRUGS.SCORE",
                 "mALC","mALC.FREQ", "mALC.EFFUNITS" ,  
                 "preterm",
                 "mLTHD", "mSCL5", "mBMIdev",
                 "fAge", "fEdu.cat", 
                 "fSMOKE", "fSMOKE.5Cs", 
                 "fALC","fALC.FREQ","fALC.TYPUNITS",
                 "fDRUGS", "fDRUGS.SCORE",
                 "mEdu.cat","mEdu.cat2","mAge.cat","mAge.cat2",
                 "mADHD","fADHD") # "mBMI" "mINC", "mSCL5", "fINC", "fAge_sq",
  
  load("data/imputed_preprocessed_cleaned_ADHDdata_pADHD_v9s.Rdata")
  
  criterion = "cADHD"
  # rename a few variables
  imputed_data = data.table(imputed_data)
  setnames(imputed_data,"mSMOK.aftPreg.DailAvg.ord","mSMOKE")
  setnames(imputed_data,"mALC.FREQ.durPreg.num","mALC.FREQ")
  setnames(imputed_data,"fSMOK.aftPreg.DailAvg.ord","fSMOKE")
  setnames(imputed_data,"fALC.FREQ.aftPreg.num","fALC.FREQ")
  setnames(imputed_data,"mLTH","mLTHD")
  setnames(imputed_data,"mSCL","mSCL5")
  
  imputed_data[,small_fga := 1-small_fga]

  imputed_data[preterm > 1, preterm := 1]
  
  imputed_data[, cSEX := as.numeric(cSEX)-1]
  
  imputed_data[,CivilStatus := as.numeric(CivilStatus)-1]
  
  
  
  if("mALC" %in% predictors) {
    imputed_data[, fALC := as.numeric(factor(fALC.EVER))-1]
    imputed_data[, mALC := as.numeric(factor(mALC.EVER))-1]
  }
  
  if("mDRUGS.SCORE" %in% predictors) {
    imputed_data[, mDRUGS.SCORE := mDRUGS]
    imputed_data[, mDRUGS := (mDRUGS.SCORE > min(imputed_data[, mDRUGS]))*1]
    imputed_data[mDRUGS > 0, mDRUGS.SCORE := scale(mDRUGS.SCORE), by = "imp"]
    imputed_data[, fDRUGS.SCORE := fDRUGS]
    imputed_data[, fDRUGS := (fDRUGS.SCORE > min(imputed_data[, fDRUGS]))*1]
    imputed_data[fDRUGS > 0, fDRUGS.SCORE := scale(fDRUGS.SCORE), by = "imp"]
  }
  
  if("mSMOKE.5Cs" %in% predictors) {
    imputed_data[, mSMOKE.5Cs := mapvalues(mSMOKE,1:6,c(0,.5,1.5,3,7,15))]
    imputed_data[, mSMOKE := (mSMOKE.5Cs > 0)*1]
    imputed_data[mSMOKE > 0, mSMOKE.5Cs := scale(mSMOKE.5Cs), by = "imp"]
    imputed_data[, fSMOKE.5Cs := mapvalues(fSMOKE,1:6,c(0,.5,1.5,3,7,15))]
    imputed_data[, fSMOKE := (fSMOKE.5Cs > 0)*1]
    imputed_data[fSMOKE > 0, fSMOKE.5Cs := scale(fSMOKE.5Cs), by = "imp"]
  }
  
  imputed_data[,fALC.TYPUNITS := mapvalues(fALC.TYPUNITS,1:7,c(0,.5,1.75,3.5,5.625,8.375,12.5))]
  imputed_data[fALC.TYPUNITS > 0, fALC := 1]
  imputed_data[, fALC.TYPUNITS := scale(fALC.TYPUNITS), by = "imp"]
  imputed_data[fALC.FREQ > 0, fALC := 1]
  imputed_data[, fALC.FREQ := scale(fALC.FREQ), by = "imp"]
  
  
  imputed_data[mALC.EFFUNITS > 8, mALC.EFFUNITS := 9]
  imputed_data[mALC == 1 & mALC.EFFUNITS == 0, mALC.EFFUNITS := 1]
  imputed_data[, mALC.EFFUNITS := scale(mALC.EFFUNITS),by = "imp"]

  imputed_data[, Age := cut(mAge,breaks = c(13,19,24,29,34,39,49),labels = c("<20","20-24","25-29","30-34","35-39","40-49"),ordered = T)]
  imputed_data[, mAge.cat := mean(mAge), by = Age]
  imputed_data[, mAge.cat := scale(mAge.cat), by = imp]
  imputed_data[, mAge.cat2 := scale(poly(mAge.cat,2)[,2],center = F), by = "imp"]
  imputed_data[, fAge := scale(fAge), by = "imp"]
  
  
  imputed_data[,Parity.cat := parity]
  
  
  imputed_data[, medu := cut(mEDU,breaks = c(.5,1.5,4.5,5.5,6.5),ordered = T)]
  imputed_data[,mEdu.cat := as.numeric(medu)-2.5]
  imputed_data[, mEdu.cat2 := scale(poly(mEdu.cat,2)[,2],center = F), by = "imp"]
  imputed_data[, fedu := cut(fEDU,breaks = c(.5,1.5,4.5,5.5,6.5),ordered = T)]
  imputed_data[,fEdu.cat := as.numeric(fedu)-2.5]
  
  imputed_data[mSCL5 > 15, mSCL5 := 15]
  imputed_data[, mSCL5 := scale(mSCL5), by = "imp"]
  imputed_data[, mLTHD := scale(mLTHD), by = "imp"]

  
  imputed_data[, mADHD := scale(mADHD), by = "imp"]
  imputed_data[, fADHD := scale(fADHD), by = "imp"]
  
  
  mBMImode = as.numeric(names(which.max(table(round(imputed_data$mBMI)))))
  imputed_data[,mBMIdev := scale((mBMI/mBMImode)), by = "imp"]

  imputed_data = imputed_data[,c(criterion,predictors,"Age","edu","parity","has_missing","group","imp","uidx","idx"),with = F]


  return(list("imputed_data" = imputed_data,
              "criterion" = criterion,
              "predictors" = predictors))
  print(table(imputed_data$preterm))
}


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}