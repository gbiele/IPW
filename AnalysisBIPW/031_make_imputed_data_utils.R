########### service functions for make_imputetd_data ##################
#######################################################################
library(data.table)
library(mi)
library(ggplot2)
library(psych)
library(ltm)
old.par = par()



############## check if imputation results are good ###################
check_mi = function(IMP){
  df = round(mipply(IMP, mean, to.matrix = TRUE), 3)
  cv.df = data.frame(apply(df,1,sd)/apply(df,1,mean))
  cv.df$var = rownames(cv.df)
  names(cv.df) = c("cv","var")
  ggplot(cv.df[cv.df$cv > .001,],aes(x = cv, y = var)) + geom_point(size = 3)
  Rhat = Rhats(IMP)
  plot(IMP)
  hist(IMP)
}


####################################################################################
#################### load and prepare imputed data (ADHD) ##########################
####################################################################################


get_impdata = function(IMP,nimps = 10, analysis_vars){
  check_missing_vars = setdiff(analysis_vars,c("mINC","fINC","HELSEREGION","imp", "uidx" ))
  impdata = complete(IMP,m = nimps)
  
  impdata.df = impdata[[1]]

  impdata.df$imp = 1
  for (i in 2:length(impdata)){
    impdata[[i]]$imp = i
    impdata.df = rbind(impdata.df,impdata[[i]])
  }
  impdata.df$imp = factor(impdata.df$imp,ordered = F)

  vars =  names(impdata.df)[-grep("missing",names(impdata.df))]
  impdata.df = impdata.df[,vars]
  
  for (v in vars[grep("FREQ",vars)]){
    tmp = impdata.df[,v]
    impdata.df[,v] = as.numeric(impdata.df[,v])
    for (L in levels(tmp)) impdata.df[tmp == L,v] = as.numeric(L)
  }
  for (v in c(vars[grep("DailAvg",vars)],"fALC.TYPUNITS")){
    impdata.df[,v] = as.numeric(impdata.df[,v])
  }
  
  # find cases with imputed values
  all_vars = names(IMP@data[[1]])
  missings = data.table(data.frame(matrix(0,ncol = length(all_vars),nrow = nrow(IMP@data[[1]]))))
  setnames(missings,names(missings),all_vars)
  missings$pattern = IMP@data[[1]]@patterns
  
  for (L in levels(IMP@data[[1]]@patterns)) {
    if (length(grep(",",L)) > 0) {
      missing_vars = strsplit(L,", ")[[1]]
    } else {
      missing_vars = L
    }
    missings[pattern == L,(missing_vars) := 1]
  }
  rm(missing_vars,L)
  
  missings = missings[,pattern := NULL]
  
  missings[,cADHDz := cADHD]
  missings[,preterm := weeks.transformed]
  missings[,weeks := weeks.transformed]
  missings[,cBMI := weight + length]
  missings[,mBMI := mWght + mTall]
  missings[,small_fga := weight + weeks]
  missings[,mDRUGS := mDU.sc.Hash + mDU.sc.Amph + mDU.sc.Ecst + mDU.sc.Coc]
  missings[,fDRUGS := fDU.sc.Hash + fDU.sc.Amph + fDU.sc.Ecst + fDU.sc.Coc]
  missings[,(names(missings)[grep("fDU.",names(missings))]) := NULL]
  missings[,(names(missings)[grep("mDU.",names(missings))]) := NULL]
  
  missings[,nothing := NULL]
  
  impdata.df$has_missing = rep(apply(missings[,check_missing_vars,with = F] > 0,1,max),nimps)
  
  return(impdata.df)  
}

########## make histogram of a number of numeric and categorical variables ##########
hists = function(df){
  par(mar=c(2,2,2,2))
  nms = names(df)
  for (b in ((1:ceiling(length(df)/25))-1)*25) {
    nr = floor(sqrt(dim(df)[2]))
    nc = ceiling(dim(df)[2]/nr)
    par(mfrow=c(nr,nc))
    for (v in nms[(1:min(c(25,length(nms)-b)))+b]) {
      y = df[,v]
      if (length(table(y))>2 & class(y) == "numeric") {
        hist(as.numeric(y),main = v,xlab = NULL,ylab = NULL) 
      } else {
        barplot(table(y),main = v,xlab = NULL,ylab = NULL)
      }
    }
  }
}



####################################################################################
############### IRT models for mothers and fathers drug use (ADHD) #################
####################################################################################
add_drugsuse = function(impdata.df){
  vars = names(impdata.df)
  dat = impdata.df[,grep("fDU",vars)]
  for (i in 1:ncol(dat)) dat[,i] = as.numeric(dat[,i])
  gpcm_model = gpcm(dat)
  fs = factor.scores.gpcm(gpcm_model)
  impdata.df$fDRUGS = NA
  for (i in 1:nrow(fs$score.dat)) {
    idx = which(rowSums(sapply(1:4,function(x) dat[,x] == fs$score.dat[i,x])) == 4)
    impdata.df$fDRUGS[idx] = fs$score.dat$z1[i]-min(fs$score.dat$z1)
  }
  
  
  dat = impdata.df[,grep("mDU",vars)]
  dat = dat[,-grep("mDU.sc.Ecst",names(dat))]
  for (i in 1:ncol(dat)) dat[,i] = as.numeric(dat[,i])
  gpcm_model = gpcm(dat)
  fs = factor.scores.gpcm(gpcm_model)
  impdata.df$mDRUGS = NA
  for (i in 1:nrow(fs$score.dat)) {
    tm = data.frame(t(matrix(rep(fs$score.dat[i,1:ncol(dat)],nrow(dat)),ncol = nrow(dat))))
    idx = which(rowSums(dat == tm) == ncol(dat))
    impdata.df$mDRUGS[idx] = fs$score.dat$z1[i]
  }
  
  rm(dat,tm)
  return(impdata.df)
}

####################################################################################
########################### recalculate some variables #############################
####################################################################################

add_vars = function(impdata.df,analysis_vars){
  ########################## small for gestational week ##############################
  
  impdata.df$small_fga = NA
  impdata.df$short_fga = NA
  impdata.df$weeks = round(48-exp(impdata.df$weeks.transformed))
  for (w in unique(impdata.df$weeks)){
    for (cs in 1:2) {
      idx = impdata.df$weeks == w & impdata.df$cSEX == cs
      if (sum(idx) > 10) {
        impdata.df$small_fga[idx] = 1*(quantile(impdata.df[idx,"weight"],.1) < impdata.df[idx,"weight"])
        impdata.df$short_fga[idx] = scale(impdata.df[idx,"length"])
      } else {
        impdata.df$small_fga[idx] = 0
        impdata.df$short_fga[idx] = 0
      }
    }
  }
  # weeks deviation from the 40th pregnancy week
  impdata.df$preterm = abs(as.numeric(cut(impdata.df$weeks,breaks = c(0,27,31,36,45)))-4)
  
  impdata.df$mBMI = impdata.df$mWght/(impdata.df$mTall/100)^2
  impdata.df$cBMI = (impdata.df$weight/1000)/(impdata.df$length/100)^2
  
  if ( "cADHDz" %in% analysis_vars) impdata.df$cADHDz = as.matrix(scale(impdata.df$cADHD))
  
  factor_vars = c("cSEX","HELSEREGION","CivilStatus","imp")
  numeric_vars = setdiff(analysis_vars,factor_vars)
  for (v in analysis_vars){
    if (v %in% factor_vars) {
      impdata.df[,v] = as.factor(impdata.df[,v])
    } else {
      impdata.df[,v] = as.numeric(impdata.df[,v])
    }
  }
  
  #if (length(grep("PP.",analysis_vars)) < 1) impdata.df = impdata.df[complete.cases(impdata.df[,analysis_vars]),]
  return(impdata.df)
}
