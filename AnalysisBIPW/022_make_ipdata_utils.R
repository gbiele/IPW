#########################################################################
#########################################################################
############### PREPARE DATA FOR IMPUTATION PROCESS #####################
#########################################################################
#########################################################################

clean_ipdata = function(ipdata){
  
  vars = names(ipdata)
  drops =  vars[c(grep("_w",vars),grep("_we",vars),grep("_bef",vars),#which(vars == "mLTHD3plus"),
              grep("Heroine",vars),grep("fALC.BINGE",vars),grep("fINC.",vars))]
  
  drop_from_drops = c("mSMOK_aftPreg_DailAvg","fSMOK_aftPreg_DailAvg","fALC_FREQ_aftPreg_num","fALC_FREQ_befPreg_num")
  drops = drops[!(drops %in% drop_from_drops)]
  ipdata = ipdata[,!(vars %in% drops)]
  vars = names(ipdata)
  
  
  # drop non-numerical variables that are count variables
  drops = c()
  for (v in vars[grep("_num",vars)]) {
    if (length(c(grep("TYPUNITS",v),grep("INC",v))) == 0) {
      drops = c(drops,-which(vars == substr(v,1,nchar(v)-4)))
    } else {
      if (v != "fALC_FREQ_aftPreg_num") drops = c(drops,-which(vars == v))
    }
  }
  ipdata = ipdata[,drops]
  
  ######## remove cases with more than 50% of missing variables ########
  x = rowSums(is.na(ipdata))
  drops = unique(which(x>(ncol(ipdata)/2)))
  ipdata = ipdata[-drops,]
  rownames(ipdata) = 1:nrow(ipdata)
  
  ######## recode smoking vars ########
  smokevars = vars[grep("DailAvg",vars)]
  for (sv in smokevars) {
    ipdata[,paste(sv,'_ord',sep = "")] = cut(ipdata[,sv],breaks = c(-1,0,1,2,4,10,20,100),labels = c("0","0.1-1","1.1-2","2.1-4","4.1-10","10.1-20",">20"))
    ipdata = ipdata[,-which(names(ipdata)==sv)]
  }
  vars = names(ipdata)
  
  
  # set correct class
  cvars = names(ipdata)[-grep("DailAvg_ord|byear|Parity|fEDU|mEDU|HELSEREGION",names(ipdata))]
  for (v in cvars){
    us = length(table(ipdata[,v]))
    if (us == 2) {
      ipdata[,v] = factor(ipdata[,v],ordered = F)
    } else if (us <= 7) {
      ipdata[,v] = factor(ipdata[,v],ordered = T)
    } else if (us > 7) {
      ipdata[,v] = as.numeric(ipdata[,v])
    }
  }

  vars = names(ipdata)
  ############### join small classes #################
  
  for (v in names(ipdata)[grep("fDU",names(ipdata))]){
    ipdata[,v] = cut(as.numeric(ipdata[,v]),breaks = c(.5:2.5,4.5),labels = c("Never","Prev", "Recent"))
  }
  for (v in names(ipdata)[grep("mDU",names(ipdata))]){
    ipdata[,v] = cut(as.numeric(ipdata[,v]),breaks = c(.5:1.5,4.5),labels = c("Never","PrevRec"))
  }
  
  non_numeric = names(ipdata)[grep("factor",lapply(ipdata,class))]
  for (iv in non_numeric) {
    y = ipdata[,iv]
    ylevels = rev(levels(y))
    for (L in ylevels) {
     #Lprop = 100*(table(y)/sum(table(y)))[names(table(y)) == L]
      #if (Lprop < .5){ #(sum(y==L,na.rm = T) < 25){
      if (sum(y == L,na.rm = t) < 500){
        #print(paste("LEVEL REMOVED FROM VARIABLE ",iv,": ",L," (",sum(y==L,na.rm = T),")",sep = ""))
        print(paste("LEVEL REMOVED FROM VARIABLE ",iv,": ",L," (only ",sum(y == L,na.rm = t)," occurences)",sep = ""))
        y[which(y==L)] = ylevels[which(L == ylevels)+1]
      }
    }
    newlabels = levels(y)[table(y) > 0]
    if (length(newlabels) >1) {
      next_to_max = newlabels[length(newlabels)-1]
      if (is.numeric(next_to_max)) {
        next_to_max = max(as.numeric(strsplit(gsub("[^[:alnum:] ]", " : ", next_to_max)," : ")[[1]]),na.rm = T)
        newlabels[length(newlabels)] = paste0(">=",next_to_max)
      }
    }
    ipdata[,iv] = factor(as.numeric(y),labels = newlabels, ordered = length(newlabels) > 2)
  }
  

  ########## adjust name of variables
  names(ipdata)[names(ipdata) == "MOBA.3"] = "MOBA3"
  names(ipdata) = gsub("_Q1","",names(ipdata))
  names(ipdata) = gsub("_",".",names(ipdata))
  for (p in seq(15,1,-1)) {
    rs = paste(".",p,sep = "")
    names(ipdata) = gsub(as.character(p),rs,names(ipdata))
  }
  
 
  ipdata$mALC.EFFUNITS[ipdata$mALC.EFFUNITS > 12] = 13
  
  ipdata$group = factor(ipdata$group,ordered = F)
  
  ipdata$mFeel = factor(ipdata$mFeel,ordered = T)
  #ipdata$fDU.sc.Hash[ipdata$fDU.sc.Hash == "LstMbfPreg"] = "Prev"
  #ipdata$fDU.sc.Hash = ordered(as.numeric(ipdata$fDU.sc.Hash),labels = c("Never","Prev"))
  #ipdata$fDU.sc.Amph[ipdata$fDU.sc.Amph == "LstMbfPreg" | ipdata$fDU.sc.Amph == "Recently"] = "Prev"
  #ipdata$fDU.sc.Amph = ordered(as.numeric(ipdata$fDU.sc.Amph),labels = c("Never","Prev"))
  
  ipdata$weeks.transformed[is.na(ipdata$weeks.transformed)] = NA
  
  ipdata$CivilStatus[ipdata$CivilStatus == "Single"] = "Other"
  ipdata$CivilStatus = factor(ipdata$CivilStatus,ordered = F)
  
  ipdata$HELSEREGION = factor(ipdata$HELSEREGION,ordered = F)
  
  return(ipdata)
  par(ask = F)
}

############################################################################### 
#################### function to make sumscores for scales #################### 
############################################################################### 

make_sum_scores = function(ipdata,scales){ 
  fADHDdiag = ipdata$fADHDdiag
  ipdata = ipdata[,names(ipdata) != "fADHDdiag"]
  for (sc in scales){
    S = ipdata[,grep(sc,names(ipdata))]
    many_missings =  which(rowMeans(is.na(S)) > .5)
    ipdata[,sc] = round(rowMeans(S,na.rm = T)*ncol(S))
    ipdata[,sc] = ipdata[,sc]-min(ipdata[,sc],na.rm = T)
    ipdata[many_missings,sc] = NA
    ipdata = ipdata[,-which(names(ipdata) %in% names(S))]
  }
  ipdata$fADHDdiag = fADHDdiag
  ipdata$mALC.RAPI = cut(ipdata$mALC.RAPI,breaks = c(-1:4,10),labels = c(paste(0:4),"5-10"))
  ipdata$mALC.TACE = cut(ipdata$mALC.TACE,breaks = c(-1,0,1,3),labels = c("0","1","2-3"))

  return(ipdata)
}