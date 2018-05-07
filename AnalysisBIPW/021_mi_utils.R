
source("022_make_ipdata_utils.R")

############################################################################### 
################################# load data ################################### 
############################################################################### 

load_RD = function(filename){

  load(filename)
  
  RD$mFeel_Useless_Q1_refl = abs(RD$mFeel_Useless_Q1-5)
  RD$mFeel_NotProud_Q1_refl = abs(RD$mFeel_NotProud_Q1-5)
  RD = RD[,-which(names(RD) == "mFeel_Useless_Q1")]
  RD = RD[,-which(names(RD) == "mFeel_NotProud_Q1")]
  
  for (i in 1:3) {
    v_name = paste("mDEP",i,sep = "")
    RD[,paste(v_name,"refl",sep = "_")] = abs(RD[,v_name]-6)
    RD = RD[,-which(names(RD) == v_name)]
  }
  
  #code to verify number of missing cADHD items: 
  nnacADHD = rowSums(is.na(RD[,grep("cADHD",names(RD))]))
  # xtabs(~byear + nnacADHD,RD)
 
  RD$length[RD$weight > 2000 & RD$length < 37.5] = NA
  RD$weight[which(RD$length > 45 & RD$weight < 1500)] = NA

  RD = RD[-which(RD$weeks<32 |
                 RD$weeks>44 |
                 RD$weight<1000 |
                 RD$length<35) ,]

    
  vnms = names(RD)
  vars = c("HELSEREGION","byear","Parity","CivilStatus","mAge","mEDU","mWght","mTall",
           "fAge","fEDU","fADHDdiag",
           "weeks_transformed","weight","length","cSEX",
                 vnms[grep("INC",vnms)],
                 vnms[grep("SMOK_EVER",vnms)],
                 vnms[grep("Preg_DailAvg",vnms)],
                 vnms[grep("ALC",vnms)],
                 vnms[grep("DU_sc_",vnms)],
                 vnms[grep("SCL",vnms)],
                 vnms[grep("Feel",vnms)],
                 vnms[grep("mDEP",vnms)],
                 vnms[grep("^mLTHD",vnms)],
                "group")
  vars = vars[-grep("fALC_TYPUNITS_",vars)]
  rm(vnms)
  RD$idx = 1:dim(RD)[1]
  vars = unique(c("idx",vars))
  return(list(RD,vars))
}


##################################################################################
############## function for data to predict ADHD scores & particip. ##############
######################### SUMSCORES BASED ON RAW DATA!!  #########################
##################################################################################
make_ipdata_scales = function(RD,vars){
  RDvars = names(RD)
  vars = c(vars,
           RDvars[grep("^mANX",RDvars)],
           RDvars[grep("^mDEP",RDvars)],
           RDvars[grep("^mADHD[0-9]",RDvars)],
           RDvars[grep("^fADHD[0-9]",RDvars)],
           RDvars[grep("^cADHD[0-9]",RDvars)])
  
  #ipdata = RD[RD$MOBA0 & RD$MOBA3,intersect(names(RD),vars)]
  ipdata = RD[RD$MOBA0,intersect(names(RD),vars)]
  
  
  names(ipdata)[names(ipdata) == "mALC_TACE_EFFUNITS"] = "mALC_EFFUNITS"
  scales = c("mALC.RAPI","mALC.TACE","mSCL","mFeel","mLTH","mADHD","fADHD","mDEP","mANX","cADHD")
  #ipdata$mLTHD3plus = (as.numeric(ipdata$mLTHD3plus)-1)*2 # added Mai 2015
  names(ipdata)[names(ipdata) == "mLTHD3plus"] = "mLTxHD3plus"
  ipdata = make_sum_scores(ipdata,scales)
  names(ipdata)[names(ipdata) == "mLTxHD3plus"] = "mLTHD3plus"
  
  ipdata = clean_ipdata(ipdata)
  
  ipdata$valid_cADHD = !is.na(ipdata$cADHD)
  
  return(ipdata)
}


##################################################################################
############################## plot included vars. ###############################
##################################################################################

hists = function(df){
  par(mar=c(2,2,2,2))
  nms = names(df)
  for (b in ((1:ceiling(length(df)/25))-1)*25) {
    par(mfrow=c(5,5))
    for (v in nms[(1:min(c(25,length(nms)-b)))+b]) {
      y = df[,v]
      if (length(table(y))>5 & class(y) == "numeric") {
        hist(as.numeric(y),main = v,xlab = NULL,ylab = NULL) 
      } else {
        barplot(table(addNA(y)),main = v,xlab = NULL,ylab = NULL)
      }
    }
  }
}

hists_mdf = function(mdf){
  library(fitdistrplus)
  par(ask=T) 
  for (k in names(mdf@index)){
    mdfclass = attr(mdf[,k],"class")[1]
    y = attr(mdf[,k],"raw_data")
    nNA = sum(is.na(y))
    varname = paste(mdfclass," ",k,", ",nNA," missing",sep = "")
    if (mdfclass == "unordered-categorical" | mdfclass == "ordered-categorical" | mdfclass == "binary") {
      t = table(addNA(y))
      names(t)[length(t)] = "NA"
      barplot(t,main = varname)
    } else if (mdfclass == "positive-continuous" | mdfclass == "continuous" | mdfclass == "nonnegative-continuous") {
      y = y[!is.na(y)]
      if (mdfclass == "nonnegative-continuous") y = y[y>0]
      fd = fitdist(y,"norm")
      h = hist(y,breaks = 25,plot = F)
      d = dnorm(h$mids,mean = fd$estimate["mean"],sd = fd$estimate["sd"])
      d = d/sum(d)*sum(h$count)
      plot(h$mids,h$counts,"s",ylim = range(c(nNA,h$count,d)),xlab = "value",main = varname,ylab = "frequency",lwd = 2)
      lines(h$mids,d,'l',col = "red",lwd = 2)
      lines(max(h$mids)-.2,nNA,'h',col = "green",lwd = 2)
      legend("topright",legend = c("observed","model"),col = c("black","red"),cex = .75,lty = 1,bty = "n")
    } else if (mdfclass == "count") {
      y = y[!is.na(y)]
      fd = fitdist(y,"nbinom")
      t = table(y)
      x =  sort(unique(y))
      d = dnbinom(x,mu = fd$estimate["mu"],size = fd$estimate["size"])
      d = d/sum(d)*sum(t)
      plot(t,ylim = range(c(t,d)),xlab = "value",main = varname,ylab = "frequency")
      segments(x,rep(0,length(x)),x,t)
      segments(x+.2,rep(0,length(x)),x+.2,d,col = "red",lwd = 2)
      lines(max(x)-.2,nNA,'h',col = "green",lwd = 2)
      legend("topright",legend = c("observed","model"),col = c("black","red"),cex = .75,lty = 1,bty = "n")
    } else {print(paste("no plot for",varname))}
    
  }
  par(ask=F) 
}


##################################################################################
############################## make mdf from raw data ############################
##################################################################################

make_mdf = function(ipdata){
  
  #mdf = mdf = missing_data.frame(ipdata[sample(nrow(ipdata),2500),],favor_positive = TRUE)
  mdf = missing_data.frame(ipdata,favor_positive = TRUE)
  
  
  count_vars = c("mADHD","cADHD","fADHD","mANX","mALC.EFFUNITS","mSCL","mDEP")
  mdf <- change(mdf, y = count_vars, what = "model", to = "zeroinfl") 
  
  continuous_vars = c("byear","mAge","fAge","mTall",
                      "weight","length","weeks.transformed")
  mdf <- change(mdf, y = continuous_vars, what = "type", to = "continuous") 
  
  pos_continuous_vars = c("mWght")
  mdf <- change(mdf, y = pos_continuous_vars, what = "type", to = "positive-continuous")
  
  # makes problems with zeroinfl model
  mdf <- change(mdf, y = c("mALC.EFFUNITS"), what = "type", to = "count")
  
  
  #irrelvant_vars = c("group","idx","valid_cADHD","mALC.EVER","fALC.EVER")
  irrelvant_vars = c("group","idx","valid_cADHD","HELSEREGION","mDEP","mANX")
  mdf <- change(mdf, y = irrelvant_vars, what = "type", to = "irrelevant") 
  
  missingness_columns =  grep("missing",colnames(mdf@X))
  for (iv in names(mdf@index)){
    dropvars = c(grep("^idx|^group",colnames(mdf@X)),grep(paste("^",iv,sep = ""),colnames(mdf@X)))
    mdf@index[[iv]] = sort(c(mdf@index[[iv]],missingness_columns))
  }
  
  return(mdf)
}


##################################################################################
################## merge different iterations of an imputation ###################
##################################################################################


merge_IMPs = function(fns){
    n.chains = length(fns)
    mdfs = vector(mode = "list",length = n.chains)
    for (k in 1:n.chains) {
      load(fns[k])
      print(paste(IMP@total_iters,"in",fns[k]))
      mdfs[[k]] = IMP@data[[1]]
      }
    IMP <- new("mi", 
               call        = IMP@call,
               data        = mdfs,
               total_iters = IMP@total_iters)
  return(IMP)
}


