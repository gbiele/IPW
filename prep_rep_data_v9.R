######################################################################################
######################################################################################
###########  COLLECT DATA FROM DIFFERENT SOURCES INTO ONE data.frame #################
######################################################################################
######################################################################################

data_root = "F:/Forskningsprosjekter/PDB 2100 - The importance of pr_/Forskningsfiler/GUBI/"

setwd("F:/Forskningsprosjekter/PDB 299 - ADHD-studien Prescho_/Forskningsfiler/GUBI/GuidoData/representativeness/IPW")
library(psych)
library(ggplot2)
library(matrixStats)
library(reshape) 
library(plyr)
library(foreign)
library(lme4)
library(MASS)
library(haven)
#ggplot(RD,aes(x = mINC, y = mEDU, group = 1)) + stat_sum(aes(size = ..n..), geom = "point") + geom_smooth(method = lm)

#####################################################
################ load MFR data ######################
#####################################################

# MFR, to get Year of Birth
MFR =read_sav("F:/Forskningsprosjekter/PDB 2100 - The importance of pr_/Forskningsfiler/GUBI/PDB2014_MFR_410_v9.sav");
# get participants in the sampling period
#MFR = MFR[MFR$FAAR > 2003 & MFR$FAAR < 2008,]
Cidx = c("PREG_ID_2014","BARN_NR","HELSEREGION","FAAR","MORS_ALDER",
         "FARS_ALDER_KAT_K8","SIVST_2","PARITET_5","EPILEPSI","MULTIVITU",
         "FOLATU","SVLEN_UL","SVLEN_SM","VEKT","LENGDE","KJONN","APGAR5",
         "C00_MALF_ALL","DODKAT", "PERINAT_DODFODT_MFR")
MFR = MFR[,Cidx]
names(MFR) = c("PID","CN","HELSEREGION","byear","mAge","fAge","MarSt","Parity","Epil","MULTIVU","FOLATU","weeks_UL","weeks_MS","weight","length","cSEX","APGAR5",
               "malformation","living","db_MFR")
MFR$cSEX = as.factor(MFR$cSEX)

MFR$living = factor(MFR$living == 0 &  MFR$db_MFR == 0, labels =c("No","Yes"))
MFR$HELSEREGION = as_factor(MFR$HELSEREGION)
contrasts(MFR$HELSEREGION) = contr.treatment(levels(MFR$HELSEREGION),base = 1)
MFR$malformation = factor(MFR$malformation,labels =c("No","Yes"))

fAge = c(17,22,27,32,37,42,47,55)
for (c in seq(8)) {MFR$fAge[MFR$fAge == c] = fAge[c]}
MFR$fAge = as.numeric(MFR$fAge)

RD = MFR
remove(MFR,Cidx)
RD$pMFR = 1

#####################################################
######## add data from 1st Questionnaire ############
#####################################################
# for eduction, income:
load("F:/Forskningsprosjekter/PDB 2100 - The importance of pr_/Forskningsfiler/GUBI/PDB2014_Q1_v9_red.Rdata");


nms = c("PID","mpreWght1", "mpreWght2","mTall", "fTall","fWght", "mEDUcomp","mEDUcur","fEDUcomp","fEDUcur","mINC","fINC",
        "mSCL_Fear_Q1", "mSCL_Nervous_Q1", "mSCL_Hopeless_Q1", "mSCL_Blue_Q1", "mSCL_Worry_Q1",
        "mFeel_Pos_Q1", "mFeel_Useless_Q1", "mFeel_NotProud_Q1", "mFeel_Valuabel_Q1",
        "mLTHD1", "mLTHD2","mLTHD3","mLTHD4","mLTHD5","mLTHD3plus",
        paste("mDU_Hash",c("Never","Prev","LstMbfPreg","DurPreg"),sep = "_"),
        paste("mDU_Amph",c("Never","Prev","LstMbfPreg","DurPreg"),sep = "_"),
        paste("mDU_Ecst",c("Never","Prev","LstMbfPreg","DurPreg"),sep = "_"),
        paste("mDU_Coc",c("Never","Prev","LstMbfPreg","DurPreg"),sep = "_"),
        paste("mDU_Heroine",c("Never","Prev","LstMbfPreg","DurPreg"),sep = "_"),
        "mALC_EVER", "mALC_FREQ_befPreg", "mALC_FREQ_durPreg", "mALC_BINGE_befPreg", "mALC_BINGE_durPreg", "mALC_TYPUNITS_befPreg", "mALC_TYPUNITS_durPreg",
        "mALC_TACE_EFFUNITS", "mALC_TACE_HurtCrit", "mALC_TACE_OughtLESS", "mALC_TACE_DrinkMorn",
        "mALC_RAPI_argued", "mALC_RAPI_SuddFound", "mALC_RAPI_Absence", "mALC_RAPI_Fainted", "mALC_RAPI_SadPeriod",
        "CivilStatus")
names(MoBa1) = nms


MoBa1$CivilStatus[which(MoBa1$CivilStatus == 0)] = NA
MoBa1$CivilStatus[which(MoBa1$CivilStatus %in% c(2,4,6))] = 4
MoBa1$CivilStatus[which(MoBa1$CivilStatus == 5)] = 3.5
MoBa1$CivilStatus[is.na(MoBa1$CivilStatus)] = NA
MoBa1$CivilStatus = factor(as.numeric(MoBa1$CivilStatus),labels = c("married","cohabitating","Single","Other"))

for (v in names(MoBa1)) MoBa1[is.na(MoBa1[,v]),v] = NA
############### organize ALC data ################################

code_alc = function(df){
  df = data.frame(df)
  vnms = names(df)
  ever_vs = vnms[grep("ALC_EVER",vnms)]
  df[which(df[,ever_vs] == 0),ever_vs] = NA
  df[,ever_vs] = factor(df[,ever_vs],labels = c("No","Yes"))
  never = df[,ever_vs] == "No"
  
  freq_vs = vnms[grep("ALC_FREQ",vnms)]
  binge_vs = vnms[grep("ALC_BINGE",vnms)]
  typun_vs = vnms[grep("ALC_TYPUNITS",vnms)]

  
  for (v in freq_vs){
    FREQ = table((df[,v]))
    names(FREQ) = names(attributes(df[,v])[["labels"]])
    df[,v][which(df[,v] == attributes(df[,v])[["labels"]]["More than 1 check box filled in"])] = NA
    
    df[,v] = max(df[,v],na.rm = T)-df[,v]+1
    df[,v] = factor(df[,v],labels = c("never","<1/month","1-3/month","1/week","2-3/week","4-5/week","6-7/week"),ordered = T)
    df[which(never),v] = "never"
    oa = attributes(df[,v])
    oa[["orig_labels_freq"]] = FREQ
    attributes(df[,v]) = oa
    
    nvar = paste(v,"_num",sep = "")
    df[,nvar] = mapvalues(as.numeric(df[,v]),1:7,c(0,.5,2,4,2.5*4,4.5*4,6.5*4))
    df[which(df[,nvar] > 0),ever_vs] = "Yes"
  }
  
  for (v in binge_vs){
    FREQ = table((df[,v]))
    names(FREQ) = names(attributes(df[,v])[["labels"]])
    df[,v][which(df[,v] == attributes(df[,v])[["labels"]]["More than 1 check box filled in"])] = NA
    
    df[,v] = max(df[,v],na.rm = T)-df[,v]
    df[,v] = factor(df[,v],labels = c("never","<1/month","1-3/month","1/week","sev.time/week"),ordered = T)
    df[which(never),v] = "never"
    oa = attributes(df[,v])
    oa[["orig_labels_freq"]] = FREQ
    attributes(df[,v]) = oa
    
    nvar = paste(v,"_num",sep = "")
    df[,nvar] = mapvalues(as.numeric(df[,v]),1:5,c(0,.5,2,4,16))
    df[which(df[,nvar] > 0),ever_vs] = "Yes"  
  }
  
  for (v in typun_vs){
    FREQ = table((df[,v]))
    names(FREQ) = names(attributes(df[,v])[["labels"]])
    df[,v][which(df[,v] == attributes(df[,v])[["labels"]]["More than 1 check box filled in"])] = NA
    
    df[,v] = max(df[,v],na.rm = T)-df[,v]+1
    df[which(never),v] = 0
    df[,v] = factor(df[,v],labels = c("0","<1","1-2","3-4","5-6","7-9",">=10"),ordered = T)
    oa = attributes(df[,v])
    oa[["orig_labels_freq"]] = FREQ
    attributes(df[,v]) = oa
    
    nvar = paste(v,"_num",sep = "")
    df[,nvar] = mapvalues(as.numeric(df[,v]),1:7,c(0,.5,1.5,3.5,5.5,8,12))
    df[which(df[,nvar] > 0),ever_vs] = "Yes"
  }
  return(df)
}

MoBa1 = code_alc(MoBa1)
MoBa1$mALC_EVER[ MoBa1$mALC_TACE_EFFUNITS > 0 ] = "Yes"
MoBa1$mALC_TACE_EFFUNITS[MoBa1$mALC_EVER == "No"] = 0


avs1 = names(MoBa1)[c(grep("mALC_TACE",names(MoBa1)),grep("mALC_RAPI",names(MoBa1)))[-1]]
for (v in avs1){MoBa1[which(MoBa1[,v] == 0),v] = NA}
MoBa1$mALC_EVER[rowSums(MoBa1[,avs1],na.rm = T) > 0] = "Yes"
never = which(MoBa1$mALC_EVER == "No")
for (a in avs1){
  MoBa1[,v][MoBa1[,v] == attributes(MoBa1[,v])[["labels"]]["More than 1 check box filled in"]] = NA
}


############### make summary score fore each drug ################################
drug_levels = c("Never","Prev","LstMbfPreg","DurPreg")
for (l in 1:length(drug_levels)){
  idx = grep(drug_levels[l],names(MoBa1))
  MoBa1[,idx] = MoBa1[,idx]*(l-1)
}
drug = c("Hash","Amph","Ecst","Coc","Heroine")
for (d in 1:length(drug)){
  idx = grep(drug[d],names(MoBa1))
  vname = paste("mDU_sc_",drug[d],sep = "")
  MoBa1[,vname] = rowMaxs(as.matrix(MoBa1[,idx]),na.rm = T)
  #MoBa1[,vname] = apply(MoBa1[,idx],1,function(x) max(x,na.rm = T))
  MoBa1[abs(MoBa1[,vname]) > 10 ,vname] = NA
  MoBa1 = MoBa1[,grep(paste("mDU_",drug[d],sep = ""),names(MoBa1))*-1]
  MoBa1[,vname] = factor(MoBa1[,vname],labels = drug_levels,ordered = T)
}


vs = c("mLTHD1", "mLTHD2","mLTHD3","mLTHD4","mLTHD5","mLTHD3plus",
       "mSCL_Fear_Q1", "mSCL_Nervous_Q1", "mSCL_Hopeless_Q1", "mSCL_Blue_Q1", "mSCL_Worry_Q1",
       "mFeel_Pos_Q1", "mFeel_Useless_Q1", "mFeel_NotProud_Q1", "mFeel_Valuabel_Q1")
for (v in vs){
  MoBa1[,v][MoBa1[,v] == attributes(MoBa1[,v])[["labels"]]["More than 1 check box filled in"]] = NA
}

MoBa1$pMoBa1 = 1
RD = merge(RD,MoBa1,by = "PID", all.x = T, all.y = T)
remove(nms,MoBa1,a,avs1,c,d,drug,drug_levels,fAge,idx,l,never,v,vname)

RD$mpreWght1[RD$mpreWght1 > 150] = NA
RD$mpreWght1[RD$mpreWght1 < 35] = NA
RD$mpreWght2[RD$mpreWght2 > 150] = NA
RD$mpreWght2[RD$mpreWght2 < 35] = NA


Q3 = read_sav("F:/Forskningsprosjekter/PDB 2100 - The importance of pr_/Forskningsfiler/GUBI/PDB2014_Q3_v9.sav");
names(Q3) = c("PID","x","y","z","mpreWght3")
Q3$mpreWght3[Q3$mpreWght3 > 150] = NA
Q3$mpreWght3[Q3$mpreWght3 < 35] = NA
RD = merge(RD,Q3[,c("PID","mpreWght3")],by = "PID", all.x = T, all.y = T)
RD$mWght = RD$mpreWght1
RD$mBMI = RD$mWght / ((RD$mTall/100)^2)
RD$mBMI = RD$mWght / ((RD$mTall/100)^2)

RD$fWght[RD$fWght > 200] = NA
RD$fWght[RD$fWght < 30] = NA

RD$mTall[RD$mTall < 125] = NA
RD$fTall[RD$fTall < 125] = NA

RD$low_b_weight = (RD$weight < 2500)*1
RD$low_b_weight = factor(RD$low_b_weight,labels = c("No","Yes"))

RD$mINC = ordered(RD$mINC,levels = 1:7, labels = c("No", "<150", "150-199", "200-299", "300-399", "400-499", ">500"))
RD$fINC = ordered(RD$fINC,levels = 1:7, labels = c("No", "<150", "150-199", "200-299", "300-399", "400-499", ">500"))
RD$mINC_num = mapvalues(as.numeric(RD$mINC),1:7,c(0,7.5,17.5,25,35,45,65))
RD$MarSt = factor(RD$MarSt,labels = c("cohabitating/married","single"))
RD$Parity = ordered(RD$Parity,labels = c("0","1","2","3","4"))
RD$MULTIVU = factor(RD$MULTIVU,labels = c("No","Yes"))
RD$FOLATU = factor(RD$FOLATU,labels = c("No","Yes"))

EDUlabels = c("basic","1-2 high school","vocational high","general high", "Bachelor","Master" )


#RD$mEDUcomp = as.numeric(cut(RD$mEDUcomp,breaks = c(0,1.5,2.5,4.5,5.5,6.5)))
#RD$mEDUcur = as.numeric(cut(RD$mEDUcur,breaks = c(0,1.5,2.5,4.5,5.5,6.5)))
#tmp = RD[,c("mEDUcomp", "mEDUcur")]
#tmp[is.na(tmp)] = 0
#RD$mEDU = apply(tmp,1,max)
#RD$mEDU[rowSums(is.na(RD[,c("mEDUcomp", "mEDUcur")]))==2] = NA
RD$mEDU = RD$mEDUcomp
idx = is.na(RD$mEDUcomp) & !is.na(RD$mEDUcur) & RD$mEDUcur != 1
RD$mEDU[idx] = RD$mEDUcur[idx] - 1
RD$mEDU = ordered(RD$mEDU,labels = EDUlabels)

tmp = cbind(mapvalues(RD$mEDUcomp,1:6,c(9,10.5,12,12,15.5,18)),
            mapvalues(RD$mEDUcur,1:6,c(8,9.75,11,11,13.75,16))) # c(8,9.75,10.5,10.5,13.75,15)
RD$mEDUyrs = rowMaxs(tmp,na.rm = T)
RD$mEDUyrs[RD$mEDUyrs<0] = rep(NA,sum(RD$mEDUyrs<0))
tgts = sort(c(9,10.5,12,12,15.5,18,8,9.75,11,11,13.75,16))
for (k in (which((RD$mAge-RD$mEDUyrs)<5))){
  RD$mEDUyrs[k] = tgts[max(which((RD$mAge[k]-tgts-5)>=0))]
}

#RD$fEDUcomp = as.numeric(cut(RD$fEDUcomp,breaks = c(0,1.5,2.5,4.5,5.5,6.5)))
#RD$fEDUcur = as.numeric(cut(RD$fEDUcur,breaks = c(0,1.5,2.5,4.5,5.5,6.5)))
#tmp = RD[,c("fEDUcomp", "fEDUcur")]
#tmp[is.na(tmp)] = 0
#RD$fEDU = apply(tmp,1,max)
#RD$fEDU[rowSums(is.na(RD[,c("fEDUcomp", "fEDUcur")]))==2] = NA
RD$fEDU = RD$fEDUcomp
idx = is.na(RD$fEDUcomp) & !is.na(RD$fEDUcur) & RD$fEDUcur != 1
RD$fEDU[idx] = RD$fEDUcur[idx] - 1
RD$fEDU = ordered(RD$fEDU,labels = EDUlabels)

tmp = cbind(mapvalues(RD$fEDUcomp,1:6,c(9,10.5,12,12,15.5,18)),
             mapvalues(RD$fEDUcur,1:6,c(8,9.75,11,11,13.75,16)))


# add variable CN to cases where MoBa1 but not MFR was available
if (length(unique(RD[which(is.na(RD$CN)),]$PID)) == length(which(is.na(RD$CN)))) {RD$CN[is.na(RD$CN)] = 1}

#####################################################
######## add data smoking data ############
#####################################################
code_smoke = function(df){
  df$SMOK_EVER[df$SMOK_EVER < 1] = NA
  df$SMOK_EVER = factor(df$SMOK_EVER,labels = c("No","Yes"))
  for (v in c("SMOK_befPreg","SMOK_aftPreg")){
    df[,v] = factor(df[,v]-1,levels = 0:2, labels = c("No","Sometimes","Daily"),ordered = T)
    idx = which(is.na(df[,v]) & df$SMOK_EVER == "No")
    df[idx,v] = "No"
  }
  
  vars = c("SMOK_befPreg_perDay","SMOK_befPreg_perWeek","SMOK_aftPreg_perDay","SMOK_aftPreg_perWeek")
  for (v in vars){
    idx = which(is.na(df[,v]) & df$SMOK_EVER == "No")
    df[idx,v] = 0
  }
  df$SMOK_befPreg_perDay[which(df$SMOK_befPreg == "No" & is.na(df$SMOK_befPreg_perDay))] = 0
  df$SMOK_aftPreg_perDay[which(df$SMOK_aftPreg == "No" & is.na(df$SMOK_aftPreg_perDay))] = 0
  
  df$SMOK_befPreg_DailAvg = df$SMOK_befPreg_perDay*.5 + df$SMOK_befPreg_perWeek/7*.5
  df$SMOK_befPreg_DailAvg[is.na(df$SMOK_befPreg_perDay)] = df$SMOK_befPreg_perWeek[is.na(df$SMOK_befPreg_perDay)]/7
  df$SMOK_befPreg_DailAvg[is.na(df$SMOK_befPreg_perWeek)] = df$SMOK_befPreg_perDay[is.na(df$SMOK_befPreg_perWeek)]
  
  df$SMOK_aftPreg_DailAvg = df$SMOK_aftPreg_perDay*.5 + df$SMOK_aftPreg_perWeek/7*.5
  df$SMOK_aftPreg_DailAvg[is.na(df$SMOK_aftPreg_perDay)] = df$SMOK_aftPreg_perWeek[is.na(df$SMOK_aftPreg_perDay)]/7
  df$SMOK_aftPreg_DailAvg[is.na(df$SMOK_aftPreg_perWeek)] = df$SMOK_aftPreg_perDay[is.na(df$SMOK_aftPreg_perWeek)]
  
  df$SMOK_aftPreg_DailAvg[which((df$SMOK_aftPreg_DailAvg > df$SMOK_befPreg_DailAvg) & is.na(df$SMOK_aftPreg_perDay))] = NA
  df$SMOK_aftPreg_DailAvg[which((df$SMOK_aftPreg_DailAvg > df$SMOK_befPreg_DailAvg) & is.na(df$SMOK_befPreg_perDay))] = NA
  df$SMOK_befPreg_DailAvg[which((df$SMOK_aftPreg_DailAvg > df$SMOK_befPreg_DailAvg) & is.na(df$SMOK_aftPreg_perDay))] = NA
  df$SMOK_befPreg_DailAvg[which((df$SMOK_aftPreg_DailAvg > df$SMOK_befPreg_DailAvg) & is.na(df$SMOK_befPreg_perDay))] = NA
  
  df$SMOK_aftPreg_DailAvg[which((df$SMOK_aftPreg_DailAvg > df$SMOK_befPreg_DailAvg) & (df$SMOK_aftPreg_perDay < df$SMOK_befPreg_perDay))] = 
    df$SMOK_aftPreg_perDay[which((df$SMOK_aftPreg_DailAvg > df$SMOK_befPreg_DailAvg) & (df$SMOK_aftPreg_perDay < df$SMOK_befPreg_perDay))]
  df$SMOK_befPreg_DailAvg[which((df$SMOK_aftPreg_DailAvg > df$SMOK_befPreg_DailAvg) & (df$SMOK_aftPreg_perDay < df$SMOK_befPreg_perDay))] = 
    df$SMOK_befPreg_perDay[which((df$SMOK_aftPreg_DailAvg > df$SMOK_befPreg_DailAvg) & (df$SMOK_aftPreg_perDay < df$SMOK_befPreg_perDay))]
  
  df$SMOK_EVER[which(df$SMOK_befPreg_DailAvg > 0 | df$SMOK_aftPreg_DailAvg > 0)] = "Yes"
  #df$SMOK_befPreg_perDay[df$SMOK_befPreg_perDay > 40] = 40
  #df = df[,-grep("Week",names(df))]
  return(df)
}

load("F:/Forskningsprosjekter/PDB 2100 - The importance of pr_/Forskningsfiler/GUBI/SMOK_PDB2014_Q1_v9.Rdata")
names(M1smoke) = c("PID", "fSMOK_befPreg","fSMOK_aftPreg","SMOK_EVER", "SMOK_aftPreg", "SMOK_aftPreg_perWeek","SMOK_aftPreg_perDay",
                   "SMOK_befPreg", "SMOK_befPreg_perWeek","SMOK_befPreg_perDay")
for (v in names(M1smoke)) M1smoke[is.na(M1smoke[,v]),v] = NA
M1smoke = data.frame(M1smoke)
M1smoke = code_smoke(M1smoke)
names(M1smoke)[-1] = paste("m",names(M1smoke)[-1],sep = "")
RD = merge(RD,M1smoke,by = "PID", all.x = T, all.y = T)
RD$mfSMOK_aftPreg = factor(RD$mfSMOK_aftPreg,levels = c(1,2), labels = c("No","Yes"))
RD$mfSMOK_befPreg = factor(RD$mfSMOK_befPreg,levels = c(1,2), labels = c("No","Yes"))


######## if possible correct MFR weight/length data with help of MoBa questionnaire 4 ################
load("F:/Forskningsprosjekter/PDB 2100 - The importance of pr_/Forskningsfiler/GUBI/PDB2014_Q4_6mnd_v9_red.sav")
for (v in names(MoBa4)) MoBa4[is.na(MoBa4[,v]),v] = NA
MoBa4 = data.frame(MoBa4)

names(MoBa4) = c("PID","CN","MBsex","MBweight","MBlength","MBweek")
MoBa4$MBlength = round(MoBa4$MBlength,digits = 0)
MoBa4$MBlength[MoBa4$MBlength > 65] = NA
MoBa4$MBlength[MoBa4$MBlength < 20] = NA
MoBa4$MBweight[MoBa4$MBweight > 6500] = NA
MoBa4$MBweight[MoBa4$MBweight < 500] = NA
RD = merge(RD,MoBa4,by = c("PID","CN"), all.x = T, all.y = T)
RD$length[is.na(RD$length)] = round(RD$MBlength[is.na(RD$length)],digits = 0)
RD$weight[is.na(RD$weight)] = round(RD$MBweight[is.na(RD$weight)],digits = 0)
RD$MBlength[is.na(RD$MBlength)] = round(RD$length[is.na(RD$MBlength)],digits = 0)
RD$MBweight[is.na(RD$MBweight)] = round(RD$weight[is.na(RD$MBweight)],digits = 0)

tdf = RD[,c("weight","length")]
#tdf[which(tdf$weight > 2000 & tdf$length < 40),] = NA
tdf = log(tdf)
mW = lm(weight ~length,data=tdf)
mL = lm(length ~weight,data=tdf)
sdW = sd(RD$weight,na.rm=T)
sdL = sd(RD$length,na.rm=T)


x = RD$length-RD$MBlength
y = RD$weight-RD$MBweight
large_devs = which((abs(x) > 5 | abs(y) > 500) & (!is.na(x) & !is.na(y)))
cols = c("weight","length","weight","MBlength","MBweight","MBlength","MBweight","length")
for (id in large_devs) {
  #print(RD[id,c("weight","length","MBweight","MBlength")])
  wlm = t(matrix(as.numeric(RD[id,cols]),nrow = 2))
  
  pW = log(wlm)
  pW[,1] = 1
  predWeight = exp(mW$coefficients %*% t(pW))
  
  pL = log(wlm[,c(2,1)])
  pL[,1] = 1
  predLength = exp(mL$coefficients %*% t(pL))
  
  pred = t(as.matrix(rbind(predWeight,predLength)))
  
  delta = pred-wlm
  normed_deviance = (delta[,2]/sdL)^2 + (delta[,1]/sdW)^2
  tmp_id = min(which(normed_deviance == min(normed_deviance)))
  
  RD$length[id] = wlm[tmp_id,2]
  RD$weight[id] = wlm[tmp_id,1]
  #print(RD[id,c("weight","length")])
}


RD$cBMI = (RD$weight/1000)/((RD$length/100)^2)
RD$cBMI[RD$cBMI < 5] = NA
RD$cBMI[RD$cBMI >25] = NA


RD$length[which(RD$length < 10)] = NA
RD$length[which(RD$length > 65)] = NA

tmp = RD[,c("weeks_UL","weeks_MS")]

RD$weeks = rowMeans(RD[,c("weeks_MS","weeks_UL")],na.rm = T)
RD$weeks_transformed = log(48-RD$weeks)
attributes(RD$weeks_transformed) = list("transform" = function(x){log(48-x)},"inverse_transform" = function(x) {48-exp(x)})

RD$Preterm = (RD$weeks < 37)*1
RD$Preterm = factor(RD$Preterm,labels = c("No","Yes"))

RD$weeks[is.na(RD$weeks)] = RD$MBweek[is.na(RD$weeks)]
RD$cSEX[is.na(RD$cSEX)] = RD$MBsex[is.na(RD$cSEX)]
RD$cSEX[!(RD$cSEX == 1 | RD$cSEX == 2)] = NA
RD$cSEX = factor(as.numeric(RD$cSEX))


#####################################################
##### add data from 6th Questionnaire (3yrs )########
#####################################################
# for included (Jan 2004 - Dec 2007), ASRS[R68], Depression [R72], Anxiety [R70]
MoBa6 = read_sav("F:/Forskningsprosjekter/PDB 2100 - The importance of pr_/Forskningsfiler/GUBI/PDB2014_Q6_3aar_v9.sav");
for (v in names(MoBa6)) MoBa6[is.na(MoBa6[,v]),v] = NA
# remove participants outside sampling period
#MoBa6 = select.cases.by.ID(PID,MoBa6)

cADHDidx = paste("GG",c(314:316,320,327,332,339:342,345),sep = "")
mADHDidx = paste("GG",503:508,sep = "")
DEPRidx = paste("GG",600:605,sep = "")
ANXIidx = paste("GG",514:521,sep = "")
Cidx = c("PREG_ID_2014","BARN_NR",mADHDidx,DEPRidx,ANXIidx,"ALDERUTFYLT_S6",cADHDidx)
MoBa6 = MoBa6[,Cidx]
names(MoBa6) = c("PID","CN",paste("mADHD",1:6,sep = ""),paste("mDEP",1:6,sep = ""),paste("mANX",1:8,sep = ""),"dAGE",paste("cADHD",1:11,sep = ""))
for (v in 3:dim(MoBa6)[2]){MoBa6[which(MoBa6[,v] ==0),v] = NA}
MoBa6$pMoBa6 = 1

RD = merge(RD,MoBa6,by = c("PID","CN"), all.x = T, all.y = T)
remove(MoBa6,cADHDidx,mADHDidx,DEPRidx,ANXIidx,Cidx)

###########################################################################
################### General Illnesses at age 3 ############################
###########################################################################

MoBa6HP =  read_sav("F:/Forskningsprosjekter/PDB 2100 - The importance of pr_/Forskningsfiler/GUBI/PDB2014_Q6_3aar_v9.sav");
cHPidx = c("PREG_ID_2014","BARN_NR",paste("GG",29:116,sep = ""))

HPS = c("Hearing","Vision","Motor","CerebralPalsy","Joints","Diabetes","ToLight","ToHeavy","Heart","Testicle",
        "Asthma","Allergi","AtopicEczema","OtherEczema","FoodAllergi","Gastrointestinal",
        "LateSpeech","Sleep","Autims","Hyperactivity","OtherBeh","OtherLT")
HPSstatus = c("No","Now","Previously","Specialist")
nms = c("PID","CN")
for (hp in HPS){nms = c(nms,paste(hp,HPSstatus,sep = "_"))}
MoBa6HP = MoBa6HP[,cHPidx]
names(MoBa6HP) = nms
for (hp in HPS){
  hpx = grep(hp,names(MoBa6HP))
  vn = paste("HP3",hp,sep = "_")
  MoBa6HP[,vn] = NA
  MoBa6HP[which(MoBa6HP[,hpx[1]] == 1),vn] = 0                              # no problem
  MoBa6HP[which(MoBa6HP[,hpx[3]] == 1 & MoBa6HP[,hpx[4]] == 1),vn] = 1      # previously / no specialist
  MoBa6HP[which(MoBa6HP[,hpx[3]] == 1 & MoBa6HP[,hpx[4]] == 2),vn] = 2      # previously / specialist
  MoBa6HP[which(MoBa6HP[,hpx[3]] == 1 & is.na(MoBa6HP[,hpx[4]])),vn] = 1.5  # previously / no repsonse to specialist questions
  MoBa6HP[which(MoBa6HP[,hpx[2]] == 1 & MoBa6HP[,hpx[4]] == 1),vn] = 3      # now / no specialist
  MoBa6HP[which(MoBa6HP[,hpx[2]] == 1 & MoBa6HP[,hpx[4]] == 2),vn] = 4      # now / specialist
  MoBa6HP[which(MoBa6HP[,hpx[2]] == 1 & is.na(MoBa6HP[,hpx[4]])),vn] = 3.5  # previously / no repsonse to specialist questions
}
MoBa6HP = MoBa6HP[,c(1,2,grep("HP3_",names(MoBa6HP)))]
#MoBa6HPl = melt(MoBa6HP,id.vars = c("PID","CN"))
#names(MoBa6HPl) = c("PID","CN","Disease","value")
#MoBa6HPl$value = factor(MoBa6HPl$value,labels = c("No","PrevNoSpec","Prev","PrevSpec","NowNoSpec","Now","NowSpec"),ordered = T)

RD = merge(RD,MoBa6HP,by = c("PID","CN"), all.x = T, all.y = T)

#####################################################
######## add data from Fathers Questionnaire ########
#####################################################


MoBaF = read_sav("F:/Forskningsprosjekter/PDB 2100 - The importance of pr_/Forskningsfiler/GUBI/PDB2014_QF_v9.sav");
#MoBaF = select.cases.by.ID(PID,MoBaF)
Cidx = c("PREG_ID_2014","FF16","FF17","FF341","FF386","FF389","FF392","FF395",
         paste("FF",535:540,sep = ""),
         paste("FF",214:221,sep = ""), #smoking
         paste("FF",222:241,sep = ""), #drugs
         paste("FF",553:556,sep = ""), #drugs - amph
         paste("FF",c(242:245,465:468,473,474),sep = "")) #alcohol

nms = c("PID","fEDUcomp","fEDUcur","fINC","fADHDdiag","fAnorex","fManDepr","fSchiz",
        paste("fANX",1:6,sep = ""),
        "fSMOK_EVER", "fSMOK_befPreg", "fSMOK_befPreg_perDay","fSMOK_befPreg_perWeek",
        "fSMOK_aftPreg", "fSMOK_aftPreg_perDay","fSMOK_aftPreg_perWeek","fSMOK_loc",
        paste("fDU_Hash",c("Never","Prev","LstMbfPreg","Recently"),sep = "_"),
        paste("fDU_Ecst",c("Never","Prev","LstMbfPreg","Recently"),sep = "_"),
        paste("fDU_Coc",c("Never","Prev","LstMbfPreg","Recently"),sep = "_"),
        paste("fDU_Heroine",c("Never","Prev","LstMbfPreg","Recently"),sep = "_"),
        paste("fDU_CentStim",c("Never","Prev","LstMbfPreg","Recently"),sep = "_"),
        paste("fDU_Amph",c("Never","Prev","LstMbfPreg","Recently"),sep = "_"),
        "fALC_EVER", "fALC_FREQ_befPreg", "fALC_FREQ_aftPreg","fALC_TYPUNITS",
        "fALC_TYPUNITS_befPreg_we","fALC_TYPUNITS_befPreg_w","fALC_TYPUNITS_aftPreg_we","fALC_TYPUNITS_aftPreg_w",
        "fALC_BINGE_befPreg","fALC_BINGE_aftPreg")
MoBaF = MoBaF[,Cidx]
names(MoBaF) = nms

for (v in names(MoBaF)) MoBaF[is.na(MoBaF[,v]),v] = NA
MoBaF = data.frame(MoBaF)

MoBaF$fINC[which(MoBaF$fINC == 0)] = NA
MoBaF$fINC = ordered(MoBaF$fINC,levels = 1:7, labels = c("No", "<150", "150-199", "200-299", "300-399", "400-499", ">500"))


# MoBaF$fEDUcomp = as.numeric(cut(MoBaF$fEDUcomp,breaks = c(0,1.5,2.5,4.5,5.5,6.5)))
# MoBaF$fEDUcur = as.numeric(cut(MoBaF$fEDUcur,breaks = c(0,1.5,2.5,4.5,5.5,6.5)))
# tmp = MoBaF[,c("fEDUcomp", "fEDUcur")]
# tmp[is.na(tmp)] = 0
# MoBaF$fEDU = apply(tmp,1,max)
# MoBaF$fEDU[rowSums(is.na(MoBaF[,c("fEDUcomp", "fEDUcur")]))==2] = NA
# MoBaF$fEDU = ordered(MoBaF$fEDU,labels = c("secondary","high1-2","high3+","bachelor","Master"))

MoBaF$fEDU = MoBaF$fEDUcomp
idx = is.na(MoBaF$fEDUcomp) & !is.na(MoBaF$fEDUcur) & MoBaF$fEDUcur != 1
MoBaF$fEDU[idx] = MoBaF$fEDUcur[idx] - 1
MoBaF$fEDU = ordered(MoBaF$fEDU,labels = EDUlabels)


for (v in c("fADHDdiag","fAnorex","fManDepr","fSchiz")){
  MoBaF[is.na(MoBaF[,v]),v] = 0
  MoBaF[,v] = factor(MoBaF[,v],labels = c("No","Yes"))
}

################## correct scale for ALC ##############

MoBaF = code_alc(MoBaF)

##################### smoking #########################
fSMOK = MoBaF[,c(1, grep("SMOK",names(MoBaF)))]
MoBaF = MoBaF[,-grep("SMOK",names(MoBaF))]
names(fSMOK) = gsub("fSMOK_","SMOK_", names(fSMOK))

fSMOK$SMOK_befPreg[fSMOK$SMOK_befPreg > 3] = NA
fSMOK$SMOK_befPreg[fSMOK$SMOK_befPreg < 1] = NA
fSMOK$SMOK_aftPreg[fSMOK$SMOK_aftPreg > 3] = NA
fSMOK$SMOK_aftPreg[fSMOK$SMOK_aftPreg < 1] = NA

M1smoke$mfSMOK_befPreg[M1smoke$mfSMOK_befPreg < 1] = NA
M1smoke$mfSMOK_aftPreg[M1smoke$mfSMOK_aftPreg < 1] = NA

#### try to add info missing from fathers questionnair by using mothers response ###
fSMOK = merge(M1smoke,fSMOK,by = "PID", all.x = T, all.y = T)
fSMOK$SMOK_EVER[which(is.na(fSMOK$SMOK_EVER) & fSMOK$mfSMOK_befPreg == 1)] = 1
fSMOK$SMOK_EVER[which(is.na(fSMOK$SMOK_EVER) & fSMOK$mfSMOK_befPreg == 2)] = 2
fSMOK$SMOK_befPreg[which(is.na(fSMOK$SMOK_befPreg) & fSMOK$mfSMOK_befPreg == 1)] = 1
fSMOK$SMOK_befPreg[which(is.na(fSMOK$SMOK_befPreg) & fSMOK$mfSMOK_befPreg == 2)] = 2
fSMOK$SMOK_aftPreg[which(is.na(fSMOK$SMOK_aftPreg) & fSMOK$mfSMOK_aftPreg == 1)] = 1
fSMOK$SMOK_aftPreg[which(is.na(fSMOK$SMOK_aftPreg) & fSMOK$mfSMOK_aftPreg == 2)] = 2
fSMOK = fSMOK[,-grep("^m",names(fSMOK))]
fSMOK = code_smoke(fSMOK)
names(fSMOK)[-1] = paste("f",names(fSMOK)[-1],sep = "")
RD = merge(RD,fSMOK,by = "PID", all.x = T, all.y = T)

############### make summary score fore each drug ################################
drug_levels = c("Never","Prev","LstMbfPreg","Recently")
for (l in 1:length(drug_levels)){
  idx = grep(drug_levels[l],names(MoBaF))
  MoBaF[,idx] = MoBaF[,idx]*(l-1)
}
# drug = c("Hash","CentStim","Ecst","Coc","Heroine")
# for (d in 1:length(drug)){
#   idx = grep(drug[d],names(MoBaF))
#   vname = paste("fDU_sc_",drug[d],sep = "")
#   MoBaF[,vname] = rowMins(MoBaF[,idx],na.rm = T)
#   MoBaF[MoBaF[,vname]>10,vname] = NA
#   MoBaF = MoBaF[,grep(paste("fDU_",drug[d],sep = ""),names(MoBaF))*-1]
#   MoBaF[,vname] = factor(MoBaF[,vname],labels = drug_levels,ordered = T)
# }

drug = c("Hash","CentStim","Ecst","Coc","Heroine","Amph")
for (d in 1:length(drug)){
  idx = grep(drug[d],names(MoBaF))
  vname = paste("fDU_sc_",drug[d],sep = "")
  MoBaF[,vname] = rowMins(as.matrix(MoBaF[,idx]),na.rm = T)
  MoBaF[is.infinite(MoBaF[,vname]),vname] = rep(NA,sum(is.infinite(MoBaF[,vname])))
  MoBaF = MoBaF[,grep(paste("fDU_",drug[d],sep = ""),names(MoBaF))*-1]
  MoBaF[,vname] = factor(MoBaF[,vname],labels = drug_levels,ordered = T)
}

# merge the 2 questions about amphetamines / central stimulants
MoBaF$fDU_sc_Amph[is.na(MoBaF$fDU_sc_Amph)] = MoBaF$fDU_sc_CentStim[is.na(MoBaF$fDU_sc_Amph)]
MoBaF = MoBaF[,-grep("fDU_sc_CentStim",names(MoBaF))]

MoBaF$pMoBaF = 1


RD = merge(RD,MoBaF,by = c("PID"), all.x = T, all.y = T)

#### harmonize data about fathers' education ###
# set value to fathers response
RD$fEDU = RD$fEDU.y
# fill missings with mothers response
RD$fEDU[is.na(RD$fEDU)] = RD$fEDU.x[is.na(RD$fEDU)]

#### harmonize data about fathers' income ###
RD$fINC = RD$fINC.y
# fill missings with mothers response
RD$fINC[is.na(RD$fINC)] = RD$fINC.x[is.na(RD$fINC)]
RD$fINC_num = mapvalues(as.numeric(RD$fINC),1:7,c(0,7.5,17.5,25,35,45,65))

remove(MoBaF,Cidx,nms)

RD$pMFR[is.na(RD$pMFR)] = 0
RD$pMoBa1[is.na(RD$pMoBa1)] = 0
RD$pMoBa6[is.na(RD$pMoBa6)] = 0
RD$pMoBaF[is.na(RD$pMoBaF)] = 0

################### miscellaneous info #####################################

ps = c("pMFR", "pMoBa1", "pMoBa6", "pMoBaF")
RD[,ps] = RD[,ps] == 1


RD$MOBA0 = RD$pMoBa1
RD$MOBA3 = RD$pMoBa6


RD$group = NA
RD$group[RD$MOBA0] = 0
RD$group[RD$MOBA3] = 1


for (v in names(RD)) {
  if (is.factor(RD[,v])) {
    if (sum(levels(RD[,v]) == "NaN") > 0) {
      new_levels = levels(RD[,v])[-which(levels(RD[,v]) == "NaN")]
      RD[!(RD[,v] %in% new_levels),v] = NA
      print(paste0(paste(new_levels, collapse = ", "),"  ::  ", paste(levels(RD[,v]), collapse = ", ")))
      RD[,v] = factor(RD[,v], levels = new_levels)
    }
  }
}

save(RD,file="RDallnew_v9.Rdata")


idx = c(which(RD$byear > 2000 & RD$living == "Yes"))
RD = RD[idx,]
rownames(RD) = 1:nrow(RD)

save(RD,file="RDallnew_v9_2001_2009.Rdata")

### Add fathers ADHD scores

load("../IPW/Data/RDallnew_v9_2001_2009.Rdata")
RD = RD[,-grep("fADHD[0-9]",names(RD))]
MoBaF = read_sav("F:/Forskningsprosjekter/PDB 2100 - The importance of pr_/Forskningsfiler/GUBI/PDB2014_QF_v9.sav");
MoBaF = MoBaF[,c("PREG_ID_2014",paste0("FF",535:540))]
names(MoBaF) = c("PID",paste0("fADHD",1:6))
MoBaF = data.frame(MoBaF)
for (k in 2:7) {
  tmp = as.numeric(MoBaF[,k]-1)
  tmp[tmp < 0] = NA
  MoBaF[,k] = tmp
}

RD = merge(RD,MoBaF,by = "PID", all.x = T, all.y = F)

v = paste0("fADHD",1:6)
fADHD = rowMeans(RD[,v],na.rm = T)
fADHD[rowSums(is.na(RD[,v])) > 2] = NA
v = paste0("mADHD",1:6)
mADHD = rowMeans(RD[,v],na.rm = T)
mADHD[rowSums(is.na(RD[,v])) > 2] = NA

save(RD,file = "../IPW/Data/RDallnew_v9_2001_2009.Rdata")
