load("data/analyses_pADHD_2SB_v9.Rdata")
library(data.table)
library(rstan)
library(tikzDevice)
my_data[cADHD > 22, cADHD := 22]
my_data_imps = my_data

rids = sample(100000:999999,1000)
k = 1

exposures = c("small_fga", "preterm", 
              "mSMOKE", "mSMOKE.5Cs", "mALC.FREQ", "mALC.EFFUNITS", 
              "fSMOKE", "fSMOKE.5Cs", "fALC.FREQ", "fALC.TYPUNITS",
              "fDRUGS", "fDRUGS.SCORE", "mDRUGS", "mDRUGS.SCORE", "mLTHD", "mSCL5")

for (cimp in 1) {
  my_data = my_data_imps[imp == as.character(cimp)]
  
  load(paste0("data/IPW_i",cimp,"_pADHD_2SB_v9.Rdata"))
  ipw_data[,ipw_idx := 1:nrow(ipw_data)]

  my_data = merge(my_data,
                  ipw_data[,c("Age", "edu", "parity","ipw_idx"), with = F],
                  by = c("Age", "edu", "parity"),all = T)

  Pparticipation = 0.231411/apply(sipw,1,median)
  my_data[, Ppart := Pparticipation[ipw_idx]]
  
  median_sipw = apply(sipw,1,median)
  
  my_data[, sipw := median_sipw[ipw_idx]]
  
  # implied conditional indendency when L is collider: D ⊥ P | L
  ms1 = summary(lm(cADHD ~ Ppart + edu + Age + poly(parity,2) , my_data))
  # implied conditional indendency when L is not collider: D ⊥ P | E, L
  ms2 = summary(lm(cADHD ~ Ppart + edu + Age + poly(parity,2) +
                     small_fga + preterm + 
                     mDRUGS + mDRUGS.SCORE + mLTHD + mSCL5 + fDRUGS + fDRUGS.SCORE ,
                   my_data))
  rbind(ms1$coefficients["Ppart",],
        ms2$coefficients["Ppart",])
  
  cors = c()
  for (ex in exposures) {
    cors = rbind(cors,
                 wtd.cor(my_data$mEdu.cat,
                         my_data[[ex]],
                         weight = my_data$sipw))
  }
  
  cors = data.table(cors)
  cors$exposure = exposures
  cors[, lb := correlation - 1.96*cors$std.err]
  cors[, ub := correlation + 1.96*cors$std.err]
  tikz(file = "LTX/figures/mEDUcors.tex",width = 13/2.54, height = 15/2.54)
  par(mar=c(3,8,.5,1), mgp=c(2,.7,0), tck=-.01)
  y = barplot(cors$correlation,
              horiz = T,
              xlab = "IPP weighted correlation",
              cex.main = 1,border = NA,
              xlim = c(min(cors$ub)*1.1,max(cors$lb)*1.1))
  axis(2, at = y, labels = gsub("_",".",exposures), las = 2, lwd = 0, cex.axis = 1)
  segments(x0 = cors$ub,
           x1 = cors$lb,
           y0 = y)
  dev.off()
}


