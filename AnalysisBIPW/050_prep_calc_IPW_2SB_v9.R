source("../analysis_utils.R")
library(rstan)
options(mc.cores = 4)
library(data.table)
# file IPW_DATA.Rdata generated in comp_SSBMOBAR.R
#sm = stan_model("stan/ipw_thetas.stan")
library(rstanarm)
library(boot)

clean_data = function(dt) {
  dt[Age == "<20", parity := 0]
  dt[Age == "<20" & (edu == "Bachelor" | edu == "Master"), edu := "High-school"]
  dt[Age == "20-24" & parity == 3, parity := 2]
  dt[Age == "20-24" & edu == "Master", edu := "Bachelor"]
  dt[Age == "20-24" & edu == "Bachelor" & parity == 2, parity := 1]
  dt[Age == "25-29" & edu == "Master" & parity > 1, parity := 1]
  dt[Age == "40-49" & edu == "Master" & parity == 3, parity := 2]
  dt[Age == "40-49" & edu == "Elementary" & parity == 0, parity := 1]
  return(dt)  
}
load("data/ssb_l.Rdata")
ssb_l_cleaned = clean_data(ssb_l)
load("data/imputed_preprocessed_cleaned_ADHDdata_pADHD_v9s.Rdata")
moba_l = imputed_data[,.(N = .N), by = c("Age","edu","parity","imp")]

moba_l$group = "moba"
ssb_l$group = "pop"
ssb_l = ssb_l[is.na(N), N := 0]

ipw_data = merge(moba_l[,list(N.moba = sum(N)), by = c("imp","Age","edu","parity")],
                 ssb_l[,list(N.ssb = sum(N)), by = c("Age","edu","parity")],
                 by = c("Age","edu","parity"),
                 all = T)

ipw_data[,P := round(N.ssb/sum(ipw_data$N.ssb,na.rm = T)*1000,digits = 1)]
ipw_data[imp == 1 & P < .5]


ipw_data[, ipw_direct := N.ssb/N.moba]
ipw_data[, sipw_direct := (N.ssb/N.moba)*sum(ipw_data$N.moba)/sum(ipw_data$N.ssb)]


ipw_data[,ageN := poly(as.numeric(ipw_data$Age), order = 3)[,1]]
ipw_data[,eduN := poly(as.numeric(ipw_data$edu), order = 2)[,1]]
ipw_data[,parityN := poly(as.numeric(ipw_data$parity), order = 2)[,1]]
ipw_data[,ageN2 := poly(as.numeric(ipw_data$Age), order = 3)[,2]]
ipw_data[,ageN3 := poly(as.numeric(ipw_data$Age), order = 3)[,3]]
ipw_data[,eduN2 := poly(as.numeric(ipw_data$edu), order = 2)[,2]]
ipw_data[,parityN2 := poly(as.numeric(ipw_data$parity), order = 2)[,2]]

ipw_data[, N.moba := round(N.moba)]
ipw_data[, N.ssb := round( N.ssb)]
ipw_data[, OP  := N.moba/N.ssb]


save(ipw_data,file = "data/IPW_DATA_pADHD_v9_2SB.Rdata")

