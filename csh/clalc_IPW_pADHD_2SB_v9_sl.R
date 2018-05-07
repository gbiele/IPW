options(mc.cores = 4)
library(data.table)
library(rstanarm)
library(boot)


args=(commandArgs(TRUE))

i = args[1]

load("/usit/abel/u1/guidopb/programming/R/representativeness/IPW/data/IPW_DATA_pADHD_v9_2SB.Rdata")

ipw_data = ipw_data[imp == i]

ipw_data[,ageNs := scale(ageN)]
ipw_data[,ageN2s := scale(ageN2)]



mx2 = stan_glmer(cbind(N.moba, N.ssb - N.moba) ~ 1 + ageNs + ageN2s + (ageNs + ageN2s | edu:parity),
                 data = ipw_data,  family = binomial(link = "logit"), iter = 6000, warmup = 4000,adapt_delta = .9999,
                 control = list(max_treedepth = 20))


k = 0
test = T
while(test) {
  pp = posterior_predict(mx2x)
  test = sum(pp == 0) > 0
  k = k + 1
  if(k == 2000) test = F
}


for (g in 1:nrow(ipw_data)) pp[,g] = pp[,g]/ipw_data$N.ssb[g]
sipw = (sum(ipw_data[,N.moba])/sum(ipw_data[,N.ssb])) / t(pp)

save(sipw,ipw_data,mx2,file = paste0("/usit/abel/u1/guidopb/programming/R/representativeness/IPW/data/IPW_i",i,"_pADHD_2SB_v9.Rdata"))
