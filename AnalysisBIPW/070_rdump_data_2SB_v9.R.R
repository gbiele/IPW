load("data/analyses_pADHD_2SB_v9.Rdata")
library(data.table)
library(rstan)
my_data[cADHD > 22, cADHD := 22]
my_data_imps = my_data

rids = sample(100000:999999,1000)
k = 1

avg_sipw = c()

for (cimp in 11:20) {
  my_data = my_data_imps[imp == as.character(cimp)]
  
  load(paste0("data/IPW_i",cimp,"_pADHD_2SB_v9.Rdata"))
  ipw_data[,ipw_idx := 1:nrow(ipw_data)]
  ipw_data[,ipw_idx := 1:nrow(ipw_data)]
  
  my_data = merge(my_data,
                  ipw_data[,c("Age", "edu", "parity","ipw_idx"), with = F],
                  by = c("Age", "edu", "parity"),all = T)

  sipw = sipw[,colSums(sipw > max(apply(sipw,1,function(x) quantile(x,.995)))) == 0]  

  sipw_fixed_sum = 
    apply(sipw[,1:2000],2, function(x1) {
      tmp = x1[my_data$ipw_idx]
      tmpw = tmp*nrow(my_data)/sum(tmp)
      return( sapply(as.list(1:nrow(sipw)), function(x2) {
        tmpw[min(which(my_data$ipw_idx == x2))]
      } ))
    })
  avg_sipw = rbind(avg_sipw,rowMeans(sipw,na.rm = T))
  for (a in 1:3) {
    f = analyses[[a]]$formula
    
    standata = list(N = nrow(my_data),
                    Y = my_data$cADHD,
                    K = length(attr(terms(f),"term.labels")),
                    K_ipw = length(attr(terms(f),"term.labels")) - ifelse(analyses[[a]]$adjust_ipw,0,1)*length(ipw_vars),
                    X = as.matrix(model.matrix(f,my_data))[,-1],
                    ipw_idx = my_data$ipw_idx,
                    ipw_weights = sipw_fixed_sum,
                    id = rids[k])
    save(standata,file = paste0("data/4Stan/a",a,"_i",cimp,"x.Rdata"))
    k = k + 1
  }
}


