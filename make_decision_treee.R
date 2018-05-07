library(data.table)
library(rpart)
cs = expand.grid("L_P" = c(F,T),
                 "CC_DL" = c(F,T),
                 "E_P" = c(F,T),
                 "E_L" = c(F,T),
                 "L_ED" = c(F,T))
cs = data.table(cs)
setkeyv(cs,c("L_P","CC_DL"))
cs[L_P == F ,analysis := "no_bias"]
cs[E_P == T & CC_DL == T & L_P == T, analysis := "AR"]
cs[E_L == T & CC_DL == T & L_P == T, analysis := "IPW"]
cs[L_ED == T & L_P == T, analysis := "MRP"]
cs[E_L == T & CC_DL == T & L_P == T & L_ED == T, analysis := "IPW"]
cs[E_P == T & CC_DL == T & L_P == T & L_ED == T, analysis := "MRP"]
cs[is.na(analysis), analysis := "no_bias"]


tfit = rpart(analysis ~ L_P + E_P + E_L + L_ED + CC_DL, data=cs,method="class",minsplit=2, minbucket=1)
plot(tfit,uniform = T)
text(tfit, use.n=TRUE, all=TRUE, cex=.8)
post(tfit, filename = "", 
     title = "SB")

cs[analysis == "MRP", analysis := "IPW"]
tfit2 = rpart(analysis ~ L_P + E_P + E_L + L_ED + CC_DL, data=cs,method="class",minsplit=2, minbucket=1)
plot(tfit2,uniform = T)
text(tfit2, use.n=TRUE, all=TRUE, cex=.8)
post(tfit2, filename = "", 
     title = "SB")
