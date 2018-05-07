library(HDInterval)
library(rstan)
library(data.table)
library(RColorBrewer)
library(bayesplot)
library(xtable)
library(gtools)

as_matrix = function(sf,var) {
  a = extract(sf,var,permuted = F);
  samples = c()
  for (chain in 1:dim(a)[2]) samples = rbind(samples,a[,chain,])
  return(samples)
}

get_stats = function(s, ROPE = NULL) {
  my_stats = paste0(ifelse(round(mean(s),digits = 2) < 0,""," "),
                    ifelse(abs(mean(s))<10,sprintf("%.2f",mean(s)),sprintf("%.0f",mean(s))),
                    " (",
                    paste(sprintf("%.2f",hdi(s,credMass = .9)),
                          collapse = ","),
                    ")")
  if (!is.null(ROPE)) my_stats = paste0(my_stats,
                                        "; ",
                                        sprintf("%.1f", logit(mean(s > min(ROPE) & s < max(ROPE)))))
  return(my_stats)
}

p_in_rope = function(s) {
  x = seq(.01,4,length = 100)
  y = sapply(x, function(x) mean(s > (x*-1) & s < x)) 
  return(list(x = x, y = y))
}

load("data/analyses_pADHD_2SB_v9.Rdata")
rm(imputed_data,ipw,missings,my_data,scaling_info, missing_birthvars,
   missing_covariates, missing_MHvars, missing_SDvars, criterion,p,
   round2targets)

pds = c()
table_datas = c()
dds = c()
ddms = c()
varss = c()
tbls = c()
tblbs = c()
plots = vector(mode = "list",length = 3)

for (a in c(1,2,3)) {
  load(paste0("samples/pADHD/a",a,".Rdata"))
  analyses[[a]]$sf = sf
  rm(sf)
  
  vars = strsplit(as.character(analyses[[a]]$formula)[[3]]," + ", fixed = T)[[1]]
  vars = vars[1:analyses[[a]]$sf@par_dims[["b_ipw"]]]
  vars = setdiff(vars,ctrl_vars)
  vars = gsub("5Cs","nCs",vars)
  n = length(vars)
  
  estim = vector(mode = "list",length = 5)
  names(estim) = c("ipw","ar","delta_ar","bias","biasm")
  estim[["ipw"]] = as_matrix(analyses[[a]]$sf,"me_ipw")[,1:n]*-22
  estim[["delta_ar"]] = as_matrix(analyses[[a]]$sf,"delta_me")[,1:n]*-22
  estim[["ar"]] = as_matrix(analyses[[a]]$sf,"me_ar")[,1:n]*-22
  estim[["bias"]] = t(t(estim[["delta_ar"]])/apply(estim[["ipw"]],2,sd))
  estim[["biasm"]] = t(t(estim[["delta_ar"]])/abs(apply(estim[["ipw"]],2,mean)))
  
  
  pd = data.table(var = vars)
  tbl = data.frame(matrix(NA,ncol = length(estim),nrow= length(vars)),row.names = vars)
  names(tbl) = names(estim)
  for (est in names(estim)) {
    pd[[paste0("m_",est)]] = colMeans(estim[[est]][,1:n])
    pd[[paste0("hdi90l_",est)]] = hdi(estim[[est]][,1:n],credMass = .9)[1,]
    pd[[paste0("hdi90u_",est)]] = hdi(estim[[est]][,1:n],credMass = .9)[2,]
    pd[[paste0("hdi50l_",est)]] = hdi(estim[[est]][,1:n],credMass = .5)[1,]
    pd[[paste0("hdi50u_",est)]] = hdi(estim[[est]][,1:n],credMass = .5)[2,]
    
    
    
    for (j in 1:n) {
      tbl[j,est] = ifelse(est %in% c("bias","biasm"),
                          get_stats(estim[[est]][,j], ROPE = c(-.5,.5)),
                          get_stats(estim[[est]][,j]))
      if (est == "bias" | est == "biasm") {
        pd[[paste0("logRR_",est)]][j] = as.numeric(strsplit(get_stats(estim[[est]][,j], ROPE = c(-.5,.5)),";")[[1]][2])
      }
    }
  }
  tblb = 
  data.frame(cbind(tbl[,1:3],
                   sapply(tbl[,4], function(x) strsplit(x,";")[[1]][1]),
                   sapply(tbl[,4], function(x) strsplit(x,";")[[1]][2]),
                   sapply(tbl[,5], function(x) strsplit(x,";")[[1]][1]),
                   sapply(tbl[,5], function(x) strsplit(x,";")[[1]][2])))
  names(tblb)[4:7] = c("bias","","biasm","")
  
  table_data = c()
  for (est in names(estim)) {
    dt = data.table(t(rbind(colMeans(estim[[est]][,1:n]),
                            apply(estim[[est]][,1:n],2,sd),
                            hdi(estim[[est]][,1:n],credMass = .9),
                            colMeans(estim[[est]][,1:n]>0))))
    dt = round(dt,digits = 3)
    setnames(dt,c("mean","sd","HDI90L","HDI90U","%>0"))
    dt[,estimate := est]
    dt[,IV := vars]
    if (is.null(dim(table_data))) {
      table_data = dt
    } else {
      table_data = rbind(table_data,dt)
    }
  }
  table_data = table_data[,c(6,7,1:5),with = F]
  table_data = table_data[estimate != "delta",]
  
  if (a == 1) {
    p_in_rope_bias = apply(estim[["bias"]],2,p_in_rope)
    p_in_rope_biasm = apply(estim[["biasm"]],2,p_in_rope)
    names(p_in_rope_bias) = names(p_in_rope_biasm) = vars
  } else {
    tmp = apply(estim[["bias"]],2,p_in_rope)
    tmpm = apply(estim[["biasm"]],2,p_in_rope)
    names(tmp) = names(tmpm) = vars
    p_in_rope_bias = c(p_in_rope_bias,tmp)
    p_in_rope_biasm = c(p_in_rope_biasm,tmpm)
  }
  
    
  dd = apply(estim$bias[,1:n],2,density)
  ddm = apply(estim$biasm[,1:n],2,density)
  
  pds = rbind(pds, pd)
  table_datas = rbind(table_datas, table_data)
  dds = c(dds,dd)
  ddms = c(ddms,ddm)
  varss = c(varss,vars)
  tbls = rbind(tbls,tbl)
  tblbs = rbind(tblbs,tblb)
}

dd = dds
ddm = ddms
pd = pds
table_data = table_datas
vars = varss
tbl = tbls
tblb = tblbs

tblx = c()
for (k in 1:ncol(tblb)) {
  tmp = tblb[,k]
  if (length(grep("\\(",tmp)) > 0) {
    tmp = gsub("^ ","",tmp)
    tblx = cbind(tblx,
                 t(sapply(tmp,function(x) strsplit(x,split = " ")[[1]])))
  } else {
    tmp = as.character(tmp)
    tmp = gsub("^ ","",tmp)
    tblx = cbind(tblx,tmp)
  }
}

row.names(tblx) = row.names(tblb)

rm(dds,ddms,pds,table_datas,varss,tbls)
n = length(vars)
pd[,y := rev(1:n)]
pd[,y1 := y -.15]
pd[,y2 := y +.15]
rownames(tblx) = gsub("small_fga","SGA",row.names(tblx) )
rownames(tblx) = gsub("SCORE","score",row.names(tblx) )
rownames(tblx) = gsub("TYPUNITS","typunits",row.names(tblx) )
rownames(tblx) = gsub("EFFUNITS","effunits",row.names(tblx) )
names(tbl) = c("IPW","AR","AR-IPW","(AR-IPW)/\\sigma[IPW]","(AR-IPW)/\\mu[IPW]")
colnames(tblx) = c("IPW","HDIIPW","AR","HDIAR","delta","HDIdelta","deltasd",
                   "HDIdeltasd","lRRdeltasd","deltamu","HDIdeltamu","lRRdeltamu")
print(xtable(tblx,caption = "Parameter estimates"),file = "LTX/tables/estimates.tex")

#####################################
############### plot ################
#####################################

library(tikzDevice)
tikz(file = "LTX/figures/estimates.tex", width = 15/2.54, height = 18/2.54)

mat = matrix(c(1,1,1,2,2),ncol = 5)
layout(mat = mat)
par (mar=c(3,8.5,2,.1), mgp=c(2,.7,0), tck=-.01)

############### plot estimates ################

plot("n",
     ylim = range(pd[,grep("y",names(pd)),with = F]),
     xlim = range(pd[,grep("m_ipw|m_ar|u_ipw|u_ar|l_ipw|l_ar",names(pd)),with = F]),
     ylab = "",
     xlab = "average marginal effect",
     yaxt = "n")
abline(v = 0, col = "grey",lty = 2)

k = 0
yos = seq(-.125,.125,length = 2)
cols = brewer.pal(3,"Dark2")[1:2]
pchs = 15:16
for (est in c("ipw","ar")) {
  k = k+1
  points(pd[[paste0("m_",est)]],pd$y + yos[k], pch = pchs[k], col = cols[k],cex = 1.5)
  segments(x0 = pd[[paste0("hdi50l_",est)]],
           x1 = pd[[paste0("hdi50u_",est)]],
           y0 = pd$y + yos[k],
           col = cols[k],lwd = 3)
  segments(x0 = pd[[paste0("hdi90l_",est)]],
           x1 = pd[[paste0("hdi90u_",est)]],
           y0 = pd$y + yos[k],
           col = cols[k],lwd = 1)
}

lbs = gsub("small_fga","SGA",pd$var)
lbs = gsub("SCORE","score",lbs)
lbs = gsub("FREQ","freq",lbs)
lbs = gsub("TYPUNITS","typunits",lbs)
axis(2,at = pd$y,labels = lbs, las = 2)

legend(x = par("usr")[1], y = 17.4,
       col = c(cols,"black","black"),
       pch = c(pchs),
       pt.cex = c(2,2,0,0),
       cex = 1.25,
       legend = c("AME$_{IPW}$","AME$_{AR}$","50\\% HDI", "90\\% HDI"),
       lwd = c(0,0,3,1),
       lty = c(1,1,1,1),
       bty = "n",
       horiz = T,
       xpd = T,x.intersp=0)

abline(h = (n-1-.5)+c(-.03,.03))
abline(h = (n-1-6-.5)+c(-.03,.03))

############### plot bias ################
sl = function(x) {return(sign(x)*sqrt(abs(x)))}
sl = function(x) {return(x)}

par (mar=c(3,.1,2,1))

# select variables for which bias is not extreme
# because ipw estimates is close to 0
vidx = which(apply(abs(pd[,grep("bias",names(pd)),with = F])<100,1,all))

tmp = sort(as.matrix(pd[, grep("bias",names(pd)), with = F]))
xlim = c(tmp[3],tail(tmp,4)[1])

xlim = range(as.matrix(pd[, grep("bias",names(pd)), with = F])[grep("fSMOKE",pd$var,invert = T),],na.rm = T)

plot(0,
     type = "n",
     ylim = range(pd[,grep("y",names(pd)),with = F]),
     xlim = xlim,
     ylab = "",
     xlab = "bias",
     yaxt = "n")
#abline(v = 0, col = "grey",lty = 2)
abline(v = c(-.5,.5), col = "red",lty = 2)
pd$y0a = pd$y + min(yos)
pd$y0b = pd$y + max(yos)
segments(x0 = sl(pd$hdi90l_bias),x1 = sl(pd$hdi90u_bias),y0 = pd$y0b, lwd = 1)
segments(x0 = sl(pd$hdi50l_bias),x1 = sl(pd$hdi50u_bias),y0 = pd$y0b, lwd = 3)
segments(x0 = sl(pd$hdi90l_biasm),x1 = sl(pd$hdi90u_biasm),y0 = pd$y0a, lwd = 1, col = "blue")
segments(x0 = sl(pd$hdi50l_biasm),x1 = sl(pd$hdi50u_biasm),y0 = pd$y0a, lwd = 3, col = "blue")


cols = c("black","black","blue")
pchs = c(17,17,18)
points(sl(pd$m_bias),pd$y0b, col = cols[1], pch = pchs[1],cex = 2)
points(sl(pd$m_biasm),pd$y0a, col = cols[3], pch = pchs[3],cex = 2)


legend(x = par("usr")[1], y = 17.4,
       col = c("white","black","blue"),
       pch = c(pchs),
       pt.cex = c(0,2,2),
       cex = 1.25,
       legend = c(" ",
                  "$\\delta/\\sigma_{IPW}$",
                  "$\\delta/\\mu_{IPW}$"),
       lwd = c(0,0,0),
       lty = c(0,0,0),
       bty = "n",
       horiz = T,
       xpd = T, text.width = 2)
abline(h = (n-1-.5)+c(-.03,.03))
abline(h = (n-1-6-.5)+c(-.03,.03))

dev.off()




#### Kruscke ROPE plots ###

plot_ropehdi = function(my_pd,add = F, ylab = T, xlab = T) {
  col = ifelse(add,"blue","black")
  idx = my_pd$x<3
  if (add == F) {
    plot(my_pd$x[idx],
         my_pd$y[idx],
         'l',
         xlab = ifelse(xlab == T, "ROPE radius",""),
         ylab = ifelse(ylab == T,"prop. posterior in ROPE",""),
         main = names(p_in_rope_bias)[k],
         ylim = c(0,1),
         col = col,
         xaxt = "n",
         yaxt = "n",
         bty = "n")
    axis(1,pos = 0)
    axis(2,pos = 0)
    lines(c(4,4),c(0,1))
    lines(c(.5,.5),c(0,1), col = "red", lty = 2)
  } else {
    lines(my_pd$x,
          my_pd$y,
          col = col)
  }
  
  idx95 = which.min(abs(my_pd$y -.90))
  idx05 = which.min(abs(my_pd$y -.10))
  points(my_pd$x[idx95],0,pch = 15, col = col)
  points(my_pd$x[idx05],0,pch = 25, col = col, bg = col)
  #lines(rep(my_pd$x[idx95],2),c(0,.9),lty = 2, col = col)
  #lines(c(0,my_pd$x[idx95]),rep(.9,2),lty = 2, col = col)
  #lines(rep(my_pd$x[idx05],2),c(1,.1),lty = 3, col = col)
  #lines(c(0,my_pd$x[idx05]),rep(.1,2),lty = 2, col = "red")
  
}

names(p_in_rope_bias) = gsub("small_fga","SGA",names(p_in_rope_bias) )
names(p_in_rope_bias) = gsub("SCORE","score",names(p_in_rope_bias) )
names(p_in_rope_bias) = gsub("TYPUNITS","typunits",names(p_in_rope_bias) )
names(p_in_rope_bias) = gsub("EFFUNITS","effunits",names(p_in_rope_bias) )
names(p_in_rope_bias) = gsub("FREQ","freq",names(p_in_rope_bias) )


tikz(file = "LTX/figures/ROPE_plots.tex", width = 15/2.54, height = 15/2.54)
par(mar=c(2,2,1.5,.01), mgp=c(1,.2,0), tck=.01, mfrow = c(4,4))
for (k in 1:length(p_in_rope_bias)) {
  plot_ropehdi(p_in_rope_bias[[k]],
               xlab = ifelse(k>12,T,F),
               ylab = ifelse(k %in% c(1,5,9,13),T,F))
  plot_ropehdi(p_in_rope_biasm[[k]],add = T)
}
dev.off()

tikz(file = "LTX/figures/logRRs.tex", width = 15/2.54, height = 18/2.54)
par (mar=c(3,8,1,1),xpd = F,mgp=c(1.75,.5,0), tck=.01, mfrow = c(1,1))
plot(0, type = "n",
     ylim = range(pd[,grep("y",names(pd)),with = F]),
     xlim = log(100)*c(-1,1),
     ylab = "",
     xlab = "log(RR) [bias in ROPE]",
     yaxt = "n", bty = "n", xaxt = "n")
rect(-log(3),par("usr")[3],log(3),par("usr")[4],col = adjustcolor("gray",alpha = .3), border = NA)
xticks = sort(c(-log(c(3,10,30,100)),log(c(3,10,30,100))))
axis(1,at = xticks, labels = exp(abs(xticks)))
abline(v = xticks,lty = 3, col = "grey")

#abline(v = 0,col = "red",lwd = 2)
mtext("bias",side = 3, line = 0,adj = .25,col = "red")
mtext("no bias",side = 3, line = 0,adj = .75,col = "green3")
par(mgp=c(1.75,.7,0))
axis(2,at = pd$y, labels = lbs, las = 2,lty = 0)
offset = .15
segments(x0 = 0,
         x1 = pd[,logRR_biasm],
         y0 = pd$y-offset,
         lwd = 12,
         lend = 1,
         col = "blue")
segments(x0 = 0,
         x1 = pd[,logRR_bias],
         y0 = pd$y+offset,
         lwd = 12,
         lend = 1,
         col = "black")
par(xpd = T)
text(mean(log(c(3,10))),3,"moderate",srt = 90)
text(mean(log(c(10,30))),3,"strong",srt = 90)
text(mean(log(c(30,100))),3,"very strong",srt = 90)
text(mean(log(c(100,exp(par("usr")[2])))),3,"extreme",srt = 90)

legend("topright",
       pch = 15,
       col = c("blue","black"),
       legend = c("$\\delta/\\mu_{IPW}$",
                  "$\\delta/\\sigma_{IPW}$"),
       box.col = "white")
dev.off()