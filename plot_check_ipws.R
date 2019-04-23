library(rstanarm)
library(tikzDevice)
library(RColorBrewer)
library(data.table)
library(Hmisc)
library(wCorr)

get_tables = function(ipw_data,weights) {
  return(
  lapply(c("parity","Age","edu"), function(v) {
    cbind(prop.table(wtd.table(ipw_data[,get(v)], weights = ipw_data$N.ssb)[[2]]),
          prop.table(wtd.table(ipw_data[,get(v)], weights = ipw_data$N.moba)[[2]]),
          prop.table(wtd.table(ipw_data[,get(v)], weights = ipw_data$N.moba*weights)[[2]]))
  })
  )
}

wc = function(x,y,wght) {
  weightedCorr(x,
               y,
               method = "Spearman",
               weights = wght)  
}

get_wcorr = function(ipw_data,weights) {
  sapply(list(ipw_data$N.ssb,ipw_data$N.moba, ipw_data$N.moba*weights), 
         function(w){
           c(wc(ipw_data$edu, ipw_data$Age, w),
             wc(ipw_data$edu, ipw_data$parity, w),
             wc(ipw_data$parity, ipw_data$Age, w))
         } )
}


cols = brewer.pal(4,"RdYlBu")

cols = adjustcolor(c("blue","green","orange","red"),alpha = .2)

pps = c()
sipws = c()

t_age = matrix(0, nrow = 6, ncol = 3)
t_parity = matrix(0, nrow = 4, ncol = 3)
t_edu = matrix(0, nrow = 4, ncol = 3)
corrs = matrix(0, nrow = 3, ncol = 3)

for (k in 1:20) {
  load(paste0("data/IPW_i",k,"_pADHD_2SB_v9.Rdata"))
  pps = rbind(pps,colMeans(posterior_predict(mx2)))
  sipw = sipw[,!apply(is.infinite(sipw),2,any)]
  sipws = rbind(sipws,rowMeans(sipw))
  
  tts = get_tables(ipw_data,sipw[,k])
  t_age = t_age + (tts[[2]]-t_age)/k
  t_parity = t_parity + (tts[[1]]-t_parity)/k
  t_edu = t_edu + (tts[[3]]-t_edu)/k
  corrs = corrs + (get_wcorr(ipw_data,sipw[,k])-corrs)/k
}


## plot balance of samples
tikz('LTX/figures/IPWbalance.tex', width=12/2.54, height=15/2.54, pointsize = 14, standAlone = F)
par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)
#layout(matrix(c(rep(1:3,3),rep(4,6)), nrow = 3))
layout(matrix(c(1,2,3,4,4),ncol = 1))
tts = list(Age = t_age, parity = t_parity, edu = t_edu)
ipw_data$parity = ordered(ipw_data$parity)
col = RColorBrewer::brewer.pal(3,"Dark2")
for (k in 1:3) {
  barplot(t(tts[[k]]),
          beside = T,
          names.arg = levels(ipw_data[,get(names(tts)[k])]),
          main = names(tts)[k],
          col = col,
          ylab = "proportion")
  
  if (k == 2)
    legend("topright",
           legend = c("Population",
                      "Unweighted MoBa",
                      "Weighted MoBa"),
           fill = col,
           bty = "n", cex = 1)
  
}
lim = range(corrs)
plot(0, type = "n",
     xlim = lim, ylim = lim, 
     xlab = "Population correlation",
     ylab = "Sample correlation")
abline(0,1, lty = 2)
points(corrs[,1],corrs[,2], pch = 16, col = col[3], cex = 2)
points(corrs[,1],corrs[,3], pch = 17, col = col[2], cex = 2)

legend("topleft",
       pch = c(17,16),
       legend = c("Unweighted MoBa",
                  "Weighted MoBa"),
       col = col[-1],
       bty = "n", cex = 1.5)
dev.off()

## plot oberved vs predicted participation rates

obs = ipw_data$N.moba / ipw_data$N.ssb
pred = colMeans(pps) / ipw_data$N.ssb

sipw = colMeans(sipws)

for (L in levels(ipw_data$edu))
  ipw_data[edu == L, col := cols[which(L == levels(ipw_data$edu))]]

tikz('LTX/figures/IPWpres.tex', width=20/2.54, height=10/2.54, pointsize = 14, standAlone = F)
par(mfrow = c(2,1))
par(mar=c(2,2,2,1), mgp=c(1,.7,0))

sobs = sqrt(obs)
spred = sqrt(pred)
lims = c(0,max(c(sobs,spred)))
plot(0, type = "n", pch = 20, col = "red", cex = .5, ylim = lims, xlim = lims,
     xlab = "predicted participation (\\%)",ylab = "observed participation (\\%)",
     xaxt  = "n",
     yaxt = "n",
     bty = "n")
lines(c(0,lims[2]),c(0,lims[2]), "l")
points(sobs,spred, cex = sqrt(ipw_data$N.ssb)/20, pch = 19,
       col = ipw_data$col)
#lbs =  c(0, .01 , .025, .05, .1, .2, .3, .4)
lbs =  c(0, .01 , .05, .1, .2, .4)
axis(1,at = sqrt(lbs), labels = lbs*100, pos = 0)
axis(2,at = sqrt(lbs), labels = lbs*100, pos = 0)


Ns = c(250,2500,10000, 40000)
ys = seq(sqrt(.001),sqrt(.075),length = 5)
points(rep(sqrt(.2),4),rep(ys[2],4), cex = sqrt(Ns)/20)
text(sqrt(.2),ys[3],
     paste0("$N_{SSB}=$",paste(Ns, collapse = ", ")),
     pos = 3)

legend(x = sqrt(.0005), y = sqrt(.4),
       pch = 16, col = cols,
       border = cols,
       legend = levels(ipw_data$edu),bty = "n", pt.cex = 2)

breaks = seq(log(min(sipws*.99)),log(max(sipws*1.01)),length = 31)
hs = vector(length = 20, mode = "list")
for (k in 1:20) hs[[k]] = hist(log(sipws[k,]),breaks = breaks, plot = F)
ylim = c(0,max(sapply(hs,function(x) max(x$counts))))


h = hist(log(sipws[1,]),
         xaxt = "n",
         main = "",
         xlab = "sIPW",
         breaks = breaks,
         ylim = ylim, yaxt = "n",
         col = adjustcolor("black", alpha = .1), border = adjustcolor("black", alpha = .05))
for (k in 2:20) hist(log(sipws[k,]),breaks = breaks, add = T,
                     col = adjustcolor("black", alpha = .05),
                     border = adjustcolor("black", alpha = .05))

lbs =  exp(h$mids)[seq(1,length(h$mids),by = 4)]
for (k in 1:length(lbs)) lbs[k] = ifelse(lbs[k] < 5, round(lbs[k],digits = 1), round(lbs[k]))
axis(1,at = log(lbs),labels  = lbs, pos = 0)
axis(2,at = seq(2,max(ylim),by = 2), pos = log(lbs)[1])

dev.off()
