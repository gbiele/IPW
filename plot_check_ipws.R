library(rstanarm)
library(tikzDevice)
library(RColorBrewer)
library(data.table)

cols = brewer.pal(4,"RdYlBu")

cols = adjustcolor(c("blue","green","orange","red"),alpha = .2)

pps = c()
sipws = c()
for (k in 1:20) {
  load(paste0("data/IPW_i",k,"_pADHD_2SB_v9.Rdata"))
  pps = rbind(pps,colMeans(posterior_predict(mx2)))
  sipw = sipw[,!apply(is.infinite(sipw),2,any)]
  sipws = rbind(sipws,rowMeans(sipw))
}



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
