library(data.table)
library(gplots)
library(RColorBrewer)
library(tikzDevice)
library(xtable)

rename_vars = function(x, ft = F) {
  x = gsub("YearsEdu", "Years\neducation", x)
  x = gsub("AgeFirstBirth", "Age at\nfirst birth", x)
  x = gsub("NumChildrBorn", "Number\nchildr. born", x)
  x = gsub("BirthWeight", "Birthweight", x)
  x = gsub("AlcUse", "Alcohol\nuse", x)
  x = gsub("CigPerDay", "Cigarettes\nper day", x)
  x = gsub("EverSmoked", "Ever\nsmoked", x)
  x = gsub("Cannabis", "Cannabis\nuse", x)
  x = gsub("SubstAbuse", "Substance\nabuse", x)
  x = gsub("Anxiety", "Anxiety\ndisorder", x)
  x = gsub("DeprSymp", "Depression\nsymptoms", x)
  x = gsub("ADHD", "ADHD", x)
  
  my_levels = c("Years\neducation", "Age at\nfirst birth", 
                "Number\nchildr. born", "Birthweight", 
                "Alcohol\nuse", "Cigarettes\nper day" ,
                "Ever\nsmoked", "Cannabis\nuse", "Substance\nabuse",
                "Anxiety\ndisorder", "Depression\nsymptoms", "ADHD")
  if(ft == T) {
    x = gsub("\n"," ", x)
    my_levels = gsub("\n"," ", my_levels)
    my_levels = my_levels[c(12,1:11)]
  }
  
  x = factor(x, levels = my_levels)
  
  return(x)
}

vn = c("p1","p2","rg","se","z","p","h2_obs","h2_obs_se","h2_in","h2_int_se","gcov_int","gcov_int_se")
fls = dir("GWASresults/my_data/rg")
for ( f in fls) {
  tmp = fread(paste0("GWASresults/my_data/rg/",f),skip = 59)
  setnames(tmp,names(tmp),vn)
  if (f == dir("GWASresults/my_data/rg")[1]) {
    rg = tmp
  } else {
    rg = rbind(rg,tmp) 
  }
    
}

rg[,p1 := gsub(".sumstats.gz","",p1)]
rg[,p2 := gsub(".sumstats.gz","",p2)]

pse = function(mstr){
  return(c(as.numeric(strsplit(strsplit(mstr,":")[[1]][2]," ")[[1]][2]),
         as.numeric(gsub("\\(|\\)","",strsplit(strsplit(mstr,":")[[1]][2]," ")[[1]][3]))))
}

geth2stats = function(txt) {
  return(data.table(
    var = sub(".h2.log","",f),
    h2 = pse(txt[26])[1],
    h2_se = pse(txt[26])[2],
    lambda_GC = as.numeric(strsplit(txt[27],":")[[1]][2]),
    mean_chi_sq = as.numeric(strsplit(txt[28],":")[[1]][2]),
    intercept = pse(txt[29])[1],
    intercept_se = pse(txt[29])[1],
    ratio = pse(txt[30])[1],
    ratio_se = pse(txt[30])[1]
  ))
}

fls = dir("GWASresults/my_data/h2")
for ( f in fls) {
  tmp = readLines(paste0("GWASresults/my_data/h2/",f))
  if (f == dir("GWASresults/my_data/h2")[1]) {
    h2 = geth2stats(tmp)
  } else {
    h2 = rbind(h2,geth2stats(tmp)) 
  }
}

rg = rg[p1 != "BWfirst" & p2 != "BWfirst"]
rg[,p1 := sub("subabu","SubstAbuse", p1)]
rg[,p2 := sub("subabu","SubstAbuse", p2)]
h2 = h2[var != "BWfirst"]
h2[,var := sub("subabu","SubstAbuse", var)]


rm(tmp,vn,f)
m_rg = matrix(NA,nrow = nrow(h2),ncol = nrow(h2))
for (v1 in 1:nrow(h2)) {
  for (v2 in 1:nrow(h2)) {
    v1n = h2$var[v1]
    v2n = h2$var[v2]
     if (v1 == v2) {
       m_rg[v1,v2] = h2[var == v1n,h2]
     } else {
       idx = apply(cbind(rg$p1 %in% c(v1n,v2n),rg$p2 %in% c(v1n,v2n)),1,all)
       if(sum(idx) == 1) {
         m_rg[v1 ,v2] = rg[idx,rg]
       } else {
         print(c(v1n,v2n))
       }
     }
  }
}

rownames(m_rg) = h2$var
colnames(m_rg) = h2$var
hm = heatmap.2(m_rg,trace = "n", col = brewer.pal(11,"RdYlGn"),breaks = 12,margins = c(10,10))

mm = ceiling(max(abs(rg$rg)*10))/10*1.05
h = hist(rg$rg,breaks = seq(-mm,mm,length = 22),plot = F)
breaks = h$breaks
n = length(h$mids)


my_rg = rg[,c("p1","p2","rg","z","p"),with  = F]
setnames(my_rg,c("p1","p2"),c("V1","V2"))
my_rg_b = my_rg
my_rg_b$V1 = my_rg$V2
my_rg_b$V2 = my_rg$V1
my_rg = rbind(my_rg,my_rg_b)
rm(my_rg_b)

h2[,z := h2/h2_se]
h2b = h2
setnames(h2b,c("var"),c("V1"))
h2[,rg := h2]
h2b[,V2 := V1]
pdata  = rbind(my_rg[,c("V1","V2","rg","z"),with = F],
               h2b[,c("V1","V2","rg","z"),with = F])

pdata$az = abs(pdata$z)
pdata$s = sqrt(pdata$az)/10


#########################################
############ plot figure ################
#########################################

res = 101
# blue red
clrs = colorRampPalette( c("blue", "white", "red"), space="rgb")(res)
max_rg = round(max(abs(pdata$rg)),digits = 1)
w = seq(from = -max_rg, to = max_rg, length = res)
dt = data.table(w, val = w)
setattr(dt, "sorted", "w") 
pdata[,rgcol := dt[J(pdata$rg),roll = "nearest"][,val]]
dt$clr = clrs
dt[,rgcol := val]
pdata = merge(pdata,dt[,c("rgcol","clr")],by = "rgcol", all.x = T, all.y = F)

# white black
bw = colorRampPalette( c("white", "black"), space="rgb")(res)
max_h2 = .5
w = seq(from = 0, to = 1, length = res)
dt = data.table(w, val = w)
setattr(dt, "sorted", "w") 
pdata[,bwcol := dt[J(pdata$rg),roll = "nearest"][,val]]
dt$bw = bw
dt[,bwcol := val]
pdata = merge(pdata,dt[,c("bwcol","bw")],by = "bwcol", all.x = T, all.y = F)
pdata[V1 == V2,clr := bw]


pdata$cex = (1-pdata$p)*4.5

pdata$V1 = rename_vars(pdata$V1)
pdata$V2 = rename_vars(pdata$V2)

pdata$x = as.numeric(pdata$V1)
pdata$y = as.numeric(factor(pdata$V2, levels = rev(levels(pdata$V2))))


nvars = length(unique(pdata$V1))

tikz("LTx/figures/rg.tex", height = 13/2.54, width = 11/2.54, standAlone = F)
#postscript("LTx/figures/eps/Fig2_rg.eps", height = 15/2.54, width = 12.9/2.54, paper = "a4", horizontal = F, pagecentre = T)

layout(matrix(c(rep(1,3),2),ncol = 1))

scale_z = max(sqrt(pdata$az))/.5
pdata$s = sqrt(pdata$az)/scale_z


par(mar = c(0,7,1,0.5))
plot(0,type = "n", xlim = c(.5,nvars+.5), ylim = c(.5,nvars+.5),
     xaxt = "n", yaxt = "n",
     ylab = "", xlab = "",bty = "n")
msz = max(sqrt(pdata$az))
for (r in 1:nrow(pdata)) {
  s = pdata$s[r]
  rect(pdata$x[r]-sqrt(4)/scale_z,
       pdata$y[r]-sqrt(4)/scale_z,
       pdata$x[r]+sqrt(4)/scale_z,
       pdata$y[r]+sqrt(4)/scale_z,
       col = "white",
       border = "gray")
  rect(pdata$x[r]-s,
       pdata$y[r]-s,
       pdata$x[r]+s,
       pdata$y[r]+s,
       col = pdata$clr[r],
       border = NA)
  if (pdata[r,V1] == pdata[r,V2]) {
    text(pdata$x[r],pdata$y[r],
         substr(sprintf("%0.2f", pdata$rg[r]),2,4),cex = .75)
  }
}
segments(x0 = rep(0,nvars)+.5,
         y0 = (1:(nvars+1))-.5,
         x1 = rep(nvars,nvars)+.5,
         col = "lightgray")
segments(y0 = rep(0,nvars)+.5,
         x0 = (1:(nvars+1))-.5,
         y1 = rep(nvars,nvars)+.5,
         col = "lightgray")

lbs = sub("_","-",levels(pdata$V1))
mtext(paste0(1:length(lbs),c(rep("\n ",11),"")),2,at = rev(1:nvars), 
      las = 2, cex = .75, adj = 0, line = 6.75)
mtext(lbs,2,at = rev(1:nvars), 
      las = 2, cex = .75, adj = 0, line = 5.25)
mtext(1:nvars,3,at = 1:nvars, cex = .75, line = -1)
#axis(2,at = 1:nvars, labels = rev(lbs),las = 2, tick = F,line = -1.5)
#lbs = sub("AlcDependence","AlcDep.",lbs)
#text(1:nvars, par("usr")[2] , labels = lbs, srt = 45, pos = 4, xpd = TRUE,offset = .1)
#text(1:nvars, par("usr")[2] , labels = lbs, srt = 360-45, pos = 2, xpd = TRUE,offset = .1)


par(mar = c(0,7,0,1))
plot(0,type = "n", xaxt = "n", yaxt = "n",bty = "n",ylab = "", xlab = "",
     xlim = c(0,10),ylim = c(-.5,2.5))

x = seq(2,10,length = length(clrs)+1)
for (i in 1:(length(x)-1)) {
  rect(x[i],1.75,x[i+1],2.25, col=clrs[i], border=NA)
}
text(0,2,"$rG$",pos = 4,cex = 1.25)
lbs = axisTicks(c(-max_rg,max_rg),log = F, nint = 6)
b = solve(matrix(c(1,1,max_rg*c(-1,1)),ncol = 2),range(x))
ticks = cbind(rep(1,length(lbs)),lbs) %*% b
segments(x0 = ticks,x1 = ticks,y0 = 2.25,y1 = 2.3)
text(ticks,
     2.25,
     labels = lbs,
     pos = 3,
     offset = 0.5)

x = seq(2,10,length = length(bw)+1)
for (i in 1:(length(x)-1)) {
  rect(x[i],1.25,x[i+1],1.75, col=bw[i], border=NA)
}
text(0,1.5,"$h^2_{SNP}$",pos = 4,cex = 1.25)
lbs = axisTicks(c(0,max_h2),log = F, nint = 6)
b = solve(matrix(c(1,1,max_h2*c(0,1)),ncol = 2),range(x))
ticks = cbind(rep(1,length(lbs)),lbs) %*% b
segments(x0 = ticks,x1 = ticks,y0 = 1.25,y1 = 1.2)
text(ticks,
     1.25,
     labels = lbs,
     pos = 1,
     offset = 0.5)


box_size = .4
zs = c(2,4,8,16,32)
r = seq(from = 2+box_size, to = 10-box_size, length = length(zs))
s = (zs^.5)/max(zs)^.5/(1/box_size)
scale4dev = diff(grconvertY(c(0,1), from = "user", to = "inches"))/diff(grconvertX(c(0,1), from = "user", to = "inches"))
yy = .45
for (k in 1:5) {
  rect(r[k]-box_size,
       yy-box_size/scale4dev,
       r[k]+box_size,
       yy+box_size/scale4dev,
       col = "white",
       border = "gray")
  rect(r[k]-s[k],
       yy-s[k]/scale4dev,
       r[k]+s[k],
       yy+s[k]/scale4dev,
       col = "gray",
       border = "gray")
  text(r[k],yy-box_size/scale4dev,zs[k],pos = 1, offset = .4)
}
text(0,yy,"$z$",pos = 4,cex = 1.25)

dev.off()

##############################################
################# make table #################
##############################################

setnames(h2,
         c("V1","h2","h2_se","z"),
         c("p2","h2p2","se(h2p2)","z(h2p2)"))
tab = rg[,c("p1","p2","rg","se","z"),with = F]
for (k in which(tab$p2 == "ADHD")) {
  tab[["p2"]][k] = tab[["p1"]][k]
  tab[["p1"]][k] = "ADHD"
}


for (k in 1:nrow(tab)) {
  p1 = tab$p1[k]
  p2 = tab$p2[k]
  if (as.numeric(rename_vars(p2,ft = T)) < as.numeric(rename_vars(p1,ft = T))) {
    tab$p1[k] = p2
    tab$p2[k] = p1
  }
}

tab = merge(tab,
            h2[,c("p2","h2p2","se(h2p2)","z(h2p2)"),with = F],
            by = "p2",
            all.x = T)
setnames(h2,
         c("p2","h2p2","se(h2p2)","z(h2p2)"),
         c("p1","h2p1","se(h2p1)","z(h2p1)"))
tab = merge(tab,
            h2[,c("p1","h2p1","se(h2p1)","z(h2p1)"),with = F],
            by = "p1",
            all.x = T)
tab[,p1 := rename_vars(p1, ft = T)]
tab[,p2 := rename_vars(p2, ft = T)]


setnames(tab,
         c("rg","se","z"),
         c("rG","se(rG)","z(rG)"))

setkeyv(tab,c("p1","p2"))

xtab = xtable(tab, 
              digits=c(0,0,0,2,2,1,2,2,1,2,2,1),
              caption = "SNP based genetic correlations and heritability estimates from publicly available GWAS summary statistics.")
print(xtab, hline.after=c(-1, 0),
      tabular.environment = "longtable",
      file = "LTX/tables/rg_h+.tex",
      include.rownames=FALSE)
