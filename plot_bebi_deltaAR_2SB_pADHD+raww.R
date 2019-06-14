setwd("F:/Forskningsprosjekter/PDB 299 - ADHD-studien Prescho_/Forskningsfiler/GUBI/GuidoData/representativeness/IPW")
library(HDInterval)
library(rstan)
library(data.table)
library(RColorBrewer)
library(bayesplot)
library(xtable)
library(gtools)

#\\resizebox{\\textwidth}{!}{
ltx_table_start = "
	
		\\begin{tabular}{rLlLlLlLlLLlL}
\\hline
& \\multicolumn{2}{c}{$IPPW$} & 
\\multicolumn{2}{c}{$AR$} &
\\multicolumn{2}{c}{$\\delta=AR-IPPW$} &
\\multicolumn{3}{c}{$\\delta/\\sigma_{IPPW}$ } &
\\multicolumn{3}{c}{$\\delta/\\mu_{IPPW}$ }  \\\\
&\\multicolumn{1}{c}{\\footnotesize m} & \\multicolumn{1}{c}{\\footnotesize HDI} &
\\multicolumn{1}{c}{\\footnotesize m} & \\multicolumn{1}{c}{\\footnotesize HDI} &
\\multicolumn{1}{c}{\\footnotesize m} & \\multicolumn{1}{c}{\\footnotesize HDI} &  
\\multicolumn{1}{c}{\\footnotesize m} & \\multicolumn{2}{c}{\\footnotesize HDI} {\\footnotesize log(RR$_b$)} &
\\multicolumn{1}{c}{\\footnotesize m} & \\multicolumn{2}{c}{\\footnotesize HDI} {\\footnotesize log(RR$_b$)} \\\\
\\hline
"
ltx_table_end = "
\\hline
\\end{tabular}
"
#}

orig_varnames = c("small_fga", "preterm",
                  "mDRUGS", "mDRUGS.SCORE", "mLTHD", "mSCL5", "fDRUGS", "fDRUGS.SCORE",
                  "mSMOKE", "mSMOKE.nCs", "mALC.FREQ", "mALC.EFFUNITS",
                  "fSMOKE", "fSMOKE.nCs", "fALC.FREQ", "fALC.TYPUNITS")

clean_varnames = c("small f.\ngest. age", "preterm",
                  "m. drug use", "m. drug use\nscore", "m. LTH", "m. SCL5",
                  "p. drug use", "p. drug use\nscore",
                  "m. smoking", "m. num.\ncigarettes", "m. alc. \nfreq.", "m. eff.\nunits alc.",
                  "p. smoking", "p. num.\ncigarettes", "p. alc. \nfreq.", "p. typ.\nunits alc.")


ft = "eps"

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
dds_ar = ddms_ar = c()
dds_raw = ddms_raw = c()
varss = c()
tbls = c()
tblbs = c()
plots = vector(mode = "list",length = 3)

for (a in c(1,2,3)) {
  load(paste0("samples/pADHD/a",a,"+raww.Rdata"))
  analyses[[a]]$sf = sf
  rm(sf)
  
  vars = strsplit(as.character(analyses[[a]]$formula)[[3]]," + ", fixed = T)[[1]]
  vars = vars[1:analyses[[a]]$sf@par_dims[["b_ipw"]]]
  vars = setdiff(vars,ctrl_vars)
  vars = gsub("5Cs","nCs",vars)
  n = length(vars)
  
  estim = vector(mode = "list",length = 9)
  names(estim) = c("ipw","ar","delta_ar","bias_ar","biasm_ar","raw","delta_raw","bias_raw","biasm_raw")
  estim[["ipw"]] = as_matrix(analyses[[a]]$sf,"me_ipw")[,1:n]*-22
  estim[["delta_ar"]] = as_matrix(analyses[[a]]$sf,"delta_ipw_ar")[,1:n]*-22
  estim[["delta_raw"]] = as_matrix(analyses[[a]]$sf,"delta_ipw_raw")[,1:n]*-22
  estim[["ar"]] = as_matrix(analyses[[a]]$sf,"me_ar")[,1:n]*-22
  estim[["raw"]] = as_matrix(analyses[[a]]$sf,"me_raw")[,1:n]*-22
  estim[["bias_ar"]] = t(t(estim[["delta_ar"]])/apply(estim[["ipw"]],2,sd))
  estim[["biasm_ar"]] = t(t(estim[["delta_ar"]])/abs(apply(estim[["ipw"]],2,mean)))
  estim[["bias_raw"]] = t(t(estim[["delta_raw"]])/apply(estim[["ipw"]],2,sd))
  estim[["biasm_raw"]] = t(t(estim[["delta_raw"]])/abs(apply(estim[["ipw"]],2,mean)))
  
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
      tbl[j,est] = ifelse(grepl("bias",est),
                          get_stats(estim[[est]][,j], ROPE = c(-.5,.5)),
                          get_stats(estim[[est]][,j]))
      if (grepl("bias",est)) {
        pd[[paste0("logRR_",est)]][j] = as.numeric(strsplit(get_stats(estim[[est]][,j], ROPE = c(-.5,.5)),";")[[1]][2])
      }
    }
  }
  tblb = 
  data.frame(cbind(tbl[,1:3],
                   sapply(tbl[,4], function(x) strsplit(x,";")[[1]][1]),
                   sapply(tbl[,4], function(x) strsplit(x,";")[[1]][2]),
                   sapply(tbl[,5], function(x) strsplit(x,";")[[1]][1]),
                   sapply(tbl[,5], function(x) strsplit(x,";")[[1]][2]),
                   tbl[,6:7],
                   sapply(tbl[,8], function(x) strsplit(x,";")[[1]][1]),
                   sapply(tbl[,8], function(x) strsplit(x,";")[[1]][2]),
                   sapply(tbl[,9], function(x) strsplit(x,";")[[1]][1]),
                   sapply(tbl[,9], function(x) strsplit(x,";")[[1]][2])))
  names(tblb)[4:7] = c("bias_ar","","biasm_ar","")
  names(tblb)[10:13] = c("bias_raw","","biasm_raw","")
  
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
  table_data = table_data[!grepl("delta",estimate),]
  
  if (a == 1) {
    p_in_rope_bias_ar = apply(estim[["bias_ar"]],2,p_in_rope)
    p_in_rope_biasm_ar = apply(estim[["biasm_ar"]],2,p_in_rope)
    p_in_rope_bias_raw = apply(estim[["bias_raw"]],2,p_in_rope)
    p_in_rope_biasm_raw = apply(estim[["biasm_raw"]],2,p_in_rope)
    names(p_in_rope_bias_ar) = names(p_in_rope_biasm_ar) = vars
    names(p_in_rope_bias_raw) = names(p_in_rope_biasm_raw) = vars
  } else {
    tmp = apply(estim[["bias_ar"]],2,p_in_rope)
    tmpm = apply(estim[["biasm_ar"]],2,p_in_rope)
    names(tmp) = names(tmpm) = vars
    p_in_rope_bias_ar = c(p_in_rope_bias_ar,tmp)
    p_in_rope_biasm_ar = c(p_in_rope_biasm_ar,tmpm)
    
    tmp = apply(estim[["bias_raw"]],2,p_in_rope)
    tmpm = apply(estim[["biasm_raw"]],2,p_in_rope)
    names(tmp) = names(tmpm) = vars
    p_in_rope_bias_raw = c(p_in_rope_bias_raw,tmp)
    p_in_rope_biasm_raw = c(p_in_rope_biasm_raw,tmpm)
  }
  
  dds_ar = c(dds_ar, apply(estim$bias_ar[,1:n],2,density))
  ddms_ar = c(ddms_ar, apply(estim$biasm_ar[,1:n],2,density))
  dds_raw = c(dds_ar, apply(estim$bias_raw[,1:n],2,density))
  ddms_raw = c(ddms_ar, apply(estim$biasm_raw[,1:n],2,density))
  
  pds = rbind(pds, pd)
  table_datas = rbind(table_datas, table_data)
  varss = c(varss,vars)
  tbls = rbind(tbls,tbl)
  tblbs = rbind(tblbs,tblb)
}

dd_ar = dds_ar
ddm_ar = ddms_ar
dd_raw = dds_raw
ddm_raw = ddms_raw
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

rm(dds_ar,ddms_ar,dds_raw,ddms_raw,pds,table_datas,varss,tbls)
n = length(vars)
pd[,y := rev(1:n)]
pd[,y1 := y -.15]
pd[,y2 := y +.15]
rownames(tblx) = gsub("small_fga","SGA",row.names(tblx) )
rownames(tblx) = gsub("SCORE","score",row.names(tblx) )
rownames(tblx) = gsub("TYPUNITS","typunits",row.names(tblx) )
rownames(tblx) = gsub("EFFUNITS","effunits",row.names(tblx) )
names(tbl) = c("IPPW","AR","AR-IPPW","(AR-IPPW)/\\sigma[IPPW]","(AR-IPPW)/\\mu[IPPW]")
statsnames = c("mean","HDI","delta","HDIdelta","deltasd","HDIdeltasd","lRRdeltasd","deltamu","HDIdeltamu","lRRdeltamu")
colnames(tblx) = c("IPPW","HDIIPPW",
                   paste0("AR-",statsnames),
                   paste0("UR-",statsnames))

#print(htmlTable(tblx), file = "LTX/tables/estimates+raww.html")
library(ReporteRs)
for (r in c("AR","UR")) {
  if (r == "AR") {
    tmptbl = tblx[,1:12]
  } else {
    tmptbl = tblx[,c(1:2,13:22)]
  }
  fn = paste0("LTX/tables/estimates-",r,".tex")
  table_data = apply(cbind(gsub("\n"," ", clean_varnames),tmptbl),
                     1,
                     function(x)
                       paste(paste0(x,
                                    collapse = " & "),
                             "\\\\"))
  cat(ltx_table_start,file = fn, append = F)
  cat(table_data, sep = "\n",
      file = fn, append = T)
  cat(ltx_table_end,file = fn, append = T)
  
  #print(xtable(tmptbl,caption = "Parameter estimates"),file = paste0("LTX/tables/estimate-",r,".tex"))
  row.names(tmptbl) = gsub("\n"," ", clean_varnames)
  colnames(tmptbl) = c("m","HDI","m","HDI","m","HDI","m","HDI","ln(RR)","m","HDI","ln(RR)")
  FT = vanilla.table(tmptbl, 
                 add.rownames = T)
  FT = addHeaderRow(FT,
                    value = c("","IPPW",r,"delta","deltasd","deltamu"),
                    colspan = c(1,2,2,2,3,3),
                    first = T, 
                    cell.properties = cellProperties(border.left = borderProperties(style = "none"),
                                                     border.right = borderProperties(style = "none"),
                                                     border.top = borderProperties(width = 2),
                                                     border.bottom = borderProperties(width = .5)))
  #FT = addHeaderRow(FT,
  #                  c("","m","HDI","m","HDI","m","HDI","m","HDI","ln(RR)","m","HDI","ln(RR)"))
  doc <- docx()
  doc <- addSection(doc, landscape = TRUE)
  doc <- addFlexTable(doc,
                      FT)
  doc <- addSection(doc)
  
  writeDoc(doc, file = paste0("LTX/tables/estimates-",r,".docx"))
}





save(tblx,pd,file = "plot_bebi_data.Rdata")

#####################################
############### plot ################
#####################################
#pdf(file = "LTX/figures/estimates.pdf",width = 17/2.54, height = 20/2.54, pointsize = 14)


library(tikzDevice)
ft ="tex"
if(ft == "tex") {
  tikz(file = "LTX/figures/estimates+raw.tex", width = 15/2.54, height = 16/2.54, standAlone = F)
} else {
  postscript("LTX/figures/eps/Fig3-estimates.eps", width = 17.4/2.54, height = 18/2.54, paper = "a4", horizontal = F, pagecentre = T)
}


mat = matrix(c(1,2,2,3,3,4,4),ncol = 7)
layout(mat = mat)

############### add labels ################

par(mar=c(3,0,1.5,0), mgp=c(2,.7,0), tck=-.01)
plot(0,type = "n",
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n",
     xlim = c(0,1),
     ylim  = range(pd[,grep("y",names(pd)),with = F]),
     bty = "n")

lvarnames = gsub("m. |p. ","",clean_varnames)
lvarnames = gsub("num.","number",lvarnames)
lvarnames = gsub("small f.","small for",lvarnames)
lvarnames = gsub("typ.\nunits alc.","typical\nunits alcoh.",lvarnames)
lvarnames = gsub("eff.\nunits alc.","effective\nunits alcoh.",lvarnames)
lvarnames = gsub("alc. \nfreq.","alcohol \nfrequency",lvarnames)

text(x = par("usr")[2],
     y = pd$y,
     labels = lvarnames,
     pos = 2
     )

segments(x0 = .25,
         y0 = c(16.5, 14.4, 10.4, 8.4, 4.4),
         y1 = c(14.6, 10.6, 8.6,  4.6, 0.6),
         col = "grey50")
text(x = 0.1, 
     y = c(15.25, 12.25, 9.25, 6.25, 2.25),
     c("child", "maternal","paternal","maternal","paternal"),
     srt = 90,pos = 3,
     col = "grey50")

############### plot estimates ################
par(mar=c(3,0,1.5,.1))

plot(0,
     type = "n",
     ylim = range(pd[,grep("y",names(pd)),with = F]),
     xlim = range(pd[,grep("delta|var|bias|^y",names(pd),invert = T),with = F]),
     ylab = "",
     xlab = "average marginal effect",
     yaxt = "n",
     bty = "n")

lines(c(0,0),c(0,max(pd$y)+.5),col = "grey",lty = 2)
#abline(v = 0, col = "grey",lty = 2)

k = 0
yos = seq(-.15,.15,length = 3)
cols = brewer.pal(3,"Dark2")[1:3]
bgs = c(cols[1],"white","white")
pchs = c(21:22,24)
for (est in c("ipw","ar","raw")) {
  k = k+1
  segments(x0 = pd[[paste0("hdi50l_",est)]],
           x1 = pd[[paste0("hdi50u_",est)]],
           y0 = pd$y + yos[k],
           col = cols[k],lwd = 3)
  segments(x0 = pd[[paste0("hdi90l_",est)]],
           x1 = pd[[paste0("hdi90u_",est)]],
           y0 = pd$y + yos[k],
           col = cols[k],lwd = 1)
  points(pd[[paste0("m_",est)]],pd$y + yos[k], pch = pchs[k], col = cols[k],cex = 1.25,bg = bgs[k])
}

lbs = gsub("small_fga","SGA",pd$var)
lbs = gsub("SCORE","score",lbs)
lbs = gsub("FREQ","freq",lbs)
lbs = gsub("TYPUNITS","typunits",lbs)

legend(mean(par("usr")[1:2]),
       par("usr")[4]+.5,
       col = cols, pt.bg = bgs,pch = c(pchs),pt.cex = 2,
       c("IPPW","AR","UR"),
       title = "AME",
       bty = "n",
       xpd = T,
       ncol = 3,
       xjust = .5,
       x.intersp = 1,
       text.width = c(1,.9,.75)/3.75)


# legend(mean(par("usr")[1:2]),
#        par("usr")[4]+.7,
#        col = cols, pt.bg = bgs,pch = c(pchs),pt.cex = 2,
#        c(ifelse(ft == "tex","$AME_{IPPW}$",expression("AME"["IPPW"])),
#          ifelse(ft == "tex","$AME_{AR}$",expression("AME"["AR"])),
#          ifelse(ft == "tex","$AME_{UR}$",expression("AME"["UR"]))),
#        bty = "n",
#        xpd = T,
#        ncol = 3,
#        xjust = .5)


legend("topright",
       lwd = c(3,1),
       c(ifelse(ft == "tex","50\\% HDI","50% HDI"),
         ifelse(ft == "tex","90\\% HDI","90% HDI")),
       bty = "n",cex = .8,
       inset = .02)

abline(h = (n-1-.5)+c(-.03,.03), col = "grey")
abline(h = (n-1-6-.5)+c(-.03,.03), col = "grey")

############### plot bias ################
sl = function(x) {return(sign(x)*sqrt(abs(x)))}
sl = function(x) {return(x)}

par (mar=c(3,.5,1.5,.1))

for (r in c("ar","raw")) {
  rr = ifelse(r == "ar","AR","UR")
  bias_cols = grep("log",grep("bias",names(pd), value = T),invert = T, value = T)
  xlim = quantile(as.matrix(pd[, bias_cols, with = F])[grep("fSMOKE",pd$var,invert = T),],c(.001,.999))
  xlim = c(-7.75,5)
  plot(0,
       type = "n",
       ylim = range(pd[,grep("y",names(pd)),with = F]),
       xlim = sl(xlim),
       xlab = paste0("standardized difference ",rr),ylab = "",
       yaxt = "n",
       xaxt = "n",
       bty = "n")
  #abline(v = 0, col = "grey",lty = 2)
  axis(1,at = seq(-7.5,7.5,by = 2.5))
  segments(x0 = c(-.5,.5),y0 = c(0,0), y1 = rep(max(pd$y),2)+.5,col = "red",lty = 2)
  
  cols = c("black","black","blue")
  pchs = c(23,23,24)
  bgs = c(cols[1:2],"white")
  
  pd$y0a = pd$y + yos[1]/3*2
  pd$y0b = pd$y - yos[1]/3*2
  segments(x0 = sl(pd[[paste0("hdi90l_bias_",r)]]),x1 = sl(pd[[paste0("hdi90u_bias_",r)]]),y0 = pd$y0b, lwd = 1, xpd = T)
  segments(x0 = sl(pd[[paste0("hdi50l_bias_",r)]]),x1 = sl(pd[[paste0("hdi50u_bias_",r)]]),y0 = pd$y0b, lwd = 3, xpd = T)
  segments(x0 = sl(pd[[paste0("hdi90l_biasm_",r)]]),x1 = sl(pd[[paste0("hdi90u_biasm_",r)]]),y0 = pd$y0a, lwd = 1, col = "blue")
  segments(x0 = sl(pd[[paste0("hdi50l_biasm_",r)]]),x1 = sl(pd[[paste0("hdi50u_biasm_",r)]]),y0 = pd$y0a, lwd = 3, col = "blue")
  
  tmp = which(pd[[paste0("hdi90l_biasm_",r)]] < par("usr")[1])
  if (length(tmp) > 0)
    arrows(x0 = par("usr")[1]*.995,x1 = par("usr")[1],y0 = pd$y0a[tmp], col = "blue", length = .1)
  tmp = which(pd[[paste0("hdi90u_biasm_",r)]] > par("usr")[2])
  if (length(tmp) > 0)
    arrows(x0 = par("usr")[2]*.995,x1 = par("usr")[2],y0 = pd$y0a[tmp], col = "blue", length = .1)
  
  tmp = which(pd[[paste0("hdi90l_bias_",r)]] < par("usr")[1])
  if (length(tmp) > 0)
    arrows(x0 = par("usr")[1]*.995,x1 = par("usr")[1],y0 = pd$y0b[tmp], col = "black", length = .1)
  tmp = which(pd[[paste0("hdi90u_bias_",r)]] > par("usr")[2])
  if (length(tmp) > 0)
    arrows(x0 = par("usr")[2]*.995,x1 = par("usr")[2],y0 = pd$y0b[tmp], col = "black", length = .1)
  
  points(sl(pd[[paste0("m_bias_",r)]]),pd$y0b, col = cols[1], pch = pchs[1],cex = 1.5, bg = bgs[1])
  points(sl(pd[[paste0("m_biasm_",r)]]),pd$y0a, col = cols[3], pch = pchs[3],cex = 1.5, bg = bgs[3])
  
  abline(h = (n-1-.5)+c(-.03,.03),xpd = F, col = "grey")
  abline(h = (n-1-6-.5)+c(-.03,.03),xpd = F, col = "grey")
  
  legend(mean(par("usr")[1:2]),
         par("usr")[4]+.5,
         col = c("black","blue"),
         pch = pchs[2:3],
         pt.bg = bgs[2:3],
         pt.lwd = 1,
         pt.cex = c(1.5,1.5),
         cex = 1.25,
         legend = c(ifelse(ft == "tex",paste0("$\\frac{\\delta_{",rr,"}}{\\sigma_{IPPW}}$"),expression(paste(delta[rr],"/",sigma["IPPW"]))),
                    ifelse(ft == "tex",paste0("$\\frac{\\delta_{",rr,"}}{\\mu_{IPPW}}$"),expression(paste(delta[rr],"/",mu["IPPW"])))),
         lwd = c(0,0),
         lty = c(0,0),
         bty = "n",
         xpd = T, text.width = 2,
         ncol = 2,
         xjust = 0.5,
         seg.len = 4,
         x.intersp = 0)
   
  
  
}


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
         ylim = c(0,1),
         col = col,
         xaxt = "n",
         yaxt = "n",
         bty = "n")
    title(names(p_in_rope_bias)[k],cex = .5)
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

for (r in c("ar", "raw")) {
  if (r == "ar") {
    p_in_rope_bias = p_in_rope_bias_ar
    p_in_rope_biasm = p_in_rope_biasm_ar
    rr = "AR"
  } else {
    p_in_rope_bias = p_in_rope_bias_raw
    p_in_rope_biasm = p_in_rope_biasm_raw
    rr = "UR"
  }
  
  names(p_in_rope_bias) = gsub("small_fga","SGA",names(p_in_rope_bias) )
  names(p_in_rope_bias) = gsub("SCORE","score",names(p_in_rope_bias) )
  names(p_in_rope_bias) = gsub("TYPUNITS","typunits",names(p_in_rope_bias) )
  names(p_in_rope_bias) = gsub("EFFUNITS","effunits",names(p_in_rope_bias) )
  names(p_in_rope_bias) = gsub("FREQ","freq",names(p_in_rope_bias) )
  
  if(ft == "tex") {
    tikz(file = paste0("LTX/figures/ROPE-plots-",rr,".tex"), width = 15/2.54, height = 15/2.54, standAlone = F)  
  } else {
    postscript(paste0("LTX/figures/eps/ROPE-plots-",rr,".eps"), width = 17.4/2.54, height = 17.4/2.54)
  }
  
  par(mar=c(2,2,1.5,.01), mgp=c(1,.2,0), tck=.01, mfrow = c(4,4))
  names(p_in_rope_bias) = gsub("\n"," ",clean_varnames)
  for (k in 1:length(p_in_rope_bias)) {
    plot_ropehdi(p_in_rope_bias[[k]],
                 xlab = ifelse(k>12,T,F),
                 ylab = ifelse(k %in% c(1,5,9,13),T,F))
    plot_ropehdi(p_in_rope_biasm[[k]],add = T)
  }
  dev.off() 
  
  
  if (ft == "tex") {
    tikz(file = paste0("LTX/figures/logRRs-",rr,".tex"), width = 15/2.54, height = 18/2.54, standAlone = F)
  } else {
    postscript(file = paste0("LTX/figures/eps/logRRs_",rr,".eps"), width = 15/2.54, height = 18/2.54)
  }
  
  
  par (mar=c(3,8,1,1),xpd = F,mgp=c(1.75,.5,0), tck=.01, mfrow = c(1,1))
  plot(0, type = "n",
       ylim = range(pd[,grep("y",names(pd)),with = F]),
       xlim = log(100)*c(-1,1),
       ylab = "",
       xlab = "log(RR) [bias in ROPE]",
       yaxt = "n", bty = "n", xaxt = "n")
  rect(-log(3),par("usr")[3],log(3),par("usr")[4],col = "gray88", border = NA)
  xticks = sort(c(-log(c(3,10,30,100)),log(c(3,10,30,100))))
  axis(1,at = xticks, labels = exp(abs(xticks))*sign(xticks))
  abline(v = xticks,lty = 3, col = "grey")
  
  #abline(v = 0,col = "red",lwd = 2)
  mtext("bias",side = 3, line = 0,adj = .25,col = "red")
  mtext("no bias",side = 3, line = 0,adj = .75,col = "green3")
  par(mgp=c(1.75,.7,0))
  axis(2,at = pd$y, labels = lbs, las = 2,lty = 0)
  offset = .15
  xs = pd[,get(paste0("logRR_biasm_",r))]
  osr = abs(xs) > par("usr")[2]
  ss = sign(xs[osr])
  xs[osr] = ss*par("usr")[2]*.95
  segments(x0 = 0,
           x1 = xs[!osr],
           y0 = pd$y[!osr]-offset,
           lwd = 12,
           lend = 1,
           col = "blue")
  if (any(osr)) {
    arrows(x0 = 0,
           x1 = xs[osr],
           y0 = pd$y[osr]-offset,
           lwd = 12,
           lend = 1,
           col = "blue") 
  }
  
  xs = pd[,get(paste0("logRR_bias_",r))]
  osr = abs(xs) > par("usr")[2]
  ss = sign(xs[osr])
  xs[osr] = ss*par("usr")[2]*.95
  segments(x0 = 0,
           x1 = xs[!osr],
           y0 = pd$y[!osr]+offset,
           lwd = 12,
           lend = 1,
           col = "black")
  if (any(osr)) {
    arrows(x0 = 0,
           x1 = xs[osr],
           y0 = pd$y[osr]+offset,
           lwd = 12,
           lend = 1)
  }
  par(xpd = T)
  text(mean(log(c(3,10))),3,"moderate",srt = 90)
  text(mean(log(c(10,30))),3,"strong",srt = 90)
  text(mean(log(c(30,100))),3,"very strong",srt = 90)
  text(mean(log(c(100,exp(par("usr")[2])))),3,"extreme",srt = 90)
  
  legend("topright",
         pch = 15,
         col = c("blue","black"),
         legend = c(ifelse(ft == "tex","$\\delta/\\sigma_{IPPW}$",expression(paste(delta,"/",sigma["IPPW"]))),
                    ifelse(ft == "tex","$\\delta/\\mu_{IPPW}$",expression(paste(delta,"/",mu["IPPW"])))),
         box.col = "white")
  dev.off()
}






