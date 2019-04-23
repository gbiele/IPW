source("load_data_v9_2SB.R")
library(xtable)

list2env(load_data(),
         envir = .GlobalEnv)

my_data = imputed_data


ctrl_vars = grep("cSEX|fEdu|fAge|mBMIdev|mADHD|fADHD",predictors,value = T)
ipw_vars = grep("mAge.cat|mEdu.cat|Parity.cat",predictors,value = T)
predictors = c(setdiff(predictors,ctrl_vars),ctrl_vars)
predictors = c(setdiff(predictors,ipw_vars),ipw_vars)


analyses = vector(mode = "list",length = 3)
# birth variables
birth_vars = grep("small_fga|preterm",predictors,value = T)
analyses[[1]] = list(formula = as.formula(paste0("cADHD ~ ",paste(c(birth_vars,ctrl_vars,ipw_vars),collapse = " + "))),
                   adjust_ipw = F)
# mental health
MH_vars = grep("DRUGS|LTH|SCL",predictors,value = T)
analyses[[2]] = list(formula = as.formula(paste0("cADHD ~ ",paste(c(MH_vars,ctrl_vars,ipw_vars),collapse = " + "))),
                   adjust_ipw = F)
# SmokeDrink_b
SD_vars = grep("SMOKE|ALC.",predictors,value = T)
analyses[[3]] = list(formula = as.formula(paste0("cADHD ~ ",paste(c(SD_vars,ctrl_vars,ipw_vars),collapse = " + "))),
                     adjust_ipw = F)

vt = data.table(Variable = c("cADHD",birth_vars,MH_vars,SD_vars,ctrl_vars, ipw_vars))
vt[Variable %in% birth_vars, Analysis := "Birth"]
vt[Variable %in% MH_vars, Analysis := "Mental Health"]
vt[Variable %in% SD_vars, Analysis := "Legal Drugs"]
vt[Variable %in% c(ctrl_vars), Analysis := "control"]
vt[Variable %in% c(ipw_vars), Analysis := "Pred.part./contr."]
setnames(vt,"Analysis","Variable Group")
docstring = c(cADHD = "child's SS of ADHD symptom ratings",
              Parity = "number of children born", 
              Parity.cat = "number of children born (0, 1, 2, >2)", 
              cSEX = "0 = boys, 1 = girls",
              cADHDz = "standarized ADHD sumscore",
              preterm = "preterm = birthweek < 37",
              small_fga = "among 5% lightes in gest. wk.",
              mAge = "mothers age (scaled)",
              mEdu = "mothers education, 4 ordered categories",
              mAge.cat = "mothers age (categorised)",
              mEDU.cat = "education, 4 ordered categories",
              mINC = "mothers income in 5 ordered categories",
              mSMOKE = "mothers smoking, 1 = yes, 0 = no" ,
              mSMOKE.5Cs = "mothers smoking, 5 cigarettes per day",
              mDRUGS = "used any drug 1 = yes, no = 0",
              mDRUGS.SCORE = "ability score from IRT model",
              mALC  = "drinking alcohol yes = 1, no = 0",
              mALC.FREQ = "drinking per week during pregancy",
              mALC.EFFUNITS = "effective units of alcohol",
              mLTHD = "SS Lifetime History of Depression (scaled)",
              mSCL5 = "SS Symptom Check list (scaled)",
              mBMI = "Body Mass Index",
              mBMIsq = "(log(mBMI)-mean(log(BMI)))^2",
              mBMIdev = "log(mBMI/mBMImode)",
              fAge = "mothers age (categorised)",
              fAge_sq = "scale(fAge)^2",
              fEdu = "fathers education, 4 ordered categories",
              fINC = "mothers income in 5 ordered categories",
              fSMOKE = "fathers smoking, 1 = yes, 0 = no" ,
              fSMOKE.5Cs = "fathers smoking, 5 cigarettes per day",
              fALC.FREQ = "drinking per week during pregancy",
              fALC  = "drinking alcohol yes = 1, no = 0",
              fALC.TYPUNITS = "typical number of units alcohol (scaled)",
              fDRUGS.SCORE = "ability score from IRT model",
              fDRUGS = "used any drug 1 = yes, no = 0",
              CivilStatus = "0=married, 1=cohabitating, 2=Other",
              mADHD = "SS mothers ADHD",
              fADHD = "SS fathers ADHD",
              HELSEREGION = "4 health regions in norway")
vt[,Variable := gsub(".cat","",Variable)]
vt[,Description := docstring[Variable], by = 1:nrow(vt)]
vt[,Variable := gsub("small_fga","SGA",Variable)]
vt = vt[,c("Variable Group", "Variable", "Description")]
vt = vt[grepl("2$",Variable) == F]
vt[, Source := "Q1"]
vt[`Variable Group` == "Birth" | Variable %in% c("cSEX","Parity"), Source := "MBRN"]
vt[grepl("^f",Variable), Source := "QF"]
vt[Variable == "ADHD", `Variable Group` := "Outcome"]
vt[Variable == "ADHD", Source := "Q6"]

for (v in unique(vt$`Variable Group`)[-1]) 
  vt[["Variable Group"]][which(vt$`Variable Group` == v)[-1]] = ""

my_caption = "Description of variables.
	\\newline Q6 = MoBa questionnaire at child age 3, Q1 =  MoBa questionnaire at pregnancy week 20, F = Moba fathers' questionnaire (week 20), MBRN = Medical Birth Registry of Norway, SS = sum score. SGA = small for gestational age. Mothers' smoking and drinking behavior refers to the first 20 weeks of the pregnancy. 
\\newline Outcome = outcome variable, Birth = Birth related exposures, Mental health = exposures related to mothers' and father's mental health, Legal Drugs = exposures from drinking and smoking, control = adjustment variables, Pred.part./ctrl. = predictors of participation (control variables in non-weighted analyses).
Variable names starting with m refer to mothers, f refers to fathers, and c refers to children. All continuous and count variables except parity were re-scaled to have a mean of zero and a standard deviation of one."
print(xtable(vt,
             caption = my_caption,
             label = "table:variables"),
      include.rownames=FALSE,
      file = "LTX/tables/variables.tex")

library(ReporteRs)
doc <- docx()
doc <- addFlexTable(doc,
                    vanilla.table(vt),
                    body.cell.props = cellProperties(padding = 0))
writeDoc(doc, file = "LTX/tables/Table2_variables.docx")


save(analyses,predictors,my_data,ctrl_vars,ipw_vars,file = "data/analyses_pADHD_2SB_v9.Rdata")

library(tikzDevice)

vars = c(birth_vars,MH_vars,SD_vars,ctrl_vars,ipw_vars)
vars = vars[-grep("2",vars)]
cm = cor(my_data[,vars,with = F])
diag(cm) = NA
cm = cm[rev(1:nrow(cm)),]
rownames(cm) = gsub("_","-",rownames(cm))
colnames(cm) = gsub("_","-",colnames(cm))
vars = gsub("_","-",vars)
max_cor = round(abs(max(cm,na.rm = T)),digits = 2)

vars = gsub("5Cs","nCs",vars)
tikz(file = "LTX/figures/covariation.tex", width = 15/2.54, height = 13/2.54)
par(mar = c(6,8,.5,.5),xpd = T)
layout(matrix(c(rep(1,10),2),nrow = 1))
image(x = 1:length(vars),
      y = 1:length(vars),
      z = t(cm),
      col = colorRampPalette(c("blue","white","red"))(101),
      zlim = max_cor*c(-1,1),
      xaxt = "n",
      yaxt = "n",
      bty = "n",
      ylab = "",
      xlab = "")

text(x = par("usr")[1]/2+1:length(vars),
     y = par("usr")[1]*.5,
     labels = vars,
     srt=45,
     pos = 2)
text(y = 1:length(vars),
     x = par("usr")[1],
     labels = rev(vars),
     pos = 2)
segments(y0 = 0.5, y1 = 28.5, x0 = c(2.5,8.5,16.5,22.5), col = "darkgray")
segments(x0 = 0.5, x1 = 28.5, y0 = 26-c(22.5,16.5,8.5,2.5), col = "darkgray")
par(mar = c(6,0,.5,3),xpd = T)
image(y = seq(-.7,.7,length = 101),
      z = matrix(seq(-max_cor, max_cor, length = 101),nrow = 1),
      col = colorRampPalette(c("blue","white","red"))(101),
      xaxt = "n",
      yaxt = "n",
      bty = "n")
axis(4,at = seq(-6,6,by = 2)/10,las = 2,line = NA,lwd = 0,lwd.ticks = 1)
text(0,-.7,pos = 1, labels = "$r$")
dev.off()


tikz(file = "LTX/figures/ADHD_sumscore.tex", width = 12/2.54, height = 7/2.54)
my_data[cADHD > 22, cADHD := 22]
par (mar=c(3,3,.0,1), mgp=c(2,.7,0), tck=-.01)
hist(my_data[imp == 1,cADHD], breaks = c(-.5:22.5), xlab = "ADHD sum score", main = "")
dev.off()