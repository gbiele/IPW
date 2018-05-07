library(data.table)
library(rstanarm)
options(mc.cores = 4)
library(xtable)

ssb = fread("../SSB/Tabeller_utdanning_2001_2009gb.csv",na.strings = "NA",header = T)
ssb[Age == "<20",Bachelor := 0]
ssb[Age == "<20",Master := 0]
ssb = ssb[Year > 2000 & Gender == "F",]
ssb[Age == "40+", Age := "40-49"]
ssb[, Age := ordered(Age,levels = c("<20", "20-24", "25-29", "30-34", "35-39", "40-49"))]
ssb[,number_children := gsub("\\+","",number_children)]
ssb[,n_unclassified := sum(.SD,na.rm = T), by = 1:nrow(ssb),.SDcols = c("Unknown","Elementary", "High-school", "Bachelor", "Master")]
ssb[, n_unclassified := total_n - n_unclassified]
#ssb[Age == "<20" & is.na(Bachelor), Bachelor := 0]




setnames(ssb,c("Unknown","Elementary", "High-school", "Bachelor", "Master"),
         paste0("edu.",c("Unknown","Elementary", "High-school", "Bachelor", "Master")))
ssb = ssb[,-c("Gender"),with = F]

ssb_l = reshape(ssb,varying = names(ssb)[grep("edu",names(ssb))], direction = "long")
setnames(ssb_l,c("time","edu"),c("edu","N"))
ssb_l[, edu := ordered(edu, levels = c("Unknown","Elementary", "High-school", "Bachelor", "Master"))]

# imputed values of cells ommitted due to small frequencies

ssb_l[, n_edu := as.numeric(edu)]
ssb_l[, n_edu := scale(n_edu)]
ssb_l[, n_children := as.numeric(number_children)]
ssb_l[, n_children := scale(n_children)]
ssb_l[, n_age := as.numeric(Age)]
ssb_l[, n_age := scale(n_age)]
ssb_l[, n_age_s := n_age^2]
ssb_l[, n_children_s := n_children^2]
ssb_l[, n_edu_s := n_edu^2]

nb_fit = stan_glm.nb(N ~ n_age*n_edu*n_children + n_age_s + n_edu_s + n_children_s, data = ssb_l)

newdata = ssb_l[is.na(N),]
newdata[,N := NULL]
p = colMeans(posterior_predict(nb_fit,newdata = newdata))
ssb_l[is.na(N), N.predicted := p]

for (k in which(rowSums(is.na(ssb)) > 0)) {
  tmp = ssb_l[number_children == ssb[k,number_children] & Year == ssb[k,Year] & Age == ssb[k,Age]]
  missings = names(ssb)[which(is.na(ssb[k,]))]
  cols = gsub("edu.","",missings)
  tmp = tmp[edu %in% cols,N.predicted]
  tmp = tmp/sum(tmp) * ssb[k, n_unclassified]
  for (e in missings) ssb[[e]][k] = tmp[e == missings]
}


ssb_l = reshape(ssb,varying = names(ssb)[grep("edu",names(ssb))], direction = "long")
setnames(ssb_l,c("time","edu"),c("edu","N"))

# set parity > 0 for mothers younger than 20 to 0,
# because this are few cases in the population
# which do not exist at all in the moba
#ssb_l[Age == "<20", number_children := "1"]
ssb_l = ssb_l[,list(N = sum(N)),
              by = c("number_children","Age","edu")]

ssb_l[, Age := ordered(as.character(Age))]
ssb_l[, parity := as.numeric(number_children) - 1]
ssb_l[, number_children := NULL]

# allocate people with unknown eduction proportioanlly to the education categories
ssb_l[edu == "Unknown", edu := NA]
for (k in which(is.na(ssb_l[,edu]))) {
  orig_counts = ssb_l[!is.na(edu) & 
                        Age == ssb_l$Age[k] &
                        parity == ssb_l$parity[k],
                      N]
  new_counts = orig_counts + orig_counts/sum(orig_counts)*ssb_l[k,N]
  ssb_l[!is.na(edu) & 
          Age == ssb_l$Age[k] &
          parity == ssb_l$parity[k],
        N := new_counts]
}
# remove rows from which peo?le were allocated to other rows
ssb_l = ssb_l[!is.na(edu)]

rm(newdata, tmp, cols, e, missings, p, k)
ssb_l[, edu := ordered(edu,levels = c("Elementary", "High-school", "Bachelor","Master"))]

# set very small groups to 0
ssb_l[N < 10, N:=0]

save(ssb_l,file = "data/ssb_l.Rdata")

clean_data = function(dt) {
  dt[Age == "<20" & (edu == "Bachelor" | edu == "Master"), edu := "High-school"]
  dt[Age == "20-24" & edu == "Master", edu := "Bachelor"]
  if (exists("mEDU",dt)) {
    dt[Age == "<20" & (edu == "High-school"), mEDU := 4]
    dt[Age == "20-24" & mEDU == 6, mEDU := 5]
  }
  dt[Age == "<20", parity := 0]
  dt[Age == "20-24" & parity == 3, parity := 2]
  dt[Age == "20-24" & edu == "Bachelor" & parity == 2, parity := 1]
  dt[Age == "25-29" & edu == "Master" & parity > 1, parity := 1]
  dt[Age == "40-49" & edu == "Master" & parity == 3, parity := 2]
  dt[Age == "40-49" & edu == "Elementary" & parity == 0, parity := 1]
  if (exists("Parity",dt)) {
    dt[Age == "<20" & (Parity != 1), Parity := 1]
    dt[Age == "20-24" & Parity == 4, Parity := 3]
    dt[Age == "20-24" & edu == "Master", edu := "Bachelor"]
    dt[Age == "20-24" & edu == "Bachelor" & Parity == 3, Parity := 2]
    dt[Age == "25-29" & edu == "Master" & Parity > 2, Parity := 2]
    dt[Age == "40-49" & edu == "Master" & Parity > 3, Parity := 3]
    dt[Age == "20-24" & edu == "High-school" & Parity > 3, Parity := 3]
    dt[Age == "40-49" & edu == "Elementary" & Parity == 1, Parity := 2]
  }
  return(dt)  
}
load("data/ssb_l.Rdata")
ssb_l_cleaned = clean_data(ssb_l)
load("data/imputed_preprocessed_ADHDdata_pADHD_v9s_50.Rdata")
imputed_data = data.table(imputed_data)
imputed_data[,parity := Parity-1]
imputed_data[parity > 3,parity := 3]
imputed_data[, Age := cut(mAge,breaks = c(13,19,24,29,34,39,60), labels = levels(ssb_l$Age))]
imputed_data[, edu := cut(mEDU,breaks = c(.5,1.5,4.5,5.5,6.5),labels = levels(ssb_l$edu),ordered = T)]
imputed_data = clean_data(imputed_data)
moba_l = imputed_data[,list(N = .N/length(unique(imputed_data$imp))), by = c("Age","edu","parity")]


# rowSums(round(xtabs(N~parity + Age ,ssb_l)/sum(ssb_l$N)*100, digits = 1))
# rowSums(round(xtabs(N~parity + Age ,moba_l)/sum(moba_l$N)*100, digits = 1))
# (round(xtabs(N~parity + Age ,moba_l)/sum(moba_l$N)*100-xtabs(N~parity + Age ,ssb_l)/sum(ssb_l$N)*100, digits = 1))
# 
# rowSums(round(xtabs(N~edu + Age ,ssb_l)/sum(ssb_l$N)*100, digits = 1))
# rowSums(round(xtabs(N~edu + Age ,moba_l)/sum(moba_l$N)*100, digits = 1))
# (round(xtabs(N~edu + Age ,moba_l)/sum(moba_l$N)*100-xtabs(N~edu + Age ,ssb_l)/sum(ssb_l$N)*100, digits = 1))

moba_l$group = "moba"
ssb_l$group = "pop"
ssb_l = ssb_l[is.na(N), N := 0]

ipw_data = merge(moba_l[,list(N.moba = sum(N)), by = c("Age","edu","parity")],
                 ssb_l[,list(N.ssb = sum(N)), by = c("Age","edu","parity")],
                 by = c("Age","edu","parity"),
                 all = T)

ipw_data[,P := round(N.ssb/sum(ipw_data$N.ssb,na.rm = T)*1000,digits = 1)]
ipw_data[P < .5]


ipw_data[, ipw_direct := N.ssb/N.moba]
ipw_data[, sipw_direct := (N.ssb/N.moba)*sum(ipw_data$N.moba)/sum(ipw_data$N.ssb)]


save(ipw_data,file = "data/IPW_DATA_pADHD_v9.Rdata")
save(imputed_data,file = "data/imputed_preprocessed_cleaned_ADHDdata_pADHD_v9s.Rdata")

data = rbind(moba_l, ssb_l)


age_edu_p = xtabs(N~edu + Age,ssb_l)/sum(ssb_l$N)
age_edu_p = rbind(age_edu_p,colSums(age_edu_p))
age_edu_p = cbind(age_edu_p,rowSums(age_edu_p))
rownames(age_edu_p)[5] = "All"
colnames(age_edu_p)[7] = "All"

age_edu_m = xtabs(N~edu + Age,moba_l)/sum(moba_l$N)
age_edu_m = rbind(age_edu_m,colSums(age_edu_m))
age_edu_m = cbind(age_edu_m,rowSums(age_edu_m))
rownames(age_edu_m)[5] = "All"
colnames(age_edu_m)[7] = "All"

age_edu = data.table(rbind(age_edu_m,
                           age_edu_p,
                           (age_edu_m*sum(moba_l$N))/(age_edu_p*sum(ssb_l$N))))
age_edu = round(age_edu*100,digits = 1)
grp = ordered(c("MoBa","Population","Coverage"),levels = c("MoBa","Population","Coverage"))
age_edu_table = cbind(data.table(group = sort(rep(grp,5)),
                                 Education = ordered(rep(rownames(age_edu_m),3),
                                                     levels = rownames(age_edu_m)),
                      age_edu))
setkeyv(age_edu_table,c("group","Education"))
age_edu_table[,group := as.character(group)]
age_edu_table$group[(1:15)[-seq(1,15,5)]] = ""
age_edu_table[["<20"]] = as.character(age_edu_table[["<20"]])
age_edu_table[["<20"]] = gsub("NaN","-",age_edu_table[["<20"]])
print.xtable(xtable(age_edu_table,
                    caption = paste0("Proportion of mothers split by age and education in study sample
                                     in study sample  (n = ",
                                     round(sum(moba_l$N)),
                                     ") and background population (Moba, n = ",
                                     round(sum(ssb_l$N)),
                                     ", as well as coverage (\\% participation) of population subgroups in the Moba). 
                                     While around 30\\% of mothers with a masters degree participated, only around 
                                     1\\% of mothers with only elementary school education or less participated."),
                    label = "table:age_edu"),
             include.rownames = F,
             file = "LTX/tables/age_edu_v9.tex")

library(vcd)
data = merge(moba_l,ssb_l,by = c("Age","edu","parity"))
setnames(data,names(data),gsub("x","moba",names(data)))
setnames(data,names(data),gsub("y","pop",names(data)))
data = data[,-grep("group",names(data)), with = F]
data = reshape(data,direction = "long", varying = c("N.moba","N.pop"), timevar = "group")
data[, participation := ifelse(group == "moba","yes","no")]
age_edu  = xtabs(N~edu + Age,ssb_l)
coverages = data.frame(xtabs(N~Age+edu,moba_l)/xtabs(N~Age+edu,ssb_l))$Freq
coverages = coverages/max(coverages,na.rm = T)*.99
coverages[is.na(coverages)] = 0.01
gp = gpar(fill = grey.colors(length(coverages)))
gp$fill = sapply(coverages, function(x) adjustcolor("black",alpha = x))
svg(filename = "LTX/figures/age_edu.svg", width = 10/2.54, height = 10/2.54, pointsize = 10)
mosaic(age_edu,gp = gp,las = 2,rot_labels=c(90,0,90,0))
dev.off()


age_parity_p = xtabs(N~parity + Age,ssb_l)/sum(ssb_l$N)
age_parity_m = xtabs(N~parity + Age,moba_l)/sum(moba_l$N)
age_parity = data.table(rbind(round(age_parity_m*100,digits = 1),
                           round(age_parity_p*100,digits = 1),
                           round(xtabs(N~parity + Age,moba_l)/xtabs(N~parity + Age,ssb_l)*100,digits = 1)))
age_parity_table = cbind(data.table(group = sort(rep(c("Moba","Pop.","xCov"),4)),
                                 parity = ordered(rep(rownames(age_parity_p),3))),
                      age_parity)
#print.xtable(xtable(age_parity_table,
#                    caption = "Mothers\\textquotesingle \\space age and parity in background population (n = 335383) and study sample (n = 47637). Moba = study sample, Pop. = population., xCov = population coverage.",
#                    label = "table:age_parity"),
#             include.rownames = F,
#             file = "4paper/LTX/tables/age_parity.tex")
