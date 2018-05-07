# see here: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation

# LD scores and data to make sniplist here: 
# https://data.broadinstitute.org/alkesgroup/LDSCORE/

sniplist = fread("GWASresults/sniplist.1000G.EUR.QC.txt")
sniplist = fread("GWASresults/w_hm3.noMHC.snplist.txt")



# cannabis use 
# http://www.nature.com/tp/journal/v6/n3/full/tp201636a.html
dt = fread("C:/Users/gubi/Desktop/GWASresults/CCI_discovery_MA_13_samples_07-18-2015_rs.txt",
           select = c("rsid","Allele1","Allele2","Effect","P-value"))
setnames(dt,"rsid","snpid")
dt = dt[snpid %in% sniplist$SNP]
setnames(dt,"P-value","pval")
setnames(dt,"Allele1","a1")
setnames(dt,"Allele2","a2")
#dt$n = round(251151/2)
dt = dt[,c("snpid","a1","a2","pval","Effect"),with = F]
fwrite(dt,file = "GWASresults/my_data/Cannabis.txt",sep = " ")

# Age at first birth
# http://www.nature.com/ng/journal/v48/n12/full/ng.3698.html
dt = fread("C:/Users/gubi/Desktop/GWASresults/AgeFirstBirth_Female.txt")
setnames(dt,"SNPID","snpid")
dt = dt[snpid %in% sniplist$SNP]
setnames(dt,"Pvalue","pval")
setnames(dt,"A1","a1")
setnames(dt,"A2","a2")
#dt$n = round(251151/2)
dt = dt[,c("snpid","a1","a2","pval","Zscore"),with = F]
fwrite(dt,file = "GWASresults/my_data/AgeFirstBirth.1000G.EUR.txt",sep = " ")

#number children born
# http://www.nature.com/ng/journal/v48/n12/full/ng.3698.html
dt = fread("C:/Users/gubi/Desktop/GWASresults/NumberChildrenEverBorn_Female.txt")
setnames(dt,"SNPID","snpid")
dt = dt[snpid %in% sniplist$SNP]
setnames(dt,"Pvalue","pval")
setnames(dt,"A1","a1")
setnames(dt,"A2","a2")
#dt$n = round(343072/2)
dt = dt[,c("snpid","a1","a2","pval","Zscore"),with = F]
fwrite(dt,file = "GWASresults/my_data/NumChildrBorn.1000G.EUR.txt",sep = " ")


# http://www.ncbi.nlm.nih.gov/pubmed/20418890
# cigarettes ber day
dt = fread("C:/Users/gubi/Desktop/GWASresults/tag.cpd.tbl")
setnames(dt,"SNP","snpid")
dt = dt[snpid %in% sniplist$SNP]
setnames(dt,"P","pval")
setnames(dt,"A1","a1")
setnames(dt,"A2","a2")
#dt$n = 74053
if (min(dt$OR,na.rm = T) < 0 ) dt[,OR := exp(OR)]
dt = dt[,c("snpid","a1","a2","pval","OR"),with = F]
fwrite(dt,file = "GWASresults/my_data/CigPerDay.1000G.EUR.txt",sep = " ")

# ever smoked
dt = fread("C:/Users/gubi/Desktop/GWASresults/tag.evrsmk.tbl")
setnames(dt,"SNP","snpid")
dt = dt[snpid %in% sniplist$SNP]
#dt$n = 74053
setnames(dt,"P","pval")
setnames(dt,"A1","a1")
setnames(dt,"A2","a2")
if (min(dt$OR,na.rm = T) < 0 ) dt[,OR := exp(OR)]
dt = dt[,c("snpid","a1","a2","pval","OR"),with = F]
fwrite(dt,file = "GWASresults/my_data/EverSmoked.1000G.EUR.txt",sep = " ")

# years education
# http://www.nature.com/nature/journal/vaop/ncurrent/full/nature17671.html
# n = 293723
dt = fread("C:/Users/gubi/Desktop/GWASresults/EduYears_Main.txt")
setnames(dt,"MarkerName","snpid")
dt = dt[snpid %in% sniplist$SNP]
setnames(dt,"Pval","pval")
setnames(dt,"A1","a1")
setnames(dt,"A2","a2")
dt = dt[,c("snpid","a1","a2","pval","Beta"),with = F]
fwrite(dt,file = "GWASresults/my_data/YearsEdu.1000G.EUR.txt",sep = " ")


# depressive symptoms
# http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3552.html#access
# n = 161460
dt = fread("C:/Users/gubi/Desktop/GWASresults/DS_Full.txt")
setnames(dt,"MarkerName","snpid")
dt = dt[snpid %in% sniplist$SNP]
setnames(dt,"Pval","pval")
setnames(dt,"A1","a1")
setnames(dt,"A2","a2")
dt = dt[,c("snpid","a1","a2","pval","Beta"),with = F]
fwrite(dt,file = "GWASresults/my_data/DeprSymp.1000G.EUR.txt",sep = " ")


# anxiety (n = 18000)
# https://www.ncbi.nlm.nih.gov/pubmed/26754954
dt = fread("C:/Users/gubi/Desktop/GWASresults/anxiety.meta.full.fs.tbl")
setnames(dt,"SNPID","snpid")
dt = dt[snpid %in% sniplist$SNP]
setnames(dt,"P.value","pval")
setnames(dt,"Allele1","a1")
setnames(dt,"Allele2","a2")
dt = dt[,c("snpid","a1","a2","pval","Effect"),with = F]
fwrite(dt,file = "GWASresults/my_data/Anxiety.1000G.EUR.txt",sep = " ")

# ADHD 
# https://www.biorxiv.org/content/early/2017/06/03/145581
# n = 55374 , 20183 ADHD cases and 35191
# downlowd kaviar data base here: http://db.systemsbiology.net/kaviar/Kaviar.downloads.html
tmp = fread("C:/Users/gubi/Desktop/Kaviar-160204-Public-hg19.vcf",
            select = c("#CHROM","POS","ID"))[ID != ".",]
setnames(tmp,"#CHROM","CHROM")
tmp = tmp[CHROM != "M" & CHROM != "X" & CHROM != "Y"]

dt = fread("C:/Users/gubi/Desktop/GWASresults/pcg_adhd_jul2017.txt")
dt[,POS := strsplit(MarkerName[1],":")[[1]][2],,by = 1:nrow(dt)]
dt[,POS := as.integer(POS)]
dt[,CHROM := strsplit(MarkerName[1],":")[[1]][1],,by = 1:nrow(dt)]
dt[,POS := as.numeric(POS)]

dt = merge(dt,tmp,by = c("CHROM","POS"),all.x = T, all.y = F)
dt = dt[!is.na(ID)]
dt = dt[SNP %in% sniplist$SNP]
setnames(dt,"SNP","snpid")


setnames(dt,"P","pval")
setnames(dt,"A1","a1")
setnames(dt,"A2","a2")
dt = dt[,c("snpid","a1","a2","pval","OR"),with = F]
fwrite(dt,file = "GWASresults/my_data/ADHD_PCG_LUNDBECK_2017.txt",sep = " ")


# Alcohol use, British BioBank
# http://www.nature.com/mp/journal/v22/n10/full/mp2017153a.html
dt = fread("C:/Users/gubi/Desktop/GWASresults/ALC_CONS_GWAS_FORUKB.txt",
           select = c("SNP","A1","A2",
                      "BETA",
                      "SE",
                      "P"))
# n = 112117
setnames(dt,"P","pval")
setnames(dt,"SNP","snpid")
setnames(dt,"BETA","beta")
setnames(dt,c("A1","A2"),c("a1","a2"))
dt = dt[,c("snpid","a1","a2","pval","beta"),with = F]
dt = dt[snpid %in% sniplist$SNP]
fwrite(dt,file = "GWASresults/my_data/Alcohol_use.txt",sep = " ")


# birth weight
# http://www.ncbi.nlm.nih.gov/pubmed/27680694
# n = 143677
dt = fread("C:/Users/gubi/Desktop/GWASresults/BW3_EUR_summary_stats.txt")
setnames(dt,"rsid","snpid")
dt = dt[snpid %in% sniplist$SNP]
setnames(dt,"p","pval")
setnames(dt,"effect_allele","a1")
setnames(dt,"other_allele","a2")
dt = dt[,c("snpid","a1","a2","pval","beta"),with = F]
fwrite(dt,file = "GWASresults/my_data/BirthWeight.1000G.EUR.txt",sep = " ")



