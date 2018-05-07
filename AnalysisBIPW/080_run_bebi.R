# getting arguments form command line: 
# https://quantitative-ecology.blogspot.no/2007/08/including-arguments-in-r-cmd-batch-mode.html

#R CMD BATCH --no-save --no-restore '--args a=1 i=1 c=1)' run_zinb.R test.out &
# R --slave --vanilla --file=run_zinb.R --args a=2 i=3 scratch=$SCRATH "s=string with spaces" > output
# Rscript myscript.R $a $i $c 


args=(commandArgs(TRUE))

m = args[1]
a = args[2]
i = args[3]
c = args[4]


library(rstan)

if (m == 1) {
  sm = readRDS("stan/bebi_logit_deltaAR_2SB.rds") 
} else {
  sm = 1
}
load(paste0("data/4Stan/a",a,"_i",i,"x.Rdata"))

setwd(system("echo $SCRATCH", intern = T))

sample_file = paste0("tmp",standata$id,".csv")

sf = sampling(sm,
              data = standata,
              iter = 750,
              warmup = 250,
              chain = 1,
              sample_file = sample_file,
              control = list(max_treedepth = 20)
)

sroot = "/cluster/home/guidopb/programming/R/representativeness/IPW/samples/"

file.copy(sample_file,
          paste0(sroot,
                 "bebi_a",a,"_i",i,"_c",c,"_",
                 ifelse(m == 1,"logit","logit"),"x.csv"),
                 overwrite = T)

save(sf,file = paste0(sroot,
                      "bebi_a",a,"_i",i,"_c",c,"_",
                      ifelse(m == 1,"logit","logit"),"x.Rdata"))
                             

