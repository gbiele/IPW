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

#sm = stan_model("zinb_deltaARCClog_2SB.stan", allow_undefined = TRUE,
#                includes = paste0('\n#include "',
#                                  file.path(getwd(), 'get_iter.hpp'), '"\n'),
#                auto_write = T)
if (m == 1) {
  sm = readRDS("../stan/zinb_deltaARCC_2SB.rds") 
} else {
  sm = readRDS("../stan/zinb_deltaARCClog_2SB.rds")
}
load(paste0("../data/pADHD_2SB/a",a,"_i",i,".Rdata"))

setwd(system("echo $SCRATCH", intern = T))

sample_file = paste0("tmp",standata$id,".csv")

sf = sampling(sm,
              data = standata,
              iter = 750,
              warmup = 250,
              chain = 1,
              init = list(list(scale = 1e-2)),
              sample_file = sample_file,
              control = list(max_treedepth = 20)
)

sroot = "/cluster/home/guidopb/programming/R/representativeness/IPW/samples/pADHD/"
file.copy(sample_file,
          paste0(sroot,
                 "samples_m1_a",a,"_i",i,"_c",c,"_v9x_",
                 ifelse(m == 1,"lin","log"),"sb.csv"),
                 overwrite = T)

save(sf,file = paste0(sroot,
                      "samples_m1_a",a,"_i",i,"_c",c,"_v9x_",
                      ifelse(m == 1,"lin","log"),"sb.Rdata"))
