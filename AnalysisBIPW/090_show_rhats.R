library(rstan)
library(bayesplot)
library(shinystan)

# number of imputations analysed
Is = 1:20

delta_plot = vector(length = 3, mode = "list")

tms = c()

for (a in 1:3) {
  rhats = c()
  sfl = vector(mode = "list", length = length(Is))
  for (i in Is) {
    sfl[[which(i == Is)]] = vector(mode = "list", length = 3)
    for (ch in 1:3) {
     fn = paste0("samples/pADHD/bebi_a",a,"_i",i,"_c",ch,"_logitx.Rdata")
     load(fn)
     sfl[[which(i == Is)]][[ch]] = sf
     tms = c(tms,sum(get_elapsed_time(sf)))
    }
    rhats = c(rhats,rhat(sflist2stanfit(sfl[[which(i == Is)]])))
  }
  hist(rhats, main = a)
  sf = sflist2stanfit(do.call(c,sfl))
  post = as.matrix(sf)
  delta_plot[[a]] = mcmc_intervals(post,regex_pars = "delta_me")
  save(sf,file = paste0("samples/pADHD/a",a,".Rdata"))
}
