library(rstan)

sm = stan_model("bebi_logit_deltaAR_2SB.stan", allow_undefined = TRUE,
                includes = paste0('\n#include "',
                                  file.path(getwd(), 'get_iter.hpp'), '"\n'),
                auto_write = T)


