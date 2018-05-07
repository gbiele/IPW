functions { 
   int get_iter(int id);
}

data { 
  int<lower=1> N;  // number of observations 
  int Y[N];        // response variable 
  int<lower=1> K;  // number of predictors
  matrix[N, K] X;  // design matrix
  int<lower=1> K_ipw; // number of population-level effects 
  int ipw_idx[N];  // weight index 
  matrix[70,2000] ipw_weights;
  int id;  
}

transformed data {
  vector<lower=0>[N] weights[1250]; 
  matrix[N,K_ipw] X_ipw;
  int didx[K_ipw,K_ipw-1];

  X_ipw = X[,1:K_ipw];
  for (k in 1:K_ipw) {
    int counter = 0;
    for (j in 1:K_ipw) {
      if (j != k) {
        counter = counter + 1;
        didx[k,counter] = j;
      }
    }
  }
  for (i in 1:1250) weights[i] = ipw_weights[ipw_idx,i];  
}

parameters { 
  vector[K_ipw] b_ipw;
  vector[K] b_ar;
  real<lower=0, upper=1> theta_ipw;
  real<lower=0, upper=1> theta_ar;
  real alpha_ipw;
  real alpha_ar;
}

transformed parameters {
  real<lower = 0> phi_ipw = (1-theta_ar)/theta_ipw;
  real<lower = 0> phi_ar = (1-theta_ar)/theta_ar;
  vector[K_ipw] me_ipw;
  vector[K_ipw] me_ar;
  vector[K_ipw] delta_me_p;
 
  for (k in 1:K_ipw) {
    vector[N] m0_ipw = X[,didx[k,]] * b_ipw[didx[k,]] + alpha_ipw;
    vector[N] m1_ipw = m0_ipw + b_ipw[k];
    vector[N] m0_ar  = X[,didx[k,]] * b_ar[didx[k,]] + alpha_ar;
    vector[N] m1_ar  = m0_ar + b_ar[k];
    me_ipw[k] = mean(inv_logit(m0_ipw) - inv_logit(m1_ipw));
    me_ar[k] = mean(inv_logit(m0_ar) - inv_logit(m1_ar));
  }
  delta_me_p = me_ipw - me_ar;
}

model { 
  b_ipw ~ normal(0,2);
  b_ar ~ normal(0,2);
  alpha_ipw ~ normal(0,2);
  alpha_ar ~ normal(0,2);
  delta_me_p ~ normal(0,.7414);
  
  { 
    vector[N] p_ar = inv_logit(X * b_ar + alpha_ar);
    vector[N] p_ipw = inv_logit(X_ipw * b_ipw + alpha_ipw);
    vector[N] lp_ipw;
    
    target += beta_binomial_lpmf(Y | 22, p_ar * phi_ar, (1-p_ar) * phi_ar);
    for (n in 1:N) {
      lp_ipw[n] = beta_binomial_lpmf(Y[n] | 22, p_ipw[n] * phi_ipw, (1- p_ipw[n]) * phi_ipw);
      }
    target += dot_product(weights[get_iter(id)+25], lp_ipw);
  }
} 


