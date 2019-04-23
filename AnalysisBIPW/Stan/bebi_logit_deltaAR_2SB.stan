functions { 
   int get_iter(int id);
}

data { 
  int<lower=1> N;      // number of observations 
  int Y[N];            // response variable 
  int<lower=1> K;      // number of predictors for AR model
  matrix[N, K] X;      // design matrix
  int<lower=1> K_raw;  // number of predictors for raw model
  int<lower=1> K_ipw;  // number for IPW model
  int ipw_idx[N];      // weight index 
  matrix[70,2000] ipw_weights;
  int id;  
}

transformed data {
  vector<lower=0>[N] weights[1250]; 
  matrix[N,K] X0[K_ipw];
  matrix[N,K] X1[K_ipw];
  matrix[N,K_ipw] X_ipw;
  matrix[N,K_raw] X_raw;
  
  X_ipw = X[,1:K_ipw];
  X_raw = X[,1:K_raw];

  for (k in 1:K_ipw) {
    X0[k] = X;
    X1[k] = X;
    for (n in 1:N) X0[k][n,k] = 0.0;
    for (n in 1:N) X1[k][n,k] = 1.0;
  }
  for (i in 1:1250) weights[i] = ipw_weights[ipw_idx,i];  
}

parameters { 
  vector[K_ipw] b_ipw;
  vector[K_ipw] delta_ar;
  vector[K_raw] delta_raw;
  vector[K-K_ipw] b_ar_r;
  real<lower=0, upper=1> theta_ipw;
  real<lower=0, upper=1> theta_ar;
  real<lower=0, upper=1> theta_raw;
  real alpha_ipw;
  real alpha_ar;
  real alpha_raw;
}

transformed parameters {
  vector[K] b_ar;
  vector[K_raw] b_raw;
  real<lower = 0> phi_ipw = (1-theta_ar)/theta_ipw;
  real<lower = 0> phi_ar = (1-theta_ar)/theta_ar;
  real<lower = 0> phi_raw = (1-theta_raw)/theta_raw;
  b_ar = append_row(b_ipw + delta_ar,b_ar_r);
  b_raw = b_ipw[1:K_raw] + delta_raw;
}

model { 
  b_ipw ~ normal(0,2);
  b_ar ~ normal(0,2);
  b_raw ~ normal(0,2);
  alpha_ipw ~ normal(0,2);
  alpha_ar ~ normal(0,2);
  alpha_raw ~ normal(0,2);
  
  { 
    vector[N] p_ar = inv_logit(X * b_ar + alpha_ar);
    vector[N] p_raw = inv_logit(X_raw * b_raw + alpha_raw);
    vector[N] p_ipw = inv_logit(X_ipw * b_ipw + alpha_ipw);
    vector[N] lp_ipw;
    
    target += beta_binomial_lpmf(Y | 22, p_ar * phi_ar, (1-p_ar) * phi_ar);
    target += beta_binomial_lpmf(Y | 22, p_raw * phi_raw, (1-p_raw) * phi_raw);
    for (n in 1:N) {
      lp_ipw[n] = beta_binomial_lpmf(Y[n] | 22, p_ipw[n] * phi_ipw, (1- p_ipw[n]) * phi_ipw);
      }
    target += dot_product(weights[get_iter(id)+25], lp_ipw);
  }
} 

generated quantities {
  vector[K_raw] me_ar;
  vector[K_raw] me_raw;
  vector[K_raw] me_ipw;
  vector[K_raw] delta_ipw_raw;
  vector[K_raw] delta_ipw_ar;
  vector[K_raw] delta_ar_raw;
  
  for (k in 1:K_raw) {
    me_ipw[k] = mean(inv_logit(X0[k][,1:K_ipw] * b_ipw + alpha_ipw) - 
                     inv_logit(X1[k][,1:K_ipw] * b_ipw + alpha_ipw));
    me_raw[k] = mean(inv_logit(X0[k][,1:K_raw] * b_raw + alpha_raw) - 
                     inv_logit(X1[k][,1:K_raw] * b_raw + alpha_raw));
    me_ar[k] = mean(inv_logit(X0[k] * b_ar + alpha_ar) - 
                    inv_logit(X1[k] * b_ar+ alpha_ar));
  }
  delta_ipw_raw = me_ipw - me_raw;
  delta_ipw_ar = me_ipw - me_ar;
  delta_ar_raw = me_ar - me_raw;
}

