// input data
data {
 //number of observations
    int<lower=0> N_obs;
    
    //number of treatment groups
    int<lower=0> n_trt;
    
    // response
    real y_obs[N_obs];
    
    //design matrix for the whole treatment
    matrix[N_obs, n_trt] trt;
    
    //time vector
    vector[N_obs] time;
}


// The parameters accepted by the model.
parameters {
  
  //parameters for the Verhulst mmodel
  //vector<lower=0>[n_trt] b;
  vector<lower=0>[n_trt] c;
  vector<lower=0>[n_trt] r; //a
  real<lower=0> N0;
  real<lower=0> sigma_N;  

}


transformed parameters {

  vector[N_obs] yhat_dens; //density mean value
  vector[N_obs] numerator;
  vector[N_obs] denominator;
  vector[N_obs] c_trt;
  vector[N_obs] r_trt;
  
  // for abundance
  r_trt = trt * r;
  c_trt = trt * c;
  
  
  //mean function for trait and density
  for (i in 1:N_obs) {
    //mean function for density
    yhat_dens[i] = N0 * exp(-c_trt[i] * (exp(r_trt[i] * time[i]) - 1));
    
  }
    
}

// The model to be estimated.

model {
  
  sigma_N ~ cauchy(0, 5); //prior for variance of abundance
  r ~ std_normal(); //rs
  //b ~ normal(0, 5); //bs
  c ~ normal(0, 5); //cs
  
  for(i in 1:N_obs)
    y_obs[i] ~ normal(yhat_dens[i], sigma_N);
  
}

generated quantities {
  
  
}