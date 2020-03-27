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
  vector<lower=0>[n_trt] K;
  vector<lower=0>[n_trt] r;
  real<lower=0> N0;
  real<lower=0> sigma_N;  

}


transformed parameters {

  vector[N_obs] yhat_dens; //density mean value
  vector[N_obs] numerator;
  vector[N_obs] denominator;
  vector[N_obs] K_trt;
  vector[N_obs] r_trt;
  
  // for abundance
  K_trt = trt * K;
  r_trt = trt * r;
  numerator = N0 * K_trt;
  
  //mean function for trait and density
  for (i in 1:N_obs) {
    //mean function for density
    denominator[i] = N0 + ((K_trt[i] - N0) * exp((-1 * r_trt[i]) * time[i]));
    yhat_dens[i] = numerator[i] / denominator[i];
    
  }
    
}

// The model to be estimated.

model {
  
  sigma_N ~ cauchy(0, 5); //prior for variance of abundance
  r ~ std_normal(); //rs
  K ~ normal(0, 5); //Ks
  
  for(i in 1:N_obs)
    y_obs[i] ~ normal(yhat_dens[i], sigma_N);
  
}

generated quantities {
  
   vector<lower=0>[n_trt] A;
  
  for(i in 1:n_trt)
    A[i] = r[i] / K[i];
  
}






