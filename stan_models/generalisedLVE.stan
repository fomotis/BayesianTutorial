//
functions {
  
  real[] biLVE(real t,  //time 
             real[] N,//state or abundance
             real[] tetha, //parameters
             real[] x_r, //data (real)
             int[] x_i //data (integer)
             ) {
               
      real dNdt[2];
      
      
      dNdt[1] =  N[1] * tetha[1] * (1 - (tetha[3]*N[1] + tetha[5]*N[2]) );
      dNdt[2] =  N[2] * tetha[2] * (1 - (tetha[4]*N[1] + tetha[6]*N[2]) );
      return dNdt;
  
  }
  
}

//data
data {
  
  int<lower=0> N_obs; //N_obs
  int<lower=0> T; //number of unique timepoints
  
  vector[2] N[N_obs]; //observed data
  int<lower=0> nsp; //number of species
  
  real t0; // first time point
  real ts[T]; // time to integrate ODE over
  vector[N_obs] time_obs; //observed time point
 
}

transformed data {
  real x_r[0];
  int x_i[0];
}

// The parameters accepted by the model. Our model
parameters {
  
  //growth rate parameters for the LVE model
  real<lower=0> r[2];
  
  //interspecific parameters
  real<lower=0> alpha12; 
  real<lower=0> alpha21; 
  
  //intraspecific parameters
  real<lower = alpha21> alpha11; 
  real<lower = alpha12> alpha22; 
  
  real<lower=0> N0[nsp]; //starting values for both 
  vector<lower=0>[nsp] sigma; //variance for the two species
  cholesky_factor_corr[nsp] L_Omega; //for the correlation matrix
  
}


transformed parameters {
  
    real N_LVE[T, nsp];
    matrix[T, nsp] N_exp;
    vector[nsp] N_hat[N_obs];
    real tetha[6];
  
    tetha[1] = r[1];
    tetha[2] = r[2];
    tetha[3] = alpha11;
    tetha[4] = alpha21;
    tetha[5] = alpha12;
    tetha[6] = alpha22;
  
    N_LVE = integrate_ode_rk45(biLVE, N0, t0, ts, tetha, x_r, x_i);
    N_exp = to_matrix(N_LVE);
    
    //setting the expected N_exp to N_LVE at the required time points
    for(k in 1:T) {
      for(j in 1:N_obs) {
        
        if(ts[k] == time_obs[j]) N_hat[j] = to_vector(N_exp[k, ]);
        
      }
    }
}


model {
  
  matrix[2, 2] L_Sigma;
  
  //priors for the LVE parameters
  
  //interspecific effects
  alpha12 ~ std_normal(); 
  alpha21 ~ std_normal();
  
  r ~ std_normal(); //growth rates
  
  //priors for starting value
  N0 ~ normal(0, 5);
  
  //priors for the variance-covariance matrix of the response in cholesky decomposition form
  L_Omega ~ lkj_corr_cholesky(4);
  L_Sigma = diag_pre_multiply(sigma, L_Omega);  //cholesky form of the variance-covariance matrix
  sigma ~ cauchy(0, 2.5);
  
  
  //and build the log-likelihood
  N ~ multi_normal_cholesky(N_hat, L_Sigma);
  
}

generated quantities {
  
  matrix[2, 2] Sigma;
  matrix[2, 2] L_Sigma;
  real cor;
  
  //reconstructing the variance-covariance matrix from the chole
  L_Sigma = diag_pre_multiply(sigma, L_Omega);
  Sigma = L_Sigma * L_Sigma';
  cor = Sigma[1,2] / pow(Sigma[1,1] * Sigma[2,2], 0.5);
  
}





