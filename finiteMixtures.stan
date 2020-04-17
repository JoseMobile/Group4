// Stan Test Model

data {
  int<lower=1> N; //sample size
  int<lower=1>K; //number of mixture 
  int<lower=1>p; //dim of a single y elem
   
  cov_matrix[p] V[N]; 
  cov_matrix[p] Omega[K]; 
  int vK[K]; 
  vector[p] Y[N]; // y1, ..., yN
}

parameters {
  simplex[K] rho;
  vector[p] mu[K];
  cov_matrix[p] sigma[K]; 
  vector[p] theta[N];
}

model {
  vector[K] log_rho = log(rho); 
  
  for(i in 1:K){
    sigma[i] ~ inv_wishart(vK[i], Omega[i]); 
  }
  
  for (i in 1:N) {
    vector[K] lps = log_rho;
    for (k in 1:K){
      lps[k] += multi_normal_lpdf(theta[i] | mu[k], sigma[k]);
    }
    target += log_sum_exp(lps);
    
    Y[i] ~ multi_normal(theta[i], V[i]);
  }
}

