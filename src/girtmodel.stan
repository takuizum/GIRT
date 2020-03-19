// IRT in Rstan runnning test

data {
  int<lower = 1> N;  //subject
  int<lower = 1> M;  //items
  int y[N,M];  // observation
}

parameters{
  
  real theta[N];
  vector<lower=0>[M] a;
  vector[M] b;
  vector<lower=0, upper=2>[N] phi;
  real mu_b;
  real<lower = 0> sigma_b;
  real mu_a;
  real<lower = 0> sigma_a;
  real mu_phi;
  real<lower = 0> sigma_phi;
}

model{
  // hyper prior
  target += normal_lpdf(mu_b | 0, 1);
  target += normal_lpdf(mu_a | 0, 1);
  target += normal_lpdf(mu_phi | 0, 1);
  target += gamma_lpdf(sigma_a | 2, 2);
  target += gamma_lpdf(sigma_b | 2, 2);
  target += gamma_lpdf(sigma_phi | 2, 2);
  // prior
  target += lognormal_lpdf(a | mu_a, sigma_a);
  target += lognormal_lpdf(phi | mu_phi, sigma_phi);
  target += normal_lpdf(theta | 0,1);
  target += normal_lpdf(b | mu_b, sigma_b);
  
  // model
  for(k in 1:M){
    for(i in 1:N){
      target += bernoulli_logit_lpmf(y[i,k] | a[k]*phi[i]/sqrt(phi[i]^2 + a[k]^2)*(theta[i]-b[k]));
      // target += bernoulli_logit_lpmf(y[i,k] | a[k]*(theta[i]-b[k]));
    }
  }
}

