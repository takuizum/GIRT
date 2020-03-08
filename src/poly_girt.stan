data {
  int<lower = 0> N; // n ob subjects
  int<lower = 0> J; // n of items
  int<lower = 2> K;  // n of categories in each item
  int<lower = 0, upper = K> Y[N, J];
}

parameters {
  vector[N] theta;
  real<lower = 0> phi[N];
  real<lower = 0> a[J];
  ordered[K-1]d[J]; // category of items
  real mu_d;
  real<lower = 0> sigma_d;
}

model {
  a ~ lognormal(0, 1);
  phi ~ lognormal(0, 1);
  theta ~ normal(0, 1);
  for(j in 1:J){
    for(k in 1:(K-1)){
      d[j, k] ~ normal(mu_d, sigma_d);
      // d[j, k] ~ normal(0,1);
    }
  }
  mu_d ~ normal(0, 5);
  sigma_d ~ cauchy(0, 5);
  for(i in 1:N){
    for(j in 1:J){
      Y[i,j] ~ ordered_logistic(1.702 * a[j]/sqrt(1 + phi[i]^2*a[j]^2) * theta[i], 1.702 * a[j]/sqrt(1 + phi[i]^2*a[j]^2)  * d[j]);
       // Y[i,j] ~ ordered_logistic(a[j] * theta[i], d[j]);  // normal GRM
    }
  }
}

