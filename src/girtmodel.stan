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
  vector<lower=0>[N] phi;
}

model{
  // 事前分布
  a ~ lognormal(1,1);
  phi ~ lognormal(0,1);
  theta ~ normal(0,1);
  b ~ normal(0,3);
  
  // model
  for(k in 1:M){
    for(i in 1:N){
      y[i,k] ~ bernoulli_logit(1.7*a[k]/sqrt(1+phi[i]^2*a[k]^2)*(theta[i]-b[k]));
      //y[i,k] ~ 1/(1+exp(-1.7*a[k]/sqrt(1+phi[i]^2*a[k]^2)*(theta[i]-b[k])));
    }
  }
}

