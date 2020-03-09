// IRT in Rstan runnning test

data {
    int<lower = 1> N;  //subject
    int<lower = 1> M;  //items
    int y[N,M];  // observation
}

parameters{
    vector[2] theta[N];
    vector[2] mu[N];
    cholesky_factor_corr[2] L[N];
    vector<lower=0>[M] a;
    vector[M] d;
}

model{
    // prior
    for(j in 1:M){
        target += lognormal_lpdf(a[j] | 0, 2);
        target += normal_lpdf(d[j] | 0, 2);
    }
    for(i in 1:N){
        target += lkj_corr_cholesky_lpdf(L[i] | 10);
        target += normal_lpdf(mu[i,1] | 0, 1);
        target += normal_lpdf(mu[i,2] | 0, 1);
        target += multi_normal_cholesky_lpdf(theta[i,] | mu[i], L[i]);
    }
    // model
    for(k in 1:M){
        for(i in 1:N){
            target += bernoulli_logit_lpmf(y[i,k] | a[k]/sqrt(1+exp(theta[i, 2])^2*a[k]^2)*theta[i, 1] + d[k]);
        }
    }
}

