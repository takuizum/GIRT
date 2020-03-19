data {
    int<lower=1> K; // number of mixture components
    int<lower=1> N; // number of data points
    int<lower=1> P; // number of items
    int y[N,P]; // observations
}

transformed data{
    real h;
    vector[K] theta; // node
    vector[K] log_pi; // weight(log)
    h = 8/(K-1.0);
    for(k in 1:K){
        theta[k] = -4+h*(k-1);
        log_pi[k] = normal_lpdf(theta[k] | 0, 1) + log(h);
    }
}

parameters {
    real<lower=0, upper=5> a[P];
    real d[P];
    //  hyper prior
    real mu_d;
    real<lower=0> sigma_d;
    real mu_a;
    real<lower=0> sigma_a;
}

model {
    // marginarization
    vector[N] log_lik;
    //   vector[N] score;
    {
        vector[2] log_prob[P,K];
        
        // log probability
        for (k in 1:K) {
            for(p in 1:P){
                real z;
                z = a[p] * theta[k] + d[p];
                log_prob[p,k, 1] = log(1 - inv_logit(p)); //  miss
                log_prob[p,k, 2] = log(inv_logit(p)); //  correct
            }
        }
        
        // likelihood
        {
            for (n in 1:N){
                vector[K] ps = log_pi;
                for (k in 1:K){
                    for(p in 1:P){
                        ps[k] = ps[k] + log_prob[p,k,y[n,p] + 1];
                    }
                }
                log_lik[n] = log_sum_exp(ps);
                // score[n] = ps'*theta;
            }
        }
    }
    // prior
    target += normal_lpdf(mu_a | 0, 1);
    target += normal_lpdf(mu_d | 0, 1);
    target += gamma_lpdf(sigma_a | 2, 2);
    target += gamma_lpdf(sigma_d | 2, 2);
    for(j in 1:P){
        target += lognormal_lpdf(a[j] | mu_a, sigma_a);
        target += normal_lpdf(d[j] | mu_d, sigma_d);
    }
    // model
    target += sum(log_lik);
}
