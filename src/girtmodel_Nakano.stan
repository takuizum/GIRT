data{
    int<lower = 1> n; 
    int<lower = 1> p; 
    int<lower = 0, upper = 1> x[n, p];
}


parameters{ 
    real<lower = 0, upper = 2> psi[p]; 
    real p2[p];
    real theta[n];
    real<lower = 0, upper = 1> phi[n];
}


model{ 
    for(i in 1:p){
        psi[i] ~ uniform(0, 2); 
        p2[i] ~ normal(0, 1);
    }
    for(j in 1:n){ 
        theta[j] ~ normal(0, 1); 
        phi[j] ~ beta(1, 4); 
        for(i in 1:p){ 
            x[j, i] ~ bernoulli_logit((1.7/sqrt(psi[i]^2+((2*phi[j])^2))) * (theta[j] - p2[i])) ;
        }
    }
}

// generated quantities{ 
//     real<lower=0> p1[p]; 
//     real log_lik[n，p]; 
//     for(i in l:p){ 
//         p1[i] = 1/psi[i]; 
//     }
//     for(j in l:n){ 
//         for(i in l:p){ 
//             log_lik[j, i] = bernoulli_logit_lpmf(x[j, i] | (1.7 * sqrt(psi[i]^2+((2*phi[j])^2))) * (theta[j] - p2[i])) 
//         }
//     }
// }