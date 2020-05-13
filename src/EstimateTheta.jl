using Optim, Statistics, Distributions, StatsFuns

function loglik(X, theta, phi, a, b; theta_prior = Normal(0, 1), phi_prior = LogNormal(0, 1))
    J = length(a)
    l = zero(typeof(theta))
    for j in 1:J
        if X[j] !== missing
            p = logistic((a[j] * phi) / sqrt(a[j] ^ 2 + phi ^ 2) * (theta - b[j]))
            l -= X[j] == 1 ? log(p) : log(1-p)
        end
    end
    return l - logpdf(theta_prior, theta) - logpdf(phi_prior, phi)
end

using RCall

R"""
library(mirt)
data <- expand.table(LSAT7)
fit <- mirt(data, 1, "2PL")
para <- coef(fit, IRTpars = TRUE, simplify = TRUE)$items
"""

@rget data
@rget para

loglik(data[1,:], 1.0, 1.0, para[:, 1], para[:, 2])
opt = optimize(par -> loglik(data[100,:], par[1], par[2], para[:, 1], para[:, 2]), [0.0, 1.0])
opt.minimizer

function EstimatePerson(X, para)
    # 
    N = size(X, 1)
    est = Matrix{Float64}(undef, N, 2)
    for i in 1:N
        #
        opt = optimize(par -> loglik(X[i,:], par[1], par[2], para[:, 1], para[:, 2]), [0.0, 1.0])
        est[i, :] = opt.minimizer[:]
        print(i, " out of ", N, "\r")
    end
    return est
end

est = EstimatePerson(data, para);
mean(est; dims = 1)
sqrt.(var(est; dims = 1))
scatter(est[:, 1], est[:, 2])

using Plots
theta = [-4:0.1:4;];
phi = [0:0.05:4;];
plot(theta, phi, [loglik(data[1,:], theta[i], phi[j], para[:, 1], para[:, 2]) for i in 1:length(theta), j in 1:length(phi)]')
plot(theta, phi, [loglik(data[1,:], theta[i], phi[j], fill(1.0, 5), para[:, 2]) for i in 1:length(theta), j in 1:length(phi)]')
