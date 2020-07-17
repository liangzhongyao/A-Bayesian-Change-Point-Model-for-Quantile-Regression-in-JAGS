jags_qr = "
model{
 for(i in 1:length(y)){
   L[i] <- step(x[i] - cp)

   mu[i] <- alpha + beta[1 + L[i]] * (x[i] - cp)
   w[i] ~ dexp(tau[1 + L[i]])
   
   me[i] <- (1 - 2 * p) / (p * (1 - p)) * w[i] + mu[i]
   pe[i] <- (p * (1 - p) * tau[1 + L[i]]) / (2 * w[i])
   y[i] ~ dnorm(me[i], pe[i])
 }

 #priors for regression
 alpha   ~ dnorm(0, 0.01)
 beta[1] ~ dnorm(0, 0.01)
 beta[2] ~ dnorm(0, 0.01)

 sigma[1] ~ dunif(0, 0.01)
 sigma[2] ~ dunif(0, 0.01)
 tau[1] <- pow(sigma[1], -2)
 tau[2] <- pow(sigma[2], -2)
 
 cp ~ dunif(0, 1.5)
}
"

library(rjags)
library(R2jags)

qt = 0.9
jags_data = list(y = y, x = x, p = qt) # import data
params = c("alpha", "beta", "sigma", "cp", "mu")
fit = jags(model.file = textConnection(jags_qr), data = jags_data, 
             n.iter = 100000, parameters.to.save = params)
