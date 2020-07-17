rm(list = ls());gc()
library(rjags)
library(R2jags)
library(tidyverse)
library(readxl)
library(quantreg)
library(gridExtra)

d = read_excel('../1.xlsx')
ggplot(d, aes(x = logCHL, y = logDMS)) +
  geom_point(alpha = 0.5)


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
qt_all = c(0.5, 0.8, 0.9, 0.95, 0.99)
p = list()

for(i in 1:5){
  qt = qt_all[i]
  jags_data = list(y = d$logDMS, x = d$logCHL, p = qt)
  params = c("alpha", "beta", "sigma", "cp", "mu")
  fit = jags(model.file = textConnection(jags_qr), data = jags_data, 
             n.iter = 100000, parameters.to.save = params)
  
  df = d %>%
    mutate(pred = fit$BUGSoutput$mean$mu)
  
  p[[i]] = ggplot(df, aes(x = logCHL)) +
    geom_point(aes(y = logDMS), alpha = 0.5) +
    geom_line(aes(y = pred), alpha = 0.5, lwd = 1.25, col = "red") +
    labs(title = paste("regression quantile = ", qt))
}

ggpubr::ggexport(filename = '../qrcp_fitted.pdf',
                 plot = p, width = 8, height = 6)

qt = 0.95
jags_data = list(y = d$logDMS, x = d$logCHL, p = qt)
params = c("alpha", "beta", "sigma", "cp", "mu")
fit = jags(model.file = textConnection(jags_qr), data = jags_data, 
           n.iter = 100000, parameters.to.save = params)

df = d %>%
  mutate(pred  = fit$BUGSoutput$mean$mu,
         upper = fit$BUGSoutput$mean$mu + fit$BUGSoutput$sd$mu,
         lower = fit$BUGSoutput$mean$mu - fit$BUGSoutput$sd$mu)

ggplot(df, aes(x = logCHL)) +
  geom_point(aes(y = logDMS), alpha = 0.5) +
  geom_line(aes(y = pred),  lwd = 1.25, col = "red") +
  geom_line(aes(y = upper), lwd = 1.25, col = "red", lty = 2) +
  geom_line(aes(y = lower), lwd = 1.25, col = "red", lty = 2) +
  labs(title = paste("regression quantile = ", qt))

boxplot(fit$BUGSoutput$sims.matrix[, 'cp'], xlab = "distribution of change point")


jags_data = list(y = d$logDMS, x = d$pH, p = qt)
params = c("alpha", "beta", "sigma", "cp", "mu")
fit = jags(model.file = textConnection(jags_qr), data = jags_data, 
           n.iter = 100000, parameters.to.save = params)
df = d %>%
  mutate(pred  = fit$BUGSoutput$mean$mu,
         upper = fit$BUGSoutput$mean$mu + fit$BUGSoutput$sd$mu,
         lower = fit$BUGSoutput$mean$mu - fit$BUGSoutput$sd$mu)
ggplot(df, aes(x = pH)) +
  geom_point(aes(y = logDMS), alpha = 0.5) +
  geom_line(aes(y = pred),  lwd = 1.25, col = "red") +
  geom_line(aes(y = upper), lwd = 1.25, col = "red", lty = 2) +
  geom_line(aes(y = lower), lwd = 1.25, col = "red", lty = 2) +
  labs(title = paste("regression quantile = ", qt))
