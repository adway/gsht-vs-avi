library(tidyverse)

num_sims = 10^3

# Scenario 1 -- Constant ATE

## Fixed data
N <- 10000
ATE <- 0.5
y_0 <- rnorm(N, 0, 1)
y_1 <- y_0 + 0.5
y_pot <- cbind(y_0, y_1)
w <- rbinom(N, 1, 0.5)
y_obs <- (1 - w) * y_pot[, 1] + w * y_pot[, 2]

# Group sequential tests
K <- 5
global_alpha <- 0.05
alphas = seq.int(from = global_alpha / K, to = global_alpha, length.out = K)

## Power for regular confidence intervals.

## Power for the asymptotic anytime valid confidence intervals (asymptotic linear models confidence sequence from Michael's paper).
g = 1

tau_hats = numeric()
s2 = numeric()
w_inv = numeric()

for (i in 2:N) {
  tau_hat = as.numeric(coef(lm(y_obs[1:i] ~ w[1:i] - 1))[1])
  resid_square = sum((y_obs[1:i] - tau_hat * w[1:i])^2)
  w_n_inv = 1 / (w[1:i] %*% w[1:i])

  s2 = c(s2, resid_square / (i - 1))
  w_inv = c(w_inv, w_n_inv)
  tau_hats = c(tau_hats, tau_hat)
}

se = sqrt(s2) * w_inv

F_n = (tau_hats / se_n)^2
Ns = seq(2:N)
nus = seq(2:N) - 1
G_n = (1 / (1 + Ns))^(1 / 2) *
  ((1 + g / (g + Ns) * (F_n / nus)) / 1 + (F_n / nus))^(-(nus + 1) / 2)

stopping_time = which(G_n < 0.05)[1] + 1
stopping_time

summary(lm(y_obs ~ w - 1))
