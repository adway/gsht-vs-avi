library(tidyverse)
library(avlm)
source("get-cutoffs.R")

# Scenario 1 -- Constant ATE

## Fixed data
N <- 4000
ATE <- 0.1
y_0 <- rnorm(N, 0, 1)
y_1 <- y_0 + ATE
y_pot <- cbind(y_0, y_1)
w <- rbinom(N, 1, 0.5)
y_obs <- (1 - w) * y_pot[, 1] + w * y_pot[, 2]

# Group sequential tests
K <- 10
global_alpha <- 0.05
alpha_cum = seq.int(
  from = global_alpha / K,
  to = global_alpha,
  length.out = K
)
alpha_spend <- rep(0.05 / K, K)

# GET CUTS

M <- 15000
cs <- cumsum(w)
idx <- seq(N / K, N, by = N / K)
nA <- cs[idx]
nB <- seq(N / K, N, by = N / K) - nA
XA <- matrix(rnorm(M * sum(nA), 0, 1), nrow = M)
XB <- matrix(rnorm(M * sum(nB), 0, 1), nrow = M)
bracket = c(0.5, 10)

T_stats <- get_t_stats(XA, XB, nA, nB, M)
cutoff_stats <- get_cuts_from_mc(T_stats, alpha_vec = alpha_spend, sides = 2)
cuts <- cutoff_stats$cuts

y_obs_t <- y_obs[w == 1]
y_obs_c <- y_obs[w == 0]
t.test(
  y_obs_t[1:nA[1]],
  y_obs_c[1:nB[1]],
  var.equal = TRUE
)$statistic

## Power for regular confidence intervals.

## Power for the asymptotic anytime valid confidence intervals (asymptotic linear models confidence sequence from Michael's paper).

classic_fit = lm(y_obs ~ w - 1)
av_fit = av(classic_fit, g = 1)
summary(av_fit)

### OLD STUFF, MANUAL IMPLEMENTATION
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

av_p = numeric()
av_p[1] = 1
for (i in 2:N) {
  classic_fit = lm(y_obs[1:i] ~ w[1:i] - 1)
  av_fit = av(classic_fit, g = 1)
  av_p[i] <- summary(av_fit)$coefficients[, 4]
}
summary(lm(y_obs ~ w - 1))

## Power from the causal estimands paper from David
