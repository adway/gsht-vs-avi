# Get the t cutoffs for a sequence of dependent t statistics using Monte Carlo.
library(tidyverse)
library(reshape2)

# To be honest about this, we need to specify the number of treatments and controls at each look. When we do the actual t test, we observe one such sequential sample and suspect it comes from a null distribution of a quasi-multivariate T with those degrees of freedom.

K <- 3
M <- 200000
nA <- rep(500, K)
nB <- rep(500, K)

XA <- matrix(rnorm(M * sum(nA), 0, 1), nrow = M)
XB <- matrix(rnorm(M * sum(nB), 0, 1), nrow = M)

alpha_spend <- rep(0.5 / K, K)
bracket = c(0.5, 10)

get_t_stats <- function(XA, XB, nA, nB) {
  K = length(nA)
  endsA <- cumsum(nA)
  endsB <- cumsum(nB)
  T_stats = matrix(NA, M, K)

  for (i in 1:K) {
    Ai = XA[, 1:endsA[i]]
    Bi = XB[, 1:endsB[i]]

    nAi = ncol(Ai)
    nBi = ncol(Bi)

    mA = rowMeans(Ai)
    mB = rowMeans(Bi)

    s2A = rowSums((Ai - mA)^2) / (nAi - 1)
    s2B = rowSums((Bi - mB)^2) / (nBi - 1)
    sp2 = ((nAi - 1) * s2A + (nBi - 1) * s2B) / (nAi + nBi - 2)
    se = sqrt(sp2 * (1 / nAi + 1 / nBi))

    T_stats[, i] = (mA - mB) / se
  }
  return
  T_stats
}

T_stats = get_t_stats(XA, XB, nA, nB)

get_cuts <- function(tmat, alpha_vec, sides = 1) {
  K <- ncol(tmat)
  cuts <- numeric(K)
  crossed <- matrix(FALSE, nrow = nrow(tmat), ncol = K)

  for (i in seq_len(K)) {
    if (i == 1) {
      survivors <- rep(TRUE, nrow(tmat))
    } else {
      survivors <- !apply(crossed[, 1:(i - 1), drop = FALSE], 1, any)
    }

    p_i <- alpha_vec[i] / (1 - sum(alpha_vec[seq_len(i - 1)]))

    if (sides == 1) {
      cuts[i] <- quantile(tmat[survivors, i], 1 - p_i)
      crossed[, i] <- tmat[, i] > cuts[i]
    } else {
      cuts[i] <- quantile(abs(tmat[survivors, i]), 1 - p_i / 2)
      crossed[, i] <- abs(tmat[, i]) > cuts[i]
    }
  }
  cuts
}

cuts <- get_cuts(T_stats, alpha_spend, sides = 1)
