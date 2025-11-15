# Get the t cutoffs for a sequence of dependent t statistics using Monte Carlo.
library(tidyverse)
library(reshape2)

# To be honest about this, we need to specify the number of treatments and controls at each look. When we do the actual t test, we observe one such sequential sample and suspect it comes from a null distribution of a quasi-multivariate T with those degrees of freedom.

K <- 5
M <- 15000
nA <- rep(50, K)
nB <- rep(50, K)

XA <- matrix(rnorm(M * sum(nA), 0, 1), nrow = M)
XB <- matrix(rnorm(M * sum(nB), 0, 1), nrow = M)

alpha_spend <- rep(0.05 / K, K)
bracket = c(0.5, 10)

get_t_stats <- function(XA, XB, nA, nB, M) {
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

get_t_stats_memory <- function(nA, nB, M) {
  K = length(nA)
  endsA = cumsum(nA)
  endsB = cumsum(nB)
  T_stats = matrix(NA, M, K)

  for (i in 1:M) {
    XA = rnorm(sum(nA), 0, 1)
    XB = rnorm(sum(nB), 0, 1)

    for (j in 1:K) {
      Ai = XA[1:endsA[j]]
      Bi = XB[1:endsB[j]]

      nAi = endsA[j]
      nBi = endsB[j]

      mA = mean(Ai)
      mB = mean(Bi)

      s2A = sum((Ai - mA)^2) / (nAi - 1)
      s2B = sum((Bi - mB)^2) / (nBi - 1)
      sp2 = ((nAi - 1) * s2A + (nBi - 1) * s2B) / (nAi + nBi - 2)
      se = sqrt(sp2 * (1 / nAi + 1 / nBi))

      T_stats[i, j] = (mA - mB) / se
    }
  }
  return(T_stats)
}

get_ipw_stats_memory <- function(nA, nB, M) {
  K = length(nA)
  endsA = cumsum(nA)
  endsB = cumsum(nB)
  IPW_stats = matrix(NA, M, K)
  for (i in 1:M) {
    XA = rnorm(sum(nA), 0, 1)
    XB = rnorm(sum(nB), 0, 1)

    for (j in 1:K) {
      Ai = XA[1:endsA[j]]
      Bi = XB[1:endsB[j]]

      nAi = endsA[j]
      nBi = endsB[j]

      mA = sum(Ai) / 0.5
      mB = sum(Bi) / 0.5

      IPW_stats[i, j] = (mA - mB) / (nAi + nBi)
    }
  }
  return(IPW_stats)
}

ipw_stats <- get_ipw_stats_memory(nA, nB, M)

T_stats <- get_t_stats(XA, XB, nA, nB, M)
T_stats <- get_t_stats_memory(nA, nB, M)

get_cuts_from_mc <- function(tmat, alpha_vec, sides = 1) {
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
      cuts[i] <- quantile(abs(tmat[survivors, i]), 1 - p_i)
      crossed[, i] <- abs(tmat[, i]) > cuts[i]
    }
  }
  return
  list("crossed" = crossed, "cuts" = cuts)
}

cutoffs <- get_cuts_from_mc(T_stats, alpha_spend, sides = 2)
cuts <- cutoffs$cuts
crossed <- cutoffs$crossed
cuts

#check for each k
num_rejects = rowSums(crossed)
check = vector()
for (i in 2:K) {
  subset_crossed = crossed[, 1:i]
  num_rejects = rowSums(subset_crossed)
  in_idx = which(num_rejects > 0)
  check[i] = length(in_idx) / length(num_rejects)
}
check

# See if the cutoffs work for a new draw of null data.

crossed_new = matrix(FALSE, nrow = nrow(ipw_stats), ncol = K)
for (i in 1:K) {
  crossed_new[, i] <- abs(ipw_stats[, i]) > cuts[i]
}

num_rejects = rowSums(crossed_new)
check = vector()
for (i in 2:K) {
  subset_crossed = crossed_new[, 1:i]
  num_rejects = rowSums(subset_crossed)
  in_idx = which(num_rejects > 0)
  check[i] = length(in_idx) / length(num_rejects)
}
check

get_cutoffs <- function(y_obs, w, K, M, alpha_spend) {
  N <- length(w)
  cs <- cumsum(w)
  idx <- seq(N / K, N, by = N / K)
  nA <- cs[idx]
  nB <- seq(N / K, N, by = N / K) - nA
  XA <- matrix(rnorm(M * sum(nA), 0, 1), nrow = M)
  XB <- matrix(rnorm(M * sum(nB), 0, 1), nrow = M)
  bracket = c(0.5, 10)

  T_stats <- get_t_stats(XA, XB, nA, nB, M)

  cutoff_stats <- get_cuts_from_mc(T_stats, alpha_spend, sides = 1)

  cuts <- cutoff_stats$cuts
  return(cuts)
}
