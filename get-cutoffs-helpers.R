# Get the t cutoffs for a sequence of dependent t statistics using Monte Carlo.
library(tidyverse)
library(reshape2)

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
      cuts[i] <- quantile(abs(tmat[survivors, i]), 1 - p_i / 2)
      crossed[, i] <- abs(tmat[, i]) > cuts[i]
    }
  }
  return
  list("crossed" = crossed, "cuts" = cuts)
}

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
