source("get-cutoffs-helpers.R")
library(rslurm)
library(avlm)
library(DBCS)
library(LambertW)

get_cuts <- function(lookup, N, K) {
  row <- lookup %>% filter(N == !!N, K == !!K)
  if (nrow(row) < 1) {
    stop("No matching parameter combination found.")
  }
  if (nrow(row) > 1) {
    message("Multiple matches; returning first.")
  }
  row$cuts[[1]]
}

sim_fun <- function(N, K, ATE, seed) {
  lookup_hyp <- readRDS("lookup-hyp.RDS") # different on the cluster. refers to a specificup file
  lookup_ipw <- readRDS("lookup-ipw.RDS")

  y_0 <- rnorm(N, 0, 1)
  y_1 <- y_0 + ATE
  y_pot <- cbind(y_0, y_1)
  w <- rbinom(N, 1, 0.5)
  y_obs <- (1 - w) * y_pot[, 1] + w * y_pot[, 2]

  cs <- cumsum(w)
  idx <- seq(N / K, N, by = N / K)
  n_treat <- cs[idx]
  n_control <- idx - n_treat

  gstT_cuts <- get_cuts(lookup_hyp, N, K)
  gstIPW_cuts <- get_cuts(lookup_ipw, N, K)

  y_t <- y_obs[w == 1]
  y_c <- y_obs[w == 0]
  # GROUP SEQUENTIAL T STATISTICS
  gstT_stats <- numeric(K)
  for (k in 1:K) {
    gstT_stats[k] <- t.test(
      y_t[1:n_treat[k]],
      y_c[1:n_control[k]],
      var.equal = TRUE
    )$statistic
  }
  gstT_reject <- abs(gstT_stats) > gstT_cuts
  gstT_stop <- if (any(gstT_reject)) which(gstT_reject)[1] else 0

  # ANYTIME VALID LINEAR MODELS TEST
  avT_pvals <- numeric()
  avT_pvals[1] <- 1
  for (i in 2:N) {
    fit <- lm(y_obs[1:i] ~ w[1:i] - 1)
    avfit <- av(fit, g = 1)
    coefmat <- summary(avfit)$coefficients
    if (!is.null(coefmat) && ncol(coefmat) >= 4 && length(coefmat[, 4]) > 0) {
      avT_pvals[i] <- coefmat[1, 4]
    } else {
      avT_pvals[i] <- 1
    }
  }

  avT_stop <- if (any(avT_pvals < 0.05)) which(avT_pvals < 0.05)[1] else 0

  # GROUP SEQUENTIAL IPW TEST
  gstIPW_stats <- numeric(K)
  for (k in 1:K) {
    Ai = y_t[1:n_treat[k]]
    Bi = y_c[1:n_control[k]]

    nAi = n_treat[k]
    nBi = n_control[k]

    mA = sum(Ai) / 0.5
    mB = sum(Bi) / 0.5

    gstIPW_stats[k] <- (mA - mB) / (nAi + nBi)
  }
  gstIPW_reject <- abs(gstIPW_stats) > gstIPW_cuts
  gstIPW_stop <- if (any(gstIPW_reject)) which(gstIPW_reject)[1] else 0

  # ANYTIME VALID IPW TEST

  df <- data.frame(w, y_obs)
  CS <- DB_CS(df, treatment = "w", response = "y_obs", t_opt = 10)
  CS_100 <- DB_CS(df, treatment = "w", response = "y_obs", t_opt = 100)
  CS_1000 <- DB_CS(df, treatment = "w", response = "y_obs", t_opt = 1000)

  tau_hat <- y_obs / 0.5 * (-1)^(1 - w)
  center <- cumsum(tau_hat) / seq_along(tau_hat)

  sig2_hat = (0.5 * y_obs + 0.5 * y_obs)^2 / 0.5^2
  Sn = cumsum(sig2_hat)
  eta = sqrt((-lamW::lambertWm1(-0.05^2 * exp(1)) - 1) / 10)
  width = 1 /
    (1:nrow(df)) *
    sqrt(((Sn * eta^2 + 1) / eta^2) * log((Sn * eta^2 + 1) / 0.05^2))

  lower = center - width
  upper = center + width

  ipw_cross <- !(lower <= 0 & upper >= 0)
  avIPW_stop_exact <- if (any(ipw_cross)) which(ipw_cross)[1] else 0

  ipw_cross <- !(CS$lower <= 0 & CS$upper >= 0)
  avIPW_stop_10 <- if (any(ipw_cross)) which(ipw_cross)[1] else 0

  ipw_cross <- !(CS_100$lower <= 0 & CS_100$upper >= 0)
  avIPW_stop_100 <- if (any(ipw_cross)) which(ipw_cross)[1] else 0

  ipw_cross <- !(CS_1000$lower <= 0 & CS_1000$upper >= 0)
  avIPW_stop_1000 <- if (any(ipw_cross)) which(ipw_cross)[1] else 0

  return(
    list(
      N = N,
      K = K,
      ATE = ATE,
      gst_star_t = gstT_stop,
      avi_star_t = avT_stop,
      gst_star_ipw = gstIPW_stop,
      avi_star_ipw = c(
        avIPW_stop_exact,
        avIPW_stop_10,
        avIPW_stop_100,
        avIPW_stop_1000
      ),
      cuts_t = gstT_cuts,
      cuts_ipw = gstIPW_cuts,
      t_stats = gstT_stats,
      ipw_stats = gstIPW_stats
    )
  )
}
