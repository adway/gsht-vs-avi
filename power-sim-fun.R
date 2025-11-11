source("get-cutoffs-helpers.R")
library(rslurm)
library(avlm)
library(DBCS)

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
  lookup <- readRDS("lookup.RDS") # different on the cluster. refers to a specificup file

  y_0 <- rnorm(N, 0, 1)
  y_1 <- y_0 + ATE
  y_pot <- cbind(y_0, y_1)
  w <- rbinom(N, 1, 0.5)
  y_obs <- (1 - w) * y_pot[, 1] + w * y_pot[, 2]

  cs <- cumsum(w)
  idx <- seq(N / K, N, by = N / K)
  nA <- cs[idx]
  nB <- seq(N / K, N, by = N / K) - nA

  cuts <- get_cuts(lookup, N, K)

  # GROUP SEQUENTIAL T STATISTICS
  y_obs_t <- y_obs[w == 1]
  y_obs_c <- y_obs[w == 0]
  gst_T <- numeric()
  for (k in 1:K) {
    gst_T[k] <- t.test(
      y_obs_t[1:nA[k]],
      y_obs_c[1:nB[k]],
      var.equal = TRUE
    )$statistic
  }
  gst_reject <- abs(gst_T) > cuts
  if (any(gst_reject)) {
    gst_K <- which(gst_reject)[1]
  } else {
    gst_K <- 0
  }

  # ANYTIME VALID LINEAR MODELS TEST
  av_p <- numeric()
  av_p[1] <- 1
  for (i in 2:N) {
    classic_fit = lm(y_obs[1:i] ~ w[1:i] - 1)
    av_fit = av(classic_fit, g = 1)
    coef_mat = summary(av_fit)$coefficients
    if (
      !is.null(coef_mat) && ncol(coef_mat) >= 4 && length(coef_mat[, 4]) > 0
    ) {
      av_p[i] <- coef_mat[1, 4]
    } else {
      av_p[i] <- 1
    }
  }
  if (any(av_p < 0.05)) {
    av_N <- which(av_p < 0.05)[1]
  } else {
    av_N <- 0
  }

  # ANYTIME VALID IPW stuff

  df <- data.frame(w, y_obs)
  treatment = "w"
  response = "y_obs"

  CS = DB_CS(df, treatment, response)

  CS_cross = !(CS$lower <= 0 & CS$upper >= 0)

  if (any(CS_cross == TRUE)) {
    avipw_N <- which(CS_cross == TRUE)[1]
  } else {
    avipw_N <- 0
  }

  return(
    list(
      N = N,
      K = K,
      ATE = ATE,
      t_star = gst_K,
      n_star = av_N,
      ipwn_star = avipw_N,
      cuts = cuts,
      t_stats = gst_T
    )
  )
}
