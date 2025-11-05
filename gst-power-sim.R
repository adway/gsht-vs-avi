source("get-cutoffs-helpers.R")
library(rslurm)
library(avlm)

sim_fun <- function(N, K, ATE, seed) {
  set.seed(seed)

  y_0 <- rnorm(N, 0, 1)
  y_1 <- y_0 + ATE
  y_pot <- cbind(y_0, y_1)
  w <- rbinom(N, 1, 0.5)
  y_obs <- (1 - w) * y_pot[, 1] + w * y_pot[, 2]

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

  # ANYTIME VALID TEST
  av_p <- numeric()
  av_p[1] <- 1
  for (i in 2:N) {
    classic_fit = lm(y_obs[1:i] ~ w[1:i] - 1)
    av_fit = av(classic_fit, g = 1)
    av_p[i] <- summary(av_fit)$coefficients[, 4]
  }
  if (any(av_p < 0.05)) {
    av_N <- which(av_p < 0.05)[1]
  } else {
    av_N <- 0
  }

  return(
    c(N = N, K = K, ATE = ATE, t_star = gst_K, n_star = av_N)
  )
}

Ns = c(500, 1000, 5000, 8000, 10000)
Ks = c(5, 10, 20, 50, 100)
ATEs = c(0.01, 0.05, 0.1, 0.5)
Ns = 500
Ks = 10
ATEs = 0.01
reps = 100

grid <- expand.grid(N = Ns, K = Ks, ATE = ATEs)
grid <- grid[rep(seq_len(nrow(grid)), each = reps), ]
grid$seed <- sample.int(1e9, nrow(grid))

print(paste("Total simulations:", nrow(grid)))

job <- slurm_apply(
  f = sim_fun,
  params = grid,
  jobname = "sim_ATE",
  nodes = 25, # configure per cluster limits
  cpus_per_node = 4, # 200 * 4 = 800 jobs in flight
  slurm_options = list(time = "1:00:00", mem = "4G")
)

print(job)
