source("get-cutoffs-helpers.R")

sim_fun <- function(N, K, seed) {
  set.seed(seed)
  M <- 15000
  nA <- rep(N / (2 * K), K)
  nB <- rep(N / (2 * K), K)

  global_alpha <- 0.05
  alpha_cum = seq.int(
    from = global_alpha / K,
    to = global_alpha,
    length.out = K
  )
  alpha_spend <- rep(global_alpha / K, K)
  bracket = c(0.5, 10)
  T_stats <- get_t_stats_memory(nA, nB, M)
  cutoff_stats <- get_cuts_from_mc(T_stats, alpha_vec = alpha_spend, sides = 2)
  cuts <- cutoff_stats$cuts

  return(list(N = N, K = K, cuts = cuts))
}
