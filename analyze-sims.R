library(tidyverse)
library(ggplot2)

files <- list.files(
  "~/data/gsht-vs-avi/ipw_cutoffs",
  pattern = "\\.RDS$",
  full.names = TRUE
)

lookup <- map_dfr(files, function(f) {
  outer <- readRDS(f) # outer is a list of runs

  # Loop over each result list inside the file
  map_dfr(outer, function(x) {
    tibble(
      N = x$N,
      K = x$K,
      cuts = list(x$cuts), # store vector in list-column
      file = basename(f) # optional traceability
    )
  })
})

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

saveRDS(lookup, "lookup-ipw.RDS")

cuts <- get_cuts(lookup, N = 500, K = 5)

lookup_long <- lookup %>%
  unnest_longer(cuts, indices_to = "look") %>% # each cutoff gets its own row
  rename(cutoff = cuts)

ggplot(lookup_long, aes(x = look, y = cutoff, color = factor(N), group = N)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  facet_wrap(~K, scales = "free") +
  labs(
    x = "Look (Interim Analysis Step)",
    y = "Cutoff Value",
    color = "N",
    title = "Sequential Cutoff Boundaries for IPW",
    subtitle = "Subset by K"
  ) +
  theme_minimal(base_size = 14)


##### ATE RESULTS
result_files <- list.files(
  "~/data/gsht-vs-avi/ate_results_with_ipw/",
  pattern = "\\.RDS$",
  full.names = TRUE
)

sim_results <- map_dfr(result_files, function(f) {
  outer <- readRDS(f) # outer is a list of runs

  # Loop over each result list inside the file
  map_dfr(outer, function(x) {
    tibble(
      N = x$N,
      K = x$K,
      cuts = list(x$cuts),
      ATE = x$ATE,
      t_star = x$t_star,
      n_star = x$n_star,
      n_ipw_star = x$ipwn_star,
      t_stats = list(x$t_stats), # store vector in list-column
      file = basename(f) # optional traceability
    )
  })
})

sim_results$t_stop <- if_else(
  sim_results$t_star == 0,
  sim_results$K,
  sim_results$t_star
)

sim_results$n_stop <- if_else(
  sim_results$n_star == 0,
  sim_results$N,
  sim_results$n_star
)

sim_results$n_ipw_stop <- if_else(
  sim_results$n_ipw_star == 0,
  sim_results$N,
  sim_results$n_ipw_star
)

sim_results$power_gsht_n <- if_else(sim_results$t_star != 0, 1, 0)
sim_results$power_avi_n <- if_else(sim_results$n_star != 0, 1, 0)
sim_results$power_avi_n_ipw <- if_else(sim_results$n_ipw_star != 0, 1, 0)

sim_results_summary <- sim_results %>%
  group_by(N, ATE, K) %>%
  summarize(
    mean_stop_gsht = mean(t_stop),
    mean_stop_avi = mean(n_stop),
    mean_stop_avi_ipw = mean(n_ipw_stop),
    # mean_power_gsht_n = mean(power_gsht_n),
    # mean_power_avi_n = mean(power_avi_n),
    # mean_power_avi_n_ipw = mean(power_avi_n_ipw)
  )

sim_results_joined <- sim_results %>%
  left_join(sim_results_summary, by = c("N", "ATE", "K"))

sim_results_joined$power_gsht_t <- if_else(
  sim_results_joined$t_stop <= sim_results_joined$mean_stop_gsht,
  1,
  0
)

sim_results_joined <- sim_results_joined %>%
  mutate(n_stop_t_scale = n_stop / (N / K))

sim_results_joined <- sim_results_joined %>%
  mutate(n_stop_t_scale_ipw = n_ipw_stop / (N / K))

sim_results_joined$power_avi_t <- if_else(
  sim_results_joined$n_stop_t_scale <= sim_results_joined$mean_stop_gsht,
  1,
  0
)

sim_results_joined$power_avi_ipw_t <- if_else(
  sim_results_joined$n_stop_t_scale_ipw <= sim_results_joined$mean_stop_gsht,
  1,
  0
)


power_summary <- sim_results_joined %>%
  group_by(N, ATE, K) %>%
  summarize(
    mean_stop_gsht = mean(t_stop),
    mean_stop_avi = mean(n_stop),
    mean_stop_avi_ipw = mean(n_ipw_stop),
    mean_stop_avi_t_scale = mean(n_stop_t_scale),
    mean_stop_avi_ipw_t_scale = mean(n_stop_t_scale_ipw),
    mean_power_gsht_n = mean(power_gsht_n),
    mean_power_avi_n = mean(power_avi_n),
    mean_power_avi_ipw_n = mean(power_avi_n_ipw),
    mean_power_gsht_t = mean(power_gsht_t),
    mean_power_avi_t = mean(power_avi_t),
    mean_power_avi_ipw_t = mean(power_avi_ipw_t)
  )

gsht_better <- power_summary[
  which(power_summary$mean_power_gsht_t > power_summary$mean_power_avi_t),
]

write_csv(power_summary, "full-power-summary.csv")
write_csv(gsht_better, "gsht-better.csv")
