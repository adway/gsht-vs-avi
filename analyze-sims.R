library(tidyverse)
library(ggplot2)

files <- list.files(
  "~/data/gsht-vs-avi/power_diffeta",
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
lookup <- readRDS("lookup-hyp.RDS")

cuts <- get_cuts(lookup, N = 1000, K = 5)

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
    title = "Sequential Cutoff Boundaries for T Test",
    subtitle = "Subset by K"
  ) +
  theme_minimal(base_size = 14)


##### ATE RESULTS
result_files <- list.files(
  "~/data/gsht-vs-avi/power_t_IPW_diffeta",
  pattern = "\\.RDS$",
  full.names = TRUE
)


### RUN THIS FOR T-TESTS
sim_results <- map_dfr(result_files, function(f) {
  outer <- readRDS(f) # outer is a list of runs

  # Loop over each result list inside the file
  map_dfr(outer, function(x) {
    tibble(
      N = x$N,
      K = x$K,
      ATE = x$ATE,
      t_star = x$gst_star_t,
      n_star = x$avi_star_t,
      t_stats = list(x$t_stats),
      t_cuts = list(x$cuts_t), # store vector in list-column
      file = basename(f) # optional traceability
    )
  })
})

### RUN THIS FOR IPW TESTS
sim_results <- map_dfr(result_files, function(f) {
  outer <- readRDS(f) # outer is a list of runs

  # Loop over each result list inside the file
  map_dfr(outer, function(x) {
    tibble(
      N = x$N,
      K = x$K,
      ATE = x$ATE,
      t_star = x$gst_star_ipw,
      n_star = x$avi_star_ipw[2],
      t_stats = list(x$ipw_stats),
      t_cuts = list(x$cuts_ipw), # store vector in list-column
      file = basename(f) # optional traceability
    )
  })
})

# sim_results <- sim_results[sim_results$N != 20000, ]

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

sim_results$power_gsht <- if_else(sim_results$t_star != 0, 1, 0)
sim_results$power_avi <- if_else(sim_results$n_star != 0, 1, 0)

sim_results_summary <- sim_results %>%
  group_by(N, ATE, K) %>%
  summarize(
    mean_stop_gsht = mean(t_stop),
    mean_stop_avi = mean(n_stop)
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
  mutate(t_stop_n_scale = mean_stop_gsht * (N / K))

# sim_results_joined <- sim_results_joined %>%
#   mutate(n_stop_t_scale_ipw = n_ipw_stop / (N / K))

sim_results_joined$power_avi_t <- if_else(
  sim_results_joined$n_stop <= sim_results_joined$t_stop_n_scale,
  1,
  0
)

# sim_results_joined$power_avi_ipw_t <- if_else(
#   sim_results_joined$n_stop_t_scale_ipw <= sim_results_joined$mean_stop_gsht,
#   1,
#   0
# )

power_summary <- sim_results_joined %>%
  group_by(N, ATE, K) %>%
  summarize(
    mean_stop_gsht = mean(t_stop),
    mean_stop_gsht_nscale = mean(t_stop_n_scale),
    mean_stop_avi = mean(n_stop),
    mean_power_gsht_n = mean(power_gsht),
    mean_power_avi_n = mean(power_avi),
    mean_power_gsht_t = mean(power_gsht_t),
    mean_power_avi_t = mean(power_avi_t)
  )

average_avi <- power_summary %>%
  group_by(N, ATE) %>% # average over K
  summarise(
    mean_power_avi_n = mean(mean_power_avi_n),
    mean_power_avi_t = mean(mean_power_avi_t),
    mean_stop_avi = mean(mean_stop_avi),
    .groups = "drop"
  )

## Make plots

ggplot() +
  # --- TEST 1 (one line per K)
  geom_line(
    data = power_summary,
    aes(
      x = N,
      y = mean_stop_gsht_nscale,
      color = as.factor(K),
      group = K
    )
  ) +
  #  --- TEST 2 (averaged over K)
  geom_line(
    data = average_avi,
    aes(x = N, y = mean_stop_avi),
    linewidth = 1.2,
    color = "black"
  ) +
  facet_wrap(~ATE, scales = "free_y") +
  labs(
    x = "Total Sample Size (N)",
    y = "Power",
    color = "K",
    title = "Mean Power by Avg. Stopping Time of GSHT",
    subtitle = expression("T-Test")
  ) +
  theme_minimal()

# gsht_better <- power_summary[
#   which(power_summary$mean_power_gsht_t > power_summary$mean_power_avi_t),
# ]

# write_csv(power_summary, "full-power-summary.csv")
# write_csv(gsht_better, "gsht-better.csv")
