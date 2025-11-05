library(tidyverse)
library(ggplot2)

files <- list.files(
  "~/data/gsht-vs-avi",
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

get_cuts(lookup, N = 5000, K = 10)

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
    title = "Sequential Cutoff Boundaries",
    subtitle = "Subset by K"
  ) +
  theme_minimal(base_size = 14)
