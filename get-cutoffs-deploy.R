library(rslurm)
library(tidyverse)
library(reshape2)

source("get-cutoffs-fun.R")

Ns = c(500, 1000, 5000, 8000, 10000)
Ks = c(5, 10, 25, 50, 100)

grid <- expand.grid(N = Ns, K = Ks)
grid <- grid[rep(seq_len(nrow(grid)), each = reps), ]
grid$seed <- sample.int(1e9, nrow(grid))

print(paste("Total simulations:", nrow(grid)))

job <- slurm_apply(
  f = sim_fun,
  params = grid,
  jobname = "get_t_cutoffs",
  nodes = 100,
  add_objects = c("get_t_stats", "get_cuts_from_mc", "get_t_stats_memory"),
  cpus_per_node = 1,
  slurm_options = list(account = "stats_dept1", time = "1:00:00", mem = "64G"),
  submit = TRUE
)

print(job)
