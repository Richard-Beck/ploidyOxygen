#!/usr/bin/env Rscript

# ---- Setup and Imports ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error: No array index provided. Usage: Rscript run_buffering_optim.R <index>")
}

array_idx <- as.integer(args[1])
set.seed(array_idx) # Reproducibility based on array index

source("code/optim_utils.R")
source("models/buffering.R")

# Ensure output directory exists
output_dir <- "data/fits/buffering"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- Load Data ----
fit_data <- readRDS("data/fit_objects/fit_data.Rds")
jobs_2N <- readRDS("data/fit_objects/jobs_2N.Rds")
jobs_4N <- readRDS("data/fit_objects/jobs_4N.Rds")

# ---- Configure Jobs ----
# Assign parameter names for root jobs
jobs_2N$start_mean_par <- ifelse(is.na(jobs_2N$parent_key), "mu_2N", NA)
jobs_2N$start_sd_par   <- ifelse(is.na(jobs_2N$parent_key), "sd_2N", NA)

jobs_4N$start_mean_par <- ifelse(is.na(jobs_4N$parent_key), "mu_4N", NA)
jobs_4N$start_sd_par   <- ifelse(is.na(jobs_4N$parent_key), "sd_4N", NA)

jobs <- rbind(jobs_2N, jobs_4N)
jobs <- jobs[order(jobs$depth), ]

# ---- Generate Initial Theta ----
# Helper for Log-Uniform sampling: exp(runif(1, log(min), log(max)))
rlunif <- function(n, min, max) {
  exp(runif(n, log(min), log(max)))
}

# Generate start values based on priors defined in comments
theta <- list(
  p_misseg = rlunif(1, 1e-6, 1e-1), # log prior
  p_wgd    = rlunif(1, 1e-6, 1e-1), # log prior
  beta     = rlunif(1, 0.01, 100),  # log prior
  n_exp    = rlunif(1, 0.1, 10),    # log prior
  smax     = runif(1, 0.1, 1),      # uniform prior
  k_o      = runif(1, 0.01, 20),    # uniform prior
  k_o_mis  = runif(1, 0.01, 20),    # uniform prior
  lam_max  = runif(1, 0.5, 2),      # uniform prior
  lam_min  = runif(1, 0.01, 1),     # uniform prior
  
  # Start Params
  mu_2N    = runif(1, 39, 49),      # uniform prior
  sd_2N    = runif(1, 1, 20),       # uniform prior
  mu_4N    = runif(1, 70, 90),      # uniform prior
  sd_4N    = runif(1, 1, 20)        # uniform prior
)

# ---- Create Config ----
sim_config <- create_sim_config(
  sim_jobs = jobs, 
  fit_data = fit_data
)

# ---- Run Optimization ----
# Optim works in log-space, so we pass log(theta)
start_time <- Sys.time()
cat(sprintf("Starting optimization for array index %d...\n", array_idx))

res <- optim(
  par = log(unlist(theta)), 
  fn = run_sim_jobs_in_memory, 
  config = sim_config
)

end_time <- Sys.time()
duration <- end_time - start_time

# ---- Save Results ----
# Back-transform parameters (exponentiate) for readability in CSV
final_params <- as.list(exp(res$par))
output_df <- data.frame(final_params)

# Add metadata
output_df$nll <- res$value
output_df$convergence <- res$convergence
output_df$seed <- array_idx
output_df$duration_mins <- as.numeric(duration, units = "mins")

output_file <- file.path(output_dir, paste0(array_idx, ".csv"))
write.csv(output_df, output_file, row.names = FALSE)

cat(sprintf("Finished. Saved to %s\n", output_file))