create_sim_config <- function(sim_jobs, fit_data, 
                              grid = list(N_min = 22L, N_max = 128L, N_unit = 22L),
                              stop_at = list(time = 60, size = 64),
                              sd_g = 0.15, eps = 1e-16) {
  
  # --- Define Likelihood Functions Internal to Config ---
  ll_g <- function(g_obs, g_model, sd_g) {
    if (is.null(g_obs) || length(g_obs) == 0L || !is.finite(g_obs)) return(0)
    dnorm(g_obs, mean = g_model, sd = sd_g, log = TRUE)
  }
  
  ll_kary <- function(kary_obs, dist_df, N_unit, eps) {
    if (is.null(kary_obs) || length(kary_obs) == 0L) return(0)
    kary_obs <- kary_obs[is.finite(kary_obs)]
    if (!length(kary_obs)) return(0)
    
    p_by_N <- dist_df$fraction
    names(p_by_N) <- as.character(dist_df$N)
    
    N_hat <- as.integer(round(kary_obs * N_unit))
    # Clamp to grid range
    range_N <- as.numeric(names(p_by_N))
    N_hat <- pmax(min(range_N), pmin(max(range_N), N_hat)) 
    
    sum(log(p_by_N[as.character(N_hat)] + eps))
  }
  
  ll_flow <- function(flow_obs, model_res) 0 # Placeholder as in original
  
  # Wrapper to score a single ID
  ll_one_id <- function(dat, model_res, grid_opts, global_opts) {
    ll_g(dat$g, model_res$growth_rate, global_opts$sd_g) +
      ll_flow(dat$flow, model_res) +
      ll_kary(dat$kary, model_res$distribution, grid_opts$N_unit, global_opts$eps)
  }
  
  # --- Return Bundle ---
  list(
    jobs = sim_jobs,
    data = fit_data,
    grid = grid,
    opts = list(stop_at = stop_at, sd_g = sd_g, eps = eps),
    funs = list(score_id = ll_one_id)
  )
}

run_sim_jobs_in_memory <- function(theta, config, keep_states = FALSE, print_ll =FALSE) {
  
  # --- 1. Unpack Config ---
  sim_jobs <- config$jobs
  fit_data <- config$data
  grid     <- config$grid
  opts     <- config$opts
  
  # Generate explicit grid vector for Gaussian calculation
  grid_vec <- grid$N_min:grid$N_max
  
  # --- 2. Parameter Transformation ---
  # Assuming theta comes in as log-values (from optim)
  theta_vals <- exp(theta)
  # Ensure theta is a named vector/list for easy lookup
  if(is.null(names(theta_vals))) stop("theta must be a named vector")
  theta_list <- as.list(theta_vals)
  
  # --- 3. Pre-calculation Cache (Physics) ---
  unique_O <- unique(sim_jobs$oxygen)
  system_cache <- vector("list", length(unique_O))
  names(system_cache) <- as.character(unique_O)
  
  for(o_lvl in unique_O) {
    # Note: We pass theta_list. The physics function should ignore 
    # extra start-dist parameters safely if it only looks up what it needs.
    sys <- build_transition_matrix(o_lvl, theta_list, grid)
    I_mat <- diag(nrow(sys$B))
    system_cache[[as.character(o_lvl)]] <- list(
      Q = sys$lambda * (sys$B - I_mat),
      grid = sys$grid,
      lambda = sys$lambda,
      B = sys$B
    )
  }
  
  # --- 4. Simulation Loop ---
  states  <- if (keep_states) vector("list", nrow(sim_jobs)) else NULL
  u_cache <- if (!keep_states) vector("list", nrow(sim_jobs)) else NULL
  if (keep_states) names(states) <- sim_jobs$sim_key
  if (!keep_states) names(u_cache) <- sim_jobs$sim_key
  
  ll_by_job <- numeric(nrow(sim_jobs))
  names(ll_by_job) <- sim_jobs$sim_key
  ll_total <- 0
  
  for (i in seq_len(nrow(sim_jobs))) {
    job <- sim_jobs[i, ]
    sim_key    <- job$sim_key
    parent_key <- job$parent_key
    O          <- job$oxygen
    
    # --- Determine Starting State ---
    start_u <- NULL
    
    if (!is.na(parent_key)) {
      # Inheritance Logic
      source_list <- if(keep_states) states else u_cache
      if (is.null(source_list[[parent_key]])) stop("Parent state missing: ", parent_key)
      
      prev_obj <- source_list[[parent_key]]
      start_u  <- if(keep_states) prev_obj$u_vec else prev_obj
      
    } else {
      # Root Logic: Parameterized Gaussian
      # We look up the specific parameter names defined in the job row
      par_mu_name <- job$start_mean_par
      par_sd_name <- job$start_sd_par
      
      if(is.null(par_mu_name) || is.null(par_sd_name)) {
        stop(paste("Root job", sim_key, "missing start_mean_par or start_sd_par columns"))
      }
      
      mu_val <- theta_list[[par_mu_name]]
      sd_val <- theta_list[[par_sd_name]]
      
      if(is.null(mu_val) || is.null(sd_val)) {
        stop(paste("Parameters", par_mu_name, "or", par_sd_name, "not found in theta"))
      }
      
      # Create Discretized Gaussian
      dist_raw <- dnorm(grid_vec, mean = mu_val, sd = sd_val)
      start_u  <- dist_raw / sum(dist_raw) # Normalize
    }
    
    # Ensure normalized (sanity check)
    if(!is.null(start_u)) start_u <- start_u / sum(start_u)
    
    # --- Run Model ---
    cached_sys <- system_cache[[as.character(O)]]
    
    res <- tryCatch({
      run_ploidy_model(
        O = O, theta = theta_list,
        start_N = NULL, # We provide start_u explicitly now
        start_u = start_u,
        grid = grid,
        stop_at = opts$stop_at,
        method = "expm",
        cached_system = cached_sys
      )
    }, error = function(e) return(NULL))
    
    if (is.null(res)) return(1e12) # Penalize failures
    
    # --- Score Data ---
    ids <- job$data_ids[[1]]
    ll_i <- 0
    for (pid in ids) {
      dat <- fit_data[[pid]]
      # Use the scoring function inside the config
      ll_i <- ll_i + config$funs$score_id(dat, res, grid, opts)
    }
    
    ll_by_job[i] <- ll_i
    ll_total <- ll_total + ll_i
    
    # --- Store State ---
    if (keep_states) {
      states[[sim_key]] <- res
    } else {
      u_cache[[sim_key]] <- res$u_vec
    }
  }
  
  out <- -ll_total
  if(print_ll){
    print(out)
    print(theta_vals)
  }
  
  if (keep_states) {
    out <- list(ll = ll_total, ll_by_job = ll_by_job)
    attr(states, "data_ids") <- setNames(sim_jobs$data_ids, sim_jobs$sim_key)
    out$states <- states
  }
  out
}