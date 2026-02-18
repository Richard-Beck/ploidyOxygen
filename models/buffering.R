suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
  library(tidyr)
  library(expm)
})

# ============================================================
# 1. THE PHYSICS ENGINE
# ============================================================

#' Build the Ploidy Transition Matrix
build_transition_matrix <- function(O, theta, grid) {
  
  N_GRID <- grid$N_min:grid$N_max
  
  # --- Fitness & Parameters ---
  
  # 1. Define a scaling factor (e.g., gamma = 0.33 for SA:V scaling)
  #scaling_factor <- (N_GRID / (2*grid$N_unit))^0.33 
  
  # 2. Scale the Half-Saturation Constant (k_o)
  # Larger cells (higher N) have higher k_o (worse affinity)
  #k_o_adjusted <- theta$k_o * scaling_factor
  
  # 3. Calculate Lambda as a Vector (one rate per ploidy)
  #lambda_val <- theta$lam_min + (theta$lam_max - theta$lam_min) * (O / (O + k_o_adjusted))
  
  lambda_val <- theta$lam_min + (theta$lam_max - theta$lam_min) * (O / (O + theta$k_o))
  p_misseg_adj <- theta$p_misseg * (1 - (O / (O + theta$k_o_mis)))
  
  # --- Internal Helper: Delta Weights ---
  .get_delta_weights <- function(N, p, b, n, sm, k_unit) {
    if (p <= 0 || N <= 0) return(c("0" = 1))
    sd <- sqrt(N * p)
    if (sd == 0) return(c("0" = 1))
    
    sN <- sm * exp(-b * ((2 * k_unit) / N)^n)
    
    z <- 9.0 
    T_range <- min(N, max(0L, ceiling(z * sd)))
    ts <- (-T_range):T_range
    out <- numeric(length(ts))
    
    for (idx in seq_along(ts)) {
      t <- ts[idx]
      ks <- seq.int(abs(t), N, by = 2)
      if (!length(ks)) next
      out[idx] <- sum(dbinom(ks, N, p) * dbinom((ks + t)/2, ks, 0.5) * (sN^ks))
    }
    names(out) <- ts
    out
  }
  
  # --- Build Matrix B ---
  B <- matrix(0, nrow = length(N_GRID), ncol = length(N_GRID), dimnames = list(N_GRID, N_GRID))
  
  for (j in seq_along(N_GRID)) {
    N <- N_GRID[j]
    weights <- .get_delta_weights(N, p_misseg_adj, theta$beta, theta$n_exp, theta$smax, grid$N_unit)
    ts <- as.integer(names(weights))
    pr <- as.numeric(weights)
    
    for (idx in seq_along(ts)) {
      t <- ts[idx]; w <- pr[idx]
      if (w == 0) next
      
      if (t == 0L) {
        B[j, j] <- B[j, j] + (1 - theta$p_wgd) * (2 * w)
      } else {
        for (Np in c(N + t, N - t)) {
          if (Np >= grid$N_min && Np <= grid$N_max) {
            i <- Np - grid$N_min + 1L
            B[i, j] <- B[i, j] + (1 - theta$p_wgd) * w
          }
        }
      }
    }
    
    Nw <- 2L * N
    if (Nw >= grid$N_min && Nw <= grid$N_max) {
      iw <- Nw - grid$N_min + 1L
      B[iw, j] <- B[iw, j] + theta$p_wgd * 1.0
    }
  }
  
  return(list(B = B, lambda = lambda_val, grid = N_GRID))
}

# ============================================================
# 2. ODE DEFINITION (For deSolve fallback)
# ============================================================

ploidy_derivs <- function(t, u, parms) {
  flux <- parms$lambda * u
  inflow <- parms$B %*% flux
  outflow <- flux
  list(as.numeric(inflow - outflow))
}

# ============================================================
# 3. UNIFIED INTERFACE (Optimized)
# ============================================================

#' Run Ploidy Evolution Simulation
#' @param cached_system Optional. A list containing pre-calculated {Q, grid, lambda} to skip matrix building.
run_ploidy_model <- function(O, theta, 
                             start_N = 44L, 
                             start_u = NULL,
                             grid = list(N_min = 22L, N_max = 200L, N_unit = 22L),
                             stop_at = list(time = 1000, size = NULL),
                             method = c("expm", "desolve"),
                             cached_system = NULL) {
  
  method <- match.arg(method)
  
  # --- 1. Physics Setup ---
  if (!is.null(cached_system)) {
    # Fast path: Use pre-calculated system
    sys <- cached_system # Expects list(Q=..., grid=..., lambda=...)
    N_GRID <- sys$grid
  } else {
    # Slow path: Build from scratch
    sys <- build_transition_matrix(O, theta, grid)
    N_GRID <- sys$grid
    # Calculate Q on the fly if needed for expm
    if (method == "expm") {
      I_mat <- diag(nrow(sys$B))
      sys$Q <- sys$lambda * (sys$B - I_mat)
    }
  }
  
  # --- 2. Initial State Setup ---
  if (!is.null(start_u)) {
    if(length(start_u) != length(N_GRID)) stop("start_u length must match grid size.")
    u0 <- start_u
  } else {
    u0 <- rep(0, length(N_GRID)); names(u0) <- N_GRID
    if(as.character(start_N) %in% names(u0)) {
      u0[as.character(start_N)] <- 1
    } else {
      stop("start_N is not within the generated grid.")
    }
  }
  
  # BUG FIX: Ensure u0 has names even if start_u lost them
  if (is.null(names(u0))) names(u0) <- as.character(N_GRID)
  
  U_initial <- sum(u0)
  
  # --- 3. Execution ---
  sol <- NULL
  
  if (method == "desolve") {
    # Note: desolve needs 'B' and 'lambda', not Q. 
    # If cached_system only has Q, this path might fail unless we cache B too.
    # Assuming standard usage here.
    
    root_func <- if (!is.null(stop_at$size)) function(t, u, pars) sum(u) - stop_at$size else NULL
    t_max <- if(!is.null(stop_at$size)) 1e5 else stop_at$time
    
    sol <- ode(
      y = u0, 
      times = seq(0, t_max, length.out = 101),
      func = ploidy_derivs,
      parms = sys, 
      method = "lsoda",
      rootfunc = root_func
    )
    
    u_final <- sol[nrow(sol), -1]
    t_elapsed <- sol[nrow(sol), "time"]
    
  } else {
    # --- Matrix Exponential Implementation ---
    Q <- sys$Q
    t_target <- stop_at$time
    
    if (!is.null(stop_at$size)) {
      # Probe Logic for Size-based stopping
      probe_t <- 1.0
      Probe_Prop <- expm::expm(Q * probe_t)
      u_probe <- as.numeric(Probe_Prop %*% u0)
      size_probe <- sum(u_probe)
      
      # Calculate effective rate
      # Protect against log(0) or negative growth
      if (size_probe <= 0 || U_initial <= 0) {
        r_eff <- 0
      } else {
        r_eff <- (log(size_probe) - log(U_initial)) / probe_t
      }
      
      if (r_eff <= 1e-9) {
        # Fallback if no growth: run to max time or 0
        t_target <- 0
      } else {
        t_target <- (log(stop_at$size) - log(U_initial)) / r_eff
        if (t_target < 0) t_target <- 0
      }
    }
    
    Propagator <- expm::expm(Q * t_target)
    
    u_final_vec <- as.numeric(Propagator %*% u0)
    names(u_final_vec) <- names(u0)
    u_final <- u_final_vec
    t_elapsed <- t_target
    
    # Synthetic history
    sol <- matrix(c(0, u0, t_elapsed, u_final), nrow = 2, byrow = TRUE)
    colnames(sol) <- c("time", names(u0))
  }
  
  # --- 4. Post-processing ---
  U_final <- sum(u_final)
  r_implied <- if(t_elapsed > 0 && U_initial > 0) (log(U_final) - log(U_initial)) / t_elapsed else 0
  
  return(list(
    time_elapsed = t_elapsed,
    final_size = U_final,
    growth_rate = r_implied,
    u_vec = as.numeric(u_final), 
    distribution = data.frame(
      N = N_GRID, 
      count = as.numeric(u_final), 
      fraction = as.numeric(u_final) / (U_final + 1e-16)
    ),
    full_history = sol
  ))
}

