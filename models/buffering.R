suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ============================================================
# INTERFACE: The Model Core
# ============================================================

#' Run Ploidy Evolution Simulation (Refactored)
#' @param O Oxygen level (scalar)
#' @param theta Named list: p_misseg, p_wgd, beta, n_exp, smax, k_o, lam_max
#' @param start_N Starting ploidy (scalar, used if start_u is NULL)
#' @param start_u Starting distribution vector (optional, for chaining)
#' @param grid List with N_min, N_max, N_unit
#' @param stop_at List with 'time' (max duration) and optional 'size' (target U)
run_ploidy_model <- function(O, theta, 
                             start_N = 44L, 
                             start_u = NULL,
                             grid = list(N_min = 22L, N_max = 200L, N_unit = 22L),
                             stop_at = list(time = 1000, size = NULL)) {
  
  N_GRID <- grid$N_min:grid$N_max
  
  # 1. Fitness: Constant across N, Michaelis-Menten function of O
  # lambda = lam_max * O / (O + k_o)
  lambda_val <- theta$lam_max * (O / (O + theta$k_o))
  lambda_vec <- rep(lambda_val, length(N_GRID))
  
  # 2. Transition Matrix B Build
  # Internal helper for missegregation survival weights
  .get_delta_weights <- function(N, p, b, n, sm, k_unit) {
    if (p <= 0 || N <= 0) return(c("0" = 1))
    sd <- sqrt(N * p)
    if (sd == 0) return(c("0" = 1))
    
    # Anchored Dosage Survival: s = smax * exp(-beta * (2k/N)^n)
    sN <- sm * exp(-b * ((2 * k_unit) / N)^n)
    
    z <- 9.0 # Bounds for numerical integration (eps_tail ~ 1e-20)
    T_range <- min(N, max(0L, ceiling(z * sd)))
    ts <- (-T_range):T_range
    out <- numeric(length(ts))
    
    for (idx in seq_along(ts)) {
      t <- ts[idx]
      ks <- seq.int(abs(t), N, by = 2)
      if (!length(ks)) next
      # Probability of k misseg events * prob of delta t * survival of k events
      out[idx] <- sum(dbinom(ks, N, p) * dbinom((ks + t)/2, ks, 0.5) * (sN^ks))
    }
    names(out) <- ts
    out
  }
  
  B <- matrix(0, nrow = length(N_GRID), ncol = length(N_GRID), dimnames = list(N_GRID, N_GRID))
  
  for (j in seq_along(N_GRID)) {
    N <- N_GRID[j]
    weights <- .get_delta_weights(N, theta$p_misseg, theta$beta, theta$n_exp, theta$smax, grid$N_unit)
    ts <- as.integer(names(weights)); pr <- as.numeric(weights)
    
    for (idx in seq_along(ts)) {
      t <- ts[idx]; w <- pr[idx]
      if (w == 0) next
      
      if (t == 0L) {
        # Normal division (2 daughters) minus WGD probability
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
    
    # Whole Genome Doubling: 1 daughter at 2N
    Nw <- 2L * N
    if (Nw >= grid$N_min && Nw <= grid$N_max) {
      iw <- Nw - grid$N_min + 1L
      B[iw, j] <- B[iw, j] + theta$p_wgd * 1.0
    }
  }
  
  # 3. Setup Initial State
  if (!is.null(start_u)) {
    if(length(start_u) != length(N_GRID)) stop("start_u length must match grid size.")
    u0 <- start_u
  } else {
    u0 <- rep(0, length(N_GRID)); names(u0) <- N_GRID
    u0[as.character(start_N)] <- 1
  }
  U_initial <- sum(u0)
  
  # 4. Termination Logic (Root finding)
  root_func <- NULL
  if (!is.null(stop_at$size)) {
    root_func <- function(t, u, pars) sum(u) - stop_at$size
  }
  
  # 5. ODE Solver
  sol <- ode(
    y = u0, 
    times = seq(0, stop_at$time, length.out = 101),
    func = function(t, u, pars) list(as.numeric(B %*% (lambda_vec * u) - (lambda_vec * u))),
    parms = NULL, 
    method = "lsoda",
    rootfunc = root_func
  )
  
  # 6. Post-processing
  u_final <- sol[nrow(sol), -1]
  U_final <- sum(u_final)
  t_elapsed <- max(sol[, "time"])
  
  # Net Growth Rate r = [ln(U_final) - ln(U_initial)] / t
  r_implied <- if(t_elapsed > 0 && U_initial > 0) (log(U_final) - log(U_initial)) / t_elapsed else 0
  
  return(list(
    time_elapsed = t_elapsed,
    final_size = U_final,
    growth_rate = r_implied,
    u_vec = as.numeric(u_final), # Raw vector for easy chaining
    distribution = data.frame(
      N = N_GRID, 
      count = as.numeric(u_final), 
      fraction = as.numeric(u_final) / (U_final + 1e-16)
    ),
    full_history = sol
  ))
}
