experiment_lineage_plot <- function(lineage_df_subset){
  library(dplyr)
  library(ape)
  library(ggplot2)
  library(ggtree)
  
  # ---- 1) edge-list -> ape::phylo (handles labeled internal nodes too) ----
  edge_df_to_phylo <- function(df, id_col="passage_id", parent_col="passage_from", edge_length=NULL) {
    df <- df %>%
      transmute(id = .data[[id_col]], parent = .data[[parent_col]]) %>%
      distinct()
    
    stopifnot(!anyDuplicated(df$id))
    stopifnot(!anyDuplicated(df$id[!is.na(df$parent)]) || TRUE) # parent can repeat (branching)
    stopifnot(!any(is.na(df$id)))
    
    # root: parent is NA or parent not present as an id
    roots <- df$id[is.na(df$parent) | !(df$parent %in% df$id)]
    if (length(roots) != 1) stop("Expected exactly one root; found: ", paste(roots, collapse=", "))
    
    nodes <- df$id
    tips  <- setdiff(nodes, df$parent[!is.na(df$parent)])
    ints  <- setdiff(nodes, tips)
    
    # ape phylo indexing: tips = 1..ntip, internal = (ntip+1)..(ntip+nnode)
    ntip <- length(tips)
    nnode <- length(ints)
    idx <- c(setNames(seq_len(ntip), tips),
             setNames(ntip + seq_len(nnode), ints))
    
    edges <- df[!is.na(df$parent)& df$parent%in%nodes,] %>%
      transmute(parent_i = unname(idx[parent]), child_i = unname(idx[id]))
    
    phy <- list(
      edge = as.matrix(edges),
      Nnode = nnode,
      tip.label = tips
    )
    class(phy) <- "phylo"
    
    # internal node labels (so you can join metadata by id everywhere)
    phy$node.label <- ints
    
    # optional edge lengths
    if (!is.null(edge_length)) {
      el <- edge_length(df)               # should return a numeric vector length nrow(edges)
      phy$edge.length <- as.numeric(el)
    } else {
      phy$edge.length <- rep(1, nrow(edges))
    }
    
    phy <- reorder.phylo(phy)
    phy
  }
  
  phy <- edge_df_to_phylo(lineage_df_subset)
  
  
  meta <- lineage_df_subset %>%
    transmute(
      label = as.character(passage_id),
      g = ifelse(is.finite(g), g, NA_real_),
      karyotyped = as.logical(karyotyped),
      has_flow = as.logical(has_flow)
    )
  
  p <- ggtree(phy, layout="circular", size=0.3)
  
  # explicit join to ggtree data
  p$data <- left_join(p$data, meta, by="label")
  
  p <- p +
    # --- karyotyped underneath ---
    geom_point2(
      data = ~ dplyr::filter(.x, karyotyped),
      aes(shape = "karyotyped"), color = "red", size = 3.2
    ) +
    geom_point2(
      data = ~ dplyr::filter(.x, has_flow),
      aes(shape = "has_flow"), color = "green", size = 3.2
    ) +
    # --- growth for all nodes ---
    geom_point2(
      aes(color = g),
      shape = 16, size = 1.4
    ) +
    scale_color_viridis_c("growth\nrate",na.value = "grey70",trans="log") +
    theme_void()+
    scale_shape_discrete("")+
    ggtitle(paste("Hypoxia experiment data:",unique(lineage_df_subset$label_value)))
  return(p)
}

library(ggplot2)
library(gridExtra)

#' Generate Enhanced Ploidy Evolution Dashboard
#' 
#' @param O Oxygen concentration (0-21 typically).
#' @param theta Named list of parameters.
#' @param start_N Starting ploidy.
#' @param times Vector of specific times to snapshot distributions.
#' @param grid List containing N_min, N_max, N_unit.
#' 
#' @return A grid of 5 ggplot objects arranged in a custom layout.
plot_ploidy_dashboard <- function(O, theta, start_N, times, 
                                  grid = list(N_min = 22L, N_max = 200L, N_unit = 22L)) {
  
  require(ggplot2)
  require(gridExtra)
  require(deSolve)
  
  # =========================================================
  # 1. PREPARE STATIC PHYSICS PLOTS
  # =========================================================
  
  # --- Plot 1: Viability (sN) vs Ploidy ---
  N_seq <- grid$N_min:grid$N_max
  viability <- theta$smax * exp(-theta$beta * ((2 * grid$N_unit) / N_seq)^theta$n_exp)
  
  df_viability <- data.frame(N = N_seq, Survival = viability)
  
  p1 <- ggplot(df_viability, aes(x = N, y = Survival)) +
    geom_line(color = "#2c3e50", size = 1.2) +
    geom_vline(xintercept = start_N, linetype = "dashed", color = "red", alpha = 0.5) +
    labs(title = "1. Viability vs Ploidy", y = "Viability (s)", x = "Ploidy (N)") +
    theme_minimal(base_size = 10)
  
  # --- Plot 2: Missegregation Rate vs Oxygen ---
  O_seq <- seq(0, max(21, O * 1.5), length.out = 100)
  # Func: p_misseg * (1 - (O / (O + k_o_mis)))
  pmis_seq <- theta$p_misseg * (1 - (O_seq / (O_seq + theta$k_o_mis)))
  
  df_pmis <- data.frame(Oxygen = O_seq, Misseg = pmis_seq)
  
  p2 <- ggplot(df_pmis, aes(x = Oxygen, y = Misseg)) +
    geom_line(color = "#e67e22", size = 1.2) +
    geom_point(aes(x = O, y = theta$p_misseg * (1 - (O / (O + theta$k_o_mis)))), 
               color = "red", size = 3) +
    geom_vline(xintercept = O, linetype = "dotted", color = "red") +
    labs(title = "2. Missegregation vs O2", y = "Prob(Misseg)", x = "Oxygen (%)") +
    theme_minimal(base_size = 10)
  
  # --- Plot 3: Cycle Rate (Lambda) vs Oxygen (NEW) ---
  # Func: lam_min + (lam_max - lam_min) * (O / (O + k_o))
  lambda_seq <- theta$lam_min + (theta$lam_max - theta$theta_min) * (O_seq / (O_seq + theta$k_o))
  # Note: Handle case where theta_min might be named differently or just subtract
  lambda_diff <- theta$lam_max - theta$lam_min
  lambda_seq <- theta$lam_min + lambda_diff * (O_seq / (O_seq + theta$k_o))
  
  current_lambda <- theta$lam_min + lambda_diff * (O / (O + theta$k_o))
  
  df_lambda <- data.frame(Oxygen = O_seq, Lambda = lambda_seq)
  
  p3 <- ggplot(df_lambda, aes(x = Oxygen, y = Lambda)) +
    geom_line(color = "#8e44ad", size = 1.2) +
    geom_point(aes(x = O, y = current_lambda), color = "red", size = 3) +
    geom_vline(xintercept = O, linetype = "dotted", color = "red") +
    labs(title = "3. Cycle Rate (Lambda) vs O2", y = "Rate (1/hr)", x = "Oxygen (%)") +
    theme_minimal(base_size = 10)
  
  # =========================================================
  # 2. RUN SIMULATION
  # =========================================================
  
  sys <- build_transition_matrix(O, theta, grid)
  
  u0 <- rep(0, length(sys$grid))
  names(u0) <- sys$grid
  if(as.character(start_N) %in% names(u0)) {
    u0[as.character(start_N)] <- 1
  } else {
    stop("start_N outside grid range")
  }
  
  # Ensure simulation runs long enough and hits specific timepoints
  t_max <- max(times)
  if(t_max == 0) t_max <- 10 # Fallback
  t_seq <- sort(unique(c(seq(0, t_max, length.out = 200), times)))
  
  sol <- deSolve::ode(
    y = u0, 
    times = t_seq, 
    func = ploidy_derivs, 
    parms = sys,
    method = "lsoda"
  )
  
  df_sol <- as.data.frame(sol)
  
  # =========================================================
  # 3. PREPARE DYNAMIC PLOTS
  # =========================================================
  
  # --- Plot 4: Net Growth Rate vs Time ---
  pop_size <- rowSums(df_sol[, -1])
  time_vec <- df_sol$time
  
  # Central Difference for smooth derivative
  growth_rates <- numeric(length(time_vec))
  n_t <- length(time_vec)
  
  for(i in 1:n_t) {
    if(i == 1) {
      dt <- time_vec[2] - time_vec[1]
      growth_rates[i] <- (log(pop_size[2]) - log(pop_size[1])) / dt
    } else if (i == n_t) {
      dt <- time_vec[n_t] - time_vec[n_t-1]
      growth_rates[i] <- (log(pop_size[n_t]) - log(pop_size[n_t-1])) / dt
    } else {
      dt <- time_vec[i+1] - time_vec[i-1]
      growth_rates[i] <- (log(pop_size[i+1]) - log(pop_size[i-1])) / dt
    }
  }
  
  # Clean up NaNs if population hits 0
  growth_rates[!is.finite(growth_rates)] <- 0
  
  df_growth <- data.frame(Time = time_vec, Rate = growth_rates)
  
  p4 <- ggplot(df_growth, aes(x = Time, y = Rate)) +
    geom_line(color = "#27ae60", size = 1) +
    geom_point(data = df_growth[df_growth$Time %in% times, ], color = "black", size = 2) +
    labs(title = "4. Pop. Growth Rate (History)", y = "d(ln N)/dt", x = "Time") +
    theme_minimal(base_size = 10)
  
  # --- Plot 5: Distributions (LARGE COLUMN) ---
  df_dist_list <- list()
  for(t_req in times) {
    # Find closest index
    idx <- which.min(abs(df_sol$time - t_req))
    row_vals <- df_sol[idx, -1]
    
    total_pop <- sum(row_vals)
    freq <- if(total_pop > 0) as.numeric(row_vals) / total_pop else as.numeric(row_vals)
    
    tmp <- data.frame(
      Time = paste0("t=", round(df_sol$time[idx], 0)),
      N = sys$grid,
      Frequency = freq
    )
    df_dist_list[[length(df_dist_list)+1]] <- tmp
  }
  
  df_dists <- do.call(rbind, df_dist_list)
  df_dists$Time <- factor(df_dists$Time, levels = unique(df_dists$Time))
  
  p5 <- ggplot(df_dists, aes(x = N, y = Frequency, fill = Time)) +
    geom_area(alpha = 0.4, position = "identity") +
    geom_line(aes(color = Time), size = 1) +
    labs(title = "5. Ploidy Distribution Evolution", x = "Ploidy (N)", y = "Frequency") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top", 
          legend.text = element_text(size=12))+
    scale_y_sqrt()
  
  # =========================================================
  # 4. ASSEMBLE DASHBOARD
  # =========================================================
  
  # Layout: Left column (plots 1,2,3,4), Right column (plot 5 full height)
  # Matrix construction: 
  # 1  5
  # 2  5
  # 3  5
  # 4  5
  
  layout_matrix <- rbind(
    c(1, 2, 5),
    c(3, 4 , 5)
  )
  
  grid.arrange(p1, p2, p3, p4, p5, layout_matrix = layout_matrix,
               widths = c(1, 1, 2))
}