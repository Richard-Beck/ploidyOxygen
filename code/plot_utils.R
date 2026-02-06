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