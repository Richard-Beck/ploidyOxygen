load_db_vars <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("Error: File '%s' does not exist.", filepath))
  }
  
  lines <- tryCatch(readLines(filepath, warn = FALSE), 
                    error = function(e) stop("Error: Failed to read file '", filepath, "'. ", e$message))
  
  if (length(lines) == 0) {
    stop(sprintf("Error: File '%s' is empty.", filepath))
  }
  
  if (any(!grepl("^[A-Za-z_][A-Za-z0-9_]*=.+$", lines))) {
    stop("Error: All lines must be in the format KEY=value, with valid variable names.")
  }
  
  kv_pairs <- strsplit(lines, "=", fixed = TRUE)
  keys <- vapply(kv_pairs, `[`, "", 1)
  values <- vapply(kv_pairs, `[`, "", 2)
  vars <- setNames(values, keys)
  
  required_keys <- c("HOST", "DBNAME", "USER", "PASSWORD")
  missing_keys <- setdiff(required_keys, names(vars))
  if (length(missing_keys) > 0) {
    stop("Error: Missing required keys: ", paste(missing_keys, collapse = ", "))
  }
  
  return(vars)
}

get_karyotyping <- function(cloneIds){
  require(DBI)
  require(RMariaDB)
  dbvars <- load_db_vars("db_creds.txt")
  db <- dbConnect(
    MariaDB(),
    host = dbvars["HOST"],
    user = dbvars["USER"],
    password = dbvars["PASSWORD"],
    dbname = dbvars["DBNAME"]
  )
  
  origin_clause <- paste0("'", cloneIds, "'", collapse = ", ")
  q <- paste0("SELECT * FROM Perspective WHERE origin IN (", origin_clause, ") AND whichPerspective='GenomePerspective'")
  
  rs <- dbSendQuery(db, q)
  res <- dbFetch(rs)
  dbClearResult(rs) 
  dbDisconnect(db)
  
  kvecs <- lapply(res$profile, function(p) {
    raw_vec <- p
    readBin(raw_vec, what = "double", n = length(raw_vec) / 8, endian = "big")
  })
  
  df <- tibble::tibble(
    id = res$origin,
    karyotype = kvecs
  )
  df
}

recover_lineage <- function(row, data) {
  current_id <- row$last
  lineage <- c(current_id)
  
  while (!is.na(current_id) && current_id != row$first) {
    parent <- data$passaged_from_id1[which(data$id == current_id)]
    if (length(parent) == 0 || is.na(parent)) {
      stop(paste("Could not trace from", row$last, "to", row$first))
    }
    lineage <- c(parent, lineage)
    current_id <- parent
  }
  
  if (current_id != row$first) {
    stop(paste("Could not find", row$first, "starting from", row$last))
  }
  
  seeding_lineage <- lineage[data$event[match(lineage, data$id)] == "seeding"]
  
  filtered_x <- do.call(rbind, lapply(seq_along(seeding_lineage), function(i) {
    dfi <- data[data$id %in% seeding_lineage[i] | data$passaged_from_id1 %in% seeding_lineage[i], ]
    dfi$adjPass <- stringr::str_pad(i, width = 2)
    dfi$date <- as.Date(dfi$date)
    dfi$num_date <- as.numeric(as.Date(dfi$date))
    dfi$num_date <- dfi$num_date - min(dfi$num_date)
    dfi$intercept <- NaN
    dfi$g <- NaN
    if (nrow(dfi) < 2) return(dfi)
    if(sum(!is.na(dfi$correctedCount))<2) return(dfi)
    fit <- lm(log(pmax(1, dfi$correctedCount)) ~ dfi$num_date)
    dfi$intercept <- exp(coef(fit)[1])
    dfi$g <- coef(fit)[2]
    dfi
  }))
  
  filtered_x$label_value    <- row$label
  filtered_x$sublabel_value <- row$label2
  return(filtered_x)
}


fill_lineage_gaps <- function(ploidy_substr,df,x){
  g <- x[x$cellLine=="SUM-159",c("id","passaged_from_id1")]
  
  
  ## perhaps split 2N 4N.
  karyotyped_lineages <- lapply(cloneIds[grepl(ploidy_substr,cloneIds)],function(id){
    ids <- c()
    while(!is.na(id)){
      ids <- c(id,ids)
      id <- g$passaged_from_id1[g$id==id]
    }
    ids
  })
  
  init_ids <- sapply(filters,'[[',"first")
  init_ids <- init_ids[grepl(ploidy_substr,init_ids)]
  
  imaged_lineages <- lapply(init_ids,function(id){
    ids <- c()
    while(!is.na(id)){
      ids <- c(id,ids)
      id <- g$passaged_from_id1[g$id==id]
    }
    ids
  })
  
  all_lineages <- c(imaged_lineages,karyotyped_lineages)
  n_lineages <- length(all_lineages)
  
  uids <- table(unlist(all_lineages))
  n_remove <- length(uids[uids==n_lineages])-1
  
  all_lineages <- lapply(all_lineages, function(i){
    i[-c(1:n_remove)]
  })
  
  x <- x[x$id%in%unlist(all_lineages),]
  x <- x[!x$id%in%df$id,]
  x$label_value <- gsub("_","",ploidy_substr)
  x
}