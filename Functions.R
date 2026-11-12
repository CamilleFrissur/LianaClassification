## Import libraries #### 
library(terra)      
library(data.table)  
library(miscTools)   
library(psych) 

# Preprocessing ####
extract_wavelength <- function(band_names) {
  
  wl <- suppressWarnings(as.numeric(band_names))
  need_regex <- is.na(wl)
  wl[need_regex] <- as.numeric(gsub(".*\\(([0-9.]+)\\).*", "\\1", band_names[need_regex]))
  
  return(wl)
}

clean_extreme_values <- function(area, valid_range = c(0, 10000)) {
  
  area_matrix <- as.matrix(area)
  
  n_below <- sum(area_matrix < valid_range[1], na.rm = TRUE)
  n_above <- sum(area_matrix > valid_range[2], na.rm = TRUE)
  n_total_out <- n_below + n_above
  n_total_values <- sum(!is.na(area_matrix))
  
  area_clamped <- clamp(area, lower = valid_range[1], upper = valid_range[2], values = TRUE)
  
  list(
    area_clean = area_clamped,
    n_below = n_below,
    n_above = n_above,
    n_total_out = n_total_out,
    perc_n_out = n_total_out/n_total_values
  )
}

clean_extreme_values_toNA <- function(area, valid_range = c(0, 10000)) {
  n_below <- global(area < valid_range[1], "sum", na.rm = TRUE)$sum
  n_above <- global(area > valid_range[2], "sum", na.rm = TRUE)$sum
  
  cat("Number of pixel out of range:", (sum(n_below)+sum(n_above)), "\n")
  
  area[area < valid_range[1]] <- NA
  area[area > valid_range[2]] <- NA
  return(area)
}

filter_good_bands <- function(area, max_na_pct_band = 0.95) {
  wl <- extract_wavelength(names(area))
  #good_bands_idx <- which(wl < 414 | wl > 447 & wl < 1300 | wl > 1500 & wl < 1800 | wl > 2000)
  good_bands_idx <- which((wl < 414) | (wl > 447 & wl < 1329.7) | (wl > 1497.9 & wl < 1774.8) | (wl > 1981 & wl < 2361.0))
  area <- area[[good_bands_idx]]
  
  area_matrix <- as.matrix(area)
  na_pct_per_band <- colSums(is.na(area_matrix)) / nrow(area_matrix)
  good_bands <- na_pct_per_band < max_na_pct_band
  
  cat("Bands kept :", sum(good_bands), "/", length(good_bands), "\n")
  cat("Bands removed: (>=", max_na_pct_band*100, "% NA) :", sum(!good_bands), "\n")
  
  area[[good_bands]]
}

data_prep_per_pixel <- function(area_filtered, label_filtered, illu_crop_area, max_na_pct = 0.05){
  area_matrix <- as.matrix(area_filtered)
  label_matrix <- as.matrix (label_filtered)[,1]
  
  valid_idx <- which(!is.na(label_matrix) & rowSums(is.na(area_matrix)) == 0)
  
  area_valid <- area_matrix[valid_idx, , drop = FALSE]
  labels_valid <- label_matrix[valid_idx]
  
  illu_matrix <- as.matrix(illu_crop_area)[, 1]
  illu_valid <- illu_matrix[valid_idx]
  
  output_list <- list(
    matrix = area_valid,
    labels = labels_valid,
    illu = illu_valid,
    valid_idx = valid_idx)
  
  return(output_list)
} 

data_prep_per_poly <- function(area_filtered, area_vector, label_name){
  
  # Mean spectrum per polygon
  poly_spectra <- extract(area_filtered, area_vector, fun = mean, na.rm = TRUE, ID = TRUE)
  poly_labels <- values(area_vector)[[label_name]][poly_spectra$ID]
  
  # Create the matrix
  poly_matrix <- as.matrix(poly_spectra[, -1])
  
  na_per_band_poly <- colSums(is.na(poly_matrix))
  good_bands_poly <- na_per_band_poly < nrow(poly_matrix)
  
  poly_valid <- poly_matrix[, good_bands_poly]
  
  # Remove polygons full of NA
  valid_poly_idx <- rowSums(is.na(poly_valid)) == 0
  
  poly_valid <- poly_valid[valid_poly_idx, ]
  label_poly_valid <- poly_labels[valid_poly_idx]
  dim(poly_valid)
  length(label_poly_valid)
  
  output_list <- list(
    matrix = poly_valid,
    labels = label_poly_valid,
    goodbands = good_bands_poly
  )
  
  return(output_list)
}

# Filters ####
filter_by_illu <- function(illu_raster, area, label, threshold){
  illu_crop_area <- crop(illu_raster, area)
  n_before <- global(!is.na(area[[1]]), "sum")$sum
  
  illu_mask <- illu_crop_area > threshold
  
  area_filtered <- mask(area, illu_mask, maskvalues = c(FALSE, NA))
  label_filtered <- mask(label, illu_mask, maskvalues = c(FALSE, NA))
  
  n_after <- global(!is.na(area_filtered[[1]]), "sum")$sum
  
  per_n_out <- n_after/n_before
  
  output_list <- list(
    area_crop = illu_crop_area,
    area_filt = area_filtered,
    label_filt = label_filtered,
    threshold_qual = per_n_out)
  
  return(output_list)
}

filter_by_vnir <- function(area, label, vnir_range = c(400, 1100)){
  
  # area_matrix <- as.matrix(area)
  # na_per_band <- colSums(is.na(area_matrix))
  # good_bands <- na_per_band < nrow(area_matrix)
  # 
  # wl <- extract_wavelength(names(area))[good_bands]
  
  wl <- extract_wavelength(names(area))
  
  vnir_idx <- which(wl >= vnir_range[1] & wl <= vnir_range[2])
  vnir_stack <- area[[vnir_idx]]
  
  vnir_mean_px <- mean(vnir_stack, na.rm = TRUE)
  vnir_mean_values <- as.vector(vnir_mean_px)
  vnir_mean_values <- vnir_mean_values[!is.na(vnir_mean_values)]
  
  threshold_value <- median(vnir_mean_values, na.rm = TRUE)
  # threshold_value <- quantile(vnir_mean_values, probs = 0.5, na.rm = TRUE)
  
  vnir_mask <- vnir_mean_px > threshold_value
  
  n_before <- global(!is.na(area[[1]]), "sum")$sum
  
  area_filtered <- mask(area, vnir_mask, maskvalues = c(FALSE, NA))
  label_filtered <- mask(label, vnir_mask, maskvalues = c(FALSE, NA))
  
  n_after <- global(!is.na(area_filtered[[1]]), "sum")$sum
  
  output_list <- list(
    threshold = threshold_value,
    area_filt = area_filtered,
    label_filt = label_filtered,
    threshold_qual = n_after / n_before
  )
  
  return(output_list)
}

filter_by_nir <- function(area, label, nir_range = c(780, 900)){
  
  # area_matrix <- as.matrix(area)
  # na_per_band <- colSums(is.na(area_matrix))
  # good_bands <- na_per_band < nrow(area_matrix)
  # 
  # wl <- extract_wavelength(names(area))[good_bands]
  
  wl <- extract_wavelength(names(area))
  
  nir_idx <- which(wl >= nir_range[1] & wl <= nir_range[2])
  
  nir_stack <- area[[nir_idx]]
  nir_mean_px <- mean(nir_stack, na.rm = TRUE)
  nir_mean_values <- as.vector(nir_mean_px)
  nir_mean_values <- nir_mean_values[!is.na(nir_mean_values)]
  
  threshold_value <- 0.20 * mean(nir_mean_values, na.rm = TRUE)
  # threshold_value <- quantile(nir_mean_values, probs = 0.2, na.rm = TRUE)
  
  nir_mask <- nir_mean_px >= threshold_value
  
  n_before <- global(!is.na(area[[1]]), "sum")$sum
  
  area_filtered <- mask(area, nir_mask, maskvalues = c(FALSE, NA))
  label_filtered <- mask(label, nir_mask, maskvalues = c(FALSE, NA))
  
  n_after <- global(!is.na(area_filtered[[1]]), "sum")$sum
  
  output_list <- list(
    threshold_value = threshold_value,
    area_filt = area_filtered,
    label_filt = label_filtered,
    threshold_qual = n_after / n_before
  )
  
  return(output_list)
}

# Metrics ####
run_chi2 <- function(pixels1, pixels2, name1, name2) {
  tab <- table(pixels1, pixels2)
  test <- chisq.test(tab)
  cat(sprintf("\n--- %s vs %s ---\n", name1, name2))
  cat(sprintf("Chi² = %.4f, df = %d, p-value = %.4f\n",
              test$statistic, test$parameter, test$p.value))
  if (test$p.value < 0.05) {
    cat("→ The filters target similar pixels (p < 0.05).\n")
  } else {
    cat("→ The filters target different pixels (p >= 0.05).\n")
  }
}

calculate_cohen_kappa <- function(r1, r2) {
  if (!identical(ext(r1), ext(r2)) || !identical(res(r1), res(r2))) {
    warning("Rasters do not have the same extent or resolution. Use of intersection.")
    r1 <- mask(r1, r2)  
    r2 <- mask(r2, r1)
  }
  
  v1 <- as.logical(as.vector(r1, mode = "integer"))
  v2 <- as.logical(as.vector(r2, mode = "integer"))
  
  # 3. Compter les combinaisons avec table() (évite les débordements)
  tab <- table(v1, v2, useNA = "no")  # useNA="no" pour ignorer les NA
  if (length(tab) < 4) {
    stop("Not enough combination to calcultate Kappa (Need of TRUE/TRUE, TRUE/FALSE, FALSE/TRUE, FALSE/FALSE).")
  }
  
  a <- tab[4]  # TRUE/TRUE
  b <- tab[3]  # TRUE/FALSE
  c <- tab[2]  # FALSE/TRUE
  d <- tab[1]  # FALSE/FALSE
  N <- a + b + c + d
  
  p_o <- (a + d) / N
  p1_true <- (a + b) / N
  p2_true <- (a + c) / N
  p1_false <- (c + d) / N
  p2_false <- (b + d) / N
  p_e <- (p1_true * p2_true) + (p1_false * p2_false)
  
  kappa <- (p_o - p_e) / (1 - p_e)
  return(kappa)
}

# Spectrum ####

sensor_to_keep <- function(matrix_valid, keep_sensor_in_overlap = c("swir", "vnir")) {
  
  keep_sensor_in_overlap <- match.arg(keep_sensor_in_overlap)
  wl <- extract_wavelength(colnames(matrix_valid))
  
  break_idx <- which(diff(wl) < 0)[1]
  
  if (is.na(break_idx))  {
    cat("No sensor overlap detected\n")
    return(matrix_valid)
  }
  vnir_idx <- 1:break_idx
  swir_idx <- (break_idx + 1):length(wl)
  
  overlap_min <- wl[swir_idx[1]]
  overlap_max <- wl[break_idx]
  
  if (keep_sensor_in_overlap == "swir") {
    drop_idx <- vnir_idx[wl[vnir_idx] >= overlap_min]   
  } else {
    drop_idx <- swir_idx[wl[swir_idx] <= overlap_max]    
  }
  
  keep_overlap <- !(seq_along(wl) %in% drop_idx)
  
  cat("Overlap min:", overlap_min, "\nOverlap max:", overlap_max, "\n")
  matrix_valid <- matrix_valid[, which(keep_overlap), drop = FALSE]}

mean_spectrum_by_group <- function(matrix_valid, labels_valid, dataset_name) {
  groups <- unique(labels_valid)
  
  wavelengths <- extract_wavelength(colnames(matrix_valid))
  
  spectra_list <- lapply(groups, function(g) {
    group_matrix <- matrix_valid[labels_valid == g, , drop = FALSE] 
    mean_spec <- colMeans(group_matrix, na.rm = TRUE)
    sd_spec <- apply(group_matrix, 2, sd, na.rm = TRUE) 
    
    data.frame(
      wavelength = wavelengths,
      mean = mean_spec,
      sd= sd_spec,
      group = paste0(dataset_name, "_", g)
    )
  })
  
  do.call(rbind, spectra_list)
}

add_segments <- function(dt, gap_nm = 8) {
  dt <- as.data.table(dt)
  dt[, segment := cumsum(c(TRUE, diff(wavelength) <= 0 | diff(wavelength) > gap_nm)), by = group]
  dt
}

med_spectrum_by_group <- function(matrix_valid, labels_valid, dataset_name) {
  groups <- unique(labels_valid)
  
  wavelengths <- extract_wavelength(colnames(matrix_valid))
  
  spectra_list <- lapply(groups, function(g) {
    group_matrix <- matrix_valid[labels_valid == g, , drop = FALSE] 
    med_spec <- colMedians(group_matrix, na.rm = TRUE)
    sd_spec <- apply(group_matrix, 2, sd, na.rm = TRUE) 
    
    data.frame(
      wavelength = wavelengths,
      med = med_spec,
      sd= sd_spec,
      group = paste0(dataset_name, "_", g)
    )
  })
  
  do.call(rbind, spectra_list)
}

# Bhattacharya ####

extract_bhatt_dist <- function(dist_raw) {
  if (is.null(dist_raw) || is.na(dist_raw)) return(NA)
  
  if (is.numeric(dist_raw) && length(dist_raw) == 1) return(as.numeric(dist_raw))
  
  if (inherits(dist_raw, "modulus")) return(as.numeric(dist_raw))
  
  if (is.list(dist_raw) && "modulus" %in% names(dist_raw)) return(as.numeric(dist_raw$modulus))
  
  tryCatch(as.numeric(dist_raw), error = function(e) NA)
}

regularize_cov <- function(sigma, epsilon = 0) {
  diag(sigma) <- diag(sigma) + epsilon * max(diag(sigma), 1)
  return(sigma)
}

bhatt_between <- function(df, label_col, class1, class2, epsilon = 0) {
  exclude_cols <- c("poly_id", "Liana", "LianaOriginal", "Type")
  numeric_cols <- sapply(df, is.numeric)
  numeric_cols[names(numeric_cols) %in% exclude_cols] <- FALSE
  
  data <- as.matrix(df[, numeric_cols, drop = FALSE])
  labels <- df[[label_col]]
  
  d1 <- data[labels == class1, ]
  d2 <- data[labels == class2, ]
  
  var_d1 <- apply(d1, 2, var)
  var_d2 <- apply(d2, 2, var)
  good_bands <- (var_d1 > 1e-6) & (var_d2 > 1e-6)
  d1 <- d1[, good_bands, drop = FALSE]
  d2 <- d2[, good_bands, drop = FALSE]
  
  if (ncol(d1) == 0) return(NA)
  
  mu1 <- colMeans(d1)
  mu2 <- colMeans(d2)
  sigma1 <- regularize_cov(cov(d1), epsilon)
  sigma2 <- regularize_cov(cov(d2), epsilon)
  
  dist_raw <- tryCatch(bhattacharyya.dist(mu1, mu2, sigma1, sigma2), error = function(e) NA)
  dist <- extract_bhatt_dist(dist_raw)
  
  if (is.na(dist)) return(NA)
  
  list(dist = dist, coeff = exp(-dist), n1 = nrow(d1), n2 = nrow(d2), p = ncol(d1))
}

bhatt_between_random <- function(df, label_col, class1, class2, n_repeats = 100, seed = 123, epsilon = 0) {
  df_filtered <- df[df[[label_col]] %in% c(class1, class2), ]
  
  exclude_cols <- c("poly_id", "Liana", "LianaOriginal", "Type")
  numeric_cols <- sapply(df_filtered, is.numeric)
  numeric_cols[names(numeric_cols) %in% exclude_cols] <- FALSE
  
  data <- as.matrix(df_filtered[, numeric_cols, drop = FALSE])
  labels <- df_filtered[[label_col]]
  
  idx_class1 <- which(labels == class1)
  idx_class2 <- which(labels == class2)
  
  distances <- numeric(n_repeats)
  
  set.seed(seed)
  for (i in 1:n_repeats) {

    sample_idx1 <- sample(idx_class1, size = floor(length(idx_class1) * 0.5), replace = FALSE)
    sample_idx2 <- sample(idx_class2, size = floor(length(idx_class2) * 0.5), replace = FALSE)
    
    d1 <- data[sample_idx1, ]
    d2 <- data[sample_idx2, ]
    
    var_d1 <- apply(d1, 2, var)
    var_d2 <- apply(d2, 2, var)
    good_bands <- (var_d1 > 1e-6) & (var_d2 > 1e-6)
    d1 <- d1[, good_bands, drop = FALSE]
    d2 <- d2[, good_bands, drop = FALSE]
    
    if (ncol(d1) == 0) {
      distances[i] <- NA
      next
    }
    
    mu1 <- colMeans(d1)
    mu2 <- colMeans(d2)
    sigma1 <- regularize_cov(cov(d1), epsilon)
    sigma2 <- regularize_cov(cov(d2), epsilon)
    
    dist_raw <- tryCatch({
      bhattacharyya.dist(mu1, mu2, sigma1, sigma2)
    }, error = function(e) NA)
    
    distances[i] <- extract_bhatt_dist(dist_raw)
  }
  
  valid_dists <- distances[!is.na(distances)]
  if (length(valid_dists) == 0) return(NA)
  
  list(
    mean_dist = mean(valid_dists),
    sd_dist = sd(valid_dists),
    min_dist = min(valid_dists),
    max_dist = max(valid_dists),
    n_repeats = n_repeats,
    n_valid = length(valid_dists),
    all_distances = distances
  )
}

bhatt_within_random <- function(df, poly_id_col, label_col, target_class,
                                n_repeats = 100, seed = 123, epsilon = 0) {
  df_class <- df[df[[label_col]] == target_class, ]
  unique_polys <- unique(df_class[[poly_id_col]])
  unique_polys <- unique_polys[!is.na(unique_polys) & unique_polys != 0]
  
  if (length(unique_polys) < 2) {
    warning(paste("Less than two poly per class", target_class))
    return(NA)
  }
  
  exclude_cols <- c(poly_id_col, label_col, "LianaOriginal", "Type")
  numeric_cols <- sapply(df_class, is.numeric)
  numeric_cols[names(numeric_cols) %in% exclude_cols] <- FALSE
  
  distances <- numeric(n_repeats)
  
  set.seed(seed)
  for (i in 1:n_repeats) {
    shuffled_polys <- sample(unique_polys)
    split_idx <- floor(length(shuffled_polys) / 2)
    group1_polys <- shuffled_polys[1:split_idx]
    group2_polys <- shuffled_polys[(split_idx + 1):length(shuffled_polys)]
    
    data1 <- df_class[df_class[[poly_id_col]] %in% group1_polys, numeric_cols, drop = FALSE]
    data2 <- df_class[df_class[[poly_id_col]] %in% group2_polys, numeric_cols, drop = FALSE]
    
    data1 <- as.matrix(data1)
    data2 <- as.matrix(data2)
    
    var1 <- apply(data1, 2, var)
    var2 <- apply(data2, 2, var)
    good_bands <- (var1 > 1e-6) & (var2 > 1e-6)
    data1 <- data1[, good_bands, drop = FALSE]
    data2 <- data2[, good_bands, drop = FALSE]
    
    if (ncol(data1) == 0) {
      distances[i] <- NA
      next
    }
    
    mu1 <- colMeans(data1)
    mu2 <- colMeans(data2)
    sigma1 <- regularize_cov(cov(data1), epsilon)
    sigma2 <- regularize_cov(cov(data2), epsilon)
    
    dist_raw <- tryCatch({
      bhattacharyya.dist(mu1, mu2, sigma1, sigma2)
    }, error = function(e) NA)
    
    distances[i] <- extract_bhatt_dist(dist_raw)
  }
  
  valid_dists <- distances[!is.na(distances)]
  if (length(valid_dists) == 0) return(NA)
  
  list(
    mean_dist = mean(valid_dists),
    sd_dist = sd(valid_dists),
    min_dist = min(valid_dists),
    max_dist = max(valid_dists),
    n_repeats = n_repeats,
    n_valid = length(valid_dists),
    all_distances = distances
  )
}

add_to_plot <- function(plot_data, group_name, result, dataset, type) {
  if (is.list(result)) {
    rbind(plot_data, data.frame(
      group = rep(group_name, length(result$all_distances)),
      distance = result$all_distances,
      dataset = rep(dataset, length(result$all_distances)),
      type = rep(type, length(result$all_distances))
      ))
  } else {
    plot_data  
  }
}

# Mahalanobis ####

# maha_dist <- function(mu1, mu2, sigma1, sigma2) {
#   sigma_pooled <- (sigma1 + sigma2) / 2
#   delta <- mu1 - mu2
#   dist <- sqrt(t(delta) %*% solve(sigma_pooled) %*% delta)
#   return(dist)
# }
# 
# maha_between <- function(df, label_col, class1, class2, epsilon = 1e-10) { 
#   exclude_cols <- c("poly_id", "Liana", "LianaOriginal", "Type")
#   numeric_cols <- sapply(df, is.numeric)
#   numeric_cols[names(numeric_cols) %in% exclude_cols] <- FALSE
#   
#   data <- as.matrix(df[, numeric_cols, drop = FALSE])
#   labels <- df[[label_col]]
#   
#   d1 <- data[labels == class1, ]
#   d2 <- data[labels == class2, ]
#   
#   var_d1 <- apply(d1, 2, var)
#   var_d2 <- apply(d2, 2, var)
#   good_bands <- (var_d1 > 1e-6) & (var_d2 > 1e-6)
#   d1 <- d1[, good_bands, drop = FALSE]
#   d2 <- d2[, good_bands, drop = FALSE]
#   
#   if (ncol(d1) == 0) return(NA)
#   
#   mu1 <- colMeans(d1)
#   mu2 <- colMeans(d2)
#   sigma1 <- regularize_cov(cov(d1), epsilon) 
#   sigma2 <- regularize_cov(cov(d2), epsilon)  
#   
#   dist_raw <- tryCatch(maha_dist(mu1, mu2, sigma1, sigma2), error = function(e) NA)
#   if (is.na(dist_raw)) return(NA)
#   
#   list(dist = dist_raw, coeff = exp(-dist_raw), n1 = nrow(d1), n2 = nrow(d2), p = ncol(d1))
# }
# 
# maha_between_random <- function(df, label_col, class1, class2,
#                                 n_repeats = 100, seed = 123, epsilon = 1e-10) { 
#   df_filtered <- df[df[[label_col]] %in% c(class1, class2), ]
#   
#   exclude_cols <- c("poly_id", "Liana", "LianaOriginal", "Type")
#   numeric_cols <- sapply(df_filtered, is.numeric)
#   numeric_cols[names(numeric_cols) %in% exclude_cols] <- FALSE
#   
#   data <- as.matrix(df_filtered[, numeric_cols, drop = FALSE])
#   labels <- df_filtered[[label_col]]
#   
#   idx_class1 <- which(labels == class1)
#   idx_class2 <- which(labels == class2)
#   
#   distances <- numeric(n_repeats)
#   set.seed(seed)
#   for (i in 1:n_repeats) {
#     sample_idx1 <- sample(idx_class1, size = floor(length(idx_class1) * 0.5), replace = FALSE)
#     sample_idx2 <- sample(idx_class2, size = floor(length(idx_class2) * 0.5), replace = FALSE)
#     
#     d1 <- data[sample_idx1, ]
#     d2 <- data[sample_idx2, ]
#     
#     var_d1 <- apply(d1, 2, var)
#     var_d2 <- apply(d2, 2, var)
#     good_bands <- (var_d1 > 1e-6) & (var_d2 > 1e-6)
#     d1 <- d1[, good_bands, drop = FALSE]
#     d2 <- d2[, good_bands, drop = FALSE]
#     
#     if (ncol(d1) == 0) {
#       distances[i] <- NA
#       next
#     }
#     
#     mu1 <- colMeans(d1)
#     mu2 <- colMeans(d2)
#     sigma1 <- regularize_cov(cov(d1), epsilon)  
#     sigma2 <- regularize_cov(cov(d2), epsilon) 
#     
#     distances[i] <- tryCatch(maha_dist(mu1, mu2, sigma1, sigma2), error = function(e) NA)
#   }
#   
#   valid_dists <- distances[!is.na(distances)]
#   if (length(valid_dists) == 0) return(NA)
#   
#   list(
#     mean_dist = mean(valid_dists),
#     sd_dist = sd(valid_dists),
#     min_dist = min(valid_dists),
#     max_dist = max(valid_dists),
#     n_repeats = n_repeats,
#     n_valid = length(valid_dists),
#     all_distances = distances
#   )
# }
# 
# 
# maha_within_random <- function(df, poly_id_col, label_col, target_class,
#                                n_repeats = 100, seed = 123, epsilon = 1e-10) { 
#   df_class <- df[df[[label_col]] == target_class & !is.na(df[[poly_id_col]]) & df[[poly_id_col]] != 0, ]
#   unique_polys <- unique(df_class[[poly_id_col]])
#   
#   if (length(unique_polys) < 2) {
#     warning(paste("Less than two poly per class", target_class))
#     return(NA)
#   }
#   
#   exclude_cols <- c(poly_id_col, label_col, "LianaOriginal", "Type")
#   numeric_cols <- sapply(df_class, is.numeric)
#   numeric_cols[names(numeric_cols) %in% exclude_cols] <- FALSE
#   all_data <- as.matrix(df_class[, numeric_cols, drop = FALSE])
#   
#   distances <- numeric(n_repeats)
#   set.seed(seed)
#   for (i in 1:n_repeats) {
#     shuffled_polys <- sample(unique_polys)
#     split_idx <- floor(length(shuffled_polys) / 2)
#     group1_polys <- shuffled_polys[1:split_idx]
#     group2_polys <- shuffled_polys[(split_idx + 1):length(shuffled_polys)]
#     
#     data1 <- all_data[df_class[[poly_id_col]] %in% group1_polys, , drop = FALSE]
#     data2 <- all_data[df_class[[poly_id_col]] %in% group2_polys, , drop = FALSE]
#     
#     var1 <- apply(data1, 2, var)
#     var2 <- apply(data2, 2, var)
#     good_bands <- (var1 > 1e-6) & (var2 > 1e-6)
#     data1 <- data1[, good_bands, drop = FALSE]
#     data2 <- data2[, good_bands, drop = FALSE]
#     
#     if (ncol(data1) == 0 || nrow(data1) < 2 || nrow(data2) < 2) {
#       distances[i] <- NA
#       next
#     }
#     
#     mu1 <- colMeans(data1)
#     mu2 <- colMeans(data2)
#     sigma1 <- regularize_cov(cov(data1), epsilon) 
#     sigma2 <- regularize_cov(cov(data2), epsilon)  
#     
#     distances[i] <- tryCatch(maha_dist(mu1, mu2, sigma1, sigma2), error = function(e) NA)
#   }
#   
#   valid_dists <- distances[!is.na(distances)]
#   if (length(valid_dists) == 0) return(NA)
#   
#   list(
#     mean_dist = mean(valid_dists),
#     sd_dist = sd(valid_dists),
#     min_dist = min(valid_dists),
#     max_dist = max(valid_dists),
#     n_repeats = n_repeats,
#     n_valid = length(valid_dists),
#     all_distances = distances
#   )
# }

maha_sq_dist <- function(mu1, mu2, sigma_inv) {
  delta <- mu1 - mu2
  as.numeric(t(delta) %*% sigma_inv %*% delta)
}

spectral_matrix <- function(df,
                            exclude_cols = c("poly_id", "Liana",
                                             "LianaOriginal", "Type")) {
  numeric_cols <- sapply(df, is.numeric)
  numeric_cols[names(numeric_cols) %in% exclude_cols] <- FALSE
  X <- as.matrix(df[, numeric_cols, drop = FALSE])
  v <- apply(X, 2, var)
  X[, v > 1e-6, drop = FALSE]
}

# Common spectral space (ACP + global Sigma)
build_global_space <- function(df,
                               use_pca = TRUE,
                               n_pcs = NULL,
                               var_threshold = 0.99,
                               epsilon = 1e-6) {
  X <- spectral_matrix(df)
  if (ncol(X) == 0) stop("Aucune bande variable détectée.")
  
  if (use_pca) {
    pca <- prcomp(X, center = TRUE, scale. = FALSE)
    if (is.null(n_pcs)) {
      cum_var <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
      n_pcs <- max(2, which(cum_var >= var_threshold)[1])
    }
    n_pcs <- min(n_pcs, ncol(X), nrow(X) - 1)
    scores <- pca$x[, 1:n_pcs, drop = FALSE]
  } else {
    pca    <- NULL
    scores <- X
  }
  
  Sigma     <- regularize_cov(cov(scores), epsilon)
  list(pca = pca, scores = scores,
       Sigma = Sigma, Sigma_inv = solve(Sigma),
       n_pcs = ncol(scores))
}

make_split <- function(idx, pid, frac = 0.5) {
  polys   <- unique(pid[idx])
  polys_a <- sample(polys, floor(length(polys) * frac))
  list(
    A = idx[pid[idx] %in% polys_a],
    B = idx[!pid[idx] %in% polys_a]
  )
}

centroid_of <- function(scores, idx) {
  colMeans(scores[idx, , drop = FALSE])
}

summarize_dists <- function(x) {
  v <- x[!is.na(x)]
  if (length(v) == 0) return(NULL)
  list(mean_dist = mean(v), sd_dist = sd(v),
       min_dist = min(v), max_dist = max(v),
       n_valid = length(v), all_distances = x)
}
maha_global_comparisons <- function(df, label_col, class1, class2,
                                    poly_id_col = "poly_id",
                                    space = NULL,
                                    n_repeats = 500,
                                    frac = 0.5,
                                    seed = 123) {
  
  if (is.null(space)) {
    space <- build_global_space(df)
  }
  scores    <- space$scores
  Sigma_inv <- space$Sigma_inv
  
  labels <- df[[label_col]]
  pid    <- df[[poly_id_col]]
  
  idx1 <- which(labels == class1 & !is.na(pid) & pid != 0)
  idx2 <- which(labels == class2 & !is.na(pid) & pid != 0)
  if (length(idx1) < 4 || length(idx2) < 4)
    stop("Trop peu d'individus (poly_id valides) dans une des classes.")
  
  set.seed(seed)
  w1 <- w2 <- bt <- rep(NA_real_, n_repeats)
  
  for (i in seq_len(n_repeats)) {
    s1 <- make_split(idx1, pid, frac)
    s2 <- make_split(idx2, pid, frac)
    
    w1[i] <- tryCatch(sqrt(maha_sq_dist(centroid_of(scores, s1$A),
                                        centroid_of(scores, s1$B), Sigma_inv)),
                      error = function(e) NA)
    w2[i] <- tryCatch(sqrt(maha_sq_dist(centroid_of(scores, s2$A),
                                        centroid_of(scores, s2$B), Sigma_inv)),
                      error = function(e) NA)
    bt[i] <- tryCatch(sqrt(maha_sq_dist(centroid_of(scores, s1$A),
                                        centroid_of(scores, s2$A), Sigma_inv)),
                      error = function(e) NA)
  }
  
  w1s <- summarize_dists(w1)
  w2s <- summarize_dists(w2)
  bts <- summarize_dists(bt)
  
  sep_z <- NA_real_
  if (!is.null(w1s) && !is.null(w2s) && !is.null(bts)) {
    pooled_sd <- sqrt((w1s$sd_dist^2 + w2s$sd_dist^2) / 2)
    if (pooled_sd > 0)
      sep_z <- (bts$mean_dist - mean(c(w1s$mean_dist, w2s$mean_dist))) / pooled_sd
  }
  
  list(within_class1 = w1s,
       within_class2 = w2s,
       between       = bts,
       sep_z         = sep_z,
       n_pcs         = space$n_pcs,
       seed          = seed)
}

maha_results_to_df <- function(res, dataset, comparison) {
  if (is.null(res$between)) {
    return(data.frame(dataset = dataset, comparison = comparison,
                      mean_dist = NA, sd_dist = NA,
                      within1 = NA, within2 = NA,
                      sep_z = NA, n_pcs = res$n_pcs))
  }
  data.frame(
    dataset    = dataset,
    comparison = comparison,
    mean_dist  = res$between$mean_dist,
    sd_dist    = res$between$sd_dist,
    within1    = res$within_class1$mean_dist,
    within2    = res$within_class2$mean_dist,
    sep_z      = res$sep_z,
    n_pcs      = res$n_pcs
  )
}

# LDA ####
stratified_group_kfold <- function(data, group_col, stratify_col, k = 5, seed = 123) {
  set.seed(seed)
  group_classes <- sapply(unique(data[[group_col]]), function(group) {
    group_data <- data[data[[group_col]] == group, ]
    names(which.max(table(group_data[[stratify_col]])))
  })
  groups_df <- data.frame(group = unique(data[[group_col]]), class = group_classes)
  folds <- vector("list", k)
  for (class in unique(groups_df$class)) {
    class_groups <- groups_df$group[groups_df$class == class]
    fold_assignments <- sample(rep(1:k, length.out = length(class_groups)))
    for (i in 1:k) folds[[i]] <- c(folds[[i]], class_groups[fold_assignments == i])
  }
  lapply(folds, function(groups) which(data[[group_col]] %in% groups))
}

calculate_f1_manual <- function(true_labels, predicted_labels) {
  cm <- table(true_labels, predicted_labels)
  f1_scores <- numeric(ncol(cm))
  names(f1_scores) <- colnames(cm)  
  for (i in 1:ncol(cm)) {
    cls <- colnames(cm)[i]          
    tp <- cm[cls, cls]              
    fp <- sum(cm[, cls]) - tp       
    fn <- sum(cm[cls, ]) - tp       
    precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    f1_scores[i] <- ifelse((precision + recall) == 0, 0, 2 * precision * recall / (precision + recall))
  }
    f1_macro <- mean(f1_scores)
  return(f1_macro)
}