# Abbreviations used #############################################

# md: manual delineation (blob dataset)
# dp: drone photos (crowns dataset)
# poly: polygons (opposite to pixel analysis)
# m : matrix

# File initialization #############################################

rm(list = ls())

## Import libraries #### 
library(tidyverse)
library(sf)
library(terra)
library(data.table)
library(stats)
library(Rtsne)
library(ggplot2)
library(paletteer)
library(uwot)
library(FactoMineR)
library(fpc)
library(MASS)
library(caret)

## Set working Directory ####
setwd ("D:/DataCamille/")


# Uploading data and cropping ##########################################

## HSI mosaic ####

# r <- rast("HSI/mosaic_2016/nou_20160921_l1c_mosaic.bsq") #2016
r <- rast("HSI/mosaic_2023/mosaic/nou_20231025_l1c_mosaic.bsq") #2023
nlyr(r) #number of spectral bands
# plot(r[[470]])
# hist(r[[470]], breaks = 100)

## Illu ####

illu <- rast("HSI/mosaic_2016/nou_20160921_ilu_mosaic.bsq")
# illu <- rast("HSI/mosaic_2023/mosaic/nou_20231025_l1c_mosaic.bsq")
nlyr(illu)
# hist(illu, breaks = 100)
# plot(illu)

## Manually delineated (MD) ####

blob <- vect("Lianes/LianaOrNotLiana_55_45.shp")
blob <- project(blob, crs(r))

area_blob <- crop(r, blob, mask=TRUE)
# plot (area_blob) # /!\ creates one plot per band
# plot(area_blob[[2]]) # the number is the rank of the band not it's value

# Taking into account the liana label
label_blob <- rasterize(blob, area_blob, field = "Type")

# Normalized data
normalized_blob <- area_blob / mean(area_blob, na.rm = TRUE)


## Crowns from drone photos (DP) ####
crowns <- vect("Couronnes/Crowns_Liana_2025_GP_TRUE.shp")
crowns <- project(crowns, crs(r))

area_crowns <- crop(r, crowns, mask=TRUE)

label_crowns <- rasterize(crowns, area_crowns, field = "LianaFullC")

normalized_crowns <- area_crowns / mean(area_crowns, na.rm =TRUE)

# Filtering pixels ##########################################################

## Cleaning extreme values ##
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



## Removing pixel based on illu ####
# Only used with the 2016 data (illu file missing fo 2023)

med_illu <- global(illu[[1]], fun = median, na.rm = TRUE)$global
perc90_illu <- global(illu[[1]], fun = quantile, probs = 0.9, na.rm = TRUE)$X90.

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

## Removing pixel based on VNIR reflectance median ####

extract_wavelength <- function(band_names) {
  
  wl <- suppressWarnings(as.numeric(band_names))
  need_regex <- is.na(wl)
  wl[need_regex] <- as.numeric(gsub(".*\\(([0-9.]+)\\).*", "\\1", band_names[need_regex]))
  
  return(wl)
}


filter_by_vnir <- function(area, label, vnir_range = c(400, 1100)){
  
  area <- clean_extreme_values(area)$area_clean
  
  area_matrix <- as.matrix(area)
  na_per_band <- colSums(is.na(area_matrix))
  good_bands <- na_per_band < nrow(area_matrix)
  
  wl <- extract_wavelength(names(area))[good_bands]
  
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


## Removing pixel based on NIR reflectance####

# names(area_crowns)
# crowns_nir <- area_crowns[[118:157]]
# mean_crowns <- mean(crowns_nir, na.rm =TRUE)
# hist(mean_crowns)

filter_by_nir <- function(area, label, nir_range = c(780, 900)){

  area <- clean_extreme_values(area)$area_clean
  

  area_matrix <- as.matrix(area)
  na_per_band <- colSums(is.na(area_matrix))
  good_bands <- na_per_band < nrow(area_matrix)

  wl <- extract_wavelength(names(area))[good_bands]
  
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



## Applying/choosing the filter #########################################################
# Comment and/or uncomment to choose the right filter

### Illu ####
# illu_crowns_out <- filter_by_illu(illu, normalized_crowns, label_crowns, med_illu)
# area_crowns_filt <- illu_crowns_out$area_filt
# label_crowns_filt <- illu_crowns_out$label_filt
# illu_crop_crowns <- illu_crowns_out$area_crop
# 
# illu_blob_out <- filter_by_illu(illu, normalized_blob, label_blob, med_illu)
# area_blob_filt <- illu_blob_out$area_filt
# label_blob_filt <- illu_blob_out$label_filt
# illu_crop_blob <- illu_blob_out$area_crop

### VNIR ####
# spec_crowns_filt_VNIR <- filter_by_vnir(area_crowns, label_crowns)
# spec_crowns_filt_VNIR$threshold_qual
# spec_crowns_filt_VNIR$threshold
# area_crowns_filt <- spec_crowns_filt_VNIR$area_filt
# label_crowns_filt<- spec_crowns_filt_VNIR$label_filt
# illu_crop_crowns <- area_crowns_filt *0
# 
# spec_blob_filt_VNIR   <- filter_by_vnir(area_blob, label_blob)
# spec_blob_filt_VNIR$threshold_qual
# spec_blob_filt_VNIR$threshold
# area_blob_filt <- spec_blob_filt_VNIR$area_filt
# label_blob_filt <- spec_blob_filt_VNIR$label_filt
# illu_crop_blob <- area_blob_filt * 0

### NIR ####
spec_crowns_filt_NIR <- filter_by_nir(area_crowns, label_crowns)
spec_crowns_filt_NIR$threshold_qual
spec_crowns_filt_NIR$threshold_value
area_crowns_filt <- spec_crowns_filt_NIR$area_filt
label_crowns_filt<- spec_crowns_filt_NIR$label_filt
spec_crowns_filt_NIR$threshold_qual  # % remaining pixels
illu_crop_crowns <-area_crowns_filt *0

spec_blob_filt_NIR <- filter_by_nir(area_blob, label_blob)
spec_blob_filt_NIR$threshold_qual
spec_blob_filt_NIR$threshold_value
area_blob_filt <- spec_blob_filt_NIR$area_filt
label_blob_filt <- spec_blob_filt_NIR$label_filt
spec_blob_filt_NIR$threshold_qual
illu_crop_blob <- area_blob_filt * 0

# Dimensionality reduction and data visualization ###############################

## Data prep ##############################################################

data_prep_per_pixel <- function(area_filtered, label_filtered, illu_crop_area, max_na_pct = 0.05){
  area_matrix <- as.matrix(area_filtered)
  label_matrix <- as.matrix (label_filtered)[,1]
  
  na_per_band <- colSums(is.na(area_matrix))
  good_bands <- na_per_band < nrow(area_matrix)
  # sum(good_bands)
  
  # na_pct_per_pixel <- rowSums(is.na(area_matrix[, good_bands])) / sum(good_bands)
  # no_na_row <- na_pct_per_pixel <= max_na_pct   
  # valid_idx <- which(!is.na(label_matrix) & no_na_row)
  # valid_idx <- which(!is.na(label_matrix) & complete.cases(area_matrix[, good_bands]))
  valid_idx <- which(!is.na(label_matrix) & rowSums(is.na(area_matrix[, good_bands])) == 0)
  # length(valid_idx)

  area_valid <- area_matrix[valid_idx, good_bands, drop =FALSE]
  labels_valid <- label_matrix[valid_idx]
  # dim(crowns_filt_valid)
  # length(labels_valid)
  
  illu_matrix <- as.matrix(illu_crop_area)[, 1]
  illu_valid <- illu_matrix[valid_idx]
  
  # if (!is.null(illu_crop_area)) {
  #   illu_matrix <- as.matrix(illu_crop_area)[, 1]
  #   illu_valid <- illu_matrix[valid_idx]
  # } else {
  #   illu_valid <- NULL
  # }
  
  output_list <- list(
    matrix = area_valid,
    labels = labels_valid,
    illu = illu_valid,
    goodbands = good_bands)
  
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

### Data prep for DP
DP_per_pixel <- data_prep_per_pixel(area_crowns_filt, label_crowns_filt, illu_crop_crowns)

crowns_filt_valid <- DP_per_pixel$matrix
illu_crowns_valid <- DP_per_pixel$illu
label_crowns_valid <- DP_per_pixel$labels


DP_per_poly <- data_prep_per_poly(area_crowns_filt, crowns, "LianaFullC")

crowns_poly_valid <- DP_per_poly$matrix
label_crowns_poly_valid <- DP_per_poly$label


### Data prep for MD
MD_per_pixel <- data_prep_per_pixel(area_blob_filt, label_blob_filt, illu_crop_blob)

blob_filt_valid <- MD_per_pixel$matrix
illu_blob_valid <- MD_per_pixel$illu
label_blob_valid <- MD_per_pixel$label


MD_per_poly <- data_prep_per_poly(area_blob_filt, blob, "Type")

blob_poly_valid <- MD_per_poly$matrix
label_blob_poly_valid <- MD_per_poly$label



## t-SNE #######################################################################

## DP per pixel
# perplexity_crown = (nrow(crowns_filt_valid)-1)/3 - 1
# perplexity_crown
set.seed(2)
tsne_out_crowns <- Rtsne(crowns_filt_valid, check_duplicates = FALSE )

tsne_plot_crowns <- data.frame(
  x = tsne_out_crowns$Y[,1],
  y = tsne_out_crowns$Y[,2],
  LianaFullCrown = as.factor(label_crowns_valid)
)

tsne_plot_crowns$illu <- illu_crowns_valid

ggplot(tsne_plot_crowns, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("0"= "olivedrab3", "1"="deeppink3","2"="purple4","3"= "orange3"),
                     labels = c("0" = "CR0", "1" = "CR1", "2"= "CR2", "3"= "CR3"))+
  labs(title= "tSNE of the HSI values of the drone photos dataset (NIR filter)" )

# Only useful if the illumination filter was used
# ggplot(tsne_plot_crowns, aes(x = x, y = y, color = illu)) +
#   geom_point() +
#   scale_color_viridis_c(name = "Illumination") +
#   labs(title = "tSNE of the HSI values of the DP dataset colored by illumination")

## For DP per polygon 
perplexity_crown_poly = (nrow(crowns_poly_valid)-1)/3 - 1
# perplexity_crown_poly
set.seed(2)
tsne_out_crowns_poly <- Rtsne(crowns_poly_valid, perplexity = perplexity_crown_poly, check_duplicates = FALSE)
tsne_plot_crowns_poly <- data.frame(
  x = tsne_out_crowns_poly$Y[,1],
  y = tsne_out_crowns_poly$Y[,2],
  LianaFullCrown = as.factor(label_crowns_poly_valid)
)

ggplot(tsne_plot_crowns_poly,aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+ #(size = 3)
  scale_color_manual(values = c("CR0"= "olivedrab3", "CR1"="deeppink3","CR2"="purple4","CR3"= "orange3"))+
  # scale_color_viridis_d(option="turbo")+
  labs(title= "tSNE of the mean spectrum per polygon of the DP dataset (NIR filter)" )

## For MD per pixel
# perplexity_blob = (nrow(blob_filt_valid)-1)/3 - 1
# perplexity_blob
set.seed(2)
tsne_out_blob <- Rtsne( blob_filt_valid, check_duplicates = FALSE)

tsne_plot_blob <- data.frame(
  x = tsne_out_blob$Y[,1],
  y = tsne_out_blob$Y[,2],
  Type = as.factor(label_blob_valid)
)

tsne_plot_blob$illu <- illu_blob_valid


ggplot(tsne_plot_blob, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "olivedrab3", "1" = "deeppink2"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "tSNE of the HSI values of the MD dataset (NIR filter)" )

# Only with illu filter
# ggplot(tsne_plot_blob, aes(x = x, y = y, color = illu)) +
#   geom_point() +
#   scale_color_viridis_c(name = "Illumination") +
#   labs(title = "tSNE of the HSI values of the MD dataset colored by illumination")

## For MD per polygon 
perplexity_blob_poly = (nrow(blob_poly_valid)-1)/3 - 1
# perplexity_blob_poly
set.seed(2)
tsne_out_blob_poly <- Rtsne(blob_poly_valid, perplexity = perplexity_blob_poly , check_duplicates = FALSE)

tsne_plot_blob_poly <- data.frame(
  x = tsne_out_blob_poly$Y[,1],
  y = tsne_out_blob_poly$Y[,2],
  Type = as.factor(label_blob_poly_valid)
)

ggplot(tsne_plot_blob_poly, aes(x = x, y = y, color = Type)) +
  geom_point() + #(size = 3)
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "olivedrab3", "1" = "deeppink2"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "tSNE of the mean spectrum per polygon of the MD dataset (NIR filter)" )


## UMAP #########################################################################

#For DP per pixel

umap_out_crowns <- umap2( seed = 2, crowns_filt_valid,  n_neighbors = min(15, nrow(crowns_filt_valid) - 1))

umap_plot_crowns <- data.frame(
  x = umap_out_crowns[,1],
  y = umap_out_crowns[,2],
  LianaFullCrown = as.factor(label_crowns_valid)
)

# umap_plot_crowns$illu <- illu_crowns_valid
  
  
ggplot(umap_plot_crowns, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("0"= "olivedrab3", "1"="deeppink3","2"="purple4","3"= "orange3"),
                     labels = c("0" = "CR0", "1" = "CR1", "2"= "CR2", "3"= "CR3"))+
  # scale_color_viridis_d(option="turbo")+
  labs(title= "UMAP of the HSI values of the DP dataset (NIR filter)" )

# Only with illu filter
# ggplot(umap_plot_crowns, aes(x = x, y = y, color = illu)) +
#   geom_point() +
#   scale_color_viridis_c(name = "Illumination") +
#   labs(title = "UMAP of the HSI values of the DP dataset colored by illumination")

#For DP per polygon
umap_out_crowns_poly <- umap2( seed = 2, crowns_poly_valid)

umap_plot_crowns_poly <- data.frame(
  x = umap_out_crowns_poly[,1],
  y = umap_out_crowns_poly[,2],
  LianaFullCrown = as.factor(label_crowns_poly_valid))

ggplot(umap_plot_crowns_poly, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("CR0"= "olivedrab3", "CR1"="deeppink3","CR2"="purple4","CR3"= "orange3"))+
  labs(title= "UMAP of the mean HSI value per polygon of the DP dataset (NIR filter)" )


#For MD per pixel
umap_out_blob <- umap2( seed = 2, blob_filt_valid)

umap_plot_blob <- data.frame(
  x = umap_out_blob[,1],
  y = umap_out_blob[,2],
  Type = as.factor(label_blob_valid)
)

umap_plot_blob$illu <- illu_blob_valid

ggplot(umap_plot_blob, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "olivedrab3", "1" = "deeppink2"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "UMAP of the HSI values of the MD dataset (NIR filter)" )

# ggplot(umap_plot_blob, aes(x = x, y = y, color = illu)) +
#   geom_point() +
#   scale_color_viridis_c(name = "Illumination") +
#   labs(title = "UMAP of the HSI values of the MD dataset colored by illumination")

#For MD per polygon
umap_out_poly <- umap2( seed = 2, blob_poly_valid)

umap_plot_poly <- data.frame(
  x = umap_out_poly[,1],
  y = umap_out_poly[,2],
  Type = as.factor(label_blob_poly_valid)
)

ggplot(umap_plot_poly, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "olivedrab3", "1" = "deeppink2"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "UMAP of the mean HSI value per polygon of the MD dataset (NIR filter)" )

# # Bhattacharyya distance ###################################################
# 
# bhatt_dist_pairs <- function(pixel_matrix, groups) {
#   
#   ids <- unique(groups)
#   
#   # calculate mu & Sigma for each group
#   group_stats <- lapply(ids, function(id) {
#     px <- pixel_matrix[groups == id, , drop = FALSE] #drop = FALSE keep the class as matrix
#     if (nrow(px) <= ncol(px)) return(NULL)  # check if n_pixels/n_bands so sigma can be inverted
#     list(mu = colMeans(px), Sigma = cov(px), n = nrow(px)) #else
#   })
#   names(group_stats) <- ids
#   group_stats <- Filter(Negate(is.null), group_stats) #remove groupe with not enough pixels , i.e for which sigma was not calculated
#   
#   valid_ids <- names(group_stats) 
#   if (length(valid_ids) < 2) { #we need at least two groups to calculate a distance
#     warning("Not enough groups to calculate a distance")
#     return(NULL)
#   }
#   
#   #Calculate the number of pairs
#   n_pairs <- length(valid_ids) * (length(valid_ids) - 1) / 2
#   
#   #Create an empty df
#   results <- data.frame(
#     group1 = character(n_pairs),
#     group2 = character(n_pairs),
#     dist_bhatt = numeric(n_pairs),
#     coeff_bhatt = numeric(n_pairs),
#     n1 = numeric(n_pairs),
#     n2 = numeric(n_pairs)
#   )
#   
#   # Check every possible pairs and calculate mu and sigma
#   row_i <- 1
#   for (i in 1:(length(valid_ids) - 1)) {
#     for (j in (i + 1):length(valid_ids)) {
#       
#       g1 <- valid_ids[i]
#       g2 <- valid_ids[j]
#       
#       s1 <- group_stats[[g1]]
#       s2 <- group_stats[[g2]]
#       
#       # Calculate the bhattacharya distance
#       d <- tryCatch(
#         bhattacharyya.dist(s1$mu, s2$mu, s1$sigma, s2$sigma),
#         error = function(e) NA)
#       
#       results$group1[row_i] <- g1
#       results$group2[row_i] <- g2
#       results$dist_bhatt[row_i] <- d
#       results$coeff_bhatt[row_i] <- exp(-d)
#       results$n1[row_i] <- s1$n
#       results$n2[row_i] <- s2$n
#       
#       row_i <- row_i + 1
#     }
#   }
#   
#   dist_df <- results
#   
#   return(dist_df)
# }
# 
# # bhatt_dist_btwn <- function(pixel_matrix, pixel_label, label1, label2) {
# #   keep <- pixel_label %in% c(label1, label2)
# #   bhatt_dist_pairs(pixel_matrix[keep, ], pixel_label[keep])
# # }
# 
# bhatt_dist_within <- function(area_filtered, vect_obj, label_name, target_label) {
#   
#   vect_obj$poly_id <- seq_len(nrow(vect_obj))
#   poly_id_r <- rasterize(vect_obj, area_filtered, field = "poly_id")
#   
#   area_matrix <- as.matrix(area_filtered)
#   poly_id_m <- as.matrix(poly_id_r)[, 1]
#   class_by_id <- values(vect_obj)[[label_name]]
#   
#   na_per_band <- colSums(is.na(area_matrix))
#   good_bands <- na_per_band < nrow(area_matrix)
#   
#   valid_idx <- which(!is.na(poly_id_m) & complete.cases(area_matrix[, good_bands]))
#   
#   pixel_matrix <- area_matrix[valid_idx, good_bands]
#   pixel_poly_id <- poly_id_m[valid_idx]
#   
#   keep_ids <- which(class_by_id == target_label)
#   keep_px <- pixel_poly_id %in% keep_ids
#   
#   bhatt_dist_pairs(pixel_matrix[keep_px, ], pixel_poly_id[keep_px])
# }
# 
# 
# # ##Distance between classes
# # 
# bhatt_dist_btwn <- function (matrix, matrix_label, label1, label2){
# 
#   # #Reduce dimensionality
#   # max_ncp <- min(nrow(matrix) - 1, ncol(matrix))
#   # res_pca <- PCA(matrix, scale.unit = TRUE, ncp = max_ncp, graph = FALSE)
#   # cum_var <- res_pca$eig[, "cumulative percentage of variance"]
#   # n_components <- which(cum_var >= 85)[1]
#   #
#   # matrix_reduced <- res_pca$ind$coord[, 1:n_components]
#   #
#   # plot(cum_var)
# 
#   matrix_reduced <- matrix
# 
#   #Separate groups
#   group0 <- matrix_reduced[matrix_label == label1, ]
#   group1 <- matrix_reduced[matrix_label == label2, ]
# 
# 
#   # Test if the distribution is Gaussian first
#   shapiro_group0 <- apply(group0, 2, function(x) shapiro.test(x)$p.value)
#   shapiro_group1 <- apply(group1, 2, function(x) shapiro.test(x)$p.value)
# 
#   # par(mfrow = c(2, 3))
#   # for (i in 1:n_components) {
#   #   qqnorm(group0[, i], main = paste("Composante", i, "- groupe 0"))
#   #   qqline(group0[, i])
#   # }
#   #
#   # par(mfrow = c(2, 3))
#   # for (i in 1:n_components) {
#   #   qqnorm(group1[, i], main = paste("Composante", i, "- groupe 1"))
#   #   qqline(group1[, i])
#   # }
# 
#   # Calculate mean vectors
#   mu0 <- colMeans(group0)
#   mu1 <- colMeans(group1)
# 
#   # Calculate covariance matrix
#   sigma0 <- cov(group0)
#   sigma1 <- cov(group1)
# 
#   
#     
#   dist_bhatt <- bhattacharyya.dist(mu0, mu1, sigma0, sigma1)
# 
#   coeff_bhatt <- exp(-dist_bhatt)
# 
#   output <- list(
#     mean0 = mu0,
#     mean1 = mu1,
#     dist_bhatt = dist_bhatt,
#     coeff_bhatt = coeff_bhatt,
#     shapiro_group0 = shapiro_group0,
#     shapiro_group1 = shapiro_group1 #,
#     # n_components = n_components
#   )
# 
#   return(output)
# }
# 
# # With DP dataset
# 
# separability_dp_between <- bhatt_dist_btwn(crowns_poly_valid[,1:135], label_crowns_poly_valid, "CR0", "CR3")
# 
# #Check if we keep the Bhattacharya distance
# separability_dp_between$shapiro_group0
# separability_dp_between$shapiro_group0 < 0.05
# 
# separability_dp_between$shapiro_group1
# separability_dp_between$shapiro_group1 < 0.05
# 
# #Results
# bhatt_dist_dp_LianaVSTrees <- separability_dp_between$dist_bhatt
# bhatt_dist_dp_LianaVSTrees
# 
# bhatt_coeff_dp_LianaVSTrees <-separability_dp_between$coeff_bhatt
# bhatt_coeff_dp_LianaVSTrees
# 
# ## With MD dataset 
# separability_md_within <- bhatt_dist_btwn(blob_poly_valid, label_blob_poly_valid, 0, 1)
# 
# separability_md_within$shapiro_group0
# separability_md_within$shapiro_group0<0.05
# 
# separability_md_within$shapiro_group1
# separability_md_within$shapiro_group1<0.05
# 
# bhatt_dist_md_LianaVSTrees <- separability_md_within$dist_bhatt
# bhatt_dist_md_LianaVSTrees
# 
# bhatt_coeff_md_LianaVSTrees <- separability_md_within$coeff_bhatt
# bhatt_coeff_md_LianaVSTrees

## Distance within classes ####

# bhatt_dist_within <- function(area_filtered, area_vector, label_name, target_label){
#
#   #Get unique id for each polygon
#   area_vector$poly_id <- seq_len(nrow(area_vector))
#
#   poly_id_r <- rasterize(area_vector, area_filtered, field = "poly_id")
#
#   area_matrix <- as.matrix(area_filtered)
#   poly_id_m <- as.matrix(poly_id_r)[, 1]
#   class_by_id <- values(area_vector)[[label_name]]
#
#   na_per_band <- colSums(is.na(area_matrix))
#   good_bands <- na_per_band < nrow(area_matrix)
#
#   valid_idx <- which(!is.na(poly_id_m) & complete.cases(area_matrix[, good_bands]))
#
#   pixel_matrix <- area_matrix[valid_idx, good_bands]
#   pixel_poly_id <- poly_id_m[valid_idx]
#
#   # Keep only the pixels of the target label
#   keep_ids <- which(class_by_id == target_label)
#   keep_px <- pixel_poly_id %in% keep_ids
#
#   pixel_matrix <- pixel_matrix[keep_px, ]
#   pixel_poly_id <- pixel_poly_id[keep_px]
#
#   # Dimensionality reduction
#   max_ncp <- min(nrow(pixel_matrix) - 1, ncol(pixel_matrix))
#   res_pca <- PCA(pixel_matrix, scale.unit = TRUE, ncp = max_ncp, graph = FALSE)
#   cum_var <- res_pca$eig[, "cumulative percentage of variance"]
#   n_components <- which(cum_var >= 85)[1]
#
#   pixel_reduced <- res_pca$ind$coord[, 1:n_components]
#
#   # Calculate mu & sigma for each polygon
#   polygon_ids <- unique(pixel_poly_id)
#
#   polygon_stats <- lapply(polygon_ids, function(id) {
#     px <- pixel_reduced[pixel_poly_id == id, , drop = FALSE]
#     if (nrow(px) <= n_components) return(NULL)  # not enough pixels to calculate sigma
#     list(mu = colMeans(px), Sigma = cov(px), n = nrow(px))
#   })
#   names(polygon_stats) <- polygon_ids
#   polygon_stats <- Filter(Negate(is.null), polygon_stats)
#
#   ids <- names(polygon_stats)
#   if (length(ids) < 2) {
#     warning("Not enough valid polygons to calculate intra-class distances.")
#     return(NULL)
#   }
#
#   # Distance by pairs
#   pairs <- combn(ids, 2, simplify = FALSE)
#
#   dist_df <- do.call(rbind, lapply(pairs, function(p) {
#     s1 <- polygon_stats[[p[1]]]
#     s2 <- polygon_stats[[p[2]]]
#     d <- tryCatch(
#       bhattacharyya.dist(s1$mu, s2$mu, s1$Sigma, s2$Sigma),
#       error = function(e) NA
#     )
#     data.frame(poly1 = p[1], poly2 = p[2], dist_bhatt = d, coeff_bhatt = exp(-d))
#   }))
#
#   output <- list(
#     distances = dist_df,
#     n_polygons = length(ids),
#     n_components = n_components,
#     mean_dist = mean(dist_df$dist_bhatt, na.rm = TRUE),
#     mean_coeff = mean(dist_df$coeff_bhatt, na.rm = TRUE)
#   )
#
#   return(output)
# }
#
# #Results
#
# bhatt_within_dp_CR0 <- bhatt_dist_within(illu_crowns_out$area_filt, crowns, "LianaFullC", "CR0")
# bhatt_within_dp_CR0$mean_dist
# bhatt_within_dp_CR0$mean_coeff
# sum(is.na(bhatt_within_dp_CR0$distances$dist_bhatt))   # number of failed pairs (singular covariance, etc.)
#
# bhatt_within_dp_CR1 <- bhatt_dist_within(illu_crowns_out$area_filt, crowns, "LianaFullC", "CR1")
# bhatt_within_dp_CR2 <- bhatt_dist_within(illu_crowns_out$area_filt, crowns, "LianaFullC", "CR2")
#
# bhatt_within_dp_CR3 <- bhatt_dist_within(illu_crowns_out$area_filt, crowns, "LianaFullC", "CR3")
# bhatt_within_dp_CR3$mean_dist
# bhatt_within_dp_CR3$mean_coeff
# sum(is.na(bhatt_within_dp_CR0$distances$dist_bhatt))
#
# bhatt_within_md_0 <- bhatt_dist_within(illu_blob_out$area_filt, blob, "Type", "0")
# bhatt_within_md_0$mean_dist
# bhatt_within_md_0$mean_coeff
# sum(is.na(bhatt_within_md_0$distances$dist_bhatt))
#
# bhatt_within_md_1 <- bhatt_dist_within(illu_blob_out$area_filt, blob, "Type", "1")
# bhatt_within_md_1$mean_dist
# bhatt_within_md_1$mean_coeff
# sum(is.na(bhatt_within_md_1$distances$dist_bhatt))
#

## Results boxplot ####

# between_df <- data.frame(
#   group = c("DP CR0 vs CR3", "MD 0 vs 1"),
#   dist_bhatt = c(bhatt_dist_dp_LianaVSTrees, bhatt_dist_md_LianaVSTrees),
#   coeff_bhatt = c(bhatt_coeff_dp_LianaVSTrees, bhatt_coeff_md_LianaVSTrees),
#   type = "Between"
# )
#
# combine_within <- function(result, group_name) {
#   df <- result$distances
#   df$group <- group_name
#   df
# }
#
# all_within <- rbind(
#   combine_within(bhatt_within_dp_CR0, "DP CR0"),
#   combine_within(bhatt_within_dp_CR1, "DP CR1"),
#   combine_within(bhatt_within_dp_CR2, "DP CR2"),
#   combine_within(bhatt_within_dp_CR3, "DP CR3"),
#   combine_within(bhatt_within_md_0, "MD 0"),
#   combine_within(bhatt_within_md_1, "MD 1")
# )
#
# ggplot(all_within, aes(x = group, y = dist_bhatt, fill = group)) +
#   geom_boxplot(outlier.shape = NA) +
#   coord_cartesian(ylim = c(0, quantile(all_within$dist_bhatt, 0.95, na.rm = TRUE))) +
#   labs(x = "Group",
#        y = "Bhattacharyya distance",
#        title = "Bhattacharyya distance within classes") +
#   theme(legend.position = "none")
#
# all_within$type <- "Within"
#
# combined <- rbind(
#   all_within[, c("group", "dist_bhatt", "coeff_bhatt", "type")],
#   between_df
# )
#
# combined$group <- factor(combined$group, levels = c(
#   "DP CR0", "DP CR1", "DP CR2", "DP CR3", "MD 0", "MD 1",
#   "DP CR0 vs CR3", "MD 0 vs 1"
# ))
#
#
# ggplot(combined, aes(x = group, y = dist_bhatt)) +
#   geom_boxplot(
#     data = subset(combined, type == "Within"),
#     aes(fill = group), outlier.shape = NA) +
#   geom_point(
#     data = subset(combined, type == "Between"),
#     color = "red", size = 3, shape = 18) +
#   geom_text(
#     data = subset(combined, type == "Between"),
#     aes(label = round(dist_bhatt, 3)),
#     color = "red", vjust = -1, size = 3.5, fontface = "bold")+
#   scale_y_log10() +
#   labs(
#     x = "Group",
#     y = "Bhattacharyya distance (log scale)",
#     title = "Bhattacharya distance : within classes (boxplots) VS between classes (red dots)") +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 45, hjust = 1))
#


# Mean spectrum comparison ###################################################

mean_spectrum_by_group <- function(matrix_valid, labels_valid, dataset_name) {
  groups <- unique(labels_valid)

  wavelengths <- extract_wavelength(colnames(matrix_valid))
  
  spectra_list <- lapply(groups, function(g) {
    mean_spec <- colMeans(matrix_valid[labels_valid == g, , drop = FALSE], na.rm = TRUE)
    data.frame(
      wavelength = wavelengths,
      reflectance = mean_spec,
      group = paste0(dataset_name, "_", g)
    )
  })

  do.call(rbind, spectra_list)
}
  
spec_crowns <- mean_spectrum_by_group(crowns_poly_valid, label_crowns_poly_valid, "DP")
spec_blob   <- mean_spectrum_by_group(blob_poly_valid, label_blob_poly_valid, "MD")

all_spectra <- rbind(spec_crowns, spec_blob)

#selected <- subset(all_spectra, group %in% c("DP_CR0", "DP_CR3", "MD_0", "MD_1") & all_spectra$wavelength<900) #VNIR only
selected <- subset(all_spectra, group %in% c("DP_CR0", "DP_CR3", "MD_0", "MD_1"))

label_config <- data.frame(
  group = c("DP_CR0", "DP_CR1", "DP_CR2", "DP_CR3", "MD_0", "MD_1"),
  display_name = c("CR0 ( No Liana DP)", "CR1 (DP)", "CR2", "CR3 (Liana DP)", "No liana (MD)", "Liana (MD)"),
  color = c("yellowgreen", "darkorange1", "purple4", "deeppink3", "darkgreen", "pink")
)

ggplot(selected, aes(x = wavelength, y = reflectance, color = group)) +
  geom_line(linewidth = 1) +
  labs(x = "Wavelength (nm)", y = "Mean reflectance",
       title = "Mean spectrum by dataset group (NIR filter)", color = "Type") +
  scale_color_manual(
    values = setNames(label_config$color, label_config$group),
    labels = setNames(label_config$display_name, label_config$group))

# Classification ##############################################################

## Linear Discriminant Analysis (LDA) #########################################


## Within DP dataset
crowns_df <- as.data.frame(crowns_filt_valid)
crowns_df$ Liana <- as.factor(label_crowns_valid)

set.seed(123)
crowns_df[,1:ncol(crowns_filt_valid)]<- scale(crowns_df[,1:ncol(crowns_filt_valid)])
sample_crowns <- sample(c(TRUE, FALSE), nrow(crowns_df), replace=TRUE, prob=c(0.7,0.3))
train_crowns <- crowns_df[sample_crowns, ]
test_crowns <- crowns_df[!sample_crowns, ] 

model_crowns <- lda(Liana ~ ., data = train_crowns)
predictions_crowns <- predict(model_crowns, test_crowns)
lda_crowns_accuracy <- mean(predictions_crowns$class==test_crowns$Liana)
lda_crowns_accuracy #0.6187719

### Stratified folded cross-validation
predictors <- crowns_df[ , -which(names(crowns_df) %in% c("Liana"))]

ctrl_crowns <- trainControl(
  method = "cv",
  number = 5,        
  savePredictions = "final"
) 

model_cv_crowns <- train(
  x = predictors,
  y = crowns_df$Liana,
  method = "lda",
  trControl = ctrl_crowns
)

cv_lda_crowns_accuracy <-model_cv_crowns$results$Accuracy 
cv_lda_crowns_accuracy 

# Predicted liana proportion
prop.table(table(model_cv_crowns$pred$pred))

## Within MD dataset
blob_df <- as.data.frame(blob_filt_valid)
blob_df$ Liana <- as.factor(label_blob_valid)

set.seed(123)
blob_df[,1:ncol(crowns_filt_valid)]<- scale(blob_df[,1:ncol(crowns_filt_valid)]) #ncol(crowns_filt_valid) should be equal to ncol(blob_filt_valid)
sample_blob <- sample(c(TRUE, FALSE), nrow(blob_df), replace=TRUE, prob=c(0.7,0.3))
train_blob <- blob_df[sample_blob, ]
test_blob <- blob_df[!sample_blob, ] 

model_blob <- lda(Liana ~ ., data = train_blob)
predictions_blob <- predict(model_blob, test_blob)
lda_blob_accuracy <- mean(predictions_blob$class==test_blob$Liana)
lda_blob_accuracy #0.8671278

### Stratified folded cross-validation
predictors <- blob_df[ , -which(names(blob_df) %in% c("Liana"))]

ctrl_blob <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "final"
) 

model_cv_blob <- train(
  x = predictors,
  y = blob_df$Liana,
  method = "lda",
  trControl = ctrl_blob
)

cv_lda_blob_accuracy <-model_cv_blob$results$Accuracy
cv_lda_blob_accuracy 

# Predicted liana proportion
prop.table(table(model_cv_blob$pred$pred))

## Between MD and DP datasets
crowns_df <- as.data.frame(crowns_filt_valid)
crowns_df$Liana <- as.character(label_crowns_valid)
crowns_df <- crowns_df[crowns_df$Liana %in% c('0','3'), ]
crowns_df$Liana[crowns_df$Liana == '3'] <- '1'
crowns_df$Liana <- as.factor(crowns_df$Liana)
type_crowns <- rep(c('crowns'), nrow(crowns_df))
crowns_df$Type <- as.factor(type_crowns)


blob_df <- as.data.frame(blob_filt_valid)
blob_df$ Liana <- as.factor(label_blob_valid)
type_blob <- rep(c('blob'), nrow(blob_df))
blob_df$Type <- as.factor(type_blob)

all_df <- rbind(crowns_df, blob_df)

set.seed(123)
all_df[,1:ncol(crowns_filt_valid)]<- scale(all_df[,1:ncol(crowns_filt_valid)])
# train_all <- all_df[all_df$Type == 'blob', ]
# test_all <- all_df[all_df$Type== 'crowns', ]
test_all <- all_df[all_df$Type == 'blob', ]
train_all <- all_df[all_df$Type== 'crowns', ]

model_all <- lda(Liana ~ . - Type, data = train_all) # remove Type from predictive values
predictions_all <- predict(model_all, test_all)
lda_all_accuracy <- mean(predictions_all$class==test_all$Liana)
lda_all_accuracy 

# Predicted liana proportion
liana_prop_pred <- sum(predictions_all$class== 1) / nrow(test_all)
noliana_prop_pred <- sum(predictions_all$class== 0) / nrow(test_all)
liana_prop_pred
noliana_prop_pred
