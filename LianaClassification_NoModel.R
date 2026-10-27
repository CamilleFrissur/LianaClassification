# Abbreviations used #############################################

# md: manual delineation (blob dataset)
## liana infestation 0 or 1
# dp: drone photos (crowns dataset)
## liana infestation from CR0 to CR3
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
library(miscTools)
library(uwot)
library(FactoMineR)
library(fpc)
library(signal)
library(MASS)
library(caret)
library(gridExtra)
library(grid)


## Set working Directory ####
setwd ("D:/DataCamille/")

# load('Classification_R/LianaClassification/myEnvironment.RData')
#ls()

# Uploading data and cropping ##########################################

## HSI mosaic ####

# r <- rast("HSI/mosaic_2016/nou_20160921_l1c_mosaic.bsq") #2016
r <- rast("HSI/mosaic_2023/mosaic/nou_20231025_l1c_mosaic.bsq") #2023
nlyr(r) #number of spectral bands
# plot(r[[470]])
# hist(r[[470]], breaks = 100)

## Illu ####

# illu <- rast("HSI/mosaic_2016/nou_20160921_ilu_mosaic.bsq")
illu <- rast("HSI/mosaic_2023/mosaic/Nou2023_Vol2_NS_mosaic_ilu.bsq")
illu <- illu$`Illumination Shadows`
nlyr(illu) # should be 1 (Illumination shadows)
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


## Crowns from drone photos (DP) ####
crowns <- vect("Couronnes/Crowns_Liana_2025_GP_TRUE.shp")
crowns <- project(crowns, crs(r))

area_crowns <- crop(r, crowns, mask=TRUE)

label_crowns <- rasterize(crowns, area_crowns, field = "LianaFullC")

## Proportion of lianas before filtering #####################################
table(values(label_crowns), useNA = "no")
prop.table(table(values(label_crowns), useNA = "no"))
table(values(crowns)[["LianaFullC"]])
prop.table(table(values(crowns)[["LianaFullC"]]))

table(values(label_blob), useNA = "no")
prop.table(table(values(label_blob), useNA = "no"))
table(values(blob)[["Type"]])
prop.table(table(values(blob)[["Type"]]))
# Filtering pixels ##########################################################
extract_wavelength <- function(band_names) {
  
  wl <- suppressWarnings(as.numeric(band_names))
  need_regex <- is.na(wl)
  wl[need_regex] <- as.numeric(gsub(".*\\(([0-9.]+)\\).*", "\\1", band_names[need_regex]))
  
  return(wl)
}


## Cleaning extreme values / NA ##
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

  area[area < valid_range[1]] <- valid_range[1]
  area[area > valid_range[2]] <- valid_range[2]
  return(area)
}


# min(area_crowns)
# max(area_crowns)
# sum(is.na(as.matrix(area_crowns)))
area_crowns <- clean_extreme_values_toNA(area_crowns)
# min(area_crowns)
# max(area_crowns)
# sum(is.na(as.matrix(area_crowns)))

# min(area_blob)
# max(area_blob)
# sum(is.na(as.matrix(area_blob)))
area_blob <- clean_extreme_values_toNA(area_blob)
# min(area_blob)
# max(area_blob)
# sum(is.na(as.matrix(area_blob)))

## Removing range & bands with more than 95% na
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

area_crowns <- filter_good_bands(area_crowns, max_na_pct_band = 0.99)
area_blob <- filter_good_bands(area_blob, max_na_pct_band = 0.99)

wl_crowns <- round(extract_wavelength(names(area_crowns)), 2)
wl_blob   <- round(extract_wavelength(names(area_blob)), 2)
common_wl <- intersect(wl_crowns, wl_blob)
common_wl <- sort(common_wl)

area_crowns <- area_crowns[[which(wl_crowns %in% common_wl)]]
area_blob   <- area_blob[[which(wl_blob %in% common_wl)]]


## Removing pixel based on illu ####
# Only used with the 2016 data (illu file missing fo 2023)

med_illu <- global(illu[[1]], fun = median, na.rm = TRUE)$global
perc90_illu <- global(illu[[1]], fun = quantile, probs = 0.9, na.rm = TRUE)$X90.
perc95_illu <- global(illu[[1]], fun = quantile, probs = 0.95, na.rm = TRUE)$X90.


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


## Removing pixel based on NIR reflectance####

# names(area_crowns)
# crowns_nir <- area_crowns[[118:157]]
# mean_crowns <- mean(crowns_nir, na.rm =TRUE)
# hist(mean_crowns)

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


## applying/choosing the filter #########################################################

# Comment and/or uncomment to choose the right filter

### Illu (median) #################
illu_crowns_out <- filter_by_illu(illu, area_crowns, label_crowns, med_illu)
illu_blob_out <- filter_by_illu(illu, area_blob, label_blob, med_illu)

# area_crowns_filt <- illu_crowns_out$area_filt
# label_crowns_filt <- illu_crowns_out$label_filt
# area_blob_filt <- illu_blob_out$area_filt
# label_blob_filt <- illu_blob_out$label_filt
# illu_crop_crowns <- illu_crowns_out$area_crop
# illu_crop_blob <- illu_blob_out$area_crop

### VNIR (mean) ####################
VNIR_crowns_out <- filter_by_vnir(area_crowns, label_crowns)
VNIR_blob_out <- filter_by_vnir(area_blob, label_blob)

# area_crowns_filt <- VNIR_crowns_out$area_filt
# label_crowns_filt <- VNIR_crowns_out$label_filt
# area_blob_filt <- VNIR_blob_out$area_filt
# label_blob_filt <- VNIR_blob_out$label_filt
# illu_crop_crowns <- area_crowns_filt *0
# illu_crop_blob <- area_blob_filt * 0


### NIR (20%) ############################
NIR_crowns_out <- filter_by_nir(area_crowns, label_crowns)
NIR_blob_out <- filter_by_nir(area_blob, label_blob)

area_crowns_filt <- NIR_crowns_out$area_filt
label_crowns_filt <- NIR_crowns_out$label_filt
area_blob_filt <- NIR_blob_out$area_filt
label_blob_filt <- NIR_blob_out$label_filt
illu_crop_crowns <- area_crowns_filt *0
illu_crop_blob <- area_blob_filt * 0

### No filter #########################
# area_crowns_filt <- area_crowns
# label_crowns_filt <- label_crowns
# area_blob_filt <- area_blob
# label_blob_filt <- label_blob
# illu_crop_crowns <- area_crowns *0
# illu_crop_blob <- area_blob * 0

## Normalization #####################
area_crowns_filt   <- area_crowns_filt / mean(area_crowns_filt, na.rm = TRUE)
area_blob_filt <- area_blob_filt / mean(area_blob_filt, na.rm = TRUE)

## Correlation (Chi2) between filters ##########################################

run_chi2 <- function(pixels1, pixels2, name1, name2) {
  tab <- table(pixels1, pixels2)
  test <- chisq.test(tab)
  cat(sprintf("\n--- %s vs %s ---\n", name1, name2))
  cat(sprintf("Chi² = %.4f, df = %d, p-value = %.4f\n",
              test$statistic, test$parameter, test$p.value))
  if (test$p.value < 0.05) {
    cat("→ The filters target different pixels (p < 0.05).\n")
  } else {
    cat("→ No significant difference (p >= 0.05).\n")
  }
}

### Blob dataset
reference_blob <- area_blob

mask_valid_blob <- !is.na(as.matrix(reference_blob[[1]]))  

pixels_illu_blob <- as.matrix(illu_blob_out$label_filt)[mask_valid_blob]
pixels_nir_blob  <- as.matrix(NIR_blob_out$label_filt)[mask_valid_blob]
pixels_vnir_blob <- as.matrix(VNIR_blob_out$label_filt)[mask_valid_blob]

stopifnot(length(pixels_illu_blob) == length(pixels_nir_blob), length(pixels_illu_blob) == length(pixels_vnir_blob))

run_chi2(pixels_illu_blob, pixels_nir_blob, "illu", "NIR")
run_chi2(pixels_nir_blob, pixels_vnir_blob, "NIR", "VNIR")
run_chi2(pixels_vnir_blob, pixels_illu_blob, "VNIR", "illu")

### Crowns dataset
reference_crowns <- area_crowns

mask_valid_crowns <- !is.na(as.matrix(reference_crowns[[1]]))  

pixels_illu_crowns <- as.matrix(illu_crowns_out$label_filt)[mask_valid_crowns]
pixels_nir_crowns  <- as.matrix(NIR_crowns_out$label_filt)[mask_valid_crowns]
pixels_vnir_crowns <- as.matrix(VNIR_crowns_out$label_filt)[mask_valid_crowns]

stopifnot(length(pixels_illu_crowns) == length(pixels_nir_crowns), length(pixels_illu_crowns) == length(pixels_vnir_crowns))

run_chi2(pixels_illu_crowns, pixels_nir_crowns, "illu", "NIR")
run_chi2(pixels_nir_crowns, pixels_vnir_crowns, "NIR", "VNIR")
run_chi2(pixels_vnir_crowns, pixels_illu_crowns, "VNIR", "illu")


# Dimensionality reduction and data visualization ###############################

## Data prep ##############################################################

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

## Proportion of lianas after data_prep #####################################
table(label_crowns_valid)
prop.table(table(label_crowns_valid))

table(label_blob_valid)
prop.table(table(label_blob_valid))

table(label_crowns_poly_valid)
prop.table(table(label_crowns_poly_valid))

table(label_blob_poly_valid)
prop.table(table(label_blob_poly_valid))

## PCA #######################################################################

## DP per pixel
pca_out_crowns <- PCA(crowns_filt_valid, scale.unit = TRUE, graph = FALSE)
var_exp_crowns <- pca_out_crowns$eig[1:2, "percentage of variance"]

pca_plot_crowns <- data.frame(
  x = pca_out_crowns$ind$coord[,1],
  y = pca_out_crowns$ind$coord[,2],
  LianaFullCrown = as.factor(label_crowns_valid)
)

ggplot(pca_plot_crowns, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("0"= "pink1", "1"="palevioletred1","2"="palevioletred3","3"= "palevioletred4"),
  #scale_color_manual(values = c("0"= "#662700", "1"="#CC5500","2"="#FF8E33","3"= "#FFCA99"),
  #scale_color_manual(values = c("0"= "#21421e", "1"="#568203","2"="#9acd32","3"= "#c5e384"),
                     
                     labels = c("0" = "CR0", "1" = "CR1", "2"= "CR2", "3"= "CR3"))+
  labs(title= "PCA of the HSI values of the drone photos dataset (Illu filter)",
       x = paste0("Dim 1 (", round(var_exp_crowns[1], 1), "%)"),
       y = paste0("Dim 2 (", round(var_exp_crowns[2], 1), "%)"))+
  theme_minimal()

# Only useful if the illumination filter was used
# pca_plot_crowns$illu <- illu_crowns_valid
# ggplot(pca_plot_crowns, aes(x = x, y = y, color = illu)) +
#   geom_point() +
#   scale_color_viridis_c(name = "Illumination") +
#   labs(title = "PCA of the HSI values of the DP dataset colored by illumination")

## For DP per polygon 
pca_out_crowns_poly <- PCA(crowns_poly_valid, scale.unit = TRUE, graph = FALSE)
var_exp_crowns_poly <- pca_out_crowns_poly$eig[1:2, "percentage of variance"]

pca_plot_crowns_poly <- data.frame(
  x = pca_out_crowns_poly$ind$coord[,1],
  y = pca_out_crowns_poly$ind$coord[,2],
  LianaFullCrown = as.factor(label_crowns_poly_valid)
)

ggplot(pca_plot_crowns_poly, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("CR0"= "pink1", "CR1"="palevioletred1","CR2"="palevioletred3","CR3"= "palevioletred4"))+
  labs(title= "PCA of the mean spectrum per polygon of the DP dataset (Illu filter)",
       x = paste0("Dim 1 (", round(var_exp_crowns_poly[1], 1), "%)"),
       y = paste0("Dim 2 (", round(var_exp_crowns_poly[2], 1), "%)"))+
  theme_minimal()

## For MD per pixel
pca_out_blob <- PCA(blob_filt_valid, scale.unit = TRUE, graph = FALSE)
var_exp_blob <- pca_out_blob$eig[1:2, "percentage of variance"]

pca_plot_blob <- data.frame(
  x = pca_out_blob$ind$coord[,1],
  y = pca_out_blob$ind$coord[,2],
  Type = as.factor(label_blob_valid)
)

ggplot(pca_plot_blob, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#ABABAB", "1" = "#001312"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "PCA of the HSI values of the MD dataset (Illu filter)",
       x = paste0("Dim 1 (", round(var_exp_blob[1], 1), "%)"),
       y = paste0("Dim 2 (", round(var_exp_blob[2], 1), "%)"))+
  theme_minimal()

# Only with illu filter
# pca_plot_blob$illu <- illu_blob_valid
# ggplot(pca_plot_blob, aes(x = x, y = y, color = illu)) +
#   geom_point() +
#   scale_color_viridis_c(name = "Illumination") +
#   labs(title = "PCA of the HSI values of the MD dataset colored by illumination")

## For MD per polygon 
pca_out_blob_poly <- PCA(blob_poly_valid, scale.unit = TRUE, graph = FALSE)
var_exp_blob_poly <- pca_out_blob_poly$eig[1:2, "percentage of variance"]

pca_plot_blob_poly <- data.frame(
  x = pca_out_blob_poly$ind$coord[,1],
  y = pca_out_blob_poly$ind$coord[,2],
  Type = as.factor(label_blob_poly_valid)
)

ggplot(pca_plot_blob_poly, aes(x = x, y = y, color = Type)) +
  geom_point() +
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#ABABAB", "1" = "#001312"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "PCA of the mean spectrum per polygon of the MD dataset (Illu filter)",
       x = paste0("Dim 1 (", round(var_exp_blob_poly[1], 1), "%)"),
       y = paste0("Dim 2 (", round(var_exp_blob_poly[2], 1), "%)"))+
  theme_minimal()


## t-SNE #######################################################################

## DP per pixel
perplexity_crown = (nrow(crowns_filt_valid)-1)/3
perplexity_crown
set.seed(2)
tsne_out_crowns <- Rtsne(crowns_filt_valid, check_duplicates = FALSE, perplexity = perplexity_crown )

tsne_plot_crowns <- data.frame(
  x = tsne_out_crowns$Y[,1],
  y = tsne_out_crowns$Y[,2],
  LianaFullCrown = as.factor(label_crowns_valid)
)


ggplot(tsne_plot_crowns, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("0"= "pink1", "1"="palevioletred1","2"="palevioletred3","3"= "palevioletred4"),
                     labels = c("0" = "CR0", "1" = "CR1", "2"= "CR2", "3"= "CR3"))+
  labs(title= "tSNE of the HSI values of the drone photos dataset (Illu filter)" )+
  theme_minimal()

# Only useful if the illumination filter was used
# tsne_plot_crowns$illu <- illu_crowns_valid
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
  scale_color_manual(values = c("CR0"= "pink1", "CR1"="palevioletred1","CR2"="palevioletred3","CR3"= "palevioletred4"))+
  # scale_color_viridis_d(option="turbo")+
  labs(title= "tSNE of the mean spectrum per polygon of the DP dataset (Illu filter)" )+
  theme_minimal()

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

ggplot(tsne_plot_blob, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#ABABAB", "1" = "#001312"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "tSNE of the HSI values of the MD dataset (Illu filter)" )+
  theme_minimal()

# Only with illu filter
# tsne_plot_blob$illu <- illu_blob_valid
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
    values = c("0" = "#ABABAB", "1" = "#001312"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "tSNE of the mean spectrum per polygon of the MD dataset (Illu filter)" )+
  theme_minimal()


## UMAP #########################################################################

#For DP per pixel

umap_out_crowns <- umap2( seed = 2, crowns_filt_valid,  n_neighbors = min(15, nrow(crowns_filt_valid) - 1))

umap_plot_crowns <- data.frame(
  x = umap_out_crowns[,1],
  y = umap_out_crowns[,2],
  LianaFullCrown = as.factor(label_crowns_valid)
)

  
  
ggplot(umap_plot_crowns, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("0"= "pink1", "1"="palevioletred1","2"="palevioletred3","3"= "palevioletred4"),
                     labels = c("0" = "CR0", "1" = "CR1", "2"= "CR2", "3"= "CR3"))+
  # scale_color_viridis_d(option="turbo")+
  labs(title= "UMAP of the HSI values of the DP dataset (Illu filter)" )+
  theme_minimal()

# Only with illu filter
# umap_plot_crowns$illu <- illu_crowns_valid
# ggplot(umap_plot_crowns, aes(x = x, y = y, color = illu)) +
  # geom_point() +
  # scale_color_viridis_c(name = "Illumination") +
  # labs(title = "UMAP of the HSI values of the DP dataset colored by illumination")

#For DP per polygon
umap_out_crowns_poly <- umap2( seed = 2, crowns_poly_valid)

umap_plot_crowns_poly <- data.frame(
  x = umap_out_crowns_poly[,1],
  y = umap_out_crowns_poly[,2],
  LianaFullCrown = as.factor(label_crowns_poly_valid))

ggplot(umap_plot_crowns_poly, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("CR0"= "pink1", "CR1"="palevioletred1","CR2"="palevioletred3","CR3"= "palevioletred4"))+
  labs(title= "UMAP of the mean HSI value per polygon of the DP dataset (Illu filter)" )+
  theme_minimal()


#For MD per pixel
umap_out_blob <- umap2( seed = 2, blob_filt_valid)

umap_plot_blob <- data.frame(
  x = umap_out_blob[,1],
  y = umap_out_blob[,2],
  Type = as.factor(label_blob_valid)
)


ggplot(umap_plot_blob, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#ABABAB", "1" = "#001312"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "UMAP of the HSI values of the MD dataset (Illu filter)" )+
  theme_minimal()

# Only with illu filter
# umap_plot_blob$illu <- illu_blob_valid
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
    values = c("0" = "#ABABAB", "1" = "#001312"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "UMAP of the mean HSI value per polygon of the MD dataset (Illu filter)" )+
  theme_minimal()

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
#   #   palevioletred1norm(group0[, i], main = paste("Composante", i, "- groupe 0"))
#   #   palevioletred1line(group0[, i])
#   # }
#   #
#   # par(mfrow = c(2, 3))
#   # for (i in 1:n_components) {
#   #   palevioletred1norm(group1[, i], main = paste("Composante", i, "- groupe 1"))
#   #   palevioletred1line(group1[, i])
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

sensor_to_keep <- function(matrix_valid, keep_sensor_in_overlap = c("swir", "vnir")) {
  
  keep_sensor_in_overlap <- match.arg(keep_sensor_in_overlap)
  wl <- extract_wavelength(colnames(matrix_valid))
  
  break_idx <- which(diff(wl) < 0)[1]       
  
  if (is.na(break_idx)) return(matrix_valid)
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
  

  matrix_valid <- matrix_valid[,which(keep_overlap), drop = FALSE]
}

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
  
crowns_poly_valid <- sensor_to_keep(crowns_poly_valid, "swir")
blob_poly_valid <- sensor_to_keep(blob_poly_valid, "swir")


spec_crowns <- mean_spectrum_by_group(crowns_poly_valid, label_crowns_poly_valid, "DP")
spec_blob   <- mean_spectrum_by_group(blob_poly_valid, label_blob_poly_valid, "MD")

all_spectra <- rbind(spec_crowns, spec_blob)

#selected <- subset(all_spectra, group %in% c("DP_CR0", "DP_CR3", "MD_0", "MD_1") & all_spectra$wavelength<900) #VNIR only
selected <- subset(all_spectra, group %in% c("DP_CR0", "DP_CR3", "MD_0", "MD_1"))



label_config <- data.frame(
  group = c("DP_CR0", "DP_CR1", "DP_CR2", "DP_CR3", "MD_0", "MD_1"),
  display_name = c("CR0 ( No Liana DP)", "CR1 (DP)", "CR2", "CR3 (Liana DP)", "MD0 (No Liana MD)", "MD1 ( Liana MD)"),
  color = c("pink1", "palevioletred1", "palevioletred3","palevioletred4", "#ABABAB", "#001312"),
  #color = c("#662700", "#CC5500", "#FF8E33","#FFCA99", "#ABABAB", "#001312"),
  #color = c("#21421e", "#568203", "#9acd32","#c5e384", "#ABABAB", "#001312"),
  
  linetype = c ("solid","solid","solid","solid", "solid", "solid")
  #linetype = c ("solid","solid","solid","solid", "solid", "solid")
  
)

add_segments <- function(dt, gap_nm = 8) {
  dt <- as.data.table(dt)
  dt[, segment := cumsum(c(TRUE, diff(wavelength) <= 0 | diff(wavelength) > gap_nm)), by = group]
  dt
}

selected <- add_segments(selected)

selected_smooth <- as.data.frame(selected)
selected_smooth$mean_smooth <- selected_smooth$mean
selected_smooth$sd_smooth   <- selected_smooth$sd

for (g in unique(selected_smooth$group)) {
  for (s in unique(selected_smooth$segment[selected_smooth$group == g])) {
    
    rows <- which(selected_smooth$group == g & selected_smooth$segment == s)
    
    if (length(rows) >= 11) {
      selected_smooth$mean_smooth[rows] <- sgolayfilt(selected_smooth$mean[rows], p = 3, n = 11)
      selected_smooth$sd_smooth[rows]   <- sgolayfilt(selected_smooth$sd[rows], p = 3, n = 11)
    }
  }
}

ggplot(selected_smooth, aes(x = wavelength)) +
  geom_ribbon(
    aes(
      ymin = mean_smooth - sd_smooth,
      ymax = mean_smooth + sd_smooth,
      group = interaction(group, segment),
      fill = group,
      colour = group
    ),
    alpha = 0.2,
    linetype = 2
  ) +
  geom_line(
    aes(
      y = mean_smooth,
      color = group,
      linetype = group,
      group = interaction(group, segment)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(label_config$color, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_linetype_manual(
    values = setNames(label_config$linetype, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_fill_manual(
    values = setNames(label_config$color, label_config$group),  
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "Mean normalized reflectance",
    title = "Spectrum by dataset group (NIR filter)",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()


ggplot(selected_smooth[selected_smooth$wavelength < 690,], aes(x = wavelength)) +
  geom_ribbon(
    aes(
      ymin = mean_smooth - sd_smooth,
      ymax = mean_smooth + sd_smooth,
      group = interaction(group, segment),
      fill = group,
      colour = group),
    alpha = 0.2,
    linetype = 2 
  ) +
  geom_line(
    aes(
      y = mean_smooth,
      color = group,
      linetype = group,
      group = interaction(group, segment)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(label_config$color, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_linetype_manual(
    values = setNames(label_config$linetype, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_fill_manual(
    values = setNames(label_config$color, label_config$group),  
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "Mean normalized reflectance",
    title = "Spectrum by dataset group in the visible region (NIR filter)",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()

ggplot(selected_smooth[selected_smooth$wavelength >700 & selected_smooth$wavelength < 1340,], aes(x = wavelength)) +
  geom_ribbon(
    aes(
      ymin = mean_smooth - sd_smooth,
      ymax = mean_smooth + sd_smooth,
      group = interaction(group, segment),
      fill = group ,
      colour = group),
    alpha = 0.2,
    linetype = 2 
  ) +
  geom_line(
    aes(
      y = mean_smooth,
      color = group,
      linetype = group,
      group = interaction(group, segment)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(label_config$color, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_linetype_manual(
    values = setNames(label_config$linetype, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_fill_manual(
    values = setNames(label_config$color, label_config$group),  
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "Mean normalized reflectance",
    title = "Spectrum by dataset group in the NIR region (NIR filter)",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()

ggplot(selected_smooth[selected_smooth$wavelength >1420 & selected_smooth$wavelength < 1790,], aes(x = wavelength)) +
  geom_ribbon(
    aes(
      ymin = mean_smooth - sd_smooth,
      ymax = mean_smooth + sd_smooth,
      group = interaction(group, segment),
      fill = group,
      colour = group),
    alpha = 0.2,
    linetype = 2 
  ) +
  geom_line(
    aes(
      y = mean_smooth,
      color = group,
      linetype = group,
      group = interaction(group, segment)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(label_config$color, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_linetype_manual(
    values = setNames(label_config$linetype, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_fill_manual(
    values = setNames(label_config$color, label_config$group),  
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "Mean normalized reflectance",
    title = "Spectrum by dataset group in the SWIR 1 region (NIR filter)",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()

## Difference in mean spectrum ####

cr0 <- subset(selected, group == "DP_CR0")
cr3 <- subset(selected, group == "DP_CR3")
md0 <- subset(selected, group == "MD_0")
md1 <- subset(selected, group == "MD_1")

diff_cr <- merge(cr0, cr3, by = "wavelength", suffixes = c("_cr0", "_cr3"))
diff_cr$mean <- diff_cr$mean_cr3 - diff_cr$mean_cr0
diff_cr$sd <- sqrt(diff_cr$sd_cr0^2 + diff_cr$sd_cr3^2) 
diff_cr$group <- "CR3 - CR0"

diff_md <- merge(md0, md1, by = "wavelength", suffixes = c("_md0", "_md1"))
diff_md$mean <- diff_md$mean_md1 - diff_md$mean_md0
diff_md$sd <- sqrt(diff_md$sd_md0^2 + diff_md$sd_md1^2)
diff_md$group <- "MD1 - MD0"

diff_spectra <- rbind(
  diff_cr[, c("wavelength", "mean", "sd", "group")],
  diff_md[, c("wavelength", "mean", "sd", "group")]
)

diff_spectra <- add_segments(diff_spectra, gap_nm = 8)

diff_spectra_smooth <- as.data.frame(diff_spectra)
diff_spectra_smooth$mean_smooth <- diff_spectra_smooth$mean
diff_spectra_smooth$sd_smooth <- diff_spectra_smooth$sd

for (g in unique(diff_spectra_smooth$group)) {
  for (s in unique(diff_spectra_smooth$segment[diff_spectra_smooth$group == g])) {
    rows <- which(diff_spectra_smooth$group == g & diff_spectra_smooth$segment == s)
    if (length(rows) >= 11) {
      diff_spectra_smooth$mean_smooth[rows] <- sgolayfilt(diff_spectra_smooth$mean[rows], p = 3, n = 11)
      diff_spectra_smooth$sd_smooth[rows]   <- sgolayfilt(diff_spectra_smooth$sd[rows], p = 3, n = 11)
    }
  }
}

diff_label_config <- data.frame(
  group = c("CR3 - CR0", "MD1 - MD0"),
  display_name = c("CR3 - CR0", "MD1 - MD0"),
  color = c("palevioletred1", "#001312"), 
  linetype = c("solid", "solid")
)

ggplot(diff_spectra_smooth, aes(x = wavelength)) +
  geom_ribbon(
    aes(
      ymin = mean_smooth - sd_smooth,
      ymax = mean_smooth + sd_smooth,
      group = interaction(group, segment),
      fill = group,
      colour = group
    ),
    alpha = 0.2,
    linetype = 2
  ) +
  geom_line(
    aes(
      y = mean_smooth,
      color = group,
      linetype = group,
      group = interaction(group, segment)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(diff_label_config$color, diff_label_config$group),
    labels = setNames(diff_label_config$display_name, diff_label_config$group)
  ) +
  scale_linetype_manual(
    values = setNames(diff_label_config$linetype, diff_label_config$group),
    labels = setNames(diff_label_config$display_name, diff_label_config$group)
  ) +
  scale_fill_manual(
    values = setNames(diff_label_config$color, diff_label_config$group),
    labels = setNames(diff_label_config$display_name, diff_label_config$group)
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "Difference in normalized reflectance",
    title = "Spectrum difference by group (CR3-CR0 and MD1-MD0, NIR filter)",
    color = "Difference (Liana - No liana)",
    linetype = "Difference (Liana - No liana)",
    fill = "Difference (Liana - No liana)"
  ) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5)  

# Visible region plot
ggplot(diff_spectra_smooth[diff_spectra_smooth$wavelength < 690, ], aes(x = wavelength)) +
  geom_ribbon(
    aes(
      ymin = mean_smooth - sd_smooth,
      ymax = mean_smooth + sd_smooth,
      group = interaction(group, segment),
      fill = group,
      colour = group
    ),
    alpha = 0.2,
    linetype = 2
  ) +
  geom_line(
    aes(
      y = mean_smooth,
      color = group,
      linetype = group,
      group = interaction(group, segment)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(diff_label_config$color, diff_label_config$group),
    labels = setNames(diff_label_config$display_name, diff_label_config$group)
  ) +
  scale_fill_manual(
    values = setNames(diff_label_config$color, diff_label_config$group),
    labels = setNames(diff_label_config$display_name, diff_label_config$group)
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "Difference in normalized reflectance",
    title = "Spectrum difference in the visible region (NIR filter)",
    color = "Difference (Liana - No Liana)",
    fill = "Difference (Liana - No Liana)"
  ) +
  theme_minimal() +
  guides(
    color = guide_legend(title = "Difference"),
    fill = "none",
    linetype = "none"
  )+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5)

# NIR plot
ggplot(diff_spectra_smooth[diff_spectra_smooth$wavelength > 700 & diff_spectra_smooth$wavelength < 1340, ], aes(x = wavelength)) +
  geom_ribbon(
    aes(
      ymin = mean_smooth - sd_smooth,
      ymax = mean_smooth + sd_smooth,
      group = interaction(group, segment),
      fill = group,
      colour = group
    ),
    alpha = 0.2,
    linetype = 2
  ) +
  geom_line(
    aes(
      y = mean_smooth,
      color = group,
      linetype = group,
      group = interaction(group, segment)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(diff_label_config$color, diff_label_config$group),
    labels = setNames(diff_label_config$display_name, diff_label_config$group)
  ) +
  scale_fill_manual(
    values = setNames(diff_label_config$color, diff_label_config$group),
    labels = setNames(diff_label_config$display_name, diff_label_config$group)
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "Difference in normalized reflectance",
    title = "Spectrum difference in the NIR region (NIR filter)",
    color = "Difference (Liana - No Liana)",
    fill = "Difference (Liana - No Liana)"
  ) +
  theme_minimal()+
  guides(
    color = guide_legend(title = "Difference"),
    fill = "none",
    linetype = "none"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5)

# SWIR1 plot
ggplot(diff_spectra_smooth[diff_spectra_smooth$wavelength > 1420 & diff_spectra_smooth$wavelength < 1790, ], aes(x = wavelength)) +
  geom_ribbon(
    aes(
      ymin = mean_smooth - sd_smooth,
      ymax = mean_smooth + sd_smooth,
      group = interaction(group, segment),
      fill = group,
      colour = group
    ),
    alpha = 0.2,
    linetype = 2
  ) +
  geom_line(
    aes(
      y = mean_smooth,
      color = group,
      linetype = group,
      group = interaction(group, segment)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(diff_label_config$color, diff_label_config$group),
    labels = setNames(diff_label_config$display_name, diff_label_config$group)
  ) +
  scale_fill_manual(
    values = setNames(diff_label_config$color, diff_label_config$group),
    labels = setNames(diff_label_config$display_name, diff_label_config$group)
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "Difference in normalized reflectance",
    title = "Spectrum difference in the SWIR 1 region (NIR filter)",
    color = "Difference (Liana - No Liana)",
    fill = "Difference (Liana - No Liana)"
  ) +
  theme_minimal() +
  guides(
    color = guide_legend(title = "Difference"),
    fill = "none",
    linetype = "none"
  )+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5)

# Median spectrum comparison ###################################################

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

spec_crowns <- med_spectrum_by_group(crowns_poly_valid, label_crowns_poly_valid, "DP")
spec_blob   <- med_spectrum_by_group(blob_poly_valid, label_blob_poly_valid, "MD")

all_spectra <- rbind(spec_crowns, spec_blob)

#selected <- subset(all_spectra, group %in% c("DP_CR0", "DP_CR3", "MD_0", "MD_1") & all_spectra$wavelength<900) #VNIR only
selected <- subset(all_spectra, group %in% c("DP_CR0", "DP_CR3", "MD_0", "MD_1"))



label_config <- data.frame(
  group = c("DP_CR0", "DP_CR1", "DP_CR2", "DP_CR3", "MD_0", "MD_1"),
  display_name = c("CR0 ( No Liana DP)", "CR1 (DP)", "CR2", "CR3 (Liana DP)", "No liana (MD)", "Liana (MD)"),
  color = c("pink1", "palevioletred1", "palevioletred3","palevioletred4", "#ABABAB", "#001312"),
  #color = c("#662700", "#CC5500", "#FF8E33","#FFCA99", "#ABABAB", "#001312"),
  #color = c("#21421e", "#568203", "#9acd32","#c5e384", "#ABABAB", "#001312"),
  
  linetype = c ("solid","solid","solid","solid", "solid", "solid")
  #linetype = c ("solid","solid","solid","solid", "solid", "solid")
  
)

add_segments <- function(dt, gap_nm = 8) {
  dt <- as.data.table(dt)
  dt[, segment := cumsum(c(TRUE, diff(wavelength) <= 0 | diff(wavelength) > gap_nm)), by = group]
  dt
}

selected <- add_segments(selected)

selected_smooth <- as.data.frame(selected)
selected_smooth$med_smooth <- selected_smooth$med
selected_smooth$sd_smooth   <- selected_smooth$sd

for (g in unique(selected_smooth$group)) {
  for (s in unique(selected_smooth$segment[selected_smooth$group == g])) {
    
    rows <- which(selected_smooth$group == g & selected_smooth$segment == s)
    
    if (length(rows) >= 11) {
      selected_smooth$med_smooth[rows] <- sgolayfilt(selected_smooth$med[rows], p = 3, n = 11)
      selected_smooth$sd_smooth[rows]   <- sgolayfilt(selected_smooth$sd[rows], p = 3, n = 11)
    }
  }
}

ggplot(selected_smooth, aes(x = wavelength)) +
  geom_ribbon(
    aes(
      ymin = med_smooth - sd_smooth,
      ymax = med_smooth + sd_smooth,
      group = interaction(group, segment),
      fill = group,
      colour = group),
    alpha = 0.2,
    linetype = 2 
  ) +
  geom_line(
    aes(
      y = med_smooth,
      color = group,
      linetype = group,
      group = interaction(group, segment)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(label_config$color, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_linetype_manual(
    values = setNames(label_config$linetype, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_fill_manual(
    values = setNames(label_config$color, label_config$group),  
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "med normalized reflectance",
    title = "Spectrum by dataset group (NIR filter)",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()

ggplot(selected_smooth, aes(x = wavelength)) +
  geom_ribbon(
    aes(
      ymin = med_smooth - sd_smooth,
      ymax = med_smooth + sd_smooth,
      group = interaction(group, segment),
      fill = group,
      colour = group),
    alpha = 0.2,
    linetype = 2 
  ) +
  geom_line(
    aes(
      y = med_smooth,
      color = group,
      linetype = group,
      group = interaction(group, segment)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(label_config$color, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_linetype_manual(
    values = setNames(label_config$linetype, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_fill_manual(
    values = setNames(label_config$color, label_config$group),  
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "med normalized reflectance",
    title = "Spectrum by dataset group (NIR filter)",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()

ggplot(selected_smooth, aes(x = wavelength)) +
  geom_ribbon(
    aes(
      ymin = med_smooth - sd_smooth,
      ymax = med_smooth + sd_smooth,
      group = interaction(group, segment),
      fill = group,
      colour = group),
    alpha = 0.2,
    linetype = 2 
  ) +
  geom_line(
    aes(
      y = med_smooth,
      color = group,
      linetype = group,
      group = interaction(group, segment)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(label_config$color, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_linetype_manual(
    values = setNames(label_config$linetype, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_fill_manual(
    values = setNames(label_config$color, label_config$group),  
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "med normalized reflectance",
    title = "Spectrum by dataset group (NIR filter)",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()

ggplot(selected_smooth[selected_smooth$wavelength < 690,], aes(x = wavelength)) +
  geom_ribbon(
    aes(
      ymin = med_smooth - sd_smooth,
      ymax = med_smooth + sd_smooth,
      group = interaction(group, segment),
      fill = group,
      colour = group),
    alpha = 0.2,
    linetype = 2 
  ) +
  geom_line(
    aes(
      y = med_smooth,
      color = group,
      linetype = group,
      group = interaction(group, segment)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(label_config$color, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_linetype_manual(
    values = setNames(label_config$linetype, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_fill_manual(
    values = setNames(label_config$color, label_config$group),  
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "med normalized reflectance",
    title = "Spectrum by dataset group in the visible region (NIR filter)",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()

ggplot(selected_smooth[selected_smooth$wavelength >700 & selected_smooth$wavelength < 1340,], aes(x = wavelength)) +
  geom_ribbon(
    aes(
      ymin = med_smooth - sd_smooth,
      ymax = med_smooth + sd_smooth,
      group = interaction(group, segment),
      fill = group,
      colour = group),
    alpha = 0.2,
    linetype = 2 
  ) +
  geom_line(
    aes(
      y = med_smooth,
      color = group,
      linetype = group,
      group = interaction(group, segment)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(label_config$color, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_linetype_manual(
    values = setNames(label_config$linetype, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_fill_manual(
    values = setNames(label_config$color, label_config$group),  
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "med normalized reflectance",
    title = "Spectrum by dataset group in the NIR region (NIR filter)",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()

ggplot(selected_smooth[selected_smooth$wavelength >1420 & selected_smooth$wavelength < 1790,], aes(x = wavelength)) +
  geom_ribbon(
    aes(
      ymin = med_smooth - sd_smooth,
      ymax = med_smooth + sd_smooth,
      group = interaction(group, segment),
      fill = group,
      colour = group),
    alpha = 0.2,
    linetype = 2 
  ) +
  geom_line(
    aes(
      y = med_smooth,
      color = group,
      linetype = group,
      group = interaction(group, segment)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = setNames(label_config$color, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_linetype_manual(
    values = setNames(label_config$linetype, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  scale_fill_manual(
    values = setNames(label_config$color, label_config$group),  
    labels = setNames(label_config$display_name, label_config$group)
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "med normalized reflectance",
    title = "Spectrum by dataset group in the SWIR 1 region (NIR filter)",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()

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

n_classes <- length(levels(train_crowns$Liana))
model_crowns <- lda(Liana ~ ., data = train_crowns, prior = rep(1/n_classes, n_classes))
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

length(levels(crowns_df$Liana))
model_cv_crowns <- train(
  x = predictors,
  y = crowns_df$Liana,
  method = "lda",
  trControl = ctrl_crowns, 
  prior = rep(1/n_classes, n_classes)
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

n_classes <- length(levels(train_blob$Liana))
model_blob <- lda(Liana ~ ., data = train_blob, prior = rep(1/n_classes, n_classes))
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

n_classes <-  length(levels(blob_df$Liana))
model_cv_blob <- train(
  x = predictors,
  y = blob_df$Liana,
  method = "lda",
  trControl = ctrl_blob,
  prior = rep(1/n_classes, n_classes)
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
train_all <- all_df[all_df$Type == 'blob', ]
test_all <- all_df[all_df$Type== 'crowns', ]
# test_all <- all_df[all_df$Type == 'blob', ]
# train_all <- all_df[all_df$Type== 'crowns', ]

n_classes <- length(levels(train_all$Liana))
model_all <- lda(Liana ~ . - Type, data = train_all, prior = rep(1/n_classes, n_classes)) # remove Type from predictive values
predictions_all <- predict(model_all, test_all)
lda_all_accuracy <- mean(predictions_all$class==test_all$Liana)
lda_all_accuracy 

# Predicted liana proportion
liana_prop_pred <- sum(predictions_all$class== 1) / nrow(test_all)
noliana_prop_pred <- sum(predictions_all$class== 0) / nrow(test_all)
liana_prop_pred
noliana_prop_pred

## LDA trained on pure labels to check % of liana predicted ##################
crowns_df_pure <- as.data.frame(crowns_filt_valid)
crowns_df_pure$LianaOriginal <- as.factor(label_crowns_valid)

crowns_df_pure$Liana <- as.character(label_crowns_valid)
crowns_df_pure$Liana[crowns_df_pure$Liana %in% c("0", "1")] <- "0"
crowns_df_pure$Liana[crowns_df_pure$Liana %in% c("2", "3")] <- "1"
crowns_df_pure$Liana <- as.factor(crowns_df_pure$Liana)

crowns_df_pure$Type <- factor("crowns", levels = c("crowns", "blob"))
crowns$poly_id <- seq_len(nrow(crowns))
poly_id_r <- rasterize(crowns, area_crowns_filt, field = "poly_id")
poly_id_m <- as.matrix(poly_id_r)[, 1]
poly_id_valid <- poly_id_m[DP_per_pixel$valid_idx]
crowns_df_pure$poly_id <- ifelse(is.na(poly_id_valid), 0, poly_id_valid)

blob_df_pure <- as.data.frame(blob_filt_valid)
blob_df_pure$Liana <- as.factor(label_blob_valid)
blob_df_pure$Type <- factor("blob", levels = c("crowns", "blob"))
blob_df_pure$LianaOriginal <- NA
blob_df_pure$poly_id <- 0

common_cols <- intersect(names(crowns_df_pure), names(blob_df_pure))
crowns_df_pure <- crowns_df_pure[, common_cols]
blob_df_pure <- blob_df_pure[, common_cols]

all_df <- rbind(crowns_df_pure, blob_df_pure)

### Trained on blob / Tested on crowns ####
set.seed(123)
band_cols <- sapply(all_df, is.numeric)
band_cols <- names(all_df)[band_cols]

all_df[, band_cols] <- scale(all_df[, band_cols])

train_all_pure <- all_df[all_df$Type == "blob", ]
test_all_pure <- all_df[all_df$Type == "crowns", ]

X_train <- train_all_pure[, band_cols]
constant_cols <- sapply(X_train, function(x) length(unique(x)) == 1)
X_train_filtered <- X_train[, !constant_cols, drop = FALSE]
band_cols_filtered <- names(X_train_filtered)

if (ncol(X_train_filtered) > 1) {
  cor_matrix <- cor(X_train_filtered, use = "complete.obs")
  high_cor_cols <- names(X_train_filtered)[
    sapply(1:ncol(X_train_filtered), function(i) any(abs(cor_matrix[i, ]) > 0.99 & cor_matrix[i, ] != 1))
  ]
  band_cols_filtered <- setdiff(band_cols_filtered, high_cor_cols)
  X_train_filtered <- X_train_filtered[, band_cols_filtered, drop = FALSE]
}

train_data_blob <- data.frame(Liana = train_all_pure$Liana, X_train_filtered)
n_classes <- length(levels(train_data_blob$Liana))
model_blob_to_crown <- lda(Liana ~ ., data = train_data_blob, tol = 1e-8)

X_test <- test_all_pure[, band_cols]
X_test_filtered <- X_test[, !constant_cols, drop = FALSE]
X_test_final <- X_test_filtered[, band_cols_filtered, drop = FALSE]
test_data_crown <- data.frame(Liana = test_all_pure$Liana, X_test_final)

predictions_blob_to_crown <- predict(model_blob_to_crown, test_data_crown)
test_all_pure$predicted_blob <- predictions_blob_to_crown$class
test_all_pure$correct_blob <- test_all_pure$predicted_blob == test_all_pure$Liana

accuracy_blob_to_crown <- mean(test_all_pure$correct_blob)
cat("Mean accuracy per crown:", round(accuracy_blob_to_crown, 3), "\n")

accuracy_per_crown_blob <- aggregate(correct_blob ~ poly_id + LianaOriginal,
                                     data = test_all_pure,
                                     FUN = mean)

names(accuracy_per_crown_blob)[3] <- "accuracy"

p1 <- ggplot(accuracy_per_crown_blob, aes(x = LianaOriginal, y = accuracy, fill = LianaOriginal)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.5) +
  scale_fill_manual(values = c("0" = "pink1", "1" = "palevioletred1",
                               "2" = "palevioletred3", "3" = "palevioletred4")) +
  labs(x = "Liana class in drone photos (DP)", y = "Mean accuracy per crown",
       title = "LDA trained on MD, tested on DP (Illu filter)") +
  theme(legend.position = "none")
p1

test_all_pure$predicted_blob_num <- as.numeric(as.character(test_all_pure$predicted_blob))

predicted_liana_per_crown_blob <- aggregate(
  predicted_blob_num ~ poly_id + LianaOriginal,
  data = test_all_pure,
  FUN = mean
)
names(predicted_liana_per_crown_blob)[3] <- "prediction"


predicted_liana_per_crown_blob$prediction <- predicted_liana_per_crown_blob$prediction * 100

p2 <- ggplot(predicted_liana_per_crown_blob, aes(x = LianaOriginal, y = prediction, fill = LianaOriginal)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.5) +
  scale_fill_manual(values = c("0" = "pink1", "1" = "palevioletred1",
                               "2" = "palevioletred3", "3" = "palevioletred4")) +
  labs(x = "Liana class in drone photos (DP)", y = "Liana predicted per crown (%)",
       title = "LDA trained on MD, tested on DP (Illu filter)") +
  theme(legend.position = "none")

p2

grid.arrange(p1, p2, ncol = 2)

### StratifiedGroupKFold ####

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

#### For crowns ####
set.seed(123)
folds_crowns <- stratified_group_kfold(crowns_df_pure, "poly_id", "Liana", k = 5)

cv_results <- list()
for (i in 1:length(folds_crowns)) {  
  test_indices <- folds_crowns[[i]]  
  train_indices <- setdiff(1:nrow(crowns_df_pure), test_indices) 
  
  train_data <- crowns_df_pure[train_indices, ]
  test_data  <- crowns_df_pure[test_indices, ]
  
  band_cols_cv <- sapply(train_data, is.numeric)
  band_cols_cv <- band_cols_cv & names(train_data) != "poly_id"
  
  train_means <- colMeans(train_data[, band_cols_cv], na.rm = TRUE)
  train_sds  <- apply(train_data[, band_cols_cv], 2, sd, na.rm = TRUE)
  
  train_data[, band_cols_cv] <- scale(train_data[, band_cols_cv], center = train_means, scale = train_sds)
  test_data[, band_cols_cv]  <- scale(test_data[, band_cols_cv], center = train_means, scale = train_sds)
  
  # Removing constant (colinear columns)
  X_train_cv <- train_data[, band_cols_cv]
  constant_cols_cv <- sapply(X_train_cv, function(x) length(unique(x)) == 1)
  X_train_filtered_cv <- X_train_cv[, !constant_cols_cv, drop = FALSE]
  
  if (ncol(X_train_filtered_cv) > 1) {
    cor_matrix_cv <- cor(X_train_filtered_cv, use = "complete.obs")
    high_cor_cols_cv <- names(X_train_filtered_cv)[
      sapply(1:ncol(X_train_filtered_cv), function(j) any(abs(cor_matrix_cv[j, ]) > 0.99 & cor_matrix_cv[j, ] != 1))
    ]
    X_train_filtered_cv <- X_train_filtered_cv[, !names(X_train_filtered_cv) %in% high_cor_cols_cv, drop = FALSE]
  }
  
  train_lda_cv <- data.frame(Liana = train_data$Liana, X_train_filtered_cv)
  n_classes <- length(levels(train_lda_cv$Liana))
  model_cv <- lda(Liana ~ ., data = train_lda_cv, tol = 1e-8, prior = rep(1/n_classes, n_classes))
  
  X_test_cv <- test_data[, band_cols_cv]
  X_test_filtered_cv <- X_test_cv[, !constant_cols_cv, drop = FALSE]
  X_test_filtered_cv <- X_test_filtered_cv[, names(X_train_filtered_cv), drop = FALSE]
  
  test_lda_cv <- data.frame(Liana = test_data$Liana, X_test_filtered_cv)
  predictions_cv <- predict(model_cv, test_lda_cv)
  
  cv_results[[i]] <- list(
    fold = i,
    accuracy = mean(predictions_cv$class == test_data$Liana),
    test_data = test_data,
    predictions = predictions_cv,
    train_indices = train_indices,
    test_indices = test_indices
  )
}

cv_accuracies <- sapply(cv_results, function(x) x$accuracy)
cat("Accuracies per fold:", paste(round(cv_accuracies, 3), collapse = ", "), "\n")
cat("Mean accuracy CV:", round(mean(cv_accuracies), 3), "\n")

all_cv_results <- do.call(rbind, lapply(cv_results, function(x) {
  n_test <- nrow(x$test_data)
  data.frame(
    fold = rep(x$fold, n_test),
    poly_id = x$test_data$poly_id,
    LianaOriginal = x$test_data$LianaOriginal,
    correct = x$predictions$class == x$test_data$Liana,
    predicted = x$predictions$class
  )
}))

all_cv_results_filtered <- subset(all_cv_results, LianaOriginal %in% c("0", "3"))

mean_accuracy_pixel_filtered <- mean(all_cv_results_filtered$correct)
mean_accuracy_pixel_filtered

accuracy_per_crown_cv <- aggregate(correct ~ poly_id + LianaOriginal + fold,
                                   data = all_cv_results,
                                   FUN = mean)
names(accuracy_per_crown_cv)[4] <- "accuracy"

accuracy_points_by_fold <- aggregate(correct ~ poly_id + LianaOriginal + fold,
                                     data = all_cv_results,
                                     FUN = mean)


accuracy_points_by_fold_filtered <- subset(accuracy_points_by_fold,LianaOriginal %in% c("0", "3"))

p3 <- ggplot() +
  geom_boxplot(
    data = accuracy_points_by_fold_filtered,
    aes(x = LianaOriginal, y = correct, color = LianaOriginal),
    fill = NA,       
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_jitter(
    data = accuracy_points_by_fold_filtered,
    aes(x = LianaOriginal, y = correct, shape = factor(fold), fill = factor(fold)),
    width = 0.15,
    alpha = 0.7,
    size = 2.5,
    color = "black",
    stroke = 0.3
  ) +
  scale_color_manual(values = c("0" = "pink1", "3" = "palevioletred4")) +
  scale_fill_manual(values = c("1" = "#B66DFF", "2" = "#FF6DB6", "3" = "#24FF24", "4" = "#6DB6FF", "5"="#DB6D00")) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +
  labs(
    x = "Liana class in drone photos (DP)",
    y = "Mean accuracy per crown",
    title = "Cross validation (StratifiedGroupKFold) on crowns (No filter)",
    caption = paste0("Overall mean accuracy = ", round(mean_accuracy_pixel_filtered * 100, 1), "%"),
    shape = "Fold",
    fill ="Fold"
  ) +
  theme(legend.position = "right")+
  theme_minimal()

p3

all_cv_results$predicted_cv_num <- as.numeric(as.character(all_cv_results$predicted))

predicted_liana_per_crown_cv <- aggregate(
  predicted_cv_num ~ poly_id + LianaOriginal + fold,
  data = all_cv_results,
  FUN = mean
)
names(predicted_liana_per_crown_cv)[4] <- "prediction"

predicted_liana_per_crown_cv$prediction <- predicted_liana_per_crown_cv$prediction * 100


p4 <- ggplot() +
  geom_boxplot(
    data = predicted_liana_per_crown_cv,
    aes(x = LianaOriginal, y = prediction, color = LianaOriginal),
    fill = NA,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_jitter(
    data = predicted_liana_per_crown_cv,
    aes(x = LianaOriginal, y = prediction, shape = factor(fold), fill = factor(fold)),
    width = 0.15,
    alpha = 0.7,
    size = 2.5,
    color = "black",
    stroke = 0.3
  ) +
  scale_color_manual(values = c("0" = "pink1", "1" = "palevioletred1", "2" = "palevioletred3", "3" = "palevioletred4")) +
  scale_fill_manual(values = c("1" = "#B66DFF", "2" = "#FF6DB6", "3" = "#24FF24", "4" = "#6DB6FF", "5"="#DB6D00")) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +
  labs(
    x = "Liana class in drone photos (DP)",
    y = "Liana predicted per crown (%)",
    title = "Mean liana prediction on all fold tested (No filter)",
    caption = paste0("Overall mean accuracy = ", round(mean_accuracy_pixel_filtered * 100, 1), "%"),
    shape = "Fold", 
    fill = "Fold"    
    )+
  theme(legend.position = "right") +
  theme_minimal()

p4

grid.arrange(p3, p4, ncol = 2)

# grid.arrange(p1, p2, p3, p4, ncol = 2)

accuracy_by_fold <- aggregate(
  correct ~ fold,
  data = subset(all_cv_results, LianaOriginal %in% c("0", "3")),
  FUN = mean
)

p4_fold_plots <- lapply(1:5, function(current_fold) {
  fold_data <- subset(predicted_liana_per_crown_cv, fold == current_fold)
  fold_accuracy <- accuracy_by_fold$correct[accuracy_by_fold$fold == current_fold]

  ggplot() +
    geom_boxplot(
      data = fold_data,
      aes(x = LianaOriginal, y = prediction, color = LianaOriginal),
      fill = NA,
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    geom_jitter(
      data = fold_data,
      aes(x = LianaOriginal, y = prediction, fill = factor(fold), shape = factor(fold)),
      width = 0.15,
      alpha = 0.7,
      size = 2.5,
      color = "black", 
      stroke = 0.3
    ) +
    scale_color_manual(values = c("0" = "pink1", "1" = "palevioletred1", "2" = "palevioletred3", "3" = "palevioletred4")) +  
    scale_fill_manual(values = c("1" = "#B66DFF", "2" = "#FF6DB6", "3" = "#24FF24", "4" = "#6DB6FF", "5" = "#DB6D00")) +  
    scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +  
    labs(
      x = "Liana class in drone photos (DP)",
      y = "Liana predicted per crown (%)",
      title = paste("Testing on fold", current_fold),
      caption = paste0("Fold mean accuracy = ", round(fold_accuracy * 100, 1), "%")
    ) +
    theme(legend.position = "none") 
    })

grid.arrange(
  grobs = c(p4_fold_plots, p4),
  ncol = 2,
  top = textGrob("Liana prediction by crown and by fold (NIR filter)",
                 gp = gpar(fontsize = 12, fontface = "bold"))
)

#### For blob ####
# nrow_blobs <- nrow(area_blob_filt)
# ncol_blobs <- ncol(area_blob_filt)
# 
# label_blob_matrix <- as.matrix(label_blob_filt)  # ✅ 2499x1979
# 
# # 2. Créer un raster binaire (1 si blob, 0 sinon)
# #    Si label_blob_filt contient déjà 0/1, on l'utilise directement :
# blob_rast <- raster::raster(label_blob_matrix)
# 
# # 3. Appliquer clump pour identifier les groupes de pixels connectés
# blob_clumped <- raster::clump(blob_rast, directions = 4)  # ✅ directions=4 pour connectivité 4-voisins
# 
# # 4. Extraire les poly_id pour les pixels valides (indices linéaires)
# poly_id_m_blob <- as.matrix(blob_clumped)  # ✅ Matrice 2499x1979 en colonne-major
# poly_id_valid_blob <- poly_id_m_blob[MD_per_pixel$valid_idx]  # ✅ Extraction par indices linéaires
# 
# # 5. Vérifier la cohérence
# if (length(poly_id_valid_blob) != length(label_blob_valid)) {
#   stop("❌ Incohérence entre poly_id_valid_blob et label_blob_valid")
# }
# 
# # 6. Créer blob_df_pure
# blob_df_pure <- as.data.frame(blob_filt_valid)
# blob_df_pure$poly_id <- poly_id_valid_blob
# blob_df_pure$Liana <- as.factor(label_blob_valid)  # ✅ Pas de transformation (0 ou 1)
# blob_df_pure$Type <- factor("blob", levels = c("crowns", "blob"))
# blob_df_pure$LianaOriginal <- NA
# 
# set.seed(123)
# folds_blob <- stratified_group_kfold(blob_df_pure, "poly_id", "Liana", k = 5)
# 
# cv_results_blob <- list()
# for (i in 1:length(folds_blob)) {  
#   test_indices <- folds_blob[[i]]  
#   train_indices <- setdiff(1:nrow(blob_df_pure), test_indices) 
#   
#   train_data <- blob_df_pure[train_indices, ]
#   test_data  <- blob_df_pure[test_indices, ]
#   
#   band_cols_cv <- sapply(train_data, is.numeric)
#   band_cols_cv <- band_cols_cv & names(train_data) != "poly_id"
#   
#   train_means <- colMeans(train_data[, band_cols_cv], na.rm = TRUE)
#   train_sds  <- apply(train_data[, band_cols_cv], 2, sd, na.rm = TRUE)
#   
#   train_data[, band_cols_cv] <- scale(train_data[, band_cols_cv], center = train_means, scale = train_sds)
#   test_data[, band_cols_cv]  <- scale(test_data[, band_cols_cv], center = train_means, scale = train_sds)
#   
#   # Removing constant (colinear columns)
#   X_train_cv <- train_data[, band_cols_cv]
#   constant_cols_cv <- sapply(X_train_cv, function(x) length(unique(x)) == 1)
#   X_train_filtered_cv <- X_train_cv[, !constant_cols_cv, drop = FALSE]
#   
#   if (ncol(X_train_filtered_cv) > 1) {
#     cor_matrix_cv <- cor(X_train_filtered_cv, use = "complete.obs")
#     high_cor_cols_cv <- names(X_train_filtered_cv)[
#       sapply(1:ncol(X_train_filtered_cv), function(j) any(abs(cor_matrix_cv[j, ]) > 0.99 & cor_matrix_cv[j, ] != 1))
#     ]
#     X_train_filtered_cv <- X_train_filtered_cv[, !names(X_train_filtered_cv) %in% high_cor_cols_cv, drop = FALSE]
#   }
#   
#   train_lda_cv <- data.frame(Liana = train_data$Liana, X_train_filtered_cv)
#   n_classes <- length(levels(train_lda_cv$Liana))
#   model_cv <- lda(Liana ~ ., data = train_lda_cv, tol = 1e-8, prior = rep(1/n_classes, n_classes))
#   
#   X_test_cv <- test_data[, band_cols_cv]
#   X_test_filtered_cv <- X_test_cv[, !constant_cols_cv, drop = FALSE]
#   X_test_filtered_cv <- X_test_filtered_cv[, names(X_train_filtered_cv), drop = FALSE]
#   
#   test_lda_cv <- data.frame(Liana = test_data$Liana, X_test_filtered_cv)
#   predictions_cv <- predict(model_cv, test_lda_cv)
#   
#   cv_results[[i]] <- list(
#     fold = i,
#     accuracy = mean(predictions_cv$class == test_data$Liana),
#     test_data = test_data,
#     predictions = predictions_cv,
#     train_indices = train_indices,
#     test_indices = test_indices
#   )
# }

#### Confusion matrix (for StratifiedGroupKFold) ##############################
all_true <- unlist(lapply(cv_results, function(x) x$test_data$Liana))
all_pred <- unlist(lapply(cv_results, function(x) x$predictions$class))

confusion_matrix_cv <- table(True = all_true, Predicted = all_pred)
print(confusion_matrix_cv)

cm_cv <- confusionMatrix(
  factor(all_pred, levels = c("0", "1")),
  factor(all_true, levels = c("0", "1")),
  positive = "1"
)
print(cm_cv) #95% CI [1] = F1-Score

confusion_df_cv <- as.data.frame(confusion_matrix_cv)
colnames(confusion_df_cv) <- c("True", "Predicted", "Freq")

ggplot(confusion_df_cv, aes(x = True, y = Predicted, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color ="white") +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(x = "True", y = "Predicted", title = "Confusion matrix for CV (Illu filter)") +
  theme_minimal()

save.image(file = file.path("Classification_R", "LianaClassification", "myEnvironment.RData"))

