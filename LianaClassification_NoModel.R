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
library(psych)
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

##Functions
source("Classification_R/LianaClassification/Functions.R")

# load('Classification_R/LianaClassification/myEnvironment.RData')


# Uploading data and cropping ##########################################

## HSI mosaic ####

# r <- rast("HSI/mosaic_2016/nou_20160921_l1c_mosaic.bsq") #2016
r <- rast("HSI/mosaic_2023/mosaic/nou_20231025_l1c_mosaic.bsq") #2023
nlyr(r) #number of spectral bands

## Illu ####
illu <- rast("HSI/mosaic_2023/mosaic/Nou2023_Vol2_NS_mosaic_ilu.bsq")
illu <- illu$`Illumination Shadows`
nlyr(illu) # should be 1 (Illumination shadows)
# hist(illu, breaks = 100)
# plot(illu)

## Manually delineated (MD) ####

blob <- vect("Lianes/LianaOrNotLiana_55_45.shp")
blob <- project(blob, crs(r))
area_blob <- crop(r, blob, mask=TRUE)
# area_blob <- mask(r, blob)
# plot (area_blob) # /!\ creates one plot per band
# plot(area_blob[[2]]) # the number is the rank of the band not it's value
label_blob <- rasterize(blob, area_blob, field = "Type") # Taking into account the liana label

## Crowns from drone photos (DP) ####
crowns <- vect("Couronnes/Crowns_Liana_2025_GP_TRUE.shp")
crowns <- project(crowns, crs(r))
area_crowns <- crop(r, crowns, mask=TRUE)
label_crowns <- rasterize(crowns, area_crowns, field = "LianaFullC")

## Proportion of lianas before filtering #####################################
table(values(crowns)[["LianaFullC"]])
prop.table(table(values(crowns)[["LianaFullC"]]))

table(values(blob)[["Type"]])
prop.table(table(values(blob)[["Type"]]))
# Filtering pixels ##########################################################

## Cleaning extreme values / NA ##
# 
# # min(area_crowns)
# # max(area_crowns)
# # sum(is.na(as.matrix(area_crowns)))
# area_crowns <- clean_extreme_values_toNA(area_crowns)
# # min(area_crowns)
# # max(area_crowns)
# # sum(is.na(as.matrix(area_crowns)))
# 
# # min(area_blob)
# # max(area_blob)
# # sum(is.na(as.matrix(area_blob)))
# area_blob <- clean_extreme_values_toNA(area_blob)
# # min(area_blob)
# # max(area_blob)
# # sum(is.na(as.matrix(area_blob)))
# 
# ## Removing range & bands with more than 95% na
# 
# area_crowns <- filter_good_bands(area_crowns, max_na_pct_band = 0.99)
# area_blob <- filter_good_bands(area_blob, max_na_pct_band = 0.99)
# 
# wl_crowns <- round(extract_wavelength(names(area_crowns)), 2)
# wl_blob   <- round(extract_wavelength(names(area_blob)), 2)
# common_wl <- intersect(wl_crowns, wl_blob)
# common_wl <- sort(common_wl)
# 
# area_crowns <- area_crowns[[which(wl_crowns %in% common_wl)]]
# area_blob   <- area_blob[[which(wl_blob %in% common_wl)]]
# 
# writeRaster(area_crowns, "area_crowns_clean.tif", overwrite = T)
# writeRaster(area_blob, "area_blob_clean.tif", overwrite = T)

#Importing cleaned rasters ####
area_crowns <- rast ("area_crowns_clean.tif")
area_blob <- rast ("area_blob_clean.tif")

## Removing pixel based on illu ####
# Only used with the 2016 data (illu file missing fo 2023)

med_illu <- global(illu[[1]], fun = median, na.rm = TRUE)$global
perc90_illu <- global(illu[[1]], fun = quantile, probs = 0.9, na.rm = TRUE)$X90.
perc95_illu <- global(illu[[1]], fun = quantile, probs = 0.95, na.rm = TRUE)$X90.

## Removing pixel based on VNIR reflectance median ####

## Removing pixel based on NIR reflectance####

# names(area_crowns)
# crowns_nir <- area_crowns[[118:157]]
# mean_crowns <- mean(crowns_nir, na.rm =TRUE)
# hist(mean_crowns)


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
m_crowns <- mean(area_crowns_filt, na.rm = TRUE)
# area_crowns_filt <- area_crowns_filt - m_crowns 
area_crowns_filt <- area_crowns_filt / m_crowns

m_blob <- mean(area_blob_filt, na.rm = TRUE)
# area_blob_filt   <- area_blob_filt - m_blob
area_blob_filt <- area_blob_filt / m_blob

# save(area_crowns_filt, area_blob_filt, file = "Classification_R/LianaClassification/areas_NIR_norm.RData")
# load("Classification_R/LianaClassification/areas_NIR_norm.RData")

## Correlation  between filters ##########################################

### Blob data prep
reference_blob <- area_blob_filt_norm
mask_valid_blob <- !is.na(reference_blob[[1]])  

pixels_illu_blob <- as.matrix(illu_blob_out$label_filt)[mask_valid_blob]
pixels_nir_blob  <- as.matrix(NIR_blob_out$label_filt)[mask_valid_blob]
pixels_vnir_blob <- as.matrix(VNIR_blob_out$label_filt)[mask_valid_blob]

stopifnot(length(pixels_illu_blob) == length(pixels_nir_blob), length(pixels_illu_blob) == length(pixels_vnir_blob))


### Crowns data prep
reference_crowns <- area_crowns_filt
mask_valid_crowns <- !is.na(as.matrix(reference_crowns[[1]]))  

pixels_illu_crowns <- as.matrix(illu_crowns_out$label_filt)[mask_valid_crowns]
pixels_nir_crowns  <- as.matrix(NIR_crowns_out$label_filt)[mask_valid_crowns]
pixels_vnir_crowns <- as.matrix(VNIR_crowns_out$label_filt)[mask_valid_crowns]

stopifnot(length(pixels_illu_crowns) == length(pixels_nir_crowns), length(pixels_illu_crowns) == length(pixels_vnir_crowns))


### Chi² ####

### Blob
run_chi2(pixels_illu_blob, pixels_nir_blob, "illu", "NIR")
run_chi2(pixels_nir_blob, pixels_vnir_blob, "NIR", "VNIR")
run_chi2(pixels_vnir_blob, pixels_illu_blob, "VNIR", "illu")

### Crowns
run_chi2(pixels_illu_crowns, pixels_nir_crowns, "illu", "NIR")
run_chi2(pixels_nir_crowns, pixels_vnir_crowns, "NIR", "VNIR")
run_chi2(pixels_vnir_crowns, pixels_illu_crowns, "VNIR", "illu")

### Cohen's Kappa ####

# BLOB 
mask_valid_blob <- !is.na(as.matrix(reference_blob[[1]])) &
  !is.na(as.matrix(illu_blob_out$label_filt)) &
  !is.na(as.matrix(NIR_blob_out$label_filt)) &
  !is.na(as.matrix(VNIR_blob_out$label_filt))

pixels_ref_blob    <- as.vector(as.matrix(reference_blob[[1]]))[mask_valid_blob]
pixels_illu_blob_med   <- as.vector(as.matrix(illu_blob_out$label_filt))[mask_valid_blob]
pixels_nir_blob_20    <- as.vector(as.matrix(NIR_blob_out$label_filt))[mask_valid_blob]
pixels_vnir_blob   <- as.vector(as.matrix(VNIR_blob_out$label_filt))[mask_valid_blob]

# CROWNS
mask_valid_crowns <- !is.na(as.matrix(reference_crowns[[1]])) &
  !is.na(as.matrix(illu_crowns_out$label_filt)) &
  !is.na(as.matrix(NIR_crowns_out$label_filt)) &
  !is.na(as.matrix(VNIR_crowns_out$label_filt))

pixels_ref_crowns  <- as.vector(as.matrix(reference_crowns[[1]]))[mask_valid_crowns]
pixels_illu_crowns_med <- as.vector(as.matrix(illu_crowns_out$label_filt))[mask_valid_crowns]
pixels_nir_crowns_20  <- as.vector(as.matrix(NIR_crowns_out$label_filt))[mask_valid_crowns]
pixels_vnir_crowns <- as.vector(as.matrix(VNIR_crowns_out$label_filt))[mask_valid_crowns]

run_kappa(pixels_illu_blob, pixels_nir_blob, "illu", "NIR")
run_kappa(pixels_nir_blob, pixels_vnir_blob, "NIR", "VNIR")
run_kappa(pixels_vnir_blob, pixels_illu_blob, "VNIR", "illu")

# Dimensionality reduction and data visualization ###############################

## Data prep ##############################################################

# load("Classification_R/LianaClassification/valid_NIR.RData")

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

# save(DP_per_pixel,DP_per_poly, MD_per_pixel, MD_per_poly, file = "Classification_R/LianaClassification/valid_NIR.RData"  )


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

## For DP per polygon 
pca_out_crowns_poly <- PCA(crowns_poly_valid, scale.unit = TRUE, graph = FALSE)
var_exp_crowns_poly <- pca_out_crowns_poly$eig[1:2, "percentage of variance"]

pca_plot_crowns_poly <- data.frame(
  x = pca_out_crowns_poly$ind$coord[,1],
  y = pca_out_crowns_poly$ind$coord[,2],
  LianaFullCrown = as.factor(label_crowns_poly_valid)
)

## For MD per pixel
pca_out_blob <- PCA(blob_filt_valid, scale.unit = TRUE, graph = FALSE)
var_exp_blob <- pca_out_blob$eig[1:2, "percentage of variance"]

pca_plot_blob <- data.frame(
  x = pca_out_blob$ind$coord[,1],
  y = pca_out_blob$ind$coord[,2],
  Type = as.factor(label_blob_valid)
)

## For MD per polygon 
pca_out_blob_poly <- PCA(blob_poly_valid, scale.unit = TRUE, graph = FALSE)
var_exp_blob_poly <- pca_out_blob_poly$eig[1:2, "percentage of variance"]

pca_plot_blob_poly <- data.frame(
  x = pca_out_blob_poly$ind$coord[,1],
  y = pca_out_blob_poly$ind$coord[,2],
  Type = as.factor(label_blob_poly_valid)
)

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

## UMAP #########################################################################

#For DP per pixel

umap_out_crowns <- umap2( seed = 2, crowns_filt_valid,  n_neighbors = min(15, nrow(crowns_filt_valid) - 1))

umap_plot_crowns <- data.frame(
  x = umap_out_crowns[,1],
  y = umap_out_crowns[,2],
  LianaFullCrown = as.factor(label_crowns_valid)
)

#For DP per polygon
umap_out_crowns_poly <- umap2( seed = 2, crowns_poly_valid)

umap_plot_crowns_poly <- data.frame(
  x = umap_out_crowns_poly[,1],
  y = umap_out_crowns_poly[,2],
  LianaFullCrown = as.factor(label_crowns_poly_valid))

#For MD per pixel
umap_out_blob <- umap2( seed = 2, blob_filt_valid)

umap_plot_blob <- data.frame(
  x = umap_out_blob[,1],
  y = umap_out_blob[,2],
  Type = as.factor(label_blob_valid)
)

#For MD per polygon
umap_out_poly <- umap2( seed = 2, blob_poly_valid)

umap_plot_poly <- data.frame(
  x = umap_out_poly[,1],
  y = umap_out_poly[,2],
  Type = as.factor(label_blob_poly_valid)
)



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
crowns_poly_valid <- sensor_to_keep(crowns_poly_valid, "swir")
blob_poly_valid <- sensor_to_keep(blob_poly_valid, "swir")


spec_crowns <- mean_spectrum_by_group(crowns_poly_valid, label_crowns_poly_valid, "DP")
spec_blob   <- mean_spectrum_by_group(blob_poly_valid, label_blob_poly_valid, "MD")

all_spectra <- rbind(spec_crowns, spec_blob)

#selected <- subset(all_spectra, group %in% c("DP_CR0", "DP_CR3", "MD_0", "MD_1") & all_spectra$wavelength<900) #VNIR only
selected <- subset(all_spectra, group %in% c("DP_CR0", "DP_CR3", "MD_0", "MD_1"))
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


# Median spectrum comparison ###################################################

spec_crowns <- med_spectrum_by_group(crowns_poly_valid, label_crowns_poly_valid, "DP")
spec_blob   <- med_spectrum_by_group(blob_poly_valid, label_blob_poly_valid, "MD")

all_spectra <- rbind(spec_crowns, spec_blob)

#med_selected <- subset(all_spectra, group %in% c("DP_CR0", "DP_CR3", "MD_0", "MD_1") & all_spectra$wavelength<900) #VNIR only
med_selected <- subset(all_spectra, group %in% c("DP_CR0", "DP_CR3", "MD_0", "MD_1"))
med_selected <- add_segments(med_selected)

med_selected_smooth <- as.data.frame(med_selected)
med_selected_smooth$med_smooth <- med_selected_smooth$med
med_selected_smooth$sd_smooth   <- med_selected_smooth$sd

for (g in unique(med_selected_smooth$group)) {
  for (s in unique(med_selected_smooth$segment[med_selected_smooth$group == g])) {
    
    rows <- which(med_selected_smooth$group == g & med_selected_smooth$segment == s)
    
    if (length(rows) >= 11) {
      med_selected_smooth$med_smooth[rows] <- sgolayfilt(med_selected_smooth$med[rows], p = 3, n = 11)
      med_selected_smooth$sd_smooth[rows]   <- sgolayfilt(med_selected_smooth$sd[rows], p = 3, n = 11)
    }
  }
}


# Linear Discriminant Analysis (LDA) #########################################

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
lda_crowns_accuracy #0.6248766 (NIR)


### Stratified folded cross-validation
predictors <- crowns_df[ , -which(names(crowns_df) %in% c("Liana"))]

ctrl_crowns <- trainControl(
   method = "cv",
   number = 5,
   savePredictions = "final" )

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

### Within MD dataset
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
lda_blob_accuracy #0.870233 (NIR)

# ### Stratified folded cross-validation
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

# ## Between MD and DP datasets
# crowns_df <- as.data.frame(crowns_filt_valid)
# crowns_df$Liana <- as.character(label_crowns_valid)
# crowns_df <- crowns_df[crowns_df$Liana %in% c('0','3'), ]
# crowns_df$Liana[crowns_df$Liana == '3'] <- '1'
# crowns_df$Liana <- as.factor(crowns_df$Liana)
# type_crowns <- rep(c('crowns'), nrow(crowns_df))
# crowns_df$Type <- as.factor(type_crowns)
# 
# blob_df <- as.data.frame(blob_filt_valid)
# blob_df$ Liana <- as.factor(label_blob_valid)
# type_blob <- rep(c('blob'), nrow(blob_df))
# blob_df$Type <- as.factor(type_blob)
# 
# all_df <- rbind(crowns_df, blob_df)
# 
# set.seed(123)
# all_df[,1:ncol(crowns_filt_valid)]<- scale(all_df[,1:ncol(crowns_filt_valid)])
# train_all <- all_df[all_df$Type == 'blob', ]
# test_all <- all_df[all_df$Type== 'crowns', ]
# # test_all <- all_df[all_df$Type == 'blob', ]
# # train_all <- all_df[all_df$Type== 'crowns', ]
# 
# n_classes <- length(levels(train_all$Liana))
# model_all <- lda(Liana ~ . - Type, data = train_all, prior = rep(1/n_classes, n_classes)) # remove Type from predictive values
# predictions_all <- predict(model_all, test_all)
# lda_all_accuracy <- mean(predictions_all$class==test_all$Liana)
# lda_all_accuracy 
# 
# # Predicted liana proportion
# liana_prop_pred <- sum(predictions_all$class== 1) / nrow(test_all)
# noliana_prop_pred <- sum(predictions_all$class== 0) / nrow(test_all)
# liana_prop_pred
# noliana_prop_pred

## LDA trained on pure labels to check % of liana predicted ####
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

### Trained on blob / Tested on crowns 
# set.seed(123)
# band_cols <- sapply(all_df, is.numeric)
# band_cols <- names(all_df)[band_cols]
# 
# all_df[, band_cols] <- scale(all_df[, band_cols])
# 
# train_all_pure <- all_df[all_df$Type == "blob", ]
# test_all_pure <- all_df[all_df$Type == "crowns", ]
# 
# X_train <- train_all_pure[, band_cols]
# constant_cols <- sapply(X_train, function(x) length(unique(x)) == 1)
# X_train_filtered <- X_train[, !constant_cols, drop = FALSE]
# band_cols_filtered <- names(X_train_filtered)
# 
# if (ncol(X_train_filtered) > 1) {
#   cor_matrix <- cor(X_train_filtered, use = "complete.obs")
#   high_cor_cols <- names(X_train_filtered)[
#     sapply(1:ncol(X_train_filtered), function(i) any(abs(cor_matrix[i, ]) > 0.99 & cor_matrix[i, ] != 1))
#   ]
#   band_cols_filtered <- setdiff(band_cols_filtered, high_cor_cols)
#   X_train_filtered <- X_train_filtered[, band_cols_filtered, drop = FALSE]
# }
# 
# train_data_blob <- data.frame(Liana = train_all_pure$Liana, X_train_filtered)
# n_classes <- length(levels(train_data_blob$Liana))
# model_blob_to_crown <- lda(Liana ~ ., data = train_data_blob, tol = 1e-8)
# 
# X_test <- test_all_pure[, band_cols]
# X_test_filtered <- X_test[, !constant_cols, drop = FALSE]
# X_test_final <- X_test_filtered[, band_cols_filtered, drop = FALSE]
# test_data_crown <- data.frame(Liana = test_all_pure$Liana, X_test_final)
# 
# predictions_blob_to_crown <- predict(model_blob_to_crown, test_data_crown)
# test_all_pure$predicted_blob <- predictions_blob_to_crown$class
# test_all_pure$correct_blob <- test_all_pure$predicted_blob == test_all_pure$Liana
# 
# accuracy_blob_to_crown <- mean(test_all_pure$correct_blob)
# cat("Mean accuracy per crown:", round(accuracy_blob_to_crown, 3), "\n")
# 
# accuracy_per_crown_blob <- aggregate(correct_blob ~ poly_id + LianaOriginal,
#                                      data = test_all_pure,
#                                      FUN = mean)
# 
# names(accuracy_per_crown_blob)[3] <- "accuracy"
# 
# p1 <- ggplot(accuracy_per_crown_blob, aes(x = LianaOriginal, y = accuracy, fill = LianaOriginal)) +
#   geom_boxplot() +
#   geom_jitter(width = 0.15, alpha = 0.5) +
#   scale_fill_manual(values = c("0" = "pink1", "1" = "palevioletred1",
#                                "2" = "palevioletred3", "3" = "palevioletred4")) +
#   labs(x = "Liana class in drone photos (DP)", y = "Mean accuracy per crown",
#        title = "LDA trained on MD, tested on DP (Illu filter)") +
#   theme(legend.position = "none")
# p1
# 
# test_all_pure$predicted_blob_num <- as.numeric(as.character(test_all_pure$predicted_blob))
# 
# predicted_liana_per_crown_blob <- aggregate(
#   predicted_blob_num ~ poly_id + LianaOriginal,
#   data = test_all_pure,
#   FUN = mean
# )
# names(predicted_liana_per_crown_blob)[3] <- "prediction"
# 
# 
# predicted_liana_per_crown_blob$prediction <- predicted_liana_per_crown_blob$prediction * 100
# 
# p2 <- ggplot(predicted_liana_per_crown_blob, aes(x = LianaOriginal, y = prediction, fill = LianaOriginal)) +
#   geom_boxplot() +
#   geom_jitter(width = 0.15, alpha = 0.5) +
#   scale_fill_manual(values = c("0" = "pink1", "1" = "palevioletred1",
#                                "2" = "palevioletred3", "3" = "palevioletred4")) +
#   labs(x = "Liana class in drone photos (DP)", y = "Liana predicted per crown (%)",
#        title = "LDA trained on MD, tested on DP (Illu filter)") +
#   theme(legend.position = "none")
# 
# p2
# 
# grid.arrange(p1, p2, ncol = 2)

### StratifiedGroupKFold trained on pure labels ####

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
  
  test_data_filtered <- test_data[test_data$Liana %in% c("0", "3"), ]
  pred_filtered <- predictions_cv$class[test_data$Liana %in% c("0", "3")]
  
  cm_filtered <- confusionMatrix(
    factor(pred_filtered, levels = c("0", "3")),
    factor(test_data_filtered$Liana, levels = c("0", "3")),
    positive = NULL)
  
  f1_filtered <- calculate_f1_manual(test_data_filtered$Liana, pred_filtered)
  
  cv_results[[i]] <- list(
    fold = i,
    accuracy = mean(predictions_cv$class == test_data$Liana),
    test_data = test_data,
    f1_filtered = f1_filtered,
    predictions = predictions_cv,
    train_indices = train_indices,
    test_indices = test_indices
  )
}

cv_accuracies <- sapply(cv_results, function(x) x$accuracy)
cv_f1 <- sapply(cv_results, function(x) x$f1_filtered)

cat("Accuracies per fold:", paste(round(cv_accuracies, 3), collapse = ", "), "\n")
cat("Mean accuracy CV:", round(mean(cv_accuracies), 3), "\n")
cat("F1-scores per fold:", paste(round(cv_f1, 3), collapse = ", "), "\n")
cat("Mean F1-score:", round(mean(cv_f1, na.rm = TRUE), 3), "\n")

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
    title = "Cross validation (StratifiedGroupKFold) on crowns (NIR filter)",
    caption = paste0("Overall mean accuracy (0&3) = ", round(mean_accuracy_pixel_filtered * 100, 1), "%"),
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
  geom_jitter(
    data = predicted_liana_per_crown_cv,
    aes(x = LianaOriginal, y = prediction, shape = factor(fold), fill = factor(fold)),
    width = 0.15,
    alpha = 0.5,
    size = 2.5,
    color = "black",
    stroke = 0.3
  ) +
  geom_boxplot(
    data = predicted_liana_per_crown_cv,
    aes(x = LianaOriginal, y = prediction, color = LianaOriginal),
    fill = NA,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("0" = "pink1", "1" = "palevioletred1", "2" = "palevioletred3", "3" = "palevioletred4")) +
  scale_fill_manual(values = c("1" = "grey", "2" = "grey", "3" = "grey", "4" = "grey", "5"="grey")) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +
  labs(
    x = "Liana class in drone photos (DP)",
    y = "Liana predicted per crown (%)",
    title = "Mean liana prediction on all fold tested (NIR filter)",
    caption = paste0("Overall mean accuracy (0&3)= ", round(mean_accuracy_pixel_filtered * 100, 1),  "% | ",
                     "Mean F1-score (0&3) = ", round(mean(cv_f1, na.rm = TRUE), 3)),
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
  fold_f1 <- cv_f1[current_fold] 


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
      caption = paste0("Fold mean accuracy = ", round(fold_accuracy * 100, 1), "% | ",
                       "Fold F1-score (0&3) = ", round(cv_f1[current_fold], 3))
    ) +
    theme(legend.position = "none") +
    theme_minimal()
    })

grid.arrange(
  grobs = c(p4_fold_plots, p4),
  ncol = 2,
  top = textGrob("Liana prediction by crown and by fold (NIR filter)",
                 gp = gpar(fontsize = 12, fontface = "bold"))
)

#### For blob ####
blob_df_pure$LianaOriginal <- blob_df_pure$Liana
blob$poly_id <- seq_len(nrow(blob))
poly_id_r <- rasterize(blob, area_blob_filt, field = "poly_id")
poly_id_m <- as.matrix(poly_id_r)[, 1]
poly_id_valid <- poly_id_m[MD_per_pixel$valid_idx]
blob_df_pure$poly_id <- ifelse(is.na(poly_id_valid), 0, poly_id_valid)
set.seed(123)
folds_blob <- stratified_group_kfold(blob_df_pure, "poly_id", "Liana", k = 5)

cv_results <- list()
for (i in 1:length(folds_blob)) {  
  test_indices <- folds_blob[[i]]  
  train_indices <- setdiff(1:nrow(blob_df_pure), test_indices) 
  
  train_data <- blob_df_pure[train_indices, ]
  test_data  <- blob_df_pure[test_indices, ]
  
  band_cols_cv <- sapply(train_data, is.numeric)
  band_cols_cv <- band_cols_cv & names(train_data) != "poly_id"
  
  train_means <- colMeans(train_data[, band_cols_cv], na.rm = TRUE)
  train_sds  <- apply(train_data[, band_cols_cv], 2, sd, na.rm = TRUE)
  
  train_data[, band_cols_cv] <- scale(train_data[, band_cols_cv], center = train_means, scale = train_sds)
  test_data[, band_cols_cv]  <- scale(test_data[, band_cols_cv], center = train_means, scale = train_sds)
  
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
  
  f1_macro <- calculate_f1_manual(test_data$Liana, predictions_cv$class)
  
  cv_results[[i]] <- list(
    fold = i,
    accuracy = mean(predictions_cv$class == test_data$Liana),
    f1_macro = f1_macro,
    test_data = test_data,
    predictions = predictions_cv,
    train_indices = train_indices,
    test_indices = test_indices
  )
}

cv_accuracies <- sapply(cv_results, function(x) x$accuracy)
cv_f1_macro <- sapply(cv_results, function(x) x$f1_macro)

cat("Accuracies per fold:", paste(round(cv_accuracies, 3), collapse = ", "), "\n")
cat("Mean accuracy CV:", round(mean(cv_accuracies), 3), "\n\n")
cat("F1-scores (macro) per fold:", paste(round(cv_f1_macro, 3), collapse = ", "), "\n")
cat("Mean F1-score (macro) CV:", round(mean(cv_f1_macro, na.rm = TRUE), 3), "\n")

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

accuracy_per_crown_cv <- aggregate(correct ~ poly_id + LianaOriginal + fold,
                                   data = all_cv_results,
                                   FUN = mean)
names(accuracy_per_crown_cv)[4] <- "accuracy"

accuracy_points_by_fold <- aggregate(correct ~ poly_id + LianaOriginal + fold,
                                     data = all_cv_results,
                                     FUN = mean)

all_cv_results$predicted_cv_num <- as.numeric(as.character(all_cv_results$predicted))

mean_accuracy_pixel <- mean(all_cv_results$correct)

predicted_liana_per_crown_cv <- aggregate(
  predicted_cv_num ~ poly_id + LianaOriginal + fold,
  data = all_cv_results,
  FUN = mean
)
names(predicted_liana_per_crown_cv)[4] <- "prediction"

predicted_liana_per_crown_cv$prediction <- predicted_liana_per_crown_cv$prediction * 100


p4_blob <- ggplot() +
  geom_jitter(
    data = predicted_liana_per_crown_cv,
    aes(x = LianaOriginal, y = prediction, shape = factor(fold), fill = factor(fold)),
    width = 0.15,
    alpha = 0.5,
    size = 2.5,
    color = "black",
    stroke = 0.3
  ) +
  geom_boxplot(
    data = predicted_liana_per_crown_cv,
    aes(x = LianaOriginal, y = prediction, color = LianaOriginal),
    fill = NA,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("0" = "#ABABAB", "1" = "#001312")) +
  scale_fill_manual(values = c("1" = "grey", "2" = "grey", "3" = "grey", "4" = "grey", "5"="grey")) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +
  labs(
    x = "Liana class in manual delineation",
    y = "Liana predicted per crown (%)",
    title = "Mean liana prediction on all fold tested (NIR filter)",
    caption = paste0("Overall mean accuracy = ", round(mean_accuracy_pixel_filtered * 100, 1), "% | ",
                     "Mean F1-score = ", round(mean(cv_f1_macro, na.rm = TRUE), 3)),
    shape = "Fold", 
    fill = "Fold"    
  )+
  theme(legend.position = "right") +
  theme_minimal()

p4_blob

accuracy_by_fold <- aggregate(
  correct ~ fold,
  data = subset(all_cv_results, LianaOriginal %in% c("0", "1")),
  FUN = mean
)

p4_fold_plots_blob <- lapply(1:5, function(current_fold) {
  fold_data <- subset(predicted_liana_per_crown_cv, fold == current_fold)
  fold_accuracy <- accuracy_by_fold$correct[accuracy_by_fold$fold == current_fold]
  fold_f1 <- cv_f1_macro[current_fold] 
  
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
    scale_color_manual(values = c("0" = "#ABABAB", "1" = "#001312")) +  
    scale_fill_manual(values = c("1" = "#B66DFF", "2" = "#FF6DB6", "3" = "#24FF24", "4" = "#6DB6FF", "5" = "#DB6D00")) +  
    scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +  
    labs(
      x = "Liana class in manual delineation",
      y = "Liana predicted per patch (%)",
      title = paste("Testing on fold", current_fold),
      caption = paste0("Fold mean accuracy = ", round(fold_accuracy * 100, 1), "% | ",
                       "Fold F1-score = ", round(fold_f1, 3))
    ) +
    theme(legend.position = "none")+
    theme_minimal()
})

grid.arrange(
  grobs = c(p4_fold_plots_blob, p4_blob),
  ncol = 2,
  top = textGrob("Liana prediction by patch and by fold (NIR filter)",
                 gp = gpar(fontsize = 12, fontface = "bold"))
)

### F1 score graphs
all_predictions <- do.call(rbind, lapply(cv_results, function(x) {
  data.frame(
    fold = rep(x$fold, nrow(x$test_data)),
    poly_id = x$test_data$poly_id,
    true = x$test_data$Liana,
    pred = x$predictions$class
  )
}))

f1_by_poly <- do.call(rbind, lapply(cv_results, function(x) {
  unique_polys <- unique(x$test_data$poly_id)  
  results <- lapply(unique_polys, function(poly) {
    poly_pixels <- x$test_data[x$test_data$poly_id == poly, ]
    poly_preds <- x$predictions$class[x$test_data$poly_id == poly]
    
    if (nrow(poly_pixels) == 0) return(NULL)  
  
    cm <- table(poly_pixels$Liana, poly_preds)
    
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
    
    data.frame(
      fold = x$fold,
      poly_id = poly,
      LianaOriginal = unique(poly_pixels$Liana), 
      f1_0 = ifelse("0" %in% names(f1_scores), f1_scores["0"], NA),
      f1_1 = ifelse("1" %in% names(f1_scores), f1_scores["1"], NA),
      f1_macro = mean(f1_scores, na.rm = TRUE) 
    )
  })
  do.call(rbind, results)  
}))

f1_by_poly$f1_mean <- f1_by_poly$f1_macro * 100 

p4_f1_blob <- ggplot() +
  geom_jitter(
    data = f1_by_poly,
    aes(x = LianaOriginal, y = f1_mean, shape = factor(fold), fill = factor(fold)),
    width = 0.15, alpha = 0.5, size = 2.5, color = "black", stroke = 0.3
  ) +
  geom_boxplot(
    data = f1_by_poly,
    aes(x = LianaOriginal, y = f1_mean, color = LianaOriginal),
    fill = NA, linewidth = 0.8, show.legend = FALSE
  ) +
  scale_color_manual(values = c("0" = "#ABABAB", "1" = "#001312")) +
  scale_fill_manual(values = c("1" = "grey", "2" = "grey", "3" = "grey", "4" = "grey", "5" = "grey")) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +
  labs(
    x = "Liana class in manual delineation",
    y = "F1-score per patch (%)",
    title = "F1-score by patch and by fold (NIR filter)",
    caption = paste0("Overall mean F1-score = ", round(mean(cv_f1_macro, na.rm = TRUE), 3)),
    shape = "Fold", fill = "Fold"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

p4_f1_fold_plots_blob <- lapply(1:5, function(current_fold) {
  fold_data <- subset(f1_by_poly, fold == current_fold)
  fold_f1 <- cv_f1_macro[current_fold]  # F1-score moyen du fold
  
  ggplot() +
    geom_boxplot(
      data = fold_data,
      aes(x = LianaOriginal, y = f1_mean, color = LianaOriginal),
      fill = NA, linewidth = 0.8, show.legend = FALSE
    ) +
    geom_jitter(
      data = fold_data,
      aes(x = LianaOriginal, y = f1_mean, fill = factor(fold), shape = factor(fold)),
      width = 0.15, alpha = 0.7, size = 2.5, color = "black", stroke = 0.3
    ) +
    scale_color_manual(values = c("0" = "#ABABAB", "1" = "#001312")) +
    scale_fill_manual(values = c("1" = "#B66DFF", "2" = "#FF6DB6", "3" = "#24FF24", "4" = "#6DB6FF", "5" = "#DB6D00")) +
    scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +
    labs(
      x = "Liana class in manual delineation",
      y = "F1-score per patch (%)",
      title = paste("F1-score on fold", current_fold),
      caption = paste0("Fold mean F1-score = ", round(fold_f1, 3))
    ) +
    theme_minimal() +
    theme(legend.position = "none")
})

grid.arrange(
  grobs = c(p4_f1_fold_plots_blob, p4_f1_blob),
  ncol = 2,
  top = textGrob("F1-score by patch and by fold (NIR filter)",
                 gp = gpar(fontsize = 12, fontface = "bold"))
)

### StratifiedGroupKFold trained on all labels ####


# Boxplot illu ####
label_blob_illu <- rasterize(blob, illu, field = "Type")       # Labels MD
combined_blob <- c(illu, label_blob_illu)
names(combined_blob) <- c("illumination", "label")
df_blob <- as.data.frame(combined_blob, na.rm = TRUE)
df_blob$label <- factor(df_blob$label,levels = c(0, 1),labels = c("MD0", "MD1"))
df_blob$type <- "MD"

label_crown_illu <- rasterize(crowns, illu, field = "LianaFullC")  # Labels DP
combined_crowns <- c(illu, label_crowns_illu)
names(combined_crowns) <- c("illumination", "label")
df_crowns <- as.data.frame(combined_crowns, na.rm = TRUE)
df_crowns$type <- "CR"

df_all <- bind_rows(df_blob, df_crowns)