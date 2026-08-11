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
library(stats)
library(Rtsne)
library(ggplot2)
library(paletteer)
library(uwot)
library(FactoMineR)
library(fpc)

## Set working Directory ####
setwd ("D:/DataCamille/")


# Uploading data and cropping ##########################################

## HSI mosaic ####
r <- rast("HSI/mosaic_2016/nou_20160921_l1c_mosaic.bsq")
# r <- rast("HSI/mosaic_2023/mosaic/nou_20231025_l1c_mosaic.bsq")
nlyr(r) #number of spectral bands

## Illu ####
illu <- rast("HSI/mosaic_2016/nou_20160921_ilu_mosaic.bsq")
#illu <- rast("HSI/mosaic_2023/mosaic/nou_20231025_l1c_mosaic.bsq")
nlyr(illu)
#hist(illu, breaks = 100)
#plot(illu)

## Manually delineated (MD) ####
blob <- vect("Lianes/LianaOrNotLiana_55_45.shp")
blob <- project(blob, crs(r))

area_blob <- crop(r, blob, mask=TRUE)
#plot (area_blob) # /!\ creates one plot per band
plot(area_blob[[3]]) # the number is the rank of the band not it's value

# Taking into account the liana label
label_blob <- rasterize(blob, area_blob, field = "Type")

# Normalized data
normalized_blob <- area_blob / mean(area_blob, na.rm = TRUE)

## Crowns from drone photos (DP) ####
crowns <- vect("Couronnes/Crowns_Liana_2025_GP.shp")
crowns <- project(crowns, crs(r))

area_crowns <- crop(r, crowns, mask=TRUE)

label_crowns <- rasterize(crowns, area_crowns, field = "LianaFullC")

normalized_crowns <- area_crowns / mean(area_crowns, na.rm =TRUE)

# Removing pixel based on illu ###########################################

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

##For DP (crowns) ####

illu_crowns_out <- filter_by_illu(illu, normalized_crowns, label_crowns, med_illu)

##For MD (blob) ####

illu_blob_out <- filter_by_illu(illu, normalized_blob, label_blob, med_illu)

# tSNE ####################################################################################

tsne_per_pixel <- function(area_filtered, label_filtered, illu_crop_area){
  area_matrix <- as.matrix(area_filtered)
  label_matrix <- as.matrix (label_filtered)[,1]
  
  na_per_band <- colSums(is.na(area_matrix))
  good_bands <- na_per_band < nrow(area_matrix)
  # sum(good_bands)
  
  valid_idx <- which(!is.na(label_matrix) & complete.cases(area_matrix[, good_bands]))
  # length(valid_idx)
  
  area_valid <- area_matrix[valid_idx, good_bands]
  labels_valid <- label_matrix[valid_idx]
  # dim(crowns_filt_valid)
  # length(labels_valid)
  
  illu_matrix <- as.matrix(illu_crop_area)[, 1]
  illu_valid <- illu_matrix[valid_idx]
  
  output_list <- list(
    matrix = area_valid,
    labels = labels_valid,
    illu = illu_valid,
    goodbands = good_bands)
  
  return(output_list)
} 

tsne_per_polygon <- function(area_filtered, area_vector, label_name, illu_raster){
  
  # Mean spectrum per polygon
  poly_spectra <- extract(area_filtered, area_vector, fun = mean, na.rm = TRUE, ID = TRUE)
  # poly_labels <- area_vector$label_name[poly_spectra$ID]
  poly_labels <- values(area_vector)[[label_name]][poly_spectra$ID]
  
  # Create the matrix
  poly_matrix <- as.matrix(poly_spectra[, -1])
  
  na_per_band_poly <- colSums(is.na(poly_matrix))
  good_bands_poly <- na_per_band_poly < nrow(poly_matrix)
  
  poly_valid <- poly_matrix[, good_bands_poly]
  
  # Remove polygons full of NA
  valid_poly_idx <- complete.cases(poly_valid)
  
  poly_valid <- poly_valid[valid_poly_idx, ]
  label_poly_valid <- poly_labels[valid_poly_idx]
  dim(poly_valid)
  length(label_poly_valid)
  
  # poly_illu <- extract(illu_raster, area_vector, fun = mean, na.rm = TRUE, ID = TRUE)
  # illu_col_name <- names(poly_illu)[2]  
  # illu_poly_valid <- poly_illu[[illu_col_name]][valid_poly_idx]
  
  output_list <- list(
    matrix = poly_valid,
    labels = label_poly_valid,
    # illu = illu_poly_valid,
    goodbands = good_bands_poly
    )
  
  return(output_list)
}


## For DP per pixel #################################################################

DP_per_pixel <- tsne_per_pixel(illu_crowns_out$area_filt, illu_crowns_out$label_filt, illu_crowns_out$area_crop)

crowns_filt_valid <- DP_per_pixel$matrix
illu_crowns_valid <- DP_per_pixel$illu
label_crowns_valid <- DP_per_pixel$labels


tsne_out <- Rtsne(seed=2, crowns_filt_valid)

tsne_plot <- data.frame(
  x = tsne_out$Y[,1],
  y = tsne_out$Y[,2],
  LianaFullCrown = as.factor(label_crowns_valid)
)

tsne_plot$illu <- illu_crowns_valid

ggplot(tsne_plot, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("0"= "olivedrab3", "1"="deeppink3","2"="purple4","3"= "orange3"),
                     labels = c("0" = "CR0", "1" = "CR1", "2"= "CR2", "3"= "CR3"))+
  labs(title= "tSNE of the HSI values of the drone photos dataset" )


ggplot(tsne_plot, aes(x = x, y = y, color = illu)) +
  geom_point() +
  scale_color_viridis_c(name = "Illumination") +
  labs(title = "tSNE of the HSI values of the DP dataset colored by illumination")

## For MD per pixel ###########################################################

MD_per_pixel <- tsne_per_pixel(illu_blob_out$area_filt, illu_blob_out$label_filt, illu_blob_out$area_crop)

blob_filt_valid <- MD_per_pixel$matrix
illu_blob_valid <- MD_per_pixel$illu
label_blob_valid <- MD_per_pixel$label


tsne_out_blob <- Rtsne(seed=2, blob_filt_valid)

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
  labs(title= "tSNE of the HSI values of the MD dataset" )

ggplot(tsne_plot_blob, aes(x = x, y = y, color = illu)) +
  geom_point() +
  scale_color_viridis_c(name = "Illumination") +
  labs(title = "tSNE of the HSI values of the MD dataset colored by illumination")


## For DP per polygon ################################################################

DP_per_poly <- tsne_per_polygon(illu_crowns_out$area_filt, crowns, "LianaFullC") #, illu)

crowns_poly_valid <- DP_per_poly$matrix
# illu_crowns_poly_valid <- DP_per_pixel$illu
label_crowns_poly_valid <- DP_per_poly$label

tsne_out_crowns_poly <- Rtsne(crowns_poly_valid, perplexity = 10, seed = 2)  #low perplexity due to few polygons

tsne_plot_crowns_poly <- data.frame(
  x = tsne_out_crowns_poly$Y[,1],
  y = tsne_out_crowns_poly$Y[,2],
  LianaFullCrown = as.factor(label_crowns_poly_valid)
)

# tsne_plot_crowns_poly$illu <- illu_crowns_poly_valid # /!\ illu_poly_valid_crowns VIDE !!!!

ggplot(tsne_plot_crowns_poly,aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+ #(size = 3)
  scale_color_manual(values = c("CR0"= "olivedrab3", "CR1"="deeppink3","CR2"="purple4","CR3"= "orange3"))+
  # scale_color_viridis_d(option="turbo")+
  labs(title= "tSNE of the mean spectrum per polygon of the DP dataset" )

# ggplot(tsne_plot_crowns_poly, aes(x = x, y = y, color = illu)) +
#   geom_point() +
#   scale_color_viridis_c(name = "Illumination") +
#   labs(title = "tSNE of of the mean spectrum per polygon of the MD dataset colored by illumination")



## For MD per polygon ###################################################################

MD_per_poly <- tsne_per_polygon(illu_blob_out$area_filt, blob, "Type")

blob_poly_valid <- MD_per_poly$matrix
#illu_blob_poly_valid <- MD_per_pixel$illu
label_blob_poly_valid <- MD_per_poly$label

# illu_poly_blob <- extract(illu, blob, fun = mean, na.rm = TRUE, ID = TRUE)
# illu_poly_valid_blob<- illu_poly_blob$illu[valid_poly_idx] # /!\ illu_poly_valid_crowns VIDE !!!!

tsne_out_blob_poly <- Rtsne(blob_poly_valid, perplexity = 10, seed = 2)  #low perplexity due to few polygons

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
  labs(title= "tSNE of the mean spectrum per polygon of the MD dataset" )


# UMAP #########################################################################

#For DP per pixel

umap_out_crowns <- umap2( seed = 2, crowns_filt_valid)

umap_plot_crowns <- data.frame(
  x = umap_out_crowns[,1],
  y = umap_out_crowns[,2],
  LianaFullCrown = as.factor(label_crowns_valid)
)

umap_plot_crowns$illu <- illu_crowns_valid
  
  
ggplot(umap_plot_crowns, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("0"= "olivedrab3", "1"="deeppink3","2"="purple4","3"= "orange3"),
                     labels = c("0" = "CR0", "1" = "CR1", "2"= "CR2", "3"= "CR3"))+
  # scale_color_viridis_d(option="turbo")+
  labs(title= "UMAP of the HSI values of the DP dataset" )

ggplot(umap_plot_crowns, aes(x = x, y = y, color = illu)) +
  geom_point() +
  scale_color_viridis_c(name = "Illumination") +
  labs(title = "UMAP of the HSI values of the DP dataset colored by illumination")

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
  labs(title= "UMAP of the HSI values of the MD dataset" )

ggplot(umap_plot_blob, aes(x = x, y = y, color = illu)) +
  geom_point() +
  scale_color_viridis_c(name = "Illumination") +
  labs(title = "UMAP of the HSI values of the MD dataset colored by illumination")

#For DP per polygon
umap_out_crowns_poly <- umap2( seed = 2, crowns_poly_valid)

umap_plot_crowns_poly <- data.frame(
  x = umap_out_crowns_poly[,1],
  y = umap_out_crowns_poly[,2],
  LianaFullCrown = as.factor(label_crowns_poly_valid)
)

ggplot(umap_plot_crowns_poly, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("CR0"= "olivedrab3", "CR1"="deeppink3","CR2"="purple4","CR3"= "orange3"))+
  labs(title= "UMAP of the mean HSI value per polygon of the DP dataset" )

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
  labs(title= "UMAP of the mean HSI value per polygon of the MD dataset" )

# Bhattacharyya distance ###################################################

##Distance between classes####

bhatt_dist_btwn <- function (matrix, matrix_label, label1, label2){

  #Reduce dimensionality
  max_ncp <- min(nrow(matrix) - 1, ncol(matrix))
  res_pca <- PCA(matrix, scale.unit = TRUE, ncp = max_ncp, graph = FALSE)
  cum_var <- res_pca$eig[, "cumulative percentage of variance"]
  n_components <- which(cum_var >= 85)[1]
  
  matrix_reduced <- res_pca$ind$coord[, 1:n_components]
  
  
  #Separate groups
  group0 <- matrix_reduced[matrix_label == label1, ]
  group1 <- matrix_reduced[matrix_label == label2, ]
  
  
  # Test if the distribution is Gaussian first
  shapiro_group0 <- apply(group0, 2, function(x) shapiro.test(x)$p.value)
  shapiro_group1 <- apply(group1, 2, function(x) shapiro.test(x)$p.value)
  
  par(mfrow = c(2, 3))
  for (i in 1:5) {
    qqnorm(group0[, i], main = paste("Composante", i, "- groupe 0"))
    qqline(group0[, i])
  }
  
  par(mfrow = c(2, 3))
  for (i in 1:5) {
    qqnorm(group1[, i], main = paste("Composante", i, "- groupe 1"))
    qqline(group1[, i])
  }
  
  # Calculate mean vectors
  mu0 <- colMeans(group0)
  mu1 <- colMeans(group1)
  
  # Calculate covariance matrix
  sigma0 <- cov(group0)
  sigma1 <- cov(group1)
  
  dist_bhatt <- bhattacharyya.dist(mu0, mu1, sigma0, sigma1)
  
  coeff_bhatt <- exp(-dist_bhatt)
  
  output <- list(
    mean0 = mu0,
    mean1 = mu1,
    dist_bhatt = dist_bhatt,
    coeff_bhatt = coeff_bhatt,
    shapiro_group0 = shapiro_group0,
    shapiro_group1 = shapiro_group1,
    n_components = n_components
  )
  
  return(output) 
}

### With DP dataset ####
separability_dp_within <- bhatt_dist_btwn(crowns_poly_valid, label_crowns_poly_valid, "CR0", "CR3")

#Check if we keep the Bhattacharya distance
separability_dp_within$shapiro_group0
separability_dp_within$shapiro_group0 < 0.05

separability_dp_within$shapiro_group1
separability_dp_within$shapiro_group1 < 0.05

#Results
bhatt_dist_dp_LianaVSTrees <- separability_dp_within$dist_bhatt
bhatt_dist_dp_LianaVSTrees

bhatt_coeff_dp_LianaVSTrees <-separability_dp_within$coeff_bhatt
bhatt_coeff_dp_LianaVSTrees

### With MD dataset ####
separability_md_within <- bhatt_dist_btwn(blob_poly_valid, label_blob_poly_valid, 0, 1)

separability_md_within$shapiro_group0
separability_md_within$shapiro_group0<0.05

separability_md_within$shapiro_group1
separability_md_within$shapiro_group1<0.05

bhatt_dist_dp_LianaVSTrees <- separability_md_within$bhatt
bhatt_dist_dp_LianaVSTrees


bhatt_dist_within <- function(area_filtered, area_vector, label_name, target_label){
  
  #Get unique id for each polygon
  area_vector$poly_id <- seq_len(nrow(area_vector))
  
  poly_id_r <- rasterize(area_vector, area_filtered, field = "poly_id")
  
  area_matrix <- as.matrix(area_filtered)
  poly_id_m <- as.matrix(poly_id_r)[, 1]
  class_by_id <- values(area_vector)[[label_name]]
  
  na_per_band <- colSums(is.na(area_matrix))
  good_bands <- na_per_band < nrow(area_matrix)
  
  valid_idx <- which(!is.na(poly_id_m) & complete.cases(area_matrix[, good_bands]))
  
  pixel_matrix <- area_matrix[valid_idx, good_bands]
  pixel_poly_id <- poly_id_m[valid_idx]
  
  # Keep only the pixels of the target label
  keep_ids <- which(class_by_id == target_label)
  keep_px <- pixel_poly_id %in% keep_ids
  
  pixel_matrix <- pixel_matrix[keep_px, ]
  pixel_poly_id <- pixel_poly_id[keep_px]
  
  # Dimensionality reduction
  max_ncp <- min(nrow(pixel_matrix) - 1, ncol(pixel_matrix))
  res_pca <- PCA(pixel_matrix, scale.unit = TRUE, ncp = max_ncp, graph = FALSE)
  cum_var <- res_pca$eig[, "cumulative percentage of variance"]
  n_components <- which(cum_var >= 85)[1]
  
  pixel_reduced <- res_pca$ind$coord[, 1:n_components]
  
  # Calculate mu & sigma for each polygon
  polygon_ids <- unique(pixel_poly_id)
  
  polygon_stats <- lapply(polygon_ids, function(id) {
    px <- pixel_reduced[pixel_poly_id == id, , drop = FALSE]
    if (nrow(px) <= n_components) return(NULL)  # pas assez de pixels pour estimer Sigma
    list(mu = colMeans(px), Sigma = cov(px), n = nrow(px))
  })
  names(polygon_stats) <- polygon_ids
  polygon_stats <- Filter(Negate(is.null), polygon_stats)
  
  ids <- names(polygon_stats)
  if (length(ids) < 2) {
    warning("Pas assez de polygones valides pour calculer des distances intra-classe.")
    return(NULL)
  }
  
  # Distance by pairs
  pairs <- combn(ids, 2, simplify = FALSE)
  
  dist_df <- do.call(rbind, lapply(pairs, function(p) {
    s1 <- polygon_stats[[p[1]]]
    s2 <- polygon_stats[[p[2]]]
    d <- tryCatch(
      bhattacharyya.dist(s1$mu, s2$mu, s1$Sigma, s2$Sigma),
      error = function(e) NA
    )
    data.frame(poly1 = p[1], poly2 = p[2], dist_bhatt = d, coeff_bhatt = exp(-d))
  }))
  
  output <- list(
    distances = dist_df,
    n_polygons = length(ids),
    n_components = n_components,
    mean_dist = mean(dist_df$dist_bhatt, na.rm = TRUE),
    mean_coeff = mean(dist_df$coeff_bhatt, na.rm = TRUE)
  )
  
  return(output)
}

#Results

bhatt_within_dp_CR0 <- bhatt_dist_within(illu_crowns_out$area_filt, crowns, "LianaFullC", "CR0")
bhatt_within_dp_CR0$mean_dist
bhatt_within_dp_CR0$mean_coeff
sum(is.na(bhatt_within_dp_CR0$distances$dist_bhatt))   # number of failed pairs (singular covariance, etc.)


bhatt_within_dp_CR3 <- bhatt_dist_within(illu_crowns_out$area_filt, crowns, "LianaFullC", "CR3")
bhatt_within_dp_CR3$mean_dist
bhatt_within_dp_CR3$mean_coeff
sum(is.na(bhatt_within_dp_CR0$distances$dist_bhatt))

bhatt_within_md_0 <- bhatt_dist_within(illu_blob_out$area_filt, blob, "Type", "0")
bhatt_within_md_0$mean_dist
bhatt_within_md_0$mean_coeff
sum(is.na(bhatt_within_md_0$distances$dist_bhatt))   

bhatt_within_md_1 <- bhatt_dist_within(illu_blob_out$area_filt, blob, "Type", "1")
bhatt_within_md_1$mean_dist
bhatt_within_md_1$mean_coeff
sum(is.na(bhatt_within_md_1$distances$dist_bhatt)) 


# Mean spectrum comparison ###################################################

mean_spectrum_by_group <- function(matrix_valid, labels_valid, wavelengths, dataset_name) {
  groups <- unique(labels_valid)
  
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

wavelengths_blob <- as.numeric(names(illu_blob_out$area_filt))[MD_per_poly$goodbands]
wavelengths_crowns <- as.numeric(names(illu_crowns_out$area_filt))[DP_per_poly$goodbands]

spec_crowns <- mean_spectrum_by_group(crowns_poly_valid, label_crowns_poly_valid, wavelengths_crowns, "DP")
spec_blob   <- mean_spectrum_by_group(blob_poly_valid, label_blob_poly_valid, wavelengths_blob, "MD")

all_spectra <- rbind(spec_crowns, spec_blob)

selected <- subset(all_spectra, group %in% c("DP_CR0", "DP_CR3", "MD_0", "MD_1"))

label_config <- data.frame(
  group = c("DP_CR0", "DP_CR1", "DP_CR2", "DP_CR3", "MD_0", "MD_1"),
  display_name = c("CR0 ( No Liana DP)", "CR1 (DP)", "CR2", "CR3 (Liana DP)", "No liana (MD)", "Liana (MD)"),
  color = c("yellowgreen", "darkorange1", "purple4", "deeppink3", "cyan3", "pink")
)

ggplot(selected, aes(x = wavelength, y = reflectance, color = group)) +
  geom_line(linewidth = 1) +
  labs(x = "Wavelength (nm)", y = "Mean reflectance", 
       title = "Mean spectrum by dataset group (DP et MD)", color = "Type") +
  scale_color_manual(
    values = setNames(label_config$color, label_config$group),
    labels = setNames(label_config$display_name, label_config$group)
    )

