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
#plot(area_blob[[3]]) # the number is the rank of the band not it's value

### Taking into account the liana label ####
label_blob <- rasterize(blob, area_blob, field = "Type")

## Crowns from drone photos (DP) ####
crowns <- vect("Couronnes/Crowns_Liana_2025_GP.shp")
crowns <- project(crowns, crs(r))

area_crowns <- crop(r, crowns, mask=TRUE)
#plot(area_crowns[[53]]) 

label_liana <- rasterize(crowns, area_crowns, field = "LianaFullC")


# Removing pixel based on illu ###########################################

#For DP crowns ####

illu_crop_crowns <- crop(illu, area_crowns)

#Number of valid pixels before taking illu into acount
n_before <- global(!is.na(area_crowns[[1]]), "sum")$sum
n_before

# Using a mask keep all pixels with illu superior to the threshold
illu_mask_050_crowns <- illu_crop_crowns > 0.05 #change the threshold based on QGIS visualization and classification results

area_crowns_filt_050 <- mask(area_crowns, illu_mask_050_crowns, maskvalues = c(FALSE, NA))
label_liana_filt_050 <- mask(label_liana, illu_mask_050_crowns, maskvalues = c(FALSE, NA))

#Number of valid pixels after taking illu into acount
n_after <- global(!is.na(area_crowns_filt_050[[1]]), "sum")$sum
n_after


#For MD crowns ####
illu_crop_blob <- crop(illu, area_blob)

illu_mask_050_blob <- illu_crop_blob > 0.05 #change the threshold based on QGIS visualization and classification results

area_blob_filt_050 <- mask(area_blob, illu_mask_050_blob, maskvalues = c(FALSE, NA))
label_blob_filt_050 <- mask(label_blob, illu_mask_050_blob, maskvalues = c(FALSE, NA))

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
  # dim(crowns_050_valid)
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

DP_per_pixel <- tsne_per_pixel(area_crowns_filt_050, label_liana_filt_050, illu_crop_crowns)

crowns_050_valid <- DP_per_pixel$matrix
illu_crowns_valid <- DP_per_pixel$illu
label_crowns_valid <- DP_per_pixel$labels

### tSNE ####

tsne_out <- Rtsne(seed=2, crowns_050_valid)

### tSNE plot ####
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

## tSNE for MD per pixel ###########################################################

### Creation of the correct matrix ####

MD_per_pixel <- tsne_per_pixel(area_blob_filt_050, label_blob_filt_050, illu_crop_blob)

blob_050_valid <- MD_per_pixel$matrix
illu_blob_valid <- MD_per_pixel$illu
label_blob_valid <- MD_per_pixel$label

### tSNE ####

tsne_out_blob <- Rtsne(seed=2, blob_050_valid)

### tSNE plot ####
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


## For DP per polygon (crown) ################################################################

DP_per_poly <- tsne_per_polygon(area_crowns_filt_050, crowns, "LianaFullC") #, illu)

crowns_poly_valid <- DP_per_poly$matrix
# illu_crowns_poly_valid <- DP_per_pixel$illu
label_crowns_poly_valid <- DP_per_poly$label

### tSNE
tsne_out_crowns_poly <- Rtsne(crowns_poly_valid, perplexity = 10, seed = 2)  #low perplexity due to few polygons

### tSNE plot
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

MD_per_poly <- tsne_per_polygon(area_blob_filt_050, blob, "Type")

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

umap_out_crowns <- umap2( seed = 2, crowns_050_valid)

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
umap_out_blob <- umap2( seed = 2, blob_050_valid)

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
  res_pca <- PCA(matrix, scale.unit = TRUE, graph = FALSE)
  cum_var <- res_pca$eig[, "cumulative percentage of variance"]
  n_components <- which(cum_var >= 85)[1]
  
  matrix_reduced <- res_pca$ind$coord[, 1:n_components]
  
  #Separate groups
  group0 <- matrix_reduced[matrix_label == label1, ]
  group1 <- matrix_reduced[matrix_label == label2, ]
  
  # Calculate mean vectors
  mu0 <- colMeans(group0)
  mu1 <- colMeans(group1)
  
  # Calculate covariance matrix
  sigma0 <- cov(group0)
  sigma1 <- cov(group1)
  
  dist_bhatt <- bhattacharyya.dist(mu0, mu1, sigma0, sigma1)
  
  output <- list(
    mean0 = mu0,
    mean1 = mu1,
    bhatt = dist_bhatt
  )
  
  return(output) 
}

### With DP dataset ####
separability_dp_within <- bhatt_dist_btwn(crowns_poly_valid, label_crowns_poly_valid, "CR0", "CR3")

bhatt_dist_dp_LianaVSTrees <- separability_dp_within$bhatt
bhatt_dist_dp_LianaVSTrees


### With MD dataset ####
separability_md_within <- bhatt_dist_btwn(blob_poly_valid, label_blob_poly_valid, 0, 1)

bhatt_dist_dp_LianaVSTrees <- separability_md_within$bhatt
bhatt_dist_dp_LianaVSTrees

## Distance within classes #### 

bhatt_dist_within <- function(shapefile, mask_illu, matrix_area_filtered, good_bands){

  shapefile$poly_ID <- seq_len(nrow(shapefile))

  label_poly_ID <- rasterize(shapefile, area, field = "poly_ID")
  label_poly_ID_filt <- mask (label_poly_ID, mask_illu, maskvalues = c(FALSE, NA))
  label_poly_ID_m <- as.matrix(label_poly_ID_filt)[, 1]

  valid_idx_poly_ID <- which(!is.na(label_poly_ID_m) & complete.cases(matrix_area_filtered[, good_bands]))

  pixel_matrix <- matrix_area_filtered[valid_idx_poly_ID, good_bands]
  pixel_ID <- label_poly_ID_m[valid_idx_poly_ID]

  res_pca <- PCA(pixel_matrix, scale.unit = TRUE, graph = FALSE)

  pca <- FactoMineR::PCA(pixel_matrix, scale.unit = TRUE, graph = FALSE)
  cum_var <- pca$eig[, "cumulative percentage of variance"]
  n_components <- which(cum_var >= 85)[1]

  pixel_reduced <- res_pca$ind$coord[, 1:n_components]
}

## Mean spectrum comparison ###################################################
# noliana_dp =crowns_050_valid[label_crowns_valid == "CR0", ]
# 
# 
# mean_noliana_dp = separability_dp_within$mean0
# mean_liana_dp = separability_dp_within$mean1
# 
# mean_noliana_md = separability_md_within$mean0
# mean_liana_md = separability_md_within$mean1

# pixel_data_crowns <- extract(area_crowns_filt_050, crowns, ID = TRUE)

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

wavelengths_blob <- as.numeric(names(area_blob_filt_050))[MD_per_poly$goodbands]
wavelengths_crowns <- as.numeric(names(area_crowns_filt_050))[DP_per_poly$goodbands]

spec_crowns <- mean_spectrum_by_group(crowns_poly_valid, label_crowns_poly_valid, wavelengths_crowns, "DP")
spec_blob   <- mean_spectrum_by_group(blob_poly_valid, label_blob_poly_valid, wavelengths_blob, "MD")

all_spectra <- rbind(spec_crowns, spec_blob)

# Garder uniquement CR0, CR3 (DP) et 0, 1 (MD)
selected <- subset(all_spectra, group %in% c("DP_0", "DP_3", "MD_0", "MD_1"))

ggplot(selected, aes(x = wavelength, y = reflectance, color = group)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Wavelength (nm)",
    y = "Mean reflectance",
    title = "Spectres moyens par groupe (DP et MD)",
    color = "Groupe"
  )
