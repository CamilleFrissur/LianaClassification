# File initialization #############################################

rm(list = ls())

## Import libraries #### 
library(tidyverse)
library(sf)
library(terra)
library(Rtsne)
library(ggplot2)
library(uwot)
library(FactoMineR)
library(fpc)

## Set working Directory ####
setwd ("D:/DataCamille/")


# Uploading data and cropping ##########################################
#r = rast("E:/DataCamille/HSI/HSI_2016_raw/NOU_FL12_20160921_182748_PS01_IMG008_atm_slice_002_x2.hyspex")

## HSI mosaic ####
r <- rast("HSI/mosaic_2016/nou_20160921_l1c_mosaic.bsq")
nlyr(r) #number of spectral bands

## Illu ####
illu <- rast("HSI/mosaic_2016/nou_20160921_ilu_mosaic.bsq")
nlyr(illu)
hist(illu, breaks = 100)
plot(illu)

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

## For DP per pixel #################################################################

### Creation of the correct matrix ####

crowns_050_m  <- as.matrix(area_crowns_filt_050)
liana_050_m <- as.matrix(label_liana_filt_050)[, 1]

na_per_band <- colSums(is.na(crowns_050_m)) # to apply before in the next version
good_bands <- na_per_band < nrow(crowns_050_m)
sum(good_bands)

valid_idx <- which(!is.na(liana_050_m) & stats::complete.cases(crowns_050_m[, good_bands]))
length(valid_idx)

crowns_050_valid <- crowns_050_m[valid_idx, good_bands]
labels_valid <- liana_050_m[valid_idx]

dim(crowns_050_valid)
length(labels_valid)

### tSNE ####

tsne_out <- Rtsne(seed=2, crowns_050_valid)

### tSNE plot ####
tsne_plot <- data.frame(
  x = tsne_out$Y[,1],
  y = tsne_out$Y[,2],
  LianaFullCrown = as.factor(labels_valid)
)

ggplot(tsne_plot, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  labs(title= "tSNE of the HSI values of the drone photos dataset" )

## tSNE for MD per pixel ###########################################################

### Creation of the correct matrix ####

blob_050_m  <- as.matrix(area_blob_filt_050)
type_050_m <- as.matrix(label_blob_filt_050)[, 1]

na_per_band_blob <- colSums(is.na(blob_050_m)) # to apply before in the next version
good_bands_blob <- na_per_band_blob < nrow(blob_050_m)
sum(good_bands_blob)

valid_idx_blob <- which(!is.na(type_050_m) & stats::complete.cases(blob_050_m[, good_bands_blob]))
length(valid_idx_blob)

blob_050_valid <- blob_050_m[valid_idx_blob, good_bands_blob]
labels_valid_blob <- type_050_m[valid_idx_blob]

dim(blob_050_valid)
length(labels_valid_blob)

### tSNE ####

tsne_out_blob <- Rtsne(seed=2, blob_050_valid)

### tSNE plot ####
tsne_plot_blob <- data.frame(
  x = tsne_out_blob$Y[,1],
  y = tsne_out_blob$Y[,2],
  Type = as.factor(labels_valid_blob)
)

ggplot(tsne_plot_blob, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "olivedrab3", "1" = "deeppink2"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "tSNE of the HSI values of the MD dataset" )


## For DP per polygon (crown) ################################################################

# Mean spectrum per polygon
crowns_poly_spectra <- extract(area_crowns_filt_050, crowns, fun = mean, na.rm = TRUE, ID = TRUE)
crowns_poly_labels <- crowns$LianaFullC[crowns_poly_spectra$ID]

# Create the matrix
crowns_poly_matrix <- as.matrix(crowns_poly_spectra[, -1])

na_per_band_crowns_poly <- colSums(is.na(crowns_poly_matrix))
good_bands_crowns_poly <- na_per_band_crowns_poly < nrow(crowns_poly_matrix)

crowns_poly_matrix_valid <- crowns_poly_matrix[, good_bands_crowns_poly]

# Remove polygons full of NA
valid_crowns_poly_idx <- stats::complete.cases(crowns_poly_matrix_valid)
crowns_poly_matrix_valid <- crowns_poly_matrix_valid[valid_crowns_poly_idx, ]
crowns_poly_labels_valid <- crowns_poly_labels[valid_crowns_poly_idx]

dim(crowns_poly_matrix_valid)
length(crowns_poly_labels_valid)

tsne_out_crowns_poly <- Rtsne(crowns_poly_matrix_valid, perplexity = 10, seed = 2)  #low perplexity due to few polygons

tsne_plot_crowns_poly <- data.frame(
  x = tsne_out_crowns_poly$Y[,1],
  y = tsne_out_crowns_poly$Y[,2],
  LianaFullCrown = as.factor(crowns_poly_labels_valid)
)

ggplot(tsne_plot_crowns_poly,aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+ #(size = 3)
  labs(title= "tSNE of the mean spectrum per polygon of the MD dataset" )


## For MD per polygon ###################################################################

# Mean spectrum per polygon
poly_spectra <- extract(area_blob_filt_050, blob, fun = mean, na.rm = TRUE, ID = TRUE)
poly_labels <- blob$Type[poly_spectra$ID]

# Create the matrix
poly_matrix <- as.matrix(poly_spectra[, -1])

na_per_band_poly <- colSums(is.na(poly_matrix))
good_bands_poly <- na_per_band_poly < nrow(poly_matrix)

poly_matrix_valid <- poly_matrix[, good_bands_poly]

# Remove polygons full of NA
valid_poly_idx <- stats::complete.cases(poly_matrix_valid)
poly_matrix_valid <- poly_matrix_valid[valid_poly_idx, ]
poly_labels_valid <- poly_labels[valid_poly_idx]

dim(poly_matrix_valid)
length(poly_labels_valid)

tsne_out_poly <- Rtsne(poly_matrix_valid, perplexity = 10, seed = 2)  #low perplexity due to few polygons

tsne_plot_poly <- data.frame(
  x = tsne_out_poly$Y[,1],
  y = tsne_out_poly$Y[,2],
  Type = as.factor(poly_labels_valid)
)

ggplot(tsne_plot_poly, aes(x = x, y = y, color = Type)) +
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
  LianaFullCrown = as.factor(labels_valid)
)

ggplot(umap_plot_crowns, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  labs(title= "UMAP of the HSI values of the DP dataset" )

#For MD per pixel
umap_out_blob <- umap2( seed = 2, blob_050_valid)

umap_plot_blob <- data.frame(
  x = umap_out_blob[,1],
  y = umap_out_blob[,2],
  Type = as.factor(labels_valid_blob)
)

ggplot(umap_plot_blob, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "olivedrab3", "1" = "deeppink2"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "UMAP of the HSI values of the MD dataset" )

#For DP per polygon
umap_out_crowns_poly <- umap2( seed = 2, crowns_poly_matrix_valid)

umap_plot_crowns_poly <- data.frame(
  x = umap_out_crowns_poly[,1],
  y = umap_out_crowns_poly[,2],
  LianaFullCrown = as.factor(crowns_poly_labels_valid)
)

ggplot(umap_plot_crowns_poly, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  labs(title= "UMAP of the mean HSI value per polygon of the DP dataset" )

#For MD per polygon
umap_out_poly <- umap2( seed = 2, poly_matrix_valid)

umap_plot_poly <- data.frame(
  x = umap_out_poly[,1],
  y = umap_out_poly[,2],
  Type = as.factor(poly_labels_valid)
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

### With DP dataset ####



### With MD dataset ####

#Reduce dimensionality
res_pca <- PCA(poly_matrix_valid, scale.unit = TRUE, graph = FALSE)

res_pca$eig[, "cumulative percentage of variance"]

n_components <- 9  # ajustez selon le seuil choisi
poly_matrix_reduced <- res_pca$ind$coord[, 1:n_components]

#Separate groups
group0 <- poly_matrix_reduced[poly_labels_valid == 0, ]
group1 <- poly_matrix_reduced[poly_labels_valid == 1, ]


# Calculate mean vectors
mu0 <- colMeans(group0)
mu1 <- colMeans(group1)

# Calculate covariance matrix
sigma0 <- cov(group0)
sigma1 <- cov(group1)

dist_bhatt <- bhattacharyya.dist(mu0, mu1, sigma0, sigma1)
dist_bhatt

## Distance within classes #### 