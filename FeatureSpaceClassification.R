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

# Uploading data and cropping ##########################################
area_blob <- rast("Classification_R/LianaClassification/FeatureSpaceData/crowns_feat.tiff")
area_crowns <- rast("Classification_R/LianaClassification/FeatureSpaceData/blob_feat.tiff")


blob <- vect("Lianes/LianaOrNotLiana_55_45.shp")
blob <- project(blob, crs(area_blob))

crowns <- vect("Couronnes/Crowns_Liana_2025_GP_TRUE.shp")
crowns <- project(crowns, crs(area_crowns))


label_blob <- rasterize(blob, area_blob, field = "Type")
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

# Pixel filtering done in Python before retrieving the feature space

# Dimensionality reduction and data visualization ###############################

## Data prep ##############################################################

data_prep_per_pixel <- function(area_rast, label, max_na_pct = 0.2){
  area_matrix <- as.matrix(area_rast)
  label_matrix <- as.matrix (label)[,1]
  
  valid_idx <- which(!is.na(label_matrix) & rowMeans(is.na(area_matrix)) <= max_na_pct)
  
  area_valid <- area_matrix[valid_idx, , drop = FALSE]
  labels_valid <- label_matrix[valid_idx]
  
  output_list <- list(
    matrix = area_valid,
    labels = labels_valid,
    valid_idx = valid_idx)
  
  return(output_list)
} 

data_prep_per_poly <- function(area_rast, area_vector, label_name, max_na_pct = 0.2){
  
  # Mean spectrum per polygon
  poly_spectra <- extract(area_rast, area_vector, fun = mean, na.rm = TRUE, ID = TRUE)
  poly_labels <- values(area_vector)[[label_name]][poly_spectra$ID]
  
  # Create the matrix
  poly_matrix <- as.matrix(poly_spectra[, -1])
  
  na_prop <- rowMeans(is.na(poly_matrix))
  
  valid_idx <- na_prop <= max_na_pct
  poly_valid <- poly_matrix[valid_idx, ]
  label_poly_valid <- poly_labels[valid_idx]
  
  good_bands_poly <- rep(TRUE, ncol(poly_matrix))
  
  
  output_list <- list(
    matrix = poly_valid,
    labels = label_poly_valid,
    goodbands = good_bands_poly
  )
  
  return(output_list)
}

### Data prep for DP
DP_per_pixel <- data_prep_per_pixel(area_crowns, label_crowns)

crowns_valid <- DP_per_pixel$matrix
label_crowns_valid <- DP_per_pixel$labels


DP_per_poly <- data_prep_per_poly(area_crowns, crowns, "LianaFullC")

crowns_poly_valid <- DP_per_poly$matrix
label_crowns_poly_valid <- DP_per_poly$label


### Data prep for MD
MD_per_pixel <- data_prep_per_pixel(area_blob, label_blob)

blob_valid <- MD_per_pixel$matrix
label_blob_valid <- MD_per_pixel$label


MD_per_poly <- data_prep_per_poly(area_blob, blob, "Type")

blob_poly_valid <- MD_per_poly$matrix
label_blob_poly_valid <- MD_per_poly$label

## Proportion of liana after data_prep #####################################
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
pca_out_crowns <- PCA(crowns_valid, scale.unit = TRUE, graph = FALSE)
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
pca_out_blob <- PCA(blob_valid, scale.unit = TRUE, graph = FALSE)
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
perplexity_crown = (nrow(crowns_valid)-1)/3
perplexity_crown
set.seed(2)
tsne_out_crowns <- Rtsne(crowns_valid, check_duplicates = FALSE, perplexity = perplexity_crown )

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
# perplexity_blob = (nrow(blob_valid)-1)/3 - 1
# perplexity_blob
set.seed(2)
tsne_out_blob <- Rtsne( blob_valid, check_duplicates = FALSE)

tsne_plot_blob <- data.frame(
  x = tsne_out_blob$Y[,1],
  y = tsne_out_blob$Y[,2],
  Type = as.factor(label_blob_valid)
)

ggplot(tsne_plot_blob, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#666", "1" = "#001312"),
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
    values = c("0" = "#666", "1" = "#001312"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "tSNE of the mean spectrum per polygon of the MD dataset (Illu filter)" )+
  theme_minimal()


## UMAP #########################################################################

#For DP per pixel

umap_out_crowns <- umap2( seed = 2, crowns_valid,  n_neighbors = min(15, nrow(crowns_valid) - 1))

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
umap_out_blob <- umap2( seed = 2, blob_valid)

umap_plot_blob <- data.frame(
  x = umap_out_blob[,1],
  y = umap_out_blob[,2],
  Type = as.factor(label_blob_valid)
)


ggplot(umap_plot_blob, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#666", "1" = "#001312"),
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
    values = c("0" = "#666", "1" = "#001312"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "UMAP of the mean HSI value per polygon of the MD dataset (Illu filter)" )+
  theme_minimal()

