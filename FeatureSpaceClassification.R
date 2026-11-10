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

source("Classification_R/LianaClassification/Functions.R")

# Uploading data and cropping ##########################################
area_blob <- rast("Classification_R/LianaClassification/FeatureSpaceData/blob_feat.tiff")
area_crowns <- rast("Classification_R/LianaClassification/FeatureSpaceData/crowns_feat.tiff")


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
  theme_minimal()+
  ylim(-20,20)


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
  theme_minimal()+
  ylim(-15,15)

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

# Classification ##############################################################

# ## Linear Discriminant Analysis (LDA) #########################################
# 
# 
# ## Within DP dataset
# crowns_df <- as.data.frame(crowns_valid)
# crowns_df$ Liana <- as.factor(label_crowns_valid)
# 
# set.seed(123)
# crowns_df[,1:ncol(crowns_valid)]<- scale(crowns_df[,1:ncol(crowns_valid)])
# sample_crowns <- sample(c(TRUE, FALSE), nrow(crowns_df), replace=TRUE, prob=c(0.7,0.3))
# train_crowns <- crowns_df[sample_crowns, ]
# test_crowns <- crowns_df[!sample_crowns, ] 
# 
# n_classes <- length(levels(train_crowns$Liana))
# model_crowns <- lda(Liana ~ ., data = train_crowns, prior = rep(1/n_classes, n_classes))
# predictions_crowns <- predict(model_crowns, test_crowns)
# lda_crowns_accuracy <- mean(predictions_crowns$class==test_crowns$Liana)
# lda_crowns_accuracy 
# 
# ### Stratified folded cross-validation
# predictors <- crowns_df[ , -which(names(crowns_df) %in% c("Liana"))]
# 
# ctrl_crowns <- trainControl(
#   method = "cv",
#   number = 5,        
#   savePredictions = "final"
# ) 
# 
# length(levels(crowns_df$Liana))
# model_cv_crowns <- train(
#   x = predictors,
#   y = crowns_df$Liana,
#   method = "lda",
#   trControl = ctrl_crowns, 
#   prior = rep(1/n_classes, n_classes)
# )
# 
# cv_lda_crowns_accuracy <-model_cv_crowns$results$Accuracy 
# cv_lda_crowns_accuracy 
# 
# # Predicted liana proportion
# prop.table(table(model_cv_crowns$pred$pred))
# 
# ## Within MD dataset
# blob_df <- as.data.frame(blob_valid)
# blob_df$ Liana <- as.factor(label_blob_valid)
# 
# set.seed(123)
# blob_df[,1:ncol(crowns_valid)]<- scale(blob_df[,1:ncol(crowns_valid)]) #ncol(crowns_valid) should be equal to ncol(blob_valid)
# sample_blob <- sample(c(TRUE, FALSE), nrow(blob_df), replace=TRUE, prob=c(0.7,0.3))
# train_blob <- blob_df[sample_blob, ]
# test_blob <- blob_df[!sample_blob, ] 
# 
# n_classes <- length(levels(train_blob$Liana))
# model_blob <- lda(Liana ~ ., data = train_blob, prior = rep(1/n_classes, n_classes))
# predictions_blob <- predict(model_blob, test_blob)
# lda_blob_accuracy <- mean(predictions_blob$class==test_blob$Liana)
# lda_blob_accuracy #0.8671278
# 
# ### Stratified folded cross-validation
# predictors <- blob_df[ , -which(names(blob_df) %in% c("Liana"))]
# 
# ctrl_blob <- trainControl(
#   method = "cv",
#   number = 5,
#   savePredictions = "final"
# ) 
# 
# n_classes <-  length(levels(blob_df$Liana))
# model_cv_blob <- train(
#   x = predictors,
#   y = blob_df$Liana,
#   method = "lda",
#   trControl = ctrl_blob,
#   prior = rep(1/n_classes, n_classes)
# )
# 
# cv_lda_blob_accuracy <-model_cv_blob$results$Accuracy
# cv_lda_blob_accuracy 
# 
# # Predicted liana proportion
# prop.table(table(model_cv_blob$pred$pred))
# 
# ## Between MD and DP datasets
# crowns_df <- as.data.frame(crowns_valid)
# crowns_df$Liana <- as.character(label_crowns_valid)
# crowns_df <- crowns_df[crowns_df$Liana %in% c('0','3'), ]
# crowns_df$Liana[crowns_df$Liana == '3'] <- '1'
# crowns_df$Liana <- as.factor(crowns_df$Liana)
# type_crowns <- rep(c('crowns'), nrow(crowns_df))
# crowns_df$Type <- as.factor(type_crowns)
# 
# 
# blob_df <- as.data.frame(blob_valid)
# blob_df$ Liana <- as.factor(label_blob_valid)
# type_blob <- rep(c('blob'), nrow(blob_df))
# blob_df$Type <- as.factor(type_blob)
# 
# all_df <- rbind(crowns_df, blob_df)
# 
# set.seed(123)
# all_df[,1:ncol(crowns_valid)]<- scale(all_df[,1:ncol(crowns_valid)])
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
# 
# ## LDA trained on pure labels to check % of liana predicted ##################
# crowns_df_pure <- as.data.frame(crowns_valid)
# crowns_df_pure$LianaOriginal <- as.factor(label_crowns_valid)
# 
# crowns_df_pure$Liana <- as.character(label_crowns_valid)
# crowns_df_pure$Liana[crowns_df_pure$Liana %in% c("0", "1")] <- "0"
# crowns_df_pure$Liana[crowns_df_pure$Liana %in% c("2", "3")] <- "1"
# crowns_df_pure$Liana <- as.factor(crowns_df_pure$Liana)
# 
# crowns_df_pure$Type <- factor("crowns", levels = c("crowns", "blob"))
# crowns$poly_id <- seq_len(nrow(crowns))
# poly_id_r <- rasterize(crowns, area_crowns, field = "poly_id")
# poly_id_m <- as.matrix(poly_id_r)[, 1]
# poly_id_valid <- poly_id_m[DP_per_pixel$valid_idx]
# crowns_df_pure$poly_id <- ifelse(is.na(poly_id_valid), 0, poly_id_valid)
# 
# blob_df_pure <- as.data.frame(blob_valid)
# blob_df_pure$Liana <- as.factor(label_blob_valid)
# blob_df_pure$Type <- factor("blob", levels = c("crowns", "blob"))
# blob_df_pure$LianaOriginal <- NA
# blob_df_pure$poly_id <- 0
# 
# common_cols <- intersect(names(crowns_df_pure), names(blob_df_pure))
# crowns_df_pure <- crowns_df_pure[, common_cols]
# blob_df_pure <- blob_df_pure[, common_cols]
# 
# all_df <- rbind(crowns_df_pure, blob_df_pure)
# 
# ### Trained on blob / Tested on crowns ####
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
# 
# ### StratifiedGroupKFold ####
# 
# stratified_group_kfold <- function(data, group_col, stratify_col, k = 5, seed = 123) {
#   set.seed(seed)
#   group_classes <- sapply(unique(data[[group_col]]), function(group) {
#     group_data <- data[data[[group_col]] == group, ]
#     names(which.max(table(group_data[[stratify_col]])))
#   })
#   groups_df <- data.frame(group = unique(data[[group_col]]), class = group_classes)
#   folds <- vector("list", k)
#   for (class in unique(groups_df$class)) {
#     class_groups <- groups_df$group[groups_df$class == class]
#     fold_assignments <- sample(rep(1:k, length.out = length(class_groups)))
#     for (i in 1:k) folds[[i]] <- c(folds[[i]], class_groups[fold_assignments == i])
#   }
#   lapply(folds, function(groups) which(data[[group_col]] %in% groups))
# }
# 
# #### For crowns
# set.seed(123)
# folds_crowns <- stratified_group_kfold(crowns_df_pure, "poly_id", "Liana", k = 5)
# 
# cv_results <- list()
# for (i in 1:length(folds_crowns)) {  
#   test_indices <- folds_crowns[[i]]  
#   train_indices <- setdiff(1:nrow(crowns_df_pure), test_indices) 
#   
#   train_data <- crowns_df_pure[train_indices, ]
#   test_data  <- crowns_df_pure[test_indices, ]
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
# 
# cv_accuracies <- sapply(cv_results, function(x) x$accuracy)
# cat("Accuracies per fold:", paste(round(cv_accuracies, 3), collapse = ", "), "\n")
# cat("Mean accuracy CV:", round(mean(cv_accuracies), 3), "\n")
# 
# all_cv_results <- do.call(rbind, lapply(cv_results, function(x) {
#   n_test <- nrow(x$test_data)
#   data.frame(
#     fold = rep(x$fold, n_test),
#     poly_id = x$test_data$poly_id,
#     LianaOriginal = x$test_data$LianaOriginal,
#     correct = x$predictions$class == x$test_data$Liana,
#     predicted = x$predictions$class
#   )
# }))
# 
# all_cv_results_filtered <- subset(all_cv_results, LianaOriginal %in% c("0", "3"))
# 
# mean_accuracy_pixel_filtered <- mean(all_cv_results_filtered$correct)
# mean_accuracy_pixel_filtered
# 
# accuracy_per_crown_cv <- aggregate(correct ~ poly_id + LianaOriginal + fold,
#                                    data = all_cv_results,
#                                    FUN = mean)
# names(accuracy_per_crown_cv)[4] <- "accuracy"
# 
# accuracy_points_by_fold <- aggregate(correct ~ poly_id + LianaOriginal + fold,
#                                      data = all_cv_results,
#                                      FUN = mean)
# 
# 
# accuracy_points_by_fold_filtered <- subset(accuracy_points_by_fold,LianaOriginal %in% c("0", "3"))
# 
# p3 <- ggplot() +
#   geom_boxplot(
#     data = accuracy_points_by_fold_filtered,
#     aes(x = LianaOriginal, y = correct, color = LianaOriginal),
#     fill = NA,       
#     linewidth = 0.8,
#     show.legend = FALSE
#   ) +
#   geom_jitter(
#     data = accuracy_points_by_fold_filtered,
#     aes(x = LianaOriginal, y = correct, shape = factor(fold), fill = factor(fold)),
#     width = 0.15,
#     alpha = 0.7,
#     size = 2.5,
#     color = "black",
#     stroke = 0.3
#   ) +
#   scale_color_manual(values = c("0" = "pink1", "3" = "palevioletred4")) +
#   scale_fill_manual(values = c("1" = "#B66DFF", "2" = "#FF6DB6", "3" = "#24FF24", "4" = "#6DB6FF", "5"="#DB6D00")) +
#   scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +
#   labs(
#     x = "Liana class in drone photos (DP)",
#     y = "Mean accuracy per crown",
#     title = "Cross validation (StratifiedGroupKFold) on crowns (No filter)",
#     caption = paste0("Overall mean accuracy = ", round(mean_accuracy_pixel_filtered * 100, 1), "%"),
#     shape = "Fold",
#     fill ="Fold"
#   ) +
#   theme(legend.position = "right")+
#   theme_minimal()
# 
# p3
# 
# all_cv_results$predicted_cv_num <- as.numeric(as.character(all_cv_results$predicted))
# 
# predicted_liana_per_crown_cv <- aggregate(
#   predicted_cv_num ~ poly_id + LianaOriginal + fold,
#   data = all_cv_results,
#   FUN = mean
# )
# names(predicted_liana_per_crown_cv)[4] <- "prediction"
# 
# predicted_liana_per_crown_cv$prediction <- predicted_liana_per_crown_cv$prediction * 100
# 
# 
# p4 <- ggplot() +
#   geom_boxplot(
#     data = predicted_liana_per_crown_cv,
#     aes(x = LianaOriginal, y = prediction, color = LianaOriginal),
#     fill = NA,
#     linewidth = 0.8,
#     show.legend = FALSE
#   ) +
#   geom_jitter(
#     data = predicted_liana_per_crown_cv,
#     aes(x = LianaOriginal, y = prediction, shape = factor(fold), fill = factor(fold)),
#     width = 0.15,
#     alpha = 0.7,
#     size = 2.5,
#     color = "black",
#     stroke = 0.3
#   ) +
#   scale_color_manual(values = c("0" = "pink1", "1" = "palevioletred1", "2" = "palevioletred3", "3" = "palevioletred4")) +
#   scale_fill_manual(values = c("1" = "#B66DFF", "2" = "#FF6DB6", "3" = "#24FF24", "4" = "#6DB6FF", "5"="#DB6D00")) +
#   scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +
#   labs(
#     x = "Liana class in drone photos (DP)",
#     y = "Liana predicted per crown (%)",
#     title = "Mean liana prediction on all fold tested (No filter)",
#     caption = paste0("Overall mean accuracy = ", round(mean_accuracy_pixel_filtered * 100, 1), "%"),
#     shape = "Fold", 
#     fill = "Fold"    
#   )+
#   theme(legend.position = "right") +
#   theme_minimal()
# 
# p4
# 
# grid.arrange(p3, p4, ncol = 2)
# 
# # grid.arrange(p1, p2, p3, p4, ncol = 2)
# 
# accuracy_by_fold <- aggregate(
#   correct ~ fold,
#   data = subset(all_cv_results, LianaOriginal %in% c("0", "3")),
#   FUN = mean
# )
# 
# p4_fold_plots <- lapply(1:5, function(current_fold) {
#   fold_data <- subset(predicted_liana_per_crown_cv, fold == current_fold)
#   fold_accuracy <- accuracy_by_fold$correct[accuracy_by_fold$fold == current_fold]
#   
#   ggplot() +
#     geom_boxplot(
#       data = fold_data,
#       aes(x = LianaOriginal, y = prediction, color = LianaOriginal),
#       fill = NA,
#       linewidth = 0.8,
#       show.legend = FALSE
#     ) +
#     geom_jitter(
#       data = fold_data,
#       aes(x = LianaOriginal, y = prediction, fill = factor(fold), shape = factor(fold)),
#       width = 0.15,
#       alpha = 0.7,
#       size = 2.5,
#       color = "black", 
#       stroke = 0.3
#     ) +
#     scale_color_manual(values = c("0" = "pink1", "1" = "palevioletred1", "2" = "palevioletred3", "3" = "palevioletred4")) +  
#     scale_fill_manual(values = c("1" = "#B66DFF", "2" = "#FF6DB6", "3" = "#24FF24", "4" = "#6DB6FF", "5" = "#DB6D00")) +  
#     scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +  
#     labs(
#       x = "Liana class in drone photos (DP)",
#       y = "Liana predicted per crown (%)",
#       title = paste("Testing on fold", current_fold),
#       caption = paste0("Fold mean accuracy = ", round(fold_accuracy * 100, 1), "%")
#     ) +
#     theme(legend.position = "none") 
# })
# 
# grid.arrange(
#   grobs = c(p4_fold_plots, p4),
#   ncol = 2,
#   top = textGrob("Liana prediction by crown and by fold (NIR filter)",
#                  gp = gpar(fontsize = 12, fontface = "bold"))
# )
# 

# Bhattacharyya distance ###################################################

#CROWNS
crowns_df <- as.data.frame(crowns_valid)
crowns_df$ Liana <- as.factor(label_crowns_valid)

crowns_df_pure <- as.data.frame(crowns_valid)
crowns_df_pure$LianaOriginal <- as.factor(label_crowns_valid)
crowns_df_pure$Liana <- as.character(label_crowns_valid)
crowns_df_pure$Liana[crowns_df_pure$Liana %in% c("0", "1")] <- "0"
crowns_df_pure$Liana[crowns_df_pure$Liana %in% c("2", "3")] <- "1"
crowns_df_pure$Liana <- as.factor(crowns_df_pure$Liana)
crowns_df_pure$Type <- factor("crowns", levels = c("crowns", "blob"))
crowns$poly_id <- seq_len(nrow(crowns))
poly_id_r <- rasterize(crowns, area_crowns, field = "poly_id")
poly_id_m <- as.matrix(poly_id_r)[, 1]
poly_id_valid <- poly_id_m[DP_per_pixel$valid_idx]
crowns_df_pure$poly_id <- ifelse(is.na(poly_id_valid), 0, poly_id_valid)

bhatt_cr0_cr3_random <- bhatt_between_random(crowns_df_pure, "LianaOriginal", "0", "3", n_repeats = 100)
bhatt_cr0_cr3_full <- bhatt_between(crowns_df_pure, "LianaOriginal", "0", "3")
bhatt_cr0_within <- bhatt_within_random(crowns_df_pure, "poly_id", "LianaOriginal", "0", n_repeats = 100)
bhatt_cr3_within <- bhatt_within_random(crowns_df_pure, "poly_id", "LianaOriginal", "3", n_repeats = 100)

# BLOB 
blob_df <- as.data.frame(blob_valid)
blob_df$ Liana <- as.factor(label_blob_valid)

blob_df_pure <- as.data.frame(blob_valid)
blob_df_pure$Liana <- as.factor(label_blob_valid)
blob_df_pure$Type <- factor("blob", levels = c("crowns", "blob"))

# blob_df_pure$LianaOriginal <- NA
# blob_df_pure$poly_id <- 0
blob_df_pure$LianaOriginal <- blob_df_pure$Liana
blob$poly_id <- seq_len(nrow(blob))
poly_id_r <- rasterize(blob, area_blob, field = "poly_id")
poly_id_m <- as.matrix(poly_id_r)[, 1]
poly_id_valid <- poly_id_m[MD_per_pixel$valid_idx]
blob_df_pure$poly_id <- ifelse(is.na(poly_id_valid), 0, poly_id_valid)

bhatt_md0_md1_random <- bhatt_between_random(blob_df_pure, "Liana", "0", "1", n_repeats = 100)
bhatt_md0_md1_full <- bhatt_between(blob_df_pure, "Liana", "0", "1")
bhatt_md0_within <- bhatt_within_random(blob_df_pure, "poly_id", "Liana", "0", n_repeats = 100)
bhatt_md1_within <- bhatt_within_random(blob_df_pure, "poly_id", "Liana", "1", n_repeats = 100)

# Blob diagnstic
var_within_0 <- apply(data_blob[blob_df_pure$Liana == "0", ], 2, var)
var_within_1 <- apply(data_blob[blob_df_pure$Liana == "1", ], 2, var)
var_between <- apply(data_blob, 2, function(x) var(c(mean(x[blob_df_pure$Liana == "0"]), mean(x[blob_df_pure$Liana == "1"]))))
 
cat("Mean variance within class 0:", mean(var_within_0), "\n")
cat("Mean variance within class 1:", mean(var_within_1), "\n")
cat("Mean variance between classes:", mean(var_between), "\n")

# Mahalanobis ####
maha_cr0_cr3_random <- maha_between_random(crowns_df_pure, "LianaOriginal", "0", "3", n_repeats = 500)
maha_cr0_cr3_full <- maha_between(crowns_df_pure, "LianaOriginal", "0", "3")
maha_cr0_within <- maha_within_random(crowns_df_pure, "poly_id", "LianaOriginal", "0", n_repeats = 500)
maha_cr3_within <- maha_within_random(crowns_df_pure, "poly_id", "LianaOriginal", "3", n_repeats = 500)

maha_md0_md1_random <- maha_between_random(blob_df_pure, "Liana", "0", "1", n_repeats = 500)
maha_md0_md1_full <- maha_between(blob_df_pure, "Liana", "0", "1")
maha_md0_within <- maha_within_random(blob_df_pure, "poly_id", "Liana", "0", n_repeats = 500)
maha_md1_within <- maha_within_random(blob_df_pure, "poly_id", "Liana", "1", n_repeats = 500)

# Linear Discriminant Analysis (LDA) #########################################

## Within DP dataset

set.seed(123)
crowns_df[,1:ncol(crowns_valid)]<- scale(crowns_df[,1:ncol(crowns_valid)])
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
set.seed(123)
blob_df[,1:ncol(crowns_valid)]<- scale(blob_df[,1:ncol(crowns_valid)]) #ncol(crowns_valid) should be equal to ncol(blob_valid)
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
# crowns_df <- as.data.frame(crowns_valid)
# crowns_df$Liana <- as.character(label_crowns_valid)
# crowns_df <- crowns_df[crowns_df$Liana %in% c('0','3'), ]
# crowns_df$Liana[crowns_df$Liana == '3'] <- '1'
# crowns_df$Liana <- as.factor(crowns_df$Liana)
# type_crowns <- rep(c('crowns'), nrow(crowns_df))
# crowns_df$Type <- as.factor(type_crowns)
# 
# blob_df <- as.data.frame(blob_valid)
# blob_df$ Liana <- as.factor(label_blob_valid)
# type_blob <- rep(c('blob'), nrow(blob_df))
# blob_df$Type <- as.factor(type_blob)
# 
# all_df <- rbind(crowns_df, blob_df)
# 
# set.seed(123)
# all_df[,1:ncol(crowns_valid)]<- scale(all_df[,1:ncol(crowns_valid)])
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

## LDA trained on pure labels to check % of liana predicted ###
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

#### For crowns
# set.seed(123)
# folds_crowns <- stratified_group_kfold(crowns_df_pure, "poly_id", "Liana", k = 5)
# 
# cv_results <- list()
# for (i in 1:length(folds_crowns)) {  
#   test_indices <- folds_crowns[[i]]  
#   train_indices <- setdiff(1:nrow(crowns_df_pure), test_indices) 
#   
#   train_data <- crowns_df_pure[train_indices, ]
#   test_data  <- crowns_df_pure[test_indices, ]
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
#   test_data_filtered <- test_data[test_data$Liana %in% c("0", "3"), ]
#   pred_filtered <- predictions_cv$class[test_data$Liana %in% c("0", "3")]
#   
#   cm_filtered <- confusionMatrix(
#     factor(pred_filtered, levels = c("0", "3")),
#     factor(test_data_filtered$Liana, levels = c("0", "3")),
#     positive = NULL)
#   
#   f1_filtered <- calculate_f1_manual(test_data_filtered$Liana, pred_filtered)
#   
#   cv_results[[i]] <- list(
#     fold = i,
#     accuracy = mean(predictions_cv$class == test_data$Liana),
#     test_data = test_data,
#     f1_filtered = f1_filtered,
#     predictions = predictions_cv,
#     train_indices = train_indices,
#     test_indices = test_indices
#   )
# }
# 
# cv_accuracies <- sapply(cv_results, function(x) x$accuracy)
# cv_f1 <- sapply(cv_results, function(x) x$f1_filtered)
# 
# cat("Accuracies per fold:", paste(round(cv_accuracies, 3), collapse = ", "), "\n")
# cat("Mean accuracy CV:", round(mean(cv_accuracies), 3), "\n")
# cat("F1-scores per fold:", paste(round(cv_f1, 3), collapse = ", "), "\n")
# cat("Mean F1-score:", round(mean(cv_f1, na.rm = TRUE), 3), "\n")
# 
# all_cv_results <- do.call(rbind, lapply(cv_results, function(x) {
#   n_test <- nrow(x$test_data)
#   data.frame(
#     fold = rep(x$fold, n_test),
#     poly_id = x$test_data$poly_id,
#     LianaOriginal = x$test_data$LianaOriginal,
#     correct = x$predictions$class == x$test_data$Liana,
#     predicted = x$predictions$class
#   )
# }))
# 
# all_cv_results_filtered <- subset(all_cv_results, LianaOriginal %in% c("0", "3"))
# 
# mean_accuracy_pixel_filtered <- mean(all_cv_results_filtered$correct)
# mean_accuracy_pixel_filtered
# 
# accuracy_per_crown_cv <- aggregate(correct ~ poly_id + LianaOriginal + fold,
#                                    data = all_cv_results,
#                                    FUN = mean)
# names(accuracy_per_crown_cv)[4] <- "accuracy"
# 
# accuracy_points_by_fold <- aggregate(correct ~ poly_id + LianaOriginal + fold,
#                                      data = all_cv_results,
#                                      FUN = mean)
# 
# 
# accuracy_points_by_fold_filtered <- subset(accuracy_points_by_fold,LianaOriginal %in% c("0", "3"))
# 
# p3 <- ggplot() +
#   geom_boxplot(
#     data = accuracy_points_by_fold_filtered,
#     aes(x = LianaOriginal, y = correct, color = LianaOriginal),
#     fill = NA,       
#     linewidth = 0.8,
#     show.legend = FALSE
#   ) +
#   geom_jitter(
#     data = accuracy_points_by_fold_filtered,
#     aes(x = LianaOriginal, y = correct, shape = factor(fold), fill = factor(fold)),
#     width = 0.15,
#     alpha = 0.7,
#     size = 2.5,
#     color = "black",
#     stroke = 0.3
#   ) +
#   scale_color_manual(values = c("0" = "pink1", "3" = "palevioletred4")) +
#   scale_fill_manual(values = c("1" = "#B66DFF", "2" = "#FF6DB6", "3" = "#24FF24", "4" = "#6DB6FF", "5"="#DB6D00")) +
#   scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +
#   labs(
#     x = "Liana class in drone photos (DP)",
#     y = "Mean accuracy per crown",
#     title = "Cross validation (StratifiedGroupKFold) on crowns (NIR filter)",
#     caption = paste0("Overall mean accuracy (0&3) = ", round(mean_accuracy_pixel_filtered * 100, 1), "%"),
#     shape = "Fold",
#     fill ="Fold"
#   ) +
#   theme(legend.position = "right")+
#   theme_minimal()
# 
# p3
# 
# all_cv_results$predicted_cv_num <- as.numeric(as.character(all_cv_results$predicted))
# 
# predicted_liana_per_crown_cv <- aggregate(
#   predicted_cv_num ~ poly_id + LianaOriginal + fold,
#   data = all_cv_results,
#   FUN = mean
# )
# names(predicted_liana_per_crown_cv)[4] <- "prediction"
# 
# predicted_liana_per_crown_cv$prediction <- predicted_liana_per_crown_cv$prediction * 100
# 
# 
# p4 <- ggplot() +
#   geom_jitter(
#     data = predicted_liana_per_crown_cv,
#     aes(x = LianaOriginal, y = prediction, shape = factor(fold), fill = factor(fold)),
#     width = 0.15,
#     alpha = 0.5,
#     size = 2.5,
#     color = "black",
#     stroke = 0.3
#   ) +
#   geom_boxplot(
#     data = predicted_liana_per_crown_cv,
#     aes(x = LianaOriginal, y = prediction, color = LianaOriginal),
#     fill = NA,
#     linewidth = 0.8,
#     show.legend = FALSE
#   ) +
#   scale_color_manual(values = c("0" = "pink1", "1" = "palevioletred1", "2" = "palevioletred3", "3" = "palevioletred4")) +
#   scale_fill_manual(values = c("1" = "grey", "2" = "grey", "3" = "grey", "4" = "grey", "5"="grey")) +
#   scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +
#   labs(
#     x = "Liana class in drone photos (DP)",
#     y = "Liana predicted per crown (%)",
#     title = "Mean liana prediction on all fold tested (NIR filter)",
#     caption = paste0("Overall mean accuracy (0&3)= ", round(mean_accuracy_pixel_filtered * 100, 1),  "% | ",
#                      "Mean F1-score (0&3) = ", round(mean(cv_f1, na.rm = TRUE), 3)),
#     shape = "Fold", 
#     fill = "Fold"    
#     )+
#   theme(legend.position = "right") +
#   theme_minimal()
# 
# p4
# 
# grid.arrange(p3, p4, ncol = 2)
# 
# # grid.arrange(p1, p2, p3, p4, ncol = 2)
# 
# accuracy_by_fold <- aggregate(
#   correct ~ fold,
#   data = subset(all_cv_results, LianaOriginal %in% c("0", "3")),
#   FUN = mean
# )
# 
# p4_fold_plots <- lapply(1:5, function(current_fold) {
#   fold_data <- subset(predicted_liana_per_crown_cv, fold == current_fold)
#   fold_accuracy <- accuracy_by_fold$correct[accuracy_by_fold$fold == current_fold]
#   fold_f1 <- cv_f1[current_fold] 
# 
# 
#   ggplot() +
#     geom_boxplot(
#       data = fold_data,
#       aes(x = LianaOriginal, y = prediction, color = LianaOriginal),
#       fill = NA,
#       linewidth = 0.8,
#       show.legend = FALSE
#     ) +
#     geom_jitter(
#       data = fold_data,
#       aes(x = LianaOriginal, y = prediction, fill = factor(fold), shape = factor(fold)),
#       width = 0.15,
#       alpha = 0.7,
#       size = 2.5,
#       color = "black", 
#       stroke = 0.3
#     ) +
#     scale_color_manual(values = c("0" = "pink1", "1" = "palevioletred1", "2" = "palevioletred3", "3" = "palevioletred4")) +  
#     scale_fill_manual(values = c("1" = "#B66DFF", "2" = "#FF6DB6", "3" = "#24FF24", "4" = "#6DB6FF", "5" = "#DB6D00")) +  
#     scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +  
#     labs(
#       x = "Liana class in drone photos (DP)",
#       y = "Liana predicted per crown (%)",
#       title = paste("Testing on fold", current_fold),
#       caption = paste0("Fold mean accuracy = ", round(fold_accuracy * 100, 1), "% | ",
#                        "Fold F1-score (0&3) = ", round(cv_f1[current_fold], 3))
#     ) +
#     theme(legend.position = "none") +
#     theme_minimal()
#     })
# 
# grid.arrange(
#   grobs = c(p4_fold_plots, p4),
#   ncol = 2,
#   top = textGrob("Liana prediction by crown and by fold (NIR filter)",
#                  gp = gpar(fontsize = 12, fontface = "bold"))
# )

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
  
  X_test_cv <- test_data[, band_cols_cv]
  X_test_filtered_cv <- X_test_cv[, names(X_test_cv) %in% names(X_train_filtered_cv), drop = FALSE]
  
  
  train_lda_cv <- data.frame(Liana = train_data$Liana, X_train_filtered_cv)
  n_classes <- length(levels(train_lda_cv$Liana))
  model_cv <- lda(Liana ~ ., data = train_lda_cv, tol = 1e-8, prior = rep(1/n_classes, n_classes))
  
  test_lda_cv <- data.frame(Liana = test_data$Liana, X_test_filtered_cv)
  predictions_cv <- predict(model_cv, test_lda_cv)
  
  temp_df <- data.frame(
    poly_id = test_data$poly_id,
    LianaOriginal = test_data$LianaOriginal,
    fold = rep(i, nrow(test_data)),
    pred_class = predictions_cv$class
  )
  
  predicted_liana_per_crown_cv_fold <- aggregate(
    pred_class ~ poly_id + LianaOriginal + fold,
    data = temp_df,
    FUN = function(x) mean(as.numeric(x))
  )
  names(predicted_liana_per_crown_cv_fold) <- c("poly_id", "LianaOriginal", "fold", "prediction")
  predicted_liana_per_crown_cv_fold$prediction <- predicted_liana_per_crown_cv_fold$prediction * 100
  
  poly_predictions <- data.frame(
    poly_id = test_data$poly_id,
    pred_pixel = predictions_cv$class
  )
  
  poly_pred_stats <- aggregate(
    pred_pixel ~ poly_id,
    data = poly_predictions,
    FUN = function(x) mean(x == "3", na.rm = TRUE)
  )
  
  poly_pred_classes <- ifelse(
    poly_pred_stats$pred_pixel <= 0.05, "0",
    ifelse(poly_pred_stats$pred_pixel <= 0.33, "1",
           ifelse(poly_pred_stats$pred_pixel <= 0.75, "2", "3"))
  )
  
  poly_true_labels <- unique(test_data[c("poly_id", "Liana")])
  
  poly_results <- data.frame(
    fold = i,
    poly_id = poly_true_labels$poly_id,
    true = poly_true_labels$Liana,
    pred = poly_pred_classes[match(poly_true_labels$poly_id, poly_pred_stats$poly_id)]
  )
  
  poly_accuracy <- mean(poly_results$pred == poly_results$true, na.rm = TRUE)
  poly_f1 <- calculate_f1_manual(poly_results$true, poly_results$pred)
  
  cv_results[[i]] <- list(
    fold = i,
    predicted_liana_per_crown_cv = predicted_liana_per_crown_cv_fold,
    poly_results = poly_results,
    poly_accuracy = poly_accuracy,
    poly_f1 = poly_f1,
    test_data = test_data,
    predictions = predictions_cv
  )
}

predicted_liana_per_crown_cv <- do.call(rbind, lapply(cv_results, function(x) x$predicted_liana_per_crown_cv))

cv_poly_accuracies <- sapply(cv_results, function(x) x$poly_accuracy)
cv_poly_f1 <- sapply(cv_results, function(x) x$poly_f1)
mean_poly_accuracy <- mean(cv_poly_accuracies)
mean_poly_f1 <- mean(cv_poly_f1, na.rm = TRUE)

accuracy_by_fold <- aggregate(
  pred == true ~ fold,
  data = do.call(rbind, lapply(cv_results, function(x) x$poly_results)),
  FUN = mean
)
names(accuracy_by_fold)[2] <- "correct"

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
  scale_fill_manual(values = c("1" = "grey", "2" = "grey", "3" = "grey", "4" = "grey", "5" = "grey")) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +
  scale_x_discrete(labels = c("0" = "CR0", "1" = "CR1", "2" = "CR2", "3" = "CR3")) +
  labs(
    x = "Liana class in crowns dataset (CR)",
    y = "Liana predicted per crown (%)",
    title = "Mean liana prediction on all fold tested (NIR filter)",
    caption = paste0(
      "Mean accuracy  = ", round(mean_poly_accuracy, 3), " | ",
      "Mean F1-score  = ", round(mean_poly_f1, 3)
    ),
    shape = "Fold", fill = "Fold"
  ) +
  theme(legend.position = "right") +
  theme_minimal()

p4

p4_fold_plots <- lapply(1:5, function(current_fold) {
  fold_data <- subset(predicted_liana_per_crown_cv, fold == current_fold)
  fold_accuracy_pixel <- mean(fold_data$prediction / 100, na.rm = TRUE)  
  fold_accuracy_poly <- cv_poly_accuracies[current_fold]  
  fold_f1_poly <- cv_poly_f1[current_fold]
  
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
    scale_fill_manual(values = c("1" = "#B66DFF", "2" = "#FF6DB6", "3" = "#24FF24", "4" = "#6DB6FF", "5" = "#DB6D00"), name ="Fold") +
    scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25), name="Fold") +
    scale_x_discrete(labels = c("0" = "CR0", "1" = "CR1", "2" = "CR2", "3" = "CR3")) +
    labs(
      x = "Liana class in crowns dataset (CR)",
      y = "Liana predicted per crown (%)",
      title = paste("Testing on fold", current_fold),
      caption = paste0(
        "Fold mean accuracy = ", round(fold_accuracy_poly, 3), " | ",
        "Fold F1-score = ", round(fold_f1_poly, 3)
      )
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

#### For blob
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
  scale_x_discrete(labels = c("0" = "MD0", "1" = "MD1")) +
  
  labs(
    x = "Liana class in manual delineation (MD)",
    y = "Liana predicted per crown (%)",
    title = "Mean liana prediction on all fold tested (NIR filter)",
    caption = paste0("Overall mean accuracy = ", round(mean_accuracy_pixel * 100, 1), "% | ",
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
    scale_x_discrete(labels = c("0" = "MD0", "1" = "MD1")) +
    labs(
      x = "Liana class in manual delineation (MD)",
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

#### CROWNS 
set.seed(123)
folds_crowns <- stratified_group_kfold(crowns_df_pure, "poly_id", "LianaOriginal", k = 5)  
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
  
  X_test_cv <- test_data[, band_cols_cv]
  X_test_filtered_cv <- X_test_cv[, names(X_test_cv) %in% names(X_train_filtered_cv), drop = FALSE]
  
  train_lda_cv <- data.frame(Liana = train_data$LianaOriginal, X_train_filtered_cv)  
  n_classes <- length(levels(train_lda_cv$Liana))
  model_cv <- lda(Liana ~ ., data = train_lda_cv, tol = 1e-8, prior = rep(1/n_classes, n_classes))
  
  test_lda_cv <- data.frame(Liana = test_data$LianaOriginal, X_test_filtered_cv)  
  predictions_cv <- predict(model_cv, test_lda_cv)
  
  poly_predictions <- data.frame(
    poly_id = test_data$poly_id,
    pred_pixel = predictions_cv$class
  )
  
  poly_pred_modes <- aggregate(
    pred_pixel ~ poly_id,
    data = poly_predictions,
    FUN = function(x) {
      as.character(names(which.max(table(x)))) 
    }
  )
  names(poly_pred_modes) <- c("poly_id", "pred_mode")
  
  poly_true_labels <- unique(test_data[c("poly_id", "LianaOriginal")])
  
  poly_results <- data.frame(
    fold = i,
    poly_id = poly_true_labels$poly_id,
    true = poly_true_labels$LianaOriginal,
    pred = poly_pred_modes$pred_mode[match(poly_true_labels$poly_id, poly_pred_modes$poly_id)]
  )
  
  poly_accuracy <- mean(poly_results$pred == poly_results$true, na.rm = TRUE)
  
  pixel_accuracy <- mean(predictions_cv$class == test_data$LianaOriginal, na.rm = TRUE)
  
  cv_results[[i]] <- list(
    fold = i,
    poly_accuracy = poly_accuracy,
    pixel_accuracy = pixel_accuracy,
    poly_results = poly_results,
    test_data = test_data,
    predictions = predictions_cv
  )
}

cv_poly_accuracies <- sapply(cv_results, function(x) x$poly_accuracy)
cv_pixel_accuracies <- sapply(cv_results, function(x) x$pixel_accuracy)

cat("Pixel-level accuracies per fold:", paste(round(cv_pixel_accuracies, 3), collapse = ", "), "\n")
cat("Mean pixel accuracy CV:", round(mean(cv_pixel_accuracies), 3), "\n\n")
cat("Polygon-level accuracies per fold:", paste(round(cv_poly_accuracies, 3), collapse = ", "), "\n")
cat("Mean polygon accuracy CV:", round(mean(cv_poly_accuracies), 3), "\n")
