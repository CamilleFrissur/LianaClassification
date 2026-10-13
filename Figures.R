library(terra)
library(ggplot2)
library(dplyr)
library(gridExtra)

## Set working Directory ####
setwd ("D:/DataCamille/")

# load('Classification_R/LianaClassification/myEnvironment.RData')
#ls()

# Uploading data and cropping ##########################################

## HSI mosaic ####
r <- rast("HSI/mosaic_2023/mosaic/nou_20231025_l1c_mosaic.bsq") #2023
nlyr(r) #number of spectral bands

## Manually delineated (MD) ####

blob <- vect("Lianes/LianaOrNotLiana_55_45.shp")
blob <- project(blob, crs(r))

area_blob <- crop(r, blob, mask=TRUE)

label_blob <- rasterize(blob, area_blob, field = "Type")


## Crowns from drone photos (DP) ####
crowns <- vect("Couronnes/Crowns_Liana_2025_GP_TRUE.shp")
crowns <- project(crowns, crs(r))

area_crowns <- crop(r, crowns, mask=TRUE)

label_crowns <- rasterize(crowns, area_crowns, field = "LianaFullC")

# Maps ################################################################

## Dataset ####
blob_df <- as.data.frame(label_blob, xy = TRUE, na.rm = TRUE) %>%mutate(dataset = "blob (MD)")

crowns_df <- as.data.frame(label_crowns, xy = TRUE, na.rm = TRUE) %>%  mutate(dataset = "crowns (DP)")

plot_df <- bind_rows
plot_df$LianaFullC <- as.factor(plot_df$LianaFullC)
plot_df$Liana <- ifelse(is.na(plot_df$Type), as.character(plot_df$LianaFullC),
                        ifelse(is.na(plot_df$LianaFullC), plot_df$Type,
                               paste0(plot_df$Type, plot_df$LianaFullC)))
ggplot() +
  geom_tile(data = blob_df, aes(x, y, fill = "blob (MD)"), alpha = 1) +
  geom_tile(data = crowns_df, aes(x, y, fill = "crowns (DP)"), alpha = 1) +
  scale_fill_manual(name = "Dataset",
                    values = c("blob (MD)" = "#001312",
                               "crowns (DP)" = "palevioletred1")) +
  theme_minimal() +
  theme(legend.position = "bottom")
## Liana labels ####

blob_df$Type <- as.factor(blob_df$Type)
crowns_df$LianaFullC <- as.factor(crowns_df$LianaFullC)

plot_blob <- ggplot() +
  geom_tile(data = blob_df, aes(x, y, fill = Type), alpha = 1) +
  scale_fill_manual(name = "Type",
                    values = c("0" = "#001312", "1" = "palevioletred1")) +
  theme_minimal() +
  theme(legend.position = "bottom")

plot_crowns <- ggplot() +
  geom_tile(data = crowns_df, aes(x, y, fill = LianaFullC), alpha = 1) +
  scale_fill_manual(name = "Liana",
                    values = c("CR0" = "#001312", "CR1" = "palevioletred1",
                               "CR2" = "blue", "CR3" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

grid.arrange(plot_blob, plot_crowns, ncol = 2)
