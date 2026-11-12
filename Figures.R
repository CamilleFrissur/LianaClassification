library(terra)
library(tidyterra)
library(maptiles)
library(ggplot2)
library(dplyr)
library(gridExtra)

## Set working Directory ####
setwd ("D:/DataCamille/")

# Maps ################################################################

## Dataset
blob_df <- as.data.frame(label_blob, xy = TRUE, na.rm = TRUE) %>%mutate(dataset = "MD")

crowns_df <- as.data.frame(label_crowns, xy = TRUE, na.rm = TRUE) %>%  mutate(dataset = "CR")

plot_df <- bind_rows (blob_df, crowns_df)
plot_df$LianaFullC <- as.factor(plot_df$LianaFullC)
plot_df$Class <- ifelse(is.na(plot_df$Type), as.character(plot_df$LianaFullC),
                        ifelse(is.na(plot_df$LianaFullC), plot_df$Type,
                               paste0(plot_df$Type, plot_df$LianaFullC)))

## Liana labels
blob_df$Type <- as.factor(blob_df$Type)
crowns_df$LianaFullC <- as.factor(crowns_df$LianaFullC)

# Retrieve topo
bbox_sf <- st_as_sfc(st_bbox(c(xmin = min(plot_df$x), ymin = min(plot_df$y),
                               xmax = max(plot_df$x), ymax = max(plot_df$y)),
                             crs = st_crs(crs(label_blob)))) |>
  st_as_sf()

tiles <- tiles <- get_tiles(bbox_sf, provider = "OpenTopoMap", crop = FALSE, zoom = 14)

tiles_proj <- project(tiles, crs(label_blob))

# Plots
plot_all <- ggplot() +
  geom_spatraster_rgb(data = tiles_proj, maxcell = Inf) +
  geom_tile(data = blob_df, aes(x, y, fill = "MD"), alpha = 1) +
  geom_tile(data = crowns_df, aes(x, y, fill = "CR"), alpha = 1) +
  scale_fill_manual(name = "Dataset",
                    values = c("MD" = "#1965FF",
                               "CR" = "#FF9932")) +
  theme_minimal() +
  theme(legend.position = "bottom")

plot_blob <- ggplot() +
  geom_tile(data = blob_df, aes(x, y, fill = Type), alpha = 1) +
  scale_fill_manual(values = c("0" = "#99EDFF", "1" = "#1965FF")) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  labs(title="b. Manually delineated (MD) dataset", fill = "Class")

plot_crowns <- ggplot() +
  geom_tile(data = crowns_df, aes(x, y, fill = LianaFullC), alpha = 1) +
  scale_fill_manual(name = "Class",
                    values = c("CR0" = "#FFEE99", "CR1" = "#FFCC65",
                               "CR2" = "#FF9932", "CR3" = "#FF6619")) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  labs(title = "a. Crowns (CR) dataset")

plot_all
grid.arrange( plot_crowns, plot_blob, ncol = 2)

# Boxplot illu ######################################################
ggplot(df_all, aes(x = label,
                   y = illumination,
                   fill = label)) +
  geom_boxplot() +
  labs(
    x = "Liana label",
    y = "Illumination",
    fill = "Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    legend.position = "top"                           
  )+
  scale_fill_manual(name = "Type",
                    values = c("MD0" = "#99EDFF",
                               "MD1" = "#1965FF",
                               "CR0"="#FFEE99", "CR1" = "#FFCC65",
                               "CR2" = "#FF9932", "CR3" = "#FF6619"))


# PCA ################################################################

## CROWNS
pca_plot_crowns_half <- pca_plot_crowns %>%
  dplyr::filter(LianaFullCrown %in% c(0, 3)) %>%
  dplyr::mutate(LianaFullCrown = droplevels(factor(LianaFullCrown)))
pca_plot_crowns_half <- pca_plot_crowns_half %>%
  sample_frac(0.5)

ggplot(pca_plot_crowns_half, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  # scale_color_manual(values = c("0"= "#FFEE99", "1"="#FFCC65","2"="#FF9932","3"= "#FF6619"),
  scale_color_manual(values = c("0"= "#FFCC65", "1"="#FFCC65","2"="#FF9932","3"= "#FF6619"),
                     #scale_color_manual(values = c("0"= "#21421e", "1"="#568203","2"="#9acd32","3"= "#c5e384"),
                     
                     labels = c("0" = "CR0", "1" = "CR1", "2"= "CR2", "3"= "CR3"))+
  labs(title= "
       PCA on CR feature space data per pixel (NIR filter)",
       x = paste0("Dim 1 (", round(var_exp_crowns[1], 1), "%)"),
       y = paste0("Dim 2 (", round(var_exp_crowns[2], 1), "%)"))+
  theme_minimal()

# # Only useful if the illumination filter was used
# pca_plot_crowns$illu <- illu_crowns_valid
# ggplot(pca_plot_crowns, aes(x = x, y = y, color = illu)) +
#   geom_point() +
#   scale_color_viridis_c(name = "Illumination") +
#   labs(title = "pca of the HSI values of the DP dataset colored by illumination")

ggplot(pca_plot_crowns_poly, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("CR0"= "#FFEE99", "CR1"="#FFCC65","CR2"="#FF9932","CR3"= "#FF6619"))+
  labs(title= "PCA on CR feature space data per crown (NIR filter)",
       x = paste0("Dim 1 (", round(var_exp_crowns_poly[1], 1), "%)"),
       y = paste0("Dim 2 (", round(var_exp_crowns_poly[2], 1), "%)"))+
  theme_minimal()

# BLOB
pca_plot_blob_half <- pca_plot_blob %>%
  sample_frac(0.33)

ggplot(pca_plot_blob_half, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#99EDFF", "1" = "#1965FF"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "PCA on MD feature space data per pixel (NIR filter)",
       x = paste0("Dim 1 (", round(var_exp_blob[1], 1), "%)"),
       y = paste0("Dim 2 (", round(var_exp_blob[2], 1), "%)"))+
  theme_minimal()+ ylim (-15,15)

# # Only with illu filter
# umap_plot_blob$illu <- illu_blob_valid
# ggplot(umap_plot_blob, aes(x = x, y = y, color = illu)) +
#   geom_point() +
#   scale_color_viridis_c(name = "Illumination") +
#   labs(title = "umap of the HSI values of the MD dataset colored by illumination")

ggplot(pca_plot_blob_poly, aes(x = x, y = y, color = Type)) +
  geom_point() +
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#99EDFF", "1" = "#1965FF"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "PCA on MD feature space data per patch (NIR filter)",
       x = paste0("Dim 1 (", round(var_exp_blob_poly[1], 1), "%)"),
       y = paste0("Dim 2 (", round(var_exp_blob_poly[2], 1), "%)"))+
  theme_minimal() + ylim (-5,10)

# tSNE ################################################################

## CROWNS

ggplot(tsne_plot_crowns_half, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("0"= "#FFEE99", "1"="#FFCC65","2"="#FF9932","3"= "#FF6619"),
                     labels = c("0" = "CR0", "1" = "CR1", "2"= "CR2", "3"= "CR3"))+
  labs(title= "tSNE of the HSI values of the drone photos dataset (Illu filter)" )+
  theme_minimal()

# Only useful if the illumination filter was used
# tsne_plot_crowns$illu <- illu_crowns_valid
# ggplot(tsne_plot_crowns, aes(x = x, y = y, color = illu)) +
#   geom_point() +
#   scale_color_viridis_c(name = "Illumination") +
#   labs(title = "tSNE of the HSI values of the DP dataset colored by illumination")

ggplot(tsne_plot_crowns_poly,aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+ #(size = 3)
  scale_color_manual(values = c("CR0"= "#FFEE99", "CR1"="#FFCC65","CR2"="#FF9932","CR3"= "#FF6619"))+
  # scale_color_viridis_d(option="turbo")+
  labs(title= "tSNE of the mean spectrum per polygon of the DP dataset (Illu filter)" )+
  theme_minimal()

## BLOB
ggplot(tsne_plot_blob, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#99EDFF", "1" = "#1965FF"),
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

ggplot(tsne_plot_blob_poly, aes(x = x, y = y, color = Type)) +
  geom_point() + #(size = 3)
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#99EDFF", "1" = "#1965FF"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "tSNE of the mean spectrum per polygon of the MD dataset (Illu filter)" )+
  theme_minimal()

## UMAP #########################################################################

## CROWNS
umap_plot_crowns_half <- umap_plot_crowns %>%
  dplyr::filter(LianaFullCrown %in% c(0, 3)) %>%
  dplyr::mutate(LianaFullCrown = droplevels(factor(LianaFullCrown)))
umap_plot_crowns_half <- umap_plot_crowns_half %>%
  sample_frac(0.5)

ggplot(umap_plot_crowns_half, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("0"= "#FFCC65", "1"="#FFCC65","2"="#FF9932","3"= "#FF6619"),
                     labels = c("0" = "CR0", "1" = "CR1", "2"= "CR2", "3"= "CR3"))+
  # scale_color_viridis_d(option="turbo")+
  labs(title= "UMAP on CR feature space data per pixel (NIR filter)" )+
  theme_minimal()

# Only with illu filter
# umap_plot_crowns$illu <- illu_crowns_valid
# ggplot(umap_plot_crowns, aes(x = x, y = y, color = illu)) +
# geom_point() +
# scale_color_viridis_c(name = "Illumination") +
# labs(title = "UMAP of the HSI values of the DP dataset colored by illumination")

ggplot(umap_plot_crowns_poly, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("CR0"= "#FFEE99", "CR1"="#FFCC65","CR2"="#FF9932","CR3"= "#FF6619"))+
  labs(title= "UMAP on CR feature space data per crown (NIR filter)" )+
  theme_minimal()

## BLOB
umap_plot_blob_half <- umap_plot_blob %>%
  sample_frac(0.33)


ggplot(umap_plot_blob_half, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#99EDFF", "1" = "#1965FF"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "UMAP on MD feature space data per pixel (NIR filter)" )+
  theme_minimal()

# Only with illu filter
# umap_plot_blob$illu <- illu_blob_valid
# ggplot(umap_plot_blob, aes(x = x, y = y, color = illu)) +
#   geom_point() +
#   scale_color_viridis_c(name = "Illumination") +
#   labs(title = "UMAP of the HSI values of the MD dataset colored by illumination")

ggplot(umap_plot_poly, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#99EDFF", "1" = "#1965FF"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "UMAP on MD feature space data per patch (NIR filter)" )+
  theme_minimal()

# Mean spectrum ################################################################
label_config <- data.frame(
  group = c("DP_CR0", "DP_CR1", "DP_CR2", "DP_CR3", "MD_0", "MD_1"),
  display_name = c("CR0 (No Liana)", "CR1 (DP)", "CR2", "CR3 (Liana)", "MD0 (No Liana)", "MD1 ( Liana)"),
  color = c("#FFEE99", "#FFCC65", "#FF9932","#FF6619", "#99EDFF", "#1965FF"),
  # color = c("#662700", "#CC5500", "#FF8E33","#FFCA99", "#99EDFF", "#1965FF"),
  # color = c("#21421e", "#568203", "#9acd32","#c5e384", "#99EDFF", "#1965FF"),
  linetype = c ("solid","solid","solid","solid", "solid", "solid")
  #linetype = c ("solid","solid","solid","solid", "solid", "solid")
  )

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
    title = "Mean spectrum by dataset class with discrepancy (NIR filter)",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()

# Visible region plot
visible <- ggplot(selected_smooth[selected_smooth$wavelength < 690,], aes(x = wavelength)) +
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
    title = "b. Visible region",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()


# NIR Region plot
NIR<- ggplot(selected_smooth[selected_smooth$wavelength >700 & selected_smooth$wavelength < 1340,], aes(x = wavelength)) +
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
    title = "Mean spectrum with discrepancy in the NIR region",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()

#SWIR1

SWIR1 <- ggplot(selected_smooth[selected_smooth$wavelength >1420 & selected_smooth$wavelength < 1790,], aes(x = wavelength)) +
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
    title = "d. SWIR 1 region",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()

grid.arrange (visible + theme(legend.position = "none"), NIR + theme(legend.position = "none"), SWIR1 + theme(legend.position = "none"), ncol=3)

## Difference in mean spectrum  ####
diff_label_config <- data.frame(
  group = c("CR3 - CR0", "MD1 - MD0"),
  display_name = c("CR3 - CR0", "MD1 - MD0"),
  color = c("#FF6619", "#1965FF"), 
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
    title = "Spectrum difference by dataset (NIR filter)",
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

# Median spectrum ##############################################################
label_config <- data.frame(
  group = c("DP_CR0", "DP_CR1", "DP_CR2", "DP_CR3", "MD_0", "MD_1"),
  display_name = c("CR0 ( No Liana DP)", "CR1 (DP)", "CR2", "CR3 (Liana DP)", "No liana (MD)", "Liana (MD)"),
  color = c("#FFEE99", "#FFCC65", "#FF9932","#FF6619", "#99EDFF", "#1965FF"),
  #color = c("#662700", "#CC5500", "#FF8E33","#FFCA99", "#99EDFF", "#1965FF"),
  #color = c("#21421e", "#568203", "#9acd32","#c5e384", "#99EDFF", "#1965FF"),
  
  linetype = c ("solid","solid","solid","solid", "solid", "solid"))
  #linetype = c ("solid","solid","solid","solid", "solid", "solid"))

  ggplot(med_selected_smooth, aes(x = wavelength)) +
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
  
  ggplot(med_selected_smooth, aes(x = wavelength)) +
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
  
  ggplot(med_selected_smooth, aes(x = wavelength)) +
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
  
  ggplot(med_selected_smooth[med_selected_smooth$wavelength < 690,], aes(x = wavelength)) +
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
  
  ggplot(med_selected_smooth[med_selected_smooth$wavelength >700 & med_selected_smooth$wavelength < 1340,], aes(x = wavelength)) +
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
  
  ggplot(med_selected_smooth[med_selected_smooth$wavelength >1420 & med_selected_smooth$wavelength < 1790,], aes(x = wavelength)) +
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
  
# Bhattacharyya ################################################################

# Initialisation
plot_bhatt_data <- data.frame(
  group = character(0),
  distance = numeric(0),
  dataset = character(0),
  type = character(0)
)

# Within
plot_bhatt_data <- add_to_plot(plot_bhatt_data, "CR0", bhatt_cr0_within, "CR", "Within_Class")
plot_bhatt_data <- add_to_plot(plot_bhatt_data, "CR3", bhatt_cr3_within, "CR", "Within_Class")
plot_bhatt_data <- add_to_plot(plot_bhatt_data, "MD0", bhatt_md0_within, "MD", "Within_Class")
plot_bhatt_data <- add_to_plot(plot_bhatt_data, "MD1", bhatt_md1_within, "MD", "Within_Class")

# Between
plot_bhatt_data <- add_to_plot(plot_bhatt_data, "CR0 vs CR3", bhatt_cr0_cr3_random, "CR", "Between_Class")
plot_bhatt_data <- add_to_plot(plot_bhatt_data, "MD0 vs MD1", bhatt_md0_md1_random, "MD", "Between_Class")

ggplot(plot_bhatt_data, aes(x = group, y = distance, fill = interaction(type, dataset))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(
    values = c(
      "Within_Class.MD" = "#1965FF",
      "Between_Class.MD" = "#FF6619",
      "Between_Class.CR" = "#FFCC65", 
      "Within_Class.CR" = "#99EDFF"),
    name = "Distance type"
  ) +
  labs(
    title = "Bhattacharyya distances distribution on HSI raw data",
    x = "Group (Dataset and distance type)",
    y = "Bhattacharyya distance"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )+
  theme_minimal()


#Mahalanobis####################################################################
# # Initialisation
# plot_maha_data <- data.frame(
#   group = character(0),
#   distance = numeric(0),
#   dataset = character(0),
#   type = character(0)
# )
# 
# # Within
# plot_maha_data <- add_to_plot(plot_maha_data, "CR0", maha_cr0_within, "CR", "Within_Class")
# plot_maha_data <- add_to_plot(plot_maha_data, "CR3", maha_cr3_within, "CR", "Within_Class")
# plot_maha_data <- add_to_plot(plot_maha_data, "MD0", maha_md0_within, "MD", "Within_Class")
# plot_maha_data <- add_to_plot(plot_maha_data, "MD1", maha_md1_within, "MD", "Within_Class")
# 
# # Between
# plot_maha_data <- add_to_plot(plot_maha_data, "CR0 vs CR3", maha_cr0_cr3_random, "CR", "Between_Class")
# plot_maha_data <- add_to_plot(plot_maha_data, "MD0 vs MD1", maha_md0_md1_random, "MD", "Between_Class")

# add_to_plot <- function(plot_df, group, res_component, dataset, type) {
#   if (is.null(res_component)) return(plot_df)   # comparaison échouée -> ignorée
#   rbind(plot_df, data.frame(
#     group    = group,
#     distance = res_component$all_distances,
#     dataset  = dataset,
#     type     = type
#   ))
# }

plot_maha_data <- data.frame(
  group    = character(0),
  distance = numeric(0),
  dataset  = character(0),
  type     = character(0)
)

# Within
plot_maha_data <- add_to_plot(plot_maha_data, "CR0", res_cr03$within_class1, "CR", "Within_Class")
plot_maha_data <- add_to_plot(plot_maha_data, "CR3", res_cr03$within_class2, "CR", "Within_Class")
plot_maha_data <- add_to_plot(plot_maha_data, "MD0", res_md01$within_class1, "MD", "Within_Class")
plot_maha_data <- add_to_plot(plot_maha_data, "MD1", res_md01$within_class2, "MD", "Within_Class")

# Between
plot_maha_data <- add_to_plot(plot_maha_data, "CR0 vs CR3", res_cr03$between, "CR", "Between_Class")
plot_maha_data <- add_to_plot(plot_maha_data, "MD0 vs MD1", res_md01$between, "MD", "Between_Class")

plot_maha_data$group <- factor(plot_maha_data$group,
                               levels = c("CR0", "CR0 vs CR3","CR3", 
                                          "MD0",  "MD0 vs MD1","MD1"))

ggplot(plot_maha_data, aes(x = group, y = distance, fill = interaction(type, dataset))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(
    values = c(
      "Within_Class.MD" = "#1965FF",
      "Between_Class.MD" = "#99EDFF",
      "Between_Class.CR" = "#FFCC65",
      "Within_Class.CR" = "#FF6619"),
    name = "Distance type"
  ) +
  labs(
    title = "Mahalanobis distances distribution on feature space data",
    x = "Group (Dataset and distance type)",
    y = "Mahalanobis distance"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )+
  theme_minimal()

# LDA ##########################################################################
# Crowns
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
  scale_color_manual(values = c("0" = "#FFEE99", "1" = "#FFCC65", "2" = "#FF9932", "3" = "#FF6619")) +
  scale_fill_manual(values = c("1" = "grey", "2" = "grey", "3" = "grey", "4" = "grey", "5" = "grey")) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +
  scale_x_discrete(labels = c("0" = "CR0", "1" = "CR1", "2" = "CR2", "3" = "CR3")) +
  labs(
    x = "Liana class in crowns dataset (CR)",
    y = "Liana predicted per crown (%)",
    title = "a. Mean liana prediction on HSI raw data (CR)",
    caption = paste0(
      "Mean accuracy (crown) = ", round(mean_poly_accuracy, 3), " | ",
      "Mean F1-score (crown) = ", round(mean_poly_f1, 3)
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
    scale_color_manual(values = c("0" = "#FFEE99", "1" = "#FFCC65", "2" = "#FF9932", "3" = "#FF6619")) +
    scale_fill_manual(values = c("1" = "#B66DFF", "2" = "#FF6DB6", "3" = "#24FF24", "4" = "#6DB6FF", "5" = "#DB6D00"), name ="Fold") +
    scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25), name="Fold") +
    scale_x_discrete(labels = c("0" = "CR0", "1" = "CR1", "2" = "CR2", "3" = "CR3")) +
    labs(
      x = "Liana class in crowns dataset (CR)",
      y = "Liana predicted per crown (%)",
      title = paste("Testing on fold", current_fold),
      caption = paste0(
        "Fold mean accuracy (crown) = ", round(fold_accuracy_poly, 3), " | ",
        "Fold F1-score (crown) = ", round(fold_f1_poly, 3)
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

# Blob
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
  scale_color_manual(values = c("0" = "#99EDFF", "1" = "#1965FF")) +
  scale_fill_manual(values = c("1" = "grey", "2" = "grey", "3" = "grey", "4" = "grey", "5"="grey")) +
  scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +
  scale_x_discrete(labels = c("0" = "MD0", "1" = "MD1")) +
  
  labs(
    x = "Liana class in manual delineation (MD)",
    y = "Liana predicted per patch (%)",
    title = "b. Mean liana prediction on HSI raw data (MD)",
    caption = paste0("Mean accuracy (patch) = ", round(mean_accuracy_poly), " | ",
                     "Mean F1-score (patch) = ", round(f1_macro_poly, 3)),
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
  fold_accuracy_poly <- accuracy_poly_by_fold$correct_majority[
    accuracy_poly_by_fold$fold == current_fold]
  fold_f1_poly <- f1_poly_by_fold[current_fold]
  
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
    scale_color_manual(values = c("0" = "#99EDFF", "1" = "#1965FF")) +
    scale_fill_manual(values = c("1" = "#B66DFF", "2" = "#FF6DB6", "3" = "#24FF24", "4" = "#6DB6FF", "5" = "#DB6D00")) +
    scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24, "5" = 25)) +
    scale_x_discrete(labels = c("0" = "MD0", "1" = "MD1")) +
    labs(
      x = "Liana class in manual delineation (MD)",
      y = "Liana predicted per patch (%)",
      title = paste("Testing on fold", current_fold),
      caption = paste0("Fold accuracy (atch) = ", round(fold_accuracy_poly) "| ",
                       "Fold F1-score (atch) = ", round(fold_f1_poly, 3))
    ) +
    theme(legend.position = "none") +
    theme_minimal()
})

grid.arrange(
  grobs = c(p4_fold_plots_blob, p4_blob),
  ncol = 2,
  top = textGrob("Liana prediction by patch and by fold (NIR filter)",
                 gp = gpar(fontsize = 12, fontface = "bold"))
)

grid.arrange(p4, p4_blob, ncol = 2)
