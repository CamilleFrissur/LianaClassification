library(terra)
library(ggplot2)
library(dplyr)
library(gridExtra)

## Set working Directory ####
setwd ("D:/DataCamille/")

# Maps ################################################################

## Dataset
blob_df <- as.data.frame(label_blob, xy = TRUE, na.rm = TRUE) %>%mutate(dataset = "blob (MD)")

crowns_df <- as.data.frame(label_crowns, xy = TRUE, na.rm = TRUE) %>%  mutate(dataset = "crowns (DP)")

plot_df <- bind_rows (blob_df, crowns_df)
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

## Liana labels

blob_df$Type <- as.factor(blob_df$Type)
crowns_df$LianaFullC <- as.factor(crowns_df$LianaFullC)

plot_blob <- ggplot() +
  geom_tile(data = blob_df, aes(x, y, fill = Type), alpha = 1) +
  scale_fill_manual(values = c("0" = "#ABABAB", "1" = "#001312")) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  labs(title="Blob dataset", fill = "Liana")

plot_crowns <- ggplot() +
  geom_tile(data = crowns_df, aes(x, y, fill = LianaFullC), alpha = 1) +
  scale_fill_manual(name = "Liana",
                    values = c("CR0" = "pink1", "CR1" = "palevioletred1",
                               "CR2" = "palevioletred3", "CR3" = "palevioletred4")) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  labs(title = "Crowns dataset")

grid.arrange(plot_blob, plot_crowns, ncol = 2)

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
                    values = c("MD0" = "#ababab",
                               "MD1" = "#001312",
                               "CR0"="pink1", "CR1" = "palevioletred1",
                               "CR2" = "palevioletred3", "CR3" = "palevioletred4"))


# PCA ################################################################

## CROWNS

ggplot(pca_plot_crowns, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  # scale_color_manual(values = c("0"= "pink1", "1"="palevioletred1","2"="palevioletred3","3"= "palevioletred4"),
  scale_color_manual(values = c("0"= "pink1", "1"="palevioletred1","2"="palevioletred3","3"= "palevioletred4"),
                     #scale_color_manual(values = c("0"= "#21421e", "1"="#568203","2"="#9acd32","3"= "#c5e384"),
                     
                     labels = c("0" = "CR0", "1" = "CR1", "2"= "CR2", "3"= "CR3"))+
  labs(title= "PCA of the HSI values of the crowns dataset (NIR filter)",
       x = paste0("Dim 1 (", round(var_exp_crowns[1], 1), "%)"),
       y = paste0("Dim 2 (", round(var_exp_crowns[2], 1), "%)"))+
  theme_minimal()

# # Only useful if the illumination filter was used
# pca_plot_crowns$illu <- illu_crowns_valid
# ggplot(pca_plot_crowns, aes(x = x, y = y, color = illu)) +
#   geom_point() +
#   scale_color_viridis_c(name = "Illumination") +
#   labs(title = "PCA of the HSI values of the DP dataset colored by illumination")

ggplot(pca_plot_crowns_poly, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("CR0"= "pink1", "CR1"="palevioletred1","CR2"="palevioletred3","CR3"= "palevioletred4"))+
  labs(title= "PCA of the mean spectrum per polygon of the crowns dataset (NIR filter)",
       x = paste0("Dim 1 (", round(var_exp_crowns_poly[1], 1), "%)"),
       y = paste0("Dim 2 (", round(var_exp_crowns_poly[2], 1), "%)"))+
  theme_minimal()

# BLOB

ggplot(pca_plot_blob, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#ABABAB", "1" = "#001312"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "PCA of the HSI values of the MD dataset (NIR filter)",
       x = paste0("Dim 1 (", round(var_exp_blob[1], 1), "%)"),
       y = paste0("Dim 2 (", round(var_exp_blob[2], 1), "%)"))+
  theme_minimal()

# # Only with illu filter
# pca_plot_blob$illu <- illu_blob_valid
# ggplot(pca_plot_blob, aes(x = x, y = y, color = illu)) +
#   geom_point() +
#   scale_color_viridis_c(name = "Illumination") +
#   labs(title = "PCA of the HSI values of the MD dataset colored by illumination")

ggplot(pca_plot_blob_poly, aes(x = x, y = y, color = Type)) +
  geom_point() +
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#ABABAB", "1" = "#001312"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "PCA of the mean spectrum per polygon of the MD dataset (NIR filter)",
       x = paste0("Dim 1 (", round(var_exp_blob_poly[1], 1), "%)"),
       y = paste0("Dim 2 (", round(var_exp_blob_poly[2], 1), "%)"))+
  theme_minimal()

# tSNE ################################################################

## CROWNS

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

ggplot(tsne_plot_crowns_poly,aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+ #(size = 3)
  scale_color_manual(values = c("CR0"= "pink1", "CR1"="palevioletred1","CR2"="palevioletred3","CR3"= "palevioletred4"))+
  # scale_color_viridis_d(option="turbo")+
  labs(title= "tSNE of the mean spectrum per polygon of the DP dataset (Illu filter)" )+
  theme_minimal()

## BLOB
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

## CROWNS
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

ggplot(umap_plot_crowns_poly, aes(x = x, y = y, color = LianaFullCrown)) +
  geom_point()+
  scale_color_manual(values = c("CR0"= "pink1", "CR1"="palevioletred1","CR2"="palevioletred3","CR3"= "palevioletred4"))+
  labs(title= "UMAP of the mean HSI value per polygon of the DP dataset (Illu filter)" )+
  theme_minimal()

## BLOB

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

ggplot(umap_plot_poly, aes(x = x, y = y, color = Type)) +
  geom_point()+
  scale_color_manual(
    name = "Presence of liana",
    values = c("0" = "#ABABAB", "1" = "#001312"),
    labels = c("0" = "No liana", "1" = "Liana")
  )+
  labs(title= "UMAP of the mean HSI value per polygon of the MD dataset (Illu filter)" )+
  theme_minimal()

# Mean spectrum ################################################################
label_config <- data.frame(
  group = c("DP_CR0", "DP_CR1", "DP_CR2", "DP_CR3", "MD_0", "MD_1"),
  display_name = c("CR0 ( No Liana DP)", "CR1 (DP)", "CR2", "CR3 (Liana DP)", "MD0 (No Liana MD)", "MD1 ( Liana MD)"),
  color = c("pink1", "palevioletred1", "palevioletred3","palevioletred4", "white", "#001312"),
  # color = c("#662700", "#CC5500", "#FF8E33","#FFCA99", "#ABABAB", "#001312"),
  # color = c("#21421e", "#568203", "#9acd32","#c5e384", "#ABABAB", "#001312"),
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
    title = "Spectrum by dataset group (NIR filter)",
    color = "Type",
    linetype = "Type",
    fill = "Type"  
  ) +
  theme_minimal()

# Visible region plot
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


# NIR Region plot
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

#SWIR1

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

## Difference in mean spectrum  ####
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

# Median spectrum ##############################################################
label_config <- data.frame(
  group = c("DP_CR0", "DP_CR1", "DP_CR2", "DP_CR3", "MD_0", "MD_1"),
  display_name = c("CR0 ( No Liana DP)", "CR1 (DP)", "CR2", "CR3 (Liana DP)", "No liana (MD)", "Liana (MD)"),
  color = c("pink1", "palevioletred1", "palevioletred3","palevioletred4", "#ABABAB", "#001312"),
  #color = c("#662700", "#CC5500", "#FF8E33","#FFCA99", "#ABABAB", "#001312"),
  #color = c("#21421e", "#568203", "#9acd32","#c5e384", "#ABABAB", "#001312"),
  
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

plot_data <- data.frame(
  group = character(0),
  distance = numeric(0),
  dataset = character(0),  # "CR" or "MD"
  type = character(0)       # "Within" or "Between"
)

#Within
plot_data <- rbind(
  plot_data,
  data.frame(
    group = rep("CR0", length(bhatt_cr0_within$all_distances)),
    distance = bhatt_cr0_within$all_distances,
    dataset = rep("CR", length(bhatt_cr0_within$all_distances)),
    type = rep("Within_Class", length(bhatt_cr0_within$all_distances))
  ),
  data.frame(
    group = rep("CR3", length(bhatt_cr3_within$all_distances)),
    distance = bhatt_cr3_within$all_distances,
    dataset = rep("CR", length(bhatt_cr3_within$all_distances)),
    type = rep("Within_Class", length(bhatt_cr3_within$all_distances))
  ),
  data.frame(
    group = rep("MD0", length(bhatt_md0_within$all_distances)),
    distance = bhatt_md0_within$all_distances,
    dataset = rep("MD", length(bhatt_md0_within$all_distances)),
    type = rep("Within_Class", length(bhatt_md0_within$all_distances))
  ),
  data.frame(
    group = rep("MD1", length(bhatt_md1_within$all_distances)),
    distance = bhatt_md1_within$all_distances,
    dataset = rep("MD", length(bhatt_md1_within$all_distances)),
    type = rep("Within_Class", length(bhatt_md1_within$all_distances))
  )
)

# Between
plot_data <- rbind(
  plot_data,
  data.frame(
    group = rep("CR0 vs CR3", length(bhatt_cr0_cr3_random$all_distances)),
    distance = bhatt_cr0_cr3_random$all_distances,
    dataset = rep("CR", length(bhatt_cr0_cr3_random$all_distances)),
    type = rep("Between_Class", length(bhatt_cr0_cr3_random$all_distances))
  ),
  data.frame(
    group = rep("MD0 vs MD1", length(bhatt_md0_md1_random$all_distances)),
    distance = bhatt_md0_md1_random$all_distances,
    dataset = rep("MD", length(bhatt_md0_md1_random$all_distances)),
    type = rep("Between_Class", length(bhatt_md0_md1_random$all_distances))
  )
)

ggplot(plot_data, aes(x = group, y = distance, fill = interaction(type, dataset))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(
    values = c(
      "Within_Class.MD" = "#1f77b4",
      "Between_Class.MD" = "#ff7f0e",
      "Between_Class.CR" = "#fdae61", 
      "Within_Class.CR" = "#98c1d9"),
    name = "Distance type"
  ) +
  labs(
    title = "Bhattacharyya distances distribution on Feature Space data",
    x = "Group (Dataset and distance type)",
    y = "Bhattacharyya distance"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
  
#Mahalanobis####################################################################
# Initialisation
plot_maha_data <- data.frame(
  group = character(0),
  distance = numeric(0),
  dataset = character(0),
  type = character(0)
)

# Within
plot_maha_data <- add_to_plot(plot_maha_data, "CR0", maha_cr0_within, "CR", "Within_Class")
plot_maha_data <- add_to_plot(plot_maha_data, "CR3", maha_cr3_within, "CR", "Within_Class")
plot_maha_data <- add_to_plot(plot_maha_data, "MD0", maha_md0_within, "MD", "Within_Class")
plot_maha_data <- add_to_plot(plot_maha_data, "MD1", maha_md1_within, "MD", "Within_Class")

# Between
plot_maha_data <- add_to_plot(plot_maha_data, "CR0 vs CR3", maha_cr0_cr3_random, "CR", "Between_Class")
plot_maha_data <- add_to_plot(plot_maha_data, "MD0 vs MD1", maha_md0_md1_random, "MD", "Between_Class")

ggplot(plot_maha_data, aes(x = group, y = distance, fill = interaction(type, dataset))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(
    values = c(
      "Within_Class.MD" = "#1f77b4",
      "Between_Class.MD" = "#ff7f0e",
      "Between_Class.CR" = "#fdae61", 
      "Within_Class.CR" = "#98c1d9"),
    name = "Distance type"
  ) +
  labs(
    title = "Mahalanobis distances distribution on Feature Space data",
    x = "Group (Dataset and distance type)",
    y = "Mahalanobis distance"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

# LDA ##########################################################################
# Crowns
  
# Blob