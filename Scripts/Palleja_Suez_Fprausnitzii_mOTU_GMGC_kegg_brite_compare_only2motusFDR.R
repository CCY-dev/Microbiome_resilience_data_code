rm(list = ls())
setwd("~/Documents/Lab/Multiple_longitudinal/Microbe_resilience/")
library(tidyverse)
library(readxl)

###### Read in the sig. GMGC associated with 108 and 110
finalSig <- read.table("20231025_Palleja_Suez_Fprausnitzii_108_110_mOTU_GMGC_either_cor_sig.txt", sep = "\t",
                        header = T)

# Read the extracted GMGC file of 108 and 110
gmgc_108_110 <- read.table("20231025_Palleja_Suez_Fprausnitzii_06108_06110_GMGC_either_sig_GMGC_extract.txt", 
                           sep = "\t", header = T)

# Bind the finalSig and gmgc_108_110 together
colnames(gmgc_108_110)[1] <- "GMGC"
gmgc_108_110_final <- left_join(finalSig, gmgc_108_110, by = "GMGC")

# Read the Kegg brite table
brite_info <- read_xlsx("/Users/Jessica/Documents/Lab/Sofia_Kegg_perl/20231018_Kegg_brite_myedit.xlsx", 
                        sheet = 2)

###### How 108 and 110 differ in their brite pathway profile?
gmgc_108_brite <- gmgc_108_110_final %>% 
  dplyr::select(c(1, 2, 3, 4, 17)) %>% 
  filter(mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06108]" &
         Spearman_p_fdr < 0.05 & 
         Spearman_rho > 0.95)

gmgc_108_brite_long <- gmgc_108_brite %>% 
  tidyr::separate_rows(BRITE, sep = ",") %>% #Unfold the cells with multiple values
  filter(BRITE != "") # Remove those blank rows

gmgc_110_brite <- gmgc_108_110_final %>% 
  dplyr::select(c(1, 2, 3, 4, 17)) %>% 
  filter(mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06110]" &
           Spearman_p_fdr < 0.05 & 
           Spearman_rho > 0.95)

gmgc_110_brite_long <- gmgc_110_brite %>% 
  tidyr::separate_rows(BRITE, sep = ",") %>% #Unfold the cells with multiple values
  filter(BRITE != "") # Remove those blank rows

all_brite_long <- rbind(gmgc_108_brite_long, gmgc_110_brite_long)
colnames(all_brite_long)[5] <- "Brite" 

all_brite_long_final <- left_join(all_brite_long, brite_info, by = "Brite") %>% 
  filter(Level1 != "") # Remove those blank rows without level info

# Plot Venn diagram of the brites of 06108 and 06110
library(ggvenn)
#pdf(file = "20231030_Fp_mOTU_significant_cor_GMGC_Kegg_brite_venn.pdf", width = 16, height = 12)

# Individaul Kegg Brite
for_venn_brite <- list(Fp_mOTU_06108 = all_brite_long_final$Brite[all_brite_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06108]"], 
                        Fp_mOTU_06110 = all_brite_long_final$Brite[all_brite_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06110]"])

ggvenn(for_venn_brite, stroke_size = 0, set_name_size = 5, stroke_color = "white",
       text_size = 1.8, auto_scale = T, show_percentage = F, show_elements = T,
       label_sep = "\n",
       fill_color = c("#FAC866", "#99C1FA")) +
  theme(plot.margin = unit(c(0,0,0, 0), "cm"), 
        axis.title = element_text(size = 6),
        title = element_text(size = 14)) +
  labs(x = "", y = "", title = "Individual Kegg Brite")

# Individaul Keg  Brite with description
for_venn_brite_description <- list(Fp_mOTU_06108 = paste0(all_brite_long_final$Brite[all_brite_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06108]"], ":", all_brite_long_final$Brite_description[all_brite_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06108]"]), 
                       Fp_mOTU_06110 = paste0(all_brite_long_final$Brite[all_brite_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06110]"], ":", all_brite_long_final$Brite_description[all_brite_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06110]"]))

ggvenn(for_venn_brite_description, stroke_size = 0, set_name_size = 5, stroke_color = "white",
       text_size = 1.8, auto_scale = T, show_percentage = F, show_elements = T,
       label_sep = "\n",
       fill_color = c("#FAC866", "#99C1FA")) +
  theme(plot.margin = unit(c(0,0,0, 0), "cm"), 
        axis.title = element_text(size = 6),
        title = element_text(size = 14)) +
  labs(x = "", y = "", title = "Individual Kegg Brite with description")

# Level 1
for_venn_level1 <- list(Fp_mOTU_06108 = all_brite_long_final$Level1[all_brite_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06108]"], 
                        Fp_mOTU_06110 = all_brite_long_final$Level1[all_brite_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06110]"])

ggvenn(for_venn_level1, stroke_size = 0, set_name_size = 5, stroke_color = "white",
       text_size = 5, auto_scale = T, show_percentage = F, 
       fill_color = c("#FAC866", "#99C1FA")) +
  theme(plot.margin = unit(c(0,0,0, 0), "cm"), 
        axis.title = element_text(size = 6),
        title = element_text(size = 14)) +
  labs(x = "", y = "", title = "Kegg Brite level 1")

# Level 2
for_venn_level2 <- list(Fp_mOTU_06108 = all_brite_long_final$Level2[all_brite_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06108]"], 
                        Fp_mOTU_06110 = all_brite_long_final$Level2[all_brite_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06110]"])

ggvenn(for_venn_level2, stroke_size = 0, set_name_size = 5, stroke_color = "white",
       text_size = 4, auto_scale = T, show_percentage = F, show_elements = T,
       label_sep = "\n",
       fill_color = c("#FAC866", "#99C1FA")) +
  theme(plot.margin = unit(c(0,0,0, 0), "cm"), 
        axis.title = element_text(size = 6),
        title = element_text(size = 14)) +
  labs(x = "", y = "",  title = "Kegg Brite level 2")

#dev.off()

#### Summarize brite up to Level1 and Level2
all_brite_summary_level1 <- all_brite_long_final %>%
  group_by(mOTU, Level1) %>% 
  summarise(Mean_rho = mean(Spearman_rho))

all_brite_summary_level2 <- all_brite_long_final %>%
  group_by(mOTU, Level2) %>% 
  summarise(Mean_rho = mean(Spearman_rho))

# Plot heatmap
#pdf(file = "20231030_Fp_mOTU_significant_cor_GMGC_brite_heatmap.pdf", width = 12, height = 8)

# Level 1
all_brite_summary_level1_wide <- pivot_wider(all_brite_summary_level1, names_from = Level1, values_from = Mean_rho) %>% 
  column_to_rownames("mOTU") %>% 
  t() %>% 
  as.data.frame()
all_brite_summary_level1_wide[is.na(all_brite_summary_level1_wide)] <- 0

pd <- as.data.frame(scale((all_brite_summary_level1_wide)))
ord <- hclust(dist(pd, method = "binary"), method = "median")$order
ord

all_brite_summary_level1_plotting <- all_brite_summary_level1
#NA_df <- data.frame(mOTU = c(rep("Faecalibacterium prausnitzii [ref_mOTU_v25_06108]", times = sum(all_brite_summary_level1_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06108]` == 0)),
#                             rep("Faecalibacterium prausnitzii [ref_mOTU_v25_06110]", times = sum(all_brite_summary_level1_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06110]` == 0))),
#                    Level1 = c(rownames(all_brite_summary_level1_wide)[all_brite_summary_level1_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06108]` == 0],
#                               rownames(all_brite_summary_level1_wide)[all_brite_summary_level1_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06110]` == 0]),
#                    Mean_rho = NA)

#all_brite_summary_level1_plotting <- rbind(all_brite_summary_level1_plotting, NA_df)
all_brite_summary_level1_plotting$Level1 <- factor(all_brite_summary_level1_plotting$Level1 , 
                                                    levels = rownames(pd)[ord])

ggplot(all_brite_summary_level1_plotting, aes(mOTU, Level1, fill = Mean_rho)) +
  geom_tile(color = "white", size = 0.5) +
  theme_light() +
  theme(plot.title = element_text(size = 18), 
        plot.subtitle = element_text(size = 10), 
        legend.position="right",
        #plot.caption = element_text(hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 14),
        panel.border = element_blank()) +                 
  guides(fill = guide_colourbar(title = "Spearman's rho",
                                raster = FALSE, nbin = 20,
                                title.position = "top",                  
                                barwidth = 0.5,                               
                                title.hjust = 0.5)) +
  scale_fill_gradient2(high = "red3", mid = "cornsilk1", midpoint = 0.9,
                       limits=c(0.9, 1), breaks=seq(0.9, 1, by = 0.05), 
                       na.value = "grey") + 
  scale_x_discrete(labels = c("F. prausnitzii [ref_mOTU_v25_06108]",
                              "F. prausnitzii [ref_mOTU_v25_06110]")) +
  ggpubr::rotate_x_text(angle = 0) +
  labs(x = "", y = "", title = "KEGG Brite level 1")

# Level2
# Clustered automatically by hclust
# Ref: https://stackoverflow.com/questions/25528059/cluster-data-in-heat-map-in-r-ggplot
all_brite_summary_level2_wide <- pivot_wider(all_brite_summary_level2, names_from = Level2, values_from = Mean_rho) %>% 
  column_to_rownames("mOTU") %>% 
  t() %>% 
  as.data.frame()
all_brite_summary_level2_wide[is.na(all_brite_summary_level2_wide)] <- 0

pd <- as.data.frame(scale((all_brite_summary_level2_wide)))
ord <- hclust(dist(pd, method = "binary"), method = "median")$order
ord

all_brite_summary_level2_plotting <- all_brite_summary_level2

#NA_df <- data.frame(mOTU = c(rep("Faecalibacterium prausnitzii [ref_mOTU_v25_06108]", times = sum(all_brite_summary_level2_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06108]` == 0)),
#                             rep("Faecalibacterium prausnitzii [ref_mOTU_v25_06110]", times = sum(all_brite_summary_level2_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06110]` == 0))),
#                    Level2 = c(rownames(all_brite_summary_level2_wide)[all_brite_summary_level2_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06108]` == 0],
#                               rownames(all_brite_summary_level2_wide)[all_brite_summary_level2_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06110]` == 0]),
#                    Mean_rho = NA)

#all_brite_summary_level2_plotting <- rbind(all_brite_summary_level2_plotting, NA_df)
all_brite_summary_level2_plotting$Level2 <- factor(all_brite_summary_level2_plotting$Level2 , levels = rownames(pd)[ord])

ggplot(all_brite_summary_level2_plotting, aes(mOTU, Level2, fill = Mean_rho)) +
  geom_tile(color = "white", size = 0.5) +
  theme_light() +
  theme(plot.title = element_text(size = 18), 
        plot.subtitle = element_text(size = 10), 
        legend.position="right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 14),
        panel.border = element_blank()) +                 
  guides(fill = guide_colourbar(title = "Spearman's rho",
                                raster = FALSE, nbin = 20,
                                title.position = "top",                  
                                barwidth = 0.5,                               
                                title.hjust = 0.5)) +
  scale_fill_gradient2(high = "red3", mid = "cornsilk1", midpoint = 0.5,
                       limits=c(0.5, 1), breaks=seq(0.5, 1, by = 0.1), 
                       na.value = "grey") +
  scale_x_discrete(labels = c("F. prausnitzii [ref_mOTU_v25_06108]",
                              "F. prausnitzii [ref_mOTU_v25_06110]")) +
  ggpubr::rotate_x_text(angle = 0) +
  labs(x = "", y = "", title = "KEGG brite level 2")

# Individual brites
all_brite_long_finalsummary <- all_brite_long_final %>% 
  mutate(Brite_description = paste0(Brite, ":", Brite_description)) %>% 
  group_by(mOTU, Brite_description) %>% 
  summarise(Mean_rho = mean(Spearman_rho))

all_brite_long_final_wide <- all_brite_long_finalsummary %>% 
  pivot_wider(names_from = Brite_description, values_from = Mean_rho) %>% 
  column_to_rownames("mOTU") %>% 
  t() %>% 
  as.data.frame()
all_brite_long_final_wide[is.na(all_brite_long_final_wide)] <- 0

pd <- as.data.frame(scale((all_brite_long_final_wide)))
ord <- hclust(dist(pd, method = "manhattan"), method = "mcquitty")$order
ord

all_brite_long_final_wide_plotting <- all_brite_long_finalsummary

NA_df <- data.frame(mOTU = c(rep("Faecalibacterium prausnitzii [ref_mOTU_v25_06108]", times = sum(all_brite_long_final_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06108]` == 0)),
                             rep("Faecalibacterium prausnitzii [ref_mOTU_v25_06110]", times = sum(all_brite_long_final_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06110]` == 0))),
                    Brite_description = c(rownames(all_brite_long_final_wide)[all_brite_long_final_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06108]` == 0],
                               rownames(all_brite_long_final_wide)[all_brite_long_final_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06110]` == 0]),
                    Mean_rho = NA)

all_brite_long_final_wide_plotting <- rbind(all_brite_long_final_wide_plotting, NA_df)
all_brite_long_final_wide_plotting$Brite_description <- factor(all_brite_long_final_wide_plotting$Brite_description , levels = rownames(pd)[ord])

ggplot(all_brite_long_final_wide_plotting, aes(mOTU, Brite_description, fill = Mean_rho)) +
  geom_tile(color = "white", size = 0.5) +
  theme_light() +
  theme(plot.title = element_text(size = 18), 
        plot.subtitle = element_text(size = 10), 
        legend.position="right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12),
        panel.border = element_blank()) +                 
  guides(fill = guide_colourbar(title = "Spearman's rho",
                                raster = FALSE, nbin = 20,
                                title.position = "top",                  
                                barwidth = 0.5,                               
                                title.hjust = 0.5)) +
  scale_fill_gradient2(high = "red3", mid = "cornsilk1", midpoint = 0.5,
                       limits=c(0.5, 1), breaks=seq(0.5, 1, by = 0.1), 
                       na.value = "grey") +
  scale_x_discrete(labels = c("F. prausnitzii [ref_mOTU_v25_06108]",
                              "F. prausnitzii [ref_mOTU_v25_06110]")) +
  ggpubr::rotate_x_text(angle = 0) +
  labs(x = "", y = "", title = "KEGG brite")


#dev.off()


