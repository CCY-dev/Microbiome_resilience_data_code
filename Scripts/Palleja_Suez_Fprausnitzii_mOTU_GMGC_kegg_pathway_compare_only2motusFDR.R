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

# Read the Kegg pathway table
pathway_info <- read_xlsx("/Users/Jessica/Documents/Lab/Sofia_Kegg_perl/20231018_Kegg_pathwayMap_myedit.xlsx", 
                        sheet = 2)

###### How 108 and 110 differ in their pathway pathway profile?
gmgc_108_pathway <- gmgc_108_110_final %>% 
  dplyr::select(c(1, 2, 3, 4, 13)) %>% 
  filter(mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06108]" &
         Spearman_p_fdr < 0.05 & 
         Spearman_rho > 0.95)

gmgc_108_pathway_long <- gmgc_108_pathway %>% 
  tidyr::separate_rows(KEGG_Pathway, sep = ",") %>% #Unfold the cells with multiple values
  filter(KEGG_Pathway != "") # Remove those blank rows

gmgc_110_pathway <- gmgc_108_110_final %>% 
  dplyr::select(c(1, 2, 3, 4, 13)) %>% 
  filter(mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06110]" &
           Spearman_p_fdr < 0.05 & 
           Spearman_rho > 0.95)

gmgc_110_pathway_long <- gmgc_110_pathway %>% 
  tidyr::separate_rows(KEGG_Pathway, sep = ",") %>% #Unfold the cells with multiple values
  filter(KEGG_Pathway != "") # Remove those blank rows

all_pathway_long <- rbind(gmgc_108_pathway_long, gmgc_110_pathway_long)
colnames(all_pathway_long)[5] <- "Map" 

all_pathway_long_final <- left_join(all_pathway_long, pathway_info, by = "Map") %>% 
  filter(Level1 != "") # Remove those blank rows without level info

# Plot Venn diagram of the pathways of 06108 and 06110
library(ggvenn)
#pdf(file = "20231030_Fp_mOTU_significant_cor_GMGC_Kegg_pathway_venn.pdf", width = 18, height = 12)

# Individaul Kegg pathway
for_venn_pathway <- list(Fp_mOTU_06108 = all_pathway_long_final$Map[all_pathway_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06108]"], 
                        Fp_mOTU_06110 = all_pathway_long_final$Map[all_pathway_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06110]"])

ggvenn(for_venn_pathway, stroke_size = 0, set_name_size = 5, stroke_color = "white",
       text_size = 1.8, auto_scale = T, show_percentage = F, show_elements = T,
       label_sep = "\n",
       fill_color = c("#FAC866", "#99C1FA")) +
  theme(plot.margin = unit(c(0,0,0, 0), "cm"), 
        axis.title = element_text(size = 6),
        title = element_text(size = 14)) +
  labs(x = "", y = "", title = "Individual KEGG pathway")

# Individaul Kegg pathway with description
for_venn_pathway_description <- list(Fp_mOTU_06108 = paste0(all_pathway_long_final$Map[all_pathway_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06108]"], ":", all_pathway_long_final$Map_description[all_pathway_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06108]"]), 
                       Fp_mOTU_06110 = paste0(all_pathway_long_final$Map[all_pathway_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06110]"], ":", all_pathway_long_final$Map_description[all_pathway_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06110]"]))

ggvenn(for_venn_pathway_description, stroke_size = 0, set_name_size = 5, stroke_color = "white",
       text_size = 1.8, auto_scale = T, show_percentage = F, show_elements = T,
       label_sep = "\n",
       fill_color = c("#FAC866", "#99C1FA")) +
  theme(plot.margin = unit(c(0,0,0, 0), "cm"), 
        axis.title = element_text(size = 6),
        title = element_text(size = 14)) +
  labs(x = "", y = "", title = "Individual KEGG pathway with description")

# Level 1
for_venn_level1 <- list(Fp_mOTU_06108 = all_pathway_long_final$Level1[all_pathway_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06108]"], 
                        Fp_mOTU_06110 = all_pathway_long_final$Level1[all_pathway_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06110]"])

ggvenn(for_venn_level1, stroke_size = 0, set_name_size = 5, stroke_color = "white",
       text_size = 5, auto_scale = T, show_percentage = F, 
       fill_color = c("#FAC866", "#99C1FA")) +
  theme(plot.margin = unit(c(0,0,0, 0), "cm"), 
        axis.title = element_text(size = 6),
        title = element_text(size = 14)) +
  labs(x = "", y = "", title = "KEGG pathway level 1")

# Level 2
for_venn_level2 <- list(Fp_mOTU_06108 = all_pathway_long_final$Level2[all_pathway_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06108]"], 
                        Fp_mOTU_06110 = all_pathway_long_final$Level2[all_pathway_long_final$mOTU == "Faecalibacterium prausnitzii [ref_mOTU_v25_06110]"])

ggvenn(for_venn_level2, stroke_size = 0, set_name_size = 5, stroke_color = "white",
       text_size = 2, auto_scale = T, show_percentage = F, show_elements = T,
       label_sep = "\n",
       fill_color = c("#FAC866", "#99C1FA")) +
  theme(plot.margin = unit(c(0,0,0, 0), "cm"), 
        axis.title = element_text(size = 6),
        title = element_text(size = 14)) +
  labs(x = "", y = "",  title = "Kegg pathway level 2")

#dev.off()

#### Summarize pathway up to Level1 and Level2
all_pathway_summary_level1 <- all_pathway_long_final %>%
  group_by(mOTU, Level1) %>% 
  summarise(Mean_rho = mean(Spearman_rho))

all_pathway_summary_level2 <- all_pathway_long_final %>%
  group_by(mOTU, Level2) %>% 
  summarise(Mean_rho = mean(Spearman_rho))

# Plot heatmap
#pdf(file = "20231030_Fp_mOTU_significant_cor_GMGC_pathway_heatmap.pdf", width = 12, height = 15)

# Level 1
all_pathway_summary_level1_wide <- pivot_wider(all_pathway_summary_level1, names_from = Level1, values_from = Mean_rho) %>% 
  column_to_rownames("mOTU") %>% 
  t() %>% 
  as.data.frame()
all_pathway_summary_level1_wide[is.na(all_pathway_summary_level1_wide)] <- 0

pd <- as.data.frame(scale((all_pathway_summary_level1_wide)))
ord <- hclust(dist(pd, method = "binary"), method = "median")$order
ord

all_pathway_summary_level1_plotting <- all_pathway_summary_level1
#NA_df <- data.frame(mOTU = c(rep("Faecalibacterium prausnitzii [ref_mOTU_v25_06108]", times = sum(all_pathway_summary_level1_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06108]` == 0)),
#                             rep("Faecalibacterium prausnitzii [ref_mOTU_v25_06110]", times = sum(all_pathway_summary_level1_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06110]` == 0))),
#                    Level1 = c(rownames(all_pathway_summary_level1_wide)[all_pathway_summary_level1_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06108]` == 0],
#                               rownames(all_pathway_summary_level1_wide)[all_pathway_summary_level1_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06110]` == 0]),
#                    Mean_rho = NA)

#all_pathway_summary_level1_plotting <- rbind(all_pathway_summary_level1_plotting, NA_df)
all_pathway_summary_level1_plotting$Level1 <- factor(all_pathway_summary_level1_plotting$Level1 , 
                                                    levels = rownames(pd)[ord])

ggplot(all_pathway_summary_level1_plotting, aes(mOTU, Level1, fill = Mean_rho)) +
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
  labs(x = "", y = "", title = "KEGG pathway level 1")

# Level2
# Clustered automatically by hclust
# Ref: https://stackoverflow.com/questions/25528059/cluster-data-in-heat-map-in-r-ggplot
all_pathway_summary_level2_wide <- pivot_wider(all_pathway_summary_level2, names_from = Level2, values_from = Mean_rho) %>% 
  column_to_rownames("mOTU") %>% 
  t() %>% 
  as.data.frame()
all_pathway_summary_level2_wide[is.na(all_pathway_summary_level2_wide)] <- 0

pd <- as.data.frame(scale((all_pathway_summary_level2_wide)))
ord <- hclust(dist(pd, method = "manhattan"), method = "median")$order
ord

all_pathway_summary_level2_plotting <- all_pathway_summary_level2

NA_df <- data.frame(mOTU = c(rep("Faecalibacterium prausnitzii [ref_mOTU_v25_06108]", times = sum(all_pathway_summary_level2_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06108]` == 0)),
                             rep("Faecalibacterium prausnitzii [ref_mOTU_v25_06110]", times = sum(all_pathway_summary_level2_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06110]` == 0))),
                    Level2 = c(rownames(all_pathway_summary_level2_wide)[all_pathway_summary_level2_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06108]` == 0],
                               rownames(all_pathway_summary_level2_wide)[all_pathway_summary_level2_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06110]` == 0]),
                    Mean_rho = NA)

all_pathway_summary_level2_plotting <- rbind(all_pathway_summary_level2_plotting, NA_df)
all_pathway_summary_level2_plotting$Level2 <- factor(all_pathway_summary_level2_plotting$Level2 , levels = rownames(pd)[ord])

ggplot(all_pathway_summary_level2_plotting, aes(mOTU, Level2, fill = Mean_rho)) +
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
  labs(x = "", y = "", title = "KEGG pathway level 2")

# Individual pathways
all_pathway_long_finalsummary <- all_pathway_long_final %>% 
  mutate(Map_description = paste0(Map, ":", Map_description)) %>% 
  group_by(mOTU, Map_description) %>% 
  summarise(Mean_rho = mean(Spearman_rho))

all_pathway_long_final_wide <- all_pathway_long_finalsummary %>% 
  pivot_wider(names_from = Map_description, values_from = Mean_rho) %>% 
  column_to_rownames("mOTU") %>% 
  t() %>% 
  as.data.frame()
all_pathway_long_final_wide[is.na(all_pathway_long_final_wide)] <- 0

pd <- as.data.frame(scale((all_pathway_long_final_wide)))
ord <- hclust(dist(pd, method = "canberra"), method = "median")$order
ord

all_pathway_long_final_wide_plotting <- all_pathway_long_finalsummary

NA_df <- data.frame(mOTU = c(rep("Faecalibacterium prausnitzii [ref_mOTU_v25_06108]", times = sum(all_pathway_long_final_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06108]` == 0)),
                             rep("Faecalibacterium prausnitzii [ref_mOTU_v25_06110]", times = sum(all_pathway_long_final_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06110]` == 0))),
                    Map_description = c(rownames(all_pathway_long_final_wide)[all_pathway_long_final_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06108]` == 0],
                               rownames(all_pathway_long_final_wide)[all_pathway_long_final_wide$`Faecalibacterium prausnitzii [ref_mOTU_v25_06110]` == 0]),
                    Mean_rho = NA)

all_pathway_long_final_wide_plotting <- rbind(all_pathway_long_final_wide_plotting, NA_df)
all_pathway_long_final_wide_plotting$Map_description <- factor(all_pathway_long_final_wide_plotting$Map_description , levels = rownames(pd)[ord])

ggplot(all_pathway_long_final_wide_plotting, aes(mOTU, Map_description, fill = Mean_rho)) +
  geom_tile(color = "white", size = 0.5) +
  theme_light() +
  theme(plot.title = element_text(size = 18), 
        plot.subtitle = element_text(size = 10), 
        legend.position="right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 8),
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
  labs(x = "", y = "", title = "KEGG pathway")


#dev.off()


