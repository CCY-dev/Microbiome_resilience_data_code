rm(list = ls())
setwd("~/Documents/Lab/Multiple_longitudinal/Microbe_resilience/")
library(tidyverse)
library(readxl)
library(ggrepel)
library(tvthemes)

###### Then investigate the GMGC associated with 108 or 110
finalSig_108_110 <- read.table("20231025_Palleja_Suez_Fprausnitzii_108_110_mOTU_GMGC_either_cor_sig.txt", sep = "\t",
                        header = T)

### Write only the unique GMGC names of finalSig_1081_110
#write.table(unique(finalSig_1081_110$GMGC), "20231025_Palleja_Suez_Fprausnitzii_06108_06110_GMGC_either_sig_names.txt", quote = F,
#            sep = "\t", row.names = F, col.names = F)

# Then take this 20231025_Palleja_Suez_Fprausnitzii_06108_06110_GMGC_either_sig_names.txt to server and extract the GMGCs
# This results in 20231025_Palleja_Suez_Fprausnitzii_06108_06110_GMGC_either_sig_GMGC_extract.txt

#### Read the extracted GMGC file of 108 and 110
gmgc_108_110 <- read.table("20231025_Palleja_Suez_Fprausnitzii_06108_06110_GMGC_either_sig_GMGC_extract.txt", 
                       sep = "\t", header = T)

gmgc_108_110_feature <- gmgc_108_110 %>% 
  dplyr::select(c(1, 10)) %>% 
  rename(GMGC = query_name) %>% 
  full_join(finalSig_108_110, by = c("GMGC")) %>% # Bind with the GMGC file
  dplyr::select(c(3, 1, 2, everything()))

gmgc_108_110_feature_long <- gmgc_108_110_feature %>% 
  tidyr::separate_rows(KEGG_Pathway, sep = ",") %>% #Unfold the cells with multiple values
  filter(KEGG_Pathway != "") %>% # Remove those blank rows
  filter(str_detect(pattern = "map", string = KEGG_Pathway))

colnames(gmgc_108_110_feature_long)[3] <- "Map" 

### Read the Kegg module table
feature_info <- read_xlsx("/Users/Jessica/Documents/Lab/Sofia_Kegg_perl/20231018_Kegg_pathwayMap_myedit.xlsx", 
                         sheet = 2)

# Bind them together
gmgc_108_110_feature_long_final <- left_join(gmgc_108_110_feature_long, feature_info, by = "Map") %>% 
  filter(Level1 != "") # Remove those blank rows without level info

##### Plot scattered plot
pdf(file = "20231028_Fp_mOTU_significant_cor_GMGC_andKegg_pathwayMap_scatter_plot.pdf", width = 12, height = 12)

# Individual features
gmgc_108_110_feature_GMGC_plotting <- gmgc_108_110_feature_long_final %>% 
  dplyr::select(c(1, 2, 5, 7)) %>% 
  group_by(GMGC, mOTU) %>% 
  summarise(Feature_mean = mean(Spearman_rho)) %>% 
  pivot_wider(names_from = 2, values_from = 3) 
colnames(gmgc_108_110_feature_GMGC_plotting)[2:3] <- c("Fp06108", "Fp06110")

ggplot(gmgc_108_110_feature_GMGC_plotting, aes(x = Fp06108, y = Fp06110)) +
  geom_point(color = "blue", alpha = 0.4, size = 0.3,
             position = position_jitter(width = 0, height = 0)) +
  theme_light() +
  theme(plot.title = element_text(size = 30), 
        plot.subtitle = element_text(size = 10), 
        legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 15)
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + # Add diagonal dashed line
  xlim(0, 1) +
  ylim(0, 1) +
  coord_fixed(ratio = 1) +
  mdthemes::md_theme_light() + # This is for plotting italic axis label
  labs(x = "*F. prausnitzii* [ref_mOTU_v25_06108]", 
       y = "*F. prausnitzii* [ref_mOTU_v25_06110]",
       title = "Individual GMGC")

#### Individual KEGG pathways
gmgc_108_110_feature_individual_plotting <- gmgc_108_110_feature_long_final %>%
  filter(Spearman_p_fdr < 0.05 & Spearman_rho > 0) %>% #Only keep those significant and positive correlations before binning into modules
  dplyr::select(c(1, 3, 5, 6:8)) %>% 
  group_by(Map, mOTU) %>% 
  summarise(Feature_mean = mean(Spearman_rho, na.rm = T)) %>% 
  pivot_wider(names_from = 2, values_from = 3) 
colnames(gmgc_108_110_feature_individual_plotting)[2:3] <- c("Fp06108", "Fp06110")

#Adding the Map info
gmgc_108_110_feature_individual_plotting <- left_join(gmgc_108_110_feature_individual_plotting, feature_info, by = "Map")

# Only label the Map that have large rho difference in the two mOTUs
gmgc_108_110_feature_individual_plotting$Fp06108[is.na(gmgc_108_110_feature_individual_plotting$Fp06108)] <- 0 # Replace the NAs with 0
gmgc_108_110_feature_individual_plotting$Fp06110[is.na(gmgc_108_110_feature_individual_plotting$Fp06110)] <- 0 # Replace the NAs with 0
gmgc_108_110_feature_individual_plotting$Difference <- abs(gmgc_108_110_feature_individual_plotting$Fp06108 - gmgc_108_110_feature_individual_plotting$Fp06110)
gmgc_108_110_feature_individual_plotting$Label <- ifelse(gmgc_108_110_feature_individual_plotting$Difference > 0.5,
                                                         yes = gmgc_108_110_feature_individual_plotting$Map_description, no = "")
gmgc_108_110_feature_individual_plotting$Color <- ifelse(gmgc_108_110_feature_individual_plotting$Difference > 0.5,
                                                         yes = "red", no = "blue")

ggplot(gmgc_108_110_feature_individual_plotting, aes(x = Fp06108, y = Fp06110)) +
  geom_point(alpha = 0.6, size = 0.4, aes(color = Color),
             position = position_jitter(width = 0, height = 0)) +
  geom_text_repel(aes(label = Label), size = 3,
            nudge_x = 0, nudge_y = 0.01, max.overlaps = 20,
            min.segment.length = 0.1) +
  scale_color_manual(values = c("blue", "red")) +
  theme_light() +
  theme(plot.title = element_text(size = 30), 
        plot.subtitle = element_text(size = 10), 
        legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 15)
  ) +    
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + # Add diagonal dashed line
  xlim(0, 0.85) +
  ylim(0, 0.85) +
  coord_fixed(ratio = 1) +
  mdthemes::md_theme_light() + # This is for plotting italic axis label
  labs(x = "*F. prausnitzii* [ref_mOTU_v25_06108]", 
       y = "*F. prausnitzii* [ref_mOTU_v25_06110]",
       title = "Individual KEGG pathway")

#### Level 1
gmgc_108_110_feature_level1_plotting <- gmgc_108_110_feature_long_final %>%
  filter(Spearman_p_fdr < 0.05 & Spearman_rho > 0) %>% #Only keep those significant and positive correlations before binning into modules
  dplyr::select(c(1, 5, 6)) %>% 
  group_by(Level1, mOTU) %>% 
  summarise(Feature_mean = mean(Spearman_rho, na.rm = T)) %>% 
  pivot_wider(names_from = 2, values_from = 3) 
colnames(gmgc_108_110_feature_level1_plotting)[2:3] <- c("Fp06108", "Fp06110")

ggplot(gmgc_108_110_feature_level1_plotting, aes(x = Fp06108, y = Fp06110)) +
  geom_point(alpha = 0.6, size = 0.4, color = "blue",
             position = position_jitter(width = 0, height = 0)) +
  geom_text_repel(aes(label = Level1), size = 2,
                  nudge_x = 0, nudge_y = 0.01, max.overlaps = 100) +
  scale_color_manual(values = c("blue", "red")) +
  theme_light() +
  theme(plot.title = element_text(size = 30), 
        plot.subtitle = element_text(size = 10), 
        legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 15)
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + # Add diagonal dashed line
  xlim(0, 1) +
  ylim(0, 1) +
  coord_fixed(ratio = 1) +
  mdthemes::md_theme_light() + # This is for plotting italic axis label
  labs(x = "*F. prausnitzii* [ref_mOTU_v25_06108]", 
       y = "*F. prausnitzii* [ref_mOTU_v25_06110]",
       title = "KEGG pathway level 1")
#### Level2
gmgc_108_110_feature_level2_plotting <- gmgc_108_110_feature_long_final %>%
  filter(Spearman_p_fdr < 0.05 & Spearman_rho > 0) %>% #Only keep those significant and positive correlations before binning into modules
  dplyr::select(c(1, 5, 7)) %>% 
  group_by(Level2, mOTU) %>% 
  summarise(Feature_mean = mean(Spearman_rho, na.rm = T)) %>% 
  pivot_wider(names_from = 2, values_from = 3) 
colnames(gmgc_108_110_feature_level2_plotting)[2:3] <- c("Fp06108", "Fp06110")

# Only label the Modules that have large rho difference in the two mOTUs
gmgc_108_110_feature_level2_plotting$Fp06108[is.na(gmgc_108_110_feature_level2_plotting$Fp06108)] <- 0 # Replace the NAs with 0
gmgc_108_110_feature_level2_plotting$Fp06110[is.na(gmgc_108_110_feature_level2_plotting$Fp06110)] <- 0 # Replace the NAs with 0
gmgc_108_110_feature_level2_plotting$Difference <- abs(gmgc_108_110_feature_level2_plotting$Fp06108 - gmgc_108_110_feature_level2_plotting$Fp06110)
gmgc_108_110_feature_level2_plotting$Label <- ifelse(gmgc_108_110_feature_level2_plotting$Difference > 0.5,
                                                         yes = gmgc_108_110_feature_level2_plotting$Level2, no = "")
gmgc_108_110_feature_level2_plotting$Color <- ifelse(gmgc_108_110_feature_level2_plotting$Difference > 0.5,
                                                         yes = "red", no = "blue")

ggplot(gmgc_108_110_feature_level2_plotting, aes(x = Fp06108, y = Fp06110)) +
  geom_point(alpha = 0.6, size = 1, aes(color = Color)
             ) +
  geom_text_repel(aes(label = Label), size = 5, 
                  max.overlaps = 100, min.segment.length = 0.1) +
  scale_color_manual(values = c("blue", "red"), ) +
  theme_light() +
  theme(plot.title = element_text(size = 30), 
        plot.subtitle = element_text(size = 10), 
        legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 15)
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + # Add diagonal dashed line
  #xlim(0.55, 0.8) +
  #ylim(0.55, 0.8) +
  coord_fixed(ratio = 1) +
  mdthemes::md_theme_light() + # This is for plotting italic axis label
  labs(x = "*F. prausnitzii* [ref_mOTU_v25_06108]", 
       y = "*F. prausnitzii* [ref_mOTU_v25_06110]",
       title = "KEGG pathway level 2")
dev.off()

