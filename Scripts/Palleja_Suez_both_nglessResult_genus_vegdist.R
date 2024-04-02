rm(list = ls())
setwd("~/Documents/Lab/Multiple_longitudinal/Microbe_resilience/")
library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(paletteer)
library(vegan)
library(cowplot)

########### Load Suez study first ###########
# Read in master table
suez_master <-read.table(file = "../Autologous_FMT/Data/Suez_stool_genus_master.txt",
                    header = T, sep = "\t", stringsAsFactors = F)
# Calculate relative abundance for Suez study
suez_abs <- suez_master[ -c(1:8)]
suez_depth <- rowSums(suez_abs, na.rm = T)
suez_genus_values <- matrix(nrow = nrow(suez_abs), ncol = ncol(suez_abs))
colnames(suez_genus_values) <- colnames(suez_abs)
for (i in 1:nrow(suez_genus_values)) {
  for (j in 1:ncol(suez_genus_values)) {
    if (suez_abs[i, j] > 0) {
      suez_genus_values[i, j] <- suez_abs[i, j]/suez_depth[i]
    } else {
      suez_genus_values[i, j] <- 0
    }
  }
}

########### Then load Palleja study ###########
# Read in master table
pal_master <- read.table("../Antibiotic_Palleja/My_ngless_result/Palleja_stool_ngelss_motus2.5_Genus_master.txt", header = T, sep = "\t")

# Calculate relative abundance for Suez study
pal_abs <- pal_master[ -c(1:3)]
pal_depth <- rowSums(pal_abs, na.rm = T)
pal_genus_values <- matrix(nrow = nrow(pal_abs), ncol = ncol(pal_abs))
colnames(pal_genus_values) <- colnames(pal_abs)
for (i in 1:nrow(pal_genus_values)) {
  for (j in 1:ncol(pal_genus_values)) {
    if (pal_abs[i, j] > 0) {
      pal_genus_values[i, j] <- pal_abs[i, j]/pal_depth[i]
    } else {
      pal_genus_values[i, j] <- 0
    }
  }
}

########### Find the common genus in the 2 studies ###########
# Check how they intersect
intersect(colnames(pal_genus_values), colnames(suez_genus_values))
setdiff(colnames(pal_genus_values), colnames(suez_genus_values))
common <- intersect(colnames(pal_genus_values), colnames(suez_genus_values))

########### Merge the two studies (baseline time point) ###########
pal_common <- pal_genus_values %>%
  as.data.frame() %>%
  dplyr::select(which(colnames(pal_genus_values) %in% common))
pal_common <- pal_common %>%
  mutate(ID = pal_master$ID, Day = pal_master$Day, Run = pal_master$run_accession, .before = 1) %>%
  filter(Day == c("0") | Day == c("4"))
pal_common <- pal_common %>%
  mutate(Study = "Palleja", .before = 1) %>% 
  mutate(Timepoint = ifelse(Day == "0", yes = "Baseline", no = "Post_antibiotics"), .before = 4)

suez_common <- suez_genus_values %>%
  as.data.frame() %>%
  dplyr::select(which(colnames(suez_genus_values)  %in% common))
suez_common <- suez_common %>%
  mutate(ID = suez_master$participant, Description =  suez_master$Description,
         Day = suez_master$day, Group = suez_master$Group, .before = 1,
         Run = suez_master$Run) %>%
  filter(Description == "naive" & Day == 7 | Description == "during_antibiotics" & Day == 6) %>%
  dplyr::select(-c(Group)) %>% 
  mutate(Study = "Suez", .before = 1)
suez_common <- suez_common %>%
  mutate(Timepoint = ifelse(Description == "naive", yes = "Baseline", no = "Post_antibiotics"), .before = 4) %>% 
  dplyr::select(-c(Description))

merged <- rbind(suez_common, pal_common)

# Only select those who have paired data at baseline and post-antibiotics
count <- merged %>% group_by(ID) %>% count() %>% filter(n == 2)
merged_paired <- merged %>% filter(ID %in% count$ID)
merged_paired_val <- merged_paired[ , -c(1:5)]

########### Distance using vegdist ###########
dist <- vegan::vegdist(merged_paired_val,  method = "euclid", na.rm = T)
dist_df <- as.data.frame(as.matrix(dist))
rownames(dist_df) <- paste0(merged_paired$ID, ".", merged_paired$Timepoint, ".", merged_paired$Run)
colnames(dist_df) <- paste0(merged_paired$ID, ".", merged_paired$Timepoint, ".", merged_paired$Run)
dist_df_long <- dist_df %>% 
  rownames_to_column("Sample1") %>% 
  pivot_longer(2:37, names_to = "Sample2", values_to = "Dist") %>% 
  filter(Dist > 0) # Remove those 0 dist ones (self-distance)

sample_data1 <- str_split_fixed(dist_df_long$Sample1, pattern = fixed("."), n = 3)
colnames(sample_data1) <- c("Sample1_id", "Sample1_time", "Sample1_run")
sample_data2 <- str_split_fixed(dist_df_long$Sample2, pattern = fixed("."), n = 3)
colnames(sample_data2) <- c("Sample2_id", "Sample2_time", "Sample2_run")

dist_df_long <- cbind(sample_data1, sample_data2, dist_df_long) 
dist_df_long_same <- subset(dist_df_long, Sample1_id == Sample2_id) %>% 
  dplyr::select(c(1:6, 9)) %>% 
  arrange(Sample1_id) %>% 
  slice(seq(1, 36, by = 2))

#write.table(dist_df_long_same, "Palleja_Suez_bothNgless_genus_euclid_vegdist.txt", sep = "\t",
#            row.names = F, col.names = T, quote = F)

# PCoA
dist <- vegan::vegdist(merged_paired_val,  method = "euclidean", na.rm = T)
all.pcoa <- stats::cmdscale(dist, k = (nrow(merged_paired_val)-1), eig=TRUE)
in_data <- as.data.frame(vegan::scores(all.pcoa, choices = c(1,2)))
in_data$Patient <- merged_paired$ID
in_data$Case <- merged_paired$Timepoint
in_data$Study <- merged_paired$Study

(g1 <- ggplot(data = in_data, aes(x = Dim1, y = Dim2, color = Study)) +
    geom_point(size = 2, alpha = 0.6) + 
    theme_classic() + 
    ggtitle("PCoA of genus level (euclidean)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_ellipse(inherit.aes = F, aes(x = Dim1, y = Dim2, color = Study)) +
    geom_line(aes(group = Patient), size = 0.2, color = "grey"))

#### Create table for machine learning
# This is distance + baseline microbe
merged_paired_baseline <- merged_paired %>% 
  filter(Timepoint == "Baseline")

dist_df_long_same_small <- dist_df_long_same %>% 
  select(1, 7)
colnames(dist_df_long_same_small)[1] <- "ID"

ML_table <- inner_join(dist_df_long_same_small, merged_paired_baseline, by = "ID") %>% 
  dplyr::select(3, everything()) %>% 
  dplyr::select(-Day) %>% 
  rename(Dist_V1V2 = Dist) # Rename colnames

# Remove those zero abundance microbes
ML_table_abundance <- ML_table %>% 
  dplyr::select(-c(1:5))
ML_table_abundance_nonzero <- ML_table_abundance %>% 
  dplyr::select(which(colSums(ML_table_abundance) > 0))

ML_table_final <- cbind(ML_table[1:5], ML_table_abundance_nonzero)

### Add shannon diversity
# Note that shannon here is different from the raw shannon in rtk result
# Here only the common genus are included
shannon <- vegan::diversity(ML_table_abundance_nonzero, index = "shannon")

ML_table_final <- ML_table_final %>% 
  mutate(Shannon = shannon, .after = 3)

#write.table(ML_table_final, "Palleja_Suez_both_nglessResult_genus_baseline_vegdist_euclid_forML.txt", sep = "\t",
#            row.names = F, col.names = T, quote = F)


### Correlation between the features and Dist_V1V2
# This is for random forest result

(cor_1 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Bifidobacterium,
                   method = "spearman"))
(g1 <- ggplot(ML_table_final, aes(x = Bifidobacterium, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Bifidobacterium at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=0.03, y=1, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_1$estimate, 3), " , p=", round(cor_1$p.value, 4))))

(cor_2 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Lactobacillus,
                   method = "spearman"))
(g2 <- ggplot(ML_table_final, aes(x = Lactobacillus, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Lactobacillus at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=0.01, y=1, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_2$estimate, 3), " , p=", round(cor_2$p.value, 4))))

(cor_3 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Anaerostipes,
                   method = "spearman"))
(g3 <- ggplot(ML_table_final, aes(x = Anaerostipes, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Anaerostipes at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=0.04, y=1, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_3$estimate, 3), " , p=", round(cor_3$p.value, 4))))

(cor_4 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Dehalococcoidales.gen..incertae.sedis,
                   method = "spearman"))
(g4 <- ggplot(ML_table_final, aes(x = Dehalococcoidales.gen..incertae.sedis, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Dehalococcoidales gen. at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=0.03, y=1, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_4$estimate, 3), " , p=", round(cor_4$p.value, 4))))

(cor_5 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Shannon,
                   method = "spearman"))
(g5 <- ggplot(ML_table_final, aes(x = Shannon, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Shannon at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=2.5, y=1, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_5$estimate, 3), " , p=", round(cor_5$p.value, 4))))

(cor_6 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Eggerthellaceae.gen..incertae.sedis,
                   method = "spearman"))
(g6 <- ggplot(ML_table_final, aes(x = Eggerthellaceae.gen..incertae.sedis, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Eggerthellaceae gen. at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=0.005, y=1, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_6$estimate, 3), " , p=", round(cor_6$p.value, 4))))

(cor_7 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Bacteroides,
                   method = "spearman"))
(g7 <- ggplot(ML_table_final, aes(x = Bacteroides, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Bacteroides at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=0.08, y=1, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_7$estimate, 3), " , p=", round(cor_7$p.value, 4))))

(cor_8 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Firmicutes.gen..incertae.sedis,
                   method = "spearman"))
(g8 <- ggplot(ML_table_final, aes(x = Firmicutes.gen..incertae.sedis, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Firmicutes gen. at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=0.08, y=1, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_8$estimate, 3), " , p=", round(cor_8$p.value, 4))))

(cor_9 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Dorea,
                   method = "spearman"))
(g9 <- ggplot(ML_table_final, aes(x = Dorea, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Dorea at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=0.02, y=1, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_9$estimate, 3), " , p=", round(cor_9$p.value, 4))))

(cor_10 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Faecalibacterium,
                   method = "spearman"))
(g10 <- ggplot(ML_table_final, aes(x = Faecalibacterium, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Faecalibacterium at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=0.07, y=1, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_10$estimate, 3), " , p=", round(cor_10$p.value, 4))))

(g1|g2)/(g3|g4)/(g5|g6) + plot_layout(guides = "collect")
ggsave("Palleja_Suez_both_ngless_result_genus_vegdist_correlation1.pdf", device = "pdf", width = 12, height = 7)

(g7|g8)/(g9|g10) + plot_layout(guides = "collect")
ggsave("Palleja_Suez_both_ngless_result_genus_vegdist_correlation2.pdf", device = "pdf", width = 12, height = 5)







