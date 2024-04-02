rm(list = ls())
setwd("~/Documents/Lab/Multiple_longitudinal/Microbe_resilience/Interpolate_missing_timepoint/")
library(tidyverse)
library(readxl)
library(paletteer)
library(vegan)

########### Load Suez study first ###########
# Read in master table (interpolated)
suez_master <-read.table(file = "Suez_ngless_rarefy190_species_interpolated.txt",
                         header = T, sep = "\t", stringsAsFactors = F)
suez_master <- suez_master %>% 
  mutate(.before = 1, Sample = paste0(Sample_id, "_", Day))
suez_values <- suez_master[ , -c(1:4)]

########### Then load Palleja study ###########
# Read in master table (interpolated)
pal_master <- read.table("Palleja_ngless_rarefy190_species_interpolated.txt",
                         header = T, sep = "\t", stringsAsFactors = F)
pal_master <- pal_master %>% 
  mutate(.before = 1, Sample = paste0(Sample_id, "_", Day))
pal_values <- pal_master[ , -c(1:4)]

########### Distance using vegdist: Suez ###########
dist_suez <- vegan::vegdist(suez_values,  method = "bray", na.rm = T)
dist_df_suez <- as.data.frame(as.matrix(dist_suez))
rownames(dist_df_suez) <- suez_master$Sample
colnames(dist_df_suez) <- suez_master$Sample
dist_df_long_suez <- dist_df_suez %>% 
  rownames_to_column("Sample1") %>% 
  pivot_longer(2:491, names_to = "Sample2", values_to = "Dist") %>% 
  filter(Dist > 0) # Remove those 0 dist ones (self-distance)

sample_data1_suez <- str_split_fixed(dist_df_long_suez$Sample1, pattern = "_", n = 2)
colnames(sample_data1_suez) <- c("Individual1", "Day1")
sample_data2_suez <- str_split_fixed(dist_df_long_suez$Sample2, pattern = "_", n = 2)
colnames(sample_data2_suez) <- c("Individual2", "Day2")

dist_df_long_bind_suez <- cbind(sample_data1_suez, sample_data2_suez, dist_df_long_suez) 

dist_df_long_bind_suez[ , c(2, 4, 7)] <- apply(dist_df_long_bind_suez[ , c(2, 4, 7)] , 2, as.numeric)

dist_df_long_bind_suez_new <- dist_df_long_bind_suez %>% 
  filter(Individual1 == Individual2 & Day1 == 0) #day 0 as reference for all samples


########### Distance using vegdist: Palleja ###########
dist_pal <- vegan::vegdist(pal_values,  method = "bray", na.rm = T)
dist_df_pal <- as.data.frame(as.matrix(dist_pal))
rownames(dist_df_pal) <- pal_master$Sample
colnames(dist_df_pal) <- pal_master$Sample
dist_df_long_pal <- dist_df_pal %>% 
  rownames_to_column("Sample1") %>% 
  pivot_longer(2:2173, names_to = "Sample2", values_to = "Dist") %>% 
  filter(Dist > 0) # Remove those 0 dist ones (self-distance)

sample_data1_pal <- str_split_fixed(dist_df_long_pal$Sample1, pattern = "_", n = 2)
colnames(sample_data1_pal) <- c("Individual1", "Day1")
sample_data2_pal <- str_split_fixed(dist_df_long_pal$Sample2, pattern = "_", n = 2)
colnames(sample_data2_pal) <- c("Individual2", "Day2")

dist_df_long_bind_pal <- cbind(sample_data1_pal, sample_data2_pal, dist_df_long_pal) 

dist_df_long_bind_pal[ , c(2, 4, 7)] <- apply(dist_df_long_bind_pal[ , c(2, 4, 7)] , 2, as.numeric)

dist_df_long_bind_pal_new <- dist_df_long_bind_pal %>% 
  filter(Individual1 == Individual2 & Day1 == 0) #day 0 as reference for all samples

#### Find out the lowest dist for each individual: Suez ####
suez_indi <- unique(dist_df_long_bind_suez_new$Individual1)

low_dist_suez <- as.data.frame(matrix(nrow = length(suez_indi),
                                      ncol = 3))
colnames(low_dist_suez) <- c("Sample_id", "recovered_day", "lowest_dist")
for (i in 1:length(suez_indi)) {
  # Suez recovery day starts on day 14
  sub <- dist_df_long_bind_suez_new %>% 
    filter(Individual1 == suez_indi[i] & Day2 >= 14)
  min_dist <- min(sub$Dist)

  low_dist_suez$Sample_id[i] <- as.numeric(unique(sub$Individual1))
  low_dist_suez$lowest_dist[i] <- min_dist
  low_dist_suez$recovered_day[i] <- sub$Day2[which(sub$Dist == min_dist)] #The day that has lowest dist
}

low_dist_suez$DistxDay <- low_dist_suez$recovered_day * low_dist_suez$lowest_dist

ggplot(low_dist_suez, aes(recovered_day, lowest_dist)) +
  geom_point() +
  theme_light()

# Merge with suez_master
suez_master_baseline <- suez_master %>% 
  filter(Day == 0)
  
suez_for_ML <- inner_join(low_dist_suez, suez_master_baseline, by = "Sample_id")

# Add scale of perturbance (BC-dist between day 0 and last day of antibiotics, day 13)
suez_perturb <- dist_df_long_bind_suez_new %>% 
  filter(Day2 == 13)

suez_for_ML <- suez_for_ML %>% 
  mutate(Perturbance_scale = suez_perturb$Dist, .after = 4)

#### Find out the lowest dist for each individual: Palleja ####
pal_indi <- unique(dist_df_long_bind_pal_new$Individual1)

low_dist_pal <- as.data.frame(matrix(nrow = length(pal_indi),
                                      ncol = 3))
colnames(low_dist_pal) <- c("Sample_id", "recovered_day", "lowest_dist")
for (i in 1:length(pal_indi)) {
  # pal recovery day starts on day 5
  sub <- dist_df_long_bind_pal_new %>% 
    filter(Individual1 == pal_indi[i] & Day2 >= 5)
  min_dist <- min(sub$Dist)
  
  low_dist_pal$Sample_id[i] <- unique(sub$Individual1)
  low_dist_pal$lowest_dist[i] <- min_dist
  low_dist_pal$recovered_day[i] <- sub$Day2[which(sub$Dist == min_dist)] #The day that has lowest dist
}

low_dist_pal$DistxDay <- low_dist_pal$recovered_day * low_dist_pal$lowest_dist

ggplot(low_dist_pal, aes(recovered_day, lowest_dist)) +
  geom_point() +
  theme_light()

# Merge with pal_master
pal_master_baseline <- pal_master %>% 
  filter(Day == 0)

pal_for_ML <- inner_join(low_dist_pal, pal_master_baseline, by = "Sample_id")

# Add scale of perturbance (BC-dist between day 0 and last day of antibiotics, day 4)
pal_perturb <- dist_df_long_bind_pal_new %>% 
  filter(Day2 == 4)

pal_for_ML <- pal_for_ML %>% 
  mutate(Perturbance_scale = pal_perturb$Dist, .after = 4)

#### Calculate the Shannon diversity #####
suez_for_ML_value <- suez_for_ML[ , 9:ncol(suez_for_ML)]
pal_for_ML_value <- pal_for_ML[ , 9:ncol(pal_for_ML)]

shannon_suez <- vegan::diversity(suez_for_ML_value, index = "shannon")
shannon_pal <- vegan::diversity(pal_for_ML_value, index = "shannon")

#### Merge the two studies' ML tables #####
common_taxa <- intersect(colnames(suez_for_ML_value), colnames(pal_for_ML_value))
suez_for_ML_value_common <- suez_for_ML_value[, colnames(suez_for_ML_value) %in% common_taxa]
pal_for_ML_value_common <- pal_for_ML_value[, colnames(pal_for_ML_value) %in% common_taxa]

suez_for_ML_common <- cbind(suez_for_ML[ , 1:8], shannon_suez, suez_for_ML_value_common)
pal_for_ML_common <- cbind(pal_for_ML[ , 1:8], shannon_pal, pal_for_ML_value_common)

colnames(suez_for_ML_common)[9] <- "Shannon"
colnames(pal_for_ML_common)[9] <- "Shannon"

merged_for_ML <- rbind(suez_for_ML_common, pal_for_ML_common)

## Remove values that are all zero
merged_for_ML_value <- merged_for_ML[ , 10:ncol(merged_for_ML)]
merged_for_ML_value_nonzero <- merged_for_ML_value[ , colSums(merged_for_ML_value) > 0]

merged_for_ML_nonzero <- cbind(merged_for_ML[ , 1:9], merged_for_ML_value_nonzero)


#write.table(merged_for_ML_nonzero, "Palleja_Suez_both_ngless_rarefy190_species_baseline_vegdist_bray_interpolated_forML.txt", sep = "\t",
#            row.names = F, col.names = T, quote = F)


#### Plot the correlations
(cor_1 <- cor.test(merged_for_ML_nonzero$DistxDay, merged_for_ML_nonzero$X_Eubacterium_.rectale._Eubacterium.rectale.CAG.36_Eubacterium.sp..41_20__Eubacterium_.rectale_,
                   method = "spearman"))
(g1 <- ggplot(merged_for_ML_nonzero, aes(x = X_Eubacterium_.rectale._Eubacterium.rectale.CAG.36_Eubacterium.sp..41_20__Eubacterium_.rectale_, y = DistxDay)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Eubacterium rectale. v.s. DistxDay") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=6, y=120, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_1$estimate, 3), " , p=", round(cor_1$p.value, 4))))

(cor_2 <- cor.test(merged_for_ML_nonzero$DistxDay, merged_for_ML_nonzero$Ruminococcus.bromii._Ruminococcus.sp..CAG.108_Ruminococcus.sp..A254.MGS.108_uncultured.Ruminococcus.sp._Ruminococcus.sp..CAG.108.related_41_35_Ruminococcus.bromii_,
                   method = "spearman"))
(g2 <- ggplot(merged_for_ML_nonzero, aes(x = Ruminococcus.bromii._Ruminococcus.sp..CAG.108_Ruminococcus.sp..A254.MGS.108_uncultured.Ruminococcus.sp._Ruminococcus.sp..CAG.108.related_41_35_Ruminococcus.bromii_, y = DistxDay)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Ruminococcus bromii v.s. DistxDay") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=3, y=130, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_2$estimate, 3), " , p=", round(cor_2$p.value, 4))))

(cor_3 <- cor.test(merged_for_ML_nonzero$DistxDay, merged_for_ML_nonzero$Clostridiales.sp.._Clostridium.sp..CAG.264_Clostridium.sp..L2.50_uncultured.Coprococcus.sp._,
                   method = "spearman"))
(g3 <- ggplot(merged_for_ML_nonzero, aes(x = Clostridiales.sp.._Clostridium.sp..CAG.264_Clostridium.sp..L2.50_uncultured.Coprococcus.sp._, y = DistxDay)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Clostridiales sp. at baseline v.s. DistxDay") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=1, y=130, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_3$estimate, 3), " , p=", round(cor_3$p.value, 4))))

(cor_4 <- cor.test(merged_for_ML_nonzero$DistxDay, merged_for_ML_nonzero$Prevotella.copri._Prevotella.copri.CAG.164_Prevotella.copri_,
                   method = "spearman"))
(g4 <- ggplot(merged_for_ML_nonzero, aes(x = Prevotella.copri._Prevotella.copri.CAG.164_Prevotella.copri_, y = DistxDay)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Prevotella copri at baseline v.s. DistxDay") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=40, y=140, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_4$estimate, 3), " , p=", round(cor_4$p.value, 4))))

(cor_5 <- cor.test(merged_for_ML_nonzero$DistxDay, merged_for_ML_nonzero$Faecalibacterium.prausnitzii._uncultured.Faecalibacterium.sp._Faecalibacterium.prausnitzii_,
                   method = "spearman"))
(g5 <- ggplot(merged_for_ML_nonzero, aes(x = Faecalibacterium.prausnitzii._uncultured.Faecalibacterium.sp._Faecalibacterium.prausnitzii_, y = DistxDay)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("F. prausnitzii uncultured at baseline v.s. DistxDay") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=6, y=140, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_5$estimate, 3), " , p=", round(cor_5$p.value, 4))))

(cor_6 <- cor.test(merged_for_ML_nonzero$DistxDay, merged_for_ML_nonzero$Perturbance_scale,
                   method = "spearman"))
(g6 <- ggplot(merged_for_ML_nonzero, aes(x = Perturbance_scale, y = DistxDay)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Perturbance_scale v.s. DistxDay") +
    theme_light() +
    labs(x = "Perturbance scale") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=0.7, y=130, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_6$estimate, 3), " , p=", round(cor_6$p.value, 4))))

(cor_7 <- cor.test(merged_for_ML_nonzero$DistxDay, merged_for_ML_nonzero$Roseburia.sp.._Roseburia.sp..CAG.45_uncultured.Roseburia.sp._,
                   method = "spearman"))
(g7 <- ggplot(merged_for_ML_nonzero, aes(x = Roseburia.sp.._Roseburia.sp..CAG.45_uncultured.Roseburia.sp._, y = DistxDay)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Roseburia sp. at baseline v.s. DistxDay") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=1.5, y=130, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_7$estimate, 3), " , p=", round(cor_7$p.value, 4))))

(cor_8 <- cor.test(merged_for_ML_nonzero$DistxDay, merged_for_ML_nonzero$Dialister.sp..CAG.486,
                   method = "spearman"))
(g8 <- ggplot(merged_for_ML_nonzero, aes(x = Dialister.sp..CAG.486, y = DistxDay)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Dialister sp. at baseline v.s. DistxDay") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=6, y=130, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_8$estimate, 3), " , p=", round(cor_8$p.value, 4))))

(cor_9 <- cor.test(merged_for_ML_nonzero$DistxDay, merged_for_ML_nonzero$Bacteroides.faecis_thetaiotaomicron._Bacteroides.sp..HMSC067B03_Bacteroides.sp..HMSC068A09_Bacteroides.sp..1_1_14_Bacteroides.faecis_Bacteroides.thetaiotaomicron_Bacteroides.sp..AR20_,
                   method = "spearman"))
(g9 <- ggplot(merged_for_ML_nonzero, aes(x = Bacteroides.faecis_thetaiotaomicron._Bacteroides.sp..HMSC067B03_Bacteroides.sp..HMSC068A09_Bacteroides.sp..1_1_14_Bacteroides.faecis_Bacteroides.thetaiotaomicron_Bacteroides.sp..AR20_, y = DistxDay)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Bacteroides faecis/thetaiotaomicron at baseline v.s. DistxDay") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=2, y=130, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_9$estimate, 3), " , p=", round(cor_9$p.value, 4))))

(cor_10 <- cor.test(merged_for_ML_nonzero$DistxDay, merged_for_ML_nonzero$Faecalibacterium.sp.._Faecalibacterium.sp..CAG.74_Faecalibacterium.sp..CAG.74_58_120_,
                    method = "spearman"))
(g10 <- ggplot(merged_for_ML_nonzero, aes(x = Faecalibacterium.sp.._Faecalibacterium.sp..CAG.74_Faecalibacterium.sp..CAG.74_58_120_, y = DistxDay)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Faecalibacterium sp. at baseline v.s. DistxDay") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=3, y=130, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_10$estimate, 3), " , p=", round(cor_10$p.value, 4))))

(g1|g2)/(g3|g4)/(g5|g6) + plot_layout(guides = "collect")
#ggsave("Palleja_Suez_rarefy190_species_vegdist_species_interpolated_correlation1.pdf", device = "pdf", width = 12, height = 7)

(g7|g8)/(g9|g10) + plot_layout(guides = "collect")
#ggsave("Palleja_Suez_rarefy190_species_vegdist_species_interpolated_correlation2.pdf", device = "pdf", width = 12, height = 5)


