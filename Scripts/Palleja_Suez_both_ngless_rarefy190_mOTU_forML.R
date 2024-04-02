rm(list = ls())
setwd("~/Documents/Lab/Multiple_longitudinal/Microbe_resilience/")
library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(paletteer)
library(vegan)

########### Load Suez study first ###########
# Read in master table (mOTU)
suez_master <-read.table(file = "../Autologous_FMT/Data/Suez_stool_mOTU_master.txt",
                    header = T, sep = "\t", stringsAsFactors = F)
suez_genus_values <- suez_master[ , -c(1:8)]


########### Then load Palleja study ###########
# Read in master table (mOTU)
pal_master <- read.table("../Antibiotic_Palleja/My_ngless_result/Palleja_stool_ngelss_motus2.5_mOTU_rarefy190_master.txt", header = T, sep = "\t")
pal_genus_values <- pal_master[ , -c(1:3)]
  
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

########### Load distance info ###########
dist_df_long_same <-read.table(file = "Palleja_Suez_bothNgless_rarefy190_genus_bray_vegdist.txt",
                         header = T, sep = "\t", stringsAsFactors = F)

#### Create table for machine learning
# This is distance + baseline microbe
merged_paired_baseline <- merged_paired %>% 
  filter(Timepoint == "Baseline")

dist_df_long_same_small <- dist_df_long_same %>% 
  select(2, 9)
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

#### Select only the ones who showed significance in 20211110_report.key p.20
ML_table_final_value <- ML_table_final[ , -c(1:5)]

Bacteroides <- ML_table_final[ , str_which(colnames(ML_table_final), "Bacteroides")]
Bifidobacterium <- ML_table_final[ , str_which(colnames(ML_table_final), "Bifidobacterium")]
Lactobacillus <- ML_table_final[ , str_which(colnames(ML_table_final), "Lactobacillus")]
Dehalococcoidales <- ML_table_final[ , str_which(colnames(ML_table_final), "Dehalococcoidales")]
Faecalibacterium <-  ML_table_final[ , str_which(colnames(ML_table_final), "Faecalibacterium")]
Eggerthellacea <- ML_table_final[ , str_which(colnames(ML_table_final), "Eggerthellacea")]

sig_mOTU <- cbind(Bacteroides, Bifidobacterium,
                  Lactobacillus, Dehalococcoidales, Faecalibacterium, 
                  Eggerthellacea)

abundance <- colSums(sig_mOTU) #Column sum
prevalence <- colSums(sig_mOTU != 0) #Count non-zero numbers

# Remove those sig mOTU with abundance = 1 or prevalence <= 2
sig_mOTU_filtered <- sig_mOTU[ , abundance > 1 & prevalence > 2]


final_table <- cbind(ML_table[1:5], sig_mOTU_filtered)


write.table(final_table, "Palleja_Suez_significant_mOTU_baseline_vegdist_bray_forML.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

# Test their correlation with distance by using model test
N <- ncol(sig_mOTU_filtered)
model_p <- as.data.frame(matrix(NA, nrow = N, ncol = 2))
colnames(model_p) <- c("Model_p", "Coefficient")
rownames(model_p) <- colnames(sig_mOTU_filtered)
for (i in 1:N) {
  a_mic <- colnames(sig_mOTU_filtered)[i]
  f1 <- as.formula(paste0("Dist_V1V2 ~ + ", a_mic))
  m1 <- lm(data = final_table, formula = f1)
  p <- summary(m1)
  model_p[i, 1] <- p$coefficients[2, 4]
  model_p[i, 2] <- p$coefficients[2, 1]
}

# Do FDR correction separately
Bacteroides_p <- model_p %>% 
  dplyr::filter(str_detect(rownames(model_p), "Bacteroides")) %>% 
  mutate(Model_q = p.adjust(Model_p, method = "fdr"))
Bifidobacterium_p <- model_p %>% 
  dplyr::filter(str_detect(rownames(model_p), "Bifidobacterium")) %>% 
  mutate(Model_q = p.adjust(Model_p, method = "fdr"))
Lactobacillus_p <- model_p %>% 
  dplyr::filter(str_detect(rownames(model_p), "Lactobacillus")) %>% 
  mutate(Model_q = p.adjust(Model_p, method = "fdr"))
Dehalococcoidales_p <- model_p %>% 
  dplyr::filter(str_detect(rownames(model_p), "Dehalococcoidales")) %>% 
  mutate(Model_q = p.adjust(Model_p, method = "fdr"))
Faecalibacterium_p <-  model_p %>% 
  dplyr::filter(str_detect(rownames(model_p), "Faecalibacterium")) %>% 
  mutate(Model_q = p.adjust(Model_p, method = "fdr"))
Eggerthellacea_p <- model_p %>% 
  dplyr::filter(str_detect(rownames(model_p), "Eggerthellacea")) %>% 
  mutate(Model_q = p.adjust(Model_p, method = "fdr"))

final_p <- rbind(Bacteroides_p, Bifidobacterium_p, Lactobacillus_p,
                 Dehalococcoidales_p, Faecalibacterium_p, Eggerthellacea_p)

# Investigate Faecalibacterium.prausnitzii..ref_mOTU_v25_06108
# Plot
library(ggpubr)
ggplot(final_table, aes(Dist_V1V2, Faecalibacterium.prausnitzii..ref_mOTU_v25_06108.)) +
  geom_point(alpha = 0.8) +
  stat_cor(method = "spearman") +
  theme_light() +
  labs(x = "V1-V2 Distance",
       title = "F. prausnitzii ref_mOTU_v25_06108", 
       y = "abundance") +
  geom_smooth(method = "lm", se = F, color = "grey", size = 0.7)

ggplot(final_table, aes(Study, Dist_V1V2)) +
  geom_boxplot(alpha = 0.8) +
  stat_compare_means(method = "wilcox.test") +
  theme_light() +
  labs(x = "V1-V2 Distance",
       y = "F. prausnitzii ref_mOTU_v25_06108")
  


# Investigate other keystone species
# Plot
library(ggpubr)
ggplot(ML_table, aes(Dist_V1V2, Akkermansia.muciniphila..ref_mOTU_v25_03591.)) +
  geom_point(alpha = 0.8) +
  stat_cor(method = "spearman") +
  theme_light() +
  labs(x = "V1-V2 Distance",
       title = "Akkermansia.muciniphila",
       y = "abundance") +
  geom_smooth(method = "lm", se = F, color = "grey", size = 0.7)


ggplot(ML_table, aes(Dist_V1V2, Methanobrevibacter.smithii..ref_mOTU_v25_03695.)) +
  geom_point(alpha = 0.8) +
  stat_cor(method = "spearman") +
  theme_light() +
  labs(x = "V1-V2 Distance",
       title = "Methanobrevibacter.smithii",
       y = "abundance") +
  geom_smooth(method = "lm", se = F, color = "grey", size = 0.7)

ggplot(ML_table, aes(Dist_V1V2, Ruminococcus.bromii..ref_mOTU_v25_00853.)) +
  geom_point(alpha = 0.8) +
  stat_cor(method = "spearman") +
  theme_light() +
  labs(x = "V1-V2 Distance",
       title = "Ruminococcus.bromii",
       y = "abundance") +
  geom_smooth(method = "lm", se = F, color = "grey", size = 0.7)
