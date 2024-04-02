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
suez_master <-read.table(file = "../Autologous_FMT/Data/Suez_stool_species_master.txt",
                    header = T, sep = "\t", stringsAsFactors = F)
suez_species_values <- suez_master[ , -c(1:8)]


########### Then load Palleja study ###########
# Read in master table (mOTU)
pal_master <- read.table("../Antibiotic_Palleja/My_ngless_result/Palleja_stool_ngelss_motus2.5_species_rarefy190_master.txt", header = T, sep = "\t")
pal_species_values <- pal_master[ , -c(1:3)]
  
########### Find the common species in the 2 studies ###########
# Check how they intersect
intersect(colnames(pal_species_values), colnames(suez_species_values))
setdiff(colnames(pal_species_values), colnames(suez_species_values))
common <- intersect(colnames(pal_species_values), colnames(suez_species_values))

########### Merge the two studies (baseline time point) ###########
pal_common <- pal_species_values %>%
  as.data.frame() %>%
  dplyr::select(which(colnames(pal_species_values) %in% common))
pal_common <- pal_common %>%
  mutate(ID = pal_master$ID, Day = pal_master$Day, Run = pal_master$run_accession, .before = 1) %>%
  filter(Day == c("0") | Day == c("4"))
pal_common <- pal_common %>%
  mutate(Study = "Palleja", .before = 1) %>% 
  mutate(Timepoint = ifelse(Day == "0", yes = "Baseline", no = "Post_antibiotics"), .before = 4)

suez_common <- suez_species_values %>%
  as.data.frame() %>%
  dplyr::select(which(colnames(suez_species_values)  %in% common))
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
dist_df_long_same <-read.table(file = "Palleja_Suez_bothNgless_rarefy190_species_bray_vegdist.txt",
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

### Add shannon diversity
# Note that shannon here is different from the raw shannon in rtk result
# Here only the common species are included
shannon <- vegan::diversity(ML_table_abundance_nonzero, index = "shannon")

ML_table_final <- ML_table_final %>% 
  mutate(Shannon = shannon, .after = 3)

#write.table(ML_table_final, "Palleja_Suez_both_ngless_rarefy190_species_baseline_vegdist_bray_forML.txt", sep = "\t",
#            row.names = F, col.names = T, quote = F)


(cor_1 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Faecalibacterium.prausnitzii._uncultured.Faecalibacterium.sp._Faecalibacterium.prausnitzii_,
                   method = "spearman"))
(g1 <- ggplot(ML_table_final, aes(x = Faecalibacterium.prausnitzii._uncultured.Faecalibacterium.sp._Faecalibacterium.prausnitzii_, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("F.prausnitzii uncultured v.s. distance") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=4.5, y=1.2, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_1$estimate, 3), " , p=", round(cor_1$p.value, 4))))

(cor_2 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Shannon,
                   method = "spearman"))
(g2 <- ggplot(ML_table_final, aes(x = Shannon, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Shannon v.s. distance") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=2.9, y=1.2, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_2$estimate, 3), " , p=", round(cor_2$p.value, 4))))

(cor_3 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Blautia.sp.._Blautia.sp..CAG.52_uncultured.Blautia.sp._,
                   method = "spearman"))
(g3 <- ggplot(ML_table_final, aes(x = Blautia.sp.._Blautia.sp..CAG.52_uncultured.Blautia.sp._, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Blautia sp at baseline v.s. distance") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=0.5, y=1.2, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_3$estimate, 3), " , p=", round(cor_3$p.value, 4))))

(cor_4 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Ruminococcus.sp..CAG.254,
                   method = "spearman"))
(g4 <- ggplot(ML_table_final, aes(x = Ruminococcus.sp..CAG.254, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Ruminococcus sp. at baseline v.s. distance") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=1.5, y=1.2, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_4$estimate, 3), " , p=", round(cor_4$p.value, 4))))

(cor_5 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Prevotella.sp..CAG.520,
                   method = "spearman"))
(g5 <- ggplot(ML_table_final, aes(x = Prevotella.sp..CAG.520, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Prevotella sp. at baseline v.s. distance") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=1, y=1.2, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_5$estimate, 3), " , p=", round(cor_5$p.value, 4))))

(cor_6 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Dehalococcoidales.species.incertae.sedis,
                   method = "spearman"))
(g6 <- ggplot(ML_table_final, aes(x = Dehalococcoidales.species.incertae.sedis, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Dehalococcoidales sp. at baseline v.s. distance") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=5, y=1.2, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_6$estimate, 3), " , p=", round(cor_6$p.value, 4))))

(cor_7 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Lactobacillus.acidophilus,
                   method = "spearman"))
(g7 <- ggplot(ML_table_final, aes(x = Lactobacillus.acidophilus, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("L. acidophilus at baseline v.s. distance") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=1.5, y=1.2, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_7$estimate, 3), " , p=", round(cor_7$p.value, 4))))

(cor_8 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Bifidobacterium.animalis,
                   method = "spearman"))
(g8 <- ggplot(ML_table_final, aes(x = Bifidobacterium.animalis, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("B. animalis at baseline v.s. distance") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=3.5, y=1.2, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_8$estimate, 3), " , p=", round(cor_8$p.value, 4))))

(cor_9 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$X_Eubacterium_.eligens._Eubacterium.eligens.CAG.72_Candidatus.Gastranaerophilales.bacterium.HUM_19__Eubacterium_.eligens_,
                   method = "spearman"))
(g9 <- ggplot(ML_table_final, aes(x = X_Eubacterium_.eligens._Eubacterium.eligens.CAG.72_Candidatus.Gastranaerophilales.bacterium.HUM_19__Eubacterium_.eligens_, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Eubacterium eligens at baseline v.s. distance") +
    theme_light() +
    labs(x = "Abundance") + 
    labs(y = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=3.5, y=1.2, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_9$estimate, 3), " , p=", round(cor_9$p.value, 4))))

(cor_10 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$Ruminococcus.species.incertae.sedis,
                    method = "spearman"))
(g10 <- ggplot(ML_table_final, aes(x = Ruminococcus.species.incertae.sedis, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("Ruminococcus sp. at baseline v.s. distance") +
    theme_light() +
    labs(x = "Abundance") + 
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=2, y=1.2, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_10$estimate, 3), " , p=", round(cor_10$p.value, 4))))

(g1|g2)/(g3|g4)/(g5|g6) + plot_layout(guides = "collect")
#ggsave("Palleja_Suez_species_vegdist_correlation1.pdf", device = "pdf", width = 12, height = 7)

(g7|g8)/(g9|g10) + plot_layout(guides = "collect")
#ggsave("Palleja_Suez_species_vegdist_correlation2.pdf", device = "pdf", width = 12, height = 5)

# Do FDR correction
cor_results <- as.data.frame(matrix(nrow = 10, ncol = 3))
colnames(cor_results) <- c("Feature", "Spearman_rho", "Spearman_p")
cor_results$Feature <- c("Faecalibacterium.prausnitzii._uncultured.Faecalibacterium.sp._Faecalibacterium.prausnitzii_",
                         "Shannon",
                         "Blautia.sp.._Blautia.sp..CAG.52_uncultured.Blautia.sp._",
                         "Ruminococcus.sp..CAG.254",
                         "Prevotella.sp..CAG.520",
                         "Dehalococcoidales.species.incertae.sedis",
                         "Lactobacillus.acidophilus",
                         "Bifidobacterium.animalis",
                         "X_Eubacterium_.eligens._Eubacterium.eligens.CAG.72_Candidatus.Gastranaerophilales.bacterium.HUM_19__Eubacterium_.eligens_",
                         "Ruminococcus.species.incertae.sedis")
cor_results$Spearman_rho <- round(c(cor_1$estimate, cor_2$estimate, cor_3$estimate,
                              cor_4$estimate, cor_5$estimate, cor_6$estimate,
                              cor_7$estimate, cor_8$estimate, cor_9$estimate,
                              cor_10$estimate), 2)
cor_results$Spearman_p <- c(cor_1$p.value, cor_2$p.value, cor_3$p.value,
                              cor_4$p.value, cor_5$p.value, cor_6$p.value,
                              cor_7$p.value, cor_8$p.value, cor_9$p.value,
                              cor_10$p.value)
cor_results$Spearman_q <- p.adjust(cor_results$Spearman_p, method = "fdr")
cor_results$Significance <- ifelse(cor_results$Spearman_q < 0.1, 
                                   yes = "Sig", no = "")

### Are the features aside from Shannon driven by Shannon?
# Test this by using model test
# Select the module in the top 10 features of random forest

top <- c("M00240", "M00373", "M00425", "M00227", "M00342",
         "M00715", "M00647", "M00704", "M00443", "M00646")


N = length(top)
model_p <- as.data.frame(matrix(NA, nrow = N, ncol = 2))
colnames(model_p) <- c("Model_p", "Coefficient")
rownames(model_p) <-  top
for (i in 1:N) {
  a_mic <- top[i]
  f1 <- as.formula(paste0("Dist_V1V2 ~ ", a_mic))
  m1 <- lm(data = ML_table_final, formula = f1)
  p <- summary(m1)
  model_p[i, 1] <- p$coefficients[2, 4]
  model_p[i, 2] <- p$coefficients[2, 1]
}

model_p$Model_q <- p.adjust(model_p$Model_p, method = "fdr")


# Below is to check the modules that significantly increased in Palleja's paper
# See if they also increased here
########### Merge the two studies (all time points, day 0 and 4) ###########
library(orddom)
library(ggpubr)

sig_Pal_paper <- c("M00227", "M00240", "M00373", "M00646", "M00647")

result_df <- as.data.frame(matrix(nrow = length(sig_Pal_paper), ncol = 2))
colnames(result_df) <- c("Wilcox_p", "Cliffs_delta")
rownames(result_df) <- sig_Pal_paper

for (i in 1:length(sig_Pal_paper)) {
  sub <- merged[ , c(3, which(colnames(merged) == sig_Pal_paper[i]))]
  colnames(sub)[2] <- "Value"
  sub1 <- sub[sub$Timepoint == "Baseline", ]
  sub2 <- sub[sub$Timepoint == "Post_antibiotics", ]
  wil_result <- wilcox.test(sub1$Value, sub2$Value)
  cliff_result <- dmes(x = sub1$Value, y = sub2$Value)
  
  result_df$Wilcox_p[i] <- wil_result$p.value
  result_df$Cliffs_delta[i] <- cliff_result$dc
}
result_df$Wicox_q <- p.adjust(result_df$Wilcox_p, method = "fdr")

ggplot(merged, aes(Timepoint, M00647)) +
  facet_grid(.~Study) +
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0.1, color = "blue", alpha = 0.5) +
  stat_compare_means(method = "wilcox.test", paired = F) +
  theme_light() +
  ggtitle("M00647")

