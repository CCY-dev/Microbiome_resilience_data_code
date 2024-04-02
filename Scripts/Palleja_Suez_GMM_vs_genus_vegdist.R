rm(list = ls())
setwd("~/Documents/Lab/Multiple_longitudinal/Microbe_resilience/")
library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(paletteer)
library(vegan)

########### Load Suez study first ###########
# Read in master table
suez_master <-read.table(file = "../Autologous_FMT/Data/Suez_stool_GMM_master.txt",
                         header = T, sep = "\t", stringsAsFactors = F)
suez_module_values <- suez_master[ , -c(1:8)]

########### Then load Palleja study ###########
# Read in master table
pal_master <- read.table("../Antibiotic_Palleja/My_ngless_result/Palleja_stool_GMM_master.txt", header = T, sep = "\t")
pal_module_values <- pal_master[ , -c(1:3)]

########### Find the common module in the 2 studies ###########
# Check how they intersect
intersect(colnames(pal_module_values), colnames(suez_module_values))
setdiff(colnames(pal_module_values), colnames(suez_module_values))
common <- intersect(colnames(pal_module_values), colnames(suez_module_values))

########### Merge the two studies (baseline time point) ###########
pal_common <- pal_module_values %>%
  as.data.frame() %>%
  dplyr::select(which(colnames(pal_module_values) %in% common))
pal_common <- pal_common %>%
  mutate(ID = pal_master$ID, Day = pal_master$Day, Run = pal_master$run_accession, .before = 1) %>%
  filter(Day == c("0") | Day == c("4"))
pal_common <- pal_common %>%
  mutate(Study = "Palleja", .before = 1) %>% 
  mutate(Timepoint = ifelse(Day == "0", yes = "Baseline", no = "Post_antibiotics"), .before = 4)

suez_common <- suez_module_values %>%
  as.data.frame() %>%
  dplyr::select(which(colnames(suez_module_values)  %in% common))
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

########### Read distance table ###########
# Here I'm using bray-curtis cuz it's count data
distance <-read.table(file = "/Users/Jessica/Documents/Lab/Multiple_longitudinal/Microbe_resilience/Palleja_Suez_bothNgless_rarefy190_genus_bray_vegdist.txt",
                         header = T, sep = "\t", stringsAsFactors = F)


#### Create table for machine learning
# This is distance + baseline microbe
merged_paired_baseline <- merged %>% 
  filter(Timepoint == "Baseline")

dist_df_long_same_small <- distance %>% 
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


#write.table(ML_table_final, "Palleja_Suez_GMM_predict_genus_vegdist_bray_forML.txt", sep = "\t",
#            row.names = F, col.names = T, quote = F)

### Correlation between the features and Dist_V1V2
# This is for random forest result

(cor_1 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$MF0099,
                   method = "spearman"))
(g1 <- ggplot(ML_table_final, aes(x = MF0099, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("MF0099 at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=2, y=1.2, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_1$estimate, 3), " , p=", round(cor_1$p.value, 4))))

(cor_2 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$MF0048,
                   method = "spearman"))
(g2 <- ggplot(ML_table_final, aes(x = MF0048, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("MF0048 at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=600, y=1.2, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_2$estimate, 3), " , p=", round(cor_2$p.value, 4))))

(cor_3 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$MF0039,
                   method = "spearman"))
(g3 <- ggplot(ML_table_final, aes(x = MF0039, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("MF0039 at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=75, y=1.2, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_3$estimate, 3), " , p=", round(cor_3$p.value, 4))))

(cor_4 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$MF0097,
                   method = "spearman"))
(g4 <- ggplot(ML_table_final, aes(x = MF0097, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("MF0097 at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=3, y=1.2, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_4$estimate, 3), " , p=", round(cor_4$p.value, 4))))

(cor_5 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$MF0116,
                   method = "spearman"))
(g5 <- ggplot(ML_table_final, aes(x = MF0116, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("MF0116 at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=55, y=1.2, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_5$estimate, 3), " , p=", round(cor_5$p.value, 4))))

(cor_6 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$MF0089,
                   method = "spearman"))
(g6 <- ggplot(ML_table_final, aes(x = MF0089, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("MF0089 at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=35, y=1.2, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_6$estimate, 3), " , p=", round(cor_6$p.value, 4))))

(cor_7 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$MF0085,
                   method = "spearman"))
(g7 <- ggplot(ML_table_final, aes(x = MF0085, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("MF0085 at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=400, y=1.2, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_7$estimate, 3), " , p=", round(cor_7$p.value, 4))))

(cor_8 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$MF0042,
                   method = "spearman"))
(g8 <- ggplot(ML_table_final, aes(x = MF0042, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("MF0042 at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=10, y=1.2, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_8$estimate, 3), " , p=", round(cor_8$p.value, 4))))

(cor_9 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$MF0022,
                   method = "spearman"))
(g9 <- ggplot(ML_table_final, aes(x = MF0022, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("MF0022 at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=32, y=1.2, size = 4,color = "red",
             label=paste0("Spearman's rho=", round(cor_9$estimate, 3), " , p=", round(cor_9$p.value, 4))))

(cor_10 <- cor.test(ML_table_final$Dist_V1V2, ML_table_final$MF0030,
                    method = "spearman"))
(g10 <- ggplot(ML_table_final, aes(x = MF0030, y = Dist_V1V2)) + 
    geom_point(alpha = 0.6, color = "cornflowerblue", size = 3, aes(shape = Study)) +
    ggtitle("MF0030 at baseline v.s. V1_V2 distance") +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    geom_smooth(method = "lm", se = F, color = "cornflowerblue", alpha = 0.8) +
    annotate(geom="text", x=280, y=1.2, size = 4,color = "black",
             label=paste0("Spearman's rho=", round(cor_10$estimate, 3), " , p=", round(cor_10$p.value, 4))))

(g1|g2)/(g3|g4)/(g5|g6) + plot_layout(guides = "collect")
#ggsave("Palleja_Suez_GMM_predict_genus_vegdist_correlation1.pdf", device = "pdf", width = 12, height = 7)

(g7|g8)/(g9|g10) + plot_layout(guides = "collect")
#ggsave("Palleja_Suez_GMM_predict_genus_vegdist_correlation2.pdf", device = "pdf", width = 12, height = 5)


#### Do FDR correction
cor_results <- as.data.frame(matrix(nrow = 10, ncol = 3))
colnames(cor_results) <- c("Feature", "Spearman_rho", "Spearman_p")
cor_results$Feature <- c("MF0099",
                         "MF0048",
                         "MF0039",
                         "MF0097",
                         "MF0116",
                         "MF0089",
                         "MF0085",
                         "MF0042",
                         "MF0022",
                         "MF0030")
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

