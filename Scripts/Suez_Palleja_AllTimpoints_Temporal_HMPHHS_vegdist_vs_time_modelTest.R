rm(list = ls())
setwd("~/Documents/Lab/Multiple_longitudinal/Microbe_resilience/Microbiome_drift/")
library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(vegan)
library(lmtest)
library(lme4)

########### Load all study dist first ###########
# Read in master table
temporal_dist <-read.table(file = "../../Temporal_variation_EMBL/analysis/Temporal_species_bray_vegdist_alien_separated.txt",
                           header = T, sep = "\t", stringsAsFactors = F)
HHS_dist <- read.table(file = "../../HealthyHumanStudy_HMP/analysis/HMPHHS_species_bray_vegdist.txt",
                       header = T, sep = "\t", stringsAsFactors = F)
pal_dist <- read.table(file = "../../Antibiotic_Palleja/Analysis_basedOn_ngless_result/Palleja_species_bray_allTimePoints_vegdist.txt",
                       header = T, sep = "\t", stringsAsFactors = F)
suez_post_anti_dist <- read.table(file = "../../Autologous_FMT/Analysis/Suez_species_after_antibiotics_bray_vegdist.txt",
                                  header = T, sep = "\t", stringsAsFactors = F)
suez_naive_dist <- read.table(file = "../../Autologous_FMT/Analysis/Suez_species_before_antibiotics_allGroups_bray_vegdist.txt",
                              header = T, sep = "\t", stringsAsFactors = F)


# Merge the distance tables
temporal_new <- temporal_dist[ , c(1, 4, 5)]
temporal_new$Study <- ifelse(temporal_new$Individual1 == "alien-antibiotics", 
                             yes = "Alien_after_antibiotics", 
                             no =  ifelse(temporal_new$Individual1 == "alien",
                                          yes = "Alien_before_antibiotics",
                                          no = "Temporal_EMBL"))

colnames(temporal_new) <- c("Individual", "Day", "Distance", "Study")

pal_dist$Study <- "Palleja"
suez_post_anti_dist$Study <- "Suez_after_antibiotics"
suez_naive_dist$Study <- "Suez_before_antibiotics"

suez_pal_merge <- rbind(pal_dist, suez_post_anti_dist, suez_naive_dist)
suez_pal_new <- suez_pal_merge[ , c(1, 4, 7, 8)]
colnames(suez_pal_new) <- c("Individual", "Day", "Distance", "Study")

merge_dist <- rbind(temporal_new, suez_pal_new)

# Bind HHS to others
HHS_new <- HHS_dist[ , c(1, 4, 5)]
HHS_new$Study <- "HMP_HHS"

colnames(HHS_new) <- c("Individual", "Day", "Distance", "Study")

merge_dist <- rbind(merge_dist, HHS_new)

# Add distance = 0 on day 0
#Count unique individuals by study
a = merge_dist %>% group_by(Study) %>% summarise(count = n_distinct(Individual))

baseline <- data.frame(Individual = unique(merge_dist$Individual),
                       Day = 0, Distance = 0, 
                       Study = c(rep("Temporal_EMBL", 6),
                                 rep("Alien_before_antibiotics", 1),
                                 rep("Alien_after_antibiotics", 1),
                                 rep("Palleja", 12),
                                 rep("Suez_after_antibiotics", 16),
                                 rep("Suez_before_antibiotics", 21),
                                 rep("HMP_HHS", 55)))

merge_dist <- rbind(baseline, merge_dist)
merge_dist <- merge_dist %>% 
  mutate(Treatment = case_when(
    Study %in% c("Temporal_EMBL", "Alien_before_antibiotics", "HMP_HHS", "Suez_before_antibiotics") ~ "Non-treated",
    Study %in% c("Suez_after_antibiotics", "Alien_after_antibiotics", "Palleja") ~ "Antibiotics-treated"
  ))

merge_dist <- merge_dist %>% 
  mutate(Treatment_binary = 
           case_when(Treatment == "Antibiotics-treated" ~ 1,
                     Treatment == "Non-treated" ~ 0))


# Remove 'Day = 0' for plotting because it was for plotting lines in another script
merge_dist_new <- merge_dist %>% 
  filter(Day != 0)

# Model test
m <- lmer(formula = Distance ~ Treatment_binary + (1|Individual),
          data = merge_dist_new, REML = T)
result <- car::Anova(m, type = 2, test.statistic = "F")
result

# Plot box plots
library(ggpubr)

merge_dist_new$Treatment <- factor(merge_dist_new$Treatment, levels = c("Non-treated", "Antibiotics-treated"))


pdf("Suez_Palleja_AllTimepoints_Temporal_HMPHHS_distance_plots_species.pdf", width = 8, height = 6)

g1 <- ggplot(merge_dist_new, aes(Treatment, Distance)) +
  geom_boxplot() +
  geom_jitter(color = "salmon", alpha = 0.4, width = 0.1, height = 0.01) +
  theme_light() +
  stat_compare_means(method = "wilcox.test", label.y = 1.2) +
  labs(x = "")
g1



merge_dist_new$Study <- factor(merge_dist_new$Study,
                           levels = c("HMP_HHS", "Temporal_EMBL",
                                      "Alien_before_antibiotics",
                                      "Suez_before_antibiotics",
                                      "Alien_after_antibiotics",
                                      "Suez_after_antibiotics",
                                      "Palleja"))

g2 <- ggplot(merge_dist_new, aes(Study, Distance)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.6, width = 0.1, height = 0.01, color = "salmon") +
  theme_light() +
  stat_compare_means(method = "kruskal.test", label.y = 1.2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") +
  labs(x = "")
g2
dev.off()
  



