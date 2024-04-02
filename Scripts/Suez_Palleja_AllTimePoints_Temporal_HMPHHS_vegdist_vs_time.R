rm(list = ls())
setwd("~/Documents/Lab/Multiple_longitudinal/Microbe_resilience/Microbiome_drift/")
library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(vegan)

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

#### Plot time v.s. distance (all samples)
merge_dist$Day <- as.numeric(merge_dist$Day)
merge_dist$Distance <- as.numeric(merge_dist$Distance)

pdf("Suez_Palleja_Temporal_HMPHHS_shift_drift_all_timepoints_species.pdf", width = 8, height = 6)
# Raw data
ggplot(merge_dist, aes(Day^(1/2), Distance)) +
  geom_point(aes(color = Study), alpha = 0.8) +
  scale_color_manual(values = c("Temporal_EMBL" = "gray54" ,
                                   "Alien_before_antibiotics" = "black",
                                   "HMP_HHS" = "darkgoldenrod1",
                                   "Suez_after_antibiotics" = "yellowgreen",
                                   "Suez_before_antibiotics" = "mediumturquoise",
                                   "Alien_after_antibiotics"  = "coral",
                                   "Palleja"= "cornflowerblue")) +
  geom_line(aes(group = Individual, color = Study), alpha = 0.5, size = 1) +
  theme_light() +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 9)) +
  guides(color = guide_legend(title="Group")) +
  labs(x = "Day (square-rooted)",
       #title = "Distance vs day",
       y = "Microbiome dissimilarity")

# Raw data & geom_smooth for group
ggplot(merge_dist, aes(Day^(1/2), Distance)) +
  geom_point(aes(color = Study), alpha = 0.4) +
  scale_color_manual(values = c("Temporal_EMBL" = "gray54" ,
                                "Alien_before_antibiotics" = "black",
                                "HMP_HHS" = "darkgoldenrod1",
                                "Suez_after_antibiotics" = "yellowgreen",
                                "Suez_before_antibiotics" = "mediumturquoise",
                                "Alien_after_antibiotics"  = "coral",
                                "Palleja"= "cornflowerblue")) +
  geom_line(aes(group = Individual, color = Study), alpha = 0.3, size = 1) +
  theme_light() +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limit=c(0, 1), oob=scales::squish) + # Limit boundary of smooth: https://stackoverflow.com/questions/26982165/ggplot-geom-smooth-exclude-negative-values
  geom_smooth(aes(group = Study, color = Study), method = "loess", se = F, lty = "D3", size = 1.5) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 9)) +  guides(color = guide_legend(title="Group")) +
  labs(x = "Day (square-rooted)",
       #title = "Distance vs day",
       y = "Microbiome dissimilarity")

# Geom_smooth for group
ggplot(merge_dist, aes(Day^(1/2), Distance)) +
  #geom_point(aes(color = Study), alpha = 0.4) +
  scale_color_manual(values = c("Temporal_EMBL" = "gray54" ,
                                "Alien_before_antibiotics" = "black",
                                "HMP_HHS" = "darkgoldenrod1",
                                "Suez_after_antibiotics" = "yellowgreen",
                                "Suez_before_antibiotics" = "mediumturquoise",
                                "Alien_after_antibiotics"  = "coral",
                                "Palleja"= "cornflowerblue")) +
  scale_fill_manual(values = c("Temporal_EMBL" = "gray54" ,
                               "Alien_before_antibiotics" = "black",
                               "HMP_HHS" = "darkgoldenrod1",
                               "Suez_after_antibiotics" = "yellowgreen",
                               "Suez_before_antibiotics" = "mediumturquoise",
                               "Alien_after_antibiotics"  = "coral",
                               "Palleja"= "cornflowerblue")) +
  #geom_line(aes(group = Individual, color = Study), alpha = 0.3, size = 0.6) +
  theme_light() +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limit=c(0, 1), oob=scales::squish) + # Limit boundary of smooth: https://stackoverflow.com/questions/26982165/ggplot-geom-smooth-exclude-negative-values
  geom_smooth(aes(group = Study, color = Study), method = "loess", se = F, lty = "D3", size = 1.5) +
  geom_ribbon(stat='smooth', method = "loess", se=TRUE, alpha=0.2, level = 0.95,
              aes(fill = Study, group = Study)) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 9)) +  guides(color = guide_legend(title="Group"), fill = "none") +
  labs(x = "Day (square-rooted)",
       #title = "Distance vs day",
       y = "Microbiome dissimilarity")

# Raw data & geom_smooth for treatment
ggplot(merge_dist, aes(Day^(1/2), Distance)) +
  geom_point(aes(color = Treatment), alpha = 0.4) +
  geom_line(aes(group = Individual, color = Treatment), alpha = 0.3, size = 1) +
  theme_light() +
  #geom_ribbon(stat='smooth', method = "loess", se=TRUE, alpha=0.2, level = 0.95,
  #            aes(color = NULL, group = Treatment)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limit=c(0, 1), oob=scales::squish) +  # Limit boundary of smooth: https://stackoverflow.com/questions/26982165/ggplot-geom-smooth-exclude-negative-values
  geom_smooth(aes(group = Treatment, color = Treatment), method = "loess", se = F, lty = "D3", size = 1.5) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 9)) +  guides(color = guide_legend(title="Group")) +
  labs(x = "Day (square-rooted)",
       #title = "Distance vs day",
       y = "Microbiome dissimilarity")

# Geom_smooth for treatment
ggplot(merge_dist, aes(Day^(1/2), Distance)) +
  #geom_point(aes(color = Treatment), alpha = 0.4) +
  #geom_line(aes(group = Individual, color = Treatment), alpha = 0.3, size = 0.6) +
  theme_light() +
  geom_ribbon(stat='smooth', method = "loess", se=TRUE, alpha=0.2, level = 0.95,
              aes(fill = Treatment, group = Treatment)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limit=c(0, 1), oob=scales::squish) +  # Limit boundary of smooth: https://stackoverflow.com/questions/26982165/ggplot-geom-smooth-exclude-negative-values
  geom_smooth(aes(group = Treatment, color = Treatment), method = "loess", se = F, lty = "D3", size = 1.5) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 9)) +  guides(color = guide_legend(title="Treatment"), fill = "none") +
  labs(x = "Day (square-rooted)",
       #title = "Distance vs day",
       y = "Microbiome dissimilarity")

# Raw data & geom_smooth for treatment
ggplot(merge_dist, aes(Day^(1/2), Distance)) +
  geom_point(aes(color = Treatment), alpha = 0.4) +
  geom_line(aes(group = Individual, color = Treatment), alpha = 0.3, size = 1) +
  theme_light() +
  #geom_ribbon(stat='smooth', method = "loess", se=TRUE, alpha=0.2, level = 0.95,
  #            aes(fill = Treatment, group = Treatment)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limit=c(0, 1), oob=scales::squish) +  # Limit boundary of smooth: https://stackoverflow.com/questions/26982165/ggplot-geom-smooth-exclude-negative-values
  #geom_smooth(aes(group = Treatment, color = Treatment), method = "loess", se = F, lty = "D3", size = 1.5) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 9)) +  guides(color = guide_legend(title="Group")) +
  labs(x = "Day (square-rooted)",
       #title = "Distance vs day",
       y = "Microbiome dissimilarity")
dev.off()


### Test the significance of difference between the two smoothed curves
##Get the smoothed curve first
# https://statisticsglobe.com/extract-stat_smooth-regression-line-fit-from-ggplot2-plot-r
# Extract information about plot
a = ggplot(merge_dist, aes(Day^(1/2), Distance)) +
  theme_light() +
  geom_ribbon(stat='smooth', method = "loess", se=TRUE, alpha=0.2, level = 0.95,
              aes(fill = Treatment, group = Treatment)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limit=c(0, 1), oob=scales::squish) +  # Limit boundary of smooth: https://stackoverflow.com/questions/26982165/ggplot-geom-smooth-exclude-negative-values
  geom_smooth(aes(group = Treatment, color = Treatment), method = "loess", se = F, lty = "D3") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(color = guide_legend(title="Treatment"), fill = "none") +
  labs(x = "Day (square-rooted)",
       title = "Distance v.s. day") 

ggp_data <- ggplot_build(a)$data[[2]] 
head(ggp_data, 15)

smoothed_curve <- ggp_data[ , c(2, 3, 5)]
colnames(smoothed_curve) <- c("Day_sqrt", "Distance", "Treatment")
smoothed_curve$Treatment <- ifelse(smoothed_curve$Treatment == 1, 
                                   yes = "antibiotics_treated", no = "non_treated")
smoothed_curve$Treatment <- factor(smoothed_curve$Treatment, levels = c("antibiotics_treated", "non_treated"))

# Fit the data (here using gaussian because Gamma or inverse.gaussian can't take data with 0 or negative value)
my_model_gaussian <- glm(Distance ~ Treatment*Day_sqrt, family = gaussian(), data = smoothed_curve)
model_summary <- summary(my_model_gaussian)
print(model_summary)

# Generate a new data frame for predictions
newdata <- expand.grid(Day_sqrt=seq(min(smoothed_curve$Day_sqrt), max(smoothed_curve$Day_sqrt), length.out=100),
                       Treatment=levels(smoothed_curve$Treatment))

# Predict using the glm model
newdata$Predicted <- predict(my_model_gaussian, newdata=newdata, type="response")

## Plot the predicted curves
# Original curve
b = ggplot(smoothed_curve, aes(Day_sqrt, Distance)) +
  theme_light() +
  geom_point(aes(color = Treatment), alpha = 0.4) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(color = guide_legend(title="Treatment"), fill = "none") +
  labs(x = "Day (square-rooted)",
       title = "Distance v.s. day") 

# Add predicted data
b + geom_line(data = newdata, aes(Day_sqrt, Predicted, color = Treatment),
              size = 1.5, linetype = 5)
#ggsave("Suez_Palleja_Temporal_HMPHHS_shift_drift_all_timepoints_species_fitted.pdf", width = 8, height = 6)

