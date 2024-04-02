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
# The suez_master here is also rarefied to 190
suez_master <-read.table(file = "../Autologous_FMT/Data/Suez_stool_species_master.txt",
                    header = T, sep = "\t", stringsAsFactors = F)
suez_species_values <- suez_master[ , -c(1:8)]


########### Then load Palleja study ###########
# Read in master table
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

########### Distance using vegdist ###########
# Here I'm using bray-curtis cuz it's count data
dist <- vegan::vegdist(merged_paired_val,  method = "bray", na.rm = T)
dist_df <- as.data.frame(as.matrix(dist))
rownames(dist_df) <- paste0(merged_paired$Study, ".", merged_paired$ID, ".", merged_paired$Timepoint, ".", merged_paired$Run)
colnames(dist_df) <- paste0(merged_paired$Study, ".", merged_paired$ID, ".", merged_paired$Timepoint, ".", merged_paired$Run)
dist_df_long <- dist_df %>% 
  rownames_to_column("Sample1") %>% 
  pivot_longer(2:37, names_to = "Sample2", values_to = "Dist") %>% 
  filter(Dist > 0) # Remove those 0 dist ones (self-distance)

sample_data1 <- str_split_fixed(dist_df_long$Sample1, pattern = fixed("."), n = 4)
colnames(sample_data1) <- c("Sample1_study", "Sample1_id", "Sample1_time", "Sample1_run")
sample_data2 <- str_split_fixed(dist_df_long$Sample2, pattern = fixed("."), n = 4)
colnames(sample_data2) <- c("Sample2_study", "Sample2_id", "Sample2_time", "Sample2_run")

dist_df_long <- cbind(sample_data1, sample_data2, dist_df_long) 
dist_df_long_same <- subset(dist_df_long, Sample1_id == Sample2_id) %>% 
  dplyr::select(c(1:8, 11)) %>% 
  arrange(Sample1_id) %>% 
  slice(seq(1, 36, by = 2))

#write.table(dist_df_long_same, "Palleja_Suez_bothNgless_rarefy190_species_bray_vegdist.txt", sep = "\t",
#            row.names = F, col.names = T, quote = F)
