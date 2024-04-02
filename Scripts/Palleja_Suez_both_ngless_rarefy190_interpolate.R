rm(list = ls())
setwd("~/Documents/Lab/Multiple_longitudinal/Microbe_resilience/Interpolate_missing_timepoint/")
library(tidyverse)
library(readxl)
library(paletteer)
library(randomcoloR)

####### Sequencing depth issue in the 2 studies
#Although the Suez sample participant 704 ERR2750079 only has depth = 68, 
# it is in the probiotics group, so not included here
# All the Suez samples used here (spontaneous group) have depth between 160 and 189
# All the Palleja samples used here have depth between 170 and 189

########### Load Suez study first ###########
# Read in master table (mOTU)
suez_master <-read.table(file = "../../Autologous_FMT/Data/Suez_stool_genus_master.txt",
                         header = T, sep = "\t", stringsAsFactors = F)
suez_genus_values <- suez_master[ , -c(1:8)]


########### Then load Palleja study ###########
# Read in master table (mOTU)
pal_master <- read.table("../../Antibiotic_Palleja/My_ngless_result/Palleja_stool_ngelss_motus2.5_genus_rarefy190_master.txt", header = T, sep = "\t")
pal_genus_values <- pal_master[ , -c(1:3)]

########### Find the common genus in the 2 studies ###########
# Check how they intersect
intersect(colnames(pal_genus_values), colnames(suez_genus_values))
setdiff(colnames(pal_genus_values), colnames(suez_genus_values))
common <- intersect(colnames(pal_genus_values), colnames(suez_genus_values))

########### Merge the two studies (all time points) ###########
# Only select 
pal_common <- pal_genus_values %>%
  as.data.frame() %>%
  dplyr::select(which(colnames(pal_genus_values) %in% common))
pal_common <- pal_common %>%
  mutate(ID = pal_master$ID, Day = pal_master$Day, Run = pal_master$run_accession, .before = 1)
pal_common <- pal_common %>%
  mutate(Study = "Palleja", .before = 1)


suez_common <- suez_genus_values %>%
  as.data.frame() %>%
  dplyr::select(which(colnames(suez_genus_values)  %in% common))
suez_common <- suez_common %>%
  mutate(ID = suez_master$participant, Description =  suez_master$Description,
         Day_raw = suez_master$day, 
         Run = suez_master$Run,
         Group = suez_master$Group,
         .before = 1) 
# Only select "Spontaneous group" in Suez and exclude those day == "NA"
suez_common <- suez_common %>% 
  filter(Group == "spontaneous reconstitution group - post antibiotics" &
           !is.na(Day_raw)) %>%
  dplyr::select(-c(Group)) %>% 
  mutate(Study = "Suez", .before = 1)

# Because Suez's day number isn't in appropriate format
# Recode Suez's day number to fit Palleja's day number format (continuous)
# naive day -1, during_antibiotics day +6, during_watchful_waiting day +13
suez_common <- suez_common %>% 
  mutate(.before = 5, 
         Day = ifelse(
           Description == "naive", yes = suez_common$Day_raw - 1,
           no = ifelse(Description == "during_antibiotics", yes = suez_common$Day_raw + 6,
                       no = suez_common$Day_raw + 13)))

suez_common <- suez_common %>%
  dplyr::select(-c(Description, Day_raw))

merged <- rbind(suez_common, pal_common)

# Remove zero-abundance taxa
merged_value <- merged[ , -c(1:4)]
nonzero_value <- merged_value[ , colSums(merged_value) > 0]

# Remove ultra-sparse taxa (>= 90% are 0)
# That is, > 10% of the value should be non-zero
nonzero_nonsparse_value <- nonzero_value[ , colSums(nonzero_value > 0) > 0.1 * nrow(nonzero_value)]

merged_new <- cbind(merged[ , 1:4], nonzero_nonsparse_value)

#### Interpolation using splinefun: Palleja ####
merged_palleja <- merged_new %>% 
  filter(Study == "Palleja")

# Interpolate the values for each individual*microbiome one by one
sample_id_pal <- unique(merged_palleja$ID)
bug_name_pal <- colnames(merged_palleja)[5:ncol(merged_palleja)]

interp_pal <- list()
for (i in 1:length(sample_id_pal)) { # Loop through sample ids
  sample_result <- list()
  for (j in 1:length(bug_name_pal)) { # Loop through taxa
    sub <- merged_palleja %>% 
      filter(ID == sample_id_pal[i]) %>% 
      dplyr::select(ID, Day, bug_name_pal[j])
    my_spline_func <- splinefun(x = sub$Day, y = sub[ , 3], method = "natural") # Natural spline
    sample_result[[j]] <- my_spline_func(0:180) #Day 0 to 180
  }
  names(sample_result) <- bug_name_pal
  sample_result_combined <- dplyr::bind_rows(sample_result) %>% #Bind all dataframes together
    mutate(.before = 1, Day = 0:180, Study = "Palleja") #Day 0 to 180
  interp_pal[[i]] <- sample_result_combined
}
names(interp_pal) <- sample_id_pal
interp_pal_combined <- dplyr::bind_rows(interp_pal, .id = "Sample_id") #Bind all dataframes together

interp_pal_combined_value <- interp_pal_combined[ , -c(1:3)]
  
# Set negative values to zero because negative abundance doesn't exist
interp_pal_combined_value[interp_pal_combined_value < 0] <- 0

# Round values so that they're still count data
interp_pal_combined_value_new <- round(interp_pal_combined_value, digits = 0)
interp_pal_combined_new <- cbind(interp_pal_combined[ , 1:3], interp_pal_combined_value_new)

#write.table(interp_pal_combined_new, "Palleja_ngless_rarefy190_genus_interpolated.txt", quote = F,
#            sep = "\t", row.names = F, col.names = T)

# Calculate relative abundance for each individual on each day
interp_pal_rel_value <- interp_pal_combined_value_new/rowSums(interp_pal_combined_value_new)
rowSums(interp_pal_rel_value) # Check if all rows sum to 1
interp_pal_rel <- cbind(interp_pal_combined[ , 1:3], interp_pal_rel_value)

#write.table(interp_pal_rel, "Palleja_ngless_rarefy190_genus_interpolated_rel_abun.txt", quote = F,
#            sep = "\t", row.names = F, col.names = T)
  
#### Interpolation using splinefun: Suez ####
merged_suez <- merged_new %>% 
  filter(Study == "Suez")

# Interpolate the values for each individual*microbiome one by one
sample_id_suez <- unique(merged_suez$ID)
bug_name_suez <- colnames(merged_suez)[5:ncol(merged_suez)]

interp_suez <- list()
for (i in 1:length(sample_id_suez)) { # Loop through sample ids
  sample_result <- list()
  for (j in 1:length(bug_name_suez)) { # Loop through taxa
    sub <- merged_suez %>% 
      filter(ID == sample_id_suez[i]) %>% 
      dplyr::select(ID, Day, bug_name_suez[j])
    my_spline_func <- splinefun(x = sub$Day, y = sub[ , 3], method = "natural") # Natural spline
    sample_result[[j]] <- my_spline_func(0:69) #Day 0 to 69
  }
  names(sample_result) <- bug_name_suez
  sample_result_combined <- dplyr::bind_rows(sample_result) %>% #Bind all dataframes together
    mutate(.before = 1, Day = 0:69, Study = "Suez") #Day 0 to 69
  interp_suez[[i]] <- sample_result_combined
}
names(interp_suez) <- sample_id_suez
interp_suez_combined <- dplyr::bind_rows(interp_suez, .id = "Sample_id") #Bind all dataframes together

interp_suez_combined_value <- interp_suez_combined[ , -c(1:3)]

# Set negative values to zero because negative abundance doesn't exist
interp_suez_combined_value[interp_suez_combined_value < 0] <- 0

# Round values so that they're still count data
interp_suez_combined_value_new <- round(interp_suez_combined_value, digits = 0)
interp_suez_combined_new <- cbind(interp_suez_combined[ , 1:3], interp_suez_combined_value_new)

#write.table(interp_suez_combined_new, "Suez_ngless_rarefy190_genus_interpolated.txt", quote = F,
#            sep = "\t", row.names = F, col.names = T)

# Calculate relative abundance for each individual on each day
interp_suez_rel_value <- interp_suez_combined_value_new/rowSums(interp_suez_combined_value_new)
rowSums(interp_suez_rel_value) # Check if all rows sum to 1
interp_suez_rel <- cbind(interp_suez_combined[ , 1:3], interp_suez_rel_value)

#write.table(interp_suez_rel, "Suez_ngless_rarefy190_genus_interpolated_rel_abun.txt", quote = F,
#            sep = "\t", row.names = F, col.names = T)

#### Merge interpolated tables (relative abundance) of Pal and Suez ###
all_interp_rel <- rbind(interp_pal_rel, interp_suez_rel)
all_interp_rel_long <- all_interp_rel %>% 
  pivot_longer(4:ncol(all_interp_rel), names_to = "Taxa", values_to = "Rel_abun")

interp_pal_rel_long <- interp_pal_rel %>% 
  pivot_longer(4:ncol(interp_pal_rel), names_to = "Taxa", values_to = "Rel_abun")

interp_suez_rel_long <- interp_suez_rel %>% 
  pivot_longer(4:ncol(interp_suez_rel), names_to = "Taxa", values_to = "Rel_abun")


#### Plot Palleja
set.seed(100)
myColor <- randomcoloR::distinctColorPalette(k = length(bug_name_pal))

ggplot(interp_pal_rel_long, aes(fill=Taxa, y=Rel_abun, x=Day)) +
  facet_wrap(.~Sample_id, nrow = 4, ncol = 3, scales='free') +
  geom_bar(position="stack", stat="identity", width = 1) +
  theme_light() +
  #scale_fill_paletteer_d("Polychrome::light", direction = -1) + #only for phylum
  scale_fill_manual(values = myColor) +
  scale_x_continuous(breaks = c(0, seq(5, 180, 5))) +
  theme(legend.position = "bottom", 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_rect(fill = 'black'))
#ggsave("Palleja_ngless_rarefy190_interpolate_natural_genus.pdf", device = "pdf", width = 35, height = 25)

#### Plot Suez
ggplot(interp_suez_rel_long, aes(fill=Taxa, y=Rel_abun, x=Day)) +
  facet_wrap(.~Sample_id, nrow = 3, ncol = 3, scales='free') +
  geom_bar(position="stack", stat="identity", width = 1) +
  theme_light() +
  #scale_fill_paletteer_d("Polychrome::light", direction = -1) +  #only for phylum
  scale_fill_manual(values = myColor) +
  scale_x_continuous(breaks = c(0, seq(5, 70, 5))) +
  theme(legend.position = "bottom", 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_rect(fill = 'black'))

#ggsave("Suez_ngless_rarefy190_interpolate_natural_genus.pdf", device = "pdf", width = 30, height = 20)











