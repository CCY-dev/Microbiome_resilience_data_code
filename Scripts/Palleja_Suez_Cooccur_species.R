rm(list = ls())
setwd("~/Documents/Lab/Multiple_longitudinal/Microbe_resilience/")
library(tidyverse)
library(readxl)
library(cooccur)

########### Load Suez study first ###########
# Read in master table
suez_master <-read.table(file = "../Autologous_FMT/Data/Suez_stool_species_master.txt",
                    header = T, sep = "\t", stringsAsFactors = F)
suez_genus_values <- suez_master[ , -c(1:8)]

########### Then load Palleja study ###########
# Read in master table
pal_master <- read.table("../Antibiotic_Palleja/My_ngless_result/Palleja_stool_ngelss_motus2.5_Species_rarefy190_master.txt", header = T, sep = "\t")
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

# Remove those zero-abundance ones
merged_paired_val_nonzero <- merged_paired_val[ , colSums(merged_paired_val) > 0]

# Remove those low-abundant ones
# Those less than 5 non-zero values are removed
unname(sort(colSums(merged_paired_val_nonzero > 0)))
merged_paired_val_nonzero_new <- merged_paired_val_nonzero[colSums(merged_paired_val_nonzero > 0) > 4]


final_merged <- cbind(merged_paired[ , 1:5], merged_paired_val_nonzero_new)

# Check F. prausnitzii uncultured prevalence
target <- final_merged %>% 
  dplyr::select(c(2:4, which(colnames(.) == "Faecalibacterium.prausnitzii._uncultured.Faecalibacterium.sp._Faecalibacterium.prausnitzii_")))
colnames(target)[4] <- "F.p.uncultured"
target_t1 <- target %>% 
  filter(Timepoint == "Baseline")
target_t2 <- target %>% 
  filter(Timepoint == "Post_antibiotics")

sum(target_t1$F.p.uncultured > 0)/length(target_t1$F.p.uncultured) # t1 prevalence too high (>0.75) for the new keystone species method (Amit et al, 2023)
sum(target_t2$F.p.uncultured > 0)/length(target_t2$F.p.uncultured) # t2 prevalence too low (<0.25) for the new keystone species method (Amit et al, 2023)

target_wide <- target %>% 
  pivot_wider(names_from = Timepoint, values_from = F.p.uncultured,
              id_cols = 1)
# By checking target_wide, 17 out of 18 individuals have their F.p.uncultured abundance changed
# before/after antibiotics. So this might be feasible to do longitudinal keystone species method (Amit et al, 2023)

### Adjust the data so it fits Cooccur function
# https://griffithdan.github.io/pages/code_and_data/cooccur.html
# https://medium.com/analytics-vidhya/how-to-create-co-occurrence-networks-with-the-r-packages-cooccur-and-visnetwork-f6e1ceb1c523
rownames(final_merged) <- final_merged$Run
input1 <- final_merged[ , -c(1:5)]
input2 <- as.data.frame(t(input1))

# Turn into presence-absence data
input3 <- ifelse(input2 > 0, yes = 1, no = 0)

# Shorten the species names
original_names <- as.data.frame(str_split_fixed(rownames(input3), pattern = fixed("."), n = 9))
original_names$V1 <- gsub("X_", "", original_names$V1)
original_names$V1 <- gsub("_", " ", original_names$V1)
original_names$V2 <- gsub("_", " ", original_names$V2)
original_names$V3 <- gsub("_", " ", original_names$V3)
original_names$V4 <- gsub("_", " ", original_names$V4)
original_names$V5 <- gsub("_", " ", original_names$V5)
original_names$V6 <- gsub("_", " ", original_names$V6)
rownames(input3) <- paste(original_names$V1, original_names$V2, original_names$V3,
                          original_names$V4, original_names$V5, original_names$V6,
                          sep = " ")

#Check if the clean name preserves all the unique species names!!
n_distinct(rownames(input3)) == n_distinct(rownames(input2)) #Yes!! :)

# Change "F. prausnitzii uncultured" to "uncultured F. prausnitzii"
which(rownames(input3) == "Faecalibacterium prausnitzii  uncultured Faecalibacterium sp  Faecalibacterium")
rownames(input3)[39] <- "uncultured Faecalibacterium prausnitzii"

co_result <- cooccur(input3, spp_names = T, thresh = T)
summary(co_result)

# Heatmap
e = plot(co_result) + 
  theme(axis.text.x = element_text(hjust = 0,
                                   vjust = 0, angle = 0,
                                   face = "italic"))
#ggsave("20240101_Palleja_Suez_Cooccur_species_matrix.pdf", plot = e, device = "pdf", width = 25, height = 25)


obs.v.exp(co_result)
b = effect.sizes(co_result)
c = pair.attributes(co_result)

# Boxplot
#fix(pair.profile) to make x axis labels vertical
d = pair.profile(co_result) + 
  theme(axis.text.x = element_text(hjust = 0,
                                   vjust = 0, angle = -90,
                                   face = "italic"))
#ggsave("2024_Palleja_Suez_Cooccur_species_boxplot.pdf", plot = d, device = "pdf", width = 20, height = 10)
 


