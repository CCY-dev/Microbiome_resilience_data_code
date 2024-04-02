rm(list = ls())
setwd("~/Documents/Lab/Multiple_longitudinal/Microbe_resilience/")
library(tidyverse)
library(readxl)
library(ggpubr)
library(patchwork)
library(paletteer)
library(ggtext)

########### Load Suez study first ###########
# Read in master table
suez_master <-read.table(file = "../Autologous_FMT/Data/Suez_stool_species_master.txt",
                    header = T, sep = "\t", stringsAsFactors = F)
suez_species_values <- suez_master[ , -c(1:8)]


########### Then load Palleja study ###########
# Read in master table
pal_master <- read.table("../Antibiotic_Palleja/My_ngless_result/Palleja_stool_ngelss_motus2.5_Species_rarefy190_master.txt", header = T, sep = "\t")

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

# Remove those zero-abundance ones
merged_paired_val_nonzero <- merged_paired_val[ , colSums(merged_paired_val) > 0]

# Remove those low-abundant ones
# Those less than 5 non-zero values are removed
unname(sort(colSums(merged_paired_val_nonzero > 0)))
merged_paired_val_nonzero_new <- merged_paired_val_nonzero[colSums(merged_paired_val_nonzero > 0) > 4]

# Test the correlation of all species v.s. all species
# To see if Faecalibacterium has higher correlation than others
library(reshape2)
library(Hmisc)
cormat <- rcorr(as.matrix(merged_paired_val_nonzero_new), type="spearman")

melted_cormat_r <- melt(cormat$r)
melted_cormat_p <- melt(cormat$P)

melted_cormat <- inner_join(melted_cormat_r, melted_cormat_p,
                            by = c("Var1", "Var2"))
colnames(melted_cormat) <- c("source", "target", "edge", "p_value")

melted_cormat <- melted_cormat %>%
  filter(source != target) %>%  # Exclude self-correlation
  mutate(q_value = p.adjust(p_value, "fdr"),
         abs_rho = abs(edge),
         dir = sign(edge))

# Filter the correlations
melted_cormat_clean <- melted_cormat %>%
  filter(q_value < 0.1 & abs_rho > 0) #q < 0.1 and rho != 0

### Clean species name. This is for plotting prettier cytoscape network for manuscript
melted_cormat_clean_edited <- melted_cormat_clean
melted_cormat_clean_edited$source <- gsub(pattern = "_", x = melted_cormat_clean_edited$source, replacement = "")
melted_cormat_clean_edited$target <- gsub(pattern = "_", x = melted_cormat_clean_edited$target, replacement = "")
melted_cormat_clean_edited$source <- gsub(pattern = "XEubacterium", x = melted_cormat_clean_edited$source, replacement = "Eubacterium")
melted_cormat_clean_edited$target <- gsub(pattern = "XEubacterium", x = melted_cormat_clean_edited$target, replacement = "Eubacterium")
melted_cormat_clean_edited$source <- gsub(pattern = "XRuminococcus", x = melted_cormat_clean_edited$source, replacement = "Ruminococcus")
melted_cormat_clean_edited$target <- gsub(pattern = "XRuminococcus", x = melted_cormat_clean_edited$target, replacement = "Ruminococcus")

## Edit source name
source_name <- str_split_fixed(melted_cormat_clean_edited$source, pattern = fixed("."), n = 4)

# Edit V4
source_name[400, 2] <- "sp Eubacterium sp CAG 76" # This is to avoid same name of different strain
source_name[source_name[ , 4] == "Oscillibacter.sp..ER4Firmicutes.bacterium.CAG.1295924", 2] <- "sp ER4" # This is to avoid same name of different strain
source_name[!source_name[ , 4] %in% c("sedis") & 
              !str_detect(string = source_name[ , 4], pattern = "CAG") &
              !source_name[ , 3] %in% c("CAG"), 4] <- ""
source_name[nchar(source_name[ , 4]) > 7 , 4] <- ""

# Edit V3
source_name[!source_name[ , 3] %in% c("incertae", "sp", "uncultured", "CAG", "Candidatus") , 3] <- ""

## Edit target name
target_name <- str_split_fixed(melted_cormat_clean_edited$target, pattern = fixed("."), n = 4)
# Edit V4
target_name[400, 2] <- "sp Eubacterium sp CAG 76" # This is to avoid same name of different strain
target_name[target_name[ , 4] == "Oscillibacter.sp..ER4Firmicutes.bacterium.CAG.1295924", 2] <- "sp ER4" # This is to avoid same name of different strain
target_name[!target_name[ , 4] %in% c("sedis") & 
              !str_detect(string = target_name[ , 4], pattern = "CAG") &
              !target_name[ , 3] %in% c("CAG"), 4] <- ""
target_name[nchar(target_name[ , 4]) > 7 , 4] <- ""
# Edit V3
target_name[!target_name[ , 3] %in% c("incertae", "sp", "uncultured", "CAG", "Candidatus") , 3] <- ""

melted_cormat_clean_edited$source <- paste(source_name[ , 1], source_name[ , 2], source_name[ , 3], source_name[ , 4], sep = " ")
melted_cormat_clean_edited$target <- paste(target_name[ , 1], target_name[ , 2], target_name[ , 3], target_name[ , 4], sep = " ")

#Check if the clean name preserves all the unique species names
n_distinct(melted_cormat_clean$source) == n_distinct(melted_cormat_clean_edited$source) 
n_distinct(melted_cormat_clean$target) == n_distinct(melted_cormat_clean_edited$target)

### Replace F. prausnitzii uncultured with uncultured F. prausnitzii
melted_cormat_clean_edited$source <- gsub(x = melted_cormat_clean_edited$source, pattern = "Faecalibacterium prausnitzii uncultured",
                                          replacement = "uncultured Faecalibacterium prausnitzii")
melted_cormat_clean_edited$target <- gsub(x = melted_cormat_clean_edited$target, pattern = "Faecalibacterium prausnitzii uncultured",
                                          replacement = "uncultured Faecalibacterium prausnitzii")

#write.table(melted_cormat_clean_edited, quote = F, sep = "\t",
#            "species_all_cor_long_clean_name.txt", row.names = F)

##### Replace the microbe in really small nodes with numbers
# The two vectors with repeated items
vector1 <- melted_cormat_clean_edited$source
vector2 <- melted_cormat_clean_edited$target

# Combine both vectors
combined_vector <- c(vector1, vector2)

# Use factor with shared levels to assign unique numbers
unique_numbers <- as.integer(factor(combined_vector, levels = unique(combined_vector)))

# Skip those microbes with large node
large_node <- c("uncultured Faecalibacterium prausnitzii ",
                "Ruminococcaceae species incertae sedis",
                "Fusicatenibacter saccharivorans  ",
                "uncultured Clostridium sp ",
                "Bifidobacterium adolescentis  ",
                "Sutterella species incertae sedis")

# Set items in large_node to NA in the combined vector
combined_vector[combined_vector %in% large_node] <- NA

# Use factor with shared levels to assign unique numbers
unique_numbers <- as.integer(factor(combined_vector, levels = unique(na.omit(combined_vector)))) #neglect the NA in level

# Assign the unique numbers back to each vector
result_vector1 <- unique_numbers[1:length(vector1)]
result_vector2 <- unique_numbers[(length(vector1) + 1):length(combined_vector)]

# Replace NA values with original names for items in large_node
result_vector1[is.na(result_vector1)] <- vector1[is.na(result_vector1)]
result_vector2[is.na(result_vector2)] <- vector2[is.na(result_vector2)]

melted_cormat_clean_edited_numbered <- melted_cormat_clean_edited
melted_cormat_clean_edited_numbered$source <- result_vector1
melted_cormat_clean_edited_numbered$target <- result_vector2

# Write a mapping table for output
mapping_tb <- data.frame(Number = c(melted_cormat_clean_edited_numbered$source, melted_cormat_clean_edited_numbered$target),
                         Microbe = c(melted_cormat_clean_edited$source, melted_cormat_clean_edited$target))
mapping_tb <- mapping_tb %>% 
  dplyr::distinct() %>% 
  filter(!Microbe %in% large_node)

#write.table(melted_cormat_clean_edited_numbered, quote = F, sep = "\t",
#            "20240101_Species_all_cor_long_clean_name_numbered.txt", row.names = F)

#write.table(mapping_tb, quote = F, sep = "\t",
#            "20240101_Species_all_cor_long_clean_name_numbered_mapping_tbl.txt", row.names = F)


##### And then import species_all_cor_long_clean_name.txt into Cytoscape

###### Plotting scattered plot: And then read in the centrality results from cytoscape
cyto_result <- read.table(file = "20240101_Species_all_cor_long_network_analyze_cleanName.csv",
                          header = T, sep = ",", stringsAsFactors = F)

cyto_result <- cyto_result %>%
  dplyr::select(name, everything()) %>%
  mutate(Txt_color = ifelse(str_detect(name, "Faecalibacterium"), "black", "black"))

# Edit the species names
high_beweenness <- cyto_result$name[cyto_result$BetweennessCentrality %in% sort(cyto_result$BetweennessCentrality, decreasing = T)[1:5]] #top 5
high_node <- cyto_result$name[cyto_result$Degree %in% sort(cyto_result$Degree, decreasing = T)[1:4]] #top3
label_for_plotting <- c(high_beweenness, high_node)
label_for_plotting <- label_for_plotting[-which(label_for_plotting == "2")]

cyto_result$name[cyto_result$name %in% label_for_plotting] <- c("uncultured *Faecalibacterium prausnitzii*", "Ruminococcaceae species incertae sedis",  
                                                                "*Fusicatenibacter saccharivorans*", "uncultured *Clostridium* sp",              
                                                                "*Bifidobacterium adolescentis*")

# Below is necessary, not duplicated code.
high_beweenness <- cyto_result$name[cyto_result$BetweennessCentrality %in% sort(cyto_result$BetweennessCentrality, decreasing = T)[1:5]] #top 5
high_node <- cyto_result$name[cyto_result$Degree %in% sort(cyto_result$Degree, decreasing = T)[1:4]] #top3
label_for_plotting <- c(high_beweenness, high_node)
label_for_plotting <- label_for_plotting[-which(label_for_plotting == "2")]

set.seed(100)
ggplot(cyto_result, aes(Degree, BetweennessCentrality), position_jitter(height = 0.001, width = 0.1)) +
  geom_jitter(alpha = 0.6, color = "cornflowerblue", size = 2.5) +
  theme_light() +
  ggtext::geom_richtext(data=subset(cyto_result, name %in% label_for_plotting),
                        size = 4, aes(label = name), 
                        fill = NA, label.color = NA, #remove label border
                        label.padding = grid::unit(rep(0, 4), "pt"), # remove padding
                        vjust="inward",hjust="inward") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  labs(x = "Node degree",
       y = "Betweenness centrality") +
  scale_color_manual(values = c("red"= "red", "black" = "black")) +
  guides(color="none")
#ggsave("20240101_Species_betweenness_centrality_allblack_cleanName.pdf", width = 8, height = 4, device = "pdf")



