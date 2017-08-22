# Load libraries
library("dplyr")
library("ggplot2")

df_irep <- read.table("results_irep.tsv", col.names=c("Sample", "Bin", "iRep"))
df_irep <- data.frame(df_irep, filtered = factor(rep(c("Yes","No"),nrow(df_irep))))

# Format sample names and bin names to more manageable state
df_irep$Sample <- gsub(pattern = "_.*", replacement = "", df_irep$Sample)
df_irep$Sample <- gsub(pattern = ".C", replacement = "", df_irep$Sample, fixed = TRUE) # This is to match with metadata 
# .C is just a indicator of which raw data was used (varying library prep for some samples)
df_irep$Bin <- gsub(".*/", "", df_irep$Bin)
df_irep$Bin <- gsub("-contigs.fa", "", df_irep$Bin)


# Read and add metadata
meta <- read.csv2("metadata.csv")
# Make additional factor for taxonomic classification of the bins
df_classification <- read.table("bin_classification.tsv", header = TRUE)
# meta$Sample_ID <- gsub(meta$Sample_ID, pattern = ".", replacement = "_", fixed = TRUE)
data_total <- left_join(df_irep, meta, by = c("Sample" = "Sample_ID"))
data_total <- left_join(data_total, df_classification, by ="Bin")
data_total$iRep <- as.numeric(as.character(data_total$iRep))

# First select only filtered iRep values (that meet all threshold criteria)
data_total_filtered <- data_total[data_total$filtered == "Yes" & is.na(data_total$iRep) == FALSE, ]
data_total_filtered <- droplevels(data_total_filtered)
data_total_unfiltered <- data_total[data_total$filtered == "No" & is.na(data_total$iRep) == FALSE, ]
data_total_unfiltered <- droplevels(data_total_unfiltered)

# Filter out candidate Limnohabitans genomes
data_total_filtered_limno <- data_total_filtered[data_total_filtered$Classification == "Comamonadaceae", ]

# Lets plot some stuff
p_limno_filtered <- ggplot(data = data_total_filtered_limno, aes(x = Season, y = iRep, fill = Season))+
  geom_point(size = 4, shape = 21)+
  scale_fill_brewer(palette = "Paired")+
  facet_grid(Site~Bin)+
  theme_bw()+
  geom_boxplot(alpha = 0.4)

png("irep_limno_filtered_season.png", width = 12, height = 5, res = 500, units = "in")
p_limno_filtered
dev.off()


# Filter out candidate Limnohabitans genomes
data_total_unfiltered_limno <- data_total_unfiltered[data_total_unfiltered$Classification == "Comamonadaceae", ]

# Lets plot some stuff
p_limno_unfiltered <- ggplot(data = data_total_unfiltered_limno, aes(x = Season, y = iRep, fill = Season))+
  geom_point(size = 4, shape = 21)+
  scale_fill_brewer(palette = "Paired")+
  facet_grid(Site~Bin)+
  theme_bw()+
  geom_boxplot(alpha = 0.4)

png("irep_limno_unfiltered_season.png", width = 15, height = 10, res = 500, units = "in")
p_limno_unfiltered
dev.off()


# Lets see for Nitrosomonas
data_total_nitro <- data_total_filtered[data_total_filtered$Classification == "Nitrosomonas", ]

p_nitro <- ggplot(data = data_total_nitro, aes(x = Site, y = iRep, fill = Season))+
  geom_point(size = 4, shape = 21)+
  scale_fill_brewer(palette = "Paired")+
  facet_grid(Site~Bin)+
  theme_bw()+
  geom_boxplot(alpha = 0.4)
