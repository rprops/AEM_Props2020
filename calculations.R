# Load libraries
library(tidyr)
library(dplyr)
library(ggplot2)

# Read data
mean_coverage <- read.table("mean_coverage_selected_final.tsv", header = TRUE)
std_coverage <- read.table("std_coverage_selected_final.tsv", header = TRUE)


# From wide to long format
mean_coverage_long <- gather(mean_coverage, Sample_ID, coverage, 
                             Fa13_BD_MLB_DN:Su13_BD_MM15_SN_C, factor_key=TRUE)
mean_coverage_long[,2] <- gsub(mean_coverage_long[,2], pattern = "_C", 
                               replacement = "")

std_coverage_long <- gather(std_coverage, Sample_ID, std_coverage, 
                            Fa13_BD_MLB_DN:Su13_BD_MM15_SN_C, 
                            factor_key=TRUE)
std_coverage_long[,2] <- gsub(std_coverage_long[,2], pattern = "_C", 
                            replacement = "")
 
coverage_data <- data.frame(mean_coverage_long, 
                            std_coverage = std_coverage_long[,3])

# Read and add metadata
meta <- read.csv2("metadata.csv")
meta$Sample_ID <- gsub(meta$Sample_ID, pattern = ".", replacement = "_", fixed = TRUE)
data_total <- left_join(coverage_data, meta, by = "Sample_ID")

# Plot abundance distributions of all bins

ggplot(data = data_total, aes(x = Season, y = log2(coverage), fill = bins))+
  geom_point(size = 4, shape = 21)+
  scale_fill_brewer(palette = "Paired")+
  facet_grid(.~Site)+
  theme_bw()+
  geom_errorbar(aes(ymin=log2(coverage-std_coverage), 
                    ymax=log2(coverage+std_coverage)), 
                    colour="black", width=.1)

ggplot(data = data_total, aes(x = Season, y = coverage, fill = bins))+
  geom_point(size = 4, shape = 21)+
  scale_fill_brewer(palette = "Paired")+
  facet_grid(.~Site)+
  theme_bw()+
  geom_errorbar(aes(ymin=coverage-std_coverage, 
                    ymax=coverage+std_coverage), 
                colour="black", width=.1)

# Add completeness estimates, GC content and genome size estimates to metadata


#### Relative abundance
# Read data
rel_abun <- read.table("abundance_selected_final.tsv", header = TRUE)


# From wide to long format
rel_abun_long <- gather(rel_abun, Sample_ID, abundance, 
                             Fa13_BD_MLB_DN:Su13_BD_MM15_SN_C, factor_key=TRUE)
rel_abun_long[,2] <- gsub(rel_abun_long[,2], pattern = "_C", 
                               replacement = "")



# Read and add metadata
meta <- read.csv2("metadata.csv")
meta$Sample_ID <- gsub(meta$Sample_ID, pattern = ".", replacement = "_", fixed = TRUE)
data_rel_abund <- left_join(rel_abun_long, meta, by = "Sample_ID")

# Plot abundance distributions of all bins

ggplot(data = data_rel_abund, aes(x = Season, y = log2(abundance), fill = Site))+
  geom_point(size = 4, shape = 21)+
  scale_fill_brewer(palette = "Paired")+
  facet_grid(.~bins)+
  theme_bw()

png("abundance.png", width = 20, height = 5, res = 500, units = "in")
ggplot(data = data_rel_abund, aes(x = Season, y = abundance, fill = Site))+
  geom_point(size = 4, shape = 21)+
  scale_fill_brewer(palette = "Paired")+
  facet_grid(.~bins)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


# Add completeness estimates, GC content and genome size estimates to metadata
"general_bins_summary_selected_final.tsv"
