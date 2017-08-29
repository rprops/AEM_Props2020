# Abundance calculations
Ruben Props  
`r format(Sys.time(), '%d %B, %Y')`  





```r
# Read data
mean_coverage <- read.table("./anvio_output/mean_coverage_selected_final.tsv", header = TRUE)
std_coverage <- read.table("./anvio_output/std_coverage_selected_final.tsv", header = TRUE)
bin_size <- read.table("./anvio_output/general_bins_summary_selected_final.tsv", header = TRUE)[, c(2,4)]
total_reads <- read.table("./anvio_output/sample_reads.tsv", header = TRUE)
read_length <- 150

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
data_total <- left_join(coverage_data, total_reads, by = c("Sample_ID" = "Bin"))
data_total <- left_join(data_total, bin_size, by = "bins")
data_total <- left_join(data_total, meta, by =  "Sample_ID")

# Calculate relative abundance of the bins
data_total$mean_rel_abundance <- 100*(data_total$coverage*data_total$bin_size)/(read_length*data_total$Total_reads)
data_total$upper_rel_abundance <- 100*((data_total$coverage+data_total$std_coverage)*data_total$bin_size)/(read_length*data_total$Total_reads)
data_total$lower_rel_abundance <- 100*((data_total$coverage-data_total$std_coverage)*data_total$bin_size)/(read_length*data_total$Total_reads)
```

# 1. Phylogenetic tree
![RAxML tree for Limnohabitans genomes](./Tree/PhyloTree_Limno.png)  

# 2. Investigate MAG- and 16S-based abundances
Formula used to calculate relative abundances:
$$Relative\ abundance =100*(\frac{mean\ coverage * bin\ size}{read\ length*total\ sample\ reads })$$


```r
# Plot abundance distributions of all bins

p_season <- ggplot(data = data_total, aes(x = bins, y = mean_rel_abundance, fill = bins))+
  geom_point(size = 4, shape = 21, alpha = 0.7)+
  geom_boxplot(alpha = 0.3)+
  scale_fill_brewer(palette = "Paired")+
  facet_grid(.~Site)+
  theme_bw()+
  geom_errorbar(aes(ymin=lower_rel_abundance, 
                    ymax=upper_rel_abundance, 
                    width=.1))+
  facet_grid(Season~.)+
  # ylim(0,1)+ 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y=element_text(size=18))+
  ylab("Mean relative abundance (%)")

p_station <- ggplot(data = data_total, aes(x = bins, y = mean_rel_abundance, fill = bins))+
  geom_point(size = 4, shape = 21, alpha = 0.7)+
  geom_boxplot(alpha = 0.3)+
  scale_fill_brewer(palette = "Paired")+
  facet_grid(.~Site)+
  theme_bw()+
  geom_errorbar(aes(ymin=lower_rel_abundance, 
                    ymax=upper_rel_abundance, 
                    width=.1))+
  facet_grid(Site~.)+
  # ylim(0,1)+ 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y=element_text(size=18))+
  ylab("Mean relative abundance (%)")

grid.arrange(p_season, p_station, nrow = 2)
```

<img src="Figures/cached/plot abundances-1.png" style="display: block; margin: auto;" />

# 3. Investigate sequence characteristics within coding DNA sequences (CDS)

```r
# First we need the files that map the gene ID to the sequence ID (linux cmd: https://github.com/rprops/MetaG_lakeMI/wiki/11.-Genome-annotation)
# These are stored in the IMG_annotation data for each genome bin

# Next, extract the %GC of each gene from the gff file
Bin_2737471681 <- extract_gc_from_gff("./IMG_annotation/IMG_2737471681/IMG_Data/133052.assembled.gff", 
                                      outputFolder = "GC_analysis")
Bin_2737471682 <- extract_gc_from_gff("./IMG_annotation/IMG_2737471682/IMG_Data/133053.assembled.gff", 
                                      outputFolder = "GC_analysis")
Bin_2737471683 <- extract_gc_from_gff("./IMG_annotation/IMG_2737471683/IMG_Data/133054.assembled.gff", 
                                      outputFolder = "GC_analysis")
Bin_2737471793 <- extract_gc_from_gff("./IMG_annotation/IMG_2737471793/IMG_Data/133647.assembled.gff", 
                                      outputFolder = "GC_analysis")
Bin_2737471794 <- extract_gc_from_gff("./IMG_annotation/IMG_2737471794/IMG_Data/133648.assembled.gff", 
                                      outputFolder = "GC_analysis")
Bin_2737471795 <- extract_gc_from_gff("./IMG_annotation/IMG_2737471795/IMG_Data/133649.assembled.gff", 
                                      outputFolder = "GC_analysis")
Bin_2737471797 <- extract_gc_from_gff("./IMG_annotation/IMG_2737471797/IMG_Data/133651.assembled.gff", 
                                      outputFolder = "GC_analysis")
Bin_2737471799 <- extract_gc_from_gff("./IMG_annotation/IMG_2737471799/IMG_Data/133653.assembled.gff", 
                                      outputFolder = "GC_analysis")
Bin_2737471802 <- extract_gc_from_gff("./IMG_annotation/IMG_2737471802/IMG_Data/133656.assembled.gff", 
                                      outputFolder = "GC_analysis")
Bin_2737471804 <- extract_gc_from_gff("./IMG_annotation/IMG_2737471804/IMG_Data/133658.assembled.gff", 
                                      outputFolder = "GC_analysis")
Bin_2737471805 <- extract_gc_from_gff("./IMG_annotation/IMG_2737471805/IMG_Data/133659.assembled.gff", 
                                      outputFolder = "GC_analysis")
Bin_2737471806 <- extract_gc_from_gff("./IMG_annotation/IMG_2737471806/IMG_Data/133660.assembled.gff", 
                                      outputFolder = "GC_analysis")


# Use these files to make dataframes mapping function (COGs) and %GC
Bin_2737471681_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_133052.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2737471681/IMG_Data/2737471681/2737471681.gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2737471681/IMG_Data/2737471681/2737471681.cog.tab.txt", 
                             gc_thresh = 0.1, output = FALSE, label = "Bin_2737471681")
```

```
## Mon Aug 28 18:39:29 2017  --- There are 2248 genes with > 0.1 %
## Mon Aug 28 18:39:29 2017  --- This is 100 % of all genes
## Mon Aug 28 18:39:29 2017  --- The 10 genes with the highest GC% are:
##                                                                                   function_id
## 2239 NADPH-dependent 2,4-dienoyl-CoA reductase, sulfur reductase, or a related oxidoreductase
## 2240                    tRNA A37 threonylcarbamoyladenosine synthetase subunit TsaC/SUA5/YrdC
## 2241                                 2-C-methyl-D-erythritol 4-phosphate cytidylyltransferase
## 2242                                                                Mg-chelatase subunit ChlD
## 2243                                                  Enoyl-CoA hydratase/carnithine racemase
## 2244                                            Sugar or nucleoside kinase, ribokinase family
## 2245                                       DNA-binding transcriptional regulator, MarR family
## 2246                    NAD(P)H-hydrate repair enzyme Nnr, NAD(P)H-hydrate dehydratase domain
## 2247                                                                              Peptidase E
## 2248                                                               Sulfite exporter TauE/SafE
##      function_name   GC
## 2239           415 67.5
## 2240           211 67.6
## 2241           230 67.6
## 2242           261 68.0
## 2243           257 68.2
## 2244           311 68.2
## 2245           126 68.3
## 2246           284 68.9
## 2247           224 69.3
## 2248           232 69.4
```

```r
Bin_2737471682_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_133053.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2737471682/IMG_Data/2737471682/2737471682.gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2737471682/IMG_Data/2737471682/2737471682.cog.tab.txt", 
                             gc_thresh = 0.1, output = FALSE, label = "Bin_2737471682")
```

```
## Mon Aug 28 18:39:29 2017  --- There are 1882 genes with > 0.1 %
## Mon Aug 28 18:39:29 2017  --- This is 100 % of all genes
## Mon Aug 28 18:39:29 2017  --- The 10 genes with the highest GC% are:
##                                                                function_id
## 1873                                   D-alanyl-D-alanine carboxypeptidase
## 1874                        Aspartate/methionine/tyrosine aminotransferase
## 1875                    DNA-binding transcriptional regulator, LysR family
## 1876                   Predicted esterase of the alpha/beta hydrolase fold
## 1877                                            Sulfite exporter TauE/SafE
## 1878                              Endonuclease YncB, thermonuclease family
## 1879    Predicted house-cleaning NTP pyrophosphatase, Maf/HAM1 superfamily
## 1880                                           L-aminopeptidase/D-esterase
## 1881 tRNA A37 threonylcarbamoyladenosine synthetase subunit TsaC/SUA5/YrdC
## 1882                                                    Glycerate-2-kinase
##      function_name   GC
## 1873           470 68.0
## 1874           393 68.1
## 1875           297 68.1
## 1876           181 68.2
## 1877           232 68.3
## 1878           192 68.5
## 1879           193 68.9
## 1880           348 69.0
## 1881           211 69.5
## 1882           422 73.0
```

```r
Bin_2737471683_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_133054.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2737471683/IMG_Data/2737471683/2737471683.gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2737471683/IMG_Data/2737471683/2737471683.cog.tab.txt", 
                             gc_thresh = 0.1, output = FALSE, label = "Bin_2737471683")
```

```
## Mon Aug 28 18:39:29 2017  --- There are 1727 genes with > 0.1 %
## Mon Aug 28 18:39:29 2017  --- This is 100 % of all genes
## Mon Aug 28 18:39:29 2017  --- The 10 genes with the highest GC% are:
##                                                                               function_id
## 1718                                            Type II secretory pathway, component PulF
## 1719                                                    Pyrroline-5-carboxylate reductase
## 1720                                                                     Shikimate kinase
## 1721                          Phosphoribosylaminoimidazole carboxylase (NCAIR synthetase)
## 1722                                                  ADP-heptose:LPS heptosyltransferase
## 1723                                             Cobalamin biosynthesis protein CobD/CbiB
## 1724 Acetolactate synthase large subunit or other thiamine pyrophosphate-requiring enzyme
## 1725                                                               Uracil-DNA glycosylase
## 1726             NAD(P)-dependent dehydrogenase, short-chain alcohol dehydrogenase family
## 1727                tRNA A37 threonylcarbamoyladenosine synthetase subunit TsaC/SUA5/YrdC
##      function_name   GC
## 1718           397 63.1
## 1719           266 63.1
## 1720           172 63.2
## 1721           375 63.7
## 1722           334 63.8
## 1723           320 63.9
## 1724           550 64.0
## 1725           202 64.2
## 1726           251 64.4
## 1727           211 64.6
```

```r
Bin_2737471793_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_133647.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2737471793/IMG_Data/2737471793/2737471793.gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2737471793/IMG_Data/2737471793/2737471793.cog.tab.txt", 
                             gc_thresh = 0.1, output = FALSE, label = "Bin_2737471793")
```

```
## Mon Aug 28 18:39:29 2017  --- There are 1476 genes with > 0.1 %
## Mon Aug 28 18:39:29 2017  --- This is 100 % of all genes
## Mon Aug 28 18:39:29 2017  --- The 10 genes with the highest GC% are:
##                                                                               function_id
## 1467 Acetolactate synthase large subunit or other thiamine pyrophosphate-requiring enzyme
## 1468                                            2-keto-3-deoxy-L-rhamnonate aldolase RhmA
## 1469            ABC-type dipeptide/oligopeptide/nickel transport system, ATPase component
## 1470   Predicted Rossmann fold nucleotide-binding protein DprA/Smf involved in DNA uptake
## 1471        NAD(P)H-dependent flavin oxidoreductase YrpB, nitropropane dioxygenase family
## 1472                                           Pimeloyl-ACP methyl ester carboxylesterase
## 1473                                    CO or xanthine dehydrogenase, FAD-binding subunit
## 1474                        tRNA A37 threonylcarbamoyladenosine modification protein TsaB
## 1475                                                            Mg-chelatase subunit ChlD
## 1476                                                               Uracil-DNA glycosylase
##      function_name   GC
## 1467           550 67.3
## 1468           255 67.4
## 1469           316 67.5
## 1470           350 67.7
## 1471           336 68.0
## 1472           282 68.0
## 1473           284 68.1
## 1474           220 68.3
## 1475           261 68.3
## 1476           202 69.2
```

```r
Bin_2737471794_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_133648.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2737471794/IMG_Data/2737471794/2737471794.gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2737471794/IMG_Data/2737471794/2737471794.cog.tab.txt", 
                             gc_thresh = 0.1, output = FALSE, label = "Bin_2737471794")
```

```
## Mon Aug 28 18:39:29 2017  --- There are 1226 genes with > 0.1 %
## Mon Aug 28 18:39:29 2017  --- This is 100 % of all genes
## Mon Aug 28 18:39:29 2017  --- The 10 genes with the highest GC% are:
##                                                                             function_id
## 1217                                     Cytochrome c-type biogenesis protein CcmH/NrfF
## 1218                                          Type II secretory pathway, component PulK
## 1219                                           Predicted ATPase with chaperone activity
## 1220                                              Geranylgeranyl pyrophosphate synthase
## 1221                                       Threonine/homoserine efflux transporter RhtA
## 1222                                                         Predicted secreted protein
## 1223                                                            Heme exporter protein D
## 1224                                          UDP-N-acetylmuramyl pentapeptide synthase
## 1225 Predicted Rossmann fold nucleotide-binding protein DprA/Smf involved in DNA uptake
## 1226                                  CO or xanthine dehydrogenase, FAD-binding subunit
##      function_name   GC
## 1217           153 66.0
## 1218           323 66.0
## 1219           490 66.1
## 1220           322 66.2
## 1221           292 66.3
## 1222           148 66.4
## 1223            67 66.7
## 1224           451 66.9
## 1225           350 67.9
## 1226           284 68.0
```

```r
Bin_2737471795_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_133649.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2737471795/IMG_Data/2737471795/2737471795.gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2737471795/IMG_Data/2737471795/2737471795.cog.tab.txt", 
                             gc_thresh = 0.1, output = FALSE, label = "Bin_2737471795")
```

```
## Mon Aug 28 18:39:29 2017  --- There are 1383 genes with > 0.1 %
## Mon Aug 28 18:39:29 2017  --- This is 100 % of all genes
## Mon Aug 28 18:39:29 2017  --- The 10 genes with the highest GC% are:
##                                                                             function_id
## 1374                                                 Glycerol-3-phosphate dehydrogenase
## 1375                                          Predicted amidophosphoribosyltransferases
## 1376 2-polyprenyl-6-methoxyphenol hydroxylase and related FAD-dependent oxidoreductases
## 1377                                 DNA-binding transcriptional regulator, LysR family
## 1378         Tat protein secretion system quality control protein TatD (DNase activity)
## 1379                                              Phosphatidylglycerophosphate synthase
## 1380                                     Methylase of polypeptide chain release factors
## 1381                                                   Methionyl-tRNA formyltransferase
## 1382                Tripartite-type tricarboxylate transporter, receptor component TctC
## 1383                                                                    Cytochrome c553
##      function_name   GC
## 1374           532 66.9
## 1375           225 67.0
## 1376           387 67.0
## 1377           297 67.1
## 1378           256 67.1
## 1379           192 67.5
## 1380           280 67.9
## 1381           307 68.0
## 1382           319 68.9
## 1383           121 69.4
```

```r
Bin_2737471797_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_133651.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2737471797/IMG_Data/2737471797/2737471797.gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2737471797/IMG_Data/2737471797/2737471797.cog.tab.txt", 
                             gc_thresh = 0.1, output = FALSE, label = "Bin_2737471797")
```

```
## Mon Aug 28 18:39:29 2017  --- There are 1076 genes with > 0.1 %
## Mon Aug 28 18:39:29 2017  --- This is 100 % of all genes
## Mon Aug 28 18:39:29 2017  --- The 10 genes with the highest GC% are:
##                                                                             function_id
## 1067              DNA-directed RNA polymerase specialized sigma subunit, sigma24 family
## 1068                                7,8-dihydro-6-hydroxymethylpterin-pyrophosphokinase
## 1069                                           tRNA A37 N6-isopentenylltransferase MiaA
## 1070                                                   FAD/FMN-containing dehydrogenase
## 1071                                            Phosphoribosyl-ATP pyrophosphohydrolase
## 1072                                                ADP-heptose:LPS heptosyltransferase
## 1073                                     UDP-N-acetylmuramoylalanine-D-glutamate ligase
## 1074         Uncharacterized conserved protein YcfJ, contains glycine zipper 2TM domain
## 1075 2-polyprenyl-6-methoxyphenol hydroxylase and related FAD-dependent oxidoreductases
## 1076                                                    Copper oxidase (laccase) domain
##      function_name   GC
## 1067           182 58.9
## 1068           160 58.9
## 1069           308 58.9
## 1070           459 59.1
## 1071            92 59.1
## 1072           334 59.5
## 1073           448 59.7
## 1074           179 59.7
## 1075           387 59.7
## 1076           249 60.0
```

```r
Bin_2737471799_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_133653.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2737471799/IMG_Data/2737471799/2737471799.gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2737471799/IMG_Data/2737471799/2737471799.cog.tab.txt", 
                             gc_thresh = 0.1, output = FALSE, label = "Bin_2737471799")
```

```
## Mon Aug 28 18:39:29 2017  --- There are 522 genes with > 0.1 %
## Mon Aug 28 18:39:29 2017  --- This is 100 % of all genes
## Mon Aug 28 18:39:29 2017  --- The 10 genes with the highest GC% are:
##                                                               function_id
## 513                                             Glutathione S-transferase
## 514                                  DNA-nicking endonuclease, Smr domain
## 515 Carbonic anhydrase or acetyltransferase, isoleucine patch superfamily
## 516                                     Pyrroline-5-carboxylate reductase
## 517                                      16S rRNA U1498 N3-methylase RsmE
## 518                                         Thiamine monophosphate kinase
## 519                      Nicotinamide mononucleotide (NMN) deamidase PncC
## 520                         Pyrimidine reductase, riboflavin biosynthesis
## 521                                                           Glucokinase
## 522                                Nucleoside-diphosphate-sugar epimerase
##     function_name   GC
## 513           211 67.0
## 514           184 67.1
## 515           176 67.2
## 516           266 67.3
## 517           246 67.4
## 518           317 67.7
## 519           162 67.7
## 520           218 67.8
## 521           320 68.0
## 522           314 69.7
```

```r
Bin_2737471802_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_133656.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2737471802/IMG_Data/2737471802/2737471802.gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2737471802/IMG_Data/2737471802/2737471802.cog.tab.txt", 
                             gc_thresh = 0.1, output = FALSE, label = "Bin_2737471802")
```

```
## Mon Aug 28 18:39:29 2017  --- There are 1691 genes with > 0.1 %
## Mon Aug 28 18:39:29 2017  --- This is 100 % of all genes
## Mon Aug 28 18:39:29 2017  --- The 10 genes with the highest GC% are:
##                                                                 function_id
## 1682                          Pyrimidine reductase, riboflavin biosynthesis
## 1683               2-C-methyl-D-erythritol 4-phosphate cytidylyltransferase
## 1684                      CO or xanthine dehydrogenase, FAD-binding subunit
## 1685    Pyrimidine deaminase domain of riboflavin biosynthesis protein RibD
## 1686                                      Phosphohistidine phosphatase SixA
## 1687                                                       Cyanate permease
## 1688                                              Mg-chelatase subunit ChlD
## 1689                                        Chromate transport protein ChrA
## 1690 Uncharacterized conserved protein YcbX, contains MOSC and Fe-S domains
## 1691                               Predicted ATPase with chaperone activity
##      function_name   GC
## 1682           218 67.8
## 1683           230 68.0
## 1684           284 68.2
## 1685           146 68.3
## 1686           163 68.3
## 1687           395 69.1
## 1688           261 69.2
## 1689           195 69.4
## 1690           270 69.6
## 1691           490 70.7
```

```r
Bin_2737471804_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_133658.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2737471804/IMG_Data/2737471804/2737471804.gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2737471804/IMG_Data/2737471804/2737471804.cog.tab.txt", 
                             gc_thresh = 0.1, output = FALSE, label = "Bin_2737471804")
```

```
## Mon Aug 28 18:39:29 2017  --- There are 749 genes with > 0.1 %
## Mon Aug 28 18:39:29 2017  --- This is 100 % of all genes
## Mon Aug 28 18:39:29 2017  --- The 10 genes with the highest GC% are:
##                                                                               function_id
## 740 Uncharacterized conserved protein, contains von Willebrand factor type A (vWA) domain
## 741                   Uncharacterized membrane protein AbrB, regulator of aidB expression
## 742                                Hydroxymethylpyrimidine/phosphomethylpyrimidine kinase
## 743                                                         Urease accessory protein UreH
## 744                                                   Rhodanese-related sulfurtransferase
## 745                                    DNA-binding transcriptional regulator, MarR family
## 746                                                 Geranylgeranyl pyrophosphate synthase
## 747                   ABC-type phosphate/phosphonate transport system, permease component
## 748                     ABC-type phosphate/phosphonate transport system, ATPase component
## 749                   ABC-type phosphate/phosphonate transport system, permease component
##     function_name   GC
## 740           395 64.2
## 741           352 64.3
## 742           263 64.7
## 743           269 65.1
## 744           110 65.1
## 745           126 65.7
## 746           322 66.3
## 747           283 66.3
## 748           258 67.7
## 749           283 67.9
```

```r
Bin_2737471805_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_133659.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2737471805/IMG_Data/2737471805/2737471805.gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2737471805/IMG_Data/2737471805/2737471805.cog.tab.txt", 
                             gc_thresh = 0.1, output = FALSE, label = "Bin_2737471805")
```

```
## Mon Aug 28 18:39:30 2017  --- There are 1172 genes with > 0.1 %
## Mon Aug 28 18:39:30 2017  --- This is 100 % of all genes
## Mon Aug 28 18:39:30 2017  --- The 10 genes with the highest GC% are:
##                                                                   function_id
## 1163           Ribulose-5-phosphate 4-epimerase/Fuculose-1-phosphate aldolase
## 1164    Asp-tRNAAsn/Glu-tRNAGln amidotransferase A subunit or related amidase
## 1165                                   Biotin-(acetyl-CoA carboxylase) ligase
## 1166                       Imidazole glycerol phosphate synthase subunit HisF
## 1167                        ABC-type sugar transport system, ATPase component
## 1168                           Trk K+ transport system, NAD-binding component
## 1169                                                Mg-chelatase subunit ChlI
## 1170 NAD(P)-dependent dehydrogenase, short-chain alcohol dehydrogenase family
## 1171                         Transcriptional regulator of heat shock response
## 1172                                                 Phosphoheptose isomerase
##      function_name   GC
## 1163           219 58.2
## 1164           475 58.2
## 1165           238 58.2
## 1166           256 58.3
## 1167           338 58.4
## 1168           225 58.5
## 1169           423 58.8
## 1170           251 58.9
## 1171           346 59.0
## 1172           176 59.2
```

```r
Bin_2737471806_gc_cog <- gc2function(seq_id_gc = "GC_analysis/seqid_GC_133660.assembled.gff.tsv", gene_id_seq_id ="./IMG_annotation/IMG_2737471806/IMG_Data/2737471806/2737471806.gene_oid_2_seq_id.txt", 
                             functions = "./IMG_annotation/IMG_2737471806/IMG_Data/2737471806/2737471806.cog.tab.txt", 
                             gc_thresh = 0.1, output = FALSE, label = "Bin_2737471806")
```

```
## Mon Aug 28 18:39:30 2017  --- There are 1584 genes with > 0.1 %
## Mon Aug 28 18:39:30 2017  --- This is 100 % of all genes
## Mon Aug 28 18:39:30 2017  --- The 10 genes with the highest GC% are:
##                                                                       function_id
## 1575                                   Predicted NAD/FAD-dependent oxidoreductase
## 1576                                                 Gamma-glutamyltranspeptidase
## 1577                                                       Acetyl esterase/lipase
## 1578                            CO or xanthine dehydrogenase, FAD-binding subunit
## 1579 3-hydroxyisobutyrate dehydrogenase or related beta-hydroxyacid dehydrogenase
## 1580                                                       Uracil-DNA glycosylase
## 1581        Asp-tRNAAsn/Glu-tRNAGln amidotransferase A subunit or related amidase
## 1582                                                   Predicted secreted protein
## 1583    ABC-type dipeptide/oligopeptide/nickel transport system, ATPase component
## 1584                tRNA A37 threonylcarbamoyladenosine modification protein TsaB
##      function_name   GC
## 1575           331 67.0
## 1576           539 67.0
## 1577           312 67.2
## 1578           284 67.3
## 1579           286 67.3
## 1580           202 67.4
## 1581           475 67.5
## 1582           148 67.8
## 1583           316 67.8
## 1584           220 69.8
```

```r
merged_gc_cog <- rbind(Bin_2737471681_gc_cog, Bin_2737471682_gc_cog, Bin_2737471683_gc_cog,
                       Bin_2737471793_gc_cog, Bin_2737471794_gc_cog, Bin_2737471795_gc_cog,
                       Bin_2737471797_gc_cog, Bin_2737471799_gc_cog, Bin_2737471802_gc_cog,
                       Bin_2737471804_gc_cog, Bin_2737471805_gc_cog, Bin_2737471806_gc_cog)
merged_gc_cog$genome_id <- as.character(merged_gc_cog$genome_id)
```

# 4. Analysis of gene length distribution

```r
p_cog_length <- easyGgplot2::ggplot2.histogram(data = merged_gc_cog, xName = 'gene_length',
                  groupName = 'genome_id', alpha = 0.5,
                  legendPosition = "top", binwidth = 0.15, addMeanLine=TRUE, meanLineColor="black",
                  meanLineType="dashed")+ theme_bw()+ ylim(0,15)+
  facet_grid(genome_id~.)+
  labs(x = "Gene length (bp)", y = "Count")+ theme(legend.position="none", axis.text.x = element_text(angle = 45))+
  ggtitle("Limnohabitans MAGs")+ xlim(0,2000)

print(p_cog_length)
```

<img src="Figures/cached/gene length analysis-1.png" style="display: block; margin: auto;" />

# 5. COG functional categories
Get COG ID to COG functional category mapping file here: ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/COG0303/cogs.csv    

The exact statistical analysis to compare genomes based on these profiles should be performed in STAMP.


```r
# Import COG mapping file
cogid_2_cogcat <- read.csv("./Mapping_files/cogid_2_cogcat.csv", sep = ",", header = FALSE, fill = TRUE,col.names = c("COG_ID", "COG_class", "function"))[, 1:2]
cogid_2_cogcat <- cogid_2_cogcat[(cogid_2_cogcat$COG_class)!="", ]
cogid_2_cogcat <- droplevels(cogid_2_cogcat)

# Read COG category file
cog_categories <- read.table("./Mapping_files/cog_categories.tsv", header = TRUE, sep = "\t")

# Merge COG metadata
cog_meta <- dplyr::left_join(cog_categories, cogid_2_cogcat, by = c("COG_class" = "COG_class"))
cog_meta <- droplevels(cog_meta)

# Merge this metadata with the genome data from before
# COGs with multiple classifications are currently still NA - work on this.
merged_gc_cog <- dplyr::left_join(merged_gc_cog, cog_meta, by = c("cog_id" = "COG_ID"))
merged_gc_cog <- merged_gc_cog[!is.na(merged_gc_cog$COG_functional_category),]

# Visualize distribution across major metabolism functional COG groups per genome.
p_cog_func_group <- ggplot(data = merged_gc_cog, aes(x=COG_functional_category, fill = COG_functional_cluster))+
  geom_bar(stat="count", width=0.7, color = "black", size = 0.75)+
  theme_bw()+
  facet_grid(genome_id~.)+
  scale_fill_brewer(palette = "Accent")+
  labs(x = "Gene length (bp)", y = "Count")+ theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1),
                                                   legend.text = element_text(size = 7))+
  ggtitle("Limnohabitans MAGs")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

print(p_cog_func_group)
```

<img src="Figures/cached/COG functional categories-1.png" style="display: block; margin: auto;" />

```r
p_cog_func_clust <- ggplot(data = merged_gc_cog, aes(x=COG_functional_cluster, fill = COG_functional_cluster))+
  geom_bar(stat="count", width=0.7, color = "black", size = 0.75)+
  theme_bw()+
  facet_grid(genome_id~.)+
  scale_fill_brewer(palette = "Accent")+
  labs(x = "Gene length (bp)", y = "Count")+ theme(legend.position="bottom",axis.text.x = element_text(angle = 90, hjust = 1),
                                                   legend.text = element_text(size = 7))+
  ggtitle("Limnohabitans MAGs")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

print(p_cog_func_clust)
```

<img src="Figures/cached/COG functional categories-2.png" style="display: block; margin: auto;" />
