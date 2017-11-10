# Limno_lakeMI
Repository for the analysis of reconstructed Limnohabitans genomes from lake Michigan

# Genome-centered metagenomics

## 1. QC and assembly
Adapter + Quality trimming (qc.sh)

Digital normalization (BBnorm.sh – downsampling over high coverage areas) + custom scripts to remove unpaired reads and interleave (max. 60x, min. 5x)

IDBA-UD for creating individual assemblies:

Non-stringent taxonomic classification using custom wrapper of `ClassifyContigNR.py` from the [DESMAN](https://github.com/chrisquince/DESMAN)
software. 

Extracting target taxonomy and format into compatible output for `Vizbin`.

## 2. Binning strategy
* Rename contigs with anvi'o script, filtering out contigs shorter than 3kb for initial binning.  

* Initial binning with `Vizbin` (default settings) overlaid by the taxonomy of interest (i.e. Beta-proteobacteria). Extract all Beta-proteobacteria bins.  

* Example vizbin plot for sample `Fa13.BD.MLB.DN`  


## 3. Binning QC and refinement




## 4. Calculation of relative abundances



