# cegwas2-nf
GWA mapping with C. elegans


## Overview of the workflow

![alt text](https://github.com/AndersenLab/cegwas2-nf/blob/master/images/Cegwas2_flow_v2.png)


## Required software packages that should be in users PATH

1. nextflow-0.31.1
1. BCFtools-v1.9
1. plink-v1.9
1. R-cegwas2-Found on GitHub
1. R-tidyverse-v1.2.1
1. R-correlateR-Found on GitHub
1. R-rrBLUP-v4.6
1. R-sommer-v3.5
1. R-RSpectra-v0.13-1

## Execution of pipeline using Nextflow
```
nextflow main.nf --traitdir=test_traits --vcf=bin/WI.20180527.impute.vcf.gz --p3d=TRUE --sthresh=BF
```

### Parameters

* `--traitdir` - is a directory that contains one file for each trait the user wants to map. The file name should correspond to the phenotype name and be in tab-delimited format (.tsv). Each phenotype file should be in the following format (replace trait_name with the phenotype of interest):

| strain | trait_name |
| --- | --- |
| JU258 | 32.73 |
| ECA640 | 34.065378 |
| ... | ... | ... |
| ECA250 | 34.096 |

* `--vcf` - is a VCF file with variant data. All strains with phenotypes should be represented in the VCF used for mapping. There should also abe a tabix-generated index file (.tbi) in the same folder as the specified VCF file that has the same name as the VCF except for the addition of the `.tbi` extension. (generated using `tabix -p vcf vcfname.vcf.gz`)

* `--p3d` - This determines what type of kinship correction to perform prior to mapping. `TRUE` corresponds to the EMMAx method and `FALSE` corresponds to the slower EMMA method. We recommend running with `--p3d=TRUE` to make sure all files of the required files are present and in the proper format, then run with `--p3d=FALSE` for a more exact mapping.

* `--sthresh` - This determines the signficance threshold required for performing post-mapping analysis of a QTL. `BF` corresponds to Bonferroni correction, `EIGEN` corresponds to correcting for the number of independent markers in your data set, and `user-specified` corresponds to a user-defined threshold, where you replace user-specified with a number. For example `--sthresh=4` will set the threshold to a `-log10(p)` value of 4. We recommend using the strict `BF` correction as a first pass to see what the resulting data looks like. If the pipeline stops at the `summarize_maps` process, no significant QTL were discovered with the input threshold. You might want to consider lowering the threshold if this occurs. 

### R scripts

* `Get_GenoMatrix_Eigen.R` - Takes a genotype matrix and chromosome name as input and identifies the number significant eigenvalues.
* `Fix_Isotype_names.R` - Take sample names present in phenotype data and changes them to isotype names found on [CeNDR](elegansvariation.org)
* `Run_Mappings.R` - Performs GWA mapping using the rrBLUP R package and the EMMA or EMMAx algorithm for kinship correction. Generates manhattan plot and phenotype by genotype plot for peak positions.
* `Summarize_Mappings.R` - Generates plot of all QTL identified in nextflow pipeline.
* `Finemap_QTL_Intervals.R` - Run EMMA/EMMAx on QTL region of interest. Generates fine map plot, colored by LD with peak QTL SNV found from genome-wide scan
* `plot_genes.R` - Runs SnpEff and generates gene plot. 

### Output 

#### Genotype_Matrix folder
* `Genotype_Matrix.tsv` - pruned LD-pruned genotype matrix used for GWAS and construction of kinship matrix
* `total_independent_tests.txt` - number of independent tests determined through spectral decomposition of the genotype matrix

#### Mappings folder

##### Data
* `traitname_processed_mapping.tsv` - Processed mapping data frame for each trait mapped
* `QTL_peaks.tsv` - List of signifcant QTL identified across all traits

##### Plots
* `traitname_manplot.pdf` - Manhattan plot for each trait that was analyzed. Two significance threshold lines are present, one for the Bonferronit corrected threshold, and another for the spectral decomposition threshold.
* `traitname_pxgplot.pdf` - Phenotype by genotype split at peak QTL positions for every significant QTL identified
* `Summarized_mappings.pdf` - A summary plot of all QTL identified

#### Fine_Mappings folder

##### Data
* `traitname_snpeff_genes.tsv` - Fine-mapping data frame for all significant QTL

##### Plots
* `traitname_qtlinterval_finemap_plot.pdf` - Fine map plot of QTL interval, colored by marker LD with the peak QTL identified from the genome-wide scan
* `traitname_qtlinterval_gene_plot.pdf` - variant annotation plot overlaid with gene CDS for QTL interval




