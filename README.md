# cegwas2-nf
GWA mapping with C. elegans


## Overview of the workflow

![alt text](https://github.com/AndersenLab/cegwas2-nf/blob/master/images/Cegwas2_flow_v2.png)

## Required software for running on QUEST

1. [nextflow-v20.0+](https://www.nextflow.io/docs/latest/getstarted.html)

*Users can either update Nextflow to the newest version to run OR load a conda environment for Nextflow v20 using the following commands:*
```
module load python/anaconda3.6
source activate /projects/b1059/software/conda_envs/nf20_env
```

## Required software for running outside of QUEST

*These packages should be in the user's PATH* 

1. [R-v3.6.0](https://www.r-project.org/)
1. [nextflow-v20.0+](https://www.nextflow.io/docs/latest/getstarted.html)
1. [BCFtools-v1.9](https://samtools.github.io/bcftools/bcftools.html)
1. [plink-v1.9](https://www.cog-genomics.org/plink2)
1. [bedtools-2.29.2](https://bedtools.readthedocs.io/en/latest/content/installation.html)
1. [R-cegwas2](https://github.com/AndersenLab/cegwas2)
1. [R-tidyverse-v1.2.1](https://www.tidyverse.org/)
1. [R-coop-0.6-2](https://cran.r-project.org/web/packages/coop/index.html)
1. [R-rrBLUP-v4.6](https://cran.r-project.org/web/packages/rrBLUP/rrBLUP.pdf)
1. [R-plotly-4.9.2](https://cran.r-project.org/web/packages/plotly/index.html)
1. [R-DT-0.12](https://cran.r-project.org/web/packages/DT/index.html)
1. [R-data.table-1.12.8](https://cran.r-project.org/web/packages/data.table/index.html)
1. [R-Rcpp-1.0.1](https://cran.r-project.org/web/packages/Rcpp/index.html)
1. [R-genetics-1.3.8.1.2](https://cran.r-project.org/web/packages/genetics/index.html)
1. [R-sommer-4.0.4](https://cran.r-project.org/web/packages/sommer/index.html)
1. [R-RSpectra-v0.13-1](https://github.com/yixuan/RSpectra)
1. [pandoc=2.12](https://pandoc.org/installing.html)
1. [R-knitr-1.28](https://cran.r-project.org/web/packages/knitr/index.html)
1. [R-rmarkdown-2.1](https://cran.r-project.org/web/packages/rmarkdown/index.html)
1. [R-cowplot-1.0.0](https://cran.r-project.org/web/packages/cowplot/index.html)
1. [R-ggbeeswarm-v0.6](https://github.com/eclarke/ggbeeswarm)

## Required data for running outside of QUEST

1. VCF(s)
   - A hard-filtered vcf containing phenotyped samples for mapping
   - A tabix-generated index hard-filtered vcf (.tbi)
   - An imputed vcf

## Testing pipeline using Nextflow
Running debug mode is a good way to quickly test if your environment is set up correctly. Entire debug run should take 2-3 minutes.
```
git clone https://github.com/AndersenLab/cegwas2-nf.git
cd cegwas2-nf
nextflow main.nf --debug
```
## Execution of pipeline using Nextflow
```
git clone https://github.com/AndersenLab/cegwas2-nf.git
cd cegwas2-nf
nextflow main.nf --traitfile <path to traitfile> --annotation bcsq [optional parameters, see below]
```

### Profiles
Users can select from a number of profiles that each run different processes for the analysis:

* `-profile standard` - This is the default profile if no profile is included. This profile runs the GWA mapping, fine mapping, burden mapping, and generates plots and dataframes in the output folder.

* `-profile manplot_only` - This profile only runs the GWA mapping and generates manhattan and phenotype-by-genotype plots for each trait. No fine mapping or burden mapping is performed. This profile is good for users with hundreds of traits.

* `-profile reports` - This profile runs all the same processes as the standard profile with an additional analysis for haplotypes and divergent regions in the QTL. In addition to the standard outputs, this profile also outputs an html markdown report for each trait with the markdown file that can be run independently on a local computer to reproduce the output files and figures. This profile works best with <30 traits.

### Parameters

* `--traitfile` - is a tab-delimited formatted (.tsv) file that contains trait information.  Each phenotype file should be in the following format (replace trait_name with the phenotype of interest):

| strain | trait_name_1 | trait_name_2 |
| --- | --- | --- |
| JU258 | 32.73 | 19.34 |
| ECA640 | 34.065378 | 12.32 |
| ... | ... | ... | 124.33 |
| ECA250 | 34.096 | 23.1 |

* `--annotation` - Users can choose between "snpeff" annotation or "bcsq" annotation. Currently, annotation is only set up with the 20210121 release.

* `nextflow main.nf --help` - will display the help message


#### Optional parameters
* `--vcf` - CeNDR release date for the VCF file with variant data. Default is "20210121". Hard-filter VCF will be used for the GWA mapping and imputed VCF will be used for fine mapping.

* `--sthresh` - This determines the signficance threshold required for performing post-mapping analysis of a QTL. `BF` corresponds to Bonferroni correction (DEFAULT), `EIGEN` corresponds to correcting for the number of independent markers in your data set, and `user-specified` corresponds to a user-defined threshold, where you replace user-specified with a number. For example `--sthresh=4` will set the threshold to a `-log10(p)` value of 4. We recommend using the strict `BF` correction as a first pass to see what the resulting data looks like. If the pipeline stops at the `summarize_maps` process, no significant QTL were discovered with the input threshold. You might want to consider lowering the threshold if this occurs. 

* `--p3d` - This determines what type of kinship correction to perform prior to mapping. `TRUE` corresponds to the EMMAx method and `FALSE` corresponds to the slower EMMA method. We recommend running with `--p3d=TRUE` to make sure all files of the required files are present and in the proper format, then run with `--p3d=FALSE` (DEFAULT) for a more exact mapping.

* `--freqUpper` - Upper bound for variant allele frequency for a variant to be considered for burden mapping. Default = 0.5

* `--minburden` - The number of strains that must share a variant for that variant to be considered for burden mapping. Default = 2

* `--refflat` - Genomic locations for genes used for burden mapping. A default generated from WS245 is provided in the repositories bin. 

* `--genes` - Genomic locations for genes formatted for plotting purposes. A default generated from WS245 is provided in the repositories bin.

### R scripts

* `Get_GenoMatrix_Eigen.R` - Takes a genotype matrix and chromosome name as input and identifies the number significant eigenvalues.
* `Fix_Isotype_names_bulk.R` - Take sample names present in phenotype data and changes them to isotype names found on [CeNDR](elegansvariation.org).
* `Run_Mappings.R` - Performs GWA mapping using the rrBLUP R package and the EMMA or EMMAx algorithm for kinship correction. Generates manhattan plot and phenotype by genotype plot for peak positions.
* `Summarize_Mappings.R` - Generates plot of all QTL identified in nextflow pipeline.
* `Finemap_QTL_Intervals.R` - Run EMMA/EMMAx on QTL region of interest. Generates fine map plot, colored by LD with peak QTL SNV found from genome-wide scan
* `plot_genes.R` - Runs SnpEff and generates gene plot. 
* `makeped.R` - Converts trait `.tsv` files to `.ped` format for burden mapping.
* `rvtest` - Executable to run burden mapping, can be found at the [RVtests homepage](https://github.com/zhanxw/rvtests)
* `plot_burden.R` - Plots the results from burden mapping.

### Output Folder Structure

```
Analysis_Results-{Date}
  |
  ├──cegwas2_report_traitname_main.html
  ├──cegwas2_report_traitname_main.Rmd
  |
  ├──Phenotypes
      ├── strain_issues.txt
      ├── pr_traitname.tsv
  ├──Genotype_Matrix
      ├── Genotype_Matrix.tsv
      ├── total_independent_tests.txt
  ├──Mappings
      ├── Data             
          ├── traitname_processed_mapping.tsv
          ├── QTL_peaks.tsv
      ├── Plots   
          ├── traitname_manplot.pdf
          ├── traitname_pxgplot.pdf
          ├── Summarized_mappings.pdf
  ├──Fine_Mappings
      ├── Data             
          ├── traitname_snpeff_genes.tsv
          ├── traitname_qtlinterval_prLD_df.tsv
      ├── Plots   
          ├── traitname_qtlinterval_finemap_plot.pdf
          ├── traitname_qtlinterval_gene_plot.pdf
  ├──Divergent_and_haplotype
      ├──all_QTL_bins.bed
      ├──all_QTL_div.bed
      ├──div_isotype_list.txt
      ├──haplotype_in_QTL_region.txt
  ├──BURDEN
      ├── VT             
          ├── Data             
              ├── traitname.VariableThresholdPrice.assoc
          ├── Plots   
              ├── traitname_VTprice.pdf
      ├── SKAT   
          ├── Data             
              ├── traitname.Skat.assoc
          ├── Plots   
              ├── traitname_SKAT.pdf

```

#### Phenotypes folder
* `strain_issues.txt` - Output of any strain names that were changed to match vcf (i.e. isotypes that are not reference strains)
* `pr_traitname.tsv` - Processed phenotype file for each trait. This is the file that goes into the mapping

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


#### BURDEN folder (Contains two subfolders VT/SKAT with the same structure)

##### Data
* `traitname.VariableThresholdPrice.assoc` - Genome-wide burden mapping result using VT price, see [RVtests homepage](https://github.com/zhanxw/rvtests)
* `traitname.Skat.assoc` - Genome-wide burden mapping result using Skat, see [RVtests homepage](https://github.com/zhanxw/rvtests)

##### Plots
* `traitname_VTprice.pdf` - Genome-wide burden mapping manhattan plot for VTprice
* `traitname_SKAT.pdf` - Genome-wide burden mapping manhattan plot for Skat
