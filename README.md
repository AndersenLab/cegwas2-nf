# cegwas2-nf
GWA mapping with C. elegans


## Overview of the workflow

![alt text](https://github.com/AndersenLab/cegwas2-nf/blob/master/images/Cegwas2_FLOWCHART.png)


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
nextflow run AndersenLab/cegwas2-nf --traitdir=test_traits --vcf=bin/WI.20180527.impute.vcf.gz --p3d=TRUE --sthresh=BF
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
