# cegwas2-nf report branch
- Will generate html report containing all plots in interactive formats for each trait and each QTL.  
- If one downloads the "Analysis-Results-data" folder, the Rmarkdown files should knit without any modifications (given all libraries are installed).
- The environment is set up via Conda and Module/R.3.6.0 on Quest.

# Important updates
### `process fix_strain_names_bulk`
- This process by default will take the strain list and convert strain names to the corresponding isotype names. Multiple strains within the same isotype will have their trait averaged. The was needed b/c (1) we used to only have variants for isotypes, and not for strains (2) one should only use 1 strain per isotype for mapping, b/c using strains very similar to each other will lead to artifacts in GWAS.
- In older version: `--fix_names=no` (or anything not “fix”) will skip the step above, and at the same time skip a outlier pruning step in cegwas2::process_phenotype. 
- **In current version, `--fix_names` is no longer used. Instead, the decision will be made based on the input vcf automatically. If the vcf name contains "20180527.imputed" (such as the default one in cegwas2-nf/bin), the old behavior will apply. Otherwise, strain names will not be converted or changed, but pruning will be done. In which case, the input vcf MUST have the same strain names as the trait file.**

### `process vcf_to_geno_matrix`
- `bcftools filter -i N_MISSING=0` will remove all variants that is `./.` (missing, unassigned) in any samples
- Hard filtered vcf has a lot of sites with `./.`, and imputation will fill these sites with `0/0` or `1/1` based on similarity of nearby sites with other strains. After the filtering step above, hard filtered vcf will have only 1/3 of the variants left comparing to imputed vcf. So we have more confidence in variants in hard filtered vcf, but much less markers will be included for mapping leading to less spatial resolution.

### `process plot_genes`
- If the input vcf doesn’t contain annotation by snpeff (which, imputed vcfs don’t, hard filtered vcfs do), cegwas2 will download the very large soft-filtered vcf from Cendr and import into R. This is likely the reason this step takes so much memory. 
- The memory allocation now dynamically depends on how many times the excecution fails. So no manual adjustment of memory for this step should be needed.
