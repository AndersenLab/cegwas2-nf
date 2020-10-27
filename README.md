# Cegwas2-nf branches
### manplot_only branch
- Use this branch if you're screening lots of traits.
- Will only run the main mapping step and write out manhattan plots, pheno-geno plots and the underlying tables.
- Errors will be ignored to let all traits finish. 
- Uses docker container andersenlab/cegwas2:latest.

### report branch
- Use this branch if you only have a few traits.
- Will generate html report containing all plots in interactive formats for each trait and each QTL.  
- If one downloads the "Analysis-Results-data" folder, the Rmarkdown files should knit without any modifications (given all libraries are installed).
- The environment is set up via Conda and module/R.3.6.0 on Quest b/c I couldn't build 1 required package into the docker container.

### master branch
- Left unchanged since 2019 Stefan's commit. In case need to reproduce results.

# Important notes and updates
### `process fix_strain_names_bulk`
- This process by default will take the strain list and convert strain names to the corresponding isotype names. Multiple strains within the same isotype will have their trait averaged. The was needed b/c (1) we used to only have variants for isotypes, and not for strains (2) one should only use 1 strain per isotype for mapping, b/c using strains very similar to each other will lead to artifacts in GWAS.
- In older version: `--fix_names=no` (or anything not “fix”) will skip the step above, and at the same time skip a outlier pruning step in cegwas2::process_phenotype. 
- **In current version, `--fix_names` is no longer used. Instead, the decision will be made based on the input vcf automatically. If the vcf name contains "20180527.imputed" (such as the default one in cegwas2-nf/bin), the old behavior will apply. Otherwise, strain names will not be converted or changed, but pruning will be done. In which case, the input vcf MUST have the same strain names as the trait file, and it's the users' responsiblity to make sure there is only 1 strain per isotype in input trait files.**
- **To implement it in your own version, copy over "bin/Fix_strain_names_bulk_new.R", and change the script in this process to be**
`Rscript --vanilla ${workflow.projectDir}/bin/Fix_Isotype_names_bulk_new.R ${phenotypes} ${params.vcf}` 
[Example here](https://github.com/AndersenLab/cegwas2-nf/commit/b54afbc2d76db20f0744fdbd11634516aa05565f).

### `process vcf_to_geno_matrix`
- `bcftools filter -i N_MISSING=0` will remove all variants that is `./.` (missing, unassigned) in any samples
- Hard filtered vcf has a lot of sites with `./.`, and imputation will fill these sites with `0/0` or `1/1` based on similarity of nearby sites with other strains. After the filtering step above, hard filtered vcf will have only 1/3 of the variants left comparing to imputed vcf. So we have more confidence in variants in hard filtered vcf, but much less markers will be included for mapping leading to less spatial resolution.

### `process plot_genes`
- If the input vcf doesn’t contain annotation by snpeff (which, imputed vcfs don’t, hard filtered vcfs do), cegwas2 will download the very large soft-filtered vcf from Cendr and import into R. This is likely the reason this step takes so much memory. 
- The memory allocation now dynamically depends on how many times the excecution fails. So no manual adjustment of memory for this step should be needed.

### Error in QTL on mtDNA
- Usually if there is a QTL on mtDNA, it will end up with a size of 1bp and throws errors in some processes. 
- Cegwas2 finds the peak marker, and defines interval of QTL by n markers to the left, and n markers to the right from the peak. on main chromosomes, cegwas2 will detect whether the marker is at end of chromosome to make sure QTL doesn't go beyond chromosome boundary. but it is not detecting it for MtDNA. Since MtDNA is so small, the peak +/- n markers usually go beyond the length of full MtDNA, in which case the QTL length will be artificially set to be 1bp.
