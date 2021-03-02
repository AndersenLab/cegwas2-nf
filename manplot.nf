#! usr/bin/env nextflow


/*
~ ~ ~ > * USER INPUT PARAMETERS 
*/
date = new Date().format( 'yyyyMMdd' )

params.traitdir  = null
params.traitfile = null
params.vcf     = null
params.p3d     = null
params.sthresh   = null
params.freqUpper = 0.05
params.minburden = 2
params.refflat   = "${workflow.projectDir}/bin/refFlat.ws245.txt"
params.genes     = "${workflow.projectDir}/bin/gene_ref_flat.Rda"
params.cendr_v   = "20180527"
params.e_mem   = "32"
params.eigen_mem = params.e_mem + " GB"
params.group_qtl = 1000
params.ci_size   = 150
params.help    = null
params.fix_names = "fix"
params.max_QTL_number = 5
params.max_QTL_size = 2000000
params.R_libpath="/projects/b1059/software/R_lib_3.6.0"



println()

/*
~ ~ ~ > * OUTPUT DIRECTORY 
*/

params.out = "Analysis_Results-${date}"

/*
~ ~ ~ > * INITIATE GENE LOCATION FILE 
*/

genes = Channel.fromPath("${params.genes}")

log.info ""
log.info "------------------------------------------"
log.info "        C. elegans GWAS pipeline "
log.info "------------------------------------------"
log.info ""


if (params.help) {
    log.info "----------------------------------------------------------------"
    log.info "                      USAGE                                     "
    log.info "----------------------------------------------------------------"
    log.info ""
    log.info "nextflow main.nf --traitfile=test_bulk --vcf=bin/WI.20180527.impute.vcf.gz --p3d=TRUE --sthresh=EIGEN # run all traits from a single file"
    log.info "nextflow main.nf --traitdir=test_bulk --p3d=TRUE --sthresh=BF # download VCF from CeNDR"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--traitfile              String                Name of file that contains phenotypes. File should be tab-delimited with the columns: strain trait1 trait2 ..."
    log.info "--vcf                    String                Name of VCF to extract variants from. There should also be a tabix-generated index file with the same name in the directory that contains the VCF. If none is provided, the pipeline will download the latest VCF from CeNDR"
    log.info "--p3d                    BOOLEAN               Set to FALSE for EMMA algortith, TRUE for EMMAx"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
    log.info "Optional arguments (General):"
    log.info "--out                    String                Name of folder that will contain the results"
    log.info "--e_mem                  String                Value that corresponds to the amount of memory to allocate for eigen decomposition of chromosomes (DEFAULT = 100)"
    log.info "--cendr_v                String                CeNDR release (DEFAULT = 20180527)"
    log.info "Optional arguments (Marker):"
    log.info "--sthresh                String                Significance threshold for QTL - Options: BF - for bonferroni correction, EIGEN - for SNV eigen value correction, or another number e.g. 4"
    log.info "--group_qtl              Integer               If two QTL are less than this distance from each other, combine the QTL into one, (DEFAULT = 1000)"
    log.info "--ci_size                Integer               Number of SNVs to the left and right of the peak marker used to define the QTL confidence interval, (DEFAULT = 150)"
    log.info "Optional arguments (Burden):"
    log.info "--freqUpper              Float                 Maximum allele frequency for a variant to be considered for burden mapping, (DEFAULT = 0.05)"
    log.info "--minburden              Interger              Minimum number of strains to have a variant for the variant to be considered for burden mapping, (DEFAULT = 2)"
    log.info "--genes                  String                refFlat file format that contains start and stop genomic coordinates for genes of interest, (DEFAULT = bin/gene_ref_flat.Rda)"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Information describing the stucture of the input files can be located in input_files/README.txt"
    log.info ""
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info ""
    log.info " Required software packages to be in users path"
    log.info "BCFtools               v1.9"
    log.info "plink                  v1.9"
    log.info "R-cegwas2              Found on GitHub"
    log.info "R-tidyverse            v1.2.1"
    log.info "R-correlateR           Found on GitHub"
    log.info "R-rrBLUP               v4.6"
    log.info "R-sommer               v3.5"
    log.info "R-RSpectra             v0.13-1"
    log.info "R-ggbeeswarm           v0.6.0"
    log.info "--------------------------------------------------------"    
    exit 1
} else {

log.info ""
log.info "Phenotype Directory                     = ${params.traitdir}"
log.info "VCF                                     = ${params.vcf}"
log.info "CeNDR Release                           = ${params.cendr_v}"
log.info "P3D                                     = ${params.p3d}"
log.info "Significance Threshold                  = ${params.sthresh}"
log.info "Max AF for Burden Mapping               = ${params.freqUpper}"
log.info "Min Strains with Variant for Burden     = ${params.minburden}"
log.info "Significance Threshold                  = ${params.sthresh}"
log.info "Gene File                               = ${params.genes}"
log.info "Result Directory                        = ${params.out}"
log.info "Eigen Memory allocation                 = ${params.eigen_mem}"
log.info "Path to R libraries.                    = ${params.R_libpath}"
log.info ""
}

/*
~ ~ ~ > * COMBINE VCF AND VCF INDEX INTO A CHANNEL
*/

if (params.vcf) {
  
  vcf = Channel.fromPath("${params.vcf}")

  vcf_index = Channel.fromPath("${params.vcf}" + ".tbi")

  vcf
    .spread(vcf_index)
    .into{vcf_to_whole_genome;
        vcf_to_fine_map;
        vcf_to_burden;
        vcf_to_query_vcf}

} else {

  process pull_vcf {

    tag {"PULLING VCF FROM CeNDR"}
    executor 'local'

    output:
      file("*.vcf.gz") into dl_vcf
      file("*.vcf.gz.tbi") into dl_vcf_index

    """
      wget https://storage.googleapis.com/elegansvariation.org/releases/${params.cendr_v}/variation/WI.${params.cendr_v}.impute.vcf.gz
      tabix -p vcf WI.${params.cendr_v}.impute.vcf.gz
    """
  }

  dl_vcf
    .spread(dl_vcf_index)
    .into{vcf_to_whole_genome;
        vcf_to_fine_map;
        vcf_to_burden;
        vcf_to_query_vcf}

}

/*
~ ~ ~ > * INITIATE MAPPING QTL GROUPING PARAMETER
*/

Channel
  .from("${params.ci_size}")
  .set{qtl_ci_size}

/*
~ ~ ~ > * INITIATE MAPPING QTL CONFIDENCE INTERVAL SIZE PARAMETER
*/

Channel
  .from("${params.group_qtl}")
  .set{qtl_snv_groupinng}

/*
~ ~ ~ > * INITIATE MAPPING METHOD CHANNEL
*/

Channel
  .from("${params.p3d}")
  .into{p3d_full;
      p3d_fine}

/*
~ ~ ~ > * INITIATE THRESHOLD CHANNEL
*/

Channel
  .from("${params.sthresh}")
  .into{sig_threshold_full;
      sig_threshold_fine}

/*
~ ~ ~ > * INITIATE PHENOTYPE CHANNEL - GENERATES A [trait_name, trait_file] TUPLE
*/



Channel
  .fromPath("${params.traitfile}")
  .set{ traits_to_strainlist }

process fix_strain_names_bulk {

  executor 'local'

  tag {"BULK TRAIT"}

  publishDir "${params.out}", mode: 'copy'

  input:
    file(phenotypes) from traits_to_strainlist

  output:
    file("pr_*.tsv") into fixed_strain_phenotypes
    file("Phenotyped_Strains.txt") into phenotyped_strains_to_analyze

  """
    # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Fix_Isotype_names_bulk.R > Fix_Isotype_names_bulk.R 

    Rscript --vanilla Fix_Isotype_names_bulk.R ${phenotypes} ${params.vcf} ${params.fix_names}
  """

}

fixed_strain_phenotypes
    .flatten()
    .map { file -> tuple(file.baseName.replaceAll(/pr_/,""), file) }
  .into{ traits_to_map;
      traits_to_burden }


phenotyped_strains_to_analyze
  .into{strain_list_genome;
      strain_list_finemap}

process vcf_to_geno_matrix {

//  executor 'local'

  publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

  cpus 1

  input:
    set file(vcf), file(index) from vcf_to_whole_genome
    file(strains) from strain_list_genome

  output:
    file("Genotype_Matrix.tsv") into geno_matrix

  """
    bcftools view -S ${strains} ${vcf} |\\
    bcftools filter -i N_MISSING=0 -Oz -o Phenotyped_Strain_VCF.vcf.gz

    tabix -p vcf Phenotyped_Strain_VCF.vcf.gz

    plink --vcf Phenotyped_Strain_VCF.vcf.gz \\
      --snps-only \\
      --biallelic-only \\
      --maf 0.05 \\
      --set-missing-var-ids @:# \\
      --indep-pairwise 50 10 0.8 \\
      --geno \\
      --allow-extra-chr

    awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
    sort -k1,1d -k2,2n > markers.txt

    bcftools query -l Phenotyped_Strain_VCF.vcf.gz |\\
    sort > sorted_samples.txt 

    bcftools view -v snps \\
    -S sorted_samples.txt \\
    -R markers.txt \\
    Phenotyped_Strain_VCF.vcf.gz |\\
    bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' |\\
    sed 's/[[# 0-9]*]//g' |\\
    sed 's/:GT//g' |\\
    sed 's/0|0/-1/g' |\\
    sed 's/1|1/1/g' |\\
    sed 's/0|1/NA/g' |\\
    sed 's/1|0/NA/g' |\\
    sed 's/.|./NA/g'  |\\
    sed 's/0\\/0/-1/g' |\\
    sed 's/1\\/1/1/g'  |\\
    sed 's/0\\/1/NA/g' |\\
    sed 's/1\\/0/NA/g' |\\
    sed 's/.\\/./NA/g' > Genotype_Matrix.tsv
  """

}

geno_matrix
  .into{eigen_gm;
      mapping_gm}


/*
============================================================
~ > *                                                  * < ~
~ ~ > *                                              * < ~ ~
~ ~ ~ > *  EIGEN DECOMPOSITION OF GENOTYPE MATRIX  * < ~ ~ ~
~ ~ > *                                              * < ~ ~
~ > *                                                  * < ~
============================================================
*/

CONTIG_LIST = ["I", "II", "III", "IV", "V", "X"]
contigs = Channel.from(CONTIG_LIST)

/*
------------ Decomposition per chromosome
*/

process chrom_eigen_variants {

  tag { CHROM }

  cpus 1
  memory params.eigen_mem

  input:
    file(genotypes) from eigen_gm
    each CHROM from contigs

  output:
    file("${CHROM}_independent_snvs.csv") into sig_snps_geno_matrix
    file(genotypes) into concat_geno


  """
    cat Genotype_Matrix.tsv |\\
    awk -v chrom="${CHROM}" '{if(\$1 == "CHROM" || \$1 == chrom) print}' > ${CHROM}_gm.tsv

    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Get_GenoMatrix_Eigen.R > Get_GenoMatrix_Eigen.R
    Rscript --vanilla Get_GenoMatrix_Eigen.R ${CHROM}_gm.tsv ${CHROM}  
  """

}

/*
------------ Sum independent tests for all chromosomes
*/

process collect_eigen_variants {

  executor 'local'

  publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

  cpus 1

  input:
    file(chrom_tests) from sig_snps_geno_matrix.collect()

  output:
    file("total_independent_tests.txt") into independent_tests

  """
    cat *independent_snvs.csv |\\
    grep -v inde |\\
    awk '{s+=\$1}END{print s}' > total_independent_tests.txt
  """

}

independent_tests
  .spread(mapping_gm)   // Genotype_Matrix.tsv
  .spread(traits_to_map)   // set val(TRAIT), file(trait_PC1.tsv) single trait file
  .spread(p3d_full)  // p3d value
  .spread(sig_threshold_full)
  .spread(qtl_snv_groupinng)
  .spread(qtl_ci_size)
  .set{mapping_data}

/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  RUN GWAS MAPPING  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/

/*
------------ Genome-wide scan
*/

process rrblup_maps {

  cpus 4

  tag { TRAIT }

  publishDir "${params.out}/Mappings/Data", mode: 'copy', pattern: "*processed_mapping.tsv"
  publishDir "${params.out}/Mappings/Plots", mode: 'copy', pattern: "*.pdf"
  

  input:
  set file("independent_snvs.csv"), file(geno), val(TRAIT), file(pheno), val(P3D), val(sig_thresh), val(qtl_grouping_size), val(qtl_ci_size) from mapping_data

  output:
  file("*raw_mapping.tsv") into raw_map
  set val(TRAIT), file(geno), file(pheno) into processed_map_to_ld
  file("*processed_mapping.tsv") into processed_map_to_summary_plot
  set val(TRAIT), file("*processed_mapping.tsv") into pr_maps_trait
  file("*.pdf") into gwas_plots

  """
  tests=`cat independent_snvs.csv | grep -v inde`

  echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Run_Mappings.R > Run_Mappings.R

    Rscript --vanilla Run_Mappings.R ${geno} ${pheno} ${task.cpus} ${P3D} \$tests ${sig_thresh} ${qtl_grouping_size} ${qtl_ci_size}

    if [ -e Rplots.pdf ]; then
        rm Rplots.pdf
    fi
  """
}

/*
------------ Generate GWAS QTL summary plot
*/

// need to run this first to find significant traits, 
// but then you lose track of it when joining channels below this process and therefore need to re run same step to find all peaks
process summarize_maps {

  publishDir "${params.out}/Mappings/Data", mode: 'copy', pattern: "*_peaks.tsv"

  input:
  file(maps) from processed_map_to_summary_plot.collect()

  output:
  file("QTL_peaks.tsv") into qtl_peaks


  """
    cat *processed_mapping.tsv |\\
    awk '\$0 !~ "\\tNA\\t" {print}' |\\
    awk '!seen[\$2,\$5,\$12,\$13,\$14]++' |\\
    awk 'NR>1{print \$5, \$2, \$12, \$13, \$14}' OFS="\\t" > QTL_peaks.tsv
    sig_maps=`wc -l QTL_peaks.tsv | cut -f1 -d' '`
    if [ \$sig_maps = 0 ]; then   # this is for all traits combined
      max_log10=`cat *processed_mapping.tsv | awk 'BEGIN {max = 0} {if (\$4>max && \$4!= "log10p") max=\$4} END {print max}'`
      echo "NO TRAITS HAD SIGNIFICANT MAPPINGS - MAXIMUM -log10p IS \$max_log10 - CONSIDER SETTING BF THRESHOLD BELOW THIS VALUE"
      exit
    fi
  """
}