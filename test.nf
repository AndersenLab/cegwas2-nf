#! usr/bin/env nextflow

nextflow.preview.dsl=2

/*
~ ~ ~ > * USER INPUT PARAMETERS 
*/
date = new Date().format( 'yyyyMMdd' )

params.traitdir  = null
params.traitfile = null
params.vcf 		 = null
params.p3d 		 = null
params.sthresh   = null
params.freqUpper = 0.05
params.minburden = 2
params.refflat   = "${workflow.projectDir}/bin/refFlat.ws245.txt"
params.genes     = "${workflow.projectDir}/bin/gene_ref_flat.Rda"
params.cendr_v   = "20180527"
params.e_mem 	 = "10"
params.eigen_mem = params.e_mem + " GB"
params.group_qtl = 1000
params.ci_size   = 150
params.fix_names = "fix"
params.help 	 = null
params.R_libpath = "/projects/b1059/software/R_lib_3.6.0"
params.burden    = true
params.finemap   = false

println()

/*
~ ~ ~ > * OUTPUT DIRECTORY 
*/

params.out = "Analysis_Results-${date}"

/*
~ ~ ~ > * INITIATE GENE LOCATION FILE 
*/


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
    log.info "--R_libpath                  String                Path to the required R libraries (DEFAULT = /projects/b1059/software/R_lib_3.4.1)"
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
~ ~ ~ > * WORKFLOW
*/
workflow {

	// VCF
	if(params.vcf) {
		vcf = Channel.fromPath("${params.vcf}")
		vcf_index = Channel.fromPath("${params.vcf}" + ".tbi")

	} else {
		vcf = pull_vcf.out.dl_vcf
		vcf_index = pull_vcf.out.dl_vcf_index
	}

	// Fix strain names
	Channel.fromPath("${params.traitfile}") | fix_strain_names_bulk
	traits_to_map = fix_strain_names_bulk.out.fixed_strain_phenotypes
			.flatten()
    		.map { file -> tuple(file.baseName.replaceAll(/pr_/,""), file) }

    // Genotype matrix
    pheno_strains = fix_strain_names_bulk.out.phenotyped_strains_to_analyze

    vcf.spread(vcf_index)
    		.combine(pheno_strains) | vcf_to_geno_matrix

    // EIGEN
	contigs = Channel.from(["I", "II", "III", "IV", "V", "X"])
	contigs.combine(vcf_to_geno_matrix.out) | chrom_eigen_variants
	chrom_eigen_variants.out.collect() | collect_eigen_variants

	// GWAS mapping
	collect_eigen_variants.out
			.spread(vcf_to_geno_matrix.out)
			.spread(traits_to_map)
			.spread(Channel.from("${params.p3d}"))
			.spread(Channel.from("${params.sthresh}"))
			.spread(Channel.from("${params.group_qtl}"))
			.spread(Channel.from("${params.ci_size}")) | rrblup_maps

	// summarize
	rrblup_maps.out.processed_map_to_summary_plot.collect() | summarize_maps

	// only continue if manplot_only = false
	if(params.finemap) {
		// LD
		rrblup_maps.out.pr_maps_trait
				.combine(vcf_to_geno_matrix.out) | LD_between_regions

		// Fine map
		summarize_maps.out.qtl_peaks
			   .splitCsv(sep: '\t')
			   .join(rrblup_maps.out.processed_map_to_ld)
			   .join(rrblup_maps.out.pr_maps_trait)
			   .spread(vcf.spread(vcf_index))
			   .spread(pheno_strains) | prep_ld_files | rrblup_fine_maps
		rrblup_fine_maps.out.prLD | concatenate_LD_per_trait

		// Plot fine map
		Channel.fromPath("${params.genes}")
				.spread(concatenate_LD_per_trait.out)
				.spread(vcf.spread(vcf_index)) | plot_genes

		// Burden mapping
		if(params.burden) {
			traits_to_map
				.spread(vcf.spread(vcf_index))
				.spread(Channel.fromPath("${params.refflat}")) | burden_mapping | plot_burden
		}
	}


	
}


process pull_vcf {

	tag {"PULLING VCF FROM CeNDR"}
	executor 'local'

	output:
		path "*.vcf.gz", emit: dl_vcf 
		path "*.vcf.gz.tbi", emit: dl_vcf_index 

	"""
		wget https://storage.googleapis.com/elegansvariation.org/releases/${params.cendr_v}/variation/WI.${params.cendr_v}.impute.vcf.gz
		tabix -p vcf WI.${params.cendr_v}.impute.vcf.gz
	"""
}

process fix_strain_names_bulk {

	executor 'local'

	tag {"BULK TRAIT"}

	publishDir "${params.out}/Phenotypes", mode: 'copy'

	input:
		file(phenotypes)

	output:
		path "pr_*.tsv", emit: fixed_strain_phenotypes 
		path "Phenotyped_Strains.txt", emit: phenotyped_strains_to_analyze 

	"""
		# add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
		echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Fix_Isotype_names_bulk.R > Fix_Isotype_names_bulk.R 

		Rscript --vanilla Fix_Isotype_names_bulk.R ${phenotypes} ${params.vcf} ${params.fix_names}
	"""

}



process vcf_to_geno_matrix {

	executor 'local'

	publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

	cpus 1

	input:
		tuple file(vcf), file(index), file(strains)

	output:
		file("Genotype_Matrix.tsv") 

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


/*
============================================================
~ > *                                                  * < ~
~ ~ > *                                              * < ~ ~
~ ~ ~ > *  EIGEN DECOMPOSITION OF GENOTYPE MATRIX  * < ~ ~ ~
~ ~ > *                                              * < ~ ~
~ > *                                                  * < ~
============================================================
*/


/*
------------ Decomposition per chromosome
*/

process chrom_eigen_variants {

	tag { CHROM }

	cpus 6
	memory params.eigen_mem

	input:
		tuple val(CHROM), file(genotypes)


	output:
		file("${CHROM}_independent_snvs.csv")


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
		file(chrom_tests) //from sig_snps_geno_matrix.collect()

	output:
		file("total_independent_tests.txt") //into independent_tests

	"""
		cat *independent_snvs.csv |\\
		grep -v inde |\\
		awk '{s+=\$1}END{print s}' > total_independent_tests.txt
	"""

}


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
	tuple file("independent_snvs.csv"), path(geno), val(TRAIT), path(pheno), val(P3D), val(sig_thresh), val(qtl_grouping_size), val(qtl_ci_size) 

	output:
	tuple TRAIT, path(geno), path(pheno), emit: "processed_map_to_ld"
	path "*processed_mapping.tsv", emit: processed_map_to_summary_plot
	tuple TRAIT, path("*processed_mapping.tsv"), emit: "pr_maps_trait"
	file("*.pdf")

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


	publishDir "${params.out}/Mappings/Plots", mode: 'copy', pattern: "*mappings.pdf"
	publishDir "${params.out}/Mappings/Data", mode: 'copy', pattern: "*_peaks.tsv"

    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }

	input:
	file(maps)

	output:
	file("*.pdf")
	path "QTL_peaks.tsv", emit: qtl_peaks
  	//val true into qtl_peaks_done


	"""
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Summarize_Mappings.R > Summarize_Mappings.R 
		
    Rscript --vanilla Summarize_Mappings.R

		cat  *processed_mapping.tsv |\\
		awk '\$0 !~ "\\tNA\\t" {print}' |\\
		awk '!seen[\$2,\$4,\$5,\$11,\$12,\$13,\$14,\$17]++' |\\
		awk 'NR>1{print \$5, \$2, \$12, \$13, \$14,\$4,\$11,\$17}' OFS="\\t" > QTL_peaks.tsv

		sig_maps=`wc -l QTL_peaks.tsv | cut -f1 -d' '`

		if [ \$sig_maps = 0 ]; then
			max_log10=`cat *processed_mapping.tsv | awk 'BEGIN {max = 0} {if (\$4>max && \$4!= "log10p") max=\$4} END {print max}'`
			echo "NO TRAITS HAD SIGNIFICANT MAPPINGS - MAXIMUM -log10p IS \$max_log10 - CONSIDER SETTING BF THRESHOLD BELOW THIS VALUE"
			exit
		fi
	"""
}

/*
------------ Generate linkage plot between QTL regions
*/

process LD_between_regions {

  tag { TRAIT }

  publishDir "${params.out}/Mappings/Data", mode: 'copy', pattern: "*LD_between_QTL_regions.tsv"

  input:
  tuple val(TRAIT), path("processed_mapping.tsv"), path("Genotype_Matrix.tsv")

  output:
  tuple val(TRAIT), path("*LD_between_QTL_regions.tsv") optional true
  //val(TRAIT) into linkage_done

  """

    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/LD_between_regions.R > LD_between_regions.R 

    Rscript --vanilla LD_between_regions.R Genotype_Matrix.tsv processed_mapping.tsv ${TRAIT}

  """
}

/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  RUN FINE MAPPING  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/

/*
------------ Extract QTL interval genotype matrix
*/

process prep_ld_files {

	tag {TRAIT}

	input:
		tuple val(TRAIT), val(CHROM), val(start_pos), val(peak_pos), val(end_pos), val(log10p), val(var_exp), val(h2), path(complete_geno), path(phenotype), path(pr_map), path(vcf), path(index), path(strains)

	output:
		tuple val(TRAIT), val(CHROM), val(start_pos), val(peak_pos), val(end_pos), val(log10p), val(var_exp), val(h2), path(complete_geno), path(phenotype), path(pr_map), path(vcf), path(index), path("*ROI_Genotype_Matrix.tsv"), path("*LD.tsv") 

	"""
		echo "HELLO"
		cat ${pr_map} |\\
		awk '\$0 !~ "\\tNA\\t" {print}' |\\
		awk '!seen[\$2,\$5,\$12,\$13,\$14]++' |\\
		awk 'NR>1{print \$2, \$5, \$12, \$13, \$14}' OFS="\\t" > ${TRAIT}_QTL_peaks.tsv

		filename='${TRAIT}_QTL_peaks.tsv'
		echo Start
		while read p; do 
		    chromosome=`echo \$p | cut -f1 -d' '`
		    trait=`echo \$p | cut -f2 -d' '`
		    start_pos=`echo \$p | cut -f3 -d' '`
		    peak_pos=`echo \$p | cut -f4 -d' '`
		    end_pos=`echo \$p | cut -f5 -d' '`
		
		cat ${phenotype} | awk '\$0 !~ "strain" {print}' | cut -f1 > phenotyped_samples.txt

		bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${vcf} \
		-S phenotyped_samples.txt |\\
		bcftools filter -i N_MISSING=0 |\\
		awk '\$0 !~ "#" {print \$1":"\$2}' > \$trait.\$chromosome.\$start_pos.\$end_pos.txt


		bcftools view --regions \$chromosome:\$start_pos-\$end_pos ${vcf} \
		-S phenotyped_samples.txt |\\
		bcftools filter -i N_MISSING=0 -Oz -o finemap.vcf.gz

		plink --vcf finemap.vcf.gz \\
			--snps-only \\
			--maf 0.05 \\
			--biallelic-only \\
			--allow-extra-chr \\
			--set-missing-var-ids @:# \\
			--geno \\
			--make-bed \\
			--recode vcf-iid bgz \\
			--extract \$trait.\$chromosome.\$start_pos.\$end_pos.txt \\
			--out \$trait.\$chromosome.\$start_pos.\$end_pos

		nsnps=`wc -l \$trait.\$chromosome.\$start_pos.\$end_pos.txt | cut -f1 -d' '`

		plink --r2 with-freqs \\
			--allow-extra-chr \\
			--snps-only \\
			--ld-window-r2 0 \\
			--ld-snp \$chromosome:\$peak_pos \\
			--ld-window \$nsnps \\
			--ld-window-kb 6000 \\
			--chr \$chromosome \\
			--out \$trait.\$chromosome:\$start_pos-\$end_pos.QTL \\
			--set-missing-var-ids @:# \\
			--vcf \$trait.\$chromosome.\$start_pos.\$end_pos.vcf.gz

		sed 's/  */\\t/g' \$trait.\$chromosome:\$start_pos-\$end_pos.QTL.ld |\\
		cut -f2-10 |\\
		sed 's/^23/X/g' | sed 's/\\t23\\t/\\tX\\t/g' > \$trait.\$chromosome.\$start_pos.\$end_pos.LD.tsv

		bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' finemap.vcf.gz |\\
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
			sed 's/.\\/./NA/g' |\\
			sed 's/^23/X/g' > \$trait.\$chromosome:\$start_pos-\$end_pos.ROI_Genotype_Matrix.tsv

		done < \$filename
	"""
}

/*
------------ Run fine mapping
*/

process rrblup_fine_maps {

  errorStrategy { task.attempt == 3 ? 'ignore' : 'retry' }

	publishDir "${params.out}/Fine_Mappings/Plots", mode: 'copy', pattern: "*_finemap_plot.pdf"
	publishDir "${params.out}/Fine_Mappings/Data", mode: 'copy', pattern: "*_prLD_df.tsv"
	tag {"${TRAIT} ${CHROM}:${start_pos}-${end_pos}"}

	input:
		tuple val(TRAIT), val(CHROM), val(start_pos), val(peak_pos), val(end_pos), val(log10p), val(var_exp), val(h2), path(complete_geno), path(phenotype), path(pr_map), path(vcf), path(index), path(roi_geno_matrix), path(roi_ld) 


	output:
		//tuple path("*pdf"), path("*prLD_df.tsv")
		tuple val(TRAIT), path(phenotype), path(roi_geno_matrix), path("*prLD_df.tsv"), emit: prLD
		tuple path("*pdf"), path("*prLD_df.tsv")

	"""
        for i in *ROI_Genotype_Matrix.tsv;
        do
        	start_pos=`echo \$i | cut -f2 -d':' | cut -f1 -d'.' | cut -f1 -d'-'`
        	end_pos=`echo \$i | cut -f2 -d':' | cut -f1 -d'.' | cut -f2 -d'-'`

        	ld_file=`ls *LD.tsv | grep "\$start_pos" | grep "\$end_pos" | tr -d '\\n'`
        	echo "\$ld_file"

          echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/Finemap_QTL_Intervals.R > Finemap_QTL_Intervals.R
          
          Rscript --vanilla Finemap_QTL_Intervals.R ${complete_geno} \$i ${pr_map} \$ld_file ${task.cpus} ${params.p3d}
        done   

	"""
}

/*
------------ Process fine mapping : concatenate QTL finemapping files for each each trait
*/

process concatenate_LD_per_trait {

	tag {TRAIT}

	executor 'local'

	input:
		tuple val(TRAIT), path(phenotype), path(roi_geno_matrix), path(processed_ld)


	output:
		set val(TRAIT), file("${TRAIT}.combined_peak_LD.tsv"), file(phenotype)

	"""
		for i in *prLD_df.tsv;
		do
		chrom=`head -2 \$i | tail -1 | cut -f2`
		start_pos=`echo \$i | rev | cut -f3 -d"_" | cut -f2 -d'.' | rev`
		end_pos=`echo \$i | rev | cut -f3 -d"_" | cut -f1 -d'.' | rev`

		echo "${TRAIT} \$chrom \$start_pos \$end_pos"

		awk 'BEGIN{OFS="\\t"};FNR == 1 {next};\$1=\$1' OFS="\\t" \$i |\
		awk -v trait="${TRAIT}" \
			-v chrom="\$chrom" \
			-v start_pos="\$start_pos" \
			-v end_pos="\$end_pos" \
			'{print \$0, trait, start_pos, end_pos}' OFS="\t" > ${TRAIT}.\$chrom.\$start_pos.\$end_pos.pr_ld.tsv
		done

		cat *pr_ld.tsv |\
		awk 'BEGIN{OFS="\\t"; print "marker", "CHROM", "POS", "log10p", "peak_marker", "peak_maf", "maf_marker_b", "ld_r2", "trait", "start_pos", "end_pos"}; 
			{print \$0}' > ${TRAIT}.combined_peak_LD.tsv
	"""

}

/*
------------ Plot Fine mapping results overlaid on gene locations
*/

process plot_genes {

	cpus 1
    memory { 32.GB * task.attempt }
    //errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
    errorStrategy 'ignore'

	tag {phenotype}

	publishDir "${params.out}/Fine_Mappings/Plots", mode: 'copy', pattern: "*_gene_plot.pdf"
	publishDir "${params.out}/Fine_Mappings/Data", mode: 'copy', pattern: "*snpeff_genes.tsv"

	input:
		tuple path(genes), val(TRAIT), path(ld), path(phenotype), path(vcf), path(vcfindex)

	output:
		tuple val(TRAIT), file("*snpeff_genes.tsv"), file("*pdf")
    	//val(TRAIT) into gene_plts_done

	"""
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/plot_genes.R > plot_genes.R
    
    Rscript --vanilla plot_genes.R ${ld} ${phenotype} ${genes} ${vcf} ${params.cendr_v}
	"""
}

/*
====================================
~ > *                          * < ~
~ ~ > *                      * < ~ ~
~ ~ ~ > *  BURDEN MAPPING  * < ~ ~ ~
~ ~ > *                      * < ~ ~
~ > *                          * < ~
====================================
*/


process burden_mapping {

	tag {TRAIT}

	publishDir "${params.out}/BURDEN/SKAT/Data", mode: 'copy', pattern: "*.Skat.assoc"
	publishDir "${params.out}/BURDEN/VT/Data", mode: 'copy', pattern: "*.VariableThresholdPrice.assoc"

	input:
		tuple val(TRAIT), file(trait_df), file(vcf), file(index), file(refflat)

	output:
		tuple val(TRAIT), file("*.Skat.assoc"), file("*.VariableThresholdPrice.assoc")

	"""
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/makeped.R > makeped.R
    
    Rscript --vanilla makeped.R ${trait_df}

		n_strains=`wc -l ${trait_df} | cut -f1 -d" "`
		min_af=`bc -l <<< "${params.minburden}/(\$n_strains-1)"`

		rvtest \\
		--pheno ${TRAIT}.ped \\
		--out ${TRAIT} \\
		--inVcf ${vcf} \\
		--freqUpper ${params.freqUpper} \\
		--freqLower \$min_af \\
		--geneFile ${refflat} \\
		--vt price \\
		--kernel skat
	"""
}

process plot_burden {

	executor 'local'

	tag {TRAIT}

	publishDir "${params.out}/BURDEN/SKAT/Plots", mode: 'copy', pattern: "*SKAT.pdf"
	publishDir "${params.out}/BURDEN/VT/Plots", mode: 'copy', pattern: "*VTprice.pdf"

	input:
		tuple val(TRAIT), file(skat), file(vt)

	output:
		tuple val(TRAIT), file("*SKAT.pdf"), file("*VTprice.pdf")

	"""
    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/plot_burden.R > plot_burden.R

    Rscript --vanilla plot_burden.R ${TRAIT} ${skat} ${vt}
	"""
}


/*
=====================================
~ > *                           * < ~
~ ~ > *                       * < ~ ~
~ ~ ~ > *  GENERATE REPORT  * < ~ ~ ~
~ ~ > *                       * < ~ ~
~ > *                           * < ~
=====================================
*/

workflow.onComplete {

    summary = """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

    """

    println summary

    def outlog = new File("${params.out}/log.txt")
    outlog.newWriter().withWriter {
        outlog << param_summary
        outlog << summary
    }

    // mail summary
    if (params.email) {
        ['mail', '-s', 'cegwas2-nf', params.email].execute() << summary
    }


}
