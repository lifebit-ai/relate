#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/relate
========================================================================================
lifebit-ai/relate Analysis Pipeline.
#### Homepage / Documentation
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    // TODO : Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run  lifebit-ai/relate --input .. -profile docker

    Mandatory arguments:
      --input [file]                  File with list of full paths to bcf files and their indexes. Bcf files can be compressed but in a readable for bcftools format.
                                      The name of the files must be consistent across files.
                                      see example:
                                      test_all_chunks_merged_norm_chr10_53607810_55447336.bcf.gz
                                      {name}_{CHR}_{START_POS}_{END_POS}.bcf.gz
                                      Consistency is important here as a variable ('region')
                                      is extracted from the filename.

      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, docker, singularity, test, awsbatch, <institute> and more

      --inputDir                      Input dir for annotation ".txt" files for each of the regions (or chunks ), comming from the metrics compoment and first aggregate annotation of SiteQC pipeline.
      --inputMichiganLDfileExclude    File with regions to be filtered out for improving quality of sites selected.
      --inputPCsancestryrelated       File with Principal Components information comming from reference resources of GEL for the inferred ancestries from the 30k dataset.
      --inputFinalPlatekeys           File with a list of platekeys to be included in the analysis , required for "create_final_king_vcf" process.
      --inputProbs200K                File with Ancestry assignments from GEL's reference resources.
      --inputUNRELATED_1KGP3          File required for the infer_ancestry process , is a two column tab separated file with platekeys on each column.
      --input1KGP3                    File required for the infer_ancestry process , is a three column tab separated file with Sample(Platekey), Family ID and Population asignments.
      --inputSuper_pop_codes          File required for the infer_ancestry process, its tab separated file with population and super population relations and related information.
      --input05both1K100K_eigenvec    File required for the infer_ancestry process, its tab separated file with eigenvectors.
      --inputGELprojection_proj_eigenvec    File required for the infer_ancestry process, its tab separated file with eigenvectors and it has "NA" values in the last column.
      --inputAncestryAssignmentProbs  File required for hwe_pruning_30k_snps process containing tab separated values for probabilities of assignments for the 31 populations code for each platekey(sample)
      --n_pca                          Number of Principal Components desired for gcta process


    Options:


    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size [str]  Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Define variables
n_pca = params.n_pca
awk_expr_create_final_king_vcf_1 = params.awk_expr_create_final_king_vcf_1

// Input list .csv file of tissues to analyse
// [chr10_52955340_55447336, test_all_chunks_merged_norm_chr10_52955340_55447336.bcf.gz, test_all_chunks_merged_norm_chr10_52955340_55447336.bcf.gz.csi]

// Define channels based on params
  Channel.fromPath(params.inputDir+'/*.txt')
                        .ifEmpty { exit 1, "Input dir for annotation txt files not found at ${params.inputDir}. Is the dir path correct?" }
                        .filter{txt -> txt =~/chr\d+/}
                        .map { txt -> ['chr'+ txt.simpleName.split('_chr').last() , txt] }
                        // Filter out chunks from chrX, chrY and chrM - they should not be analysed in Ancestry and Relatedness pipeline
                        .filter { it[0] =~ /chr[^XYM]/ }
                        .set { ch_bcftools_site_metrics_subcols }

  Channel.fromPath(params.inputFinalPlatekeys)
                        .ifEmpty { exit 1, "Input file with samples and platekeys data not found at ${params.inputFinalPlatekeys}. Is the file path correct?" }
                        .set { ch_inputFinalPlatekeys }
  Channel.fromPath(params.inputUNRELATED_1KGP3).set { ch_inputUNRELATED_1KGP3 }
  Channel.fromPath(params.input1KGP3).set { ch_input1KGP3 }
  Channel.fromPath(params.inputSuper_pop_codes).set { ch_inputSuper_pop_codes }
  Channel.fromPath(params.input05both1K100K_eigenvec).set { ch_input05both1K100K_eigenvec }
  Channel.fromPath(params.inputGELprojection_proj_eigenvec).set { ch_GELprojection_proj_eigenvec }

  Channel.fromPath(params.inputPCsancestryrelated)
                        .ifEmpty { exit 1, "Input file with Michigan LD data not found at ${params.inputPCsancestryrelated}. Is the file path correct?" }
                        .set { ch_inputPCsancestryrelated }

  Channel.fromPath(params.inputAncestryAssignmentProbs)
                        .ifEmpty { exit 1, "Input file with Michigan LD data not found at ${params.inputAncestryAssignmentProbs}. Is the file path correct?" }
                        .set { ch_inputAncestryAssignmentProbs }

  Channel.fromPath(params.inputMichiganLDfileExclude)
                        .ifEmpty { exit 1, "Input file with Michigan LD for excluding regions  not found at ${params.inputMichiganLDfileExclude}. Is the file path correct?" }
                        .set { ch_inputMichiganLDfileExclude }
if (params.input.endsWith(".csv")) {

  Channel.fromPath(params.input)
                        .ifEmpty { exit 1, "Input .csv list of input tissues not found at ${params.input}. Is the file path correct?" }
                        .splitCsv(sep: ',',  skip: 1)
                        .map { bcf, index -> ['chr'+file(bcf).simpleName.split('_chr').last() , file(bcf), file(index)] }
                        .filter{bcf -> bcf =~/chr\d+/}
                        // Filter out chunks from chrX, chrY and chrM - they should not be analysed in Ancestry and Relatedness pipeline
                        .filter { it[0] =~ /chr[^XYM]/ }
                        .set { ch_bcfs }

}


// Plink files for mend_err_p* processes


// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
// TODO nf-core: Report custom parameters here
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"




 /* STEP_17
 * STEP - sort_compress: Sort and compress site metric data for KING step
 */
process sort_compress {
    publishDir "${params.outdir}/bcftools_site_metrics_subcols/", mode: params.publish_dir_mode

    input:
    tuple val(region), file(bcftools_site_metrics_subcols) from ch_bcftools_site_metrics_subcols

    output:
    tuple val(region), file("BCFtools_site_metrics_SUBCOLS${region}_sorted.txt.gz"), file("BCFtools_site_metrics_SUBCOLS${region}_sorted.txt.gz.tbi") into ch_sort_compress

    script:
    """
    sort -k2 -n ${bcftools_site_metrics_subcols} > BCFtools_site_metrics_SUBCOLS${region}_sorted.txt
    bgzip -f BCFtools_site_metrics_SUBCOLS${region}_sorted.txt && \
    tabix -s1 -b2 -e2 BCFtools_site_metrics_SUBCOLS${region}_sorted.txt.gz
    """
}
//  KING WORKFLOW

/* STEP_18
 * STEP - filter_regions: Produce BCFs of our data filtered to sites pass sites
 */

// Bcf and metrics channels are joined by region value, that carries exact identity of the bcf file and
// its metrics file, to ensure bcf files are filtered by correct metrics. Joining is done by first element
// in each tuple that carries the region value.
ch_bcf_and_metrics_joined =
    ch_bcfs.join(ch_sort_compress)

process filter_regions {
    publishDir "${params.outdir}/regionsFiltered/", mode: params.publish_dir_mode

    input:
    tuple val(region), file(bcf), file(index), file(site_metrics_file), file(site_metrics_file_index) from ch_bcf_and_metrics_joined

    output:
    tuple val(region), file("${region}_regionsFiltered.bcf") into ch_regions_filtered

    script:
    """
    bcftools view ${bcf} \
    -T ${site_metrics_file}  \
    -Ob \
    -o ${region}_regionsFiltered.bcf
    """
}


process further_filtering {
    publishDir "${params.outdir}/further_filtering/", mode: params.publish_dir_mode

    input:
    tuple val(region), file(bcf_filtered) from ch_regions_filtered

    output:
    tuple val(region), file("MichiganLD_regionsFiltered_${region}.bcf"), file("MAF_filtered_1kp3intersect_${region}.txt") into ch_further_filtering

    script:
    """
    bcftools view ${bcf_filtered} \
    -i 'INFO/OLD_MULTIALLELIC="." & INFO/OLD_CLUMPED="."' \
    -v snps  | \
    bcftools annotate \
    --set-id '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC' | \
    bcftools +fill-tags -Ob \
    -o MichiganLD_regionsFiltered_${region}.bcf \
    -- -t MAF
    #Produce filtered txt file
    bcftools query MichiganLD_regionsFiltered_${region}.bcf \
    -i 'MAF[0]>0.01' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%MAF\\n' | \
    awk -F "\t" '{ if((\$5 == "G" && \$6 == "C") || (\$6 == "G" && \$5 == "C") || (\$5 == "A" && \$6 == "T") || (\$6 == "A" && \$5 == "T")) {next} { print \$0} }' \
    > MAF_filtered_1kp3intersect_${region}.txt
    """
}

// /* STEP_20
//  * STEP - create_final_king_vcf: Produce new BCF just with filtered sites
//  */
process create_final_king_vcf {
    publishDir "${params.outdir}/create_final_king_vcf/", mode: params.publish_dir_mode

    input:
    tuple val(region), file(filtered_bcf), file(maf_filtered_variants) from ch_further_filtering
    each file(agg_samples_txt) from ch_inputFinalPlatekeys

    output:
    tuple val(region), file("${region}_filtered.vcf.gz"), file("${region}_filtered.vcf.gz.tbi") into ch_create_final_king_vcf

    script:
    """
    #Now filter down our file to just samples we want in our GRM. This removes any withdrawals that we learned of during the process of aggregation
    #Store the header
    bcftools view \
    -S ${agg_samples_txt} \
    --force-samples \
    -h ${filtered_bcf} \
    > ${region}_filtered.vcf

    #Then match against all variant cols in our subsetted bcf to our maf filtered, intersected sites and only print those that are in the variant file.
    #Then append this to the stored header, SNPRelate needs vcfs so leave as is
    bcftools view \
    -H MichiganLD_regionsFiltered_${region}.bcf \
    -S ${agg_samples_txt} \
    --force-samples \
    | awk -F '\t' '${awk_expr_create_final_king_vcf_1}' ${maf_filtered_variants} - >> ${region}_filtered.vcf
    bgzip ${region}_filtered.vcf
    tabix ${region}_filtered.vcf.gz
    """
}

/* STEP_21
 * STEP - concat_king_vcf: Concatenate compressed vcfs to per chromosome files
 */

// For concatenation we need vcf files to be grouped by chromosome.
ch_vcfs_groupped_by_chr =
    ch_create_final_king_vcf
        // Change region name to only contain chromosome number.
        .map{ region, vcf, index -> [region.split("_")[0].replace("chr",""), vcf, index] }
        // Group the tuples by chromosome number (by default the first element).
        .groupTuple()
        // The resulting channel has tuples of files coming from same chromome and should look as following:
        // [10, [chr10_10000_20000_filtered.vcf.gz, chr10_20001_30000_filtered.vcf.gz], [chr10_10000_20000_filtered.vcf.gz.tbi, chr10_20001_30000_filtered.vcf.gz.tbi]]
        // [22, [chr22_10000_20000_filtered.vcf.gz, chr22_20001_30000_filtered.vcf.gz], [chr22_10000_20000_filtered.vcf.gz.tbi, chr22_20001_30000_filtered.vcf.gz.tbi]]


process concat_king_vcf {
    publishDir "${params.outdir}/concat_king_vcf/", mode: params.publish_dir_mode

    input:
    tuple val(chr), file(vcf_files), file(vcf_file_indexes) from ch_vcfs_groupped_by_chr

    output:
    tuple val(chr), file("chrom${chr}_merged_filtered.vcf.gz"), file("chrom${chr}_merged_filtered.vcf.gz.tbi") into ch_vcfs_per_chromosome

    script:
    """
    bcftools concat ${vcf_files} \
    -Oz \
    -o chrom${chr}_merged_filtered.vcf.gz && \

    tabix chrom${chr}_merged_filtered.vcf.gz
    """
}

// /* STEP_22
//  * STEP - make_bed_all: Make BED files for 1000KGP3 intersected vcfs
//  */

process make_bed_all {
    publishDir "${params.outdir}/make_bed_all/", mode: params.publish_dir_mode

    input:
    tuple val(chr), file(vcf), file(index) from ch_vcfs_per_chromosome
    each file(michiganld_exclude_regions_file) from ch_inputMichiganLDfileExclude

    output:
    tuple val(chr), file("BED_${chr}.bed"), file("BED_${chr}.bim"), file("BED_${chr}.fam") into ch_make_bed_all

    script:
    """
    stringQuery='#-\$r/\$a-.-.'
    plink2 --vcf ${vcf} \
    --make-bed \
    --vcf-half-call m \
    --set-missing-var-ids chr@:\$stringQuery \
    --new-id-max-allele-len 60 missing \
    --exclude range ${michiganld_exclude_regions_file} \
    --double-id \
    --real-ref-alleles \
    --allow-extra-chr \
    --threads ${task.cpus} \
    --out BED_${chr}
    """

}
/* STEP_23
 * STEP - ld_bed: LD prune SNPs
 */

 process ld_bed {
    publishDir "${params.outdir}/ld_bed/", mode: params.publish_dir_mode

    input:
    tuple val(chr), file(bed), file(bim), file(fam) from ch_make_bed_all

    output:
    file("BED_LDpruned_${chr}*") into ch_ld_bed

    script:
    plink_base = bed.baseName
    """
    #Not considering founders in this as all of our SNPs are common
    plink  \
    --keep-allele-order \
    --bfile ${plink_base} \
    --indep-pairwise 500kb 1 0.1 \
    --threads ${task.cpus} \
    --out BED_LD_${chr}

    #Now that we have our correct list of SNPs (prune.in), filter the original
    #bed file to just these sites
    plink \
    --make-bed \
    --bfile ${plink_base} \
    --keep-allele-order \
    --extract BED_LD_${chr}.prune.in \
    --double-id \
    --allow-extra-chr \
    --threads ${task.cpus} \
    --out BED_LDpruned_${chr}
    """
}

/* STEP_24
 * STEP - merge_autosomes: Merge autosomes to genome wide BED files
 */

process merge_autosomes {
    publishDir "${params.outdir}/merge_autosomes/", mode: params.publish_dir_mode

    input:
    file(bed_ld_files) from ch_ld_bed.collect()

    output:
    tuple file("autosomes_LD_pruned_1kgp3Intersect.bed"),
          file("autosomes_LD_pruned_1kgp3Intersect.bim"),
          file("autosomes_LD_pruned_1kgp3Intersect.fam"),
          file("autosomes_LD_pruned_1kgp3Intersect.nosex") into (ch_merged_autosomes_hwe_pruning_30k_snps,
                                                                 ch_merged_autosomes_king_coefficients)

    script:
    """
    for i in {1..22}; do if [ -f "BED_LDpruned_\$i.bed" ]; then echo BED_LDpruned_\$i >> mergelist.txt; fi ;done
    plink --merge-list mergelist.txt \
    --make-bed \
    --out "autosomes_LD_pruned_1kgp3Intersect"
    rm mergelist.txt
    """
}

/* STEP_25
 *Purpose: produce a first pass HWE filter.
 */
process hwe_pruning_30k_snps {
    publishDir "${params.outdir}/hwe_pruning_30k_snps/", mode: params.publish_dir_mode

    input:
    tuple file(bed), file(bim), file(fam), file(nosex) from ch_merged_autosomes_hwe_pruning_30k_snps
    file(ancestry_assignment_probs) from ch_inputAncestryAssignmentProbs
    file(pc_sancestry_related) from ch_inputPCsancestryrelated

    output:
    file("hwe1e-5_superpops_195ksnps") into ch_hwe_pruning_30k_snps

    script:
    plink_base = bed.baseName
    """
    hwe_pops.R --ancestry_assignment_probs='${ancestry_assignment_probs}' \
               --pc_sancestry_related='${pc_sancestry_related}'

    for pop in AFR EUR SAS EAS; do
        echo \${pop}
        awk '{print \$1"\\t"\$1}' \${pop}pop.txt > \${pop}keep

        plink \
        --keep-allele-order \
        --make-bed \
        --bfile ${plink_base} \
        --out \${pop}

        plink --bfile \${pop} --hardy midp --out \${pop} --nonfounders
    done

    #Combine the HWE and produce a list of pass
    hwe_produce_pass.R
    """
}

/* STEP_26
 */
process king_coefficients{
    publishDir "${params.outdir}/king_coefficients/", mode: params.publish_dir_mode

    input:
    tuple file(bed), file(bim), file(fam), file(nosex) from ch_merged_autosomes_king_coefficients
    file(significant_superpops_snps) from ch_hwe_pruning_30k_snps

    output:
    tuple file("${bed.baseName}_unrelated.bed"),
          file("${bed.baseName}_unrelated.bim"),
          file("${bed.baseName}_unrelated.fam") into ch_king_coefficients_unrelated
    tuple file("${bed.baseName}_related.bed"),
          file("${bed.baseName}_related.bim"),
          file("${bed.baseName}_related.fam") into ch_king_coefficients_related
  	file("${bed.baseName}_triangle_HWE1_5.king.cutoff.in.id") into ch_unrelatedlist

    script:
    plink_base = bed.baseName
    """
    plink2 \
    --bfile ${plink_base} \
    --extract ${significant_superpops_snps} \
    --make-king triangle bin \
    --out "${plink_base}_triangle_HWE1e_5" \
    --thread-num ${task.cpus}
    echo "done1"

    plink2 \
    --bfile ${plink_base} \
    --king-cutoff "${plink_base}_triangle_HWE1e_5" 0.0442 && \
    mv "plink2.king.cutoff.in.id" "${plink_base}_triangle_HWE1_5.king.cutoff.in.id" && \
    mv "plink2.king.cutoff.out.id" "${plink_base}_triangle_HWE1_5.king.cutoff.out.id"

    plink2 \
    --bfile ${plink_base} \
    --make-bed \
    --keep "${plink_base}_triangle_HWE1_5.king.cutoff.in.id" \
    --out "${plink_base}_unrelated"
    echo "done2"

    plink2 \
    --bfile ${plink_base} \
    --make-bed \
    --remove "${plink_base}_triangle_HWE1_5.king.cutoff.in.id" \
    --out "${plink_base}_related"
    echo "done3"
    """
}
/* STEP_27a
 */
process gcta{
    publishDir "${params.outdir}/gcta/", mode: params.publish_dir_mode

    input:
    tuple file(bed), file(bim), file(fam) from ch_king_coefficients_unrelated

    output:
    tuple file("${bed.baseName}.eigenval"),
          file("${bed.baseName}.eigenvec"),
          file("${bed.baseName}.eigenvec.PROJ.eigenvec"),
          file("${bed.baseName}.grm.N.bin"),
          file("${bed.baseName}.grm.bin"),
          file("${bed.baseName}.grm.id") into ch_gcta

    script:
    plink_base = bed.baseName
    """
    gcta64 --bfile ${plink_base} \
    --make-grm-bin \
    --thread-num ${task.cpus} \
    --out ${plink_base}

    gcta64 --grm ${plink_base} \
    --pca ${n_pca} \
    --out ${plink_base} \
    --thread-num ${task.cpus}

    gcta64 --bfile ${plink_base} \
    --pc-loading ${plink_base} \
    --out ${plink_base} \
    --thread-num ${task.cpus}

    awk 'BEGIN{OFS="    "}{print \$0, "NA"}' "${plink_base}.eigenvec" > "${plink_base}.eigenvec.PROJ.eigenvec"
    """
}
/* STEP_27
 */
process infer_ancestry{
    publishDir "${params.outdir}/infer_ancestry/", mode: params.publish_dir_mode

    input:
    file(kgp3_sample_table) from ch_input1KGP3
    file(super_pop_codes) from ch_inputSuper_pop_codes
    file(kgp3_unrel) from ch_inputUNRELATED_1KGP3
    file(eigenvec) from ch_input05both1K100K_eigenvec
    file(projections) from ch_GELprojection_proj_eigenvec
    tuple file("autosomes_LD_pruned_1kgp3Intersect_unrelated.eigenval"), file("autosomes_LD_pruned_1kgp3Intersect_unrelated.eigenvec"), file("autosomes_LD_pruned_1kgp3Intersect_unrelated.eigenvec.PROJ.eigenvec"), file("autosomes_LD_pruned_1kgp3Intersect_unrelated.grm.N.bin"), file("autosomes_LD_pruned_1kgp3Intersect_unrelated.grm.bin"), file("autosomes_LD_pruned_1kgp3Intersect_unrelated.grm.id") from ch_gcta

    output:
    file("predicted_ancestries.tsv") into ch_infer_ancestry
    file("results.RDS") into ch_infer_ancestry_2

    script:
    """
    mkdir Ancestries
    infer_ancestry.R  \
    ${kgp3_sample_table} \
    ${super_pop_codes} \
    ${kgp3_unrel} \
    ${eigenvec} \
    ${projections} \
    "Ancestries"
    """
}

/*
 * Completion notification
 */
workflow.onComplete {

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[lifebit-ai/annotate]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        log.info "-${c_purple}[lifebit-ai/annotate]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """
    """.stripIndent()
}