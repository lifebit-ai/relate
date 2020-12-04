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
        --input [file]                  File with list of full paths to bcf files and their indexes.
                                        Bcf files can be compressed but in a readable for bcftools format.
                                        Example:
                                        #-----my_bcf_files_list.csv-----------#
                                        | bcf,index                           |
                                        | <file1.bcf>,<file1.bcf.idx>         |
                                        | <file2.bcf.gz>,<file2.bcf.gz.csi>   |
                                        | <file3.bcf.bgz>,<file3.bcf.bgz.tbx> |
                                        #-------------------------------------#
                                        The name of the files must be consistent across files
                                        and follow a specific pattern:
                                        {name}_{CHR}_{START_POS}_{END_POS}.bcf.gz
                                        Example:
                                        test_all_chunks_merged_norm_chr10_53607810_55447336.bcf.gz
                                        Consistency is important here as a variable ('region')
                                        is extracted from the filename.

        --included_samples [file]       File with a list of participant IDs to be included in the analysis.

        --siteqc_results_dir [dir]      The results directory of SiteQC pipeline, that has all metrics files
                                        generated from bcf/vcf files provided with --input option. Only final siteqc
                                        annotation files from sub-folder "Annotate" will be used in this pipeline.

        --michigan_ld_file [file]       File with regions to be filtered out for improving quality of sites selected.

        --ancestry_probs [file]         File required for hwe_pruning_30k_snps process containing tab separated values
                                        for probabilities of assignments for the 31 populations code for each participant.

        --pcs_ancestry [file]           File with Principal Components information comming from reference resources of GEL
                                        for the inferred ancestries from the 30k dataset.

      #Example files temporarily required:

        --unrelated_1kgp3 [file]        File required for the infer_ancestry process, is a two column tab separated file
                                        with platekeys on each column.
        --samplelist_1kgp3 [file]       File required for the infer_ancestry process, is a three column tab separated file
                                        with Sample(Platekey), Family ID and Population asignments.
        --example_eigenvec [file]       File required for the infer_ancestry process, its tab separated file with eigenvectors.
        --example_proj_eigenvec [file]  File required for the infer_ancestry process, its tab separated file with eigenvectors.
                                        and it has "NA" values in the last column.


    Optional arguments:
        --super_pop_codes               Tab separated file with population and super population relations and related information.
                                        Required for the infer_ancestry process. It is a general reference file.
                                        By default the following file is used:
                                        s3://lifebit-featured-datasets/projects/gel/siteqc/super_pop_codes.tsv

        --n_pca [integer]               Number of Principal Components desired for gcta process. (Default: 20)


    Other options:
        -profile [str]                  Configuration profile to use. Can use multiple (comma separated).
                                        Available: standard, test. (Default: standard)
        --outdir [file]                 The output directory where the results will be saved.
        --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move. (Default: copy)
        -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    """.stripIndent()
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

/*
 * Check all important required inputs
 */

// Check if user provided input csv file containing paths to bcf files and their indexes
if (!params.input) exit 1, "The list of input bcf/vcf files was not provided. \nPlease specify a csv file containing paths to bcf/vcf files and their indexes with --input [file] option. \nUse --help option for more information."
if (!params.included_samples) exit 1, "The list of participant IDs was not provided. \nPlease specify a text file containing participant IDs (platekeys) with --included_samples [file] option. \nUse --help option for more information."
if (!params.siteqc_results_dir) exit 1, "The SiteQC results directory was not provided. \nPlease specify a path to SiteQC results directory with --siteqc_results_dir [path] option. \nUse --help option for more information."
if (!params.michigan_ld_file) exit 1, "The file specifying genomic regions to exclude from Relatedness and Ancestry inference was not provided. \nPlease specify such file with --michigan_ld_file [file] option. \nUse --help option for more information."
if (!params.ancestry_probs) exit 1, "The file with Ancestry assignment probabilities was not provided. \nPlease specify such file with --ancestry_probs [file] option. \nUse --help option for more information."
if (!params.pcs_ancestry) exit 1, "The files with Ancestry Principal Components was not provided. \nPlease specify such file with --pcs_ancestry [file] option. \nUse --help option for more information."


// Define channels based on params

// Input list .csv file of tissues to analyse
// [chr10_52955340_55447336, test_all_chunks_merged_norm_chr10_52955340_55447336.bcf.gz, test_all_chunks_merged_norm_chr10_52955340_55447336.bcf.gz.csi]
ch_bcfs = Channel.fromPath(params.input)
            .ifEmpty { exit 1, "Input .csv list of input tissues not found at ${params.input}. Is the file path correct?" }
            .splitCsv(sep: ',',  skip: 1)
            .filter{bcf -> bcf =~/chr\d+/} //only autosmes are analyzed in Ancestry and Relatedness pipeline.
            .map { bcf, index -> ['chr'+file(bcf).simpleName.split('_chr').last() , file(bcf), file(index)] }

ch_included_samples = Channel.fromPath(params.included_samples)
            .ifEmpty { exit 1, "Input file with samples and platekeys data not found at ${params.included_samples}. Is the file path correct?" }

ch_bcftools_site_metrics_subcols = Channel.fromPath(params.siteqc_results_dir + '/Annotation/*.txt')
            .ifEmpty { exit 1, "Input annotation txt files not found at ${params.siteqc_results_dir}/Annotation. Is the dir path correct?" }
            .filter{txt -> txt =~/chr\d+/} //only autosmes are analyzed in Ancestry and Relatedness pipeline.
            .map { txt -> ['chr'+ txt.simpleName.split('_chr').last() , txt] }

ch_super_pop_codes = Channel.fromPath(params.super_pop_codes)
            .ifEmpty { exit 1, "Input file with superpopulation codes was not found at ${params.super_pop_codes}. Is the file path correct?" }

ch_michigan_ld_file = Channel.fromPath(params.michigan_ld_file)
            .ifEmpty { exit 1, "Input file with Michigan LD for excluding regions  not found at ${params.michigan_ld_file}. Is the file path correct?" }

ch_ancestry_probs = Channel.fromPath(params.ancestry_probs)
            .ifEmpty { exit 1, "Input file with ancestry assignment probabilities not found at ${params.ancestry_probs}. Is the file path correct?" }

ch_pcs_ancestry = Channel.fromPath(params.pcs_ancestry)
            .ifEmpty { exit 1, "Input file with ancestry Principal Components not found at ${params.pcs_ancestry}. Is the file path correct?" }


ch_unrelated_1kgp3 = Channel.fromPath(params.unrelated_1kgp3)
ch_samplelist_1kgp3 = Channel.fromPath(params.samplelist_1kgp3)
ch_example_eigenvec = Channel.fromPath(params.example_eigenvec)
ch_example_proj_eigenvec = Channel.fromPath(params.example_proj_eigenvec)



// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
// TODO nf-core: Report custom parameters here
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Input bcf list']   = params.input
summary['Included samples'] = params.included_samples
summary['Siteqc results dir']=params.siteqc_results_dir
summary['Michigan LD file'] = params.michigan_ld_file
summary['Ancestry probs']   = params.ancestry_probs
summary['PCs ancestry']     = params.pcs_ancestry
summary['Super pop codes']  = params.super_pop_codes
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
    publishDir "${params.outdir}/regions_filtered/", mode: params.publish_dir_mode

    input:
    tuple val(region), file(bcf), file(index), file(site_metrics_file), file(site_metrics_file_index) from ch_bcf_and_metrics_joined

    output:
    tuple val(region), file("${region}_regions_filtered.bcf") into ch_regions_filtered

    script:
    """
    bcftools view ${bcf} \
    -T ${site_metrics_file}  \
    -Ob \
    -o ${region}_regions_filtered.bcf
    """
}


process further_filtering {
    publishDir "${params.outdir}/further_filtering/", mode: params.publish_dir_mode

    input:
    tuple val(region), file(bcf_filtered) from ch_regions_filtered

    output:
    tuple val(region), file("MichiganLD_regions_filtered_${region}.bcf"), file("MAF_filtered_1kp3intersect_${region}.txt") into ch_further_filtering

    script:
    """
    bcftools view ${bcf_filtered} \
    -i 'INFO/OLD_MULTIALLELIC="." & INFO/OLD_CLUMPED="."' \
    -v snps  | \
    bcftools annotate \
    --set-id '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC' | \
    bcftools +fill-tags -Ob \
    -o MichiganLD_regions_filtered_${region}.bcf \
    -- -t MAF
    #Produce filtered txt file
    bcftools query MichiganLD_regions_filtered_${region}.bcf \
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
    each file(agg_samples_txt) from ch_included_samples

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
    -H ${filtered_bcf} \
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
    each file(michiganld_exclude_regions_file) from ch_michigan_ld_file

    output:
    tuple val(chr), file("BED_${chr}.bed"), file("BED_${chr}.bim"), file("BED_${chr}.fam") into ch_make_bed_all

    script:
    """
    string_query='#-\$r/\$a-.-.'
    plink2 --vcf ${vcf} \
    --make-bed \
    --vcf-half-call m \
    --set-missing-var-ids chr@:\$string_query \
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
    file(ancestry_probs) from ch_ancestry_probs
    file(pc_sancestry_related) from ch_pcs_ancestry

    output:
    file("hwe1e-5_superpops_195ksnps") into ch_hwe_pruning_30k_snps

    script:
    plink_base = bed.baseName
    """
    hwe_pops.R --ancestry_assignment_probs='${ancestry_probs}' \
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
    file(kgp3_sample_table) from ch_samplelist_1kgp3
    file(super_pop_codes) from ch_super_pop_codes
    file(kgp3_unrel) from ch_unrelated_1kgp3
    file(eigenvec) from ch_example_eigenvec
    file(projections) from ch_example_proj_eigenvec
    tuple file("autosomes_LD_pruned_1kgp3Intersect_unrelated.eigenval"), file("autosomes_LD_pruned_1kgp3Intersect_unrelated.eigenvec"), file("autosomes_LD_pruned_1kgp3Intersect_unrelated.eigenvec.PROJ.eigenvec"), file("autosomes_LD_pruned_1kgp3Intersect_unrelated.grm.N.bin"), file("autosomes_LD_pruned_1kgp3Intersect_unrelated.grm.bin"), file("autosomes_LD_pruned_1kgp3Intersect_unrelated.grm.id") from ch_gcta

    output:
    file("predicted_ancestries.tsv") into ch_infer_ancestry
    file("results.RDS") into ch_infer_ancestry_r_data

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