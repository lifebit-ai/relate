#!/usr/bin/env bash
/**
 * @ Author: Daniel Rhodes
 * @ Create Time: 2020-09-15 11:44:20
 * @ Description: Run funcs from siteQC_funcs.sh
 */


#Similar to siteQC_funcs.sh, this is the script meant to be a cleanedup version 
#of the aggv2 run. This is not the script actually used to generate the aggv2 (see run_site_metrics.sh)

set -u
### Config ###
#Define the following when changing locations and datasets.
#Make sure input, resources and out are directories that exist
export wd='/re_gecip/shared_allGeCIPs/drhodes/covid_siteqc/testing/src/'
export input='/re_gecip/BRS/covid_dragen_aggregation/workflow_output/covid_dragen_aggregation/temp_output/'
export out="${wd//src/output}"
export resources="${wd//src/resources}"
export bedfile="${resources}autosomes_names.txt"
export mainname="all_chunks_merged_norm"
export sampleList="final_samplelist.txt"
export michiganList="MichiganLD_liftover_exclude_regions_PARSED.txt"
export snps30k="aggV2_R9_M30K_1KGP3_ancestry_assignment_probs.tsv"
export pcs30k="aggV2_bedmerge_30KSNPs_labkeyV9_08062020_update_PCsancestryrelated.tsv"
export 1kgp3_sample_table="1KGP3.sample_table"
export super_pop_codes="super_pop_codes.tsv"
export 1kgp3_unrel="1KGP3_30K_unrel_autosomes.fam"

###

#Contains the site QC functions
source siteQC_funcs.sh

echo 'Loading modules for Helix HPC'
    export bcftoolsLoad='bio/BCFtools/1.10.2-foss-2018b'
    export RLoad='lang/R/3.6.0-foss-2019a'
    export plinkLoad='bio/PLINK/1.9b_4.1-x86_64'
    export plink2Load='bio/PLINK/2.00-devel-20200409-x86_64'


mkdir -p ${wd}/logs
#### Start ####
## Run tasks that need only be done once across all chunks ##

# Trio define
bsub -J "triodefine" -q short -P bio -U as1 -cwd ${wd} -e logs/triodefine.err_%J -o logs/triodefine.out_%J  triodata_define

#Complete sites
bsub -J "completeSites" -q long -P bio -U as1 -cwd ${wd} -e logs/completeSites.err_%J -o logs/completeSites.out_%J completeSites

#Pull 1KGP3 sites, 1 job per autosome
bsub -J "pull_1KGPSites[1-22]" -q long -P bio -U as1 -cwd ${wd} -e logs/pull_1KGPSites.err_%J -o logs/pull_1KGPSites.out_%J pull_1KGPSites

##

## The following must be run across all chunks ##
bsub -J "startFile[1-1292]" -q long -P bio -U as1 -cwd ${wd} -e logs/startFile.err_%J -o logs/startFile.out_%J startFile
bsub -J "missingness1[1-1292]" -q long -P bio -U as1 -cwd ${wd} -e logs/missingness1.err_%J -o logs/missingness1.out_%J missingness1
bsub -J "missingness2[1-1292]" -q long -P bio -U as1 -cwd ${wd} -e logs/missingness2.err_%J -o logs/missingness2.out_%J missingness2
bsub -J "medianCovAll[1-1292]" -q long -P bio -U as1 -cwd ${wd} -e logs/medianCovAll.err_%J -o logs/medianCovAll.out_%J medianCovAll
bsub -J "medianCovNonMiss[1-1292]" -q long -P bio -U as1 -cwd ${wd} -e logs/medianCovNonMiss.err_%J -o logs/medianCovNonMiss.out_%J medianCovNonMiss
bsub -J "medianGQ[1-1292]" -q long -P bio -U as1 -cwd ${wd} -e logs/medianGQ.err_%J -o logs/medianGQ.out_%J medianGQ
bsub -J "ABRatioP1[1-1292]" -q long -P bio -U as1 -cwd ${wd} -e logs/ABRatioP1.err_%J -o logs/ABRatioP1.out_%J ABRatioP1
bsub -J "ABRatioP2[1-1292]" -q long -P bio -U as1 -cwd ${wd} -e logs/ABRatioP2.err_%J -o logs/ABRatioP2.out_%J ABRatioP2
##

#Upon completion of the above, aggregate the results to create annotation file
bsub -J "aggregateAnnotation[1-1292]" -q long -P bio -U as1 -cwd ${wd} -e logs/aggregateAnnotation.err_%J -o logs/aggregateAnnotation.out_%J aggregateAnnotation

#Now we have the files containing high threshold pass sites
#Now we want to find the HQ SNPs
#First convert to BED format
bsub -J "bedfile[1-1292]" -q long -P bio -U as1 -cwd ${wd} -e logs/bedfile.err_%J -o logs/bedfile.out_%J tobed_file
bsub -J "hq_extract[1-1292]" -q long -P bio -U as1 -cwd ${wd} -e logs/hq_extract.err_%J -o logs/hq_extract.out_%J hq_extract
