#!/bin/bash


############## CONFIG ##############
#Filepaths in here must be hardcoded due to restrictions in passing arguments to script when
#submitting as jobs
#
#NOTE - BSUB wd header is hardcoded, CWD will have to be altered in run_site_metrics.sh header
#     - Please include trailing '/' in dirpaths
#
#Input dir containing all the initially aggregated chunks
#input='/re_gecip/BRS/release_8_aggregation/workflow_output/data/'
input='/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_chr18/subsample_60k_match/'
#
#Output - Where you want output docs to go. 
out='/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_chr18_subsample/'
#
#Resources
resources='/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/resources/'
#This should contain the following
# 1. List of trios - see trio_define.R for how this is used
# 2. UPD cases - these are filtered out of trio data
# 3. VCF header file named vcfheader.hdr
#These will be added to in the generation of .keep and .fam file derived from the list of trios
#
#Trio data file
triofile="extended_interpreted_trios_grch38.tsv"
#
#Samples to be included
includedSamples="innerjoin_samples.txt"
#
#The chr1_pos1_pos2 for each chunk
chr_pos_pos="chunks.txt"
#
#everything that comes before chr_pos1_pos2 in the filename
mainname="all_chunks_merged_norm"
#
#Define working directory (for job submission)
wd='/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/src/'
#
#Limit to specific chrom - mostly for testing
specChrom='18'
#
#Path to R packages (see note below)
rlibs='/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/resources/R_packages/'

#Sample lists for XX/XY QC
xx='xx_females_illumina_ploidy_samples_40740.tsv'
xy='xy_males_illumina_ploidy_samples_35924.tsv'

######## NOTE ########
# For all of this to work, certain R packages have to be installed on the HPC.
# Currently I am storing these in resources/R_packages
# SNPRelate
# gdsfmt
# GENESIS
####################################