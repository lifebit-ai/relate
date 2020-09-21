#!/bin/bash

#BSUB -P bio
#BSUB -q pipeline
#BSUB -cwd /re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/src_actual/
#BSUB -o logs/site_metrics.out.%J.%I
#BSUB -e logs/site_metrics.err.%J.%I
#BSUB -J "site_metrics_master"


### Start ###
echo 'Please note cwd is defined in script header for LSF'

#Check we have access to functions file
if [ ! -f sample_metrics.sh ]; then
    echo 'Check loc of functions file (site_metrics.sh) and running from src' >&2
    exit 1
fi

if [ ! -f config.sh ]; then
    echo 'Check loc of config file (config.sh) and running from src' >&2
    exit 1
fi

if [ ! -w logs ]; then
    echo 'Please create 'logs' file in src/' >&2
    exit 1
fi
savepath=$LD_LIBRARY_PATH
source config.sh
source sample_metrics.sh

## Export config vars ##
export input
export out
export resources
export triodata="${resources}/${triofile}"
#Samples to be included
export aggregateSamples="${resources}/${includedSamples}"
#The chr1_pos1_pos2 for each chunk
export bedfile="${resources}/${chr_pos_pos}"
#everything that comes before chr_pos1_pos2 in the filename
export mainname
export wd
##

## Test inputs ##
if [ ! -d ${input} ]; then
    echo 'Check input directory exists' >&2
    exit 1
fi

if [ ! -f "${bedfile}" ]; then
    echo 'Check filepath to .bed file' >&2
    exit 1
fi

if [ ! -f sample_metrics.sh ]; then
    echo 'File "sample_metric_funcs.sh" not found, aborting.' >&2
    exit 1
fi

if [ ! -r ${out} ]; then
    echo 'Please ensure input is a dir and is readable.' >&2
    exit 1
fi

if [ ! -r ${resources} ]; then
    echo 'Please ensure input is a dir and is readable.' >&2
    exit 1
fi

if [ ! -w ${out} ]; then
    echo 'Please ensure output is a dir and is writeable.' >&2
    exit 1
fi

if [ ! -s ${triodata} ]; then
    echo 'Please ensure triodata file exists in the resources dir.' >&2
    exit 1
fi

if [ ! -s ${aggregateSamples} ]; then
    echo 'Please ensure aggregate samples file exists in the resources dir.' >&2
    exit 1
fi

## End test ##

## Work start ##
#If we have defined a specific chromosome, filter to these only
if [ ! -z "$specChrom" ]; then
    regions=$(grep -v .csi ${bedfile} | grep "chr${specChrom}"  | grep -oP "(?<=${mainname}_)chr\d+_\d+_\d+(?=\.bcf$)")
else
    regions=$(grep -v .csi ${bedfile} | grep -oP "(?<=${mainname}_)chr\d+_\d+_\d+(?=\.bcf$)")
fi
echo 'Generating .keep and .fam file'
export -f triodata_define
bsub -J "triodefine" -q short -P bio -U as1 -cwd ${wd} -e logs/triodefine.err_%J -o logs/triodefine.out${i}_%J  triodata_define


#Export functions as environment variables so they can be submitted as jobs
export -f startFile
export -f missingness1
export -f missingness2
export -f medianCoverageAll
export -f ABRatioP1
export -f ABRatioP2
export -f completeSites
export -f medianGQ
export -f medianCoverageNonMiss


# export -f makeBed

# export -f aggregateAnnotation

### MENDEL ERRORS ###
export -f MendErrp1
export -f MendErrp2
export -f MendelDist
export -f MendErrp3
export LD_LIBRARY_PATH=$savepath

for i in $regions; do
    export i
    #Parse infile paths
    export infile=${input}/${mainname}_${i}.bcf
    #Mendel errors
    bsub -J "mendelErrp1${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/mendelErrp1.err${i}_%J -o logs/mendelErrp1.out${i}_%J  MendErrp1
    bsub -J "mendelErrp2${i}" -q long -P bio -U as1 -w "done(mendelErrp1${i})" -cwd ${wd} -e logs/mendelErrp2.err${i}_%J -o logs/mendelErrp2.out${i}_%J  MendErrp2
  
done

#Then we need to find the 'good families'. 

bsub -J "mendelDist" -q long -P bio -U as1 -cwd ${wd} -e logs/mendelDist.err_%J -o logs/mendelDist.out_%J  MendelDist

#Then we calculate the site mendel errors for this set only
for i in $regions; do
    export i
    #Parse infile paths
    export infile=${input}/${mainname}_${i}.bcf
    bsub -J "MendErrSites_${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/MendErrSites_err${i}_%J -o logs/MendErrSites.out${i}_%J MendErrp3
done

###
echo 'Calculating site metrics...\n'
#explicit rerun of failed chunks for 'chr8_141801176_145138636'

for i in $regions; do
    export i
    #Parse infile paths
    export infile=${input}/${mainname}_${i}.bcf
    #Begin jobs
    
    #Write backbone
    bsub -J "startFile${i}" -q short -P bio -U as1 -cwd ${wd} -e logs/startFile.err${i}_%J -o logs/Startfile.out${i}_%J startFile
    
    #Missingness
    bsub -J "missingSite1_${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/missingness1.err${i}_%J -o logs/missingness1.out${i}_%J missingness1
    bsub -J "missingSite2_${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/missingness2.err${i}_%J -o logs/missingness2.out${i}_%J missingness2
    
    #Coverage
    bsub -J "medianCovNonMiss_${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/meanCoverageNonMiss.err${i}_%J -o logs/meanCoverageNonMiss.out${i}_%J medianCoverageNonMiss
    bsub -J "medianCoverageAll_${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/medianCoverageAll.err${i}_%J -o logs/medianCoverageAll.out${i}_%J medianCoverageAll
    
    #AB ratio
    bsub -J "ABRatio1${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/ABratio1.err${i}_%J -o logs/ABratio1.out${i}_%J ABRatioP1
    bsub -J "ABRatio2${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/ABratio2.err${i}_%J -o logs/ABratio2.out${i}_%J ABRatioP2
    
    #GQ
    bsub -J "medianGQ_${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/medianGQ.err${i}_%J -o logs/medianGQ.out${i}_%J medianGQ
    
    #Complete sites
    bsub -J "completeSites_${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/completeSites.err${i}_%J -o logs/completeSites.out${i}_%J completeSites
    
    
done

echo 'Entering relatedness inference workflow...\n'

#This is the second run at filtering
export -f pullAC
export -f pull_1KGPSites
export -f aggregateAnnotation
export -f regions_filter
export -f sort_compress
export -f final_KING_BCF
export king='T'
export LD_LIBRARY_PATH=$savepath
for i in $regions; do
    export i
    #Parse infile paths
    export infile=${input}/${mainname}_${i}.bcf
    #bsub -q long -P bio -U as1 -e logs/pullAC.err${i}_%J -o logs/pullAC.out${i}_%J  pullAC
    #bsub -J pull1kgp3${i} -q long -P bio -U as1 -e logs/pull1kgp3.err${i}_%J -o logs/pull1kgp3.out${i}_%J  pull_1KGPSites
    #bsub -J annotate${i} -q long -P bio -U as1 -e logs/annotate.err${i}_%J -o logs/annotate.out${i}_%J  aggregateAnnotation
    #bsub -J compress${i} -q long -P bio -U as1 -e logs/sort_compress.err${i}_%J -o logs/sort_compress.out${i}_%J  sort_compress
    #bsub -J regionsfilter${i} -q  long -P bio -U as1 -e logs/regionsFilter.err${i}_%J -o logs/regionsFilter.out${i}_%J  regions_filter
    #bsub -J furtherfilter${i} -q long -P bio -U as1 -e logs/furtherFilter.err${i}_%J -o logs/furtherFilter.out${i}_%J  further_filtering
    bsub -J final_filter${i} -q long -P bio -U as1 -e logs/finalFilter.err${i}_%J -o logs/finalFilter.out${i}_%J  final_KING_BCF
done
unset king

export -f concat_KING_VCF
export -f makeBedAll
export -f LD_BED
for i in {1..22}; do
    export i
    #bsub -J concatVCF${i} -q long -P bio -U as1 -e logs/concatVCF.err${i}_%J -o logs/concatVCF.out${i}_%J  concat_KING_VCF
    #bsub -J BEDcreate${i} -q long -P bio -U as1 -e logs/BEDcreate.err${i}_%J -o logs/BEDcreate.out${i}_%J  makeBedAll
    bsub -J LDprune${i} -q long -P bio -U as1 -e logs/LDprune.err${i}_%J -o logs/LDprune.out${i}_%J  LD_BED
done
export -f merge_autosomes
export -f king_coefficients
bsub -J merge_autos -q long -P bio -U as1 -e logs/mergeAutos.err${i}_%J -o logs/mergeAutos.out${i}_%J  merge_autosomes
bsub -J king_coefficients -q long -P bio -U as1 -e logs/king.err${i}_%J -o logs/king.out${i}_%J  king_coefficients
#bsub -J pcAIR -q long -P bio -U as1 -e logs/pcair.err_%J -o logs/pcair.out_%J 'module load $RLoad; echo $RLoad; script pc_air.R'

#Now let's do the HWE for each of the chunks
export LD_LIBRARY_PATH=$savepath
export -f  p_hwe
for i in $regions; do
    export i
    #Parse infile paths
    export infile=${input}/${mainname}_${i}.bcf
    bsub -J p_hwe${i} -q long -P bio -U as1 -e logs/hwe.err${i}_%J -o logs/her.out${i}_%J p_hwe
done


### SEX QC ###
export LD_LIBRARY_PATH=$savepath
export sexChrom='true'
export -f completeSitesSexChrom
Sregions=$(grep -v .csi ${bedfile} | grep -oP "(?<=${mainname}_)chr[X]_\d+_\d+(?=\.bcf$)")
for i in $Sregions; do
    export xx #defined in config, files to xx and xy data
    export xy
    export i
    #Parse infile paths
    export infile=${input}/${mainname}_${i}.bcf
    #Begin jobs
    
    #Write backbone
    bsub -J "startFile${i}" -q short -P bio -U as1 -cwd ${wd} -e logs/startFile.err${i}_%J -o logs/Startfile.out${i}_%J startFile
    
    #Missingness
    bsub -J "missingSite1_${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/missingness1.err${i}_%J -o logs/missingness1.out${i}_%J missingness1
    bsub -J "missingSite2_${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/missingness2.err${i}_%J -o logs/missingness2.out${i}_%J missingness2
    
    #Coverage
    bsub -J "medianCovNonMiss_${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/meanCoverageNonMiss.err${i}_%J -o logs/meanCoverageNonMiss.out${i}_%J medianCoverageNonMiss
    bsub -J "medianCoverageAll_${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/medianCoverageAll.err${i}_%J -o logs/medianCoverageAll.out${i}_%J medianCoverageAll
    
    #AB ratio
    bsub -J "ABRatio1${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/ABratio1.err${i}_%J -o logs/ABratio1.out${i}_%J ABRatioP1
    bsub -J "ABRatio2${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/ABratio2.err${i}_%J -o logs/ABratio2.out${i}_%J ABRatioP2
    
    #GQ
    bsub -J "medianGQ_${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/medianGQ.err${i}_%J -o logs/medianGQ.out${i}_%J medianGQ
done

export -f sexChromAnnotation
for i in $Sregions; do
    export i
    bsub -J sexAnnotate${i} -q long -P bio -U as1 -e logs/sexAnnotate.err${i}_%J -o logs/sexAnnotate.out${i}_%J  sexChromAnnotation
done


#Run check for annotation files, then run final annotation
export -f check_annotation_length
export -f endAggregateAnnotation
for i in $regions; do
    export i
    #Parse infile paths
    export infile=${input}/${mainname}_${i}.bcf
    #bsub -J checking_${i} -q long -P bio -U as1 -e logs/checking.err${i}_%J -o logs/checking.out${i}_%J check_annotation_length
    bsub -J final_annotate_${i} -q long -P bio -U as1 -e logs/final.err${i}_%J -o logs/final.out${i}_%J endAggregateAnnotation
done
#Run the header script in sample_metrics.sh

#Now annotate the vcfs




### VEP ANNOTATION ###
#Run me on pegasus for now
export -f VEPannotate
for i in $regions; do
    export i
    export infile=${input}/${mainname}_${i}.bcf
    bsub -J "VEP_annotate${i}" -R "span[hosts=1]" -q long -P bio -W 336:0 -cwd ${wd} -o logs/VEP_annotate.out${i}_%J -e logs/VEP_annotate.err${i}_%J VEPannotate
done
#Now just add the splice AI data


# ## End work ##
# ### END ###