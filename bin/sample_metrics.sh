#!/bin/bash
/**
 * @ Author: Daniel Rhodes
 * @ Create Time: 2020-09-15 11:16:12
 * @ Description: Site QC metrics for the aggV2
 */


#This is the version being created as a cleaned up version of aggV2 run. Actual original code is found in sample_metrics.sh, which
#serves as an indicator of how the code was actually run. This is an effort to tidy up the code an optimise.
#TODO
#Optimise median calcs by removing sed command and incorporating into awk

#Site QC functions

#Define trios
triodata_define(){  
    ###############
    ##  Purpose: Create a .fam file of confirmed trios in our samples
    ##  Input: List of trios, samples included in agg
    ##  Output: .fam and .keep files
    ###############
    module load $RLoad
    Rscript trio_define.R $triodata $aggregateSamples $triodata.fam $triodata.keep
}

completeSites(){
    ###############
    ##  Purpose: Make sure the number of samples is listed in resources
    ##  Input: BCF
    ##  Output: Txt file with N samples
    ###############
    module load $bcftoolsLoad
    export infile=`sed -n "1p" $bedfile`
    #All this will do is create a file containing the number of samples to be used
    #by the annotation script
    if [ ! -f "${resources}/N_samples" ]; then
        bcftools query -l ${input}${infile} | wc -l > ${resources}/N_samples
    fi
}

startFile(){
    ###############
    ##  Purpose: Create a backbone of IDs for other data to be joined to
    ##  Input: BCF
    ##  Output: txt file with IDs
    ###############
    module load $bcftoolsLoad
    mkdir -p ${out}startfile
    echo 'Creating backbone file'
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi

    # Query chrom and pos for annotation file
    if "$sexChrom"; then
        #XY
        bcftools query -f '%CHROM %POS %REF %ALT %INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC\n' ${input}{${infile} -S ${resources}${xy} --output ${out}startfile/start_file_${i}_XY
        #XX
        bcftools query -f '%CHROM %POS %REF %ALT %INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC\n' ${input}${infile} -S ${resources}${xx} --output ${out}startfile/start_file_${i}_XX
    else
        #Autosomes
        bcftools query -f '%CHROM %POS %REF %ALT %INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC\n' ${input}${infile} --output ${out}startfile/start_file_${i}
    fi
}

missingness1(){
    ###############
    ##  Purpose: Missingness step 1, count fully missing GTs
    ##  Input: BCF
    ##  Output: Txt file with ID and count
    ###############
    module load $bcftoolsLoad
    echo 'Calculating missing sites'
    mkdir -p ${out}missing
    module load $bcftoolsLoad

    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
    if "$sexChrom"; then
        #XY
        bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -i 'GT="./." & FORMAT/DP=0'  -S ${resources}${xy} \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing1_${i}_XY
        #XX
        bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -i 'GT="./." & FORMAT/DP=0'  -S ${resources}${xx} \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing1_${i}_XX
    else
        #autosomes
        bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -i 'GT="./." & FORMAT/DP=0' \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing1_${i}
    fi
}

missingness2(){
    ###############
    ##  Purpose: Count complete GTs only
    ##  Input: BCF
    ##  Output: Txt file with ID and count
    ###############
    module load $bcftoolsLoad
    echo 'Calculating missing sites'
    mkdir -p ${out}missing
    module load $bcftoolsLoad

    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
    if "$sexChrom"; then
        #XY
        bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -e 'GT~"\."' -S ${resources}${xy} \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing2/missing2_${i}_XY
        #XX
        bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -e 'GT~"\."' -S ${resources}${xx} \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing2_${i}_XX
    else
        #Autosomes
        bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -e 'GT~"\."' \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing2_${i}
    fi
}

medianCovAll(){
    ###############
    ##  Purpose: Produce median value for depth across all GT
    ##  Input: BCF
    ##  Output: Txt file with ID and median depth
    ###############
    module load $bcftoolsLoad
    mkdir -p ${out}/medianCovAll
    echo 'Calculating median depth for all GT...'
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi

    if "$sexChrom"; then
        #XY
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' ${input}${infile} -S ${resources}${xy} |  \
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
                print $1"\t"median}' > ${out}medianCovAll/medianCov_${i}
        #XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' ${input}${infile} -S ${resources}${xx}|  \
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
                print $1"\t"median}' > ${out}medianCovAll/medianCov_${i}
    else
        #Autosomes
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' ${input}${infile} |  \
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
                print $1"\t"median}' > ${out}medianCovAll/medianCov_${i}
    fi
}

medianCovNonMiss(){
    module load $bcftoolsLoad
    mkdir -p ${out}/medianCovNonMiss
    echo 'Calculating median depth for non missing GT...'
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi

    if "$sexChrom"; then
        #XY
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' -e 'GT~"\."' ${input}${infile} -S ${resources}${xy} |  \
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
                print $1"\t"median}' > ${out}medianCovNonMiss/medianCovNonMiss_${i}
        #XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' -e 'GT~"\."' ${input}${infile} -S ${resources}${xx}|  \
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
                print $1"\t"median}' > ${out}medianCovNonMiss/medianCovNonMiss_${i}
    else
        #Autosomes
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' -e 'GT~"\."' ${input}${infile} |  \
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
                print $1"\t"median}' > ${out}medianCovNonMiss/medianCovNonMiss_${i}
    fi
}

medianGQ(){
    ###############
    ##  Purpose: Calculate median GQ
    ##  Input: BCF
    ##  Output: txt file - ID and median
    ###############
    module load $bcftoolsLoad
    echo 'Calculating median GQ for non missing GT...'
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
    mkdir -p ${out}/medianGQ

    if "$sexChrom"; then
        #XY
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GQ ] \n' -e 'GT~"\."' ${input}${infile} -S ${resources}${xy} |\
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2; if(median > 99) {median=99};
        print $1"\t"median}' > ${out}medianGQ/medianGQ_${i}_XY
        #XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GQ ] \n' -e 'GT~"\."' ${input}${infile} -S ${resources}${xx} |\
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2; if(median > 99) {median=99};
                print $1"\t"median}' > ${out}medianGQ/medianGQ_${i}_XX
    else
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GQ ] \n' -e 'GT~"\."' ${input}${infile} |\
        sed s/[[:space:]]\\./\ 0/g | \
        awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2; if(median > 99) {median=99};
                print $1"\t"median}' > ${out}medianGQ/medianGQ_${i}
    fi
}

ABRatioP1(){
    ###############
    ##  Purpose: AB ratio calculation - number of hets passing binomial test (reads supporting het call)
    ##  Input: BCF
    ##  Output: Txt file with Nhets that pass
    ###############
    module load $bcftoolsLoad
    echo 'AB ratio part 1'
    mkdir -p ${out}AB_hetPass
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi

    if "$sexChrom"; then
        #We only calculate AB ratio for XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n'  -S ${resources}${xx} \
        -i 'GT="het" & binom(FMT/AD) > 0.01' ${input}${infile} | \
        awk '{print $1"\t"NF -1}' > ${out}AB_hetPass/hetPass_${i}_XX
    else
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n' \
        -i 'GT="het" & binom(FMT/AD) > 0.01' ${input}${infile} | \
        awk '{print $1"\t"NF -1}' > ${out}AB_hetPass/hetPass_${i}
    fi
}

ABRatioP2(){
    ###############
    ##  Purpose: Number of het GTs for p2 AB ratio
    ##  Input: BCF
    ##  Output: txt file with ID and N hets
    ###############
    module load $bcftoolsLoad
    echo 'AB ratio part 2'
    mkdir -p ${out}AB_hetAll
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi

    if "$sexChrom"; then
        #We only calculate AB ratio for XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n' -i 'GT="het"' ${input}${infile} -S ${resources}${xx} | \
        awk '{print $1"\t"NF-1}' > ${out}AB_hetAll/hetAll_${i}_XX
    else
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n' -i 'GT="het"' ${input}${infile}  | \
        awk '{print $1"\t"NF-1}' > ${out}AB_hetAll/hetAll_${i}
    fi
}

pull_1KGPSites(){
    ###############
    ##  Purpose: Pull sites from 1000KGP3
    ##  Input: 1000KGP3 genotypes vcf
    ##  Output: Txt file
    ###############
    mkdir -p ${out}/1KGP3_intersect/
    chr="${LSB_JOBINDEX}"
    kgpin=$(ls -d "/public_data_resources/1000-genomes/20130502_GRCh38/"*.vcf.gz | \
    grep ${chr}_ | grep -v genotypes )
    
    zcat ${kgpin} | awk '/^##/ {next} { print $1"\t"$2"\t"$4"\t"$5}'  > tmp_1kgp_${chr}.txt
}

aggregateAnnotation(){
    ###############
    ##  Purpose: Annotate and make pass/fail. If king set to T in env, print subset of cols
    ##  Input: All the outputs of step 1 metrics
    ##  Output: A single text file containing annotated variants
    ###############
    #TODO - add option for chrom X here
    module load $RLoad
    
    mkdir -p ${out}Annotation
    #Location for summary data
    mkdir -p ${out}Annotation/Summary_stats #loc hardcoded in R script
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$7) {print $5"_"$6"_"$7}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
    
    export R_LIBS="$(dirname `which R`)/../lib64/R/library" #Set me to mitigate library mismatch on helix
    #Remember to set path var for .libpath
    Rscript annotatePerChunk.R \
    ${i} \
    ${out}startfile/start_file_${i} \
    ${out}missing/missing1_${i} \
    ${out}missing/missing2_${i} \
    ${out}medianCovAll/medianCov_${i} \
    ${out}medianCovNonMiss/medianCovNonMiss_${i} \
    ${out}medianGQ/medianGQ_${i} \
    ${out}AB_hetAll/hetAll_${i} \
    ${out}AB_hetPass/hetPass_${i} \
    ${out}MendelErrSites/MendErr_${i}.lmendel \
    ${out}Annotation \
    ${resources}/N_samples \
    ${out}AC_counts/${i}_AC \
    tmp_1kgp_${chr}.txt
}


### KING WORKFLOW ###

sort_compress(){
    ###############
    ##  Purpose: Sort and compress site metric data for KING step
    ##  Input: Text file containing annotated variants
    ##  Output: bgzipped and tabixed metrics information
    ###############

    autosomes=`grep -v chrX ${bedfile} | grep -v chrY | grep -v chrX | grep -v chrM`
    export infile=`echo $autosomes | cut -d' ' -f ${LSB_JOBINDEX}`
    #Only autosomes
    export i=`echo $infile | awk -F "_" 'sub(/.bcf/,"",$7) {print $5"_"$6"_"$7}'`

    module load $bcftoolsLoad
    sort -k2 -n ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}.txt > \
    ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt
    bgzip -f ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt && \
    tabix -s1 -b2 -e2 ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt.gz
    rm ${out}Annotation/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt
}

regions_filter(){
    ###############
    ##  Purpose: Produce BCFs of our data filtered to sites pass sites
    ##  Input: BCF, site metrics KING sites
    ##  Output: filtered bcf
    ###############
    module load $bcftoolsLoad
    mkdir -p ${out}AnnotatedVCFs/regionsFiltered
    autosomes=`grep -v chrX ${bedfile} | grep -v chrY | grep -v chrX | grep -v chrM`
    export infile=`echo $autosomes | cut -d' ' -f ${LSB_JOBINDEX}`
    #Only autosomes
    export i=`echo $infile | awk -F "_" 'sub(/.bcf/,"",$7) {print $5"_"$6"_"$7}'`

    bcftools view ${input}${infile} \
    -T ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt.gz  \
    -Ob \
    -o ${out}AnnotatedVCFs/regionsFiltered/${i}_regionsFiltered.bcf
}

#Regions_filter and further_filtering split for speed - split bcftools runs into multiple ops

further_filtering(){
    ###############
    ##  Purpose: Second stage filtering to give biallelic SNPs intersected with 1000KGP3 with MAF > 0.01
    ##  Input: Filtered BCFs, Michigan high LD pos
    ##  Output: BCF file, and txt file with positions
    ###############
    #Previously used to filter out Michigan high LD sites in this func, this now occurs
    #downstream using plink
    module load $bcftoolsLoad
    autosomes=`grep -v chrX ${bedfile} | grep -v chrY | grep -v chrX | grep -v chrM`
    export infile=`echo $autosomes | cut -d' ' -f ${LSB_JOBINDEX}`
    #Only autosomes
    export i=`echo $infile | awk -F "_" 'sub(/.bcf/,"",$7) {print $5"_"$6"_"$7}'`

    bcftools view ${out}AnnotatedVCFs/regionsFiltered/${i}_regionsFiltered.bcf \
    -i 'INFO/OLD_MULTIALLELIC="." & INFO/OLD_CLUMPED="."' \
    -v snps  | \
    bcftools annotate \
    --set-id '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC' | \
    bcftools +fill-tags -Ob \
    -o ${out}AnnotatedVCFs/regionsFiltered/MichiganLD_regionsFiltered_${i}.bcf \
    -- -t MAF
    #Produce filtered txt file
    bcftools query ${out}AnnotatedVCFs/regionsFiltered/MichiganLD_regionsFiltered_${i}.bcf \
    -i 'MAF[0]>0.01' -f '%CHROM\t%POS\t%REF\t%ALT\t%MAF\n' | \
    awk -F "\t" '{ if(($5 == "G" && $6 == "C") || ($6 == "G" && $5 == "C") || ($5 == "A" && $6 == "T") || ($6 == "A" && $5 == "T")) {next} { print $0} }' \
    > ${out}AnnotatedVCFs/MAF_filtered_1kp3intersect_${i}.txt
}

final_KING_BCF(){
    ###############
    ##  Purpose: Produce new BCF just with filtered sites
    ##  Input: Txt file of samples in agg, filtered regions bcf, intersected sites txt file
    ##  Output: Filtered compressed and indexed vcf
    ###############
    mkdir -p ${out}/KING
    autosomes=`grep -v chrX ${bedfile} | grep -v chrY | grep -v chrX | grep -v chrM`
    export infile=`echo $autosomes | cut -d' ' -f ${LSB_JOBINDEX}`
    #Only autosomes
    export i=`echo $infile | awk -F "_" 'sub(/.bcf/,"",$7) {print $5"_"$6"_"$7}'`
    #Now filter down our file to just samples we want in our GRM. This removes any withdrawals that we learned of during the process of aggregation
    #Store the header
    bcftools view \
    -S ${resources}${sampleList} \
    --force-samples \
    -h ${out}AnnotatedVCFs/regionsFiltered/MichiganLD_regionsFiltered_${i}.bcf \
    > ${out}KING/${i}_filtered.vcf
    
    #Then match against all variant cols in our subsetted bcf to our maf filtered, intersected sites and only print those that are in the variant file.
    #Then append this to the stored header, SNPRelate needs vcfs so leave as is
    bcftools view \
    -H ${out}AnnotatedVCFs/regionsFiltered/MichiganLD_regionsFiltered_${i}.bcf \
    -S ${resources}${sampleList} \
    --force-samples \
    | awk -F '\t' 'NR==FNR{c[$1$2$3$4]++;next}; c[$1$2$4$5] > 0' ${out}AnnotatedVCFs/MAF_filtered_1kp3intersect_${i}.txt - >> ${out}KING/${i}_filtered.vcf
    bgzip ${out}KING/${i}_filtered.vcf
    tabix ${out}KING/${i}_filtered.vcf.gz
}

concat_KING_VCF(){
    ###############
    ##  Purpose: Concatenate compressed vcfs to per chromosome files
    ##  Input: Compressed vcfs
    ##  Output: 1 per chrom compressed vcf
    ###############
    module load $bcftoolsLoad
    i="${LSB_JOBINDEX}"
    mkdir -p ${out}perChrom_KING

    find ${out}KING -type f -name "chr${i}_*.vcf.gz" > tmp.files_chrom${i}.txt
    bcftools concat \
    -f tmp.files_chrom${i}.txt \
    -Oz \
    -o ${out}perChrom_KING/chrom${i}_merged_filtered.vcf.gz && \
    tabix ${out}perChrom_KING/chrom${i}_merged_filtered.vcf.gz && \
    rm tmp.files_chrom${i}.txt
}

# POSSIBLE CHECK - can read file and check for ACTG combos that should have been filtered
makeBedAll(){
    ###############
    ##  Purpose: Make BED files for 1000KGP3 intersected vcfs
    ##  Input: per chrom vcf.gz
    ##  Output: BED files per chrom
    ###############
    module load $plink2Load
    module load $bcftoolsLoad
    echo 'Creating bed file'
    mkdir -p ${out}BEDref
    module load $plink2Load
    autosomes=`grep -v chrX ${bedfile} | grep -v chrY | grep -v chrX | grep -v chrM`
    export infile=`echo $autosomes | cut -d' ' -f ${LSB_JOBINDEX}`
    #Only autosomes
    export i=`echo $infile | awk -F "_" 'sub(/.bcf/,"",$7) {print $5"_"$6"_"$7}'`

    mkdir -p ${out}BED
    plink2 --vcf ${out}perChrom_KING/chrom${i}_merged_filtered.vcf.gz \
    --make-bed \
    --vcf-half-call m \
    --set-missing-var-ids chr@:#-\$r/\$a-.-. \
    --new-id-max-allele-len 60 missing \
    --exclude range ${resources}${michiganList} \
    --double-id \
    --real-ref-alleles \
    --allow-extra-chr \
    --threads 30 \
    --out ${out}BED/BED_${i}
}

LD_BED(){
    ###############
    ##  Purpose: LD prune SNPs
    ##  Input: Previoulsy produced BED, BIM, FAM files
    ##  Output: BED without high LD SNPs
    ###############
    #Not considering founders in this as all of our SNPs are common
    module load $plinkLoad
    plink  \
    --keep-allele-order \
    --bfile ${out}BED/BED_${i} \
    --indep-pairwise 500kb 1 0.1 \
    --threads 30 \
    --out ${out}BED/BED_LD_${i}
    
    #Now that we have our correct list of SNPs (prune.in), filter the original
    #bed file to just these sites
    plink \
    --make-bed \
    --bfile ${out}BED/BED_${i} \
    --keep-allele-order \
    --extract ${out}/BED/BED_LD_${i}.prune.in \
    --double-id \
    --allow-extra-chr \
    --threads 30 \
    --out ${out}BED/BED_LDpruned_${i}
}

merge_autosomes(){
    ###############
    ##  Purpose: Merge autosomes to genome wide BED files
    ##  Input: per chrom bed files (LD pruned)
    ##  Output: genome wide BED files
    ###############
    module load $plinkLoad
    for i in {1..22}; do
        echo ${out}BED/BED_LDpruned_$i >> mergelist.txt
    done
    plink --merge-list mergelist.txt \
    --make-bed \
    --out ${out}BED/autosomes_LD_pruned_1kgp3Intersect
    rm mergelist.txt
}

hwe_pruning_30k_snps()){
    ###############
    ##  Purpose: produce a first pass HWE filter. 
    ##  Input: LD pruned snps, bed files
    ##  Output: per pop unrelated bed files, list of hwe filtered snps
    ###############
    # Note the input file for unrels (third supplied argument to R) MUST have a column
    # 'unrelated_set' with either 1 or 0 denoting unrelated (1) or related (0)
    # Currently I have no provision to test for this
    module load $RLoad
    module load $plinkLoad
    R -e 'library(data.table);
    library(dplyr);
    args <- commandArgs(trailingOnly = T);
    setwd(args[1]);
    dat <- fread(args[2]) %>% as_tibble();
    unrels <- fread(args[3]) %>% as_tibble() %>% filter(unrelated_set == 1);
    dat <- dat %>% filter(plate_key %in% unrels$plate_key);
    for(col in c("AFR","EUR","SAS","EAS")){dat[dat[col]>0.8,c("plate_key",col)] %>% write.table(paste0(col,"pop.txt"), quote = F, row.names=F, sep = "\t")}' \
    --args \
    ${out}BED/ \
    ${resources}${snps30k} \
    ${resources}${pcs30k}

    bedmain="${out}/BED/autosomes_LD_pruned_1kgp3Intersect"
    for pop in AFR EUR SAS EAS; do
        echo ${pop}
        awk '{print $1"\t"$1}' ${pop}pop.txt > ${pop}keep
        plink \
        --keep ${pop}keep \
        --keep-allele-order \
        --make-bed \
        --bfile ${bedmain} \
        --out ${pop}
    
        plink --bfile ${pop} --hardy midp --out ${pop} --nonfounders
    done


    #Combine the HWE and produce a list of pass 
    R -e 'library(data.table);
    library(dplyr);
    args <- commandArgs(trailingOnly = T);
    setwd(args[1]);
    dat <- lapply(c("EUR.hwe","AFR.hwe", "SAS.hwe", "EAS.hwe"),fread);
    names(dat) <- c("EUR.hwe","AFR.hwe", "SAS.hwe", "EAS.hwe");
    dat <- dat %>% bind_rows(.id="id");
    write.table(dat, "combinedHWE.txt", row.names = F, quote = F)
    #Create set that is just SNPS that are >1e-5 in all pops
    dat %>% filter(P >1e-5) %>% group_by(SNP) %>% count() %>% filter(n==4) %>% select(SNP) %>% distinct() %>%
    write.table("hwe1e-5_superpops_195ksnps", sep="\t", row.names = F, quote = F)
    ' --args ${out}BED/ 
    #If we wanted to give the user options to change the MAF cutoff, we would need to expose this option in the filter
}

king_coefficients(){
    module load $plink2Load
    mkdir -p ${out}KING/matrix
    plink2 \
    --bfile ${out}BED/autosomes_LD_pruned_1kgp3Intersect \
    --extract ${out}BED/hwe10e-6_superpops_195ksnps \
    --make-king triangle bin \
    --out ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_triangle_HWE1e_5 \
    --thread-num 30

    #Filter the file
    plink2 --bfile ${out}BED/autosomes_LD_pruned_1kgp3Intersect_triangle_HWE1e_5 \
    --make-bed \
    --keep ${out}KING/matrix/plink2.king.cutoff.in.id \
    --out ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_unrelated

    #Also produce a related set
    #Filter the file
    plink2 --bfile ${out}BED/autosomes_LD_pruned_1kgp3Intersect_triangle_HWE1e_5 \
    --make-bed \
    --remove ${out}KING/matrix/plink2.king.cutoff.in.id \
    --out ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_related
}

##
# Run Thanos PCs script
##

##
# Run the Ancestry inference script
##
infer_ancestry(){
    module load $RLoad
    mkdir -p ${out}Ancestries
    export R_LIBS="$(dirname `which R`)/../lib64/R/library" #Set me to mitigate library mismatch on helix
    Rscript infer_ancestry.R  \
    ${resources}${1kgp3_sample_table} \
    ${resources}${super_pop_codes} \
    ${resources}${1kgp3_unrel} \
    #Need to think about how these will actually link up
    "/re_gecip/BRS/thanos/ethnicities_aggV2/aggV2_ethn_200Ksites/1KGP3_intersectGEL_200Kset_perpopHWE1e-6_unrel_maf0.05both1K100K.eigenvec" \
    "/re_gecip/BRS/thanos/ethnicities_aggV2/aggV2_ethn_200Ksites/1KGP3_intersectGEL_200Kset_perpopHWE1e-6_unrel_maf0.05both1K100K_GELprojection.proj.eigenvec" \
    ${out}Ancestries
}

## END OF ANCESTRY INFERENCE/KING WORKFLOW ##


p_hwe(){
    module load $bcftoolsLoad
    module load $plinkLoad
    mkdir -p ${out}HWE
    #Take the unrelated set, split by ancestry
    ####
    # What is needed for this?
    # Check the dat and unrel vars in the R script below. the dat is inferred ancestries per sample.
    # Unrels is given by Thanos
    # Once these are changed the script will run fine as is
    # You can run this interactively as computational overheads are v low
    ####
    R -e 'library(data.table)
    library(dplyr)
         dat <- fread("/re_gecip/shared_allGeCIPs/drhodes/covid_siteqc/testing/output/ancestry/maf5_ancestries.txt") %>% as_tibble()
         unrel <- fread("/re_gecip/BRS/genomicc_data/genomicc_WGS/relatedness/aggCOVID1_N2008_HQSNPs_MAF5_65K.king.cutoff.in.id") %>% as_tibble()
        dat <- dat %>% filter(mafList.pcsgel.Sample %in% unrel$IID)
        for(col in c("AFR","EUR","SAS","EAS")){dat[dat[col]>0.8,c("mafList.pcsgel.Sample",col)] %>%
        write.table(paste0("/re_gecip/shared_allGeCIPs/drhodes/covid_siteqc/testing/output/ancestry/",col,"_unrelated_pop.txt"), quote = F, row.names=F, sep = "\t")};'


        for pop in AFR EUR SAS EAS; do
            
            echo -e "Calculating HWE on ${pop} samples..."
            plink --bfile /re_gecip/BRS/genomicc_data/genomicc_WGS/HG_SNPs/aggCOVID1_N2008_HQSNPs_MAF5_65K \
            --hardy midp \
            --keep ${out}/ancestry/${pop}_unrelated_pop.keep \
            --double-id \
            --allow-extra-chr \
            --out ${out}/HWE/${i}_$pop
        done

    #Combine the HWE and produce a list of pass 
    cd ${wd}../outputs/HWE
    R -e 'library(data.table);
    library(dplyr);
    dat <- lapply(c(list.files(pattern = "*.hwe")),fread);
    names(dat) <- c("EUR.hwe","AFR.hwe", "SAS.hwe", "EAS.hwe");
    dat <- dat %>% bind_rows(.id="id");
    write.table(dat, "combinedHWE.txt", row.names = F, quote = F)'
    R -e 'library(dplyr); library(data.table);
    dat <- fread("combinedHWE.txt") %>% as_tibble()
     #Create set that is just SNPS that are >1e-5 in all pops
    dat %>% filter(P >1e-5) %>% group_by(SNP) %>% count() %>% filter(n==4) %>% select(SNP) %>% distinct() %>%
    write.table("hwe1e-5_superpops_high_conf_snps", sep="\t", row.names = F, quote = F)
    '
}


#Export functions

#Step 1 metric calculations
export -f startFile
export -f missingness1
export -f missingness2
export -f medianCovAll
export -f ABRatioP1
export -f ABRatioP2
export -f completeSites
export -f medianGQ
export -f medianCovNonMiss
export -f aggregateAnnotation
export -f pull_1KGPSites

#King workflow



#Step 2 metric calculations

#Aggregation




