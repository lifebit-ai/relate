#!/bin/bash

# @Author: daniel
# @Date:   2019-11-08T10:17:07+00:00
# @Email:  daniel.rhodes@genomicsengland.co.uk
# @Project: Aggregate_QC
# @Last modified by:   daniel
# @Last modified time: 2019-11-12T09:08:32+00:00


# All within based on code by MW - see run.sh for details on script origins

#Modules, specified according to which cluster it is being run on
if lsid | grep -q pegasus; then
    echo 'Loading modules for Pegasus HPC'
    export vcftoolsLoad='vcftools/0.1.15'
    export bcftoolsLoad='bcftools/1.10.2'
    export plinkLoad='PLINK/1.90'
    export RLoad='R/3.6.0'
else
    echo 'Loading modules for Helix HPC'
    export bcftoolsLoad='bio/BCFtools/1.10.2-foss-2018b'
    export RLoad='lang/R/3.6.0-foss-2019a'
    export plinkLoad='bio/PLINK/1.9b_4.1-x86_64'
    export plink2Load='bio/PLINK/2.00-devel-20200409-x86_64'
    export kingLoad='bio/KING/2.2.4'
    
fi

#TODO
#Need to add the new list from Chris of samples to be included
#This needs to be used for final output, and for creation of bed files
#### Define funcs ####
triodata_define(){  
    ###############
    ##  Purpose: Create a .fam file of confirmed trios in our samples
    ##  Input: List of trios, samples included in agg
    ##  Output: .fam and .keep files
    ###############
    module load $RLoad
    Rscript trio_define.R $triodata $aggregateSamples $triodata.fam $triodata.keep
}

startFile(){
    ###############
    ##  Purpose: Create a backbone of IDs for other data to be joined to
    ##  Input: BCF
    ##  Output: txt file with IDs
    ###############
    module load $bcftoolsLoad
    echo 'Creating backbone file'
    # Query chrom and pos for annotation file
    mkdir -p ${out}startfile
    if "$sexChrom"; then
        #XY
        bcftools query -f '%CHROM %POS %REF %ALT %INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC\n' ${infile} -S ${resources}${xy} --output ${out}startfile/start_file_${i}_XY
        #XX
        bcftools query -f '%CHROM %POS %REF %ALT %INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC\n' ${infile} -S ${resources}${xx} --output ${out}startfile/start_file_${i}_XX
    else
        #Autosomes
        bcftools query -f '%CHROM %POS %REF %ALT %INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC\n' ${infile} --output ${out}startfile/start_file_${i}
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
    mkdir -p ${out}missing2
    
    if "$sexChrom"; then
        #XY
        bcftools query $infile -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -i 'GT="./." & FORMAT/DP=0'  -S ${resources}${xy} \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing2/missing1_${i}_XY
        #XX
        bcftools query $infile -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -i 'GT="./." & FORMAT/DP=0'  -S ${resources}${xx} \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing2/missing1_${i}_XX
    else
        #autosomes
        bcftools query $infile -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -i 'GT="./." & FORMAT/DP=0' \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing2/missing1_${i}
    fi
}

missingness2(){
    ###############
    ##  Purpose: Count complete GTs only
    ##  Input: BCF
    ##  Output: Txt file with ID and count
    ###############
    module load $bcftoolsLoad
    echo 'Calculating missing sites (complete only)'
    mkdir -p ${out}missing2
    if "$sexChrom"; then
        #XY
        bcftools query $infile -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -e 'GT~"\."' -S ${resources}${xy} \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing2/missing2_${i}_XY
        #XX
        bcftools query $infile -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -e 'GT~"\."' -S ${resources}${xx} \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing2/missing2_${i}_XX
    else
        #Autosomes
        bcftools query $infile -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -e 'GT~"\."' \
        | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing2/missing2_${i}
    fi
}

completeSites(){
    ###############
    ##  Purpose: Make sure the number of samples is listed in resources
    ##  Input: BCF
    ##  Output: Txt file with N samples
    ###############
    module load $bcftoolsLoad
    #All this will do is create a file containing the number of samples to be used
    #by the annotation script
    if [ ! -f "${resources}/N_samples" ]; then
        bsub -q short -P bio -e logs/n_samples_err%J -o logs/n_samples_out%J 'module load $bcftoolsLoad; bcftools query -l ${infile} | wc -l > ${resources}/N_samples'
    fi
}

completeSitesSexChrom(){
    ###############
    ##  Purpose: Same as completeSites, but splitting by sex
    ##  Input: BCF
    ##  Output: XX and XY number of samples
    ###############
    #All this will do is create a file containing the number of samples to be used
    #by the annotation script
    if [ ! -f "${resources}/N_samplesXX" ]; then
        bsub -q short -P bio -e logs/n_samplesXX_err%J -o logs/n_samplesXX_out%J 'module load $bcftoolsLoad; bcftools query -l ${infile} -S ${resources}/xx_females_illumina_ploidy_samples_40740.tsv | wc -l > ${resources}/N_samplesXX'
    fi
    if [ ! -f "${resources}/N_samplesXY" ]; then
        bsub -q short -P bio -e logs/n_samplesXY_err%J -o logs/n_samplesXY_out%J 'module load $bcftoolsLoad; bcftools query -l ${infile} -S ${resources}/xy_males_illumina_ploidy_samples_35924.tsv | wc -l > ${resources}/N_samplesXY'
    fi
}

medianCoverageAll(){
    ###############
    ##  Purpose: Produce median value for depth across all GT
    ##  Input: BCF
    ##  Output: Txt file with ID and median depth
    ###############
    module load $bcftoolsLoad
    mkdir -p ${out}medianCoverageAll
    if "$sexChrom"; then
        #XY
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' $infile -S ${resources}${xy} |\
        awk '{sum=0; n=split($0,a)-1; for(i=2;i<=n;i++) sum+=a[i]; asort(a);
                median=n%2?a[n/2+1]:(a[n/2]+a[n/2+1])/2;
        print $1"\t"median}' > ${out}medianCoverageAll/medianCoverageAll${i}_XY
        #XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' $infile -S ${resources}${xx} |\
        awk '{sum=0; n=split($0,a)-1; for(i=2;i<=n;i++) sum+=a[i]; asort(a);
                median=n%2?a[n/2+1]:(a[n/2]+a[n/2+1])/2;
        print $1"\t"median}' > ${out}medianCoverageAll/medianCoverageAll${i}_XX
    else
        #autosomes
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' $infile |\
        awk '{sum=0; n=split($0,a)-1; for(i=2;i<=n;i++) sum+=a[i]; asort(a);
                median=n%2?a[n/2+1]:(a[n/2]+a[n/2+1])/2;
        print $1"\t"median}' > ${out}medianCoverageAll/medianCoverageAll${i}
    fi
}

medianCoverageNonMiss(){
    ###############
    ##  Purpose: Median coverage for fully present genotypes
    ##  Input: BCF
    ##  Output: Txt file with ID and median
    ###############
    module load $bcftoolsLoad
    #Allele specific, coverage of all non-missing genotypes
    #The NF-1 in sum/=NF-1 accounts for us treating the CHROM:POS-REF/ALT as a value of 0
    mkdir -p ${out}medianCoverageNonMiss
    if "$sexChrom"; then
        #XY
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' -e 'GT~"\."' $infile -S ${resources}${xy} | \
        awk '{sum=0; n=split($0,a)-1; for(i=2;i<=n;i++) sum+=a[i]; asort(a);
            median=n%2?a[n/2+1]:(a[n/2]+a[n/2+1])/2;
        print $1"\t"median}' > \
        ${out}medianCoverageNonMiss/medianNonMiss_depth_${i}_XY
        #XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' -e 'GT~"\."' $infile -S ${resources}${xx} | \
        awk '{sum=0; n=split($0,a)-1; for(i=2;i<=n;i++) sum+=a[i]; asort(a);
            median=n%2?a[n/2+1]:(a[n/2]+a[n/2+1])/2;
        print $1"\t"median}' > \
        ${out}medianCoverageNonMiss/medianNonMiss_depth_${i}_XX
    else
        #autosomes
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' -e 'GT~"\."' $infile | \
        awk '{sum=0; n=split($0,a)-1; for(i=2;i<=n;i++) sum+=a[i]; asort(a);
            median=n%2?a[n/2+1]:(a[n/2]+a[n/2+1])/2;
        print $1"\t"median}' > \
        ${out}medianCoverageNonMiss/medianNonMiss_depth_${i}
    fi
}

medianGQ(){
    ###############
    ##  Purpose: Calculate median GQ
    ##  Input: BCF
    ##  Output: txt file - ID and median
    ###############
    module load $bcftoolsLoad
    mkdir -p ${out}medianGQ
    if "$sexChrom"; then
        #XY
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GQ ] \n' -e 'GT~"\."' $infile -S ${resources}${xy} |\
        awk '{sum=0; n=split($0,a)-1; for(i=2;i<=n;i++) sum+=a[i]; asort(a);
                median=n%2?a[n/2+1]:(a[n/2]+a[n/2+1])/2; if(median > 99) {median=99};
        print $1"\t"median}' > ${out}medianGQ/medianGQ_${i}_XY
        #XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GQ ] \n' -e 'GT~"\."' $infile -S ${resources}${xx} |\
        awk '{sum=0; n=split($0,a)-1; for(i=2;i<=n;i++) sum+=a[i]; asort(a);
                median=n%2?a[n/2+1]:(a[n/2]+a[n/2+1])/2; if(median > 99) {median=99};
        print $1"\t"median}' > ${out}medianGQ/medianGQ_${i}_XX
    else
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GQ ] \n' -e 'GT~"\."' $infile |\
        awk '{sum=0; n=split($0,a)-1; for(i=2;i<=n;i++) sum+=a[i]; asort(a);
                median=n%2?a[n/2+1]:(a[n/2]+a[n/2+1])/2; if(median > 99) {median=99};
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
    
    if "$sexChrom"; then
        #We only calculate AB ratio for XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n'  -S ${resources}${xx} \
        -i 'GT="het" & binom(FMT/AD) > 0.01' ${infile} | \
        awk '{print $1"\t"NF -1}' > ${out}AB_hetPass/hetPass_${i}_XX
    else
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n' \
        -i 'GT="het" & binom(FMT/AD) > 0.01' ${infile} | \
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
    if "$sexChrom"; then
        #We only calculate AB ratio for XX
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n' -i 'GT="het"' ${infile} -S ${resources}${xx} | \
        awk '{print $1"\t"NF-1}' > ${out}AB_hetAll/hetAll_${i}_XX
    else
        bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n' -i 'GT="het"' ${infile}  | \
        awk '{print $1"\t"NF-1}' > ${out}AB_hetAll/hetAll_${i}
    fi
}

MendErrp1(){
    ###############
    ##  Purpose: Create a bed file for the mendel error calcs
    ##  Input: Input BCF
    ##  Output: BED related files
    ###############
    module load $bcftoolsLoad
    module load $plink2Load
    mkdir -p ${out}MendelErr
    #Need to generate binary files based on just the subset of samples we want
    plink2 --bcf ${infile} \
    --keep $triodata.keep \
    --make-bed \
    --vcf-half-call m \
    --set-missing-var-ids @:#,\$r,\$a \
    --new-id-max-allele-len 60 missing\
    --double-id \
    --real-ref-alleles \
    --allow-extra-chr \
    --out ${out}MendelErr/BED_trio_${i}
}

pullAC(){
    ###############
    ##  Purpose: Pull AC from all files and store for addition to site metrics
    ##  Input: BCF
    ##  Output: Txt file with AC
    ###############
    module load $bcftoolsLoad
    mkdir -p ${out}AC_counts
    bcftools query -f \
    '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC %INFO/AC \n' \
    ${infile}  > ${out}AC_counts/${i}_AC
}

pull_1KGPSites(){
    ###############
    ##  Purpose: Pull sites from 1000KGP3
    ##  Input: 1000KGP3 genotypes vcf
    ##  Output: Txt file
    ###############
    chr="${i%%_*}"
    kgpin=$(ls -d "/public_data_resources/1000-genomes/20130502_GRCh38/"*.vcf.gz | \
    grep ${chr}_ | grep -v genotypes )
    
    zcat ${kgpin} | awk '/^##/ {next} { print $1"\t"$2"\t"$4"\t"$5}'  > tmp_1kgp${i}.txt
}

MendErrp2(){
    ###############
    ##  Purpose: Calculate mendelian errors
    ##  Input: BED related files, our own .fam file
    ##  Output: lmendel, fmendel, imendel
    ###############
    module load $plinkLoad
    #Calc mendel errors using these
    #Half calls are set to missing
    #Note we use a predefined .fam file
    plink --bed ${out}MendelErr/BED_trio_${i}.bed \
    --bim ${out}MendelErr/BED_trio_${i}.bim \
    --fam $triodata.fam \
    --allow-extra-chr \
    --allow-no-sex \
    --mendel summaries-only \
    --out ${out}MendelErr/MendErr_${i}
}

MendelDist(){
    ###############
    ##  Purpose: Summary stats and good families for Mendel errors
    ##  Input: .fmendel files for all families
    ##  Output: summary stats, MendelFamilies_4SD.csv for good families
    ###############
    module load $RLoad
    #Plot the family wise mendel dists
    Rscript mendelerror_family_plotting.R ${out}MendelErr
    #Following this, it is necessary to decide what the cut off will be.
    #Then produce a list of FIDs that will be KEPT. This is passed to plink for
    #mendel error round 2.
}

MendErrp3(){
    ###############
    ##  Purpose: Calculate mendel errors on just good families
    ##  Input: BED files, list of good families (.fam file)
    ##  Output: mendel error summaries
    ###############
    module load $plinkLoad
    mkdir -p ${out}MendelErrSites
    #Calc mendel errors using these
    #Half calls are set to missing
    plink --bed ${out}MendelErr/BED_trio_${i}.bed \
    --bim ${out}MendelErr/BED_trio_${i}.bim \
    --fam $triodata.fam \
    --allow-extra-chr \
    --allow-no-sex \
    --keep-fam ${out}MendelErr/MendelFamilies_4SD.fam \
    --mendel summaries-only --out ${out}MendelErrSites/MendErr_${i}
    #Check what we need to do with the families file for this
}

#Just before this step is where i want a checklist, that all chunks are completed

aggregateAnnotation(){
    ###############
    ##  Purpose: Annotate and make pass/fail. If king set to T in env, print subset of cols
    ##  Input: All the outputs of step 1 metrics
    ##  Output: A single text file containing annotated variants
    ###############
    module load $RLoad
    #Some of the below will need to be changed when we finalise filepaths etc for data
    mkdir -p ${out}Annotation_newtest
    Rscript annotatePerChunk.R \
    ${i} \
    ${out}startfile/start_file_${i} \
    ${out}missing2/missing1_${i} \
    ${out}missing2/missing2_${i} \
    ${out}medianCoverageAll/medianCoverageAll${i} \
    ${out}medianCoverageNonMiss/medianNonMiss_depth_${i} \
    ${out}medianGQ/medianGQ_${i} \
    ${out}AB_hetAll/hetAll_${i} \
    ${out}AB_hetPass/hetPass_${i} \
    ${out}MendelErrSites/MendErr_${i}.lmendel \
    ${out}Annotation_newtest \
    ${resources}/N_samples \
    ${out}AC_counts/${i}_AC \
    tmp_1kgp${i}.txt #this will be deleted in next step
}

sexChromAnnotation(){
    ###############
    ##  Purpose: Annotate and make pass/fail. If king set to T in env, print subset of cols
    ##  Input: All the outputs of step 1 metrics
    ##  Output: A single text file containing annotated variants
    ###############
    module load $RLoad
    mkdir -p ${out}Annotation_newtest/sexQC/

    mkdir -p ${out}Annotation_newtest
    Rscript sexChromAnnotation.R \
    ${i} \
    ${out}startfile/start_file_${i} \
    ${out}missing2/missing1_${i} \
    ${out}missing2/missing2_${i} \
    ${out}medianCoverageAll/medianCoverageAll${i} \
    ${out}medianCoverageNonMiss/medianNonMiss_depth_${i} \
    ${out}medianGQ/medianGQ_${i} \
    ${out}AB_hetAll/hetAll_${i} \
    ${out}AB_hetPass/hetPass_${i} \
    ignore \
    ${out}Annotation_newtest/sexQC/ \
    ${resources}/N_samples
}

check_annotation_length(){
    #Check the input files to ensure they are the same length
    #Have to store as variables due to redirect issues
    touch checking_annotations.txt
    backbone=`wc -l ${out}startfile/start_file_${i}`
    files="${out}missing2/missing1_${i}
    ${out}missing2/missing2_${i}
    ${out}medianCoverageAll/medianCoverageAll${i}
    ${out}medianCoverageNonMiss/medianNonMiss_depth_${i}
    ${out}medianGQ/medianGQ_${i}
    ${out}AB_hetAll/hetAll_${i}
    ${out}MendelErrSites/MendErr_${i}.lmendel
    ${out}AC_counts/${i}_AC
    ${out}HWE/${i}_AFR.hwe
    ${out}HWE/${i}_EUR.hwe
    ${out}HWE/${i}_EAS.hwe
    ${out}HWE/${i}_SAS.hwe"
    for f in $files; do
        #Check the file exists, if not note this
        if [ ! -f $f ]; then
            echo -e "$f \t not found \t ${i}" >> checking_annotations.txt
        fi
        
        #Check file lengths
        f2=`wc -l $f`
        if [ "$f2" != "$backbone" ]; then
            echo -e "${f2} \t ${f} \t line length \t ${backbone} \t${i}" >> checking_annotations.txt
        fi
    done
}

sort_compress(){
    ###############
    ##  Purpose: Sort and compress site metric data for KING step
    ##  Input: Text file containing annotated variants
    ##  Output: bgzipped and tabixed metrics information
    ###############
    #Necessary step for makegVCF
    module load $bcftoolsLoad
    sort -k2 -n ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}.txt > \
    ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt
    bgzip -f ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt && \
    tabix -s1 -b2 -e2 ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt.gz
    #rm ${out}Annotation/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt &&
    #rm tmp_1kgp${i}.txt
}

duplicate_check(){
    ###############
    ##  Purpose: Checks for presence of MNPs by checking duplicate CHROM POS REF ALT
    ##  Input: startfile backbone
    ##  Output: TSV with duplicates
    ###############
    #NOTE:
    #This will not need to be run as part of the main pipeline. This is an ancillary function.
    #
    module load $RLoad
    #We know there are duplicate CHROM POS REF ALTs due to MNPs, let's double check exactly what the state of play is
    #Take all startfiles and AC counts, print out all the duplicate lines to one file
    tee duplication_check.R <<EOF
        library(data.table); library(dplyr); library(magrittr);
        files <- list.files('/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_actual/startfile', full.names = T)
        files <- files[!grepl('_X.$|_X.$', files)]
        output <- lapply(files, function(x){
            pos <- gsub(".*_chr","chr",x)
            f <- fread(x) %>% as_tibble()
            f <- f %>% group_by(V1, V2, V3, V4) %>% filter(n() > 1) %>% ungroup()
            f %<>% mutate(id = paste0(V1,':', V2, '-', V3,'/',V4,'-', V5))
            ac <- fread(paste0('/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_actual/AC_counts/',pos,'_AC')) %>% as_tibble()
            outter <- left_join(f, ac, by=c('id'='V1'))
            names(outter) <- c('CHROM','POS','REF','ALT','oldMulti_oldClump','uniqID','AC')
            return(outter)
        })
        out <- bind_rows(output)
        out %>% fwrite('MNP_variants.tsv', sep = '\t')
EOF
chmod +x duplication_check.R
bsub -q long -P bio -U as1 -e $HOME/duplicate.err -o $HOME/duplicate.out 'Rscript duplication_check.R'
}

### KING WORKFLOW ###
regions_filter(){
    ###############
    ##  Purpose: Produce BCFs of our data filtered to sites pass sites
    ##  Input: BCF, site metrics KING sites
    ##  Output: filtered bcf
    ###############
    module load $bcftoolsLoad
    mkdir -p ${out}AnnotatedVCFs/regionsFiltered
    bcftools view ${infile} \
    -T ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt.gz  \
    -Ob \
    -o ${out}AnnotatedVCFs/regionsFiltered/${i}_regionsFiltered.bcf
}

further_filtering(){
    ###############
    ##  Purpose: Second stage filtering to give biallelic SNPs intersected with 1000KGP3 with MAF > 0.01
    ##  Input: Filtered BCFs, Michigan high LD pos
    ##  Output: BCF file, and txt file with positions
    ###############
    module load $bcftoolsLoad
    #Note that the michigan sites filtering here is incorrect. Doing a workaround downstream
    bcftools view ${out}AnnotatedVCFs/regionsFiltered/${i}_regionsFiltered.bcf \
    -i 'INFO/OLD_MULTIALLELIC="." & INFO/OLD_CLUMPED="."' \
    -T ^${resources}/MichiganLD_liftover_exclude_regions.txt \
    -v snps  | \
    bcftools annotate \
    --set-id '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC' | \
    bcftools +fill-tags -Ob \
    -o ${out}AnnotatedVCFs/regionsFiltered/MichiganLD_regionsFiltered_${i}.bcf \
    -- -t MAF
    #Produce filtered txt file
    bcftools query ${out}AnnotatedVCFs/regionsFiltered/MichiganLD_regionsFiltered_${i}.bcf \
    -i 'MAF[0]>0.01' -f '%CHROM\t%POS\t%REF\t%ALT\t%MAF\n' | \
    awk -F "\t" '{ if(($3 == "G" && $4 == "C") || ($3 == "A" && $4 == "T")) {next} { print $0} }' \
    > ${out}AnnotatedVCFs/MAF_filtered_1kp3intersect_${i}.txt
}

final_KING_BCF(){
    ###############
    ##  Purpose: Produce new BCF just with filtered sites
    ##  Input: Txt file of samples in agg, filtered regions bcf, intersected sites txt file
    ##  Output: Filtered compressed and indexed vcf
    ###############
    mkdir -p ${out}/KING
    #Now filter down our file to just samples we want in our GRM. This removes any withdrawals that we learned of during the process of aggregation
    #Store the header
    bcftools view \
    -S ${resources}78389_final_platekeys_agg_v9.txt \
    --force-samples \
    -h ${out}AnnotatedVCFs/regionsFiltered/MichiganLD_regionsFiltered_${i}.bcf \
    > ${out}KING/${i}_filtered.vcf
    
    #Then match against all variant cols in our subsetted bcf to our maf filtered, intersected sites and only print those that are in the variant file.
    #Then append this to the stored header, SNPRelate needs vcfs so leave as is
    bcftools view \
    -H ${out}AnnotatedVCFs/regionsFiltered/MichiganLD_regionsFiltered_${i}.bcf \
    -S ${resources}78389_final_platekeys_agg_v9.txt \
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
    mkdir -p ${out}perChrom_KING
    find ${out}KING -type f -name "chr${i}_*.vcf.gz" > tmp.files_chrom${i}.txt
    bcftools concat \
    -f tmp.files_chrom${i}.txt \
    -Oz \
    -o ${out}perChrom_KING/chrom${i}_merged_filtered.vcf.gz && \
    tabix ${out}perChrom_KING/chrom${i}_merged_filtered.vcf.gz && \
    rm tmp.files_chrom${i}.txt
}

makeBedAll(){
    ###############
    ##  Purpose: Make BED files for 1000KGP3 intersected vcfs
    ##  Input: per chrom vcf.gz
    ##  Output: BED files per chrom
    ###############
    module load $plinkLoad
    module load $bcftoolsLoad
    echo 'Creating bed file'
    mkdir -p ${out}BEDref
    
    bcftools view ${out}perChrom_KING/chrom${i}_merged_filtered.vcf.gz \
    -Ov |\
    plink --vcf /dev/stdin \
    --vcf-half-call m \
    --double-id \
    --make-bed \
    --real-ref-alleles \
    --allow-extra-chr \
    --out ${out}BEDref/BED_${i}
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
    --exclude range ${resources}/MichiganLD_liftover_exclude_regions_PARSED.txt \
    --keep-allele-order \
    --bfile ${out}BED/BED_${i} \
    --indep-pairwise 500kb 1 0.1 \
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

hwe_pruning_30k_data(){
    #Aim of this function is to produce a first pass HWE filter. 
    #We use:
    #The 195k SNPs from above
    #The intersection bfiles (on all 80k)
    #Then we make BED files of unrelated individuals for each superpop (using only unrelated samples from 30k)
    #We do this using the inferred ancestries from the 30k

module load lang/R/3.6.0-foss-2019a
#
R -e 'library(data.table);
 library(dplyr);
  dat <- fread("aggV2_R9_M30K_1KGP3_ancestry_assignment_probs.tsv") %>% as_tibble();
  unrels <- fread("/re_gecip/BRS/thanos/aggV2_bedmerge_30KSNPs_labkeyV9_08062020_update_PCsancestryrelated.tsv") %>% as_tibble() %>% filter(unrelated_set == 1);
  dat <- dat %>% filter(plate_key %in% unrels$plate_key);
for(col in c("AFR","EUR","SAS","EAS")){dat[dat[col]>0.8,c("plate_key",col)] %>% write.table(paste0(col,"pop.txt"), quote = F, row.names=F, sep = "\t")}'

#Run on both 200k and 30k SNP set
module load bio/PLINK/1.9b_4.1-x86_64
#Need to do in unrelateds
#So find the 30k unrelated set

bedmain="${out}/BED/autosomes_LD_pruned_1kgp3Intersect"
for pop in AFR EUR SAS EAS; do
    echo ${pop}
    awk '{print $1"\t"$1}' ${pop}pop.txt > ${pop}keep
    plink \
    --keep ${pop}keep \
    --make-bed \
    --bfile ${bedmain} \
    --out ${pop}
    
    plink --bfile ${pop} --hardy --out ${pop} --nonfounders
done

#Combine the HWE and produce a list of pass 
R -e 'library(data.table);
library(dplyr);
dat <- lapply(c("EUR.hwe","AFR.hwe", "SAS.hwe", "EAS.hwe"),fread);
names(dat) <- c("EUR.hwe","AFR.hwe", "SAS.hwe", "EAS.hwe");
dat <- dat %>% bind_rows(.id="id");
write.table(dat, "combinedHWE.txt", row.names = F, quote = F)'
R -e 'library(dplyr); library(data.table);
    dat <- fread("combinedHWE.txt") %>% as_tibble()
     #Create set that is just SNPS that are >10e-6 in all pops
    dat %>% filter(P >10e-6) %>% group_by(SNP) %>% count() %>% filter(n==4) %>% select(SNP) %>% distinct() %>%
    write.table("hwe10e-6_superpops_195ksnps", sep="\t", row.names = F, quote = F)
    '
R -e 'library(dplyr); library(data.table);
    dat <- fread("combinedHWE.txt") %>% as_tibble()
     #Create set that is just SNPS that are >10e-2 in all pops
    dat %>% filter(P >10e-2) %>% group_by(SNP) %>% count() %>% filter(n==4) %>% select(SNP) %>% distinct() %>%
    write.table("hwe10e-2_superpops_195ksnps", sep="\t", row.names = F, quote = F)
    '
}

king_coefficients(){
    module load $plink2Load
    mkdir -p ${out}KING/matrix
    plink2 --bfile \
    ${out}BED/autosomes_LD_pruned_1kgp3Intersect \
    --make-king square \
    --out \
    ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect \
    --thread-num 30
}
#writing two customs just to do the above, but with the variable cut-off SNP list
#bsub -q long -P bio -U as1 -e $HOME/30k10-2.err 'module load $plink2Load; plink2 --bfile ${out}BED/autosomes_LD_pruned_1kgp3Intersect --extract /re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_actual/hwe_prune/hwe10e-2_superpops_195ksnps --make-king triangle bin --out ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_triangle_HWE10_2'
#bsub -q long -P bio -U as1 -e $HOME/30k10-10.err 'module load $plink2Load; plink2 --bfile ${out}BED/autosomes_LD_pruned_1kgp3Intersect --extract /re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_actual/hwe_prune/hwe10e-6_superpops_195ksnps --make-king triangle bin --out ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_triangle_HWE10_6'


king_coefficients_Alternate(){
    #The main difference for this is that we are aiming to do all the pcAIR
    #using other tools. Therefore the output needs to be different
    module load $plink2Load
    mkdir -p ${out}KING/matrix
    plink2 --bfile \
    ${out}BED/autosomes_LD_pruned_1kgp3Intersect \
    --make-king triangle bin \
    --out \
    ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_triangle \
    --thread-num 30
}
pcair_alternate(){
    #This isn't actually intended to run as a function, it is just to stop stuff running
    #when sourcing this file that we wrap it in a function
    #Alternate approach to producing the PC-relate info
    module load bio/FlashPCA2/2.0-foss_2018b
    module load bio/PLINK/1.9b_4.1-x86_64

    #Create partition using plink
    module load $plink2Load
    plink2 --bfile ${out}BED/autosomes_LD_pruned_1kgp3Intersect \
       --king-cutoff ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_triangle 0.0442 && \
       mv plink2.king.cutoff.in.id autosomes_LD_pruned_1kgp3Intersect.king.cutoff.in.id && \
       mv plink2.king.cutoff.out.id autosomes_LD_pruned_1kgp3Intersect.king.cutoff.out.id


    ##Partitions for the alternate triangles
    plink2 --bfile ${out}BED/autosomes_LD_pruned_1kgp3Intersect \
       --king-cutoff ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_triangle_HWE10_2 0.0442 && \
       mv plink2.king.cutoff.in.id  autosomes_LD_pruned_1kgp3Intersect_triangle_HWE10_2.king.cutoff.in.id && \
       mv plink2.king.cutoff.out.id  autosomes_LD_pruned_1kgp3Intersect_triangle_HWE10_2.king.cutoff.out.id
    plink2 --bfile ${out}BED/autosomes_LD_pruned_1kgp3Intersect \
       --king-cutoff ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_triangle_HWE10_6 0.0442 && \
       mv plink2.king.cutoff.in.id  autosomes_LD_pruned_1kgp3Intersect_triangle_HWE10_6.king.cutoff.in.id && \
       mv plink2.king.cutoff.out.id  autosomes_LD_pruned_1kgp3Intersect_triangle_HWE10_6.king.cutoff.out.id

    #Let's now have a look at how muhc overlaps in the kinship based on the HWE cutoffs
    R -e 'library(data.table); library(dplyr); library(magrittr);
        dat <- fread("autosomes_LD_pruned_1kgp3Intersect.king.cutoff.in.id") %>% as_tibble();
        hwe2 <- fread("autosomes_LD_pruned_1kgp3Intersect_triangle_HWE10_2.king.cutoff.in.id") %>% as_tibble();
        hwe6 <- fread("autosomes_LD_pruned_1kgp3Intersect_triangle_HWE10_6.king.cutoff.in.id") %>% as_tibble();
        dat <- bind_rows(dat, hwe2, hwe6, .id="id");
        dat %>% group_by(id) %>% summarise(n()); 
        dat %>% group_by(IID) %>% summarise(n=n()) %>% count(n)  '



    #Filter the file
    plink2 --bfile ${out}BED/autosomes_LD_pruned_1kgp3Intersect \
    --make-bed \
    --keep ${out}KING/matrix/plink2.king.cutoff.in.id \
    --out ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_unrelated

    #Also produce a related set
    #Filter the file
    plink2 --bfile ${out}BED/autosomes_LD_pruned_1kgp3Intersect \
    --make-bed \
    --remove ${out}KING/matrix/plink2.king.cutoff.in.id \
    --out ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_related

}

prep_hwe(){
    ###############
    ##  Purpose: Produce list of samples by super pop (probability > 0.8)
    ##  Input: Ancestry assignments
    ##  Output: Sample lists per pop with probability score
    ###############
    module load $RLoad
    tee prep_hwe.R <<EOF  
    library(data.table)
    library(dplyr)
         dat <- fread("/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_actual/Ancestries/aggV2_ancestry_assignment_probs_1KGP3_200K.tsv") %>% as_tibble()
        for(col in c("AFR","EUR","SAS","EAS")){dat[dat[col]>0.8,c("Sample",col)] %>%
        write.table(paste0('/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_actual/Ancestries/',col,"pop.txt"), quote = F, row.names=F, sep = "\t")};
EOF
    chmod +x prep_hwe.R
    Rscript prep_hwe.R
    #Produce an unrelated version of each of these lists too
    for pop in AFR EUR SAS EAS; do
        echo -e "Running ${pop}..."
        awk 'NR==FNR{a[$2]; next} ($1) in a' ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_triangle_HWE10_6.king.cutoff.in.id ${out}/Ancestries/${pop}pop.txt > \
        ${out}/Ancestries/${pop}_unrelated.txt
        awk '{print $1"\t"$1}' ${out}/Ancestries/${pop}_unrelated.txt > ${out}/Ancestries/${pop}_unrelated.keep
    done
    wc -l ${out}/Ancestries/*_unrelated.txt
}

p_hwe(){
    module load $bcftoolsLoad
    module load $plinkLoad
    mkdir -p ${out}HWE
    #Taking the files produced in prep_hwe, run HWE on each subset of samples for each file.
    #It would be faste rto parallelise across files too, but the super-pops apart from EUR shouldn't
    #take too long
        for pop in AFR EUR SAS EAS; do
            echo -e "Calculating HWE on ${pop} samples..."
            plink --bcf ${infile} \
            --hardy midp \
            --keep ${out}/Ancestries/${pop}_unrelated.keep \
            --double-id \
            --allow-extra-chr \
            --out ${out}/HWE/${i}_$pop
        done
}

endAggregateAnnotation(){
    ###############
    ##  Purpose: Annotate and make pass/fail. print subset of cols
    ##  Input: All the outputs of siteQC metrics
    ##  Output: A single text file containing annotated variants
    ###############
    module load $RLoad
    export final='TRUE'
    #Some of the below will need to be changed when we finalise filepaths etc for data
    mkdir -p ${out}Annotation_final
    Rscript annotatePerChunk.R \
    ${i} \
    ${out}startfile/start_file_${i} \
    ${out}missing2/missing1_${i} \
    ${out}missing2/missing2_${i} \
    ${out}medianCoverageAll/medianCoverageAll${i} \
    ${out}medianCoverageNonMiss/medianNonMiss_depth_${i} \
    ${out}medianGQ/medianGQ_${i} \
    ${out}AB_hetAll/hetAll_${i} \
    ${out}AB_hetPass/hetPass_${i} \
    ${out}MendelErrSites/MendErr_${i}.lmendel \
    ${out}Annotation_final \
    ${resources}/N_samples \
    ${out}AC_counts/${i}_AC \
    ignore #unused arg in this argument
}

#Most of the following functions down to make_header() are checks run on the data, not directly needed for the running of the pipeline

linecheck(){
    #Check that the end aggregation files have the same number of lines as the input bcfs
    module load $bcftoolsLoad
    mkdir -p ${out}lineCounts
    l1=`bcftools query -f '%POS\n' ${infile} | wc -l`
    l2=`grep -v \# ${out}Annotation_final/BCFtools_site_metrics_${i}.txt | wc -l `
    printf "$l1\t$l2\t${i}\n" >> ${out}lineCounts/summary.txt
    #For each line now check that the number of output FILTER flags is the correct number
    #Do this by binding onto the summary file for the chunks and N flags
}

sublinecheck(){
    export -f linecheck
    for i in $regions; do
    export i
    #Parse infile paths
    export infile=${input}/${mainname}_${i}.bcf
    bsub -J "linecheck${i}" -q long -P bio -U as1 -cwd ${wd} -e logs/linecheck.err${i}_%J -o logs/linecheck.out${i}_%J  linecheck
    done
}

gnomad_checks(){
#gnomad checks

# cd /public_data_resources/gnomad/v3
##Identify all frequencies for failures
# bcftools view gnomad.genomes.r3.0.sites.vcf.bgz | grep -v \\#\\# | awk '{print $7}' | sort | uniq -c | sort -nr
gnomad_chrom(){
    bcftools view /public_data_resources/gnomad/v3/gnomad.genomes.r3.0.sites.vcf.bgz -t chr${i} | grep -v \\#\\# | awk '{print $7}' |\
    sort | uniq -c | sort -nr > $HOME/summary_chr${i}.txt
}
export -f gnomad_chrom
gnomad_biallelic_snps(){ 
    module load bio/BCFtools/1.10.2-foss-2018b
    bcftools view /public_data_resources/gnomad/v3/gnomad.genomes.r3.0.sites.vcf.bgz -m2 -M2 -v snps | grep -v \\#\\# > $HOME/ALL_snps.txt
}
export -f gnomad_biallelic_snps
#Cut this file down to what we need and then split it into individual files
#cut -f 1,2,3,4,5,6,7 ALL_snps.txt > all_snps_cpraF.txt | awk '{print>$1}' 


bsub -J "gnomad_biallelic" -q long -P bio -U as1 -e $HOME/gnomad_biallelic.err${i}_%J -o $HOME/gnomad_bilallelic.out${i}_%J gnomad_biallelic_snps 
gnomad_all_snps(){ 
    module load bio/BCFtools/1.10.2-foss-2018b
    bcftools view /public_data_resources/gnomad/v3/gnomad.genomes.r3.0.sites.vcf.bgz -v snps | grep -v \\#\\# | awk '{print $7}' |\
    sort | uniq -c | sort -nr > $HOME/summary_all_snps.txt
}
export -f gnomad_all_snps
bsub -J "gnomad_all_snps" -q long -P bio -U as1 -e $HOME/gnomad_allSNPs.err_%J -o $HOME/gnomad_allSNPS.out_%J gnomad_all_snps 
#Per chrom pass rate
for i in {1..22}; do
    export i
    bsub -J "gnomad_chrom${i}" -q long -P bio -U as1 -e $HOME/gnomad_chrom${i}.err_%J -o $HOME/gnomad_chrom${i}.out_%J gnomad_chrom 
    done


#Intersect for 100 random chunks what passes in ours and gnomad
passfail_gnomad(){
    module load $bcftoolsLoad
    #shuf -n 100 ${bedfile} > gnomad_intersect_check_chunks.txt
    mkdir -p gnomad_intersect
    #Which chrom are we on
    infile=`sed -n "${LSB_JOBINDEX}p" gnomad_intersect_check_chunks.txt`
    chr=`echo ${infile} | awk -F '[/_]' '{print $5}'`
    bcftools query /public_data_resources/gnomad/v3/gnomad.genomes.r3.0.sites.vcf.bgz -r ${chr} -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\n' > gnomad_intersect/gnomad_${infile}.txt
}
export -f passfail_gnomad
 #bsub  -J "gnomad_passfail[1-100]" -q long -P bio -cwd ${wd} -U as1 -e $HOME/gnomad_intersect.err${f}_%J -o $HOME/gnomad_intersect.out${f}_%J  passfail_gnomad
#WHEN ANNOTATING - USE THREAD ARGUMENTS AND SUBMIT WITH MULTIPLE CORES


make_header(){
    #Need to correct paths for this
    [[ ! -f all_seen_flags.txt ]] && awk '{print $1}' all_summary_filter.txt | sort | uniq  > all_seen_flags.txt
    #Start with header from input file
    #Remember that we need to use the header with removed samples
    bcftools view -h all_chunks_merged_norm_chr11_52570777_53572719.bcf -S \
    ${resources}78389_final_platekeys_agg_v9.txt --force-samples > orig.hdr
    #Clean out INFO lines, these wil be repopulated
    grep -v \#\#INFO orig.hdr > orig.hdr2 && mv orig.hdr2 orig.hdr
    
    #Now we need to insert our header lines
    tee headerbuild.R << EOF
    #This script runs interactively, doesn't work otherwise, not sure why
        library(data.table);library(dplyr);library(magrittr);library(stringr);
        #Read in the header options
        toAdd <- fread("all_seen_flags.txt", header = F) %>% as_tibble() %>% 
        filter(V1 != 'PASS');
        starter <- '##FILTER=<ID=';
        mid <- ',Description="';
        end <- '.">'
        #Construct descriptions. We won't be all that descriptive for the filter field
        #more so for the info field. Just state failures for filter
        #Probably should add what the failing values are. Maybe only include these for the single instance onces
        #to save on space
        singlecases <- setNames(
            as.list(
                c('Fails depth, median depth < 10',
                'Fails GQ, median GQ <15',
                'Fails missingness, missingness > 0.05',
                'Fails completeGTRatio, compelete sites <0.5',
                'Fails AB ratio, AB ratio of (het sites passing binomial distribution with p-value < 0.01 / all het sites) < 0.25',
                'Fails phwe_eur, site is out of HWE (<10e-6) for inferred unrelated inferred eur superpopulation'
                )),
                c('depth',
                'GQ',
                'missingness',
                'completeGTRatio',
                'ABratio',
                'phwe_eur')
            )

        singlecases <- paste0(starter, names(singlecases), mid, singlecases, end)


    #Now do all the other filter combinations
        toAdd %<>% mutate(V2 = case_when(str_count(V1, ':') == 1 ~ str_replace(V1, ':',' and '),
                               str_count(V1,':') > 1 ~ str_replace_all(V1,':',', '),
                                TRUE ~ V1))
        toAdd %<>% mutate(V2 = ifelse(str_count(V2,',') > 1, sub(".([^,]*)$", " and\\1", V2), V2))


        toAdd %<>% mutate(V3 = paste0(starter,V1, mid, 'Fails ', V2, end))

        #Now sort out the info
        infostart <- '##INFO=<ID='
        infomid <- ',Number=.,Type=Float,Description="'
        infoend <- '">'
    
        infos <- setNames(
            as.list(
                    c('Median depth (taken from the DP FORMAT field) of all samples. Used for filter flag depth.',
                      'Median depth (taken from the DP FORMAT field) from samples with non-missing genotypes.',
                      'Median genotype quality(taken from the GQ FORMAT field) from samples with non-missing genotypes. Used for filter flag GQ.',
                      "Percent of fully missing genotypes (GT = './.' and FORMAT/DP = 0)", #v important to have inverted commas this way round here
                      'The ratio of complete sites/total number of samples',
                      'For each het call, a binomial test is conducted for reads supporting the ref and alt alleles. AB ratio is the hets showing imbalance (p<0.01) divided by the total number of hets.',
                      'The number of Mendel Errors at this site, calculated from confirmed trios',
                      'HWE mid p-value in inferred unrelated inferred afr superpop.',
                      'HWE mid p-value in inferred unrelated inferred amr superpop.',
                      'HWE mid p-value in inferred unrelated inferred eas superpop.',
                      'HWE mid p-value in inferred unrelated inferred eur superpop.',
                      'HWE mid p-value in inferred unrelated inferred sas superpop.'
                      )),
                    c("medianDepthAll",
                    "medianDepthNonMiss",
                    "medianGQ",
                    "missingness",
                    "completeGTRatio",
                    "ABratio",
                    "MendelSite",
                    "phwe_afr",
                    "phwe_amr",
                    "phwe_eas",
                    "phwe_eur",
                    "phwe_sas"
                    )
                )
            #Now build this info full string
            infosout <- paste0(infostart, names(infos), infomid, infos, infoend )

    #Now print this out to file, and add it to the rest of the header
    names(singlecases <- NULL)
    d <- c(singlecases, toAdd$V3, infosout) %>% unlist() %>% as.data.frame() 
    fwrite(d, 'additional_header.txt', quote = F, col.names = F)

EOF
Rscript headerbuild.R
}


annotate_bcfs(){
#Version to generalise
#First prepare annotation file
    module load $BCFtoolsLoad

    autosomes=`grep -v chrX ${bedfile} | grep -v chrY | grep -v chrX | grep -v chrM`
    export infile=`echo $autosomes | cut -d' ' -f ${LSB_JOBINDEX}`

    #Only autosomes
    export i=`echo $infile | awk -F "_" 'sub(/.bcf/,"",$7) {print $5"_"$6"_"$7}'`

    #Change the ABratio name for the files
    sed -i -e '1s/AB_Ratio/ABratio/' ${out}Annotation_final/BCFtools_site_metrics_${i}.txt
    #Compress and tabix files
    bgzip -f ${out}Annotation_final/BCFtools_site_metrics_${i}.txt
    tabix -f -s1 -b2 -e2 ${out}Annotation_final/BCFtools_site_metrics_${i}.txt.gz 

    #Now use final sample list to annotate bcf
    bcftools view -S samplesList.txt ${input}${infile} -Ou |\
    bcftools annotate  \
    -x FILTER,^INFO/OLD_MULTIALLELIC,^INFO/OLD_CLUMPED \
    -a  ${out}Annotation_final/BCFtools_site_metrics_${i}.txt.gz \
    -h additional_header.txt \
    --threads 30 \
    -c CHROM,POS,REF,ALT,missingness,medianDepthAll,medianDepthNonMiss,medianGQ,completeGTRatio,MendelSite,ABratio,phwe_afr,phwe_eur,phwe_eas,phwe_sas,FILTER | \
    bcftools +fill-tags \
    -o /gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data/gel_mainProgramme_aggV2_${i}.vcf.gz \
    -Oz \
    --threads 30 \
    -- -d -t AC,AC_Hom,AC_Het,AC_Hemi,AN 

    bcftools index /gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data/gel_mainProgramme_aggV2_${i}.vcf.gz
}

run_annotate(){
    #ALREADY RUN 1-130
    bsub -J "BCFAnnotate[131-1292]" -q long -P bio -U as1 -cwd ${wd} -o logs/BCF_anotate.out_%J -e logs/BCF_annotate.err_%J annotate_bcfs
}


#### VEP ANNOTATION ####

VEPannotate_missing(){
    ###############
    ##  Purpose: Functional annotation of input bcf
    ##  Input: Input bcf
    ##  Output: Annotated vcf.gz
    ###############
    mkdir -p ${out}/VEP_annotation2
    module load bcftools/1.10.2
    module load vep/98
    export PERL5LIB=$PERL5LIB:${LOFTEE38}
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.bcf/,"",$7) {print $5"_"$6"_"$7}'`
    export outfile=${out}/VEP_annotation2/VEP_annotation_${i}.vcf.gz
    #Not running --everything, as won't give proper gnomAD and TOPMed output
    #Note --nearest transcript requires additional perl module
    bcftools view ${input}${infile} -G |  \
     bcftools annotate -x ^INFO/OLD_MULTIALLELIC,INFO/OLD_CLUMPED -Ov | \
    vep --cache \
    --offline \
    --format vcf \
    --vcf \
    --assembly GRCh38 \
    --dir_cache /tools/apps/vep/98/ensembl-vep/.vep \
    --cache_version 98 \
    --verbose \
    --species homo_sapiens \
    --no_stats \
    --fasta /public_data_resources/reference/GRCh38/GRCh38Decoy_no_alt.fa \
    --sift b \
    --polyphen b \
    --ccds \
    --uniprot \
    --hgvs \
    --symbol \
    --numbers \
    --domains \
    --regulatory \
    --canonical \
    --protein \
    --biotype \
    --uniprot \
    --tsl \
    --appris \
    --gene_phenotype \
    --af \
    --af_1kg \
    --af_esp \
    --max_af \
    --pubmed \
    --variant_class \
    --mane \
    --overlaps \
    --custom /public_data_resources/gnomad/v3/gnomad.genomes.r3.0.sites.vcf.bgz,gnomADg,vcf,exact,0,AF,AF_afr,AF_amr,AF_asj,AF_eas,AF_sas,AF_fin,AF_nfe,AF_oth,AF_ami,AF_male,AF_female \
    --custom /public_data_resources/TOPMed/allele_frequencies/bravo-dbsnp-all.vcf.gz,topmedg,vcf,exact,0,AF,SVM \
    --custom /public_data_resources/phylop100way/hg38.phyloP100way.bw,PhyloP,bigwig \
    --custom /public_data_resources/vep_resources/Build-38/gerp_conservation_scores.homo_sapiens.GRCh38.bw,GERP,bigwig \
    --custom /public_data_resources/clinvar/20190219/clinvar/vcf_GRCh38/clinvar_20190219.vcf.gz,ClinVar,vcf,exact,0,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNVC,CLNVCSO,CLNVI \
    --plugin LoF,loftee_path:${LOFTEE38},human_ancestor_fa:${LOFTEE38HA},gerp_bigwig:${LOFTEE38GERP},conservation_file:${LOFTEE38SQL} \
    --plugin CADD,/public_data_resources/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz \
    --plugin SpliceRegion \
    --plugin SpliceAI,snv=/public_data_resources/SpliceAI/Predicting_splicing_from_primary_sequence-66029966/genome_scores_v1.3/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/public_data_resources/SpliceAI/Predicting_splicing_from_primary_sequence-66029966/genome_scores_v1.3/spliceai_scores.raw.indel.hg38.vcf.gz \
    --compress_output bgzip \
    --force_overwrite \
    --fork 4 \
    --output_file ${outfile}
}

missingVEPs(){
    bsub -J "vepAnnotate[1-1371]%250" -q pipeline -P bio -cwd ${wd} -o logs/vepmiss_annotate.out_%J -e logs/vepmiss_annotate.err_%J VEPannotate_missing
}



# #Workflow for stuff that did't work - namely pc-relate
# #Now use flashpca
#     flashpca --bfile ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_unrelated \
#      _flashpca.txt
# #This outputs eigenvalues + vectors, PCs and PVE
# #Once for unrelated
#  bsub -J 'flashpca' -P bio -q long -U as1 -cwd ${wd} -e logs/flashpca.err_%J -o logs/flaspca.out_%J 'module load bio/FlashPCA2/2.0-foss_2018b; flashpca --bfile ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_unrelated -f _unrelated.txt --outload loadings.txt --outmeansd meansd.txt'
# #Then project unrelated
# bsub -J 'flashpca' -P bio -q long -U as1 -cwd ${wd} -e logs/flashpca.err_%J -o logs/flaspca.out_%J 'module load bio/FlashPCA2/2.0-foss_2018b; flashpca --bfile ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_related --inmeansd meansd.txt --outproj projections.txt --inload loadings.txt -f _related.txt'
# }


# #Now we feed this to pcrelate
# #Going to just write the R script here
# #This currently doesn't work due to significant issues with PCrelate implementation. See pcrelate related R scripts
# pcrelate_run(){
#     module load lang/R/3.6.0-foss-2019a
# R -e 'source("/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/src/project_setup.R");
# library("gdsfmt");
# library("SeqArray");
# library("GENESIS");
# library("GWASTools");
# library("SNPRelate");
# setwd("/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_actual/KING/matrix");
# genodat <- GdsGenotypeReader(filename = "autosomes_merged_genos.gds");
# genoData <- GenotypeData(genodat);
# genoiter <- GenotypeBlockIterator(genoData);
# pcs_related <- fread("pcs_related.txt");
# pcs_unrelated <- fread("pcs_unrelated.txt");
# pcs <- pcs_related %>% bind_rows(pcs_unrelated);
# pc <- data.matrix(pcs[,3:ncol(pcs)]);
# rownames(pc) <- pcs$IID
# unrels <-fread("plink2.king.cutoff.in.id")

# pcRELATE <- pcrelate(genoiter,
#                      pcs = pc, 
#                      training.set = unrels$IID,
#                      ibd.probs = FALSE,
#                      sample.block.size = 10000,
#                      sample.
#                      )

# saveRDS(pcRELATE, "pcrelate1.RDS")
# '
# }
# #bsub -J 'pcrelate2' -P bio -q long -U as1 -cwd ${wd} -e logs/pcrelate.err_%J -o logs/pcrelate.out_%J pcrelate_run

# hwe_prune_preking(){
#     mkdir -p ${out}hwe_prune
#     module load bio/PLINK/1.9b_4.1-x86_64
#     module load lang/R/3.6.0-foss-2019a

#     tee hwe_1000kgp3_unrels.R << EOF
#         library(data.table);
#         library(dplyr);
#         #Read in the 1000KGP3 data and produce output of just unrelated by super pop
#         dat <- fread("/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_actual/hwe_prune/1KGP3.codes") %>% as_tibble();
#         unrels <- fread("/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_actual/hwe_prune/UNRELATED_1KGP3", header = F) %>% as_tibble();
#         d <- unrels %>% left_join(select(dat, Sample, Super_Population), by = c('V1'='Sample'));
#         for(col in c("AFR","EUR","SAS","EAS")){
#             #All we need to output is the sample list for each pop
#             d$V1[d$Super_Population == col] %>% write.table(paste0("/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_actual/hwe_prune/1000KGP3_",col,"_pop.txt"), quote = F, row.names=F, sep = "\t", col.names = F)}
# EOF
# Rscript hwe_1000kgp3_unrels.R #this isn't running via R script, works interactively

# #Now we subset our data to just the SNPs we want, and the unrels we want for our HWE
# for pop in AFR EUR SAS EAS; do
#     echo "Working on ${pop}..."
#     plink \
#     --keep ${out}/hwe_prune/${pop}keep \
#     --extract ${out}BED/autosomes_LD_pruned_1kgp3Intersect.bim
#     --make-bed \
#     --real-ref-alleles \
#     --bfile ${out}/hwe_prune/aggV2_bedmerge_30KSNPs_labkeyV9_unrelatedKING \
#     --out ${out}/hwe_prune/${pop}

#     #Now using these keep files, 
#     plink --bfile ${out}/hwe_prune/aggV2_bedmerge_30KSNPs_labkeyV9_unrelatedKING --hardy midp --out ${out}/hwe_prune/${pop} --nonfounders
# done
### END ###