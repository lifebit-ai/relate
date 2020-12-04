#!/usr/bin/env Rscript

# Helper script to combine the HWE files and produce a list of passed
# variants for hwe_pruning_30k_snps process.

############################## ARGUMENTS SECTION #############################
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | "help" %in% args ) {
  cat("
      The helper R Script hwe_produce_pass.R
      Mandatory arguments:
            Script args are hardcoded beacuse it always expects fixed set of input files that
            is always produced by hwe_pops.R script that is run just before current script
            in the hwe_pruning_30k_snps process.

         --help                      - you are reading it


      The typical command for running the script is as follows:

          ./hwe_produce_pass.R

      \n")

  q(save="no")
}


############################## LIBRARIES SECTION #############################

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))


# ############################### SCRIPT SECTION ###############################

#Combine the HWE
dat <- lapply(c("EUR.hwe","AFR.hwe", "SAS.hwe", "EAS.hwe"),fread);
names(dat) <- c("EUR.hwe","AFR.hwe", "SAS.hwe", "EAS.hwe");
dat <- dat %>% bind_rows(.id="id");
write.table(dat, "combinedHWE.txt", row.names = F, quote = F)

#Create set that is just SNPS that are >1e-5 in all pops
dat %>%
    filter(P >1e-5) %>%
    group_by(SNP)   %>%
    count()         %>%
    filter(n==4)    %>%
    select(SNP)     %>%
    distinct()      %>%
    write.table("hwe1e-5_superpops_195ksnps", row.names = F, quote = F)