#!/usr/bin/env Rscript

# Helper script to generate population files for hwe_pruning_30k_snps process.

############################# ARGUMENTS SECTION #############################
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) | (length(args) == 1) ) {
  cat("
      The helper R Script hwe_pops.R

      Mandatory arguments:
        --ancestry_assignment_probs='path'    - The path to the ancestry assignment probs file.
        --pc_sancestry_related='path'         - The path to the ancestry Principal Components (PCs) file.

      Optional arguments:
        --help                                - You are reading it.

     Usage:
          The typical command for running the script is as follows:

          ./hwe_pops.R --ancestry_assignment_probs='aggV2_ancestry_assignment_probs.tsv' \n
                --pc_sancestry_related='aggV2_PCsancestryrelated.tsv'

      \n")

  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")

argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)


############################# LIBRARIES SECTION #############################

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))

####################### VARIABLES REASSIGNMENT SECTION ######################

# Facilitates testing and protects from wh-spaces, irregular chars

# required
ancestry_assignment_probs       <- args$ancestry_assignment_probs
pc_sancestry_related            <- args$pc_sancestry_related

cat("\n")
cat("ARGUMENTS SUMMARY")
cat("\n")
cat("ancestry_assignment_probs: ", ancestry_assignment_probs, "\n",sep="")
cat("pc_sancestry_related  : ", pc_sancestry_related,   "\n",sep="")


############################## SCRIPT SECTION ###############################

dat <- fread(ancestry_assignment_probs) %>% as_tibble();
unrels <- fread(pc_sancestry_related) %>% as_tibble() %>% filter(unrelated_set == 1);
dat <- dat %>% filter(plate_key %in% unrels$plate_key);

for (col in c("AFR","EUR","SAS","EAS")) {

    dat[dat[col]>0.8,c("plate_key",col)] %>%
    write.table(paste0(col,"pop.txt"), quote = F, row.names=F)

}