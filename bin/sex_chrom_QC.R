

# Based on MWs scripts:
# https://cnfl.extge.co.uk/display/REBS/Chromosome+Y+and+X+QC
# Expects input filename, file path, and output path
# This will run over each of the sex chrom chunks

# Setup -------------------------------------------------------------------
source('project_setup.R')


args <- commandArgs(TRUE)

infile <- commandArgs[1]
fileloc <- commandArgs[2]
output <- commandArgs[3]



# Define input filenames --------------------------------------------------
lmiss <- file.path(fileloc, paste0("male_missing_",infile,".lmiss"))
dp <- file.path(fileloc, paste0("male_depth_", infile, ".ldepth.mean"))
gq <- file.path(fileloc, paste0("male_GQ_mean_", infile))

# Start -------------------------------------------------------------------

chr <- strsplit(as.character(infile),"_", fixed = TRUE )[[1]][1]
print(c(chr))
cat("this is", chr)

start <- fread(file.path(fileloc, "male_start_file", infile))

  

