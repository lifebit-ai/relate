#!/usr/bin/env Rscript

# @Author: daniel 
# @Date: 2019-11-08T10:44:09+00:00
# @Email:daniel.rhodes@genomicsengland.co.uk
# @Project: Aggregate_QC  
# @Last modified by:daniel 

#Really this script mostly worries about fam file generation, specifically for
#trios for calc of mendel errors


#Takes working directory, input file and output filename as args

# Setup -------------------------------------------------------------------
args <- commandArgs(trailingOnly = T)
input <- as.character(args[1])
inaggregate <- as.character(args[2])
famout <- as.character(args[3])
keepout <- as.character(args[4])
#Or locally
#setwd('Documents/AggregateVCF/variant_QC/GEL/Data/')
source("../src/project_setup.R")
### 


# Functions ---------------------------------------------------------------
prep_fam <- function(dat, extendedTrio = F){
  #If true, extended trios allows multiple lines for parents in more than one
  #trio (i.e. mother, father proband: mother, father unaffected sibling)
  #This can be a problem for some calcs for PLINK.
  #If F , set to keep proband containing trios
  tmp <- dat %>% 
    distinct() %>%
    group_map( ~ .x %>% 
                 mutate(PAT_p =
                          ifelse(grepl('^Father$',
                                       trio),
                                 sampleId, NA),
                        MAT_p = 
                          ifelse(grepl('^Mother$', 
                                       trio),
                                 sampleId, NA)) %>%
                 #From here pad the PAT_p and MAT_p so proband info is filled,
                 #and then revert the data to 0 in mother and father
                 group_by(family_id) %>% 
                 tidyr::fill(PAT_p, .direction = 'up') %>%
                 tidyr::fill(MAT_p, .direction = 'up') %>%
                 mutate(PAT_p = 
                          ifelse(trio %in% c('Father', 'Mother'),
                                 as.character(0), PAT_p),
                        MAT_p = 
                          ifelse(trio %in% c('Father', 'Mother'),
                                 as.character(0), MAT_p)) %>% 
                 dplyr::select(family_id, 
                               sampleId,
                               PAT_p, MAT_p, everything()) %>%
                 arrange(family_id, factor(trio, levels = c('Offspring',
                                                            'Mother',
                                                            'Father'))) %>%
                 ungroup() %>%
               mutate(sex = ifelse(sex == 'FEMALE', 2, 1)) %>%
                 mutate(case = -99)
    )
  tmp <- tmp[[1]] %>% 
    group_by(family_id) %>% 
    filter(n() > 2) %>%
    ungroup()#we don't want duos, but we want to be more specific in 
  #filtering trios
  if(!extendedTrio){
    #Only keep one of the offspring, preferentially keep proband
    tmp %<>% 
      mutate(offspring = case_when(PAT_p == 0 & sex == 1 ~ 1,
                                   MAT_p == 0 & sex == 2 ~ 2,
                                   TRUE ~ 0)) %>%
      mutate(isProband = ifelse(isProband == 'TRUE','A','B')) %>% #distinct sorts alphabetically   
      arrange(family_id, offspring, isProband) %>% 
      distinct(family_id, offspring, .keep_all = T ) %>%
      group_by(family_id) %>%
      filter(n() > 2) %>% #we have to do this again because we may have trios that 
    #were composed of parent offspring offpspring, so will now be duo
      select(-offspring, -isProband) %>%
      ungroup()
  }
  tmp %<>% select(-trio)
  return(tmp)
}

# Start -------------------------------------------------------------------
#Start from trio list supplied by Chris
dat <- fread(input) %>%
  as_tibble()

#Filter to samples actually in the aggregate
present <- fread(inaggregate, header = F) %>%
  as_tibble()

#Grep out the name etc. then filter and check for complete trios  
present %<>% 
  mutate(V1 = basename(V1)) %>%
  mutate(V1 = gsub('\\.genome.*','', V1))
cat('Assumes stable file naming convention to pull platekey\n')

#Let's also remove any that are UPD cases
cat('Filtering out UPD cases. Check that this is an updated set on running.\n')
cat('Currently hardcoded in trio_define.R\n')
#Info from https://cnfl.extge.co.uk/display/BTS/Uniparental+disomy+cases+detected 
updfile <- file.path(dirname(famout), 'UPD_cases_14012020.csv')
upd <- fread(updfile, header = F) %>% as_tibble()
present %<>% filter(!V1 %in% upd$V1)

cat('Dimensions prior to filtering based on aggregate membership: ', dim(dat), '\n')
dat %<>% filter(sampleId %in% present$V1)  
cat('Dimensions post filtering based on aggregate membership: ', dim(dat), '\n')

#Only keep whole trios
dat %<>% 
  group_by(family_id) %>%
  filter(n() > 2) %>%
  ungroup()
cat('Dimensions post filtering to complete trios: ', dim(dat), '\n')

#Let's produce a FAM file out of this
fam <- dat  %>%
  dplyr::select(family_id,
                sampleId,
                trio,
                sex,
                isProband) %>% 
  as_tibble() %>%
  prep_fam()
cat('Total number of trios: ', nrow(fam)/3, '\n')
#Write to file
fam %>%
  write.table(famout, row.names = F, col.names = F, quote = F)

fam %>% 
  select(family_id, sampleId) %>% 
  mutate(family_id = sampleId) %>%
  write.table(keepout, row.names = F, col.names = F, quote = F)

### END ###