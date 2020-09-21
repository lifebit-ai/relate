#!/usr/bin/env Rscript


# Setup -------------------------------------------------------------------
# all files will be explicilty passed to this script. 
source("../src_actual/project_setup.R")

args <- commandArgs(trailingOnly = T)


argList <- c('pos',
             'start',
             'missing1',
             'missing2',
             'coverageAll',
             'coverageNonMiss',
             'quality',
             'allHets',
             'passHets',
             'mendel',
             'outdir',
             'nsamp',
             'acCount',
             'kgp3txt')


allArgs <- list()
allArgs <- lapply(1:length(argList),
                  function(i){
                    allArgs <- as.character(args[i])
                  })
names(allArgs) <- argList
#Also adding the hardy location here
allArgs$hardy <- '/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_actual/HWE/'
# If king is defined as T in bash env only output pass col
env <- Sys.getenv()
if(any(grepl('^king$', names(env), ignore.case = F))){
  king <- env[grep('king', names(env))][[1]]
} else {
  king <- F
}

#If final is included, then pull hwe data too
if(any(grepl('^final$', names(env), ignore.case = F))){
  final <- env[grep('final', names(env))][[1]]
} else {
  final <- F
}

cat(king)
cat(final)
### 


# Funcs -------------------------------------------------------------------

startProc <- function(x){
  dat <- fread(x$start) %>% 
    as_tibble() %>% 
    mutate(ID = paste0(V1, ":", V2, "-", V3, "/", V4, "-", V5)) %>%
    select(ID, everything()) %>%
    setNames(c("ID","#CHROM","POS","REF","ALT", "OC_OM")) %>%
    select(-OC_OM)
  return(dat)
}
# Start file --------------------------------------------------------------
# dat <- fread(argList$start) %>% 
#   as_tibble() %>% 
#   mutate(ID = paste0(V1, ":", V2, "-", V3, "/", V4, "-", V5)) %>%
#   select(ID, everything()) %>%
#   setNames(c("ID","#CHROM","POS","REF","ALT", "OC_OM")) %>%
#   select(-OC_OM)

#cat('Start file done...\n')

# Missingness -------------------------------------------------------------
missProc <- function(dat, x){
  miss1 <- fread(x$missing1) %>% 
    as_tibble() %>%
    setNames(c("ID", "miss1")) 
  
  miss2 <- fread(x$missing2) %>% as_tibble() %>% setNames(c("ID", "miss2"))
  
  miss <- full_join(miss1, miss2, by="ID") %>%
    mutate(miss1 = replace_na(miss1, 0),
           miss2 = replace_na(miss2, 0)) %>%
    mutate(missingness = miss1/(miss2 + miss1)) %>%
    select(-miss1)
  
  dat %<>% left_join(miss, by="ID") %>%
    mutate(missingness = replace_na(missingness, 0)) #sites as NAs here should
  #be set to fully missing
}

# miss1 <- fread(missing1) %>% 
#   as_tibble() %>%
#   setNames(c("ID", "miss1")) 
# 
# miss2 <- fread(missing2) %>% as_tibble() %>% setNames(c("ID", "miss2"))
# 
# miss <- full_join(miss1, miss2, by="ID") %>%
#   mutate(miss1 = replace_na(miss1, 0),
#          miss2 = replace_na(miss2, 0)) %>%
#   mutate(missingness = miss1/(miss2 + miss1)) %>%
#   select(-miss1)
# 
# dat %<>% left_join(miss, by="ID") %>%
#   mutate(missingness = replace_na(missingness, 0)) #sites as NAs here should
# #be set to fully missing

#cat('Missingness done...\n')
#rm(miss, miss1) #we want to keep miss2 for now

# Coverage ----------------------------------------------------------------
covProc <- function(dat, x){
  covAll <- fread(x$coverageAll) %>% 
    as_tibble() %>%
    setNames(c("ID", "medianCovAll"))
  covNonMiss <- fread(x$coverageNonMiss) %>%
    as_tibble() %>%
    setNames(c("ID", "medianCovNonMiss"))
  
  dat %<>% 
    left_join(covAll, by = "ID") %>%
    left_join(covNonMiss, by = "ID") %>%
    mutate(medianCovNonMiss = replace_na(medianCovNonMiss, 0))#sites as NAs here should
  #be set to 0 coverage
}

# covAll <- fread(coverageAll) %>% 
#   as_tibble() %>%
#   setNames(c("ID", "medianCovAll"))
# covNonMiss <- fread(coverageNonMiss) %>%
#   as_tibble() %>%
#   setNames(c("ID", "medianCovNonMiss"))
# 
# dat %<>% 
#   left_join(covAll, by = "ID") %>%
#   left_join(covNonMiss, by = "ID") %>%
#   mutate(medianCovNonMiss = replace_na(medianCovNonMiss, 0))#sites as NAs here should
# #be set to 0 coverage

#cat('Coverage done...\n')
#rm(covAll, covNonMiss)

# GQ ----------------------------------------------------------------------
gqProc <- function(dat, x){
  GQ <- fread(x$quality) %>%
    as_tibble() %>%
    setNames(c("ID", "medianGQ"))
  
  dat %<>% left_join(GQ, by="ID") %>%
    mutate(medianGQ = replace_na(medianGQ, 0))#sites as NAs here should
  #be set to 0 coverage
}

# GQ <- fread(quality) %>% as_tibble() %>% setNames(c("ID", "medianGQ"))
# 
# dat %<>% left_join(GQ, by="ID") %>%
#   mutate(medianGQ = replace_na(medianGQ, 0))#sites as NAs here should
# #be set to 0 coverage

#cat('GQ done...\n')
#rm(GQ)

# Allelic imbalance -------------------------------------------------------

abProc <- function(dat, x){
  hetAll <- fread(x$allHets) %>%
    as_tibble() %>%
    rename('ID'='V1', 'AB_hetAll'='V2')
  
  hetPass <- fread(x$passHets) %>% 
    as_tibble() %>%
    as_tibble() %>%
    rename('ID'='V1', 'AB_hetPass'='V2')
  
  #We expect fewer passes than all hets, so NAs will be introduced,
  #set these to 0. Essentially this means those without AB scores can still 
  #pass if other flags are fine
  ab <- hetAll %>%
    left_join(hetPass, by='ID') %>%
    mutate(AB_hetPass = replace_na(AB_hetPass, 0)) %>%
    mutate(AB_Ratio = AB_hetPass/AB_hetAll) %>%
    select(ID, AB_Ratio)
  
  dat %<>% 
    left_join(ab, by = "ID") %>%
    mutate(AB_Ratio = replace_na(AB_Ratio, 99)) #scores of 99 will be changed downstream
}


# hetAll <- fread(allHets) %>%
#   as_tibble() %>%
#   rename('ID'='V1', 'AB_hetAll'='V2')
# 
# hetPass <- fread(passHets) %>% 
#   as_tibble() %>%
#   as_tibble() %>%
#   rename('ID'='V1', 'AB_hetPass'='V2')
# 
# #We expect fewer passes than all hets, so NAs will be introduced,
# #set these to 0
# ab <- hetAll %>%
#   left_join(hetPass, by='ID') %>%
#   mutate(AB_hetPass = replace_na(AB_hetPass, 0)) %>%
#   mutate(AB_Ratio = AB_hetPass/AB_hetAll) %>%
#   select(ID, AB_Ratio)
# 
# dat %<>% 
#   left_join(ab, by = "ID") %>%
#   mutate(AB_Ratio = replace_na(AB_Ratio, 99))
#cat('AB done...\n')
#rm(hetAll, hetPass, ab)


# Complete sites ----------------------------------------------------------
completeProc <- function(dat, x){
  nsamples <- fread(x$nsamp) %>%     
    as_tibble()
  dat %<>% mutate( miss2 = replace_na(miss2, 0),
                   completeSites = 1-((nsamples$V1-miss2)/nsamples$V1))
}

# #Here we use the number of samples and the second missingness file to produce
# #a percentage of sites that are complete for a row
# 
# 
# nsamples <- fread(nsamp) %>% as_tibble()
# dat %<>%
#   mutate( miss2 = replace_na(miss2, 0),
#           completeSites = 1-((nsamples$V1-miss2)/nsamples$V1))


# MendelErrors ------------------------------------------------------------
#Mendel errors are being moved from this section, they are not important for 
#pass filters
mendelProc <- function(dat, x){
  me <- fread(x$mendel) %>%
    as_tibble() 
  if(nrow(me) != nrow(dat)){
    stop('Mendel error site file length does not match backbone file')
  } else {
    dat$MendelSite <- me$N
  }
  return(dat)
}


# HWE ---------------------------------------------------------------------
hweProc <- function(dat){
  hardyfiles <- paste0(c('AFR','EUR','EAS','SAS'))
  hwe <- lapply(hardyfiles,
                function(pop) {
                  f <- fread(file.path(allArgs$hardy,
                                       paste0(allArgs$pos,
                                              '_',
                                              pop,
                                              '.hwe'))) %>%
                    as_tibble()
                  names(f) <- c(names(f)[-ncol(f)], paste0('phwe_',tolower(pop)))
                  f <- f[,ncol(f)]
                })
  #Remember format of plink files switches ref and alt, aside from that the
  #vars should be in the same order
  out <- bind_cols(hwe)
  dat %<>% bind_cols(out)
  return(dat)
}

pull_data <- function(){
  #Runs each of the functions above
  dat <- startProc(allArgs)
  cat('Startfile done ...\n')
  dat %<>% missProc(allArgs)
  cat('Missingness done...\n')
  dat %<>% covProc(allArgs)
  cat('Coverage done...\n')
  dat %<>% gqProc(allArgs)
  cat('GQ done...\n')
  dat %<>% completeProc(allArgs)
  cat('Completeness done...\n')
  dat %<>% mendelProc(allArgs)
  cat('Mendel sites done...\n')
  dat %<>% abProc(allArgs)
  cat('Data collated\n')
  return(dat)
}

king_parse <- function(dat){
  #Now let's filter down to biallelic SNPs only
  dat %<>% 
    filter(grepl('chr\\d+:\\d+-[[:alpha:]]/[[:alpha:]]-\\.-\\.',ID))
  
  #For the king run we want AC count and 1kgp3 sites included
  ## AC Count 
  ac <- fread(allArgs$acCount) %>% 
    as_tibble() %>%
    setNames(c("ID", "AC"))
  dat %<>% left_join(ac, by=c('ID')) 
  
  ## 1KGP3 sites
  kgp3 <- fread(allArgs$kgp3txt, header = T) %>% 
    as_tibble() %>%
    mutate(ID=paste0('chr',`#CHROM`,':',POS,'-',REF,'/',ALT))
  
  dat %<>% 
    mutate(tmpID=paste0(`#CHROM`,':',POS,'-',REF,'/',ALT)) %>%
    mutate(kgp3 = ifelse(tmpID %in% kgp3$ID, 1, 0)) %>% 
    filter(kgp3 == 1 & AC >= 20)
  return(dat)
}

king_filter <- function(dat){
  #If we are doing this for KING, we want to filter out all irrelevant SNPs
  #These are sites with AC >=20, 
  dat %<>% 
    mutate(FILTER = 
             ifelse((missingness <= 0.01 & 
                       medianCovAll >= 30 & 
                       medianGQ >= 30 &
                       AB_Ratio >= 0.90 &
                       completeSites >= 0.90), 
                    "PASS", 'NA'),
           FILTER = ifelse(missingness > 0.01, paste0(FILTER, ':missingness'), FILTER),
           FILTER = ifelse(medianCovAll < 30, paste0(FILTER, ':depth'), FILTER),
           FILTER = ifelse(AB_Ratio < 0.90, paste0(FILTER, ':ABratio'), FILTER),
           FILTER = ifelse(completeSites < 0.9, paste0(FILTER, ':completeGTRatio'), FILTER),
           FILTER = ifelse(medianGQ < 30, paste0(FILTER, ':GQ'), FILTER)) %>%
    mutate(FILTER =
             ifelse(grepl('^NA:', FILTER), str_sub(FILTER, 4), FILTER)) 
  return(dat)
}

standard_filter <- function(dat){
  dat %<>%  
    mutate(FILTER = 
             ifelse((missingness <= 0.05 & 
                       medianCovAll >= 10 & 
                       medianGQ >= 15 &
                       completeSites >= 0.5 &
                       AB_Ratio >= 0.25 & 
                       phwe_eur >= 10e-6), 
                    "PASS", 'NA'),
           FILTER = ifelse(missingness > 0.05, paste0(FILTER, ':missingness'), FILTER),
           FILTER = ifelse(medianCovAll < 10, paste0(FILTER, ':depth'), FILTER),
           FILTER = ifelse(AB_Ratio < 0.25, paste0(FILTER, ':ABratio'), FILTER),
           FILTER = ifelse(completeSites < 0.5, paste0(FILTER, ':completeGTRatio'), FILTER),
           FILTER = ifelse(medianGQ < 15, paste0(FILTER, ':GQ'), FILTER),
           FILTER = ifelse(phwe_eur < 10e-6, paste0(FILTER, ':phwe_eur'), FILTER)) %>%
    mutate(FILTER =
             ifelse(grepl('^NA:', FILTER), str_sub(FILTER, 4), FILTER))
  return(dat)
}


data_clean <- function(dat){
  cat('Cleaning data...\n')
  dat %<>% 
    mutate_at(.vars = vars(-ID,
                           -`#CHROM`,
                           -POS,
                           -REF,
                           -ALT,
                           -FILTER),
              .funs = list( ~ gsub("(\\.?<![0-9])0+", "", 
                                   round(x = ., digits = 3),
                                   perl = TRUE)))  %>% 
    mutate_at(.vars = vars(),
              .funs = list(~ replace_na(., '.')))
  
  dat %<>% select(-ID)
  
  #Scrub the 99 AB_ratio scores, revert to '.'
  dat %<>% mutate(AB_Ratio = ifelse(AB_Ratio == 99, '.', AB_Ratio))
  dat %<>% select(-'miss2')
  return(dat)
}

# Start -------------------------------------------------------------------

dat <- pull_data()

#Add HWE data if necessary
if(final){ dat %<>% hweProc }

dat[,5:ncol(dat)] %>% summarise_all(list( ~(sum(is.na(.))))) %>%
  as.data.frame()



# Populate VCF FILTER field -----------------------------------------------
#This operates in two streams, if it is for king, we are using more stringent
#filters.
if(king){
  dat %<>% king_parse %>% 
    king_filter()
} else {
  dat %<>% standard_filter()
}
cat('VCF filter field populated...\n')

if(final){
  #Tidy up the data a bit - will require the right columns for final data output
  dat %<>% data_clean()
  #Redefine some of the names to keep things uniform
  dat %<>% rename('medianDepthAll'='medianCovAll',
                  'medianDepthNonMiss'='medianCovNonMiss',
                  'completeGTRatio'='completeSites')
}



#Write data to file
if(king){
  filename <- paste0(allArgs$outdir, '/BCFtools_site_metrics_SUBCOLS', allArgs$pos,'.txt')
  dat %<>% select('#CHROM',
                  'POS',
                  'REF',
                  'ALT',
                  'FILTER') %>%
    filter(FILTER == 'PASS') %>%
    write.table(filename,
                sep ='\t',
                col.names = T,
                row.names = F, 
                quote = F)
  
} else{
  filename <- paste0(allArgs$outdir, '/BCFtools_site_metrics_', allArgs$pos,'.txt')
  dat %>% 
    write.table(filename,
                sep ='\t',
                col.names = T,
                row.names = F, 
                quote = F)
}
cat('END')
### END ###