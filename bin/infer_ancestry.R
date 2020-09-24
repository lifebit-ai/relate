#Produce RF predicting the pops based on projected SNPs on 1000kgp3

# Setup -------------------------------------------------------------------
lapply(c("data.table", "tidyverse", "magrittr"), library, character.only = T)
args <- commandArgs(trailingOnly = T)

argsList <- c('kgp3_sample_table',
              'super_pop_codes',
              'unrel_30k',
              'eigenvec',
              'projections',
              'outdir')
names(args) <= argsList
# Functions ---------------------------------------------------------------
data_setup <- function(mafList, phen, k1g){
  #Set up sample names
  colnames(mafList$pcsk1g)[1]="Sample"
  colnames(k1g)[1]="Sample"
  mafList$k1g <- k1g
  colnames(mafList$pcsk1g)[1]="Sample"
  colnames(mafList$pcsk1g)[2:21]=paste0("PC",1:20)
  
  colnames(mafList$pcsgel)[1]="Sample"
  colnames(mafList$pcsgel)[2:21]=paste0("PC",1:20)
  mafList$train_data<-merge(phen,mafList$k1g,by="Sample")
  mafList$train_data<-merge(mafList$train_data,mafList$pcsk1g,by="Sample")
  mafList <- lapply(mafList, as_tibble)
  return(mafList)
}


grid_search_parameters <- function(train_data, pop){
  #Function runs grid search for PCs and Ntrees as defined in the code below.
  #Supply pop as either Super_population, or Population to train on super or sub pops
  res <- list(ntrees = c(1,10,20,30,40,50,100,200,300,400,500,1000),
              npcs = c(1:20)) %>%
    cross_df() %>%
    mutate(err.rate = NA)
  
  
  for(i in 1:nrow(res)){
    #Define our variables for each run
    ntrees <- res$ntrees[i]
    npcs <- res$npcs[i] 
    rfdat=train_data[,1:(5+npcs)] %>% as_tibble()
    dat <- rfdat[,6:ncol(rfdat)]
    y = select(rfdat, !!pop) %>%
      mutate(inpop = as.factor(!!as.symbol(pop)))
    
    #Train random forest algorithm on individual population labels
    fit_pop=randomForest(dat, y=y$inpop, ntree=ntrees)
    
    #Store the best OOB for each model
    res$err.rate[i] <- min(fit_pop$err.rate[,1])
  }
  res %<>% arrange(err.rate)
  return(res)
}

fit_model <- function(mafList, nPCs, nTrees, pop, ...){
  rfdat= mafList$train_data[,1:(5+nPCs)] #the plus indicates the number of PCs included
  y <- as.factor(pull(rfdat, pop))
  #Train random forest algorithm on individual population labels
  fit_pop <- randomForest(rfdat[,6:ncol(rfdat)],
                          y= y,
                          ntree= 400,
                          keep.inbag = T,
                          ...)
  return(fit_pop)
}

pred_labels <- function(mafList, model){
  out <- paste0('pred_GEL_',gsub('_fit','',model))
  mafList[[out]] <- predict(mafList[[model]],mafList$pcsgel,type="prob")
  mafList[[out]] <- data.frame(mafList$pcsgel$Sample,mafList[[out]]) %>% as_tibble()
  return(mafList)
}


sub_pop_confmat <- function(mafList, super_pop){
  #Only meant to work on the subpops
  #Lets see how we do in our sub pop when we compress to super pops
  #Take the confusion matrix, rename the variables and sum the matrix by groups
  confmat <- mafList$sub_fit$confusion
  
  #Drop the error col
  confmat <- confmat[,-ncol(confmat)]
  superlabs <- super_pop %>% select(Population, Super_Population)
  labs <- colnames(confmat) %>%
    tibble::enframe(name = NULL) %>%
    left_join(superlabs, by=c('value'='Population'))
  colnames(confmat) <- rownames(confmat) <- labs$Super_Population
  tmp <- rowsum(confmat, row.names(confmat))
  tmp <- tmp %>% t()
  tmp <- rowsum(tmp, row.names(tmp))
  mafList[['adjusted_confmat_sub']] <- t(tmp)
  return(mafList)
}

assign_pops <- function(mafList, phen, super_pop){
  spops=unique(phen$Super_Population)
  pop_prob<-list()
  
  #For each super-population sum the probabilities of assignement to each
  #subpopulation for pop_prob and get the super-pop predicted values for super_prob.
  for (lpop in 1:length(spops)){
    focal_pop=sort(
      match(
        unlist(
          super_pop[which(super_pop$Super_Population==spops[lpop]),1]),colnames(mafList$pred_GEL_sub)
      )
    )
    pop_prob[[lpop]]=rowSums(mafList$pred_GEL_sub[,focal_pop])
    super_prob=mafList$pred_GEL_super[,which(colnames(mafList$pred_GEL_super)==spops[lpop])]
  }
  
  #Bind and format final tables
  superpops_probs=data.frame(mafList$pcsgel$Sample,do.call("cbind",pop_prob))
  colnames(superpops_probs)=c("Sample",spops)
  colnames(mafList$pred_GEL_sub)[1]="Sample"
  
  #merge individual pop probabilites and super-pop probabilities
  ancestries=merge(superpops_probs, mafList$pred_GEL_sub,by="Sample") %>%
    as_tibble()
  return(ancestries)
}

# Start -------------------------------------------------------------------

#Read population labels
indiv_pop<-fread(args$kgp3_sample_table)
super_pop<-fread(args$super_pop_codes)
phen<-merge(indiv_pop,super_pop,by="Population")

#Read plate_keys of 1KGP3 unrelated individuals
k1g<-fread(args$unrel_30k)[,1]

#Read 1KGP3 PCs and their projection to aggV2 with 200k SNPs
maf1 <- list()
maf1$pcsk1g<-fread(args$eigenvec)[,-2]
maf1$pcsgel<-fread(args$projections)[,-c(2,23)]

#Prep data
maf1 %<>% data_setup(phen = phen, k1g = k1g)


# Random forest -----------------------------------------------------------
#Train on subpops
library(randomForest)
set.seed(123)

#Lets see how many PCs and how many trees we should be using.

#Run search across both sets
maf1$pop_super <- grid_search_parameters(train_dat = maf1$train_data, pop = "Super_Population")
maf1$pop_sub <- grid_search_parameters(train_dat = maf1$train_data , pop = 'Population')

##
# Manually inspect the error rates to decide the parameters to use
# Remember the method used is pretty short-hand so interpretation may be required
##

#Fit the models for both sets
maf1$super_fit <- fit_model(mafList = maf1, nPCs = 10, nTrees = 1000, pop = "Super_Population")
maf1$sub_fit <-  fit_model(mafList = maf1, nPCs = 20, nTrees = 1000, pop = "Population")


#Predict labels on our data

maf1 %<>% pred_labels(model = 'super_fit')
maf1 %<>% pred_labels(model = 'sub_fit')


#Create confmatrix for the sub-pop to super-population
maf1 %<>% sub_pop_confmat(super_pop)


#Assign those ancestries
maf1$ancestries <- assign_pops(maf1, phen, super_pop)


#Write entire object to file
maf1 %>% saveRDS(file = file.path(args$outdir, 'results.RDS' ))

#Output just the predictions
maf1$ancestries %>% fwrite(file.path(args$outdir,'predicted_ancestries.tsv', sep = '\t'))



#### END ###
