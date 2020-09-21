
#module load lang/R/3.6.0-foss-2019a

lapply(c("data.table", "tidyverse", "magrittr"), library, character.only = T)

#Produce RF predicting the pops based on projected SNPs on 1000kgp3

#Read population labels
indiv_pop<-fread("/re_gecip/BRS/thanos/ethnicities_aggV2/1KGP3.sample_table")
super_pop<-fread("/re_gecip/BRS/thanos/ethnicities_aggV2/super_pop_codes.tsv")
phen<-merge(indiv_pop,super_pop,by="Population")

#Read plate_keys of 1KGP3 unrelated individuals
k1g<-fread("/re_gecip/BRS/thanos/ethnicities_aggV2/1KGP3_30K_unrel_autosomes.fam")[,1]

#Read 1KGP3 PCs and their projection to aggV2 with 200k SNPs
pcsk1g<-fread("/re_gecip/BRS/thanos/ethnicities_aggV2/aggV2_ethn_200Ksites/hwe10e-6_superpops_195ksnps_1KGP3samples_unrel.eigenvec")[,-2]
pcsgel<-fread("/re_gecip/BRS/thanos/ethnicities_aggV2/aggV2_ethn_200Ksites/autosomes_LD_pruned_1kgp3Intersect_ALL_hwe10e-6.proj.eigenvec")[,-c(2,23)]

#Set up sample names
colnames(k1g)[1]="Sample"

colnames(pcsk1g)[1]="Sample"
colnames(pcsk1g)[2:21]=paste0("PC",1:20)

colnames(pcsgel)[1]="Sample"
colnames(pcsgel)[2:21]=paste0("PC",1:20)

#Create training data
train_data<-merge(phen,k1g,by="Sample")
train_data<-merge(train_data,pcsk1g,by="Sample")

## RANDOM FOREST
#Train on subpops
library(randomForest)
set.seed(123)



#Lets see how many PCs and how many trees we should be using.

grid_search_parameters <- function(pop){
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
  return(res)
}

#Run search
pop_super <- grid_search_parameters(pop = "Super_Population")
pop_sub <- grid_search_parameters(pop = 'Population')

#The above searches suggest that the best performing models are 
pop_super %>% arrange(err.rate)
pop_sub %>% arrange(err.rate)

###Based on the above, going to use the following
# For the sub-pop prediction, 20PCs, and 400 trees
# For the super-pop prediction, 8PCs, and 400

rfdat=train_data[,1:(5+20)] #the plus indicates the number of PCs included
#Train random forest algorithm on individual population labels
fit_pop=randomForest(rfdat[,6:ncol(rfdat)],
                     y=as.factor(rfdat$Population),
                     ntree=400,
                     keep.inbag = T)


rfdat=train_data[,1:(5+8)]
#Train random forest algorithm on individual population labels
fit_super=randomForest(rfdat[,6:ncol(rfdat)],
                       y=as.factor(rfdat$Super_Population),
                       ntree=500,
                       keep.inbag = T)



#Lets see how we do in our sub pop when we compress to super pops
#Take the confusion matrix, rename the variables and sum the matrix by groups
confmat <- fit_pop$confusion
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
adjusted_confmat <- t(tmp)


#Use fitted pop model to predict on projected GEL PCs
pred_GEL_pop = predict(fit_pop,pcsgel,type="prob")
pred_GEL_pop=data.frame(pcsgel$Sample,pred_GEL_pop)


#Use fitted model to predict on projected GEL PCs
pred_GEL_super = predict(fit_super,pcsgel,type="prob")
pred_GEL_super=data.frame(pcsgel$Sample,pred_GEL_super)

pop_prob<-list()
spops=unique(phen$Super_Population)

#pdf("singlepop_ancestry_vs_superanc.pdf",height=6,width=9)
#par(mfrow=c(2,3))

assign_pops <- function(x){
  #For each super-population sum the probabilities of assignement to each subpopulation for pop_prob and get the super-pop predicted values for super_prob.
  for (lpop in 1:length(spops)){
    focal_pop=sort(match(unlist(super_pop[which(super_pop$Super_Population==spops[lpop]),1]),colnames(pred_GEL_pop)))
    pop_prob[[lpop]]=rowSums(pred_GEL_pop[,focal_pop])
    super_prob=pred_GEL_super[,which(colnames(pred_GEL_super)==spops[lpop])]
    
    #plot(super_prob,pop_prob[[lpop]],xlab="Super pop trained probability",ylab="Sum prob of single pops",main=paste0(spops[lpop],",cor=",round(cor(super_prob,pop_prob[[lpop]]),2)),pch=20)
    #abline(0,1)
  }
  #dev.off()
  
  #Bind and format final tables
  superpops_probs=data.frame(pcsgel$Sample,do.call("cbind",pop_prob))
  colnames(superpops_probs)=c("Sample",spops)
  colnames(pred_GEL_pop)[1]="Sample"
  
  #merge individual pop probabilites and super-pop probabilities
  ancestries=merge(superpops_probs,pred_GEL_pop,by="Sample")
}
#write-out inferred ancestries
write.table(ancestries,"aggV2_ancestry_assignment_probs_1KGP3_200K.tsv",row.names=F,quote=F)


#################
###Further investigation

lapply(c('ggplot2','dplyr','magrittr','randomForest','data.table'), 
       library, as.character = T)
##### Checking to see what the ancestry assignment is for the HWE-2 + HWE-10 data based on varying thresholds of probability
fit_model <- function(pcsgel, pcsk1g){
  #Read population labels
  indiv_pop<-fread("/re_gecip/BRS/thanos/ethnicities_aggV2/1KGP3.sample_table")
  super_pop<-fread("/re_gecip/BRS/thanos/ethnicities_aggV2/super_pop_codes.tsv")
  phen<-merge(indiv_pop,super_pop,by="Population")
  
  #Read plate_keys of 1KGP3 unrelated individuals
  k1g<-fread("/re_gecip/BRS/thanos/ethnicities_aggV2/1KGP3_30K_unrel_autosomes.fam")[,1]
  
  #Set up sample names
  colnames(k1g)[1]="Sample"
  
  colnames(pcsk1g)[1]="Sample"
  colnames(pcsk1g)[2:21]=paste0("PC",1:20)
  
  colnames(pcsgel)[1]="Sample"
  colnames(pcsgel)[2:21]=paste0("PC",1:20)
  
  #Create training data
  train_data<-merge(phen,k1g,by="Sample")
  train_data<-merge(train_data,pcsk1g,by="Sample")
  #First fit the model
  #Checking on super populations with 6 pcs and 400 trees
  rfdat=train_data[,1:(5+6)] #the plus indicates the number of PCs included
  #Train random forest algorithm on individual population labels
  fit_super=randomForest(rfdat[,6:ncol(rfdat)],
                         y=as.factor(rfdat$Super_Population),
                         ntree=400,
                         keep.inbag = T)
  l <- list()
  l$model <- fit_super
  l$rfdat <- train_data[,1:(5+6)]
  l$pcsgel <- pcsgel
  l$pcsk1g <- pcsk1g
  return(l)
}

pcsk1ghwe2<-fread("/re_gecip/BRS/thanos/ethnicities_aggV2/aggV2_ethn_200Ksites/hwe10e-2_superpops_195ksnps_1KGP3samples_unrel.eigenvec")[,-2]
pcsgelhwe2 <- fread("/re_gecip/BRS/thanos/ethnicities_aggV2/aggV2_ethn_200Ksites/autosomes_LD_pruned_1kgp3Intersect_ALL_hwe10e-2.proj.eigenvec")[,-c(2,23)]
pcsk1ghwe6<-fread("/re_gecip/BRS/thanos/ethnicities_aggV2/aggV2_ethn_200Ksites/hwe10e-6_superpops_195ksnps_1KGP3samples_unrel.eigenvec")[,-2]
pcsgelhwe6 <- fread("/re_gecip/BRS/thanos/ethnicities_aggV2/aggV2_ethn_200Ksites/autosomes_LD_pruned_1kgp3Intersect_ALL_hwe10e-6.proj.eigenvec")[,-c(2,23)]

modhwe2 <- fit_model(pcsgelhwe2, pcsk1ghwe2)
modhwe6 <- fit_model(pcsgelhwe6, pcsk1ghwe6)
#Save the models and associated datat
modhwe2 %>% saveRDS('/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_actual/Ancestries/InvestigateHWE/hwe2_pcs6_trees400.RDS')
modhwe6 %>% saveRDS('/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_actual/Ancestries//InvestigateHWE/hwe6_pcs6_trees400.RDS')


#Now check the labels
pred_labels <- function(model){
  pred_GEL_super = predict(model$model,model$pcsgel,type="prob")
  pred_GEL_super=data.frame(model$pcsgel$Sample,pred_GEL_super) %>% as_tibble()
}

#How are the ancestry assignments across them
maxpop <- function(labs){
  preds <- c('AFR','AMR','EAS','EUR','SAS')
  labs %<>%
    mutate(ethnmax = apply(.[,preds], 1,
                           function(x) names(x)[which.max(x)]))
  
  labs %<>% 
    mutate(ethn0.8 =  
             ifelse(apply(.[,preds], 1, function(x) max(x)) > 0.8, 1, 0),
           ethn0.5 = 
                    ifelse(apply(.[,preds], 1, function(x) max(x)) > 0.5, 1, 0))
  return(labs)
}
labhwe2 <- pred_labels(modhwe2) %>% maxpop()
labhwe6 <- pred_labels(modhwe6) %>% maxpop()


out <- left_join(labhwe2, labhwe6, by = 'model.pcsgel.Sample', suffix = c(".hwe2",".hwe6"))
out %>% fwrite('/re_gecip/shared_allGeCIPs/drhodes/Aggregation_79k/out_actual/Ancestries/InvestigateHWE/anc_probs_hwe2_hwe6_multiProbThreshold.tsv',
               sep='\t')

count(out, ethn0.8, ethnmax) %>% filter(ethn0.8 ==1)

#Lets also just look at the concordance with the 30k data
#To do this I have to repeat the above but for the 30k
pcsk1g30<-fread("/re_gecip/BRS/thanos/ethnicities_aggV2/1KGP3_30K_unrel_autosomes.eigenvec")[,-2]
pcsgel30<-fread("/re_gecip/BRS/thanos/ethnicities_aggV2/aggV2_1KGP3_30K_projection.proj.eigenvec")[,-c(2,23)]
k30model <- fit_model(pcsgel30, pcsk1g30)
k30labs <- pred_labels(k30model) %>% maxpop()

tmp <- labhwe6 %>% left_join(k30labs, by='model.pcsgel.Sample', suffix = c('hwe6','k30'))

tmp %>% filter(ethnmaxhwe6 == ethnmaxk30 & ethn0.8hwe6 == ethn0.8k30)






#Compare the probabilities for the train/pred super vs train/pred sub
subsuper <- ancestries %>% select(Sample, AFR, SAS, EAS, EUR, AMR) %>% 
  left_join(pred_GEL_super, by=c('Sample'='pcsgel.Sample'),
            suffix = c('.sub','.super'))
library(ggplot)
library(patchwork)
#Plot some things
df <- subsuper
plots <- lapply(c('AFR','AMR','EUR','SAS','EAS'), function(x) {
  popsub <- paste0(x, '.sub')
  popsuper <- paste0(x,'.super')
  df %>% 
    ggplot(aes_string(popsub, popsuper)) +
    geom_point(alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = 0.8, colour = 'red', linetype = 'dashed') +
    geom_hline(yintercept = 0.8, colour = 'red', linetype = 'dashed') +
    theme_minimal()
})

probs <- c('AFR','SAS','EAS','AMR','EUR')
subprobs <- paste0(probs,'.sub')
superprobs <- paste0(probs, '.super')

df %<>%
  mutate(ethnmax_super = apply(.[,superprobs], 1,
                               function(x) names(x)[which.max(x)]))

df %<>% 
  mutate(ethn0.8_super = 
           ifelse(apply(.[,superprobs], 1, function(x) max(x)) > 0.8, 1, 0))

df %<>%
  mutate(ethnmax_sub = apply(.[,subprobs], 1,
                             function(x) names(x)[which.max(x)])) 
df %<>% 
  mutate(ethn0.8_sub = 
           ifelse(apply(.[,subprobs], 1, function(x) max(x)) > 0.8, 1, 0))

df %<>% mutate(matching = ifelse(
  gsub('\\..*','',ethnmax_sub) == 
    gsub('\\..*','',ethnmax_super), 1, 0 ))
df %>% count(matching) %>% mutate(perc = n/sum(n) *100)
#So we end up with a 95.6% match between our super and sub pop assignments.

#Lets look at those that don't match
df %>% 
  filter(matching == 0) %>% 
  ggplot(aes(EUR.sub, EUR.super)) + 
  geom_point(alpha = 0.3) +
  xlim(c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_vline(xintercept = 0.8, colour = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 0.8, colour = 'red', linetype = 'dashed') +
  theme_minimal() 



#Now we need to bring the 30k SNPs and probabilities into this
k30 <- fread('aggV2_ancestry_assignment_probs_1KGP3_30K.tsv') %>% as_tibble()


#These are sub population based super pop estimates, let's compare the 30k snp sub pops agains 200k subpops
df2 <- df %>% select(Sample, AFR.sub, SAS.sub, EAS.sub, AMR.sub, EUR.sub) %>% 
  left_join(select(k30, AFR, SAS, EAS, AMR, EUR, Sample), by='Sample')

names(df2) <- c('sample', 
                paste0(probs, '.200ksub'), paste0(probs, '.30ksub'))

plots <- lapply(probs,  function(x) {
  popsub <- paste0(x, '.200ksub')
  popsuper <- paste0(x,'.30ksub')
  df2 %>% 
    ggplot(aes_string(popsub, popsuper)) +
    geom_point(alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = 0.8, colour = 'red', linetype = 'dashed') +
    geom_hline(yintercept = 0.8, colour = 'red', linetype = 'dashed') +
    theme_minimal()
})

probs30k <- paste0(probs,'.30ksub')
probs200k <- paste0(probs, '.200ksub')

df2 %<>%
  mutate(ethnmax_30k = apply(.[,probs30k], 1,
                             function(x) names(x)[which.max(x)]))

df2 %<>% 
  mutate(ethn0.8_30k = 
           ifelse(apply(.[,probs30k], 1, function(x) max(x)) > 0.8, 1, 0))

df2 %<>%
  mutate(ethnmax_200k = apply(.[,probs200k], 1,
                              function(x) names(x)[which.max(x)])) 
df2 %<>% 
  mutate(ethn0.8_200k = 
           ifelse(apply(.[,probs200k], 1, function(x) max(x)) > 0.8, 1, 0))

df2 %<>% mutate(matching = ifelse(
  gsub('\\..*','',ethnmax_30k) == 
    gsub('\\..*','',ethnmax_200k), 1, 0 ))
df2 %>% count(matching) %>% mutate(perc = n/sum(n) *100)
