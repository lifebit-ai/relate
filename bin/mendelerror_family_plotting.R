#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(stringr)

args <- commandArgs(trailingOnly = T)

options(bitmaptype="cairo")
indir <- list.files(path = args[1],
                    pattern = '*.fmendel',
                    full.names = T)
outfile <- file.path(args[1], 'SD_mendel_family.png')
fmendel <- lapply(indir, fread) 
fmendel %<>% 
  bind_rows() %>% 
  as_tibble() %>%
  group_by(FID, PAT, MAT, CHLD) %>%
  summarise(N=sum(N))


#Plot distribution
summary(fmendel$N) %>%
  capture.output(file=file.path(args[1], 'FamilyWise.summarystats'))

sdv <- 4*sd(fmendel$N, na.rm = T)
mn <- mean(fmendel$N, na.rm = T)

fmendel %>% 
  filter(N < mn+sdv & N > mn+(-1*sdv)) %>%
  ungroup() %>%
  select(FID) %>%
  distinct() %>%
    write.csv(file.path(args[1],'MendelFamilies_4SD.fam'),
            row.names = F,
            quote = F)
md <- median(fmendel$N, na.rm = T)
sd <- sd(fmendel$N)
png(filename = outfile, 
    width = 6, height = 6,
    units = 'in', res = 300)
fmendel %>% ggplot(aes(N)) +
  geom_density(fill = 'grey35', alpha = 0.7) + 
  geom_vline(xintercept = mn,
             linetype = 'dashed',
             colour = "#29c440") +
  geom_vline(xintercept = md,
             linetype = 'dashed',
             colour = "#29c48b") +
  geom_vline(xintercept = sd*2,
           linetype = 'dotted',
           colour = "#a4dd29") +
  geom_vline(xintercept = sd*-2,
             linetype = 'dotted',
             colour = "#a4dd29") +
  geom_vline(xintercept = sd*3,
             linetype = 'dotted',
             colour = "#dd8329") +
  geom_vline(xintercept = sd*-3,
           linetype = 'dotted',
           colour = "#dd8329") +
  geom_vline(xintercept = sd*4,
             linetype = 'dotted',
             colour = "#dd4a29") +
  geom_vline(xintercept = sd*-4,
             linetype = 'dotted',
             colour = "#dd4a29")
dev.off()


### END ###
