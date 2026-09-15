# AUTHOR : Pascale Lemieux
# extract read count from aggregated fasta file

#import packages and custom functions

library(tidyverse)
library(Biostrings)
library(magrittr)
library(ggpubr)
library(ggsci)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory with the raw data
setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/')


# read reference file 
file_dir <- 
  read_csv('ref_file.csv')[, -c(1:2)]

all_ref<-read_csv('all_reference.csv')

ref_condition <- 
  read_csv('sample_description.csv')
  

# import of fasta file
# import all reads using custom function
all_sample <- tibble()
for(f in file_dir$Agg_file){
  reads <- import_reads(sample_name = f, suffix = 'aggregated_July2025/')
  reads$file <- f
  all_sample <- bind_rows(all_sample, reads)
}

# add reads from May sequencing
file_dir_may <- 
  read_csv('ref_file_May2025.csv')[, -c(1:2)]

for(f in file_dir_may$Agg_file){
  reads <- import_reads(sample_name = f, suffix = 'aggregated_May2025/')
  reads$file <- f
  all_sample <- bind_rows(all_sample, reads)
}

# dissect sample name to identify sample conditions
conditions <- matrix(unlist(strsplit(all_sample$file, '_')), byrow = T, ncol = 4)[, 1:3]

all_sample <- bind_cols(all_sample, conditions)
colnames(all_sample)[5:7] <- c('sample', 'timepoint', 'condition')

all_sample$number <- 
  as.numeric(all_sample$number)

# vizualize read count distribution
ggplot(all_sample)+
  geom_histogram(aes(x = number))+
  scale_x_log10(breaks = c(1, 10, 50, 100, 150, 200, 1000, 1500, 2000))+
  ylab('number of sequences')+
  xlab('read count')+
  #xlim(c(1, 1e7))+
  t+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# merge read counts with sample description
all_sample <- 
merge(all_sample, 
      ref_condition, 
      by.x = c('sample', 'condition', 'timepoint'), 
      by.y = c('Sample', 'Condition', 'Timepoint'))


# verify number of reads post aggregation
sum(all_sample$number)
#   174'311'643 July , 214'593'201 total reads

# sum the read counts per sequence
# sum the reads from both sequencing runs
all_sample%<>%
  dplyr::group_by(sequence, sample, timepoint, condition, Assay, Replicate, Pool)%>%
  dplyr::summarise(number_both = sum(number))%>%
  ungroup()%>%
  unique()

sum(all_sample$number_both)

# merge with dataframe with the reference sequences

all_ref%>%
  select(family, aa_seq, dna_seq)%>%
  unique()->x

reads_info <- 
merge(all_sample, 
      x, 
      by.x = c('sequence'), 
      by.y = c('dna_seq'), 
      all =F)


# verify number of reads post merging, not that much lost
sum(reads_info$number, na.rm  =TRUE)
# 203'438'538

# save reads info to file
write_csv(reads_info, '~/PL_projects/PL_papers/PPI_optimization_paper/data/reads_info_both_PCA2.csv')

# retrieve reads info
reads_info <- read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/reads_info_both_PCA2.csv')

# vizualize read count distribution associated with a sample
ggplot(reads_info)+
  geom_histogram(aes(x = number_both))+
  scale_x_log10(breaks = c(1, 10, 50, 100, 150, 200, 1000, 1500, 2000))+
  ylab('number of sequences')+
  xlab('read count')+
  #xlim(c(1, 1e7))+
  t+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

reads_info%>%
  dplyr::group_by(sample, timepoint, condition, family)%>%
  dplyr::mutate(read_family = sum(number_both, na.rm = TRUE))%>%
  select(timepoint, condition, read_family, family, Pool, Assay, Replicate, sample)%>%
  unique()-> family_reads

sum(family_reads$read_family)

family_reads[!grepl(pattern = 'PRM_0', family_reads$family), 'family'] <- 'reference'

family_reads$family <- 
  gsub(pattern = 'PRM_0', '', family_reads$family)

family_reads$sample_name <- str_c(family_reads$Assay, family_reads$Pool,
 family_reads$condition, family_reads$timepoint, family_reads$Replicate,sep= '_')

# vizualise read count per sample
family_reads%>%
dplyr::filter(family != 'reference')%>%
ggplot()+
  geom_col(aes(x = sample_name, y = read_family, fill = family), color = 'black')+
  scale_fill_manual(values = c("#0B0405FF", "#2E1E3CFF", "#413D7BFF", "#37659EFF", "#348FA7FF", "#40B7ADFF", "#8AD9B1FF"))+
  ylab('read count')+
  t+
  theme(axis.text.x = element_text(angle =45, vjust = 1, hjust = 1))+
  guides(fill = guide_legend(title = 'peptide family', position = 'bottom', nrow = 1))->read_count_per_family

saveRDS(read_count_per_family, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/read_count_per_sample_family.rds',)



# check coverage (read count relative to the expected number of difference sequence)
# per sample
lib_description <- 
  read_csv('library_description.csv')

merge(family_reads, 
      lib_description[, c(1,2,4)], 
      by = c('Pool', 'family'))->coverage_read


coverage_read%<>%
  dplyr::group_by(sample_name)%>%
  dplyr::mutate(diversity_sample = sum(diversity))%>%
  dplyr::ungroup()%>%
  dplyr::group_by(sample_name, family)%>%
  dplyr::mutate(diversity_family = sum(diversity))%>%
  dplyr::mutate(coverage_family = read_family/diversity_family)%>%
  dplyr::ungroup()



ggplot(coverage_read)+
    geom_tile(aes(x = family, y = sample_name, fill = coverage_family))+
    geom_point(data =coverage_read[coverage_read$timepoint == 'T0' & coverage_read$coverage_family < 100, ], 
               aes(x = family, y = sample_name), color ='red', shape =4)+
  scale_fill_viridis_c(option = 'G', direction = 1)+
  ylab('read count')+
  t+
  theme(axis.text.x = element_text(angle =45, vjust = 1, hjust = 1), 
        legend.position = 'bottom')

  coverage_read%<>%
  dplyr::group_by(sample, condition, timepoint)%>%
  dplyr::mutate(coverage_sample = sum(read_family)/diversity_sample)

#ggsave('./figures/coverage_per_sample_both.png',width = 7, height = 13)

# vizualise coverage for T0 samples
coverage_read%>%
  dplyr::group_by(sample_name, timepoint)%>%
  dplyr::summarise(coverage_sample = sum(coverage_family, na.rm = TRUE))%>%
  ungroup()->coverage_sample

  ggplot(coverage_sample[coverage_sample$timepoint == 'T0', ])+
    geom_col(aes(y = sample_name, 
                 x = coverage_sample))+
    geom_vline(xintercept = 100, color = 'red', linetype = 'dashed')+
    ylab('samples')+
    xlab('coverage per sample (T0)')+
    scale_x_log10(breaks = c(1, 10, 50, 100, 150, 200, 1000, 1500, 2000))+
    t+
    theme(axis.text.x = element_text(angle =45, vjust = 1, hjust = 1))->coverage_t0

saveRDS(coverage_t0, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/coverage_per_sample_T0_both.rds')  
  
ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/PCA2/coverage_per_sample_T0_both.png', width = 7, height = 5)  
         