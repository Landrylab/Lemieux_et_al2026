# AUTHOR : Pascale Lemieux
# extract read count from aggregated fasta file

#import packages and functions

library(tidyverse)
library(Biostrings)
library(magrittr)
library(ggpubr)
library(ggsci)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory with the raw data
setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/')


# read reference file 
file_dir <- 
  read_csv('ref_file.csv')[, -c(1:2)]

ref_seq <- 
  read_csv('df_variants.csv')

spiked_seq <- 
  read_csv('spiked_seq.csv')
# select PRM in the screen
PRM <- c('PRM_0152', 'PRM_0246','PRM_0250', 'PRM_0299', 'PRM_0363', 'PRM_0366', 'PRM_0385')

ref_seq%<>%
 filter(family %in% PRM)

colnames(ref_seq)[2] <- 'aa_seq'
ref_seq$dna_seq <- str_c(ref_seq$dna_seq, 'TAA')

spiked_seq$dna_seq <- str_c(spiked_seq$dna_seq, 'TAA')
colnames(spiked_seq)[1] <- 'family'
spiked_seq$codon1 <- spiked_seq$family
spiked_seq$codon2 <- spiked_seq$family

all_ref <- bind_rows(ref_seq[, c('family','aa_seq', 'dna_seq', 'codon1', 'codon2')], 
                     spiked_seq[, c('family','aa_seq', 'dna_seq','codon1', 'codon2')])


ref_condition <- 
expand.grid('timepoint' = c('T0', 'T1', 'T2'), 
            'replicate' = c('R1', 'R2', 'R3'),
            'assay' = c('PPI', 'Avail'))

all_ref <- 
expand_grid(ref_condition, all_ref)



# import of fasta file
# import all reads using custom function

all_sample <- tibble()
for(f in file_dir$Agg_file){
  reads <- import_reads(sample_name = f, suffix = 'aggregated/')
  reads$sample <- f
  all_sample <- bind_rows(all_sample, reads)
}

# dissect sample name to identify sample conditions
conditions <- matrix(unlist(strsplit(all_sample$sample, '_')), byrow = T, ncol = 4)[, 1:3]

all_sample <- bind_cols(all_sample, conditions)
colnames(all_sample)[5:7] <- c('assay', 'timepoint', 'replicate')

all_sample$number <- 
  as.numeric(all_sample$number)

# verify number of reads post aggregation
sum(all_sample$number)
#   21583074

# merge with dataframe with the reference sequences

reads_info <- 
merge(all_sample, 
      all_ref, 
      by.x = c('sequence', 'timepoint', 'replicate', 'assay'), 
      by.y = c('dna_seq','timepoint', 'replicate', 'assay'), 
      all.y = TRUE)

reads_info$number <- as.numeric(reads_info$number)
# verify number of reads post merging, not that much lost
sum(reads_info$number, na.rm  =T)
# 20322511

reads_info <- 
merge(reads_info, 
      ref_seq, 
      by.x = c('sequence', 'aa_seq', 'family','aa_seq', 'codon1', 'codon2'), 
      by.y = c('dna_seq', 'aa_seq', 'family','aa_seq', 'codon1', 'codon2'), 
      all = TRUE)

# save reads info to file
write_csv(reads_info, '~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/reads_info.csv')

reads_info <- read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/reads_info.csv')

# Extract read count per sample and per PRM family
reads_info%<>%
  dplyr::group_by(sample)%>%
  dplyr::mutate(read_sample = sum(number))

reads_info[!(grepl(pattern = 'PRM', reads_info$family)), 'family'] <- 'reference'

reads_info%>%
  dplyr::group_by(sample, family)%>%
  dplyr::mutate(read_family = sum(number)) -> sum_reads

sum_reads <- 
unique(sum_reads[, c(3,11,20:21)])

#sum_reads[!(grepl(pattern = 'PRM', sum_reads$family)), 'family'] <- 'reference'
sum_reads <- sum_reads[!(is.na(sum_reads$read_family)), ]

sum_reads$sample <- 
gsub(pattern = '_agg.fasta', replacement = '', sum_reads$sample, fixed = T)

sum_reads$family <- 
gsub(pattern = 'PRM_0', '', sum_reads$family)

FigS4A <- 
ggplot(sum_reads)+
  geom_col(aes(x = sample, y = read_family, fill = family), color = 'black')+
  #scale_y_continuous(transform = 'log2')+
  #scale_fill_frontiers(alpha = 0.9)+
  scale_fill_viridis_d(option = 'G', direction = -1)+
  ylab('read count')+
  t+
  theme(axis.text.x = element_text(angle =45, vjust = 1, hjust = 1))+
  guides(fill = guide_legend(title = 'peptide family', ))

saveRDS(FigS4A, '~/PL_projects/PL_papers/PPI_optimization_paper/figures/FigS4A.rds')

sum_reads%>%
  dplyr::group_by(sample)%>%
  dplyr::mutate(sum_read = sum(read_family))


sum_all <- unique(sum_reads[, c(2,3)])

sum(sum_all$read_sample)
# 20M passed all the steps

# from fastq 47 259 242 raw reads R1 and R2
# from fastq trimmed 47259242 reads R1 and R2  - > 23.5M
# from fasta merged 21681666 reads (lost only 2M)

# 1 450 553 from NA values, not significant

# Not that much lost, 23M in fastq ended up with 20M. 
