# AUTHOR : Pascale Lemieux
# analysis of sequence frequency changes 

#import packages and functions

library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggsci)

source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# set working directory
setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/')

#import dataset
reads_info <- read_csv('reads_info.csv')
od_values <- read_csv('od_values.csv')


# compute number of generation per sample x timepoint
od_values$gen_T1 <- 
unlist(lapply(od_values$OD_end_T1, get_generation))

od_values$gen_T1_T2 <- 
  unlist(lapply(od_values$OD_end_T2, get_generation))

od_values$gen_T2 <- od_values$gen_T1 + od_values$gen_T1_T2

gen_values <- 
pivot_longer(od_values[,c(1,2,9,11)], 
             names_to = 'timepoint', 
             cols = c(gen_T1, gen_T2))

gen_values$timepoint <- 
gsub(pattern = 'gen_', replacement = '', gen_values$timepoint)

colnames(gen_values)[4] <- 'n_gen'

# replace na values for read count by 0
reads_info$number <- 
  replace_na(reads_info$number, 0)

# add 1 to the read count of each sequence

reads_info%<>%
  mutate(cor_count = number+1)

# long to wide format
seq_wide <- 
pivot_wider(reads_info[, -c(9,11)], id_cols = c('sequence', 'aa_seq', 'codon1', 'codon2', 'deg1', 'deg2', 'pos1', 'pos2', 'max1', 'max2', 'family'), 
            names_from = c('assay', 'timepoint', 'replicate'), values_from = 'cor_count')

# select for sequences with > 20 reads at T0 in all replicates
seq_wide%>%
  select(PPI_T0_R1, PPI_T0_R2, PPI_T0_R3, sequence)%>%
  filter(PPI_T0_R1 > 20 & PPI_T0_R2 > 20 & PPI_T0_R3 > 20)-> PPI_sequence

seq_wide%>%
  select(Avail_T0_R1, Avail_T0_R2, Avail_T0_R3, sequence)%>%
  filter(Avail_T0_R1 > 20 & Avail_T0_R2 > 20 & Avail_T0_R3 > 20)-> Avail_sequence

col <- grep(colnames(seq_wide), pattern = 'PPI')

PPI_data <- 
  seq_wide[, c(1:11, col)]
PPI_data <- 
  PPI_data[PPI_data$sequence %in% PPI_sequence$sequence, ]


col <- grep(colnames(seq_wide), pattern = 'Avail')
Avail_data <- 
  seq_wide[, c(1:11, col)]

Avail_data <- 
Avail_data[Avail_data$sequence %in% Avail_sequence$sequence, ]

# compute frequency changes in Avail assay
# get total read count per sample
tot_reads <- colSums(Avail_data[, 12:20])

# compute read frequency in each sample
Avail_data%>%
  mutate(freq_T0_R1 = Avail_T0_R1/tot_reads[1], 
         freq_T1_R3 = Avail_T1_R3/tot_reads[2],
         freq_T0_R2 = Avail_T0_R2/tot_reads[3],
         freq_T0_R3 = Avail_T0_R3/tot_reads[4], 
         freq_T1_R2 = Avail_T1_R2/tot_reads[5], 
         freq_T2_R1 = Avail_T2_R1/tot_reads[6], 
         freq_T2_R2 = Avail_T2_R2/tot_reads[7], 
         freq_T2_R3 = Avail_T2_R3/tot_reads[8],
         freq_T1_R1 = Avail_T1_R1/tot_reads[9]) -> Avail_data

# compute log2 fold change
Avail_data%>%
  mutate(delta_T0_T1_R1 = log2(freq_T1_R1) - log2(freq_T0_R1), 
         delta_T0_T2_R1 = log2(freq_T2_R1) - log2(freq_T0_R1), 
         delta_T0_T1_R2 = log2(freq_T1_R2) - log2(freq_T0_R2),
         delta_T0_T2_R2 = log2(freq_T2_R2) - log2(freq_T0_R2), 
         delta_T0_T1_R3 = log2(freq_T1_R3) - log2(freq_T0_R3), 
         delta_T0_T2_R3 = log2(freq_T2_R3) - log2(freq_T0_R3)) -> freq_Avail



# compute frequency changes PPI assay
# get total read count per sample
tot_reads <- colSums(PPI_data[, 12:20])
# compute read frequency in each sample
PPI_data%>%
  mutate(freq_T0_R1 = PPI_T0_R1/tot_reads[1], 
         freq_T1_R3 = PPI_T1_R3/tot_reads[2],
         freq_T0_R2 = PPI_T0_R2/tot_reads[3],
         freq_T0_R3 = PPI_T0_R3/tot_reads[4], 
         freq_T1_R2 = PPI_T1_R2/tot_reads[5], 
         freq_T2_R1 = PPI_T2_R1/tot_reads[6], 
         freq_T2_R2 = PPI_T2_R2/tot_reads[7], 
         freq_T2_R3 = PPI_T2_R3/tot_reads[8],
         freq_T1_R1 = PPI_T1_R1/tot_reads[9]) -> PPI_data

# compute log2 fold change
PPI_data%>%
  mutate(delta_T0_T1_R1 = log2(freq_T1_R1) - log2(freq_T0_R1), 
         delta_T0_T2_R1 = log2(freq_T2_R1) - log2(freq_T0_R1), 
         delta_T0_T1_R2 = log2(freq_T1_R2) - log2(freq_T0_R2),
         delta_T0_T2_R2 = log2(freq_T2_R2) - log2(freq_T0_R2), 
         delta_T0_T1_R3 = log2(freq_T1_R3) - log2(freq_T0_R3), 
         delta_T0_T2_R3 = log2(freq_T2_R3) - log2(freq_T0_R3)) -> freq_PPI

# pivot back to long with replicates

long_delta_avail <- 
pivot_longer(freq_Avail[, c(1:11, 30:35)], cols = c(12:17), names_to = c('calcul', 'timepoint_start', 'timepoint_end', 'replicate'),
             names_sep = '_', values_to = 'delta', )


long_delta_ppi <- 
  pivot_longer(freq_PPI[, c(1:11, 30:35)], cols = c(12:17), names_to = c('calcul', 'timepoint_start', 'timepoint_end', 'replicate'),
               names_sep = '_', values_to = 'delta', )

# associate generation number with log2 fold change
long_delta_avail <- 
merge(long_delta_avail, 
      gen_values[gen_values$assay == 'Avail', ], 
      by.x = c('replicate', 'timepoint_end'), 
      by.y = c('replicate', 'timepoint'))

# normalize by n_gen
long_delta_avail%<>%
  mutate(sel_coef = delta/n_gen)

# associate generation number with log2 fold change
long_delta_ppi <- 
  merge(long_delta_ppi, 
        gen_values[gen_values$assay == 'PPI', ], 
        by.x = c('replicate', 'timepoint_end'), 
        by.y = c('replicate', 'timepoint'))

# normalize by n_gen
long_delta_ppi%<>%
  mutate(sel_coef = delta/n_gen)

# total variants per assay 43080
nrow(long_delta_avail)/43080
# 0.9286908
nrow(long_delta_ppi)/43080
# 0.9551532

write_csv(long_delta_avail, '~/PL_projects/PL_papers/PPI_optimization_paper/data/sel_coef_availability.csv')
write_csv(long_delta_ppi, '~/PL_projects/PL_papers/PPI_optimization_paper/data/sel_coef_ppi.csv')

