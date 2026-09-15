# script to prepare the codon optimized dna sequence for ordering
# add homology region to the destination vector to peptides and PBDs dna sequences
# using the output of Codon_optimization.ipynb
library(tidyverse)
library(pkgmaker)

##################################
##     Variables to modify      ##
# directory of the repository with the output files of Codon_optimization.ipynb
# and with the output file of peptide_generation.R
path = '~/PL_projects/PL_papers/PPI_optimization_paper/code/library_design/'
output_dir = './'

## Homology sequences to add ##
# if amplification with primers
# Use the peptide sequence + homology as A_F & A_R is Cyct_fusion (atgcgtacacgcgtctgtac)
# with Tm = 57
D3linker_homo = 'GGAGGTGGAGCTAGC'
D3cyct_homo = 'TAAaagcttattagttatgt'


# Amplification of PBD GeneFragment or gBlock
# Use A_F : aggttttgggacgctcgaag & A_R : GGTCAGGTGGAGGCGGATCC
# with Tm = 57
D12linker_homo = 'AGGCGGGTCTGGAGGAGGTGGGTCAGGTGGAGGCGGATCC'
D12cyct_homo = 'TAAactagtcaaattaaagccttcgagcgtcccaaaacct'


# choice for the first few PBD to order based on their Primary ID
PBD_order <- paste0('PRM_0', c('250', '304', '131', '299', '152', '366', '385', '025', '061', '246', '214', '363'))

# sequence dividing the PBD sequences (ordered in pair) in the gBlock 
inter_seq <- 'GTGGATCCAC'

##################################

setwd(path)

# read the codon optimized sequences of interest
pep <- read_csv(paste0(path, 'peptide_deg_DNA.csv'))
domain <- read_csv(paste0(path, 'domain_DNA.csv'))

# read the table with the library information
ref <- read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/code/library_design/library_09102024.csv')

# merge tables together to have one ligne per domain-peptide combination

ref <- 
merge(ref,
  pep[, c(2,3)], 
  by.x = 'peptide_sequence', 
  by.y = 'aa_seq')

ref <- 
  merge(ref, 
        domain[, c(2,3)], 
        by.x = 'Domain_Sequence', 
        by.y = 'aa_seq', 
        suffixes = c('.pep', '.domain'))


table(ref[, c('Primary_ID', 'method')])



# add homology sequences to the dna sequence of the peptides and domains

ref$dna_seq.pep_homology <- str_c(D3linker_homo, ref$dna_seq, D3cyct_homo)
ref$dna_seq.domain_homology <- str_c(D12linker_homo, ref$dna_seq, D12cyct_homo)

# write file with information associated to the peptide-PBD pairs and homology region added to the dna sequences 
write_csv(ref, paste0(path, 'homology_library22042025.csv'))


# We order each domain in pair to reduce synthesis cost. 
# To do so, we verify that the ends (20pb) of each sequence has at least 
# 10 pb of differences. This way, we make sure that the PCR amplification
# from the gene fragments is specific.

to_order <- unique(ref$dna_seq.domain)
split_seq <- strsplit( to_order, split = '')

# select the first and last 20pb of each PBD dna sequence
l_bar = 20
barcode_start <- vector(mode = 'character', length = 27)
barcode_end <- vector(mode = 'character', length = 27)

for(i in 1:length(split_seq)){
  barcode_start[i] <- str_c(split_seq[[i]][1:l_bar], collapse = '')
  end <- length(split_seq[[i]])
  n_end <- end - l_bar
  barcode_end[i] <- str_c(split_seq[[i]][n_end:end], collapse = '')
}

to_order <- 
  bind_cols('seq' = to_order, 'barcode_start' = barcode_start, 'barcode_end' = barcode_end, id = unique(ref$Primary_ID))

# create the object to save the comparison information
identity_start <- matrix(nrow = length(barcode_end), ncol = length(barcode_end)) 
identity_end <- matrix(nrow = length(barcode_end), ncol = length(barcode_end)) 

# Loop through the starts and ends to compare if they have
# more than 10 pb od difference
for(i in 1:nrow(to_order)){
  # starts
  seq_barcode <- to_order$barcode_start[i]
  x <- lapply(to_order$barcode_start, str_diff, y=seq_barcode)
  x <- lapply(x, length)
  identity_start[, i] <-unlist(x)
  
  # ends
  seq_barcode <- to_order$barcode_end[i]
  x <- lapply(to_order$barcode_end, str_diff, y=seq_barcode)
  x <- lapply(x, length)
  identity_end[, i] <-unlist(x)
}

# Boolean matrix, TRUE = more than 10pb of difference between starts and ends 
# of the 2 sequences compared 
diff_10 <- (identity_end>=10 & identity_start>=10)
row.names(diff_10) <- to_order$id
colnames(diff_10)<- to_order$id
diff_10

to_order$dna_lenght <-  nchar(unlist(to_order[, 'seq']))


# choice for the first few PBD to order
diff_10[rownames(diff_10) %in% PBD_order, colnames(diff_10) %in% PBD_order]


# get the %gc for  
get_gc <- 
function(seq_str){
  gc <- sum(table(seq_str)[c(2, 4)])
  all <- sum(table(seq_str)[c(1:4)])
  per <- gc/all*100
  return(per)
}

seq_list <- strsplit(domain$dna_seq, '')
lapply(seq_list, get_gc)


# manually identify pairs that could be order together use the dna lenght to order <500pb gene fragment
pairs <-
  list(c('PRM_0366', 'PRM_0025'), 
       c('PRM_0385', 'PRM_0061'), 
       c('PRM_0363', 'PRM_0246'), 
       c('PRM_0214', 'PRM_0131'), 
       c('PRM_0250', 'PRM_0304'), 
       c('PRM_0152', 'PRM_0299'))

frequency(strsplit(domain$dna_seq, ''))

to_order$dna_lenght <-  nchar(unlist(to_order[, 'seq']))
to_order[to_order$id %in% unlist(pairs), c('id', 'dna_lenght')]

# combine the pairs in one sequence and write ordering file
PBD_order <- tibble()
for (i in 1:length(pairs)){
  PBD_order[i, 'name'] <- 
   str_c(gsub(pattern = 'PRM_', replacement = '', pairs[[i]]), collapse = '_')
  
  s1 <- unique(ref[ref$Primary_ID == pairs[[i]][1], 'dna_seq.domain'])
  s2 <-  unique(ref[ref$Primary_ID == pairs[[i]][2], 'dna_seq.domain'])
  PBD_order[i, 'dna_sequence'] <- 
    str_c(c(s1, s2), collapse = inter_seq)
}

seq_list <- strsplit(PBD_order$dna_sequence, '')
lapply(seq_list, get_gc)


write_csv(PBD_order, paste0(output_dir, 'PBD_order09102024.csv'))

# generate primer to amplify the selected PBD dna sequence

# don't forget to reverse complement the end 
D12linker_homo = tolower('AGGCGGGTCTGGAGGAGGTGGGTCAGGTGGAGGCGGATCC')
D12cyct_homo = tolower('aggttttgggacgctcgaaggctttaatttgactagtTTA')

library(seqinr)
to_order$O_F <- str_c(D12linker_homo, to_order$barcode_start)

O_R <- strsplit(to_order$barcode_end, '')
O_R <- lapply(O_R, comp, forceToLower = FALSE)
O_R <- lapply(O_R, rev)
O_R <- lapply(O_R, str_c, collapse = '')

to_order$O_R <- str_c(D12cyct_homo, unlist(O_R))
sub_order <- to_order[to_order$id %in% unlist(pairs), ]
sub_order$id <- 
factor(sub_order$id, 
       levels = unlist(pairs), 
       labels = unlist(pairs))

oligo <- 
pivot_longer(sub_order[, -c(1:3)], 
             cols = c('O_F', 'O_R'), 
             names_to = 'O', 
             values_to = 'sequence')

oligo$name <- paste(oligo$id, oligo$O, 'for insertion downstream of D12 linker of pGD110')

oligo$id<- 
  factor(oligo$id, 
         levels = unlist(pairs), 
         labels = unlist(pairs))

oligo <- 
oligo[order(oligo$id), ]

# prepare oligos for peptide cloning
oligo_peptide <- ref[ref$Primary_ID %in% unlist(pairs),  c('Primary_ID', 'dna_seq.pep_homology', 'method')]

# order sequences by method and by domain id

oligo_peptide$Primary_ID<- 
  factor(oligo_peptide$Primary_ID, 
         levels = unlist(pairs), 
         labels = unlist(pairs))

colnames(oligo_peptide)[2] <- 'sequence'

exp_oligo <- oligo_peptide[oligo_peptide$method %in% c('NGS', 'HAL', 'PWM_top'),]
deg_oligo <- oligo_peptide[oligo_peptide$method == 'deg',]
exp_oligo<- exp_oligo[order(exp_oligo$Primary_ID), ]

# create unique name for each oligo
name_oligo <- table(exp_oligo$Primary_ID)
id <- names(name_oligo)
x <- vector(mode = 'character')

for (i in 1:length(id)){
  x <- c(x, paste0(id[i],'_',1:name_oligo[i]))
}
exp_oligo$name <- x

write_csv(exp_oligo, 'id_peptide.csv')

exp_oligo <- exp_oligo[order(exp_oligo$Primary_ID), ]

exp_oligo$name <- 
  str_c('O_F for peptide ', exp_oligo$name, ' insertion downstream of D3 linker of pGD110')

deg_oligo$name <- 
  str_c('O_F for degenerate peptide ', deg_oligo$Primary_ID, ' insertion downstream of D3 linker of pGD110')
deg_oligo <- deg_oligo[order(deg_oligo$Primary_ID), ]

# combine oligo for domain amplification and peptide insertion

oligo <- 
bind_rows(oligo[, c('sequence', 'name')], 
          exp_oligo[ ,c('sequence', 'name')], 
          deg_oligo[ ,c('sequence', 'name')])

# write file with all oligos
write_csv(oligo, 
          'oligo_order16102024.csv')

sum(nchar(exp_oligo$sequence)>=90)
sum(nchar(deg_oligo$sequence)>=90)
