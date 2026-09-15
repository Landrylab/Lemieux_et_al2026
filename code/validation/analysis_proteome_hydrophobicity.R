# AUTHOR : Pascale Lemieux
# compare the mVenus and D12 hydrophobicity vs the general hydrophobicity of cytosolic protein.
source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

# seq mVenus-linker-DHFR12
mVenusD12<-'MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKLICTTGKLPVPWPTLVTTLGYGVQCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGGVQLADHYQQNTPIGDGPVLLPDNHYLSYQSKLSKDPNEKRDHMVLLEFVTAAGITHGMDELYKGGGGSGGGGSMVRPLNCIVAVSQNMGIGKNGDYPWPPLRNESKYFQRMTTTSSVEGKQNLVIMGRKTWFSIPEKNRPLKDRINIVLSRELKEPPRGAHFLAKSLDDALRLIEQPELGT*'

# seq mVenus
mVenus<-'MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKLICTTGKLPVPWPTLVTTLGYGVQCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGGVQLADHYQQNTPIGDGPVLLPDNHYLSYQSKLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK'

# seq linker-D12
D12<-'GGGGSGGGGSMVRPLNCIVAVSQNMGIGKNGDYPWPPLRNESKYFQRMTTTSSVEGKQNLVIMGRKTWFSIPEKNRPLKDRINIVLSRELKEPPRGAHFLAKSLDDALRLIEQPELGT*'

# import  cornell hydrophobicity scale
hydrophobicity<-read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/hydrophobicity.csv')
colnames(hydrophobicity)[1]<-'AA'

sticky<-hydrophobicity%>%
  dplyr::select(AA, Cornette)

colnames(sticky)[2]<-'sticky_score'

aa_seq<-c(mVenusD12, mVenus, D12)

sticky_seq<-tibble('id' = c('mVenusD12', 'mVenus', 'D12'), aa_seq)

# modify the aa sequence based on the stop codon
aa_seq <- 
  strsplit(aa_seq, '*', fixed = T)

aa_seq <- 
  lapply(aa_seq, `[[`, 1)
aa_seq <- 
  strsplit(unlist(aa_seq), '')


for (i in 1:length(aa_seq)) {
  if (length(aa_seq[[i]]) == 0){
    sticky_seq[i, 'sticky_score'] <- 0
    next
  }
  sticky_seq[i, 'sticky_score_abs'] <- get_stickyness(aa_seq[[i]], sticky, absolute = TRUE)
  sticky_seq[i, 'sticky_score_rel'] <- get_stickyness(aa_seq[[i]], sticky, absolute = FALSE)
}

# compute stickiness absolute and relative for all protein of S cerevisiae

library(seqinr)

proteome_yeast<-
read.fasta('~/PL_projects/PL_papers/PPI_optimization_paper/data/orf_trans.fasta')

sticky$AA<-
  tolower(sticky$AA)

aa_proteome<-
  unlist(lapply(proteome_yeast, str_c, collapse = ''))

proteome_sticky<-tibble('gene'=names(proteome_yeast), 'aa_seq' = aa_proteome)

for (i in 1:length(aa_proteome)) {
  if (length(proteome_yeast[[i]]) == 0){
    proteome_sticky[i, 'sticky_score'] <- 0
    next
  }
  proteome_sticky[i, 'sticky_score_abs'] <- get_stickyness(proteome_yeast[[i]], sticky, absolute = TRUE)
  proteome_sticky[i, 'sticky_score_rel'] <- get_stickyness(proteome_yeast[[i]], sticky, absolute = FALSE)
}

write_csv(proteome_sticky, '~/PL_projects/PL_papers/PPI_optimization_paper/data/hydrophobicity_proteome.csv')

proteome_sticky<-
  read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/hydrophobicity_proteome.csv')

#relative stickiness
ggplot()+
  geom_density(data = proteome_sticky, 
                 aes(x = sticky_score_rel), bins = 50, fill= '#373737')+
  geom_vline(data = sticky_seq[sticky_seq$id =='mVenusD12',], 
             aes(xintercept = sticky_score_rel), color = '#fde333ff')+
  t+
  xlab('yeast proteome\nrelative hydrophobicity')+
  theme(legend.position = 'bottom')->YFP_sticky

ggsave(YFP_sticky, filename='~/PL_projects/PL_papers/PPI_optimization_paper/figures/validation/YFP_hydro.svg', 
       height = 4, width =7)

# absolute stickiness
ggplot()+
  geom_histogram(data = proteome_sticky, 
                 aes(x = sticky_score_rel), bins = 50)+
  geom_vline(data = sticky_seq, 
             aes(xintercept = sticky_score_rel, color = id))+
  t+
  theme(legend.position = 'bottom')




