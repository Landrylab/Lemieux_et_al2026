# AUTHOR : Pascale Lemieux
# compare degronopedia to peptide sequences screened with DHFR PCA

setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/')
source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')


degron<-
readxl::read_xlsx('./validation/DEGRONOPEDIA_degron_motifs.xlsx', sheet = 2)


# filter the degron database for yeast

degron%>%
  filter(grepl(pattern = 'S.cerevisiae', Organism))->degron_yeast

# import median scores of the July2025 experiment (single mutant availability)
specificity <- 
  read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/signif_score_welch20_filter.csv')

avail_screen2 <- subset(specificity, pool == 'P1' & condition == 'MTX')

# import the scores from the February2025 experiment (double mutant availability)
avail <- 
  read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA1/avail_signif_score.csv')

avail%<>%
  subset(family != 'reference')%>%
  mutate(family = as.numeric(family))


colnames(avail)[11]<-'med_norm'


# select unique aa sequences from the July 2025 screen
avail_screen2%>%
  select(aa_seq, family, condition, med_norm, med_sel, pos, deg, peptide_max, side)%>%
  subset(family != 'reference')%>%
  dplyr::mutate(family= as.numeric(family))%>%
  unique()-> avail_screen2_u

# combine both experiment for single mutants
all_avail_single <- 
  bind_rows(avail_screen2_u, 
            avail)

# compare degron to peptide sequences tested in availability

degron_scan<-tibble()
for (deg_seq in degron_yeast$Degron_regex) {
  

 pep_seq<-grep(deg_seq, all_avail_single$aa_seq)
 pep_seq<-all_avail_single$aa_seq[pep_seq]
 
 sub_deg<-
 tibble('aa_seq' = pep_seq, 
        'degron' = deg_seq)
 
 degron_scan<-
  bind_rows(degron_scan, sub_deg)
}

degron_scan%<>%
  left_join(degron_yeast, join_by(degron == 'Degron_regex'))

all_avail_single<-
  left_join(all_avail_single, degron_scan)

table(all_avail_single$Degron)

all_avail_single%>%
  filter(!is.na(degron))%>%
  select(aa_seq, family, degron, med_norm)%>%
  unique()->unique_deg

# 467/3466
all_avail_single$stop<-grepl(pattern = '*', all_avail_single$aa_seq, fixed = TRUE)

all_avail_single%>%
  filter(!stop)%>%
ggplot()+
  #facet_grid(cols = vars(Degron_location))+
  geom_jitter(data = all_avail_single[!is.na(all_avail_single$degron) & !all_avail_single$stop, ], 
             aes(x=as.factor(family), y=med_norm, color = Degron_location), width = 0.1, shape =1)+
  geom_violin(aes(as.factor(family), y=med_norm), fill ='transparent', quantile.colour = 'black', quantile.linetype = 'dashed')+
  ylab('Avail. score')+
  xlab('peptide family')+
  scale_color_manual(values = c('#a7ea52ff', '#ff8021ff'))+
  t+
  guides(color = guide_legend(title = 'degron type', position = 'bottom'))


ggsave('~/PL_projects/PL_papers/PPI_optimization_paper/figures/supplementary/FigS_degron.png', 
       width = 7, height = 4)

all_avail_single%>%
  filter(!stop)%>%
  filter(!is.na(Degron))%>%
  dplyr::group_by(family)%>%
  summarise(n(), med =median(med_norm))


all_avail_single%>%
  filter(!stop)%>%
  filter(!is.na(Degron))%>%
  filter(family==363)%>%
  dplyr::mutate(test = med_norm<0.0198)%>%
  select(test)%>%sum()
  










