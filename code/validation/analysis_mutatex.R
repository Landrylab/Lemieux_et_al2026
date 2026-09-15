# 

setwd('~/PL_projects/PL_papers/PPI_optimization_paper/data/model_docking/results/interface_ddgs/')

dir_list = dir(pattern="3*")

models_info <- matrix(unlist(strsplit(dir_list[1:15], '_')), nrow = 15, byrow = TRUE)

models_info<-
tibble::tibble(
PBD = gsub(models_info[, 1], pattern = '-w-pep', replacement = ''), 
rep = models_info[,2])


mut_list<-read_table('~/PL_projects/PL_papers/PPI_optimization_paper/data/model_docking/mutation_list.txt', col_names = 'AA_mut')


all_PBD<-tibble::tibble()

for(j in 1:length(dir_list[1:15])){

pbd<-
  unlist(strsplit(dir_list[j], '-'))[1] 

setwd(paste0('~/PL_projects/PL_papers/PPI_optimization_paper/data/model_docking/results/interface_ddgs/', dir_list[j], '/A-B/'))

file = dir(pattern="*")

all_ddG<-tibble::tibble()

for (i in 1:length(file)) {
 
   pos_mut<-unlist(strsplit(file[i], split = '*'))
  
    aa<-pos_mut[1]
    chain<-pos_mut[2]
    pos<-str_c(pos_mut[3:length(pos_mut)], collapse = '')
    
    
  ddG_table<-read_table(file = file[i],skip = 1, col_names = c('avg', 'sd', 'min', 'max'))

  bind_cols('chain' =chain, 
            'position' = pos, 
            'aa_ref' = aa, 
            mut_list, ddG_table)->ddG_table
  
  all_ddG<-
    bind_rows(all_ddG, ddG_table)
  
}

all_ddG<-bind_cols('PBD'=pbd, 'rep' = models_info[j, 'rep'], all_ddG)

all_PBD<-
  bind_rows(all_PBD, all_ddG)

}


source('~/PL_projects/PL_papers/PPI_optimization_paper/code/functions.R')

all_PBD$position<-as.numeric(all_PBD$position)

ddG_summary<-
data_summary(all_PBD, 'avg', c('PBD', 'chain', 'position', 'AA_mut'))
  
ggplot(ddG_summary)+
  geom_density(aes(sd), color = 'red')+
  geom_density(aes(avg))+
  xlim(-2.5, 2.5)


ddG_summary%>%
  dplyr::mutate(top_avg = case_when(avg > 2.5 ~2.5, 
                                    avg < -2.5 ~ -2.5, 
                                    (2.5>= avg) & (avg >=-2.5) ~ avg))->vis_ddG


write_csv(all_PBD, '~/PL_projects/PL_papers/PPI_optimization_paper/data/model_docking/ddG_interface.csv')


ggplot(vis_ddG)+
  facet_grid(cols = vars(chain), 
             rows = vars(PBD), scales = 'free', space = 'free')+
  geom_tile(aes(x = position, y = AA_mut, fill = top_avg))+
  #geom_point(aes(x = position, y = aa_ref), color = 'red', shape = 4, size = 1)+
  scale_fill_viridis_c(limits = c(-2.5, 2.5))+
  theme_classic2()+
  guides(fill = guide_colorbar())


PCA_score<-
read_csv('~/PL_projects/PL_papers/PPI_optimization_paper/data/PCA2/signif_score_welch20_filter.csv')


PCA_score%<>%
  dplyr::mutate(paired = case_when((pool == 'P4' & family == '363' | 
                                      pool == 'P5' & family == '366' | 
                                      pool == 'P6' & family == '385') ~ TRUE  , 
                                   (pool == 'P4' & family != '363' | 
                                      pool == 'P5' & family != '366' | 
                                      pool == 'P6' & family != '385') ~ FALSE), 
                double_mutant = case_when(is.na(max)~TRUE,
                                          !is.na(max)~ FALSE))


PCA_score%>%
  dplyr::filter(paired & condition == 'MTX' & !double_mutant)%>%
  dplyr::filter(!grepl('*', aa_seq, fixed = TRUE))->sub_signif

ddG_summary%>%
  dplyr::filter(chain == 'B')%>%
  right_join(sub_signif, join_by(position ==pos, AA_mut == deg, PBD == family))->comp_ddG_PCA


comp_ddG_PCA%>%
  dplyr::filter(peptide_max)%>%
  dplyr::select(PBD, med_norm)%>%
  unique()->ref_pep

comp_ddG_PCA%>%
ggplot()+
  facet_grid(vars(PBD), scales = 'free')+
  geom_point(aes(x = med_norm, y = avg, color = as.factor(position)))+
  geom_vline(data =ref_pep, aes(xintercept = med_norm))+
  stat_cor(aes(x = med_norm, y = avg), method = 'spearman')+
  scale_color_viridis_d()+
  #ylim(-2.5, 5)+
  xlab('PPI score')+
  ylab('ddG interface')+
  theme_classic2()







