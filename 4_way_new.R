
# libraries ---------------------------------------------------------------


# load libraries
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(PFun)
library(forcats)
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(gtools)
library(broom)
library(here)
library(openxlsx)
library(coin)
library(RVAideMemoire)
library(colorspace)
library(plotly)
library(ggpubr)
library(patchwork)
library(here)
library(cowplot)

# options(width = 150)

# 4-way plot --------------------------------------------------------------

library(gridExtra)
adjustments = as_tibble(read.csv('Contrast_adjustments.csv'))

###
### 4-way plot
###

# is a hell of a code, it's super messy, but it works...

worm.sum[worm.sum$Strain == 'CRP', ]$Strain = 'crp'
worm.sum$Strain = as.factor(worm.sum$Strain)



# print the list of variables 
for (name in unique(worm.sum$Strain)){
  print(name)
  print(unique((worm.sum %>% filter(Strain == name))$Experiment))
}

# select only one experiment and strain
# contrasts and other main variables of this part

# specify labels to plot 
labels_4way = c("Uracil_N","Cytidine_N","Uridine_N","Cytidine-3'-monophosphate",
                "Cytidine- 2',3'-cyclic monophosphate","Uridine-3'-monophosphate",
                "Uridine-2',3'-cyclic-monophosphate","Cytidine-2'-monophosphate",
                "Cytidine-5'-monophosphate","Uridine-2'-monophosphate",
                "Uridine-5'-monophosphate","Uridine","D-Galactose",
                "Glycerol","D-Sorbitol","D-Trehalose","Dulcitol",
                "Maltose","D-Mannose",'alpha-D-Glucose')

labels_4way_nogluc = c("Uracil_N","Cytidine_N","Uridine_N","Cytidine-3'-monophosphate",
                "Cytidine- 2',3'-cyclic monophosphate","Uridine-3'-monophosphate",
                "Uridine-2',3'-cyclic-monophosphate","Cytidine-2'-monophosphate",
                "Cytidine-5'-monophosphate","Uridine-2'-monophosphate",
                "Uridine-5'-monophosphate","Uridine","D-Galactose",
                "Glycerol","D-Sorbitol","D-Trehalose","Dulcitol",
                "Maltose","D-Mannose")

# specify nucleotide labels
labels_nuc = c("Uracil_N","Cytidine_N","Uridine_N","Cytidine-3'-monophosphate",
               "Cytidine- 2',3'-cyclic monophosphate","Uridine-3'-monophosphate",
               "Uridine-2',3'-cyclic-monophosphate","Cytidine-2'-monophosphate",
               "Cytidine-5'-monophosphate","Uridine-2'-monophosphate",
               "Uridine-5'-monophosphate","Uridine")
# specify sugar labels
labels_sug = c("D-Galactose","Glycerol","D-Sorbitol",
               "D-Trehalose","Dulcitol","Maltose","D-Mannose")



# BW ----------------------------------------------------------------------




# select only one experiment and strain
# contrasts and other main variables of this part

strain = 'BW'
experiment = "T5"


str_C = paste(strain, '_C', sep = '')
str_T = paste(strain, '_T', sep = '')
str_CT = paste(strain, '_5FU', sep = '')
alln = paste(strain, '_5FU_Alln', sep = '')
adjs = paste(strain, '_5FU_adj', sep = '')


# for example, to mimic Pov's code
worm.sum %>% filter(Experiment == experiment, Strain == strain)

worm.sum %>% 
  write_csv(here('Summary', 'worm_dev_scores.csv'))

worm.old = worm.sum %>%
  filter(Experiment == experiment & Strain == strain) %>% 
  ungroup %>%
  mutate(Description = 'C. elegans development at 5uM 5-FU',
         Contrast = 'Ce_Dev5',   # C. elegans development with 5FU
         Contrast_type = 'Treatment',
         FDR = NA) %>%
  select(Description, Contrast, Contrast_type, Plate, Well, logFC = Median, SE = SD, FDR)

# check that we filtered stuff properly
worm.old %>%
  group_by(Plate) %>%
  summarise(N = n())


# Join bacterial and worm results (change contrasts)
jointresults = allresults$results %>%
  filter((Contrast %in% c(str_C, str_T, str_CT, alln))) %>%
  mutate(Contrast = fct_recode(Contrast, adjs = alln)) %>% 
  select(Description, Contrast, Contrast_type, Plate, Well, logFC, SE, FDR) %>%
  bind_rows(worm.old) %>%
  left_join(info) %>%
  select(Description:Well,Index:KEGG_ID,logFC:FDR)


jointcast = jointresults %>%
  select(Contrast, Plate, Well, Index, Metabolite, MetaboliteU, EcoCycID, KEGG_ID, logFC, SE, FDR) %>%
  gather(Stat, Value, logFC, SE, FDR) %>%
  unite(CS, Contrast, Stat) %>%
  spread(CS, Value) %>%
  rename(Ce_Dev5_Median = Ce_Dev5_logFC, Ce_Dev5_SD = Ce_Dev5_SE) %>%
  select(-Ce_Dev5_FDR)


# write.csv(jointresults,paste(odir,'/Ecoli_results_All_As_old_screen.csv', sep = ''),row.names = FALSE)
# write.csv(jointcast,paste(odir,'/Ecoli_results_sidebyside_As_old_screen.csv', sep = ''),row.names = FALSE)


# Multiplex 
jointresults.multi = multiplex(jointresults,c("Plate","Well","Index","Metabolite","MetaboliteU","EcoCycID","KEGG_ID"))
jointresults.multi2 = multiplex(jointresults,c("Plate","Well","Index","Metabolite","MetaboliteU","EcoCycID","KEGG_ID"),2)




########################
NGMa = adjustments[adjustments$Contrast == alln, ]$a
NGMb = adjustments[adjustments$Contrast == alln, ]$b



# gradcolours = c('#71B83B','yellow','orange','red')
gradcolours = c('red','orange', 'yellow','#71B83B') ## new colors for the paper

# gradcolours = c("#9A1324", "#D53E24", "#F3E930", "#74B533")

total = jointresults.multi %>%
  filter(x_Contrast == str_C & y_Contrast == str_T & z_Contrast == 'Ce_Dev5')

total = total %>% 
  mutate(score = (abs(y_logFC - (x_logFC-NGMa))*z_logFC),
         score2 = abs(y_logFC - (x_logFC-NGMa))) %>% 
  arrange(score) %>% 
  mutate(Index = factor(Index, levels = Index)) %>% 
  drop_na(z_logFC)

total = total %>% 
  mutate(new_labels = case_when(MetaboliteU %in% labels_4way ~ MetaboliteU,
                                TRUE ~ ''),
         label_size = case_when(z_logFC == 4 ~ 3,
                                TRUE ~ 2),
         label_color = case_when(MetaboliteU %in% labels_nuc ~ 'Nucleotide',
                                 MetaboliteU %in%  labels_sug ~ 'Sugars',
                                 TRUE ~ '')
  ) 


BW.scores = total %>% 
  select(Plate, Well, MetaboliteU, x_Contrast, z_logFC, score, score2)


### BW plot ####
BW.plot = total %>% 
  select(-x_Description:-x_FDR,
         -y_Description:-y_FDR,
         -z_Description:-z_Contrast_type) %>% 
  ggplot(aes(x = reorder(MetaboliteU, score), y = score)) +
  geom_point(aes(fill = z_logFC, size = z_logFC),pch=21) +
  geom_hline(aes(yintercept = 1), alpha = 0.9, color = 'grey') +
  scale_fill_gradientn(colours = gradcolours,
                       breaks = c(1,2,3,4), limits = c(1,4), 
                       guide = "legend", name = 'C. elegans\nphenotype') +
  scale_color_manual(values = c('black',
                                '#C70B00', # nucleotides 
                                '#310CB3' # sugars
  )) +
  geom_text_repel(aes(x = MetaboliteU, y = score, 
                      label= new_labels,
                      size = label_size,
                      color = factor(label_color)),
                  box.padding = unit(0.6, 'lines'),
                  segment.alpha = 0.4,
                  max.overlaps = 50) +
  # geom_text_repel(aes(x = MetaboliteU, y = score, label = ifelse(z_logFC >= 4, MetaboliteU, '')), 
  #                 box.padding = unit(0.6, "lines"), segment.alpha = 0.4) + 
  labs(x = 'Metabolite',
       y = 'Composed score',
       title = 'BW25113 - Wild type') +
  coord_cartesian(xlim = c(0,400)) +
  # ylim(0, 15) +
  scale_size(guide='none') +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size=18),
        legend.text = element_text(size=18)
        ) +
  guides(color = 'none',fill = guide_legend(override.aes = list(size=8)))


names(total)

total %>% 
  select(-x_Description:-x_FDR,
         -y_Description:-y_FDR,
         -z_Description:-z_Contrast_type, 
         -z_SE, -z_FDR) %>% 
  rename(worm_score = z_logFC) %>% 
  write_csv(here('Summary', '4way_BW.csv'))




# pyrE --------------------------------------------------------------------


strain = 'pyrE'
experiment = "T5"


str_C = paste(strain, '_C', sep = '')
str_T = paste(strain, '_T', sep = '')
str_CT = paste(strain, '_5FU', sep = '')
alln = paste(strain, '_5FU_Alln', sep = '')
adjs = paste(strain, '_5FU_adj', sep = '')


# for example, to mimic Pov's code
worm.sum %>% filter(Experiment == experiment, Strain == strain)

worm.old = worm.sum %>%
  filter(Experiment == experiment & Strain == strain) %>% 
  ungroup %>%
  mutate(Description = 'C. elegans development at 5uM 5-FU',
         Contrast = 'Ce_Dev5',   # C. elegans developement with 5FU
         Contrast_type = 'Treatment',
         FDR = NA) %>%
  select(Description, Contrast, Contrast_type, Plate, Well, logFC = Median, SE = SD, FDR)

# check that we filtered stuff properly
worm.old %>%
  group_by(Plate) %>%
  summarise(N = n())


# Join bacterial and worm results (change contrasts)
jointresults = allresults$results %>%
  filter((Contrast %in% c(str_C, str_T, str_CT, alln))) %>%
  mutate(Contrast = fct_recode(Contrast, adjs = alln)) %>% 
  select(Description, Contrast, Contrast_type, Plate, Well, logFC, SE, FDR) %>%
  bind_rows(worm.old) %>%
  left_join(info) %>%
  select(Description:Well,Index:KEGG_ID,logFC:FDR)


jointcast = jointresults %>%
  select(Contrast, Plate, Well, Index, Metabolite, MetaboliteU, EcoCycID, KEGG_ID, logFC, SE, FDR) %>%
  gather(Stat, Value, logFC, SE, FDR) %>%
  unite(CS, Contrast, Stat) %>%
  spread(CS, Value) %>%
  rename(Ce_Dev5_Median = Ce_Dev5_logFC, Ce_Dev5_SD = Ce_Dev5_SE) %>%
  select(-Ce_Dev5_FDR)


# write.csv(jointresults,paste(odir,'/Ecoli_results_All_As_old_screen.csv', sep = ''),row.names = FALSE)
# write.csv(jointcast,paste(odir,'/Ecoli_results_sidebyside_As_old_screen.csv', sep = ''),row.names = FALSE)


# Multiplex 
jointresults.multi = multiplex(jointresults,c("Plate","Well","Index","Metabolite","MetaboliteU","EcoCycID","KEGG_ID"))
jointresults.multi2 = multiplex(jointresults,c("Plate","Well","Index","Metabolite","MetaboliteU","EcoCycID","KEGG_ID"),2)




########################
NGMa = adjustments[adjustments$Contrast == alln, ]$a
NGMb = adjustments[adjustments$Contrast == alln, ]$b
r2 = adjustments[adjustments$Contrast == alln, ]$r2

sr = sqrt((1-r2)/(379-2))

gradcolours = c('#71B83B','yellow','orange','red')





total = jointresults.multi %>%
  filter(x_Contrast == str_C & y_Contrast == str_T & z_Contrast == 'Ce_Dev5')

total = total %>% 
  mutate(score = (abs(y_logFC - (x_logFC-NGMa))*z_logFC),
         score2 = abs(y_logFC - (x_logFC-NGMa))) %>% 
  arrange(score) %>% 
  mutate(Index = factor(Index, levels = Index)) %>% 
  drop_na(z_logFC)

total = total %>% 
  mutate(new_labels = case_when(MetaboliteU %in% labels_4way ~ MetaboliteU,
                                TRUE ~ ''),
         label_size = case_when(z_logFC == 4 ~ 3,
                                TRUE ~ 2),
         label_color = case_when(MetaboliteU %in% labels_nuc ~ 'Nucleotide',
                                 MetaboliteU %in%  labels_sug ~ 'Sugars',
                                 TRUE ~ '')
         ) 


pyrE.scores = total %>% 
  select(Plate, Well, MetaboliteU, x_Contrast, z_logFC, score, score2)

### pyrE plot ####
pyrE.plot = total %>% 
  ggplot(aes(x = reorder(MetaboliteU, score), y = score)) +
  geom_point(aes(fill = z_logFC, size = z_logFC),pch=21) +
  geom_hline(aes(yintercept = 1), alpha = 0.9, color = 'grey') +
  scale_fill_gradientn(colours = gradcolours,
                         breaks = c(1,2,3,4), limits = c(1,4), 
                         guide = "legend", name = 'C. elegans\nphenotype') +
  scale_color_manual(values = c('black',
                                '#C70B00', # nucleotides 
                                '#310CB3' # sugars
                                )) + 
  geom_text_repel(aes(x = MetaboliteU, y = score, 
                      label= new_labels, 
                      size = label_size,
                      color = factor(label_color)),
                  box.padding = unit(0.6, 'lines'),
                  segment.alpha = 0.4,
                  max.overlaps = 50) +
  # geom_text_repel(aes(x = MetaboliteU, y = score, 
  #                     label = ifelse(z_logFC >= 4, MetaboliteU, '')), 
  #                 box.padding = unit(0.6, "lines"), 
  #                 segment.alpha = 0.4,
  #                 max.overlaps = 50) +
  labs(x = 'Metabolite',
       y = 'Composed score',
       title = bquote(~Delta*pyrE ~ "mutant")) +
  coord_cartesian(xlim = c(0,400)) +
  # ylim(0, 15) +
  scale_size(guide='none') +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(size=16, face="bold", color ='black'),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size=18),
        legend.text = element_text(size=18)) +
  guides(color = 'none',
         fill = guide_legend(override.aes = list(size=8)))


total %>% 
  select(-x_Description:-x_FDR,
         -y_Description:-y_FDR,
         -z_Description:-z_Contrast_type, 
         -z_SE, -z_FDR) %>% 
  rename(worm_score = z_logFC) %>% 
  write_csv(here('Summary', '4way_pyrE.csv'))





# TM ----------------------------------------------------------------------

# select only one experiment and strain
# contrasts and other main variables of this part

strain = 'TM'
experiment = "T250"


str_C = paste(strain, '_C', sep = '')
str_T = paste(strain, '_T', sep = '')
str_CT = paste(strain, '_5FU', sep = '')
alln = paste(strain, '_5FU_Alln', sep = '')
adjs = paste(strain, '_5FU_adj', sep = '')


# for example, to mimic Pov's code
worm.sum %>% filter(Experiment == experiment, Strain == strain)

worm.old = worm.sum %>%
  filter(Experiment == experiment & Strain == strain) %>% 
  ungroup %>%
  mutate(Description = 'C. elegans development at 5uM 5-FU',
         Contrast = 'Ce_Dev5',   # C. elegans developement with 5FU
         Contrast_type = 'Treatment',
         FDR = NA) %>%
  select(Description, Contrast, Contrast_type, Plate, Well, logFC = Median, SE = SD, FDR)

# check that we filtered stuff properly
worm.old %>%
  group_by(Plate) %>%
  summarise(N = n())


# Join bacterial and worm results (change contrasts)
jointresults = allresults$results %>%
  filter((Contrast %in% c(str_C, str_T, str_CT, alln))) %>%
  mutate(Contrast = fct_recode(Contrast, adjs = alln)) %>% 
  select(Description, Contrast, Contrast_type, Plate, Well, logFC, SE, FDR) %>%
  bind_rows(worm.old) %>%
  left_join(info) %>%
  select(Description:Well,Index:KEGG_ID,logFC:FDR)


jointcast = jointresults %>%
  select(Contrast, Plate, Well, Index, Metabolite, MetaboliteU, EcoCycID, KEGG_ID, logFC, SE, FDR) %>%
  gather(Stat, Value, logFC, SE, FDR) %>%
  unite(CS, Contrast, Stat) %>%
  spread(CS, Value) %>%
  rename(Ce_Dev5_Median = Ce_Dev5_logFC, Ce_Dev5_SD = Ce_Dev5_SE) %>%
  select(-Ce_Dev5_FDR)


# write.csv(jointresults,paste(odir,'/Ecoli_results_All_As_old_screen.csv', sep = ''),row.names = FALSE)
# write.csv(jointcast,paste(odir,'/Ecoli_results_sidebyside_As_old_screen.csv', sep = ''),row.names = FALSE)


# Multiplex 
jointresults.multi = multiplex(jointresults,c("Plate","Well","Index","Metabolite","MetaboliteU","EcoCycID","KEGG_ID"))
jointresults.multi2 = multiplex(jointresults,c("Plate","Well","Index","Metabolite","MetaboliteU","EcoCycID","KEGG_ID"),2)




########################
NGMa = adjustments[adjustments$Contrast == alln, ]$a
NGMb = adjustments[adjustments$Contrast == alln, ]$b
r2 = adjustments[adjustments$Contrast == alln, ]$r2

# sr = sqrt((1-r2)/(379-2))

sr = 1/r2

gradcolours = c('#71B83B','yellow','orange','red')




total = jointresults.multi %>%
  filter(x_Contrast == str_C & y_Contrast == str_T & z_Contrast == 'Ce_Dev5')

total = total %>% 
  mutate(score = (abs(y_logFC - (x_logFC-NGMa))*z_logFC),
         score2 = abs(y_logFC - (x_logFC-NGMa))) %>% 
  arrange(score) %>% 
  mutate(Index = factor(Index, levels = Index)) %>% 
  drop_na(z_logFC)

total = total %>% 
  mutate(new_labels = case_when(MetaboliteU %in% labels_4way ~ MetaboliteU,
                                TRUE ~ ''),
         label_size = case_when(z_logFC >= 3.5 ~ 3,
                                TRUE ~ 2),
         label_color = case_when(MetaboliteU %in% labels_nuc ~ 'Nucleotide',
                                 MetaboliteU %in%  labels_sug ~ 'Sugars',
                                 TRUE ~ '')
  ) 

TM.scores = total %>% 
  select(Plate, Well, MetaboliteU, x_Contrast, z_logFC, score, score2)


### TM plot ####
TM.plot = total %>% 
  ggplot(aes(x = reorder(MetaboliteU, score), y = score)) +
  geom_point(aes(fill = z_logFC, size = z_logFC),pch=21) +
  # annotate("rect", xmin = 0, xmax = 400, ymin = 1 - sr, ymax = 1 + sr, alpha = .2) +
  geom_hline(aes(yintercept = 1), alpha = 0.9, color = 'grey') +
  scale_fill_gradientn(colours = gradcolours,
                       breaks = c(1,2,3,4), limits = c(1,4), 
                       guide = "legend", name = 'C. elegans\nphenotype') +
  scale_color_manual(values = c('black',
                                '#C70B00', # nucleotides 
                                '#310CB3' # sugars
  )) +
  geom_text_repel(aes(x = MetaboliteU, y = score, 
                      label= new_labels,
                      size = label_size,
                      color = factor(label_color)),
                  box.padding = unit(0.6, 'lines'),
                  segment.alpha = 0.4,
                  max.overlaps = 800) +
  # geom_text_repel(aes(x = MetaboliteU, y = score, label = ifelse(z_logFC >= 4, MetaboliteU, '')), 
  #                 box.padding = unit(0.6, "lines"), segment.alpha = 0.4) + 
  labs(x = 'Metabolite',
       y = 'Composed score',
       title = bquote(~Delta*upp*Delta*udp*Delta*udk ~ "mutant")) +
  coord_cartesian(xlim = c(0,400)) +
  # ylim(0, 15) +
  scale_size(guide='none') +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_text(size=16, face="bold", color ='black'),
        plot.title = element_text(size=18),
        legend.text = element_text(size=18)
        ) +
  guides(color = 'none',
         fill = guide_legend(override.aes = list(size=8)))



total %>% 
  select(-x_Description:-x_FDR,
         -y_Description:-y_FDR,
         -z_Description:-z_Contrast_type, 
         -z_SE, -z_FDR) %>% 
  rename(worm_score = z_logFC) %>% 
  write_csv(here('Summary', '4way_TM.csv'))


# 
# 
# 
# library(ggpubr)
# ggarrange(BW.plot, pyrE.plot, TM.plot, nrow = 3)
# p4 = ggarrange(BW.plot, pyrE.plot, TM.plot, nrow = 3)
# 
# 
# ggsave(p4, file = here('Summary', '4wayScreening_score.pdf'),
#        height = 10, width = 13)
# 



p4 = BW.plot / pyrE.plot / TM.plot +
  plot_layout(guides = 'collect')

p4

ggsave(p4, file = here('Summary', '4wayScreening_score_v2.pdf'),
       height = 10, width = 13)










# PAPER version ------------------------------------------------------


met_dat = function(strain='BW',experiment='T5') {

    str_C = paste(strain, '_C', sep = '')
    str_T = paste(strain, '_T', sep = '')
    str_CT = paste(strain, '_5FU', sep = '')
    alln = paste(strain, '_5FU_Alln', sep = '')
    adjs = paste(strain, '_5FU_adj', sep = '')
    
    
    # for example, to mimic Pov's code
    worm.sum %>% filter(Experiment == experiment, Strain == strain)
    
    worm.old = worm.sum %>%
      filter(Experiment == experiment & Strain == strain) %>% 
      ungroup %>%
      mutate(Description = 'C. elegans development at 5uM 5-FU',
             Contrast = 'Ce_Dev5',   # C. elegans developement with 5FU
             Contrast_type = 'Treatment',
             FDR = NA) %>%
      select(Description, Contrast, Contrast_type, Plate, Well, logFC = Median, SE = SD, FDR)
    
    # Join bacterial and worm results (change contrasts)
    jointresults = allresults$results %>%
      filter((Contrast %in% c(str_C, str_T, str_CT, alln))) %>%
      mutate(Contrast = fct_recode(Contrast, adjs = alln)) %>% 
      select(Description, Contrast, Contrast_type, Plate, Well, logFC, SE, FDR) %>%
      bind_rows(worm.old) %>%
      left_join(info) %>%
      select(Description:Well,Index:KEGG_ID,logFC:FDR)
    
    
    jointcast = jointresults %>%
      select(Contrast, Plate, Well, Index, Metabolite, MetaboliteU, EcoCycID, KEGG_ID, logFC, SE, FDR) %>%
      gather(Stat, Value, logFC, SE, FDR) %>%
      unite(CS, Contrast, Stat) %>%
      spread(CS, Value) %>%
      rename(Ce_Dev5_Median = Ce_Dev5_logFC, Ce_Dev5_SD = Ce_Dev5_SE) %>%
      select(-Ce_Dev5_FDR)
    
    
    # write.csv(jointresults,paste(odir,'/Ecoli_results_All_As_old_screen.csv', sep = ''),row.names = FALSE)
    # write.csv(jointcast,paste(odir,'/Ecoli_results_sidebyside_As_old_screen.csv', sep = ''),row.names = FALSE)
    
    
    # Multiplex 
    jointresults.multi = multiplex(jointresults,c("Plate","Well","Index","Metabolite","MetaboliteU","EcoCycID","KEGG_ID"))
    jointresults.multi2 = multiplex(jointresults,c("Plate","Well","Index","Metabolite","MetaboliteU","EcoCycID","KEGG_ID"),2)
    
    
    NGMa = adjustments[adjustments$Contrast == alln, ]$a
    NGMb = adjustments[adjustments$Contrast == alln, ]$b
    
    
    gradcolours = c('#71B83B','yellow','orange','red')
    
    
    total = jointresults.multi %>%
      filter(x_Contrast == str_C & y_Contrast == str_T & z_Contrast == 'Ce_Dev5')
    
    total = total %>% 
      mutate(score = (abs(y_logFC - (x_logFC-NGMa))*z_logFC),
             score2 = abs(y_logFC - (x_logFC-NGMa))) %>% 
      arrange(score) %>% 
      mutate(Index = factor(Index, levels = Index)) %>% 
      drop_na(z_logFC)
    
    total = total %>% 
      mutate(new_labels = case_when(MetaboliteU %in% labels_4way_nogluc ~ MetaboliteU,
                                    TRUE ~ ''),
             points_alpha = case_when(MetaboliteU %in% labels_4way ~ 1,
                                      TRUE ~ 0.15),
             label_size = case_when(z_logFC == 4 ~ 3,
                                    TRUE ~ 2),
             label_color = case_when(MetaboliteU %in% labels_nuc ~ 'Nucleotide',
                                     MetaboliteU %in%  labels_sug ~ 'Sugars',
                                     TRUE ~ ''),
             new_labels = case_when(new_labels == "Cytidine-3'-monophosphate" ~ "3'-CMP",
                                    new_labels == "Cytidine- 2',3'-cyclic monophosphate" ~ "2',3'-cCMP",
                                    new_labels == "Uridine-3'-monophosphate" ~ "3'-UMP",
                                    new_labels == "Uridine-2',3'-cyclic-monophosphate" ~ "2',3'-cUMP",
                                    new_labels == "Cytidine-2'-monophosphate" ~ "2'-CMP",
                                    new_labels == "Cytidine-5'-monophosphate" ~ "5'-CMP",
                                    new_labels == "Uridine-2'-monophosphate" ~ "2'-UMP",
                                    new_labels == "Uridine-5'-monophosphate" ~ "5'-UMP",
                                    new_labels == "Cytidine-3'-monophosphate" ~ 'cCMP',
                                    new_labels == "Cytidine-3'-monophosphate" ~ 'cCMP',
                                    new_labels == "Cytidine-3'-monophosphate" ~ 'cCMP',
                                    TRUE ~ new_labels)
      ) 


return(total)

}


scatter.4way= function(data.4way, padding = 0.2) {
  data.4way %>% 
    ggplot(aes(x = reorder(MetaboliteU, score2), y = score2)) +
    geom_hline(aes(yintercept = 1), alpha = 0.4, color = 'grey') +
    geom_point(aes(fill = z_logFC, alpha = points_alpha), pch=21, size = 5) +
    # annotate("rect", xmin = 0, xmax = 400, ymin = 1 - sr, ymax = 1 + sr, alpha = .2) +
    scale_fill_gradientn(colours = gradcolours,
                         breaks = c(1,2,3,4), limits = c(1,4), 
                         guide = "legend", name = 'Developmental\nscore') +
    scale_color_manual(values = c('black',
                                  '#C70B00', # nucleotides 
                                  '#310CB3' # sugars
    )) +
    geom_text_repel(aes(x = MetaboliteU, y = score2, 
                        label= new_labels),
                    box.padding = padding,
                    # nudge_y = 0.6,
                    # nudge_x = -10,
                    # segment.curvature = -0.1,
                    # segment.ncp = 3,
                    # segment.angle = 20,
                    size = 3.5,
                    segment.alpha = 0.4,
                    max.overlaps = 100) +
    # geom_text_repel(aes(x = MetaboliteU, y = score2, label = ifelse(z_logFC >= 4, MetaboliteU, '')), 
    #                 box.padding = unit(0.6, "lines"), segment.alpha = 0.4) + 
    labs(x = 'Metabolite',
         y = 'Normalised \nbacterial growth'
    ) +
    coord_cartesian(xlim = c(-10,415),
                    ylim = c(-0.5,4.5)) +
    # ylim(0, 15) +
    scale_size(guide="none") +
    theme_classic() +
    theme(
      # axis.title.y=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x = element_text(size=16, face="bold", color ='black'),
      axis.title.y = element_text(size=16, face="bold", color ='black'),
      plot.title = element_text(size=18),
      legend.text = element_text(size=18)
    ) +
    guides(color = FALSE,
           alpha = FALSE,
           fill = guide_legend(override.aes = list(size=8)))
  
}


scatter.4way.TM= function(data.4way, padding = 0.2, nudge_y = 1.5) {
  data.4way %>% 
    ggplot(aes(x = reorder(MetaboliteU, score2), y = score2)) +
    geom_hline(aes(yintercept = 1), alpha = 0.4, color = 'grey') +
    geom_point(aes(fill = z_logFC, alpha = points_alpha), pch=21, size = 5) +
    # annotate("rect", xmin = 0, xmax = 400, ymin = 1 - sr, ymax = 1 + sr, alpha = .2) +
    scale_fill_gradientn(colours = gradcolours,
                         breaks = c(1,2,3,4), limits = c(1,4), 
                         guide = "legend", name = 'Developmental\nscore') +
    scale_color_manual(values = c('black',
                                  '#C70B00', # nucleotides 
                                  '#310CB3' # sugars
    )) +
    geom_text_repel(aes(x = MetaboliteU, y = score2, 
                        label= new_labels),
                    box.padding = padding,
                    nudge_y = nudge_y,
                    direction     = "y",
                    # nudge_y = 0.6,
                    # nudge_x = -10,
                    # segment.curvature = -0.1,
                    # segment.ncp = 3,
                    # segment.angle = 20,
                    size = 3.5,
                    segment.alpha = 0.4,
                    max.overlaps = 100) +
    # geom_text_repel(aes(x = MetaboliteU, y = score2, label = ifelse(z_logFC >= 4, MetaboliteU, '')), 
    #                 box.padding = unit(0.6, "lines"), segment.alpha = 0.4) + 
    labs(x = 'Metabolite',
         y = 'Normalised \nbacterial growth'
    ) +
    coord_cartesian(xlim = c(-10,415),
                    ylim = c(-0.5,4.5)) +
    # ylim(0, 15) +
    scale_size(guide="none") +
    theme_classic() +
    theme(
      # axis.title.y=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x = element_text(size=16, face="bold", color ='black'),
      axis.title.y = element_text(size=16, face="bold", color ='black'),
      plot.title = element_text(size=18),
      legend.text = element_text(size=18)
    ) +
    guides(color = FALSE,
           alpha = FALSE,
           fill = guide_legend(override.aes = list(size=8)))
  
}


### Plot new versions ####

# BW
BW_total = met_dat(strain = 'BW', experiment = 'T5') 


BW.plot = BW_total %>%
  mutate(glucose_lab = case_when(MetaboliteU == 'alpha-D-Glucose' ~ 'Glucose'),
         .before = "Index") %>% 
  scatter.4way(padding = 0.2) + 
  labs(title = 'BW25113 - Wild type',
       x = NULL, y = NULL) +
  geom_text_repel(aes(label = glucose_lab),
                   box.padding = 0.6,
                  nudge_y = 1, 
                   segment.alpha = 0.4,
                   max.overlaps = 100,
                  size = 5)


# pyrE
pyrE_total = met_dat(strain = 'pyrE', experiment = 'T5') 

pyrE.plot = pyrE_total %>% 
  mutate(glucose_lab = case_when(MetaboliteU == 'alpha-D-Glucose' ~ 'Glucose'),
         .before = "Index") %>% 
  scatter.4way(padding = 0.3) + 
  labs(title = bquote(~Delta*italic(pyrE) ~ "mutant"),
       x = NULL) +
  geom_text_repel(aes(label = glucose_lab),
                  box.padding = 0.6,
                  nudge_y = 0.9, 
                  segment.alpha = 0.4,
                  max.overlaps = 100,
                  size = 5)

# TM
TM_total = met_dat(strain = 'TM', experiment = 'T250') 

TM.plot = TM_total %>% 
  mutate(glucose_lab = case_when(MetaboliteU == 'alpha-D-Glucose' ~ 'Glucose'),
         .before = "Index") %>% 
  scatter.4way.TM(padding = 0.3, nudge_y = 0) +
  labs(title = bquote(~Delta*italic(upp)*Delta*italic(udp)*Delta*italic(udk) ~ "mutant"),
       y = NULL) +
  geom_text_repel(aes(label = glucose_lab),
                  box.padding = 0.6,
                  nudge_y = 1.6, 
                  segment.alpha = 0.4,
                  max.overlaps = 100,
                  size = 5)




p4 = BW.plot / pyrE.plot / TM.plot +
  plot_layout(guides = 'collect')

p4 

# ggsave(p4, file = here('Summary', '4wayScreening_score_2.pdf'),
#        height = 10, width = 13)
# 
# 
# ggsave(p4, file = here('Summary', '4wayScreening_score_3.pdf'),
#        height = 8, width = 10)


ggsave(p4, file = here('FIGURES', '4wayScreening_composite.pdf'),
       height = 10, width = 13)




# check how many metabolites go from 3/4 to 1/2 in pyrE from BW
BW.scores %>% 
  select(Plate, Well, MetaboliteU, BW_worm = z_logFC) %>% 
  left_join(pyrE.scores %>% 
              select(Plate, Well, MetaboliteU, pyrE_worm = z_logFC)
) %>% 
  mutate(big_change = case_when(BW_worm %in% c(1,1.5,2,2.5) & pyrE_worm %in% c(3,3.5,4) ~ 'UP',
                           BW_worm %in%  c(3,3.5,4) & pyrE_worm %in% c(1,1.5,2,2.5) ~ 'DOWN',
                           TRUE ~ 'NO CHANGE'),
         small_change = case_when(BW_worm %in% c(1,1.5) & pyrE_worm %in% c(2,2.5) ~ 'UP',
                                  BW_worm %in%  c(2,2.5) & pyrE_worm %in% c(1,1.5) ~ 'DOWN',
                                  TRUE ~ 'DONT CARE')
         )  %>% view
 





# SUPP scatterplots -------------------------------------------------------

### BW cases ----------------

NGMa = adjustments[adjustments$Contrast == "BW_5FU_Alln", ]$a %>% round(2)
NGMb = adjustments[adjustments$Contrast == "BW_5FU_Alln", ]$b %>% round(2)


#### 0uM ---------

BW_total = met_dat(strain = 'BW', experiment = 'C0') 

BW_total %>%
  mutate(new_labels = case_when(MetaboliteU == 'alpha-D-Glucose' ~ 'Glucose',
                                TRUE ~ new_labels),
         .before = "Index") %>% 
  ggplot(aes(y = y_logFC, x = x_logFC, fill = z_logFC)) +
  geom_abline(intercept = 0, slope = 1, alpha = 1, 
              color = 'grey', 
              linetype = 'longdash') +
  geom_abline(aes(intercept = NGMb, slope = NGMa), 
              alpha = 0.8, color = 'red') +
  geom_errorbarh(aes(y = y_logFC, xmax = x_logFC + x_SE, 
                     xmin = x_logFC - x_SE), height = 0, 
                 alpha = 0.3, color = 'grey50') +
  geom_errorbar(aes(x = x_logFC, ymax = y_logFC + y_SE, 
                    ymin = y_logFC - y_SE), 
                width = 0, alpha = 0.3, 
                color = 'grey50') +
  geom_point(shape = 21,
             size = 4, alpha = 0.7) +
  # geom_text_repel(aes(label= new_labels),
  #                 box.padding = 0.4,
  #                 size = 3.5,
  #                 segment.alpha = 0.4,
  #                 max.overlaps = 100) +
  scale_fill_gradientn(colours = gradcolours,
                         breaks = c(1,2,3,4), 
                         limits = c(1,4), 
                         guide = "legend", 
                         name = 'C. elegans\nphenotype') +
  annotate("text", x = -3.5, y = 2,
           label = glue::glue("y = {NGMa}x + {NGMb}"),
           size=6) +
  theme_light()

  

ggsave(file = here('FIGURES', "supp", '4way_BW_0uM.pdf'),
       height = 10, width = 13)


#### 1.5 uM ---------


BW_total = met_dat(strain = 'BW', experiment = 'T1.5') 

BW_total %>%
  mutate(new_labels = case_when(MetaboliteU == 'alpha-D-Glucose' ~ 'Glucose',
                                TRUE ~ new_labels),
         .before = "Index") %>% 
  ggplot(aes(y = y_logFC, x = x_logFC, fill = z_logFC)) +
  geom_abline(intercept = 0, slope = 1, alpha = 1, 
              color = 'grey', 
              linetype = 'longdash') +
  geom_abline(aes(intercept = NGMb, slope = NGMa), 
              alpha = 0.8, color = 'red') +
  geom_errorbarh(aes(y = y_logFC, xmax = x_logFC + x_SE, 
                     xmin = x_logFC - x_SE), height = 0, 
                 alpha = 0.3, color = 'grey50') +
  geom_errorbar(aes(x = x_logFC, ymax = y_logFC + y_SE, 
                    ymin = y_logFC - y_SE), 
                width = 0, alpha = 0.3, 
                color = 'grey50') +
  geom_point(shape = 21,
             size = 4, alpha = 0.7) +
  # geom_text_repel(aes(label= new_labels),
  #                 box.padding = 0.4,
  #                 size = 3.5,
  #                 segment.alpha = 0.4,
  #                 max.overlaps = 100) +
  scale_fill_gradientn(colours = gradcolours,
                       breaks = c(1,2,3,4), 
                       limits = c(1,4), 
                       guide = "legend", 
                       name = 'C. elegans\nphenotype') +
  annotate("text", x = -3.5, y = 2,
           label = glue::glue("y = {NGMa}x + {NGMb}"),
           size=6) +
  theme_light()



ggsave(file = here('FIGURES', "supp", '4way_BW_1.5uM.pdf'),
       height = 10, width = 13)



#### 5 uM ---------


BW_total = met_dat(strain = 'BW', experiment = 'T5') 

BW_total %>%
  mutate(new_labels = case_when(MetaboliteU == 'alpha-D-Glucose' ~ 'Glucose',
                                TRUE ~ new_labels),
         .before = "Index") %>% 
  ggplot(aes(y = y_logFC, x = x_logFC, fill = z_logFC)) +
  geom_abline(intercept = 0, slope = 1, alpha = 1, 
              color = 'grey', 
              linetype = 'longdash') +
  geom_abline(aes(intercept = NGMb, slope = NGMa), 
              alpha = 0.8, color = 'red') +
  geom_errorbarh(aes(y = y_logFC, xmax = x_logFC + x_SE, 
                     xmin = x_logFC - x_SE), height = 0, 
                 alpha = 0.3, color = 'grey50') +
  geom_errorbar(aes(x = x_logFC, ymax = y_logFC + y_SE, 
                    ymin = y_logFC - y_SE), 
                width = 0, alpha = 0.3, 
                color = 'grey50') +
  geom_point(shape = 21,
             size = 4, alpha = 0.7) +
  # geom_text_repel(aes(label= new_labels),
  #                 box.padding = 0.4,
  #                 size = 3.5,
  #                 segment.alpha = 0.4,
  #                 max.overlaps = 100) +
  scale_fill_gradientn(colours = gradcolours,
                       breaks = c(1,2,3,4), 
                       limits = c(1,4), 
                       guide = "legend", 
                       name = 'C. elegans\nphenotype') +
  annotate("text", x = -3.5, y = 2,
           label = glue::glue("y = {NGMa}x + {NGMb}"),
           size=6) +
  theme_light()



ggsave(file = here('FIGURES', "supp", '4way_BW_5uM.pdf'),
       height = 10, width = 13)


### pyrE cases ---------

NGMa = adjustments[adjustments$Contrast == "pyrE_5FU_Alln", ]$a %>% round(2)
NGMb = adjustments[adjustments$Contrast == "pyrE_5FU_Alln", ]$b %>% round(2)

#### 0uM ---------

BW_total = met_dat(strain = 'pyrE', experiment = 'T0') 
 
BW_total %>%
  mutate(new_labels = case_when(MetaboliteU == 'alpha-D-Glucose' ~ 'Glucose',
                                TRUE ~ new_labels),
         .before = "Index") %>% 
  ggplot(aes(y = y_logFC, x = x_logFC, fill = z_logFC)) +
  geom_abline(intercept = 0, slope = 1, alpha = 1, 
              color = 'grey', 
              linetype = 'longdash') +
  geom_abline(aes(intercept = NGMb, slope = NGMa), 
              alpha = 0.8, color = 'red') +
  geom_errorbarh(aes(y = y_logFC, xmax = x_logFC + x_SE, 
                     xmin = x_logFC - x_SE), height = 0, 
                 alpha = 0.3, color = 'grey50') +
  geom_errorbar(aes(x = x_logFC, ymax = y_logFC + y_SE, 
                    ymin = y_logFC - y_SE), 
                width = 0, alpha = 0.3, 
                color = 'grey50') +
  geom_point(shape = 21,
             size = 4, alpha = 0.7) +
  # geom_text_repel(aes(label= new_labels),
  #                 box.padding = 0.4,
  #                 size = 3.5,
  #                 segment.alpha = 0.4,
  #                 max.overlaps = 100) +
  scale_fill_gradientn(colours = gradcolours,
                       breaks = c(1,2,3,4), 
                       limits = c(1,4), 
                       guide = "legend", 
                       name = 'C. elegans\nphenotype') +
  annotate("text", x = -3.5, y = 2,
           label = glue::glue("y = {NGMa}x + {NGMb}"),
           size=6) +
  theme_light()



ggsave(file = here('FIGURES', "supp", '4way_pyrE_0uM.pdf'),
       height = 10, width = 13)



#### 1.5 uM ---------

BW_total = met_dat(strain = 'pyrE', experiment = 'T1.5') 

BW_total %>%
  mutate(new_labels = case_when(MetaboliteU == 'alpha-D-Glucose' ~ 'Glucose',
                                TRUE ~ new_labels),
         .before = "Index") %>% 
  ggplot(aes(y = y_logFC, x = x_logFC, fill = z_logFC)) +
  geom_abline(intercept = 0, slope = 1, alpha = 1, 
              color = 'grey', 
              linetype = 'longdash') +
  geom_abline(aes(intercept = NGMb, slope = NGMa), 
              alpha = 0.8, color = 'red') +
  geom_errorbarh(aes(y = y_logFC, xmax = x_logFC + x_SE, 
                     xmin = x_logFC - x_SE), height = 0, 
                 alpha = 0.3, color = 'grey50') +
  geom_errorbar(aes(x = x_logFC, ymax = y_logFC + y_SE, 
                    ymin = y_logFC - y_SE), 
                width = 0, alpha = 0.3, 
                color = 'grey50') +
  geom_point(shape = 21,
             size = 4, alpha = 0.7) +
  # geom_text_repel(aes(label= new_labels),
  #                 box.padding = 0.4,
  #                 size = 3.5,
  #                 segment.alpha = 0.4,
  #                 max.overlaps = 100) +
  scale_fill_gradientn(colours = gradcolours,
                       breaks = c(1,2,3,4), 
                       limits = c(1,4), 
                       guide = "legend", 
                       name = 'C. elegans\nphenotype') +
  annotate("text", x = -3.5, y = 2,
           label = glue::glue("y = {NGMa}x + {NGMb}"),
           size=6) +
  theme_light()



ggsave(file = here('FIGURES', "supp", '4way_pyrE_1.5uM.pdf'),
       height = 10, width = 13)



#### 5 uM ---------

BW_total = met_dat(strain = 'pyrE', experiment = 'T5') 

BW_total %>%
  mutate(new_labels = case_when(MetaboliteU == 'alpha-D-Glucose' ~ 'Glucose',
                                TRUE ~ new_labels),
         .before = "Index") %>% 
  ggplot(aes(y = y_logFC, x = x_logFC, fill = z_logFC)) +
  geom_abline(intercept = 0, slope = 1, alpha = 1, 
              color = 'grey', 
              linetype = 'longdash') +
  geom_abline(aes(intercept = NGMb, slope = NGMa), 
              alpha = 0.8, color = 'red') +
  geom_errorbarh(aes(y = y_logFC, xmax = x_logFC + x_SE, 
                     xmin = x_logFC - x_SE), height = 0, 
                 alpha = 0.3, color = 'grey50') +
  geom_errorbar(aes(x = x_logFC, ymax = y_logFC + y_SE, 
                    ymin = y_logFC - y_SE), 
                width = 0, alpha = 0.3, 
                color = 'grey50') +
  geom_point(shape = 21,
             size = 4, alpha = 0.7) +
  # geom_text_repel(aes(label= new_labels),
  #                 box.padding = 0.4,
  #                 size = 3.5,
  #                 segment.alpha = 0.4,
  #                 max.overlaps = 100) +
  scale_fill_gradientn(colours = gradcolours,
                       breaks = c(1,2,3,4), 
                       limits = c(1,4), 
                       guide = "legend", 
                       name = 'C. elegans\nphenotype') +
  annotate("text", x = -3.5, y = 2,
           label = glue::glue("y = {NGMa}x + {NGMb}"),
           size=6) +
  theme_light()



ggsave(file = here('FIGURES', "supp", '4way_pyrE_5uM.pdf'),
       height = 10, width = 13)


### TM cases ---------

NGMa = adjustments[adjustments$Contrast == "TM_5FU_Alln", ]$a %>% round(2)
NGMb = adjustments[adjustments$Contrast == "TM_5FU_Alln", ]$b %>% round(2)

#### 0uM ---------

BW_total = met_dat(strain = 'TM', experiment = 'T0') 

BW_total %>%
  mutate(new_labels = case_when(MetaboliteU == 'alpha-D-Glucose' ~ 'Glucose',
                                TRUE ~ new_labels),
         .before = "Index") %>% 
  ggplot(aes(y = y_logFC, x = x_logFC, fill = z_logFC)) +
  geom_abline(intercept = 0, slope = 1, alpha = 1, 
              color = 'grey', 
              linetype = 'longdash') +
  geom_abline(aes(intercept = NGMb, slope = NGMa), 
              alpha = 0.8, color = 'red') +
  geom_errorbarh(aes(y = y_logFC, xmax = x_logFC + x_SE, 
                     xmin = x_logFC - x_SE), height = 0, 
                 alpha = 0.3, color = 'grey50') +
  geom_errorbar(aes(x = x_logFC, ymax = y_logFC + y_SE, 
                    ymin = y_logFC - y_SE), 
                width = 0, alpha = 0.3, 
                color = 'grey50') +
  geom_point(shape = 21,
             size = 4, alpha = 0.7) +
  # geom_text_repel(aes(label= new_labels),
  #                 box.padding = 0.4,
  #                 size = 3.5,
  #                 segment.alpha = 0.4,
  #                 max.overlaps = 100) +
  scale_fill_gradientn(colours = gradcolours,
                       breaks = c(1,2,3,4), 
                       limits = c(1,4), 
                       guide = "legend", 
                       name = 'C. elegans\nphenotype') +
  annotate("text", x = -3.5, y = 1,
           label = glue::glue("y = {NGMa}x + {NGMb}"),
           size=6) +
  theme_light()



ggsave(file = here('FIGURES', "supp", '4way_TM_0uM.pdf'),
       height = 10, width = 13)



#### 250uM ---------

BW_total = met_dat(strain = 'TM', experiment = 'T250') 

BW_total %>%
  mutate(new_labels = case_when(MetaboliteU == 'alpha-D-Glucose' ~ 'Glucose',
                                TRUE ~ new_labels),
         .before = "Index") %>% 
  ggplot(aes(y = y_logFC, x = x_logFC, fill = z_logFC)) +
  geom_abline(intercept = 0, slope = 1, alpha = 1, 
              color = 'grey', 
              linetype = 'longdash') +
  geom_abline(aes(intercept = NGMb, slope = NGMa), 
              alpha = 0.8, color = 'red') +
  geom_errorbarh(aes(y = y_logFC, xmax = x_logFC + x_SE, 
                     xmin = x_logFC - x_SE), height = 0, 
                 alpha = 0.3, color = 'grey50') +
  geom_errorbar(aes(x = x_logFC, ymax = y_logFC + y_SE, 
                    ymin = y_logFC - y_SE), 
                width = 0, alpha = 0.3, 
                color = 'grey50') +
  geom_point(shape = 21,
             size = 4, alpha = 0.7) +
  # geom_text_repel(aes(label= new_labels),
  #                 box.padding = 0.4,
  #                 size = 3.5,
  #                 segment.alpha = 0.4,
  #                 max.overlaps = 100) +
  scale_fill_gradientn(colours = gradcolours,
                       breaks = c(1,2,3,4), 
                       limits = c(1,4), 
                       guide = "legend", 
                       name = 'C. elegans\nphenotype') +
  annotate("text", x = -3.5, y = 1,
           label = glue::glue("y = {NGMa}x + {NGMb}"),
           size=6) +
  theme_light()



ggsave(file = here('FIGURES', "supp", '4way_TM_250uM.pdf'),
       height = 10, width = 13)








# save TABLES for PAPER ---------------------------------------------------


## BW ------


bw_c = met_dat(strain = 'BW', experiment = 'C0') %>% 
  select(-Index, 
         -x_Description:-x_Contrast_type, 
         -y_Description:-y_Contrast_type,
         -z_Description:-z_Contrast_type,
         -z_FDR:-label_color, -z_SE) %>% 
  rename(
    log2FC_control = x_logFC, log2FC_control_SE = x_SE,
    FDR_control = x_FDR, log2FC_5FU = y_logFC,
    log2FC_5FU_SE = y_SE, FDR_5FU = y_FDR,
    worm_score_control = z_logFC
  ) %>% 
  mutate(description = "BW25113", .before = "Plate")

bw_15 = met_dat(strain = 'BW', experiment = 'T1.5') %>% 
  select(-Index, 
         -x_Description:-x_Contrast_type, 
         -y_Description:-y_Contrast_type,
         -z_Description:-z_Contrast_type,
         -z_FDR:-label_color, -z_SE) %>% 
  rename(
    log2FC_control = x_logFC, log2FC_control_SE = x_SE,
    FDR_control = x_FDR, log2FC_5FU = y_logFC,
    log2FC_5FU_SE = y_SE, FDR_5FU = y_FDR,
    worm_score_1.5 = z_logFC
  ) %>% 
  mutate(description = "BW25113", .before = "Plate")

bw_5 = met_dat(strain = 'BW', experiment = 'T5') %>% 
  select(-Index, 
         -x_Description:-x_Contrast_type, 
         -y_Description:-y_Contrast_type,
         -z_Description:-z_Contrast_type,
         -z_FDR:-label_color, -z_SE) %>% 
  rename(
    log2FC_control = x_logFC, log2FC_control_SE = x_SE,
    FDR_control = x_FDR, log2FC_5FU = y_logFC,
    log2FC_5FU_SE = y_SE, FDR_5FU = y_FDR,
    worm_score_5 = z_logFC
  ) %>% 
  mutate(description = "BW25113", .before = "Plate")


bw_c %>% left_join(bw_15) %>% left_join(bw_5) %>% 
  write_csv("FIGURES/DATATABLES/4way_BW.csv")


## pyrE ------


bw_c = met_dat(strain = 'pyrE', experiment = 'T0') %>% 
  select(-Index, 
         -x_Description:-x_Contrast_type, 
         -y_Description:-y_Contrast_type,
         -z_Description:-z_Contrast_type,
         -z_FDR:-label_color, -z_SE) %>% 
  rename(
    log2FC_control = x_logFC, log2FC_control_SE = x_SE,
    FDR_control = x_FDR, log2FC_5FU = y_logFC,
    log2FC_5FU_SE = y_SE, FDR_5FU = y_FDR,
    worm_score_control = z_logFC
  ) %>% 
  mutate(description = "pyrE", .before = "Plate")

bw_15 = met_dat(strain = 'pyrE', experiment = 'T1.5') %>% 
  select(-Index, 
         -x_Description:-x_Contrast_type, 
         -y_Description:-y_Contrast_type,
         -z_Description:-z_Contrast_type,
         -z_FDR:-label_color, -z_SE) %>% 
  rename(
    log2FC_control = x_logFC, log2FC_control_SE = x_SE,
    FDR_control = x_FDR, log2FC_5FU = y_logFC,
    log2FC_5FU_SE = y_SE, FDR_5FU = y_FDR,
    worm_score_1.5 = z_logFC
  ) %>% 
  mutate(description = "pyrE", .before = "Plate")

bw_5 = met_dat(strain = 'pyrE', experiment = 'T5') %>% 
  select(-Index, 
         -x_Description:-x_Contrast_type, 
         -y_Description:-y_Contrast_type,
         -z_Description:-z_Contrast_type,
         -z_FDR:-label_color, -z_SE) %>% 
  rename(
    log2FC_control = x_logFC, log2FC_control_SE = x_SE,
    FDR_control = x_FDR, log2FC_5FU = y_logFC,
    log2FC_5FU_SE = y_SE, FDR_5FU = y_FDR,
    worm_score_5 = z_logFC
  ) %>% 
  mutate(description = "pyrE", .before = "Plate")


bw_c %>% left_join(bw_15) %>% left_join(bw_5) %>% 
  write_csv("FIGURES/DATATABLES/4way_pyrE.csv")

## TM ------


bw_c = met_dat(strain = 'TM', experiment = 'T0') %>% 
  select(-Index, 
         -x_Description:-x_Contrast_type, 
         -y_Description:-y_Contrast_type,
         -z_Description:-z_Contrast_type,
         -z_FDR:-label_color, -z_SE) %>% 
  rename(
    log2FC_control = x_logFC, log2FC_control_SE = x_SE,
    FDR_control = x_FDR, log2FC_5FU = y_logFC,
    log2FC_5FU_SE = y_SE, FDR_5FU = y_FDR,
    worm_score_control = z_logFC
  ) %>% 
  mutate(description = "TM", .before = "Plate")

bw_15 = met_dat(strain = 'TM', experiment = 'T250') %>% 
  select(-Index, 
         -x_Description:-x_Contrast_type, 
         -y_Description:-y_Contrast_type,
         -z_Description:-z_Contrast_type,
         -z_FDR:-label_color, -z_SE) %>% 
  rename(
    log2FC_control = x_logFC, log2FC_control_SE = x_SE,
    FDR_control = x_FDR, log2FC_5FU = y_logFC,
    log2FC_5FU_SE = y_SE, FDR_5FU = y_FDR,
    worm_score_250 = z_logFC
  ) %>% 
  mutate(description = "TM", .before = "Plate")


bw_c %>% 
  ggplot(aes(x = log2FC_control, y = log2FC_5FU)) +
  geom_point()


bw_c %>% left_join(bw_15)  %>% 
  write_csv("FIGURES/DATATABLES/4way_TM.csv")

# ternary plot ------------------------------------------------------------



remove = c('Tween 20', 'Tween 40', 'Tween 80')

pyrE.scores
BW.scores
TM.scores

scores = BW.scores %>% 
  bind_rows(pyrE.scores, TM.scores) %>%
  filter(!(MetaboliteU %in% remove)) %>% 
  arrange(x_Contrast,Plate, MetaboliteU)

# BW_C, pyrE_C, TM_C

BW = scores %>% filter(x_Contrast == 'BW_C') %>% select(score) %>% t %>% as.vector
pyrE = scores %>% filter(x_Contrast == 'pyrE_C') %>% select(score) %>% t %>% as.vector 
TM = scores %>% filter(x_Contrast == 'TM_C') %>% select(score) %>% t %>% as.vector

wBW = scores %>% filter(x_Contrast == 'BW_C') %>% select(worm = z_logFC) %>% t %>% as.vector
wpyrE = scores %>% filter(x_Contrast == 'pyrE_C') %>% select(worm = z_logFC) %>% t %>% as.vector
wTM = scores %>% filter(x_Contrast == 'TM_C') %>% select(worm = z_logFC) %>% t %>% as.vector



# # normalise between 0 and 1
# BW = scores %>% filter(x_Contrast == 'BW_C') %>% mutate(score3 = ((score - min(score))/(max(score)- min(score)))) %>% select(score3) %>% t %>% as.vector
# pyrE = scores %>% filter(x_Contrast == 'pyrE_C') %>% mutate(score3 = ((score - min(score))/(max(score)- min(score)))) %>% select(score3) %>% t %>% as.vector 
# TM = scores %>% filter(x_Contrast == 'TM_C') %>% mutate(score3 = ((score - min(score))/(max(score)- min(score)))) %>% select(score3) %>% t %>% as.vector


label = scores %>% filter(x_Contrast == 'BW_C') %>% select(MetaboliteU) %>% t %>% as.character

df = data.frame(BW, pyrE, TM, wBW, wpyrE, wTM, label) 

thr = 4
df = df %>% mutate(Size = case_when(wBW >= thr | wpyrE >= thr | wTM >= thr ~ 3,
                               wBW < thr | wpyrE < thr | wTM < thr ~ 1)) %>% 
  mutate(Phenotype = case_when(wTM == 4 ~ 'TM',
                               wBW == 4 & wpyrE == 4 ~ 'BW and pyrE',
                               wpyrE == 4 & wBW != 4 ~ 'pyrE',
                               wpyrE != 4 & wBW == 4 ~ 'BW',
                               TRUE ~ 'Other'),
         Phenotype = as.factor(Phenotype)) %>% 
  as_tibble

df_small = df %>% 
  filter(Size == 1)

df_big = df %>% 
  filter(Size == 5)

# axis layout
axis <- function(title) {
  list(
    title = title,
    titlefont = list(
      size = 10
    ),
    tickfont = list(
      size = 15
    ),
    tickcolor = 'rgba(0,0,0,0)',
    ticklen = 5
  )
}


fig = df_small %>% plot_ly()
# add small points
fig = fig %>% add_trace(
  type = 'scatterternary',
  mode = 'markers',
  a = ~BW,
  b = ~pyrE,
  c = ~TM,
  text = ~label,
  name = 'Other',
  marker = list( 
    size = 4,
    color = '#0CEB96',
    sizemode = 'diameter',
    opacity = 0.5
  )
)

fig = fig %>% add_trace(
  data = df_big,
  type = 'scatterternary',
  # mode = 'markers',
  mode = 'text',
  a = ~BW,
  b = ~pyrE,
  c = ~TM,
  text = ~label,
  textposition = "top right",
  color = ~Phenotype,
  marker = list( 
    size = 10,
    # color = '#DB392E',
    sizemode = 'diameter'
  )
)

fig <- fig %>% layout(
  ternary = list(
    sum = 100,
    aaxis = list(title = 'BW', min=0.1),
    baxis = list(title = 'pyrE', min=0.1),
    caxis = list(title = 'TM', min = 0)
  )
)

fig

# to save files from plot_ly
server = orca_serve()
server$export(fig,here('Summary',"4wayScreening_score_ternary.pdf"))


df %>% View


# version without irrelevant points

df.reduced = df %>% filter(size == 14) %>% 
  mutate(Phenotype = case_when(wTM == 4 ~ 'TM',
                               wBW == 4 & wpyrE == 4 ~ 'BW and pyrE',
                               wpyrE == 4 & wBW != 4 ~ 'pyrE',
                               wpyrE != 4 & wBW == 4 ~ 'BW'))


fig <- df.reduced %>% plot_ly()
fig <- fig %>% add_trace(
  type = 'scatterternary',
  mode = 'markers',
  a = ~BW,
  b = ~pyrE,
  c = ~TM,
  text = ~label,
  color = ~Phenotype,
  marker = list( 
    # symbol = 1,
    
    size = ~size,
    line = list('width' = 2)
  )
)
fig <- fig %>% layout(
  # title = "Simple Ternary Plot with Markers",
  ternary = list(
      sum = 100,
      aaxis = list(title = 'BW', min=0.3),
      baxis = list(title = 'pyrE', min=0.1),
      caxis = list(title = 'TM', min = 0)
    )
)

fig





### ggtern ####

# versio without nucleotide classes

library(ggtern)

df %>% 
  mutate(Phenotype = factor(Phenotype, levels = c('BW','BW and pyrE',
                                                  'pyrE','TM','Other'))) %>% 
  mutate(label = case_when(Phenotype == 'Other' ~ '',
                           TRUE ~ label)) %>% 
  ggtern(aes(pyrE,BW,TM, color = Phenotype, size = Size,
             alpha = Size)) +
  # stat_density_tern(geom = 'polygon',
  #                   bdl = 0.01,
  #                   n = 100,
  #                   aes(fill=..level..,
  #                       alpha=..level../5),
  #                   weight = 1,
  #                   base = "ilr") +
  geom_point() +
  theme_rgbw() +
  labs(title = "4-way screen")    +
  scale_color_manual(values = c('#DB1D51', # red
                                '#F28413', # orange
                                '#D6D613', # yellow
                                '#5913F2', # dark blue
                                '#0CEB96'  # light blue
                                )) + 
  scale_fill_gradient(low = "blue",high = "red") +
  scale_size(range = c(1, 3)) +
  scale_alpha(range = c(0.4,1)) + 
  guides(fill = "none", alpha = "none", size = 'none') +
  annotate(geom  = 'text',
                    # nuc  # sug  #pyrE  #BW
           x     = c(0.5  , 0.2 , 0.2 , 0.2),
           y     = c(0.5  , 0.2 , 0.45 , 0.8),
           z     = c(0.1 , 1 , 0.3 , 0.3),
           # angle = c(0,30,60),
           # vjust = c(1.5,0.5,-0.5),
           label = c("Pyrimidines","Sugars","Thymidine","Uridine"),
           color = c("#F28413","#5913F2",'#B0B000','#DB1D51')) +
  theme_nomask() +
  labs(color = 'Rescued\nphenotype') + 
  theme(legend.key = element_rect(fill = NA, color = NA))


ggsave(file = here('Summary', '4wayScreening_ggtern.pdf'),
       height = 9, width = 12)



### ggtern with classes ####

df_class = EC_classes %>% rename(label = MetaboliteU) %>%
  select(label, EcoCyc_Classes) %>%
  full_join(df)

# classes from enrichment
tm_cats = c('|Alcohols|','|Carbohydrates|','|Glycans|',
            '|Hexitols|','|Sugar−alcohols|','|Sugar|')

bw_cats = c('|All−Nucleosides|','|Nucleosides|','|Nucleotides|',
            '|Organic−heterocyclic−compound|','|Organonitrogen−Heterocyclic−Compounds|',
            '|Pyrimidine−Nucleosides|','|Pyrimidine−ribonucleotides|',
            '|Pyrimidines|')

total_cats = c(tm_cats,bw_cats)

# clean the dataset and create two categories for sugars and nucleotides, 
# drop na values, and tidy everything a bit more
df_class = df_class %>%
  mutate(
    Phenotype = factor(Phenotype, levels = c('BW','BW and pyrE',
                                             'pyrE','TM','Other')),
    EcoCyc_Classes = case_when(EcoCyc_Classes %in% tm_cats & Phenotype != 'Other' ~ 'Sugars',
                                    EcoCyc_Classes %in% bw_cats & Phenotype != 'Other' ~ 'Nucleotides',
                                    TRUE ~ '')) %>%
  distinct(label,.keep_all = T) %>%
  mutate(label = case_when(Phenotype == 'Other' ~ '',
                           TRUE ~ label)) %>%
  drop_na(BW)

# fix table
df_class = df_class %>%
  mutate(EcoCyc_Classes = case_when(
    Phenotype %in% c('BW and pyrE','BW','pyrE') ~ 'Nucleotides',
    Phenotype == 'TM' ~ 'Sugars',
    TRUE ~ ''))

# Plot with density according to metabolites

df_class %>%
    filter(EcoCyc_Classes %in% c('Sugars','Nucleotides')) %>%
    ggtern(aes(pyrE,BW,TM,
               color = EcoCyc_Classes)) +
    stat_density_tern(geom = 'polygon',
                      size = 0.5,
                      bdl = 0.05,
                      n = 100,
                      aes(alpha=..level..,
                          color = EcoCyc_Classes),
                      weight = 1,
                      base = "ilr") + 
  geom_point(data = df_class, pch = 21, 
             color = 'grey70',
             aes(group = Phenotype, 
                 fill = Phenotype, 
                 size=Size)) +
  scale_fill_manual(
    values = c('#DB1D51', # red
               '#F28413', # orange
               '#D6D613', # yellow
               '#5913F2', # dark blue
               '#0CEB96'  # light blue
    )
  ) +
  scale_color_manual(
    values = c(
      '#FA754B', # orange
      '#3D9EF0'# blue
      )
  ) + 
  theme_rgbw() +
  scale_size(range = c(1, 5)) +
  guides(alpha = "none", size = 'none') +
  theme_nomask() +
  labs(color = 'Rescued\nphenotype') + 
  theme(legend.key = element_rect(fill = NA, color = NA))

ggsave(file = here('Summary', '4wayScreening_ggtern_class_density_v1.pdf'),
       height = 9, width = 12)




### ggtern with class filling ####

# Plot with density according to metabolites

df_class %>%
  filter(EcoCyc_Classes %in% c('Sugars','Nucleotides')) %>%
  ggtern(aes(pyrE,BW,TM,
             color = EcoCyc_Classes)) +
  stat_density_tern(geom = 'polygon',
                    size = 0.5,
                    bdl = 0.058,
                    n = 100,
                    aes(alpha=..level..,
                        fill = EcoCyc_Classes),
                    weight = 1,
                    base = "ilr") + 
  geom_point(data = df_class, pch = 21, 
             color = 'grey70',
             aes(group = Phenotype, 
                 fill = Phenotype, 
                 size=Size)) +
  scale_fill_manual(
    values = c('#DB1D51', # red
               '#FA4F22', # orange
               '#FA754B',  # orange CLASS
               '#0CEB96', # light blue
               '#D6D613', # yellow
               '#3D9EF0', # blue CLASS
               '#5913F2' # dark blue
               
    )
  ) +
  scale_color_manual(
    values = c(
      '#FA754B', # orange
      '#3D9EF0'  # blue
    )
  ) + 
  theme_rgbw() +
  scale_size(range = c(1, 5)) +
  guides(alpha = "none", size = 'none') +
  theme_nomask() +
  labs(color = 'Rescued\nphenotype') + 
  theme(legend.key = element_rect(fill = NA, color = NA))


ggsave(file = here('Summary', '4wayScreening_ggtern_class_density_v4.pdf'),
       height = 9, width = 12)




df %>% 
  select(-wBW:-wpyrE, -Size) %>%
  rename(Metabolite = label) %>% 
  write_csv("FIGURES/DATATABLES/bact_scores_ternary.csv")
  



# more simpler 4-way plots ------------------------------------------------


strain = 'BW'
experiment = "T5"

bw_worms = worm.sum %>% filter(Strain == strain, Experiment == experiment) %>% 
  select(-Experiment, -Mean:-SD)

bw_scores =
  results %>% 
  filter(Contrast %in% c('BW_5FU', 'BW_C', 'BW_T')) %>% 
  select(Contrast, Strain:FDR, -SE, -p.value, -FDR, -t.value) %>% 
  pivot_wider(names_from = Contrast, values_from = logFC) %>% 
  left_join(bw_worms) %>% 
  mutate(slope = adjustments %>% 
           filter(Contrast == 'BW_5FU_Alln') %>% 
           pull(a)) %>%
  mutate(BW_adj = (BW_C - slope),
         score = abs(BW_T - BW_adj))





bw_scores %>%
  ggplot(aes(x = fct_reorder(MetaboliteU, score), y = score)) +
  geom_point(aes(fill = Median), shape = 21) +
  scale_fill_gradientn(colours = gradcolours,
                       breaks = c(1,2,3,4), limits = c(1,4),
                       guide = "legend", name = 'C. elegans\nphenotype') +
  theme_cowplot(14)







# PCA plot Paper ----------------------------------------------------------


pca_b_data = data.b %>%
  filter(Strain %in% c("BW", "pyrE", "TM") & 
           Metabolite != 'Negative Control')

# some data frame transformations
pca_b_data = pca_b_data %>%
  select(SampleID, MetaboliteU, logAUC) %>%
  spread(MetaboliteU, logAUC) %>%
  data.frame(check.names = F)

rownames(pca_b_data) = pca_b_data[,1]; pca_b_data = pca_b_data[,-1]

# lets compute the PCA
res.pca = PCA(pca_b_data, scale.unit = TRUE, ncp = 5, graph = F)

# metadata 
meta_var = bioinfo %>% filter(Strain %in%  c("BW", "pyrE", "TM"))

# extract info about the individuals
ind = get_pca_ind(res.pca)
ind_df = data.frame(ind$coord[,1], ind$coord[,2], ind$coord[,3], meta_var$Type, 
                    meta_var$Strain)

colnames(ind_df) = c('Dim1', 'Dim2', 'Dim3', 'Type', 'Strain')

# get ellipses based on the correlation
getellipse = function(x, y, sc = 1) {
  as.data.frame(ellipse::ellipse(cor(x, y),
                                 scale = c(sd(x) * sc, sd(y) * sc),
                                 centre = c(mean(x), mean(y))))
}

# make a data frame from ellipses
ell = ind_df %>% group_by(Type, Strain) %>% do(getellipse(.$Dim1, .$Dim2, 1)) %>% data.frame

# plot!
ggplot(ind_df, aes(x = Dim1, y = Dim2, color = Strain, group = interaction(Type, Strain))) + 
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell, aes(x = x, y = y, group = interaction(Type, Strain), linetype = Type), size = 1) +
  geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Type, Strain), 
                               linetype = Type, fill = Strain), size = 1, alpha = 0.3) +
  xlab(paste("PC1 - ", round(res.pca$eig[1,2], 1), " % of variance", sep = "")) + 
  ylab(paste("PC2 - ", round(res.pca$eig[2,2], 1), " % of variance", sep = "")) +
  cowplot::theme_cowplot(14)

# save the file
ggsave(file = here('FIGURES', "PCA_main_fig2.pdf"),
       width = 100, height = 80, units = 'mm', 
       scale = 2, device = cairo_pdf, family = "Arial")










