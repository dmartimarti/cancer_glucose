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

options(width = 170)

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

pyrE.scores = total %>% 
  select(Plate, Well, MetaboliteU, x_Contrast, z_logFC, score, score2)


pyrE.plot = total %>% 
  ggplot(aes(x = reorder(MetaboliteU, score), y = score)) +
  geom_point(aes(colour = z_logFC, size = z_logFC)) +
  geom_hline(aes(yintercept = 1), alpha = 0.9, color = 'grey') +
  scale_colour_gradientn(colours = gradcolours,
                         breaks = c(1,2,3,4), limits = c(1,4), 
                         guide = "legend", name = 'C. elegans\nphenotype') +
  geom_text_repel(aes(x = MetaboliteU, y = score, label = ifelse(z_logFC >= 4, MetaboliteU, '')), 
                  box.padding = unit(0.6, "lines"), segment.alpha = 0.4) + 
  labs(x = 'Metabolite',
       y = 'Composed score') +
  coord_cartesian(xlim = c(0,400)) +
  scale_size(guide=FALSE) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 



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



gradcolours = c('#71B83B','yellow','orange','red')



total = jointresults.multi %>%
  filter(x_Contrast == str_C & y_Contrast == str_T & z_Contrast == 'Ce_Dev5')

total = total %>% 
  mutate(score = (abs(y_logFC - (x_logFC-NGMa))*z_logFC),
         score2 = abs(y_logFC - (x_logFC-NGMa))) %>% 
  arrange(score) %>% 
  mutate(Index = factor(Index, levels = Index)) %>% 
  drop_na(z_logFC)

BW.scores = total %>% 
  select(Plate, Well, MetaboliteU, x_Contrast, z_logFC, score, score2)



BW.plot = total %>% 
  ggplot(aes(x = reorder(MetaboliteU, score), y = score)) +
  geom_point(aes(colour = z_logFC, size = z_logFC)) +
  geom_hline(aes(yintercept = 1), alpha = 0.9, color = 'grey') +
  scale_colour_gradientn(colours = gradcolours,
                         breaks = c(1,2,3,4), limits = c(1,4), 
                         guide = "legend", name = 'C. elegans\nphenotype') +
  geom_text_repel(aes(x = MetaboliteU, y = score, label = ifelse(z_logFC >= 4, MetaboliteU, '')), 
                  box.padding = unit(0.6, "lines"), segment.alpha = 0.4) + 
  labs(x = 'Metabolite',
       y = 'Composed score') +
  coord_cartesian(xlim = c(0,400)) +
  scale_size(guide=FALSE) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 





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

TM.scores = total %>% 
  select(Plate, Well, MetaboliteU, x_Contrast, z_logFC, score, score2)



TM.plot = total %>% 
  ggplot(aes(x = reorder(MetaboliteU, score), y = score)) +
  geom_point(aes(colour = z_logFC, size = z_logFC)) +
  # annotate("rect", xmin = 0, xmax = 400, ymin = 1 - sr, ymax = 1 + sr, alpha = .2) +
  geom_hline(aes(yintercept = 1), alpha = 0.9, color = 'grey') +
  scale_colour_gradientn(colours = gradcolours,
                         breaks = c(1,2,3,4), limits = c(1,4), 
                         guide = "legend", name = 'C. elegans\nphenotype') +
  geom_text_repel(aes(x = MetaboliteU, y = score, label = ifelse(z_logFC >= 4, MetaboliteU, '')), 
                  box.padding = unit(0.6, "lines"), segment.alpha = 0.4) + 
  labs(x = 'Metabolite',
       y = 'Composed score') +
  coord_cartesian(xlim = c(0,400)) +
  scale_size(guide=FALSE) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 




library(ggpubr)
ggarrange(BW.plot, pyrE.plot, TM.plot, nrow = 3)
p4 = ggarrange(BW.plot, pyrE.plot, TM.plot, nrow = 3)


ggsave(p4, file = here('Summary', '4wayScreening_score.pdf'),
       height = 10, width = 13)








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



