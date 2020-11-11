
## Master script for Leo's plots, thesis edition

# load libraries ####
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



options(width = 170)


# sub library -------------------------------------------------------------

# glucose

# read pathways
library(readxl)
library(sda)

pathways = read_excel("Sub_library/big_screen/EcoCyc_pathways_KeioSublibrary_08-12-18.xlsx") %>% 
  select(Genes, Pathway = KEGG3)

# read and transform data
glucose = read_xlsx('Sub_library/Summary_Keio Sublibrary_Glucose Supp_10mM_N2 worms.xlsx', sheet = 'Summary_good') 

glucose = glucose %>%
  gather(Drug, Score, `0`, `1`, `2.5`, `5`) %>% # make table long
  unite(ID, Genes, Supplement, Replicate, remove = FALSE) %>% 
  mutate(Replicate = as.numeric(Replicate),
         Supplement_mM = as.factor(Supplement_mM),
         Supplement = as.factor(Supplement),
         Genes = as.factor(Genes),
         Pathway = as.factor(Pathway),
         ID = as.factor(ID),
         Drug = as.factor(Drug)) %>%
  select(Genes, Well, ID, Pathway, Replicate, Supplement, Supplement_mM, Drug, Score)


# table with genotypes and pathways

# pathways = glucose %>% select(Genes, Pathway) %>% unique

# data summary
glu.sum = glucose %>%
  group_by(Supplement, Supplement_mM, Drug, Genes) %>%
  summarise(Median_Score = median(Score, na.rm = TRUE),
            MAD = mad(Score, na.rm = TRUE),
            Mean = mean(Score, na.rm = TRUE),
            SD = sd(Score, na.rm = TRUE)) %>%
  mutate(BW_Score = Median_Score[Genes == 'BW'],
         BW_Mean = Mean[Genes == 'BW']) %>%
  ungroup %>%
  group_by(Genes, Drug) %>%
  ungroup %>%
  mutate(BW_norm = Median_Score - BW_Score,
         BW_Mean_norm = Mean - BW_Mean)


### SCATTER PLOT
# original plot
drug = 5
pos = position_jitter(width = 0.05, height = 0.05, seed = 1) # to plot names in jitter positions
glu.sum %>% 
  filter(Drug == drug) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>%
  unite(Supp, Supplement, Supplement_mM) %>%
  spread(Supp, BW_norm) %>%
  left_join(pathways) %>% 
  ggplot(aes(x = Glucose_0, y = Glucose_10)) + 
  geom_hline(yintercept = 0, colour = 'grey30') +
  geom_vline(xintercept = 0, colour = 'grey30') +
  geom_point(aes(colour = Pathway),position = pos, size = 2) + 
  # scale_color_manual(values = grad) + 
  # scale_fill_manual(values = grad) +
  # geom_text_repel(aes(label = ifelse(Genes == 'ppnP', as.character(Genes), '')), position = pos) + # this point went rogue
  geom_text_repel(aes(label = Genes), position = pos) +
  labs(title = expression(paste("5FU + Glucose effect on ", italic('C. elegans'), " N2 phenotype", sep = '')),
       x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' N2 phenotype', sep = ' ')),
       y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' N2 phenotype ' , bold('(Glucose)'), sep = ' '))) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.text = element_text(size = 6)) + 
  guides(colour = guide_legend(override.aes = list(size = 4))) # make lengend points larger
# 
# quartz.save(file = here('Summary', paste0('Scatter_sub_lib_glucose_',as.character(drug),'uM.pdf')),
#             type = 'pdf', dpi = 300, height = 10, width = 12)



# centroids


centr = glu.sum %>% 
  filter(Drug == drug) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>%
  unite(Supp, Supplement, Supplement_mM) %>%
  spread(Supp, BW_norm) %>%
  left_join(pathways) %>% 
  group_by(Pathway) %>% 
  summarise(Glucose_0_centroid = mean(Glucose_0),
            sd_0 = sd(Glucose_0),
            Glucose_10_centroid = mean(Glucose_10),
            sd_10 = sd(Glucose_10)) 

glu.sum %>% 
  filter(Drug == drug) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>%
  unite(Supp, Supplement, Supplement_mM) %>%
  spread(Supp, BW_norm) %>%
  left_join(pathways) %>% 
  left_join(centr) %>% 
  ggplot(aes(x = Glucose_0, y = Glucose_10)) + 
  geom_hline(yintercept = 0, colour = 'grey30') +
  geom_vline(xintercept = 0, colour = 'grey30') +
  geom_point(aes(x = Glucose_0_centroid, y = Glucose_10_centroid, color = Pathway), size = 12, alpha = 0.1) +
  geom_point(aes(colour = Pathway),position = pos, size = 2) + 
  geom_text_repel(aes(label = Genes), position = pos) +
  labs(title = expression(paste("5FU + Glucose effect on ", italic('C. elegans'), " N2 phenotype", sep = '')),
       x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' N2 phenotype', sep = ' ')),
       y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' N2 phenotype ' , bold('(Glucose)'), sep = ' '))) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.text = element_text(size = 6)) + 
  guides(colour = guide_legend(override.aes = list(size = 4))) # make lengend points larger

ggsave(file = here('Summary', 'Dev.assay.centroids_genes.pdf'),
       height = 6, width = 7)

# scatter plot for centroids

centr %>% 
  ggplot(aes(x = Glucose_0_centroid, y = Glucose_10_centroid)) + 
  geom_hline(yintercept = 0, colour = 'grey30') +
  geom_vline(xintercept = 0, colour = 'grey30') +
  geom_point(aes(colour = Pathway), size = 7) + 
  # geom_errorbar(aes(xmin = Glucose_0_centroid - sd_0, xmax = Glucose_0_centroid + sd_0)) +
  # geom_errorbar(aes(ymin = Glucose_10_centroid - sd_10, ymax = Glucose_10_centroid + sd_10)) +
  # geom_text_repel(aes(label = Genes), position = pos) +
  labs(title = expression(paste("5FU + Glucose effect on ", italic('C. elegans'), " N2 phenotype", sep = '')),
       x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' N2 phenotype', sep = ' ')),
       y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' N2 phenotype ' , bold('(Glucose)'), sep = ' '))) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.text = element_text(size = 6)) + 
  guides(colour = guide_legend(override.aes = list(size = 4))) # make lengend points larger
ggsave(file = here('Summary', 'Dev.assay.centroids.pdf'),
       height = 6, width = 7)
# histogram of centroids and sd

glu.sum %>% 
  filter(Drug == drug) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>%
  unite(Supp, Supplement, Supplement_mM) %>%
  # spread(Supp, BW_norm) %>%
  left_join(pathways) %>% 
  group_by(Pathway, Supp) %>% 
  summarise(Mean = mean(BW_norm),
            SD = sd(BW_norm)) %>% 
  filter(Pathway != 'Control') %>% 
  ggplot(aes(x = reorder(Pathway, Mean), y = Mean, fill = Supp, group = Supp)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = (Mean) - (SD), ymax = (Mean) + (SD)), position = position_dodge(), alpha = 1, color = 'grey50') +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



## reorganization of resistance/sensitivity


glu.sum %>% 
  filter(Drug == drug) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>%
  unite(Supp, Supplement, Supplement_mM) %>%
  # spread(Supp, BW_norm) %>%
  left_join(pathways) %>% 
  mutate(resistance = case_when(BW_norm >= 1 ~ 'Resistant',
                                BW_norm <= -1 ~ 'Sensitive',
                                BW_norm < 1 & BW_norm > -1 ~ 'Neutral')) %>% 
  count(Pathway, Supp, resistance) %>%
  ggplot(aes(x = Pathway, y = n, group = interaction(Supp, resistance), fill = interaction(resistance, Supp))) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






thr = 1
# with total n genes
glu.sum %>% 
  filter(Drug == drug) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>%
  unite(Supp, Supplement, Supplement_mM) %>%
  # spread(Supp, BW_norm) %>%
  left_join(pathways) %>% 
  mutate(resistance = case_when(BW_norm >=  thr ~ 'Resistant',
                                BW_norm <= -thr ~ 'Sensitive',
                                BW_norm < thr & BW_norm > -thr ~ 'Neutral')) %>% 
  count(Pathway, Supp, resistance) %>%
  group_by(Pathway, Supp) %>% 
  mutate(Total = sum(n),
         Prop = n/Total * 100) %>% 
  filter(Pathway != 'Control') %>% 
  ggplot(aes(x = Pathway, y = n, group = resistance, fill = resistance)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  facet_wrap(~Supp, ncol = 1) +
  theme_light() +
  labs(x = 'Pathways',
       y = 'Number of genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = here('Summary', 'Dev.assay.glucose_barplot_total.pdf'),
       height = 6, width = 7)

# with percentages
glu.sum %>% 
  filter(Drug == drug) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>%
  unite(Supp, Supplement, Supplement_mM) %>%
  # spread(Supp, BW_norm) %>%
  left_join(pathways) %>% 
  mutate(resistance = case_when(BW_norm >= thr ~ 'Resistant',
                                BW_norm <= -thr ~ 'Sensitive',
                                BW_norm < thr & BW_norm > -thr ~ 'Neutral')) %>% 
  count(Pathway, Supp, resistance) %>%
  group_by(Pathway, Supp) %>% 
  mutate(Total = sum(n),
         Prop = n/Total * 100) %>% 
  filter(Pathway != 'Control') %>% 
  ggplot(aes(x = Pathway, y = Prop, group = resistance, fill = resistance)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  facet_wrap(~Supp, ncol = 1) +
  theme_light() +
  labs(x = 'Pathways',
       y = 'Proportion of genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(file = here('Summary', 'Dev.assay.glucose_barplot_percentage.pdf'),
       height = 6, width = 7)










