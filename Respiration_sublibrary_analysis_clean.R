# load libraries
library(tidyverse)
library(readxl)
library(ggrepel)
library(gtools)
library(broom)
library(here)
library(openxlsx)
library(ggpubr)
library(cowplot)

theme_set(theme_cowplot(15))



# data --------------------------------------------------------------------

# data from my original experiment, in March - 2020
resp

# data from Tanara's first experiment (August - 2021)
respT

# data from Tanara's second experiment (December - 2021)
respTnew


# join all datasets

joint_resp =resp %>% 
  bind_rows(respT %>% 
  mutate(Replicate = Replicate + 3)) %>% 
  bind_rows(respTnew %>% 
  mutate(Replicate = Replicate + 6))

remove = c('crp_prpB:K','crp_pyrE:K','gltA_gpt','gltA_hpt',
           'gpt_hpt','hpt_guaA','pyrE_gpt','pyrE_guaA','pyrE_hpt',
           'upp_udp_p(prpB)')

# remove genes
joint_resp = joint_resp %>% 
  filter(!(Genes %in% remove))


joint_resp.sum = joint_resp %>%
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

drug = 5
pos = position_jitter(width = 0.05, height = 0.05, seed = 1) # to plot names in jitter positions
joint_resp.sum %>% 
  filter(Drug == drug) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>%
  unite(Supp, Supplement, Supplement_mM) %>%
  spread(Supp, BW_norm) %>%
  ggplot(aes(x = Glucose_0, y = Glucose_10)) + 
  geom_hline(yintercept = 0, colour = 'grey30') +
  geom_vline(xintercept = 0, colour = 'grey30') +
  geom_point(position = pos, size = 2) + 
  geom_text_repel(aes(label = Genes), position = pos, max.overlaps = 100) +
  labs(title = expression(paste("5FU + Glucose effect on ", italic('C. elegans'),
                                " N2 phenotype", sep = '')),
       x = "Normalised median scores of *C. elegans* N2 phenotype",
       y = "Normalised median scores of *C. elegans* N2 phenotype with **Glucose**") +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.text = element_text(size = 6),
        axis.title.y = ggtext::element_markdown(),
        axis.title.x = ggtext::element_markdown()) + 
  guides(colour = guide_legend(override.aes = list(size = 4))) 


ggsave(file = here('Summary/resp_sublibrary', 'Scatter_sub_lib_JOINT_DATASETS_5uM.pdf'),
       height = 10, width = 12)



#### density plots for pathways ####

ecocyc_zip = ecocyc %>% 
  group_by(gene) %>% 
  mutate(pathway = str_c(pathway, collapse = ';')) %>% 
  distinct(gene, .keep_all = TRUE) %>% ungroup


# get the genes that are related to the relevant pathways
tca_genes = ecocyc %>% 
  filter(pathway %in% c(
    'superpathway of glyoxylate bypass and TCA',
    'superpathway of glycolysis, pyruvate dehydrogenase, TCA, and glyoxylate bypass',
    'TCA cycle I (prokaryotic)')) %>% 
  group_by(gene) %>% 
  mutate(pathway = str_c(pathway, collapse = ';')) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  mutate(pathway = 'TCA cycle')


pyr_genes = ecocyc %>% 
  filter(pathway %in% c('pyrimidine ribonucleotides de novo biosynthesis',
                        'superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis',
                        'superpathway of histidine, purine, and pyrimidine biosynthesis')) %>% 
  group_by(gene) %>% 
  mutate(pathway = str_c(pathway, collapse = ';')) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  mutate(pathway = 'Pyridimidine biosynthesis')


tca_genes$gene %in% pyr_genes$gene
pyr_genes$gene %in% tca_genes$gene

interesting_genes = tca_genes %>% 
  bind_rows(pyr_genes) %>% 
  rename(Genes = gene)

# plot with density for the two pathways
joint_resp.sum %>% 
  filter(Drug == drug) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>%
  unite(Supp, Supplement, Supplement_mM) %>%
  spread(Supp, BW_norm) %>%
  filter(!(abs(Glucose_10) < 0.6 & Glucose_0 < 1)) %>%
  left_join(interesting_genes ) %>%
  drop_na(pathway) %>%
  ggplot(aes(x = Glucose_0, y = Glucose_10)) + 
  geom_density_2d(aes(color = pathway)) +
  geom_hline(yintercept = 0, colour = 'grey30') +
  geom_vline(xintercept = 0, colour = 'grey30') +
  geom_point(position = pos, size = 2, aes(color = pathway)) 
  


# quadrants enrichment ----------------------------------------------------

# Using the merged data table, I'll separate the genes in different 
# spaces having at least a difference of 1 respect to the control:
#
# 1. upper-left, == resistant with glucose, normal without 
# 2. upper-right, == resistant in every condition
# 3. right, == resistant without glucose
# 4. bottom-right == sensitive with glucose, resistant without glu
# 5. bottom-left == sensitive in both conditions
# 6. bottom part == all sensitive to glucose
# 7. upper part == all resistant to glucose


drug = 5
joint_resp.wide = joint_resp.sum %>% 
  filter(Drug == drug) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>%
  unite(Supp, Supplement, Supplement_mM) %>%
  spread(Supp, BW_norm)


# get gene groups, modify thresholds if necessary
group_1 = joint_resp.wide %>% 
  filter(Glucose_0 <= 0.5, Glucose_10 > 0.5) 

group_2 = joint_resp.wide %>% 
  filter(Glucose_0 >= 1, Glucose_10 > 0.5) 

group_3 = joint_resp.wide %>% 
  filter(Glucose_0 >= 1) %>% 
  filter( Glucose_10 >= -0.5)

group_4 = joint_resp.wide %>% 
  filter(Glucose_0 >= 1, Glucose_10 < -0.5) 

group_5 = joint_resp.wide %>% 
  filter(Glucose_0 <= 0.5, Glucose_10 < -0.5) 

group_6 = joint_resp.wide %>% 
  filter(Glucose_10 < -0.5) 

group_7 = joint_resp.wide %>% 
  filter(Glucose_10 > 0.5) 


# calculate enrichment for the different sets
genes = group_1$Genes
g1.enrich = enrich(genes, db = ecocyc) %>% mutate(direction = 'group_1')

genes = group_2$Genes
g2.enrich = enrich(genes, db = ecocyc) %>% mutate(direction = 'group_2')

genes = group_3$Genes
g3.enrich = enrich(genes, db = ecocyc) %>% mutate(direction = 'group_3')

genes = group_4$Genes
g4.enrich = enrich(genes, db = ecocyc) %>% mutate(direction = 'group_4')

genes = group_5$Genes
g5.enrich = enrich(genes, db = ecocyc) %>% mutate(direction = 'group_5')

genes = group_6$Genes
g6.enrich = enrich(genes, db = ecocyc) %>% mutate(direction = 'group_6')

genes = group_7$Genes
g7.enrich = enrich(genes, db = ecocyc) %>% mutate(direction = 'group_7')


merge_enrich = g1.enrich %>% 
  bind_rows(g2.enrich, g3.enrich, g4.enrich, 
            g5.enrich, g6.enrich, g7.enrich) %>% 
  filter(pval <= 0.1)


list_of_datasets = list('Enrich resp sublibrary' = merge_enrich)

write.xlsx(list_of_datasets, 
           here('summary/resp_sublibrary','Ecocyc_enrichment_Resp_sublibrary.xlsx'), 
           colNames = T, rowNames = T) 


# get significant pathways and add a factor multiplier
sig_paths = merge_enrich %>% 
  filter(fdr < 0.05) %>%
  mutate(sig_factor = case_when(fdr.stars == '*' ~ 1,
                                fdr.stars == '**' ~ 2,
                                fdr.stars == '***' ~ 3))

 


expanded = merge_enrich %>% 
  filter(fdr <= 0.05) %>% 
  expand(categories) %>% 
  mutate(Glucose_0 = 0, Glucose_10 = 0)


sig_paths %>% 
  filter(categories %in% 'NADH to cytochrome bd oxidase electron transfer I')


# get scores depending on the quadrant from where it's been enriched 
paths_scores = sig_paths %>% 
  mutate(
    Glucose_10 = case_when(
      direction %in% c('group_4','group_5',
                       'group_6') ~ -1,
      direction %in% c('group_1','group_2',
                       'group_7') ~ 1,
      TRUE ~ 0),
    Glucose_0 = case_when(
      direction %in% c('group_3','group_4','group_2') ~ 1,
      TRUE ~ 0
    )
    ) %>% 
  group_by(categories) %>% 
  summarise(Glucose_0 = sum(Glucose_0),
            Glucose_10 = sum(Glucose_10))
  


# test plot
paths_scores %>% 
  filter(str_detect(categories,'TCA|pyrimidine')) %>% 
  mutate(categories = str_wrap(categories, width = 15)) %>% 
  ggplot(aes(x = Glucose_10, y = Glucose_0)) +
  geom_text(aes(label = categories), position = position_jitter(width = 0.2, height = 0.4)) 


