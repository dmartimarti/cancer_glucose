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
                       'group_6') ~ sig_factor * -1,
      direction %in% c('group_1','group_2',
                       'group_7') ~ sig_factor * 1,
      TRUE ~ 0),
    Glucose_0 = case_when(
      direction %in% c('group_3','group_4','group_2') ~ sig_factor * 1,
      TRUE ~ 0
    )
    ) %>% 
  group_by(categories) %>% 
  summarise(Glucose_0 = sum(Glucose_0),
            Glucose_10 = sum(Glucose_10))
  


# Plot pathways -----------------------------------------------------------



# test plot with all pathways
paths_scores %>% 
  # filter(str_detect(categories,'TCA|pyrimidine')) %>% 
  mutate(categories = str_wrap(categories, width = 25)) %>% 
  ggplot(aes(x = Glucose_10, y = Glucose_0)) +
  geom_text(aes(label = categories), position = position_jitter(width = 0.2, height = 0.4)) 


select_paths = c('superpathway of histidine, purine, and pyrimidine biosynthesis',
                 'superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis',
                 # 'superpathway of pyrimidine ribonucleotides de novo biosynthesis',
                 'superpathway of histidine, purine, and pyrimidine biosynthesis',
                 'UMP biosynthesis I',
                 'superpathway of glycolysis, pyruvate dehydrogenase, TCA, and glyoxylate bypass',
                 'TCA cycle I (prokaryotic)',
                 'superpathway of purine nucleotides de novo biosynthesis II',
                 'NADH to cytochrome bo oxidase electron transfer I',
                 'NADH to cytochrome bd oxidase electron transfer I',
                 'pentose phosphate pathway',
                 'pentose phosphate pathway (non-oxidative branch)')


# getting there
paths_scores %>% 
  filter(categories %in% select_paths) %>% 
  # filter(str_detect(categories,'TCA|pyrimidine')) %>% 
  mutate(categories = str_wrap(categories, width = 20)) %>% 
  ggplot(aes(x = Glucose_10, y = Glucose_0)) +
  # geom_label(aes(label = categories), 
  #            position = position_jitter(width = 0.2, height = 0.4)) +
  geom_label_repel(aes(label = categories), 
                   min.segment.length = Inf,
                   arrow = NULL)

paths_scores_tidy = paths_scores %>% 
  filter(categories %in% select_paths) %>% 
  # filter(str_detect(categories,'TCA|pyrimidine')) %>% 
  mutate(categories = str_wrap(categories, width = 20))


expand_grid(Glucose_0 = seq(from = 0, by = 0.01, to = 3),
            Glucose_10 = seq(from = -6, by = 0.01, to = 2)) %>% 
  mutate(fill = Glucose_0 * Glucose_10 + (Glucose_10)) %>% 
  ggplot(aes(y = Glucose_0,
             x = Glucose_10)) +
  geom_tile(aes(fill = fill),
            width=0.05) +
  scale_fill_gradient2() +
  geom_label_repel(
    data = paths_scores_tidy,
    aes(x = Glucose_10, 
        y = Glucose_0,
        label = categories),
    min.segment.length = Inf,
    alpha = 0.3,
    seed = 1234) +
  geom_label_repel(
    data = paths_scores_tidy,
    aes(x = Glucose_10, 
        y = Glucose_0,
        label = categories),
    min.segment.length = Inf,
    label.size = NA,
    fill = NA,
    seed = 1234) +
  guides(fill = 'none') +
  labs(
    x = 'Sensitivity to 5FU + Glucose score',
    y = 'Sensitivity to 5FU score'
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

ggsave(here('summary/resp_sublibrary', 'Pathways_sensitivity.pdf'),
       height = 9, width = 11)



# Fixing the background colors --------------------------------------------


background = expand_grid(Glucose_0 = seq(from = 0, by = 0.01, to = 3),
            Glucose_10 = seq(from = -6, by = 0.01, to = 2)) %>% 
  mutate(fill = Glucose_0 * Glucose_10 + (Glucose_10)) 



min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}



# this code will do a modified version of min/max norm
vect = background$Glucose_10
new_vect = c()
for (i in vect){
  if (i < 0){
    x = (-(abs(i)/abs(min(vect))))
    new_vect = c(new_vect,x)
  } else if(i > 0) {
    x = (i/max(vect))
    new_vect = c(new_vect,x)
  } else if(i == 0){
    x = 0
    new_vect = c(new_vect,x)
  }
}

background['Glucose_10_fill'] = new_vect

background = background %>% 
  mutate(Glucose_0_fill = min_max_norm(Glucose_0))


background %>% 
  mutate(new_fill = (Glucose_0_fill * Glucose_10_fill) + Glucose_10_fill) %>% 
  ggplot(aes(y = Glucose_0,
             x = Glucose_10)) +
  geom_tile(aes(fill = new_fill),
            alpha = 1,
            width = 0.05) +
  scale_fill_gradient2() +
  geom_label_repel(
    data = paths_scores_tidy,
    aes(x = Glucose_10, 
        y = Glucose_0,
        label = categories),
    min.segment.length = Inf,
    alpha = 0.3,
    seed = 1234) +
  geom_label_repel(
    data = paths_scores_tidy,
    aes(x = Glucose_10, 
        y = Glucose_0,
        label = categories),
    min.segment.length = Inf,
    label.size = NA,
    fill = NA,
    seed = 1234) +
  guides(fill = 'none') +
  labs(
    x = 'Sensitivity to 5FU + Glucose score',
    y = 'Sensitivity to 5FU score'
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))


ggsave(here('summary/resp_sublibrary', 'Pathways_sensitivity_v2.pdf'),
       height = 9, width = 11)



#### version 3 ####

glu_targets = c('TCA cycle I\n(prokaryotic)',
'superpathway of\nglycolysis, pyruvate\ndehydrogenase,\nTCA, and glyoxylate\nbypass')

pyr_targets = c('superpathway\nof pyrimidine\ndeoxyribonucleotides\nde novo biosynthesis')

paths_scores_tidy_filter = paths_scores_tidy %>% 
  mutate(color = case_when(categories %in% glu_targets ~ 'blue',
                           categories %in% pyr_targets ~ 'red',
                           TRUE ~ 'grey30'),
         size = case_when(categories %in% glu_targets | categories %in% pyr_targets ~ 5,
                          TRUE ~ 3)) %>% 
  filter(!(categories %in% c('superpathway of\nhistidine, purine,\nand pyrimidine\nbiosynthesis',
                           'UMP biosynthesis I')))


background %>% 
  # mutate(new_fill = (Glucose_0 * Glucose_10) ) %>% 
  mutate(new_fill = (Glucose_0_fill * Glucose_10_fill) ) %>% 
  ggplot(aes(y = Glucose_0,
             x = Glucose_10)) +
  geom_tile(aes(fill = new_fill),
            alpha = 1,
            width = .05) +
  scale_fill_gradient2(
    high = 'red',
    mid = 'grey95',
    low = 'blue') +
  geom_label_repel(
    data = paths_scores_tidy_filter,
    aes(x = Glucose_10, 
        y = Glucose_0,
        label = categories,
        size = size),
    min.segment.length = Inf,
    alpha = 0.3,
    seed = 1234) +
  geom_label_repel(
    data = paths_scores_tidy_filter,
    aes(x = Glucose_10, 
        y = Glucose_0,
        label = categories,
        size = size,
        color = color),
    min.segment.length = Inf,
    label.size = NA,
    fill = NA,
    seed = 1234) +
  guides(fill = 'none',
         size = 'none',
         color = 'none') +
  labs(
    x = 'Sensitivity to 5FU + Glucose score',
    y = 'Sensitivity to 5FU score'
  ) +
  scale_size(range = c(3,6)) +
  scale_color_manual(values = c('darkred', 'grey50', 'darkblue')) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_cowplot(20)



ggsave(here('summary/resp_sublibrary', 'Pathways_sensitivity_v3.pdf'),
       height = 9, width = 11)





background %>% 
  mutate(new_fill = (Glucose_0_fill * Glucose_10_fill) ) %>% 
  filter(Glucose_10 >= 0) %>% 
  ggplot(aes(y = Glucose_0,
             x = Glucose_10)) +
  geom_tile(aes(fill = new_fill),
            alpha = 1) +
  scale_fill_gradient2(
    high = 'red',
    mid = 'grey95',
    low = 'blue')


