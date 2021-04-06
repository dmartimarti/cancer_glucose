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

# options(width = 150)

theme_set(theme_classic())


# Read data ##################

# read and transform data
resp = read_xlsx('Sub_library/Resp_hits_200317.xlsx', sheet = 'Summary_good') 

# remove bad data
resp = resp %>% 
  filter(Genes != 'EMPTY') %>% 
  filter(Genes != 'tyrA') # gene didn't have worms for some samples

resp = resp %>%
  gather(Drug, Score, `0`, `5`) %>% # make table long
  unite(ID, Genes, Supplement, Replicate, remove = FALSE) %>% 
  mutate(Replicate = as.numeric(Replicate),
         Supplement_mM = as.factor(Supplement_mM),
         Supplement = as.factor(Supplement),
         Genes = as.factor(Genes),
         ID = as.factor(ID),
         Drug = as.factor(Drug)) %>%
  select(Genes, Well, ID, Replicate, Supplement, Supplement_mM, Drug, Score)


# table with genotypes and pathways

# data summary
resp.sum = resp %>%
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
resp.sum %>% 
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


ggsave(file = here('Summary', 'Scatter_sub_lib_RESP_5uM.pdf'),
       height = 10, width = 12)




# compare old and new glu -------------------------------------------------


library(patchwork)


drug = 5
pos = position_jitter(width = 0.05, height = 0.05, seed = 1) # to plot names in jitter positions
new = resp.sum %>% 
  filter(Drug == drug) %>% 
  filter(Genes %in% unique(as.character(glu.sum$Genes))) %>% 
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
  guides(colour = guide_legend(override.aes = list(size = 4))) # make lengend points larger

old = glu.sum %>% 
  filter(Drug == drug) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>%
  unite(Supp, Supplement, Supplement_mM) %>%
  spread(Supp, BW_norm) %>%
  left_join(pathways) %>% 
  ggplot(aes(x = Glucose_0, y = Glucose_10)) + 
  geom_hline(yintercept = 0, colour = 'grey30') +
  geom_vline(xintercept = 0, colour = 'grey30') +
  geom_point(aes(colour = Pathway),position = pos, size = 2) + 
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
  guides(colour = guide_legend(override.aes = list(size = 4))) # make lengend points larger


new + old

ggsave(file = here('Summary', 'Scatter_sub_COMPARISON.pdf'),
       height = 10, width = 24)






# merge both datasets -----------------------------------------------------


glu_res = glucose %>% select(-Pathway) %>%
  filter(Drug == 5) %>%
  rbind(resp %>% 
          mutate(Replicate = Replicate + 10)
        )

# data summary
glu_res.sum = glu_res %>%
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
glu_res.sum %>% 
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




ggsave(file = here('Summary', 'Scatter_sub_lib_RESP_5uM_MERGED.pdf'),
       height = 10, width = 12)







# quadrants enrichment ----------------------------------------------------

# Using the merged data table, I'll separate the genes in 5 different 
# spaces: 
# 1. upper-left, == resistant with glu, normal without 
# 2. upper-right, == resistant in every condition
# 3. middle-right, == resistant without glu, normal with glu
# 4. bottom-right == sensitive with glucose, resistant without glu
# 5. bottom-left == sensitive in both conditions
drug = 5
glu_res.wide = glu_res.sum %>% 
  filter(Drug == drug) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>%
  unite(Supp, Supplement, Supplement_mM) %>%
  spread(Supp, BW_norm)

# get gene groups, modify thresholds if necessary
group_1 = glu_res.wide %>% 
  filter(Glucose_0 <= 0.5, Glucose_10 >= 0.5) 

group_2 = glu_res.wide %>% 
  filter(Glucose_0 <= 0.5, Glucose_10 <= -0.5) 

group_3 = glu_res.wide %>% 
  filter(Glucose_0 >= 1.5, Glucose_10 <= 0.5) %>% 
  filter( Glucose_10 >= -0.5)

group_4 = glu_res.wide %>% 
  filter(Glucose_0 >= 1.5, Glucose_10 >= 0.5) 

group_5 = glu_res.wide %>% 
  filter(Glucose_0 >= 1.5, Glucose_10 <= -0.5) 



# let's do enrichment from those
# KEGG
# GO function
# EcoCyc 

# read the data table with gene info

ecocyc = read_csv("EcoCyc_genes_pathways.csv")

##### Hypergeometric function ########


# hypergeometric test function

enrich = function(gene, db){
  # initiate variables
  pval = c()
  m_total = c()
  x_total = c()
  k_total = c()
  gene_in_cat = c()
  db = as.data.frame(db)
  cats = unique(db[,2])
  
  for (cat in cats){
    subcat = db[db[,2] == cat,]
    N = (db %>% distinct(.[,1]) %>% count())$n
    m = dim(subcat)[1]
    n = N - m
    x = sum(gene %in% subcat[,1])
    k = sum(gene %in% db[,1]) # genes with at least 1 annotation!
    p = phyper(q=x-1, m=m, n=n, k=k, lower.tail=FALSE)
    
    # save variables
    m_total = c(m_total, m)
    x_total = c(x_total, x)
    k_total = c(k_total, k)
    gene_in_cat = c(gene_in_cat)
    pval = c(pval, p)
  }
  
  # build the table
  table = tibble(categories = cats, N = N, elm_in_cat = m_total, gene_in_cat = x_total, k_tot = k, pval = pval) %>% 
    mutate(
      p.stars = gtools::stars.pval(pval),
      fdr = p.adjust(pval, method = 'fdr'),
      fdr.stars = gtools::stars.pval(fdr)
    ) %>% 
    arrange(pval)
  
  return(table)
  
}

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


merge_enrich = g1.enrich %>% 
  bind_rows(g2.enrich,g3.enrich,g4.enrich,g5.enrich) %>% 
  filter(pval <= 0.1)


list_of_datasets = list('Enrichment respiration sublibrary' = merge_enrich)

write.xlsx(list_of_datasets, 
           here('summary','Ecocyc_enrichment_Resp_sublibrary.xlsx'), 
           colNames = T, rowNames = T) 


# plot

expanded = merge_enrich %>% 
  filter(fdr <= 0.1) %>% 
  select(categories, fdr, fdr.stars, direction) %>% 
  expand(categories,direction)

expanded = expanded %>% left_join(merge_enrich %>% 
                         filter(fdr <= 0.05) %>% 
                         select(categories, fdr, fdr.stars, direction)) %>% 
  replace_na(list(fdr = 1, p.stars = ''))







# plots -------------------------------------------------------------------

removals = c('superpathway of sulfate assimilation and cysteine biosynthesis',
             'sulfate activation for sulfonation',
             'pyruvate decarboxylation to acetyl CoA I',
             'inosine-5\'-phosphate biosynthesis I',
             'guanosine ribonucleotides de novo biosynthesis',
             'citrate degradation', 'adenosine nucleotides degradation II',
             '2-oxoglutarate decarboxylation to succinyl-CoA')

# plot KEIO
expanded %>% 
  filter(!(categories %in% removals)) %>% 
  mutate(Category = cut(fdr , 
                        breaks=c(-Inf, 0.001,  0.01, 0.05, 0.1, Inf), 
                        labels=c("<0.001","<0.01","<0.05", '<0.1', 'NS') )) %>% 
  mutate(direction = case_when(direction == 'group_1' ~ 'Upper-Left',
                               direction == 'group_2' ~ 'Bottom-Left',
                               direction == 'group_3' ~ 'Mid-right',
                               direction == 'group_4' ~ 'Upper-Right',
                               direction == 'group_5' ~ 'Bottom-Right')) %>% 
  ggplot(aes(y = categories, x = direction, fill = Category)) +
  geom_tile() +
  labs(
    x = 'Quadrant',
    y = 'Terms'
    ) +
  # Case for all fdr values
  # scale_fill_manual(values = c('#2229F0', '#2292F0', '#22E4F0', '#C7F0F0', 'white'),
  #                   name = 'FDR') +
  # case if one case is missing
  scale_fill_manual(values = c('#2229F0', '#2292F0', '#22E4F0',  'white'),
                    name = 'FDR') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here('summary', 'EcoCyc_enrich_resp_sublib.pdf'), height =  7, width = 13)


### clean pathways to include them in the scatterplot

pathway_select = c("chorismate biosynthesis I",
                   "folate transformations III (E. coli)",
                   "glycine cleavage",
                   "superpathway of glycolysis, pyruvate dehydrogenase, TCA, and glyoxylate bypass",
                   "superpathway of histidine, purine, and pyrimidine biosynthesis",
                   "superpathway of pyridoxal 5'-phosphate biosynthesis and salvage",
                   "superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis",
                   "UMP biosynthesis I")


# set coordinates for categories
cat_coord = expanded %>% 
  filter(categories %in%  pathway_select) %>% 
  filter(fdr <= 0.05) %>% 
  mutate(Glucose_0 = case_when(direction == 'group_2' ~ 0.25,
                               direction == 'group_3' ~ 1.75,
                               direction == 'group_4' ~ 1.75,
                               direction == 'group_5' ~ 1.75),
         Glucose_10 = case_when(direction == 'group_2' ~ -0.75,
                                direction == 'group_3' ~ 0,
                                direction == 'group_4' ~ 0.75,
                                direction == 'group_5' ~ -0.75))






# dont like this version
glu_res.wide %>% 
  bind_rows(cat_coord %>% select(Glucose_0,Glucose_10,categories)) %>%
  ggplot(aes(x = Glucose_0, y = Glucose_10)) + 
  geom_hline(yintercept = 0, colour = 'grey30') +
  geom_vline(xintercept = 0, colour = 'grey30') +
  geom_point(position = pos, size = 2) + 
  # geom_text_repel(aes(label = Genes), position = pos, max.overlaps = 100) +
  geom_text_repel(aes(label = categories), position = pos, max.overlaps = 100) +
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





genes_to_include = factor(c('upp', 'guaB', 'guaA', 'gcvP', 'gcvH', 'gcvT', 
                     'nuoI', 'purH',  'pdxB', 'pdxA', 'gpp', 'glyA',
                     'gltA', 'atpH', 'sucD', 'atpA','folM',
                     'lpd', 'gpt', 'nuoK', 'nuoL', 'purK', 'purN', 'pyrE',
                     'pyrD', 'pyrB', 'holC', 'sucA', 'purL','pyrF'))

glu_res.wide %>% 
  mutate(
    Genes = as.character(Genes),
    Genes = case_when(Genes %in% genes_to_include ~ Genes,
                           TRUE ~ '')) %>% 
  ggplot(aes(x = Glucose_0, y = Glucose_10)) + 
  annotate("rect", xmin = 1.3, xmax = 2.2, ymin = -0.5, ymax = 0.5, 
           alpha = .3, fill = 'blue') +
  annotate("rect", xmin = 1.3, xmax = 2.2, ymin = 0.5, ymax = 1.2, 
           alpha = .3, fill = 'green') +
  annotate("rect", xmin = 1.3, xmax = 2.2, ymin = -1.2, ymax = -0.5, 
           alpha = .3, fill = 'red') +
  annotate("rect", xmin = -0.2, xmax = 0.65, ymin = -1.2, ymax = -0.45, 
           alpha = .3, fill = 'yellow') +
  annotate("label", x = c(0.2, 1.65, 1.6,1.65), 
           y = c(-0.65, -0.75,-0.2,0.8), 
           label = c("UMP biosynthesis \nPyrimidine (deoxy)ribonucleotide\n de novo synthesis \nChorismate synthesis",
                     'Glycolisis, Pyruvate and \nTCA superpathway \nPyridoxal 5\'-phosphate biosynthesis',
                     "Glycine cleavage \nFolate transformations III",
                     "Glycine cleavage \nFolate transformations III \nPurine de novo biosynthesis II"),
           alpha = 0.5) +
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


ggsave(file = here('Summary', 'Scatter_sub_lib_with_Enrichment.pdf'),
       height = 10, width = 12)



# small summary of the above figure: I merged Leo's and my dataset into one,
# then I extracted the genes from each one of the important places in the map,
# and did a hypergeometric test with EcoCyc categories. I summarised the cats 
# in the most useful/informative, and represented them in the plot











# quadrants 2 enrichment ----------------------------------------------------

# Using the merged data table, I'll separate the genes in 5 different 
# spaces: 
# 1. upper == resistant with glu
# 2. right == resistant in normal conditions
# 3. bottom == sensitive with glu
# 4. center == slightly resistant in normal conditions

drug = 5
glu_res.wide = glu_res.sum %>% 
  filter(Drug == drug) %>% 
  select(Supplement, Supplement_mM, Genes, BW_norm) %>%
  unite(Supp, Supplement, Supplement_mM) %>%
  spread(Supp, BW_norm)

# get gene groups, modify thresholds if necessary
group_1 = glu_res.wide %>% 
  filter(Glucose_10 >= 0.5) 

group_2 = glu_res.wide %>% 
  filter(Glucose_0 >= 1.5) 

group_3 = glu_res.wide %>% 
  filter(Glucose_10 <= -0.5)

group_4 = glu_res.wide %>% 
  filter(Glucose_0 >= 0.5, Glucose_0 <= 1.5) 




# let's do enrichment from those
# KEGG
# GO function
# EcoCyc 

# read the data table with gene info

ecocyc = read_csv("EcoCyc_genes_pathways.csv")

# calculate enrichment for the different sets
genes = group_1$Genes
g1.enrich = enrich(genes, db = ecocyc) %>% mutate(direction = 'group_1')

genes = group_2$Genes
g2.enrich = enrich(genes, db = ecocyc) %>% mutate(direction = 'group_2')

genes = group_3$Genes
g3.enrich = enrich(genes, db = ecocyc) %>% mutate(direction = 'group_3')

genes = group_4$Genes
g4.enrich = enrich(genes, db = ecocyc) %>% mutate(direction = 'group_4')


merge_enrich = g1.enrich %>% 
  bind_rows(g2.enrich,g3.enrich,g4.enrich) %>% 
  filter(pval <= 0.1)


list_of_datasets = list('Enrichment respiration sublibrary' = merge_enrich)

write.xlsx(list_of_datasets, 
           here('summary','Ecocyc_enrichment_Resp_sublibrary_v2.xlsx'), 
           colNames = T, rowNames = T) 


# plot

expanded = merge_enrich %>% 
  filter(fdr <= 0.05) %>% 
  select(categories, fdr, fdr.stars, direction) %>% 
  expand(categories,direction)

expanded = expanded %>% left_join(merge_enrich %>% 
                                    filter(fdr <= 0.05) %>% 
                                    select(categories, fdr, fdr.stars, direction)) %>% 
  replace_na(list(fdr = 1, p.stars = ''))







# plots -------------------------------------------------------------------

removals = c('superpathway of sulfate assimilation and cysteine biosynthesis',
             'sulfate activation for sulfonation',
             'pyruvate decarboxylation to acetyl CoA I',
             'inosine-5\'-phosphate biosynthesis I',
             'guanosine ribonucleotides de novo biosynthesis',
             'citrate degradation', 'adenosine nucleotides degradation II',
             '2-oxoglutarate decarboxylation to succinyl-CoA')

# 1. upper == resistant with glu
# 2. right == resistant in normal conditions
# 3. bottom == sensitive with glu
# 4. center == slightly resistant in normal conditions


# plot enrichment
expanded %>% 
  # filter(!(categories %in% removals)) %>% 
  mutate(Category = cut(fdr , 
                        breaks=c(-Inf, 0.001,  0.01, 0.05, 0.1, Inf), 
                        labels=c("<0.001","<0.01","<0.05", '<0.1', 'NS') )) %>% 
  mutate(direction = case_when(direction == 'group_1' ~ 'Resistant glucose',
                               direction == 'group_2' ~ 'Resistant normal',
                               direction == 'group_3' ~ 'Sensitive glucose',
                               direction == 'group_4' ~ '(slightly) resistant normal')) %>% 
  ggplot(aes(y = categories, x = direction, fill = Category)) +
  geom_tile() +
  labs(
    x = 'Quadrant',
    y = 'Terms'
  ) +
  # Case for all fdr values
  # scale_fill_manual(values = c('#2229F0', '#2292F0', '#22E4F0', '#C7F0F0', 'white'),
  #                   name = 'FDR') +
  # case if one case is missing
  scale_fill_manual(values = c('#2229F0', '#2292F0', '#22E4F0',  'white'),
                    name = 'FDR') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here('summary', 'EcoCyc_enrich_resp_sublib_v2.pdf'), height =  10, width = 13)







