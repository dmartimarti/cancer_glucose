
# libraries ---------------------------------------------------------------

library(tidyverse)
library(readr)
library(readxl)
library(here)
library(glue)
library(ComplexHeatmap)
library(rstatix)
library(cowplot)
library(extrafont)
# font_import()

theme_set(theme_cowplot(15))


metab = read_excel("Rdata.xlsx", sheet = "Rdata")



names(metab)

# make the tible longer

metab_clean = metab %>% 
  pivot_longer(cols = `(R)C(S)S-alliin`:xanthosine,
               names_to = 'Metabolite', values_to = 'value') %>% 
  separate(Sample, into = c('Genotype', 'Condition', 'Time'), sep = '_') %>% 
  separate(Code, into = c('ID','Sample','experiment'), sep = '_') %>% 
  mutate(experiment = case_when(is.na(experiment) ~ 1,
                                TRUE ~ 2),
         Genotype = factor(Genotype, levels = c('p53', 'WT')),
         Condition = factor(Condition, levels = c('Micit', 'Control')),
         # DATA IMPUTATION: 0 for 1s
         value = case_when(value == 0 ~ 1,
                           TRUE ~ value)) %>% 
  select(-ID,-Sample,-Time, Genotype, Condition, experiment, Metabolite, value)


### metaboanalyst csv ####
# create a version to be used in metaboanalyst
metab %>% 
  separate(Code, into = c('ID','Name','experiment'), sep = '_') %>% 
  mutate(experiment = case_when(is.na(experiment) ~ 1,
                                TRUE ~ 2)) %>% 
  separate(Sample, into = c('Genotype', 'Condition', 'Time'), sep = '_') %>% 
  filter(experiment == 1) %>% 
  select(-ID:-Name, -experiment,-Time) %>% 
  unite(Sample, Genotype:Condition, remove=F) %>% 
  filter(Genotype == 'WT') %>% 
  group_by(Genotype, Condition) %>% 
  mutate(ID = seq(1,n(),1), .before=Genotype) %>% 
  arrange(Genotype, Condition) %>% 
  select(-Sample) %>% 
  unite(Sample, Condition,ID, sep = '_', remove=F) %>% 
  ungroup %>% 
  select(-ID, -Genotype) %>% 
  write_csv(here('metaboanalyst_WT.csv'))

metab %>% 
  separate(Code, into = c('ID','Name','experiment'), sep = '_') %>% 
  mutate(experiment = case_when(is.na(experiment) ~ 1,
                                TRUE ~ 2)) %>% 
  separate(Sample, into = c('Genotype', 'Condition', 'Time'), sep = '_') %>% 
  filter(experiment == 1) %>% 
  select(-ID:-Name, -experiment,-Time) %>% 
  unite(Sample, Genotype:Condition, remove=F) %>% 
  filter(Genotype == 'p53') %>% 
  group_by(Genotype, Condition) %>% 
  mutate(ID = seq(1,n(),1), .before=Genotype) %>% 
  arrange(Genotype, Condition) %>% 
  select(-Sample) %>% 
  unite(Sample, Condition:ID, sep = '_', remove=F) %>% 
  write_csv(here('metaboanalyst_p53.csv'))



### multiomics ####

metab %>% 
  separate(Code, into = c('ID','Name','experiment'), sep = '_') %>% 
  mutate(experiment = case_when(is.na(experiment) ~ 1,
                                TRUE ~ 2)) %>% 
  separate(Sample, into = c('Genotype', 'Condition', 'Time'), sep = '_') %>% 
  filter(experiment == 1) %>% 
  select(-ID:-Name, -experiment,-Time) %>% 
  unite(Sample, Genotype:Condition, remove=F) %>% 
  filter(Genotype == 'WT') %>% 
  group_by(Genotype, Condition) %>% 
  mutate(ID = seq(1,n(),1), .before=Genotype) %>% 
  arrange(Genotype, Condition) %>% 
  select(-Sample) %>% 
  unite(Sample, Condition,ID, sep = '_', remove=F) %>% 
  ungroup %>% 
  select(-ID, -Genotype) %>% 
  pivot_longer(cols = `(R)C(S)S-alliin`:xanthosine,
               names_to = 'Metabolite', values_to = 'value') %>% 
  select(-Condition) %>% 
  pivot_wider(names_from = 'Sample', values_from = value) %>% 
  write_csv('metab_multiomics.csv')

# the original tubes were split into two different tubes, so I will separate 
# the data into two different experiments as well. I'll do the stats per 
# separate and only take the metabs that are significant in both experiments

## boxplots ####

# change the name of some metabolites that give problems with glue and filenames

metab_plots = metab_clean %>%
  mutate(Metabolite = str_replace_all(Metabolite, '/', ','))


query_metab = '2,3-Bisphosphoglyceric acid'

metab_plots %>% 
  unite(Sample, Genotype, Condition, sep = ', ', remove = FALSE) %>% 
  filter(Metabolite == query_metab) %>% 
  ggplot(aes(x = Sample, y = value, fill = Sample)) +
  geom_boxplot() + 
  labs(
    x = 'Sample',
    y = 'Normalised ion intensity',
    title = query_metab,
  ) +
  scale_fill_manual(
    values = c('#2661E0', # blue
               '#E0CD1D', # yellow
               '#25CC89', # green
               '#E34F32'  # orange
    )
  ) +
  geom_point(position = position_jitter(width = 0.1)) + 
  facet_wrap(~experiment) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(here('summary', 'individual_plots_split', 
            glue::glue('{query_metab}_boxplot.pdf')),
       height = 7, width = 11)  


metab_list = metab_plots %>% distinct(Metabolite) %>% pull(Metabolite)

# save plots with experiments splits
for (metabolite in metab_list){
  
  cat(glue::glue('Plotting metabolite {metabolite} \n\n'))
  
  metab_plots %>% 
    unite(Sample, Genotype, Condition, sep = ', ', remove = FALSE) %>% 
    filter(Metabolite == metabolite) %>% 
    ggplot(aes(x = Sample, y = value, fill = Sample)) +
    geom_boxplot() + 
    labs(
      x = 'Sample',
      y = 'Normalised ion intensity',
      title = metabolite,
    ) +
    scale_fill_manual(
      values = c('#2661E0', # blue
                 '#E0CD1D', # yellow
                 '#25CC89', # green
                 '#E34F32'  # orange
      )
    ) +
    geom_point(position = position_jitter(width = 0.1)) + 
    facet_wrap(~experiment) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(here('summary', 'individual_plots_split', 
              glue::glue('{metabolite}_boxplot.pdf')),
         height = 7, width = 11)  
}


# save plots all together

for (metabolite in metab_list){
  
  cat(glue::glue('Plotting metabolite {metabolite} \n\n'))
  
  metab_plots %>% 
    unite(Sample, Genotype, Condition, sep = ', ', remove = FALSE) %>% 
    filter(Metabolite == metabolite) %>% 
    ggplot(aes(x = Sample, y = value, fill = Sample)) +
    geom_boxplot() + 
    labs(
      x = 'Sample',
      y = 'Normalised ion intensity',
      title = metabolite,
    ) +
    scale_fill_manual(
      values = c('#2661E0', # blue
                 '#E0CD1D', # yellow
                 '#25CC89', # green
                 '#E34F32'  # orange
      )
    ) +
    geom_point(position = position_jitter(width = 0.1)) + 
    # facet_wrap(~experiment) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(here('summary', 'individual_plots', 
              glue::glue('{metabolite}_boxplot.pdf')),
         height = 7, width = 9)  
}






# stats -------------------------------------------------------------------

# WT micit vs control

wt_stats = metab_clean %>% 
  filter(Genotype == 'WT') %>% 
  mutate(value = log2(value)) %>%
  # mutate(Condition = as.factor) %>% 
  group_by(Metabolite, experiment) %>% 
  t_test(value ~ Condition, p.adjust.method = "fdr", detailed = T) 


weird_mets_wt = wt_stats %>% 
  filter(p < 0.5) %>% 
  group_by(Metabolite) %>% 
  count() %>% 
  arrange(n) %>% 
  filter(n == 1) %>% pull(Metabolite)

wt_stats = wt_stats %>% 
  filter(experiment == 1) %>% 
  arrange(estimate) %>% 
  select(-experiment, -`.y.`, -n1, -n2, -estimate1, -estimate2,
         -conf.low:-alternative) %>% 
  rename(log2FC = estimate) %>% 
  mutate(p.stars = gtools::stars.pval(p)) 

wt_stats %>% 
  write_csv(here('summary', 'WT_stats.csv'))



wt_stats %>% 
  filter(log2FC < 10) %>% 
  mutate(
         sig = case_when(abs(log2FC) > 0.5 & p < 0.05 ~ "significant",
                         TRUE ~ NA),
         label = case_when(sig == 'significant' ~ Metabolite,
                           TRUE ~ NA)) %>% 
  ggplot(aes(x = log2FC, y = -log10(p), fill = sig)) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_vline(xintercept = -0.5,  linetype = 'dashed') +
  geom_vline(xintercept = 0.5,  linetype = 'dashed') +
  geom_point(shape = 21, size = 3, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = label),
                           box.padding = 0.1) 

ggsave("summary/volcano_S6D.pdf", height = 10, width = 11)
  



# p53 micit vs control

p53_stats = metab_clean %>% 
  filter(Genotype == 'p53') %>% 
  mutate(value = log2(value)) %>% 
  # mutate(Condition = as.factor) %>% 
  group_by(Metabolite, experiment) %>% 
  t_test(value ~ Condition, p.adjust.method = "fdr", detailed = T) 


weird_mets_p53 = p53_stats %>% 
  filter(p < 0.5) %>% 
  group_by(Metabolite) %>% 
  count() %>% 
  arrange(n) %>% 
  filter(n == 1) %>% pull(Metabolite)

p53_stats = p53_stats %>% 
  filter(experiment == 1) %>% 
  arrange(statistic) %>% 
  select(-experiment, -`.y.`, -n1, -n2, -estimate1, -estimate2,
         -conf.low:-alternative) %>% 
  rename(log2FC = estimate) %>% 
  mutate(p.stars = gtools::stars.pval(p)) 

p53_stats %>% 
  write_csv(here('summary', 'p53_stats.csv'))


 
## get metabolite lists ####

# WT
wt_stats %>% 
  filter(p < 0.05, log2FC < 0) %>%
  select(Metabolite) %>% 
  write.table(here("summary", "metab_lists", "wt_down.txt"), 
              col.names = F, row.names = F, quote = F)

wt_stats %>% 
  filter(p < 0.05, log2FC > 0) %>%
  select(Metabolite) %>% 
  write.table(here("summary", "metab_lists", "wt_up.txt"), 
              col.names = F, row.names = F, quote = F)
# p53
p53_stats %>% 
  filter(p < 0.05, log2FC < 0) %>%
  select(Metabolite) %>% 
  write.table(here("summary", "metab_lists", "p53_down.txt"), 
              col.names = F, row.names = F, quote = F)

p53_stats %>% 
  filter(p < 0.05, log2FC > 0) %>%
  select(Metabolite) %>% 
  write.table(here("summary", "metab_lists", "p53_up.txt"), 
              col.names = F, row.names = F, quote = F)



# plot enrich -------------------------------------------------------------


read_enrich = function(data = "wt_down") {
  enrich = read_csv(
    here("summary", "metaboanalyst", "enrich", 
         data, "msea_ora_result.csv")) %>% 
    rename(Pathway = `...1`,
           p = `Raw p`,
           Holm_p = `Holm p`) %>% 
    mutate(dataset = data) %>% 
    separate(dataset, into = c('genotype', 'direction'), sep='_')
  
  return(enrich)
}

wt_down_enrich = read_enrich('wt_down')
wt_up_enrich = read_enrich('wt_up')
p53_down_enrich = read_enrich('p53_down')
p53_up_enrich = read_enrich('p53_up')


enrich = wt_down_enrich %>% 
  bind_rows(wt_up_enrich,
            p53_down_enrich,
            p53_up_enrich)

enrich %>% 
  write_csv(here("summary", "metaboanalyst", "enrich",
                 "enrichment_stats.csv"))
enrich %>% 
  write_csv(here("summary", 
                 "enrichment_stats.csv"))


paths = c('Amino Sugar Metabolism','Citric Acid Cycle','Alanine Metabolism',
          'Butyrate Metabolism', 'Folate Metabolism', 'Gluconeogenesis',
          'Glucose-Alanine Cycle', 'Glutamate Metabolism', 
          'Glutathione Metabolism', 'Glycolisis', 
          'Mitochondrial Electron Transport Chain', 'Purine Metabolism',
          'Pyrimidine Metabolism', 'Pyruvate Metabolism', 'Urea Cycle',
          'Marburg Effect', 'Glycerol Phosphate Shuttle'
          )

# version with p-values
enrich %>% 
  filter(p < 0.05) %>%
  mutate(logpval = -log10(p)) %>% 
  filter(logpval > 2) %>%
  filter(Pathway %in% paths) %>% 
  unite(sample, genotype, direction, remove = F) %>% 
  # mutate(Pathway = str_wrap(Pathway, width = 25)) %>% 
  mutate(sample = case_when(sample == 'p53_down' ~ 'p53\nDown',
                            sample == 'p53_up' ~ 'p53\nUp',
                            sample == 'wt_down' ~ 'WT\nDown',
                            sample == 'wt_up' ~ 'WT\nUp')) %>% 
  mutate(sample = factor(sample, levels = c('WT\nUp','WT\nDown',
                                            'p53\nUp','p53\nDown'))) %>% 
  filter(sample %in% c('WT\nUp', 'WT\nDown')) %>% 
  ggplot(aes(x = sample, y = Pathway, fill = logpval,
             size = logpval, color = logpval)) +
  # geom_tile() +
  geom_point() +
  labs(
    x = 'Pathway',
    y = 'Sample'
  ) +
  # guides(color = guide_legend(title='-log10(p-value)'),
  #        # size = 'none',
  #        fill = 'none') +
  scale_y_discrete(limits=rev) +
  # scale_fill_gradient(low = 'white', high = 'red') +
  scale_color_gradient(low = 'white', high = 'red', guide = 'legend') +
  scale_size_continuous(range = c(1, 7)) +
  guides(color=guide_legend(title = '-log10\n(p-value)'), 
         size = guide_legend(title = '-log10\n(p-value)'),
         fill='none') 

ggsave(here("summary", "metaboanalyst", "enrich", "Enrich_heatmap_pval_cutoff2.pdf"),
       height = 10, width = 8)


ggsave(here("summary", "metaboanalyst", "enrich", "Enrich_heatmap_poster.pdf"),
       height = 5.5, width = 7)


ggsave(here("summary", "metaboanalyst", "enrich", "Enrich_heatmap_poster2.pdf"),
       height = 5.5, width = 5)


### version for presentation ####


# version with p-values
enrich %>% 
  filter(p < 0.05) %>%
  mutate(logpval = -log10(p)) %>% 
  mutate(categories = cut(logpval, 
                          breaks = c(-Inf, 2, 4, Inf),
                          labels = c('2', '4', '6'))) %>% 
  unite(sample, genotype, direction, remove = F) %>% 
  # mutate(Pathway = str_wrap(Pathway, width = 25)) %>% 
  mutate(sample = case_when(sample == 'p53_down' ~ 'p53\nDown',
                            sample == 'p53_up' ~ 'p53\nUp',
                            sample == 'wt_down' ~ 'WT\nDown',
                            sample == 'wt_up' ~ 'WT\nUp')) %>% 
  ggplot(aes(x = sample, y = Pathway, fill = categories,
             size = categories,
             color = categories)) +
  # geom_tile() +
  geom_point() +
  labs(
    x = 'Pathway',
    y = 'Sample'
  ) +
  # guides(color = guide_legend(title='-log10(p-value)'),
  #        # size = 'none',
  #        fill = 'none') +
  scale_y_discrete(limits=rev) +
  # scale_fill_gradient(low = 'white', high = 'red') +
  scale_color_manual(values = c( '#E6C7B1','#E6945A', '#E6660C')) +
  # scale_color_gradient(low = 'white', high = 'red', guide = 'legend') +
  # scale_size_continuous(range = c(1, 7)) +
  guides(color=guide_legend(title = '-log10(p-value)'), 
         size = guide_legend(title = '-log10(p-value)'),
         fill='none')

ggsave(here("summary", "metaboanalyst", "enrich", "Enrich_heatmap_pval.pdf"),
       height = 10, width = 10)



# paper versions ----------------------------------------------------------
library(extrafont)
# version with p-values
enrich %>% 
  filter(p < 0.05) %>%
  mutate(logpval = -log10(p)) %>% 
  mutate(categories = cut(logpval, 
                          breaks = c(-Inf, 2, 4, Inf),
                          labels = c('2', '4', '6'))) %>% 
  unite(sample, genotype, direction, remove = F) %>% 
  # mutate(Pathway = str_wrap(Pathway, width = 25)) %>% 
  mutate(sample = case_when(sample == 'p53_down' ~ 'p53\nDown',
                            sample == 'p53_up' ~ 'p53\nUp',
                            sample == 'wt_down' ~ 'WT\nDown',
                            sample == 'wt_up' ~ 'WT\nUp')) %>% 
  filter(sample %in% c("WT\nDown", "WT\nUp")) %>% 
  ggplot(aes(x = sample, y = Pathway, fill = categories,
             size = categories,
             color = categories)) +
  # geom_tile() +
  geom_point() +
  labs(
    x = 'Pathway',
    y = 'Sample'
  ) +
  # guides(color = guide_legend(title='-log10(p-value)'),
  #        # size = 'none',
  #        fill = 'none') +
  scale_y_discrete(limits=rev) +
  # scale_fill_gradient(low = 'white', high = 'red') +
  scale_color_manual(values = c( '#E6C7B1','#E6945A', '#E6660C')) +
  # scale_color_gradient(low = 'white', high = 'red', guide = 'legend') +
  # scale_size_continuous(range = c(1, 7)) +
  guides(color=guide_legend(title = '-log10(p-value)'), 
         size = guide_legend(title = '-log10(p-value)'),
         fill='none') +
  theme_cowplot(font_family = "Arial")

ggsave(here("summary", "metaboanalyst", "enrich", "Enrich_heatmap_pval_WT.pdf"),
       height = 9, width = 7)


enrich %>% 
  filter(p < 0.05) %>%
  mutate(logpval = -log10(p)) %>% 
  mutate(categories = cut(logpval, 
                          breaks = c(-Inf, 2, 4, Inf),
                          labels = c('2', '4', '6'))) %>% 
  unite(sample, genotype, direction, remove = F) %>% 
  # mutate(Pathway = str_wrap(Pathway, width = 25)) %>% 
  mutate(sample = case_when(sample == 'p53_down' ~ 'p53\nDown',
                            sample == 'p53_up' ~ 'p53\nUp',
                            sample == 'wt_down' ~ 'WT\nDown',
                            sample == 'wt_up' ~ 'WT\nUp')) %>% 
  filter(sample %in% c("WT\nDown", "WT\nUp")) %>% 
  write_csv("summary/metaboanalyst/enrich/Enrich_heatmap_pval_WT.csv")


enrich %>% 
  filter(p < 0.05) %>%
  mutate(logpval = -log10(p)) %>% 
  mutate(categories = cut(logpval, 
                          breaks = c(-Inf, 2, 4, Inf),
                          labels = c('2', '4', '6'))) %>% 
  unite(sample, genotype, direction, remove = F) %>% 
  # mutate(Pathway = str_wrap(Pathway, width = 25)) %>% 
  mutate(sample = case_when(sample == 'p53_down' ~ 'p53\nDown',
                            sample == 'p53_up' ~ 'p53\nUp',
                            sample == 'wt_down' ~ 'WT\nDown',
                            sample == 'wt_up' ~ 'WT\nUp')) %>% 
  filter(sample %in% c("p53\nDown", "p53\nUp")) %>% 
  ggplot(aes(x = sample, y = Pathway, fill = categories,
             size = categories,
             color = categories)) +
  # geom_tile() +
  geom_point() +
  labs(
    x = 'Pathway',
    y = 'Sample'
  ) +
  # guides(color = guide_legend(title='-log10(p-value)'),
  #        # size = 'none',
  #        fill = 'none') +
  scale_y_discrete(limits=rev) +
  # scale_fill_gradient(low = 'white', high = 'red') +
  scale_color_manual(values = c( '#E6C7B1','#E6945A', '#E6660C')) +
  # scale_color_gradient(low = 'white', high = 'red', guide = 'legend') +
  # scale_size_continuous(range = c(1, 7)) +
  guides(color=guide_legend(title = '-log10(p-value)'), 
         size = guide_legend(title = '-log10(p-value)'),
         fill='none') +
  theme_cowplot(font_family = "Arial")

ggsave(here("summary", "metaboanalyst", "enrich", "Enrich_heatmap_pval_p53.pdf"),
       height = 10, width = 10)





enrich %>% 
  filter(FDR < 0.05) %>%
  mutate(logpval = -log10(FDR)) %>% 
  unite(sample, genotype, direction, remove = F) %>% 
  mutate(sample = case_when(sample == 'p53_down' ~ 'p53\nDown',
                            sample == 'p53_up' ~ 'p53\nUp',
                            sample == 'wt_down' ~ 'WT\nDown',
                            sample == 'wt_up' ~ 'WT\nUp')) %>% 
  ggplot(aes(x = sample, y = Pathway, fill = logpval)) +
  # geom_tile() +
  geom_point(aes(size = logpval, color = logpval)) +
  labs(
    x = 'Pathway',
    y = 'Sample'
  ) +
  # guides(color = guide_legend(title='-log10(FDR)'),
  #        size = 'none',
  #        fill = 'none') +
  scale_y_discrete(limits=rev) +
  scale_color_gradient(low = 'white', high = 'red')+
  scale_size_continuous(range = c(1, 7)) +
  guides(color=guide_legend(title = '-log10(p-value)'), 
         size = guide_legend(title = '-log10(p-value)'),
         fill='none')

ggsave(here("summary", "metaboanalyst", "enrich", "Enrich_heatmap_FDR.pdf"),
       height = 10, width = 9)



# boxplots with p-vals####

# change the name of some metabolites that give problems with glue and filenames

metab_plots = metab_clean %>% 
  filter(experiment == 1) %>% 
  mutate(value = log2(value)) %>% 
  left_join(wt_stats %>% 
              select(Metabolite, 
                     log2FC_wt = log2FC, 
                     p.stars_wt = p.stars)) %>% 
  left_join(p53_stats %>% 
              select(Metabolite, 
                     log2FC_p53 = log2FC, 
                     p.stars_p53 = p.stars))

minmax_vals = metab_plots %>% 
  group_by(Metabolite) %>% 
  summarise(min_val = min(value),
            max_val = max(value))

max_val = 24.32724
min_val = 19.90562

metab_plots %>% 
  unite(Sample, Genotype, Condition, sep = ', ', remove = FALSE) %>% 
  filter(Metabolite == query_metab) %>% 
  ggplot(aes(x = Sample, y = value, fill = Sample)) +
  geom_boxplot() + 
  labs(
    x = 'Sample',
    y = 'Normalised ion intensity (log2)',
    title = query_metab,
  ) +
  scale_fill_manual(
    values = c('#2661E0', # blue
               '#E0CD1D', # yellow
               '#25CC89', # green
               '#E34F32'  # orange
    )
  ) +
  geom_text(aes(label = p.stars_p53), 
            y = max_val * 1.05, x = 1.5, size = 5) +
  geom_text(aes(label = paste('log2FC: ',round(log2FC_p53,2))), 
            y = max_val * 1.03, x = 1.5, size = 4) +
  geom_text(aes(label = p.stars_wt), 
            y = max_val * 1.05, x = 3.5, size = 5) +
  geom_text(aes(label = paste('log2FC: ',round(log2FC_wt,2))), 
            y = max_val * 1.03, x = 3.5, size = 4) +
  ylim(min_val * 0.98, max_val * 1.05) +
  guides(fill = 'none') +
  geom_point(position = position_jitter(width = 0.1)) + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(here('summary', 'individual_plots_stats', 
            glue::glue('{query_metab}_boxplot.pdf')),
       height = 7, width = 8)  




metab_list = metab_plots %>% distinct(Metabolite) %>% pull(Metabolite)

for (metabolite in metab_list) {
  temp = metab_plots %>% 
    unite(Sample, Genotype, Condition, sep = ', ', remove = FALSE) %>% 
    filter(Metabolite == metabolite) 
  
  min_val = round(min(temp$value),2) 
  max_val = round(max(temp$value),2)
  
  cat(glue::glue('Plotting metabolite {metabolite} with min {min_val} and max {max_val}\n\n'))
  
  temp %>% 
    ggplot(aes(x = Sample, y = value, fill = Sample)) +
    geom_boxplot() + 
    labs(
      x = 'Sample',
      y = 'Normalised ion intensity (log2)',
      title = metabolite,
    ) +
    scale_fill_manual(
      values = c('#2661E0', # blue
                 '#E0CD1D', # yellow
                 '#25CC89', # green
                 '#E34F32'  # orange
      )
    ) +
    geom_text(aes(label = p.stars_p53), 
              y = max_val * 1.05, x = 1.5, size = 5) +
    geom_text(aes(label = paste('log2FC: ',round(log2FC_p53,2))), 
              y = max_val * 1.03, x = 1.5, size = 4) +
    geom_text(aes(label = p.stars_wt), 
              y = max_val * 1.05, x = 3.5, size = 5) +
    geom_text(aes(label = paste('log2FC: ',round(log2FC_wt,2))), 
              y = max_val * 1.03, x = 3.5, size = 4) +
    ylim(min_val * 0.97, max_val * 1.06) +
    guides(fill = 'none') +
    geom_point(position = position_jitter(width = 0.08)) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(here('summary', 'individual_plots_stats', 
              glue::glue('{metabolite}_boxplot.pdf')),
         height = 7, width = 8)
  
}




# plot for presentation ---------------------------------------------------

## Heatmap for selected paths ####

# Filipe wants a plot with all the metabolites from the Purine, Pyrimidine and
# TCA/Glycolysis pathways, p53 and WT
# I went to Metaboanalyst and parsed all the metabolites to get the metabolites
# within each pathway so I can match them

# read pathway data
library(readxl)
library(ComplexHeatmap)

paths = read_excel("pyr_pur_TCA_paths.xlsx")

paths = paths %>% 
  distinct(Pathway, Metabolite) 
# %>%   mutate(Metabolite = tolower(Metabolite))



metab_paths = metab_clean %>% 
  # mutate(Metabolite = tolower(Metabolite)) %>% 
  filter(Metabolite %in% paths$Metabolite) %>% 
  left_join(paths) %>% 
  drop_na() 
# %>% 
#   distinct(experiment, Genotype, 
#            Condition, Metabolite, Pathway, .keep_all = T)



metab_paths %>% 
  filter(Metabolite == 'ADP') %>% 
  filter(Pathway == 'Purine')

# test that we are plotting them right
metab_paths %>% 
  filter(Pathway == 'Pyrimidine') %>% 
  mutate(log_val = log2(value)) %>% 
  unite(cond, Condition, Genotype, remove = F) %>% 
  ggplot(aes(y = Metabolite, x = log_val, fill = cond)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.0), alpha = 0.4)


# create a matrix with purine compounds
pur_mat = metab_paths %>% 
  filter(Pathway == 'Purine') %>% 
  unite(Sample, Genotype, Condition, sep = ', ', remove = FALSE) %>% 
  mutate(Sample = factor(Sample, levels = c('WT, Control',
                                            'WT, Micit',
                                            'p53, Control',
                                            'p53, Micit'))) %>% 
  mutate(logval = log2(value)) %>% 
  group_by(Metabolite, Sample) %>%
  summarise(Mean = mean(logval)) %>% 
  mutate(
    z_score = (Mean - mean(Mean)) / sd(Mean)) %>% 
  select(Sample, Metabolite, z_score) %>% 
  pivot_wider(names_from = Sample, values_from = z_score) %>% 
  column_to_rownames('Metabolite') %>% 
  as.matrix()


# create a matrix with purine compounds
pyr_mat = metab_paths %>% 
  filter(Pathway == 'Pyrimidine') %>% 
  unite(Sample, Genotype, Condition, sep = ', ', remove = FALSE) %>% 
  mutate(Sample = factor(Sample, levels = c('WT, Control',
                                            'WT, Micit',
                                            'p53, Control',
                                            'p53, Micit'))) %>% 
  mutate(logval = log2(value)) %>% 
  group_by(Metabolite, Sample) %>%
  summarise(Mean = mean(logval)) %>% 
  mutate(
    z_score = (Mean - mean(Mean)) / sd(Mean)) %>% 
  select(Sample, Metabolite, z_score) %>% 
  pivot_wider(names_from = Sample, values_from = z_score) %>% 
  column_to_rownames('Metabolite') %>% 
  as.matrix()

# create a matrix with TCA compounds
tca_mat = metab_paths %>% 
  filter(Pathway == 'TCA/Glycolysis') %>% 
  unite(Sample, Genotype, Condition, sep = ', ', remove = FALSE) %>% 
  mutate(Sample = factor(Sample, levels = c('WT, Control',
                                            'WT, Micit',
                                            'p53, Control',
                                            'p53, Micit'))) %>% 
  mutate(logval = log2(value)) %>% 
  group_by(Metabolite, Sample) %>%
  summarise(Mean = mean(logval)) %>% 
  mutate(
    z_score = (Mean - mean(Mean)) / sd(Mean)) %>% 
  select(Sample, Metabolite, z_score) %>% 
  pivot_wider(names_from = Sample, values_from = z_score) %>% 
  column_to_rownames('Metabolite') %>% 
  as.matrix()


ha = HeatmapAnnotation(
  Condition = c(rep(c('Control','Micit'),2)), 
  Genotype = c('WT', 'WT', 'p53', 'p53'),
  col = list(
    Condition = c('Control' = '#25CC89', 'Micit' = '#E34F32' ),
    Genotype = c('WT' = '#2661E0', 'p53' = '#E0CD1D')
  )
)

ha_half = HeatmapAnnotation(
  Condition = c(rep(c('Control','Micit'),1)), 
  Genotype = c('WT', 'WT'),
  col = list(
    Condition = c('Control' = '#25CC89', 'Micit' = '#E34F32'),
    Genotype = c('WT' = '#2661E0')
  )
)

Heatmap(pur_mat[,1:2],
        name = 'Z-score',
        cluster_columns = FALSE,
        column_names_rot = 0,
        column_names_side = "top",
        column_names_centered = T,
        show_column_names = T,
        top_annotation = ha_half) 

quartz.save(file = here('summary', 'heatmap_purines_metabolomics.pdf'),
            type = 'pdf', dpi = 300, height = 7, width = 5)

Heatmap(pyr_mat[,1:2],
        name = 'Z-score',
        cluster_columns = FALSE,
        column_names_rot = 0,
        column_names_side = "top",
        column_names_centered = T,
        top_annotation = ha_half,
        show_column_names = T)

quartz.save(file = here('summary', 'heatmap_pyrimidines_metabolomics.pdf'),
            type = 'pdf', dpi = 300, height = 7, width = 6)


Heatmap(tca_mat,
        name = 'Z-score',
        cluster_columns = FALSE,
        column_names_rot = 0,
        column_names_side = "top",
        column_names_centered = T,
        top_annotation = ha,
        show_column_names = T)

quartz.save(file = here('summary', 'heatmap_TCA_metabolomics.pdf'),
            type = 'pdf', dpi = 300, height = 7, width = 7)

# ht_list = h1 %v% h2 %v% h3
# draw(ht_list)

## landscape mode -------
ha_r = rowAnnotation(
  Condition = c(rep(c('Control','Micit'),2)), 
  Genotype = c('WT', 'WT', 'p53', 'p53'),
  col = list(
    Condition = c('Control' = '#25CC89', 'Micit' = '#E34F32' ),
    Genotype = c('WT' = '#2661E0', 'p53' = '#E0CD1D')
  )
)

Heatmap(t(tca_mat),
        name = 'Z-score',
        cluster_columns = TRUE,
        cluster_rows = FALSE,
        column_names_rot = 90,
        column_names_side = "bottom",
        column_names_centered = T,
        # top_annotation = ha,
        left_annotation = ha_r,
        show_column_names = T)

quartz.save(file = here('summary', 'heatmap_TCA_metabolomics_landscape.pdf'),
            type = 'pdf', dpi = 300, height = 5, width = 7)


tca_mat %>% 
  as_tibble(rownames = "Metabolite") %>% 
  write_csv(here('summary', 'heatmap_TCA_metabolomics_landscape.csv'))


# Paper figures -----------------------------------------------------------

## pyrimidine & purine heatmaps ----------------------

# create a matrix with purine compounds
pyrpur_mat = metab_paths %>% 
  filter(Pathway %in% c('Purine', 'Pyrimidine')) %>% 
  unite(Sample, Genotype, Condition, sep = ', ', remove = FALSE) %>% 
  mutate(Sample = factor(Sample, levels = c('WT, Control',
                                            'WT, Micit',
                                            'p53, Control',
                                            'p53, Micit'))) %>% 
  arrange(Pathway) %>% 
  # filter(Genotype == 'WT') %>% 
  mutate(logval = log2(value)) %>% 
  group_by(Metabolite, Sample) %>%
  summarise(Mean = mean(logval)) %>% 
  mutate(
    z_score = (Mean - mean(Mean)) / sd(Mean)) %>% 
  select(Sample, Metabolite, z_score) %>% 
  pivot_wider(names_from = Sample, values_from = z_score) %>% 
  column_to_rownames('Metabolite') %>% 
  as.matrix()

# without z-score, just log2 values
pyrpur_mat = metab_paths %>%
  filter(Pathway %in% c('Purine', 'Pyrimidine')) %>%
  unite(Sample, Genotype, Condition, sep = ', ', remove = FALSE) %>%
  mutate(Sample = factor(Sample, levels = c('WT, Control',
                                            'WT, Micit',
                                            'p53, Control',
                                            'p53, Micit'))) %>%
  arrange(Pathway) %>%
  filter(Genotype == 'WT') %>%
  mutate(logval = log2(value)) %>%
  select(Sample, Metabolite, logval) %>%
  group_by(Metabolite, Sample) %>%
  summarise(Mean = mean(logval)) %>%
  pivot_wider(names_from = Sample, values_from = Mean) %>%
  column_to_rownames('Metabolite') %>%
  as.matrix()




Heatmap(pyrpur_mat[,1:2],
        name = 'Z-score',
        cluster_columns = FALSE,
        column_names_rot = 0,
        column_names_side = "top",
        column_names_centered = T,
        # top_annotation = ha_half,
        show_column_names = T) 

quartz.save(file = here('summary', 'heatmap_purines_metabolomics.pdf'),
            type = 'pdf', dpi = 300, height = 7, width = 5)



purine_order = c("ADP-ribose", "fumarate", "NAD", "NADH", "NADP", "GTP",
                 "ATP", "NADPH", "FAD", "GDP", "ADP", "adenosine",
                 "adenine", "dATP", "glutamate", "deoxyadenosine",
                 "urate", "GMP", "AMP", "xanthosine", "xanthine", 
                 "hypoxanthine", "inosine", "IMP", "aspartate", 
                 "glutamine", "glycine")

pyr_order = c("glutamine", "NADP", "ATP", "NADPH", "FAD", 
              "deoxycytidine", "orotate", "ADP", "UDP",
              "dTDP", "CMP", "deoxyuridine", "UMP", "uridine",
              "cytidine", "N-carbamoylaspartate")




### dot plot ----------

metab_paths %>%
  filter(Pathway %in% c('Purine', 'Pyrimidine')) %>%
  unite(Sample, Genotype, Condition, sep = ', ', remove = FALSE) %>%
  mutate(Sample = factor(Sample, levels = c('WT, Control',
                                            'WT, Micit',
                                            'p53, Control',
                                            'p53, Micit'))) %>%
  arrange(Pathway) %>%
  filter(Genotype == 'WT') %>%
  mutate(logval = log2(value)) %>%
  select(Sample, Metabolite, logval) %>%
  group_by(Metabolite, Sample) %>%
  summarise(Mean = mean(logval)) %>% 
  pivot_wider(names_from = Sample, values_from = Mean) %>% 
  mutate(diff = `WT, Micit` - `WT, Control`) %>% 
  ggplot(aes(x = diff, y = fct_reorder(Metabolite, diff))) +
  geom_vline(xintercept = 0, colour = 'grey60', alpha = 0.8) +
  geom_point(aes(color = diff), size = 4) +
  scale_colour_gradient2() +
  labs(
    x = "log2 Fold Change (Treatment vs Control)",
    y = NULL,
    color = "log2FC"
  ) +
  theme_cowplot(15, font_family = "Arial")

ggsave(here('summary', "dotplot_pyr_pur_log2FC.pdf"),
       height = 7, width = 5.5)


metab_paths %>%
  filter(Pathway %in% c('Purine', 'Pyrimidine')) %>%
  unite(Sample, Genotype, Condition, sep = ', ', remove = FALSE) %>%
  mutate(Sample = factor(Sample, levels = c('WT, Control',
                                            'WT, Micit',
                                            'p53, Control',
                                            'p53, Micit'))) %>%
  arrange(Pathway) %>%
  filter(Genotype == 'WT') %>%
  mutate(logval = log2(value)) %>%
  select(Sample, Metabolite, logval) %>%
  group_by(Metabolite, Sample) %>%
  summarise(Mean = mean(logval)) %>% 
  pivot_wider(names_from = Sample, values_from = Mean) %>% 
  mutate(diff = `WT, Micit` - `WT, Control`) %>% 
  write_csv("summary/dotplot_pyr_pur_log2FC_data.csv")


## selected boxplots #####

# helper function
plot_metab_boxplot = function(query_metab =  presentation_metabs[1]){
  
  metab_plots %>% 
    unite(Sample, Genotype, Condition, sep = ', ', remove = FALSE) %>% 
    filter(Metabolite %in% query_metab) %>% 
    mutate(Condition = factor(Condition, levels = c('Control', 'Micit')),
           Genotype = factor(Genotype, levels = c('WT', 'p53'))) %>% 
    ggplot(aes(x = Condition, y = value, fill = Condition)) +
    geom_boxplot(show.legend = F) + 
    labs(
      x = 'Sample',
      y = NULL,
      title = query_metab
    ) +
    scale_fill_manual(
      values = c(
        '#25CC89', # green
        '#E34F32'  # orange
      )
    ) +
    geom_point(position = position_jitter(width = 0.1),
               show.legend = F) + 
    facet_wrap(~Genotype) +
    # facet_grid(~Genotype*Metabolite) +
    theme_cowplot(15) +
    theme(plot.title = element_text(hjust = 0.5))
}


plot_metab_boxplot(query_metab = presentation_metabs[3])

presentation_metabs = c('NAD', 'NADH', 'glycerol 3-phosphate')

for (metab in presentation_metabs){
  
  plot_metab_boxplot(query_metab = metab)
  
  ggsave(here('summary', 'selected_boxplots', 
              glue::glue('{metab}_boxplot.pdf')),
         height = 5, width = 7)  
}



# all 3 plots together because why not

library(ggpubr)

pl1 = plot_metab_boxplot(query_metab = presentation_metabs[1]) +
  labs(y = 'Normalised ion intensity (log2)')
pl2 = plot_metab_boxplot(query_metab = presentation_metabs[2])
pl3 = plot_metab_boxplot(query_metab = presentation_metabs[3]) 

figure = ggarrange(pl1, pl2, pl3,
          ncol = 3, nrow = 1)
figure

ggsave(here('summary', 'selected_boxplots', 
            'selected_presentation_boxplot.pdf'),
       height = 5, width = 13)  

## stats for boxplots ####

# WT micit vs control

wt_select = metab_clean %>% 
  filter(Metabolite %in% presentation_metabs) %>% 
  filter(Genotype == 'WT') %>% 
  mutate(value = log2(value)) %>%
  # mutate(Condition = as.factor) %>% 
  group_by(Metabolite) %>% 
  t_test(value ~ Condition, p.adjust.method = "fdr", detailed = T) %>% 
  mutate(test = 'WT',.before = Metabolite)

p_select = metab_clean %>% 
  filter(Metabolite %in% presentation_metabs) %>% 
  filter(Genotype == 'p53') %>% 
  mutate(value = log2(value)) %>%
  # mutate(Condition = as.factor) %>% 
  group_by(Metabolite) %>% 
  t_test(value ~ Condition, p.adjust.method = "fdr", detailed = T)  %>% 
  mutate(test = 'p53',.before = Metabolite)


interaction_terms = metab_clean %>% 
  filter(Metabolite %in% presentation_metabs) %>% 
  # filter(Genotype == 'p53') %>% 
  mutate(value = log2(value)) %>%
  group_by(Metabolite) %>% 
  anova_test(value ~ Condition * Genotype, detailed = T) 


wt_select %>% 
  bind_rows(p_select) %>% 
  write_csv(here('summary', 'selected_boxplots',
                 'stats_t_test.csv'))


interaction_terms %>% 
  write_csv(here('summary', 'selected_boxplots',
                 'interaction_terms.csv'))




