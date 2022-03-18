
# libraries ---------------------------------------------------------------

library(tidyverse)
library(readr)
library(readxl)
library(here)
library(glue)
library(ComplexHeatmap)
library(rstatix)
library(cowplot)


theme_set(theme_cowplot(15))


# read the data ----------------------
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



# version with p-values
enrich %>% 
  filter(p < 0.05) %>%
  mutate(logpval = -log10(p)) %>% 
  unite(sample, genotype, direction, remove = F) %>% 
  # mutate(Pathway = str_wrap(Pathway, width = 25)) %>% 
  mutate(sample = case_when(sample == 'p53_down' ~ 'p53\nDown',
                            sample == 'p53_up' ~ 'p53\nUp',
                            sample == 'wt_down' ~ 'WT\nDown',
                            sample == 'wt_up' ~ 'WT\nUp')) %>% 
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
  guides(color=guide_legend(title = '-log10(p-value)'), 
         size = guide_legend(title = '-log10(p-value)'),
         fill='none')

ggsave(here("summary", "metaboanalyst", "enrich", "Enrich_heatmap_pval.pdf"),
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









