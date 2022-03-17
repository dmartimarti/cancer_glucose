
# libraries ---------------------------------------------------------------

library(tidyverse)
library(readxl)
library(here)
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
         Genotype = factor(Genotype, levels = c('WT', 'p53')),
         Condition = factor(Condition, levels = c('Control', 'Micit')),
         # DATA IMPUTATION: 0 for 1s
         alue = case_when(value == 0 ~ 1,
                           TRUE ~ value)) %>% 
  select(-ID,-Sample,-Time, Genotype, Condition, experiment, Metabolite, value)


# the original tubes were split into two different tubes, so I will separate 
# the data into two different experiments as well. I'll do the stats per 
# separate and only take the metabs that are significant in both experiments

## boxplots ####

# change the name of some metabolites that give problems with glue and filenames

metab_plots = metab_clean %>% 
  mutate(Metabolite = str_replace_all(Metabolite, '/', ','))


query_metab = '2,3-phosphoglycerate'

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
    # facet_wrap(~experiment) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(here('summary', 'individual_plots', 
              glue::glue('{metabolite}_boxplot.pdf')),
         height = 7, width = 9)  
}



# stats -------------------------------------------------------------------

# WT micit vs control
metab_clean








