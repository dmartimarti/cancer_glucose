library(tidyverse)
library(readxl)
library(cowplot)


# read data ---------------------------------------------------------------


bca = read_excel("metabolomics/BCA_producersnonproducers.xlsx", 
                                        sheet = "Micit_norm")

bact_class = read_excel("metabolomics/BCA_producersnonproducers.xlsx", 
                         sheet = "list_bugs")

# plot by sample #####
## pellet ####
bca %>% 
  filter(Sample == 'Pellet') %>% 
  filter(Bacteria != 'gltA prpB double mutant') %>% 
  mutate(Producers = factor(Producers, levels = c('Producer',
                                                  'Non-producer'))) %>%
  ggplot(aes(y = fct_reorder(Bacteria, Phylum), 
             x = Micit, fill = Phylum)) +
  geom_bar(stat = 'identity') +
  facet_grid(Producers ~ ., 
             scales = "free_y",
             space = "free") +
  labs(
    x = 'Micit',
    y = 'Bacteria'
  ) +
  scale_fill_manual(values = c(
    '#02D616',
    '#2D5AD6',
    '#D9024B',
    '#DBA10D'
    )) +
  theme_half_open(13) +
  scale_y_discrete(limits=rev) 

ggsave('metabolomics/plots/barplot_pellet.pdf',
       height = 7, width = 7)



# plot by phyla -----------------------------------------------------------

# plot by sample #####
## pellet ####
bca %>% 
  # filter(Sample == 'Pellet') %>% 
  filter(Bacteria != 'gltA prpB double mutant') %>% 
  mutate(Sample = str_replace(Sample, 'Pellet', 'Endogenous'),
         Sample = str_replace(Sample, 'Supernatant', 'Excreted')) %>% 
  ggplot(aes(x = Phylum, 
             y = Micit, fill = Phylum)) +
  geom_boxplot(show.legend = F, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(),
             show.legend = F) +
  facet_wrap(~ Sample) +
  labs(
    x = NULL,
    y = 'Micit concentration'
  ) +
  scale_fill_manual(values = c(
    '#F9B91A', # actinobac
    '#D9024B', # bacteroidetes
    '#6C2AF1', # Firmicutes
    '#00B9EA'  # proteobact
  )) +
  theme_half_open() +
  theme(
    axis.text.x =  element_text(angle = 45, hjust = 1)
  )

ggsave('metabolomics/plots/boxplot_phyla.pdf',
       height = 5, width = 5.6)

## supernatant####
bca %>% 
  filter(Sample == 'Supernatant') %>% 
  filter(Bacteria != 'gltA prpB double mutant') %>% 
  mutate(Producers = factor(Producers, levels = c('Producer',
                                                  'Non-producer'))) %>% 
   ggplot(aes(y = fct_reorder(Bacteria, Phylum), 
             x = Micit, fill = Phylum)) +
  geom_bar(stat = 'identity') +
  facet_grid(Producers ~ ., 
             scales = "free_y",
             space = "free") +
  labs(
    x = 'Micit',
    y = 'Bacteria'
  ) +
  scale_fill_manual(values = c(
    '#02D616',
    '#2D5AD6',
    '#D9024B',
    '#DBA10D'
  )) +
  theme_half_open(13) +
  scale_y_discrete(limits=rev)

ggsave('metabolomics/plots/barplot_supernatant.pdf',
       height = 7, width = 7)





bca %>% 
  full_join(bact_class)


bca %>% 
  left_join(bact_class) %>% 
  filter(Sample == 'Pellet') %>% 
  drop_na(Class) %>% 
  ggplot(aes(y = Bacteria, x = Micit, 
             fill = Sample)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~Class) + 
  theme_cowplot(15) 


bca %>% 
  left_join(bact_class) %>% 
  filter(Sample == 'Supernatant') %>% 
  drop_na(Class) %>% 
  ggplot(aes(y = Bacteria, x = Micit, 
             fill = Sample)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~Class) + 
  theme_cowplot(15) 
  