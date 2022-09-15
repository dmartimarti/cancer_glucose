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
  scale_y_discrete(limits=rev) +
  # panel_border() +p
  theme(strip.text.y = element_text(angle = 0))

ggsave('metabolomics/plots/barplot_pellet.pdf',
       height = 7, width = 7)

## supernatant####
bca %>% 
  filter(Sample == 'Supernatant') %>% 
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
  scale_y_discrete(limits=rev) +
  # panel_border() +p
  theme(strip.text.y = element_text(angle = 0))

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
  