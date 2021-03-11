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


### SCATTER PLOT ####

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
  labs(title = expression(paste("5FU + Glucose effect on ", italic('C. elegans'), " N2 phenotype", sep = '')),
       x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' N2 phenotype', sep = ' ')),
       y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' N2 phenotype ' , bold('(Glucose)'), sep = ' '))) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.text = element_text(size = 6)) + 
  guides(colour = guide_legend(override.aes = list(size = 4))) # make lengend points larger

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
  # labs(title = expression(paste("5FU + Glucose effect on ", italic('C. elegans'), " N2 phenotype", sep = '')),
  #      x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' N2 phenotype', sep = ' ')),
  #      y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' N2 phenotype ' , bold('(Glucose)'), sep = ' '))) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.text = element_text(size = 6)) + 
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
  # scale_color_manual(values = grad) + 
  # scale_fill_manual(values = grad) +
  # geom_text_repel(aes(label = ifelse(Genes == 'ppnP', as.character(Genes), '')), position = pos) + # this point went rogue
  geom_text_repel(aes(label = Genes), position = pos, max.overlaps = 100) +
  # labs(title = expression(paste("5FU + Glucose effect on ", italic('C. elegans'), " N2 phenotype", sep = '')),
  #      x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' N2 phenotype', sep = ' ')),
  #      y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' N2 phenotype ' , bold('(Glucose)'), sep = ' '))) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.text = element_text(size = 6)) + 
  guides(colour = guide_legend(override.aes = list(size = 4))) # make lengend points larger


new + old

ggsave(file = here('Summary', 'Scatter_sub_COMPARISON.pdf'),
       height = 10, width = 24)
