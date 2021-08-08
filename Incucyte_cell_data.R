library(tidyverse)
library(readr)
library(readxl)
library(here)
library(broom)
library(cowplot)

theme_set(theme_classic())

# read data ---------------------------------------------------------------


data = read_excel("cell_summary.xlsx")

data = data %>% drop_na(confluence)

data = data %>% 
  mutate(IC = factor(IC),
         Micit_mM = as.factor(Micit_mM),
         cell = as.factor(cell))

# correct by starting everything from 0
data = data %>% 
  group_by(cell, IC, Micit_mM, replicate, biorep) %>% 
  mutate(confluence2 = confluence - min(confluence))


# this is a second run of incucyte sent by Tanara on 16/07/2021
data2 = read_excel('Run2_cells.xlsx') %>% 
  drop_na(confluence) %>% 
  mutate(IC = factor(IC),
         Micit_mM = as.factor(Micit_mM),
         cell = as.factor(cell)) %>% 
  group_by(cell, IC, Micit_mM, replicate, biorep) %>% 
  mutate(confluence2 = confluence - min(confluence))



# Growth data -------------------------------------------------------------



### summarise data ####

# summarise things
data.sum = data %>% 
  group_by(cell, IC, biorep, Micit_mM) %>% 
  # mutate(confluence = confluence - min(confluence)) %>% # remove min value from each cell
  group_by(cell, IC, elapsed, biorep, Micit_mM) %>% 
  summarise(Mean = mean(confluence2),
            SD = sd(confluence2))


data2.sum = data2 %>% 
  group_by(cell, IC, biorep, Micit_mM) %>% 
  group_by(cell, IC, elapsed, biorep, Micit_mM) %>% 
  summarise(Mean = mean(confluence2),
            SD = sd(confluence2))

data.sum %>% 
  filter(cell == 'HCT116') %>% 
  ggplot(aes(x = elapsed, y = Mean, color = Micit_mM, fill = Micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.1)+
  geom_line() +
  scale_x_continuous(breaks = seq(0, 335, by = 60)) +
  theme(legend.position = "top", strip.text = element_text(size = 13)) +
  facet_wrap(~IC*biorep*cell)

data2.sum %>% 
  filter(cell == 'SK-CO-1') %>% 
  ggplot(aes(x = elapsed, y = Mean, color = Micit_mM, fill = Micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.1)+
  geom_line() +
  scale_x_continuous(breaks = seq(0, 335, by = 60)) +
  theme(legend.position = "top", strip.text = element_text(size = 13)) +
  facet_wrap(~IC*biorep*cell)


## there is a bad rep, which one
# 
# data %>% 
#   filter(cell == 'LS1123') %>% 
#   filter(IC == 1000) %>% 
#   ggplot(aes(x = elapsed, y = confluence, color = Micit_mM)) +
#   geom_line(size =1.2) +
#   facet_wrap(~replicate*Micit_mM*biorep)
# 
# 
# 
# data.sum.filt = data %>% 
#   filter(!(replicate == 2 & IC == 1000 & Micit_mM == 0 & biorep == 1 & cell == 'LS1123')) %>% 
#   group_by(cell, IC, elapsed, biorep, Micit_mM) %>% 
#   summarise(Mean = mean(confluence),
#             SD = sd(confluence))

# data.sum.filt %>% 
#   # filter(Micit_mM == 10) %>% 
#   ggplot(aes(x = elapsed, y = Mean, color = Micit_mM, fill = Micit_mM)) +
#   geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2)+
#   geom_line() +
#   labs(title = 'LS1123',
#        y = 'Mean confluence (+- SD)',
#        x = 'Time (h)') +
#   scale_x_continuous(breaks = seq(0, 335, by = 60)) +
#   theme(legend.position = "top", strip.text = element_text(size = 13)) +
#   # facet_grid(rows = vars(biorep), cols = vars(IC)) 
#   scale_fill_viridis_d() + 
#   scale_color_viridis_d() +
#   facet_wrap(~IC*biorep*cell)

# ggsave(here('summary', 'growth_curves_bioreps.pdf'), 
#        height = 9, 
#        width = 12)



## summary stats but merging bioreps

# data.sum.filt_2 = data %>% 
#   filter(!(replicate == 2 & IC == 1000 & Micit_mM == 0 & biorep == 1 & cell == 'LS1123')) %>% 
#   group_by(cell, IC, elapsed, Micit_mM) %>% 
#   summarise(Mean = mean(confluence),
#             SD = sd(confluence))
# 
# 
# data.sum.filt_2 %>% 
#   # filter(Micit_mM == 10) %>% 
#   ggplot(aes(x = elapsed, y = Mean, color = Micit_mM, fill = Micit_mM)) +
#   geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2)+
#   geom_line() +
#   labs(title = 'LS1123',
#        y = 'Mean confluence (+- SD)',
#        x = 'Time (h)') +
#   scale_x_continuous(breaks = seq(0, 335, by = 60)) +
#   theme(legend.position = "top", strip.text = element_text(size = 13)) +
#   facet_wrap(~IC*cell)

## it's bad


## create directories for plotting


dir_name = 'summary/growth_curves_bioreps'

ifelse(!dir.exists(dir_name),dir.create(dir_name),print('Folder exists!'))



### plot individual growth curves ####

cells = unique(data$cell)
for (line in cells){

  data.sum %>% 
    filter(cell == line) %>% 
    ggplot(aes(x = elapsed, y = Mean, color = Micit_mM, fill = Micit_mM)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2)+
    geom_line() +
    labs(y = 'Mean confluence (+- SD)',
         x = 'Time (h)') +
    scale_x_continuous(breaks = seq(0, 335, by = 60)) +
    theme(legend.position = "top", strip.text = element_text(size = 13)) +
    scale_fill_viridis_d() + 
    scale_color_viridis_d() +
    facet_wrap(~IC*biorep)
  
  ggsave(here('summary/growth_curves_bioreps/', paste0(line,'growth_curves_bioreps.pdf')), 
         height = 9, 
         width = 12)
  
}

# RUN 2

dir_name = 'summary/growth_curves_bioreps_run2'

ifelse(!dir.exists(dir_name),dir.create(dir_name),print('Folder exists!'))

cells = unique(data2$cell)
for (line in cells){
  
  data2.sum %>% 
    filter(cell == line) %>% 
    ggplot(aes(x = elapsed, y = Mean, color = Micit_mM, fill = Micit_mM)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2)+
    geom_line() +
    labs(y = 'Mean confluence (+- SD)',
         x = 'Time (h)') +
    scale_x_continuous(breaks = seq(0, 335, by = 60)) +
    theme(legend.position = "top", strip.text = element_text(size = 13)) +
    scale_fill_viridis_d() + 
    scale_color_viridis_d() +
    facet_wrap(~IC*biorep)
  
  ggsave(here(dir_name, paste0(line,'growth_curves_bioreps.pdf')), 
         height = 9, 
         width = 12)
  
}



# AUC calculation ---------------------------------------------------------

library(MESS)

data_auc = data %>% 
  mutate(auc = auc(elapsed,confluence2)) %>% 
  distinct(auc,.keep_all = T) %>% 
  ungroup %>%
  select(cell:biorep,Micit_mM,auc)


# let's grab an example

cell_ex = 'HCT116'

data.sum %>% 
  filter(cell == cell_ex) %>% 
  ggplot(aes(x = elapsed, y = Mean, color = Micit_mM, fill = Micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2)+
  geom_line() +
  scale_fill_viridis_d() + 
  scale_color_viridis_d() +
  scale_x_continuous(breaks = seq(0, 335, by = 60)) +
  theme(legend.position = "top", strip.text = element_text(size = 13)) +
  facet_wrap(~IC*biorep*cell)


test = data_auc %>% filter(cell == cell_ex, IC == 250, biorep == 6)

test %>% ggplot(aes(x = Micit_mM, y = auc, fill = Micit_mM)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  scale_fill_viridis_d()


model = aov(auc ~ Micit_mM, data = test)

tidy(TukeyHSD(model))



## SECOND RUN

data2 %>% group_by(cell, IC, Micit_mM, replicate, biorep) %>% 
  summarise(N = n())

data2_auc = data2 %>% 
  mutate(auc = auc(elapsed,confluence2)) %>% 
  distinct(auc,.keep_all = T) %>% 
  ungroup %>%
  select(cell:biorep,Micit_mM,auc)


# let's grab an example

cell_ex = 'LoVo'

data2.sum %>% 
  filter(cell == cell_ex) %>% 
  ggplot(aes(x = elapsed, y = Mean, color = Micit_mM, fill = Micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2)+
  geom_line() +
  scale_fill_viridis_d() + 
  scale_color_viridis_d() +
  scale_x_continuous(breaks = seq(0, 335, by = 60)) +
  theme(legend.position = "top", strip.text = element_text(size = 13)) +
  facet_wrap(~IC*biorep*cell)


test = data_auc %>% filter(cell == cell_ex, biorep == 1)

test %>% ggplot(aes(x = Micit_mM, y = auc, fill = Micit_mM)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  scale_fill_viridis_d()


model = aov(auc ~ Micit_mM, data = test)

tidy(TukeyHSD(model))
  

# multi univariate stats --------------------------------------------------


library(openxlsx)
# library(multcomp)
# library(furrr)
# library(tictoc)

# 
# plan('multisession')
# plan('sequential')
# 
# availableCores()

statsR = data_auc %>% 
  group_by(cell,biorep,IC) %>% 
  nest() %>% 
  mutate( model = map(data, ~aov(auc~Micit_mM, data = .x)) ,
          tukey = map(.x = model, ~TukeyHSD(.x)),
          tukey = map(tukey, tidy))  %>% 
  select(-model,-data) %>% 
  unnest(cols = c(tukey)) %>% 
  mutate(p.stars = gtools::stars.pval(adj.p.value))

# the file is saved later in the script
# write.xlsx(statsR, here('summary','multi_univariate_stats.xlsx'))



statsR_run2 = data2_auc %>% 
  group_by(cell,biorep,IC) %>% 
  nest() %>% 
  mutate( model = map(data, ~aov(auc~Micit_mM, data = .x)) ,
          tukey = map(.x = model, ~TukeyHSD(.x)),
          tukey = map(tukey, tidy))  %>% 
  select(-model,-data) %>% 
  unnest(cols = c(tukey)) %>% 
  mutate(p.stars = gtools::stars.pval(adj.p.value))


statsR_global_run2 = data2_auc %>% 
  group_by(cell,IC) %>% 
  nest() %>% 
  mutate( model = map(data, ~aov(auc~Micit_mM, data = .x)) ,
          tukey = map(.x = model, ~TukeyHSD(.x)),
          tukey = map(tukey, tidy))  %>% 
  select(-model,-data) %>% 
  unnest(cols = c(tukey)) %>% 
  mutate(p.stars = gtools::stars.pval(adj.p.value))


list_of_tables = list(
  'bioreps' = statsR_run2,
  'merged_bioreps' = statsR_global_run2
)


# the file is saved later in the script
write.xlsx(list_of_tables, here('summary','multi_univariate_stats_RUN2.xlsx'))




### bioreps summary ####

# summarise things
data.sum.bioreps = data %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  # 
  mutate(elapsed = round(elapsed,0)) %>% 
  summarise(Mean = mean(confluence2),
            SD = sd(confluence2),
            SEM = SD/sqrt(n()))


data.sum.bioreps %>% 
  filter(cell == 'HCT116') %>% 
  ggplot(aes(x = elapsed, y = Mean, color = Micit_mM, fill = Micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SEM, ymax = Mean + SEM), color = NA, alpha = 0.1)+
  geom_line() +
  scale_x_continuous(breaks = seq(0, 335, by = 60)) +
  labs(y = 'Mean confluence (+- SEM)',
       x = 'Time (h)') +
  scale_fill_viridis_d() + 
  scale_color_viridis_d() +
  theme(legend.position = "top", strip.text = element_text(size = 13)) +
  facet_wrap(~IC*cell)



cells = unique(data$cell)
for (line in cells){
  
  data.sum.bioreps %>% 
    filter(cell == line) %>% 
    ggplot(aes(x = elapsed, y = Mean, color = Micit_mM, fill = Micit_mM)) +
    geom_ribbon(aes(ymin = Mean - SEM, ymax = Mean + SEM), color = NA, alpha = 0.2)+
    geom_line() +
    labs(y = 'Mean confluence (+- SEM)',
         x = 'Time (h)') +
    scale_x_continuous(breaks = seq(0, 335, by = 60)) +
    theme(legend.position = "top", strip.text = element_text(size = 13)) +
    scale_fill_viridis_d() + 
    scale_color_viridis_d() +
    facet_wrap(~IC*cell)
  
  ggsave(here('summary/growth_curves/', paste0(line,'_growth_curves.pdf')), 
         height = 9, 
         width = 12)
  
}
  


### selected cell lines ####

# HCT116 has a biorep that is slightly shorter than the others, so filter the 
# max time of it to the others and let's see

data %>% 
  filter(cell == 'HCT116', IC == 250, biorep == 2) %>% 
  select(elapsed) %>% 
  arrange(desc(elapsed))


# filter data for selected cell lines
selected = data %>% filter(cell == 'HCT116', IC == 250, elapsed < 195,
                biorep %in% c(2,3,4,6)) %>% 
  bind_rows(data  %>% filter(cell == 'HCT116', IC == 500, elapsed < 195,
                             biorep %in% c(2,3))) %>% 
  bind_rows(data  %>% filter(cell == 'LoVo')) %>% 
  bind_rows(data  %>% filter(cell == 'SW948', IC == 1000, elapsed < 200)) %>% 
  bind_rows(data  %>% filter(cell == 'SK-CO-1', elapsed < 270)) %>% 
  bind_rows(data  %>% filter(cell == 'DLD-1', IC == 1000)) %>% 
  bind_rows(data  %>% filter(cell == 'LS123', elapsed < 270))
  


# summarise things
data.sum.select = selected %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  # 
  mutate(elapsed = round(elapsed,0)) %>% 
  summarise(Mean = mean(confluence2),
            SD = sd(confluence2),
            SEM = SD/sqrt(n()))



cells = unique(selected$cell)
for (line in cells){
  
  data.sum.select %>% 
    filter(cell == line) %>% 
    ggplot(aes(x = elapsed, y = Mean, color = Micit_mM, fill = Micit_mM)) +
    geom_ribbon(aes(ymin = Mean - SEM, ymax = Mean + SEM), color = NA, alpha = 0.2)+
    geom_line() +
    labs(y = 'Mean confluence (+- SEM)',
         x = 'Time (h)') +
    scale_x_continuous(breaks = seq(0, 335, by = 60)) +
    theme(legend.position = "top", strip.text = element_text(size = 13)) +
    scale_fill_viridis_d() + 
    scale_color_viridis_d() +
    facet_wrap(~IC*cell)
  
  ggsave(here('summary/growth_curves_selected/', paste0(line,'_growth_curves.pdf')), 
         height = 9, 
         width = 12)
  
}

#### summary for the paper

data.sum.select %>% 
  mutate(selected = case_when(cell == 'HCT116' & IC == 500 ~ 'Yes',
                              cell == 'DLD-1' ~ 'Yes',
                              cell == 'LoVo' & IC == 1000 ~ 'Yes',
                              TRUE ~ 'No')) %>% 
  filter(selected == 'Yes') %>%
  ggplot(aes(x = elapsed, y = Mean, color = Micit_mM, fill = Micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SEM, ymax = Mean + SEM), color = NA, alpha = 0.2)+
  geom_line() +
  labs(y = 'Mean confluence (+- SEM)',
       x = 'Time (h)') +
  scale_x_continuous(breaks = seq(0, 335, by = 60)) +
  theme(legend.position = "top", strip.text = element_text(size = 13)) +
  scale_fill_viridis_d() + 
  scale_color_viridis_d() +
  facet_wrap(~cell, scales = 'free_y') +
  theme_cowplot(19)

ggsave(here('summary/growth_curves_selected/', 'selected_paper.pdf'), 
       height = 7, 
       width = 12)


# normalise AUC -----------------------------------------------------------


auc_viab = data_auc %>% 
  # filter(Micit_mM == 0) %>% 
  group_by(cell, IC, biorep) %>% 
  mutate(Viability = auc/auc[Micit_mM == 0]) %>% 
  ungroup
  
  
auc_viab %>% 
  filter(cell == "CCD841_CoN", IC == "500", biorep == 1)

data_auc %>% 
  filter(cell == "CCD841_CoN", IC == "500", biorep == 1)


# stats by biorep
statsR_viab = auc_viab %>% 
  group_by(cell,biorep,IC) %>% 
  nest() %>% 
  mutate( model = map(data, ~aov(Viability~Micit_mM, data = .x)) ,
          tukey = map(.x = model, ~TukeyHSD(.x)),
          tukey = map(tukey, tidy))  %>% 
  select(-model,-data) %>% 
  unnest(cols = c(tukey)) %>% 
  mutate(p.stars = gtools::stars.pval(adj.p.value))



# stats without bioreps
statsR_viab_NOBIOREPS = auc_viab %>% 
  group_by(cell,IC) %>% 
  nest() %>% 
  mutate( model = map(data, ~aov(Viability~Micit_mM, data = .x)) ,
          tukey = map(.x = model, ~TukeyHSD(.x)),
          tukey = map(tukey, tidy))  %>% 
  select(-model,-data) %>% 
  unnest(cols = c(tukey)) %>% 
  mutate(p.stars = gtools::stars.pval(adj.p.value))


list_of_tables = list(
  'stats_AUC' = statsR,
  'stats viability (bioreps)' = statsR_viab,
  'stats viability' = statsR_viab_NOBIOREPS
)

write.xlsx(list_of_tables, here('summary','multi_univariate_stats.xlsx'))



# summarise data

library(ggpubr)

auc_viab.sum = auc_viab %>% 
  group_by(cell, IC, Micit_mM) %>% 
  summarise(Mean = mean(Viability, na.rm = T),
            SD = sd(Viability, na.rm = T),
            SEM = SD / sqrt(n()))




auc_viab %>% 
  filter(cell == 'HCT116') %>% 
  ggplot(aes(x = Micit_mM, y =  Viability, color = Micit_mM, fill = Micit_mM)) +
  # geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  stat_summary(fun.data = "mean_cl_boot", size = 1)  +
  labs(y = 'Normalized confluence',
       x = '2-methylisocitrate (mM)') +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  facet_wrap(~IC*biorep) 
  
  





cells = unique(auc_viab$cell)
for (line in cells){
  
  
  auc_viab %>% 
    filter(cell == line) %>% 
    ggplot(aes(x = Micit_mM, y =  Viability, color = Micit_mM, fill = Micit_mM)) +
    # geom_boxplot() +
    geom_point(position = position_jitterdodge()) +
    stat_summary(fun.data = "mean_cl_boot", size = 1)  +
    labs(y = 'Normalized confluence',
         x = '2-methylisocitrate (mM)',
         color = '2-methylisocitrate (mM)',
         fill = '2-methylisocitrate (mM)') +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    facet_wrap(~IC) 
  
  
  ggsave(here('summary/Viability_plots/', paste0(line,'_viability.pdf')), 
         height = 9, 
         width = 12)
  
}


# plot with bioreps separately

cells = unique(auc_viab$cell)
for (line in cells){

  auc_viab %>% 
    filter(cell == line) %>% 
    ggplot(aes(x = Micit_mM, y =  Viability, color = Micit_mM, fill = Micit_mM)) +
    # geom_boxplot() +
    geom_point(position = position_jitterdodge()) +
    stat_summary(fun.data = "mean_cl_boot", size = 1)  +
    labs(y = 'Normalized confluence',
         x = '2-methylisocitrate (mM)',
         color = '2-methylisocitrate (mM)',
         fill = '2-methylisocitrate (mM)') +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    facet_wrap(~IC*biorep) 
  
  
  ggsave(here('summary/Viability_plots/bioreps/', paste0(line,'_viability.pdf')), 
         height = 9, 
         width = 12)
  
}




# viab plots paper --------------------------------------------------------


auc_viab %>% 
  mutate(selected = case_when(cell == 'HCT116' & IC == 500 ~ 'Yes',
                              cell == 'DLD-1' ~ 'Yes',
                              cell == 'LoVo' & IC == 1000 ~ 'Yes',
                              TRUE ~ 'No')) %>% 
  filter(selected == 'Yes') %>% 
  # filter(Micit_mM != 0) %>% 
  ggplot(aes(x = Micit_mM, y =  Viability,  fill = Micit_mM)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  stat_summary(fun.data = "mean_cl_boot", size = 1)  +
  labs(y = 'Normalized confluence',
       x = '2-methylisocitrate (mM)',
       color = '2-methylisocitrate (mM)',
       fill = '2-methylisocitrate (mM)') +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  facet_wrap(~cell) +
  theme_cowplot(18) +
  theme(legend.position="top")


ggsave(here('summary/Viability_plots/', 'selected_paper.pdf'), 
       height = 7, 
       width = 12)

