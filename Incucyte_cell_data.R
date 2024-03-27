library(tidyverse)
library(readr)
library(readxl)
library(here)
library(broom)
library(cowplot)
library(ComplexHeatmap)
library(extrafont)

theme_set(theme_cowplot(15))

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
  filter(cell == 'U2-OS') %>% 
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

# Viability ---------------------------------------------------------------



data_auc %>% 
  filter(cell == 'CCD841_CoN') %>% 
  filter(!(IC == 500 & biorep == 1)) %>% 
  group_by(cell, IC, biorep) %>% 
  mutate(Viability = auc/auc[Micit_mM == 0]) 


auc_viab = data_auc %>% 
  # filter(Micit_mM == 0) %>% 
  filter(!(cell == 'CCD841_CoN' & IC == 500 & biorep == 1)) %>% 
  # filter(!(cell == 'U2-OS' & IC == 1000 & biorep == 2)) %>% 
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
  
  



stats4plot = statsR_viab_NOBIOREPS %>% 
  filter(contrast %in% c('5-1', '1-0','5-0','10-0')) %>%
  mutate(estimate = case_when(contrast == '5-1' ~ 0,
                              TRUE ~ estimate),
         p.stars = case_when(contrast == '5-1' ~ '',
                             TRUE ~ p.stars)) %>% 
  mutate(Micit_mM = case_when(contrast == '1-0' ~ 1,
                              contrast == '5-0' ~ 5,
                              contrast == '10-0' ~ 10,
                              contrast == '5-1' ~ 0)) %>% 
  mutate(Micit_mM = factor(Micit_mM, levels = c(0,1,5,10)))


auc_viab %>% 
  left_join(stats4plot) %>% 
  filter(cell == 'HCT116') %>% 
  ggplot(aes(x = Micit_mM, y =  Viability, color = Micit_mM, fill = Micit_mM)) +
  geom_point(position = position_jitterdodge()) +
  stat_summary(fun.data = "mean_cl_boot", size = 1)  +
  geom_text(aes(label = p.stars), 
            y = 1.2,
            color = 'black',
            size = 6) +
  ylim(0,1.3) +
  labs(y = 'Normalized confluence',
       x = '2-methylisocitrate (mM)',
       color = '2-methylisocitrate (mM)',
       fill = '2-methylisocitrate (mM)') +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  facet_wrap(~IC) 




cells = unique(auc_viab$cell)
for (line in cells){
  
  
  auc_viab %>% 
    left_join(stats4plot) %>% 
    filter(cell == line) %>% 
    ggplot(aes(x = Micit_mM, y =  Viability, color = Micit_mM, fill = Micit_mM)) +
    # geom_boxplot() +
    geom_point(position = position_jitterdodge()) +
    stat_summary(fun.data = "mean_cl_boot", size = 1)  +
    geom_text(aes(label = p.stars), 
              y = 1.2,
              color = 'black',
              size = 6) +
    ylim(0,1.3) +
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


# presentation plots ------------------------------------------------------


select_viab = auc_viab %>% 
  filter(cell %in% c('HCT116', 'LoVo', 'DLD-1')) %>% 
  filter((IC == 250 & cell == 'HCT116') |
           (IC == 1000 & cell == 'LoVo') |
           (IC == 1000 & cell == 'DLD-1'))

select_viab %>% 
  left_join(stats4plot %>% ungroup) %>% 
  ggplot(aes(x = Micit_mM, y =  Viability, color = Micit_mM, fill = Micit_mM)) +
  # geom_boxplot() +
  geom_point(position = position_jitterdodge(),
             show.legend = FALSE) +
  stat_summary(fun.data = "mean_cl_boot", size = 1, show.legend = F)  +
  geom_text(aes(label = p.stars), 
            y = 1.2,
            color = 'black',
            size = 6) +
  ylim(0,1.3) +
  labs(y = 'Normalized confluence',
       x = '2-methylisocitrate (mM)',
       color = '2-methylisocitrate (mM)',
       fill = '2-methylisocitrate (mM)') +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  facet_wrap(~cell) 

ggsave(here('summary/Viability_plots/presentation/', 
            'selected_viability.pdf'), 
       height = 5, 
       width = 9)



select_viab %>% 
  group_by(cell, IC, Micit_mM) %>% 
  summarise(Viab_mean = mean(Viability),
            Viab_sd = sd(Viability)) %>% 
  left_join(stats4plot %>% ungroup) %>% 
  select(cell, IC, Micit_mM, Viab_mean, Viab_sd, p.stars) %>% 
  mutate(lwr = Viab_mean - Viab_sd,
         upr = Viab_mean + Viab_sd) %>% 
  # mutate(Micit_mM = as.factor(as.character(Micit_mM))) %>%
  ungroup %>% 
  ggplot(aes(x = Micit_mM, y = Viab_mean, 
             group = cell)) +
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr,
                  fill = cell),
              alpha = 0.6) +
  geom_line(color = 'black') +
  geom_point(show.legend = F, size = 3) +
  scale_fill_viridis_d(end = 0.7) +
  labs(y = 'Cell viability ',
       x = '2-methylisocitrate (mM)',
       fill = 'Cell line') +
  facet_wrap(~cell)




select_viab %>% 
  left_join(stats4plot %>% ungroup) %>% 
  ggplot(aes(x = Micit_mM, y =  Viability, color = Micit_mM, fill = Micit_mM)) +
  geom_point(position = position_jitterdodge(),
             show.legend = FALSE, size = 1.9) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.4, 
               show.legend = F, color = 'black',
               position = position_nudge(x = 0.24))  +
  geom_text(aes(label = p.stars), 
            y = 1.2,
            color = 'black',
            size = 6) +
  ylim(0,1.3) +
  labs(y = 'Normalized confluence',
       x = '2-methylisocitrate (mM)',
       color = '2-methylisocitrate (mM)',
       fill = '2-methylisocitrate (mM)') +
  scale_fill_viridis_d() +
  scale_color_viridis_d(begin = 0.24, end = 0.94) +
  facet_wrap(~cell) 


ggsave(here('summary/Viability_plots/presentation/', 
            'selected_viability.pdf'), 
       height = 4, 
       width = 6)







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
  left_join(stats4plot) %>%
  mutate(selected = case_when(cell == 'HCT116' & IC == 500 ~ 'Yes',
                              cell == 'DLD-1' ~ 'Yes',
                              cell == 'LoVo' & IC == 1000 ~ 'Yes',
                              TRUE ~ 'No')) %>% 
  filter(selected == 'Yes') %>% 
  # filter(Micit_mM != 0) %>% 
  ggplot(aes(x = Micit_mM, y =  Viability,  fill = Micit_mM)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  stat_summary(fun.data = "mean_cl_boot", size = 1)  +
  geom_text(aes(label = p.stars), 
            y = 1.25,
            color = 'red',
            size = 7) +
  ylim(0.4,1.4) +
  labs(y = 'Normalized confluence',
       x = '2-methylisocitrate (mM)',
       color = '2-methylisocitrate (mM)',
       fill = '2-methylisocitrate (mM)') +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  facet_wrap(~cell) +
  guides(fill = 'none') +
  theme_cowplot(18) +
  theme(legend.position="top")


ggsave(here('summary/Viability_plots/', 'selected_paper.pdf'), 
       height = 7, 
       width = 12)





# heatmap ======------------------------

# remove MDA-MB-453 because it didn't grow well
cell_metadata = read_xlsx("cell_line_metadata.xlsx") %>% 
  rename(cell = `Cell line`) %>% 
  filter(cell != "MDA-MB-453")

cell_order = c("CCD841_CoN", "CCD-18Co", "HCT116", "SW1417", 
               "LoVo", "HT29", "SW837", "LS123", "MCF7",
               "RKO", "T-47D", "SW48", "T84", "SK-CO-1", 
               "U2-OS", "Hs 578T", 
               "DLD-1",  "MDA-MB-231", "SW948", "LS411N")

tissue_data = cell_metadata %>% 
  filter(cell %in% cell_order) %>% 
  select(cell, Tissue) %>% 
  column_to_rownames("cell")

tissue_data = tissue_data[cell_order,]



ha = HeatmapAnnotation(Tissue = tissue_data,
                       col = list(Tissue = c("Colon" = "#F0C73C", 
                                          "Breast" = "#EA43F0", 
                                          "Cecum" = "#3CF0CE",
                                          "Bone" = "#4B4BB5")))



# stats by biorep
auc_viab_stats = auc_viab %>% 
  filter(cell != 'HCT116_lowG') %>% 
  filter(cell %in% cell_order) %>% 
  # group_by(cell, IC, biorep, Micit_mM) %>%
  # summarise(Viability = mean(Viability)) %>%
  group_by(cell) %>% 
  nest() %>% 
  mutate( model = map(data, ~aov(Viability~Micit_mM, data = .x)) ,
          tukey = map(.x = model, ~TukeyHSD(.x)),
          tukey = map(tukey, tidy))  %>% 
  select(-model,-data) %>% 
  unnest(cols = c(tukey)) %>% 
  mutate(p.stars = gtools::stars.pval(adj.p.value)) 

auc_stats_matrix = auc_viab_stats %>% 
  filter(str_detect(contrast, "-0")) %>% 
  separate(contrast, into = c("Micit_mM", "control"), sep = "-") %>% 
  filter(cell %in% cell_order) %>% 
  select(cell, Micit_mM, p.stars) %>% 
  pivot_wider(names_from = Micit_mM, values_from = p.stars) %>% 
  mutate(`0` = "", .before = `1`) %>% 
  column_to_rownames("cell") %>% 
  as.matrix()
  
auc_stats_matrix = t(auc_stats_matrix[cell_order,])

cell_matrix = auc_viab %>% 
  filter(!(str_detect(cell, "ENEM"))) %>% 
  filter(!(str_detect(cell, "lowG"))) %>% 
  filter(cell %in% cell_order) %>% 
  group_by(cell, Micit_mM) %>% 
  summarise(mean_viab = mean(Viability)) %>% 
  ungroup %>% 
  pivot_wider(names_from = Micit_mM, values_from = mean_viab) %>% 
  column_to_rownames("cell") %>% 
  as.matrix()

cell_matrix = t(cell_matrix[cell_order,])

library(circlize)

col_fun = colorRamp2(c(0.4, 1, 1.4), c("red", "white", "green"))

pushViewport(viewport(gp = gpar(fontfamily = "Arial")))
ht = cell_matrix %>%
  Heatmap(cluster_columns = F,
          cluster_rows = F,
          col = col_fun,
          name = "Mean cell\nviability (%)", 
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%s", auc_stats_matrix[i, j]), 
                      x, y, 
                      gp = gpar(fontsize = 10))},
          top_annotation = ha,
          column_names_side = "top")
draw(ht, newpage = FALSE)

dev.copy2pdf(device = cairo_pdf,
             file = "summary/heatmap_viability_all.pdf",
             width = 10, height = 6, useDingbats = FALSE)



cell_matrix %>% t %>% 
  Heatmap(cluster_columns = T,
          cluster_rows = T,
          col = col_fun,
          name = "Mean cell\nviability (%)")
          # cell_fun = function(j, i, x, y, width, height, fill) {
          #   grid.text(sprintf("%s", auc_stats_matrix[i, j]), 
          #             x, y, 
          #             gp = gpar(fontsize = 10))},)


auc_viab %>% 
  filter(!(str_detect(cell, "ENEM"))) %>% 
  filter(!(str_detect(cell, "lowG"))) %>% 
  filter(cell %in% cell_order) %>% 
  group_by(cell, Micit_mM) %>% 
  summarise(mean_viab = mean(Viability)) %>% 
  ungroup %>% 
  filter(Micit_mM %in% c(0,10)) %>% 
  pivot_wider(names_from = Micit_mM, values_from = mean_viab) %>% 
  ggplot(aes(x = `0`, y = `10`)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(
    x = "Viability at control",
    y = "Viability at 10 mM Micit"
  ) +
  ggpubr::stat_cor(method = "pearson") 

  


# supplementary cell growth  ----------------------------------------------


# cell lines and conditions to select (initial conditions - biorep)
# CCD-18Co - 1000 (2, 3)
# CCD841_Co - 1000 (3,4)
# DLD-1 - 1000 (2)
# HCT116 - 500 (3, 4)
# Hs 578T - 500 (2)
# HT29 - 1000 (2)
# LoVo - 100- (2,3)
# LS123 - 1000 (1,2)
# LS411N - 1000 (1,2)
# MCF7 - 1000 (1)
# RKO - 250 (1)
# SK-CO-1 - 1000 (1,2)
# SW48 - 2000 (1)
# SW837 - 1000 (1)
# SW948 - 1000 (1)
# SW1417 - 1000 (1)
# T84 - 1000 (2)

cell_data = data

cell1 = cell_data %>% 
  filter(cell == "CCD-18Co", biorep %in% c(2,3), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell2 = cell_data %>% 
  filter(cell == "CCD841_CoN", biorep %in% c(4,3), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell3 = cell_data %>% 
  filter(cell == "DLD-1", biorep %in% c(2), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell4 = cell_data %>% 
  filter(cell == "HCT116", biorep %in% c(3,4), IC == 500) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell5 = cell_data %>% 
  filter(cell == "Hs 578T", biorep %in% c(2), IC == 500) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell6 = cell_data %>% 
  filter(cell == "HT29", biorep %in% c(2), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell7 = cell_data %>% 
  filter(cell == "LoVo", biorep %in% c(2,3), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell8 = cell_data %>% 
  filter(cell == "LS123", biorep %in% c(1), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell9 = cell_data %>% 
  filter(cell == "LS411N", biorep %in% c(2), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell10 = cell_data %>% 
  filter(cell == "MCF7", biorep %in% c(1), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell11 = cell_data %>% 
  filter(cell == "RKO", biorep %in% c(1), IC == 250) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell12 = cell_data %>% 
  filter(cell == "SK-CO-1", biorep %in% c(1), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell13 = cell_data %>% 
  filter(cell == "SW48", biorep %in% c(1), IC == 2000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell14 = cell_data %>% 
  filter(cell == "SW837", biorep %in% c(1), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell15 = cell_data %>% 
  filter(cell == "SW948", biorep %in% c(1), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell16 = cell_data %>% 
  filter(cell == "SW1417", biorep %in% c(1), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell17 = cell_data %>% 
  filter(cell == "T84", biorep %in% c(2), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell18 = cell_data %>% 
  filter(cell == "ZR-75-1", biorep %in% c(1), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell19 = cell_data %>% 
  filter(cell == "MDA-MB-453", biorep %in% c(1), IC == 500) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell20 = cell_data %>% 
  filter(cell == "MDA-MB-231", biorep %in% c(1), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell21 = cell_data %>% 
  filter(cell == "T-47D", biorep %in% c(1), IC == 2000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell22 = cell_data %>% 
  filter(cell == "SK-BR-3", biorep %in% c(1), IC == 500) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup

cell23 = cell_data %>% 
  filter(cell == "U2-OS", biorep %in% c(1), IC == 1000) %>% 
  group_by(cell, IC, elapsed, Micit_mM) %>% 
  summarise(Mean = mean(confluence2), SD = sd(confluence2)) %>% ungroup


cells_paper = cell1 %>% 
  bind_rows(cell2,cell3,cell4,cell5,cell6,cell7,cell8,cell9,cell10,cell11,
            cell12,cell13,cell14,cell15,cell16,cell17,cell20,
            cell21, cell23)


cell_order 

cells_paper %>% 
  filter(elapsed < 300) %>% 
  ggplot(aes(x = elapsed, y = Mean, color = Micit_mM, fill = Micit_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), 
              color = NA, alpha = 0.2)+
  geom_line() +
  labs(y = 'Mean confluence (+- SD)',
       x = 'Time (h)') +
  scale_x_continuous(breaks = seq(0, 240, by = 60),
                     limits = c(0, 240)) +
  # xlim(0, 240) +
  theme_cowplot(15,
                font_family = 'Arial') +
  # background_grid() +
  panel_border() +
  theme(legend.position = "bottom", 
        strip.text = element_text(size = 13)) +
  scale_fill_viridis_d(option = "inferno", end = 0.9, direction = -1) + 
  scale_color_viridis_d(option = "inferno", end = 0.9, direction = -1) +
  facet_wrap(~cell, nrow = 4)

ggsave('summary/growth_supplementary_figure.pdf', height = 11, width = 14)



# save data tables PAPER --------------------------------------------------



auc_viab %>% 
  filter(cell %in% cell_order) %>% 
  write_csv("summary/TABLES/auc_viability.csv")




auc_viab_stats %>% 
  filter(cell %in% cell_order) %>% 
  filter(str_detect(contrast, "-0")) %>% 
  separate(contrast, into = c("Micit_mM", "control"), sep = "-") %>% 
  filter(cell %in% cell_order) %>% ungroup %>% 
  write_csv("summary/TABLES/auc_viability_stats.csv")


cells_paper %>% 
  write_csv("summary/TABLES/growth_summary_stats.csv")
  



# max growth rate ---------------------------------------------------------



library(growthrates)
# helper function
easyfit = function(df, h){
  fit = fit_easylinear(df$elapsed, df$Mean, h = h)
  return(coef(fit))
}





cells_paper %>% 
  group_by(cell) %>% 
  slice_max(Mean) %>% arrange(elapsed)

test = cells_paper %>% 
  filter(cell == 'HCT116', Micit_mM == 0, elapsed > 10)

fit_easylinear(test$elapsed, test$Mean, h = 10)

easyfit(test, h = 10)
fit_names = names(easyfit(test, h = 10))


mumax = cells_paper %>% 
  # remove first timepoints
  filter(elapsed > 10) %>% 
  # filter(Cell == 'WT', Condition == 'Control') %>% 
  group_by(cell, Micit_mM) %>%
  nest() %>%
  # modify the value of h to take a larger or smaller time window
  summarise(numax = map(data, easyfit, 30)) %>% 
  unnest_wider(numax) %>% 
  ungroup

mumax %>% filter(Micit_mM == 0) %>% arrange(desc(mumax)) %>% 
  select(-Micit_mM) %>% write_csv("summary/TABLES/mumax_cells.csv")


mumax %>% 
  # filter(cell == 'HCT116') %>% 
  mutate(micit = as.numeric(as.vector(Micit_mM))) %>% 
  ggplot(aes(x = micit, y = mumax)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = "pearson") +
  facet_wrap(~cell, scale = 'free_y')



mumax %>% 
  filter(!(cell %in% c('CCD-18Co',
                     'CCD841_CoN'))) %>% 
  mutate(micit = as.numeric(as.vector(Micit_mM))) %>% 
  filter(Micit_mM %in% c(0, 10)) %>% 
  select(cell, Micit_mM, mumax) %>% 
  pivot_wider(values_from = mumax, names_from = Micit_mM) %>% 
  ggplot(aes(x = `0`, y = `10`)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(
    x = "Max GR at control",
    y = "Max GR at 10 mM Micit"
  ) +
  ggpubr::stat_cor(method = "pearson") 
  

mumax %>% 
  filter(!(cell %in% c('CCD-18Co',
                       'CCD841_CoN'))) %>% 
  mutate(micit = as.numeric(as.vector(Micit_mM))) %>% 
  filter(Micit_mM %in% c(0, 10)) %>% 
  select(cell, Micit_mM, mumax) %>% 
  pivot_wider(values_from = mumax, names_from = Micit_mM) %>% 
  mutate(diff = `10` - `0`) %>% 
  arrange(diff) %>% 
  ggplot(aes(x = diff, y = fct_reorder(cell, diff))) +
  geom_point(aes(fill = diff), shape = 21, size = 2) +
  labs(y = "cell line", 
       x = "Difference of mumax (treatment vs control)") %>% 
  theme_cowplot(15)



mumax %>% 
  filter(!(cell %in% c('CCD-18Co',
                       'CCD841_CoN',
                       'LS123'))) %>%
  mutate(micit = as.numeric(as.vector(Micit_mM))) %>% 
  filter(Micit_mM %in% c(0, 10)) %>% 
  select(cell, Micit_mM, mumax) %>% 
  pivot_wider(values_from = mumax, names_from = Micit_mM) %>% 
  mutate(diff = `10` - `0`) %>% 
  arrange(diff) %>% 
  ggplot(aes(x = `0`, y = diff)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(
    x = "Max GR at control",
    y = "Difference of mumax (treatment vs control)"
  ) +
  ggrepel::geom_text_repel(aes(label = cell)) +
  ggpubr::stat_cor(method = "pearson") 


