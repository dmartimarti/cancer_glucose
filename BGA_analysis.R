# BGA analysis

library(tidyverse)
library(readxl)
library(ggrepel)
library(glue)


options(width = 220)

# my own library of functions
source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')


odir <- 'Plots'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

# directory 
# "/Users/dmarti14/Documents/MRC_Postdoc/Projects/Leo/BGA"


# Get timeseries data
# It will take a while
time.data = read_csv('Output/Timeseries.csv', quote = "\"") %>%
	filter(Data == '595nm_f') %>%
	select(-c(Variable1, Replicate_x)) %>%
	gather(Time_s, OD, matches('\\d')) %>%
	filter(!is.na(Media)) %>% # Remove empty values if there are missmatches
	rename(Str = Metformin_mM) %>%
	mutate(Str = as.character(Str), # Change strain namings
           Time_s = as.numeric(Time_s),
           Time_h = Time_s/3600,
           Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
           Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
           Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
           Col = factor(Col, levels = LETTERS[1:8])) %>%
 	select(c(-File, -Pattern, -Strain)) %>%
 	rename(Strain = Str, Replicate = Replicate_y) %>%
 	mutate_at(c('Well', 'Media', 'Strain'), as.factor)



tsum = time.data %>%
	group_by(Strain, Media, Time_h) %>%
	summarise(Mean = mean(OD), SD = sd(OD), SE = SD/sqrt(length(OD))) %>%
	ungroup 


## plotting all temporal series in a grid
lvls = naturalsort::naturalsort(unique(tsum$Well))

tsum %>%
  # filter(Well == 'A1') %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Media, color = Media)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Media") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  facet_wrap(vars(Strain)) +
  # scale_colour_manual(name = legend_lab, values = colours) +
  # scale_fill_manual(name = legend_lab, values = colours) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50")) 



ggsave(file = paste(odir,'/Bacterial_growth.pdf',sep = ''),
       width = 80, height = 80, units = 'mm', scale = 2, device = 'pdf')





#### second data set from Leo at BGA_10_01_2019


# Get timeseries data
# It will take a while
time.data = read_csv('BGA_10_01_2019/Output/Timeseries.csv', quote = "\"") %>%
	filter(Data == '595nm_f') %>%
	select(-c(Replicate_x)) %>%
	rename(FU = `5FU`) %>%
	gather(Time_s, OD, matches('\\d')) %>%
	filter(!is.na(Media)) %>% # Remove empty values if there are missmatches
	rename(Str = Metformin_mM) %>%
	mutate(Str = as.character(Str), # Change strain namings
           Time_s = as.numeric(Time_s),
           Time_h = Time_s/3600,
           Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
           Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
           Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
           Col = factor(Col, levels = LETTERS[1:8])) %>%
 	select(c(-File, -Pattern, -Media)) %>%
 	rename(Strain = Str, Replicate = Replicate_y) %>%
 	mutate_at(c('Well', 'Strain'), as.factor)




tsum = time.data %>%
	group_by(Strain, Background, FU, Time_h) %>%
	summarise(Mean = mean(OD), SD = sd(OD), SE = SD/sqrt(length(OD))) %>%
	ungroup 


## plotting all temporal series in a grid
lvls = naturalsort::naturalsort(unique(tsum$Well))

tsum %>%
  # filter(Well == 'A1') %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Background, color = Background)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Background") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  facet_wrap(vars(Strain, FU), ncol = 6) +
  # scale_colour_manual(name = legend_lab, values = colours) +
  # scale_fill_manual(name = legend_lab, values = colours) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50")) 


ggsave(file = paste(odir,'/Bacterial_growth.pdf',sep = ''),
       width = 80, height = 80, units = 'mm', scale = 2, device = 'pdf')






##### 
## new bacterial growth assay from 19/03/2019


setwd('/Users/dmarti14/Documents/MRC_Postdoc/Projects/Leo/BGA/19_03_19')




# Get timeseries data

time.data = read_csv('Output/Timeseries.csv', quote = "\"") %>%
  filter(Data == '595nm_f') %>%
  select(-c(Variable1, Replicate_x)) %>%
  gather(Time_s, OD, matches('\\d')) %>%
  filter(!is.na(Media)) %>% # Remove empty values if there are missmatches
  rename(Str = Metformin_mM) %>%
  mutate(Str = as.character(Str), # Change strain namings
           Time_s = as.numeric(Time_s),
           Time_h = Time_s/3600,
           Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
           Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
           Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
           Col = factor(Col, levels = LETTERS[1:8])) %>%
  select(c(-File, -Pattern, -Strain)) %>%
  rename(Strain = Str, Replicate = Replicate_y) %>%
  mutate_at(c('Well', 'Media', 'Strain'), as.factor)



tsum = time.data %>%
  group_by(Strain, Media, Time_h) %>%
  summarise(Mean = mean(OD), SD = sd(OD), SE = SD/sqrt(length(OD))) %>%
  ungroup 



## plotting all temporal series in a grid
lvls = naturalsort::naturalsort(unique(time.data$Well))

tsum %>%
  # filter(Well == 'A1') %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Media, color = Media)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Media") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  facet_wrap(vars(Strain)) +
  # scale_colour_manual(name = legend_lab, values = colours) +
  # scale_fill_manual(name = legend_lab, values = colours) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50")) 



ggsave(file = paste(odir,'/BGA_BW_pyrE_gltA_TM.pdf',sep = ''),
       width = 80, height = 80, units = 'mm', scale = 2, device = 'pdf')




##### 
## new bacterial growth assay from 29/03/2019


setwd('/Users/dmarti14/Documents/MRC_Postdoc/Projects/Leo/BGA/29_03_19')




# Get timeseries data

time.data = read_csv('Output/Timeseries.csv', quote = "\"") %>%
  filter(Data == '595nm_f') %>%
  select(-c(Variable1, Replicate_x)) %>%
  gather(Time_s, OD, matches('\\d')) %>%
  filter(!is.na(Media)) %>% # Remove empty values if there are missmatches
  rename(Str = Metformin_mM) %>%
  mutate(Str = as.character(Str), # Change strain namings
           Time_s = as.numeric(Time_s),
           Time_h = Time_s/3600,
           Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
           Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
           Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
           Col = factor(Col, levels = LETTERS[1:8])) %>%
  select(c(-File, -Pattern, -Strain)) %>%
  rename(Strain = Str, Replicate = Replicate_y) %>%
  mutate_at(c('Well', 'Media', 'Strain'), as.factor)



tsum = time.data %>%
  group_by(Strain, Media, Time_h) %>%
  summarise(Mean = mean(OD), SD = sd(OD), SE = SD/sqrt(length(OD))) %>%
  ungroup 



## plotting all temporal series in a grid
lvls = naturalsort::naturalsort(unique(time.data$Well))

tsum %>%
  # filter(Well == 'A1') %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Media, color = Media)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Media") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  facet_wrap(vars(Strain)) +
  # scale_colour_manual(name = legend_lab, values = colours) +
  # scale_fill_manual(name = legend_lab, values = colours) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50")) 



ggsave(file = paste(odir,'/BGA_BW_pyrE_gltA_TM.pdf',sep = ''),
       width = 100, height = 100, units = 'mm', scale = 2, device = 'pdf')









##### 
## new bacterial growth assay from 03/04/2019


setwd('/Users/dmarti14/Documents/MRC_Postdoc/Projects/Leo/BGA/03_04_19')




# Get timeseries data

time.data = read_csv('Output/Timeseries.csv', quote = "\"") %>%
  filter(Data == '595nm_f') %>%
  select(-c(Variable1, Replicate_x)) %>%
  gather(Time_s, OD, matches('\\d')) %>%
  filter(!is.na(Media)) %>% # Remove empty values if there are missmatches
  rename(Str = Metformin_mM) %>%
  mutate(Str = as.character(Str), # Change strain namings
           Time_s = as.numeric(Time_s),
           Time_h = Time_s/3600,
           Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
           Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
           Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
           Col = factor(Col, levels = LETTERS[1:8])) %>%
  select(c(-File, -Pattern, -Strain)) %>%
  rename(Strain = Str, Replicate = Replicate_y) %>%
  mutate_at(c('Well', 'Media', 'Strain'), as.factor)



tsum = time.data %>%
  group_by(Strain, Media, Time_h) %>%
  summarise(Mean = mean(OD), SD = sd(OD), SE = SD/sqrt(length(OD))) %>%
  ungroup 



## plotting all temporal series in a grid
lvls = naturalsort::naturalsort(unique(time.data$Well))

tsum %>%
  # filter(Well == 'A1') %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Strain, color = Strain)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Strain") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  facet_wrap(vars(Media)) +
  # scale_colour_manual(name = legend_lab, values = colours) +
  # scale_fill_manual(name = legend_lab, values = colours) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50")) 



ggsave(file = paste(odir,'/03_04_BGA_BW_pyrE_gltA_TM_strains.pdf',sep = ''),
       width = 100, height = 100, units = 'mm', scale = 2, device = 'pdf')









####################################
## new bacterial growth assay from 12/04/2019
## This time is with 4 different media, and 5 different drug concentrations
####################################



setwd("/Users/dmarti14/Documents/MRC_Postdoc/Projects/Leo/BGA/12_04_19_fourMedia")

# libraries

library(tidyverse)
library(readxl)
library(ggrepel)
library(here)


options(width = 220)

# my own library of functions
source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')


odir <- 'Plots'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



# Get timeseries data
# It will take a while
time.data = read_csv('Output/Timeseries.csv', quote = "\"") %>%
  filter(Data == '595nm_f') %>%
  gather(Time_s, OD, matches('\\d')) %>%
  select(-c(Strain)) %>%
  filter(!is.na(Media)) %>% # Remove empty values if there are missmatches
  rename(FU_uM = Metformin_mM,
         Strain = Media,
         Media = Var) %>%
  mutate( 
           Time_s = as.numeric(Time_s),
           Time_h = Time_s/3600,
           Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
           Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
           Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
           Col = factor(Col, levels = LETTERS[1:8])) %>%
  select(c(-File, -Pattern, -Replicate_x, -Reader)) %>%
  rename(Replicate = Replicate_y) %>%
  mutate_at(c('Well', 'Media', 'Strain'), as.factor)



tsum = time.data %>%
  group_by(Strain, Media, FU_uM, Time_h) %>%
  summarise(Mean = mean(OD), SD = sd(OD), SE = SD/sqrt(length(OD))) %>%
  ungroup 



# this piece of code separates the strains in three groups: 1 = (OP50, C_aq), 2 = (BW), 3 = (ndk, pyrE, TM)
tsum = tsum %>% 
  mutate(Group = ifelse(Strain %in% c('OP50p', 'C_aquatica'), 1, 
    ifelse(Strain %in% c('ndk', 'pyrE', 'TM'), 3, 2)))

# plot different drug concentrations and groups
tsum %>%
  filter(FU_uM == 50, Group >= 2) %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Strain, color = Strain)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Strain") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  facet_wrap(vars(Media)) +
  # scale_colour_manual(name = legend_lab, values = colours) +
  # scale_fill_manual(name = legend_lab, values = colours) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50")) 



ggsave(file = here('Plots', 'BGA_50um_Group2.pdf'),
       width = 80, height = 80, units = 'mm', scale = 2, device = 'pdf')





# plot different media with all drug concentrations 
tsum %>%
  filter(Media == 'Soy', Group >= 2) %>%
  mutate(FU_uM = as.factor(FU_uM)) %>%
  ggplot(aes(x = Time_h, y = Mean, fill = FU_uM, color = FU_uM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "FU_uM") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  facet_wrap(vars(Strain)) +
  # scale_colour_manual(name = legend_lab, values = colours) +
  # scale_fill_manual(name = legend_lab, values = colours) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50")) 

ggsave(file = here('Plots', 'BGA_Soy_Group2.pdf'),
       width = 80, height = 80, units = 'mm', scale = 2, device = 'pdf')










####################################
## new bacterial growth assay from 18/04/2019
## This time is with 2 different media, and 5 different drug concentrations
####################################

# working directory: /Users/dmarti14/Documents/MRC_Postdoc/Projects/Leo/BGA/04_18_19

# BGA analysis

library(tidyverse)
library(readxl)
library(ggrepel)
library(here)


options(width = 220)

# my own library of functions
source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')


odir <- 'Plots'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")





# Get timeseries data
# It will take a while
time.data = read_csv('Output/Timeseries.csv', quote = "\"") %>%
    filter(Data == '595nm_f') %>%
    gather(Time_s, OD, matches('\\d')) %>%
    select(c(-File, -Pattern, -Replicate_x, -Reader)) %>%
    filter(!is.na(Media)) %>% # Remove empty values if there are missmatches
    rename(FU_uM = Metformin_mM,
         Strain = Media,
         Media = Var) %>%
    mutate(Time_s = as.numeric(Time_s),
         Time_h = Time_s/3600,
         Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
         Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
         Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
         Col = factor(Col, levels = LETTERS[1:8])) %>%
    rename(Replicate = Replicate_y) %>%
    mutate_at(c('Well', 'Media', 'Strain'), as.factor) 


tsum = time.data %>%
    group_by(Strain, Media, FU_uM, Time_h) %>%
    summarise(Mean = mean(OD), 
              SD = sd(OD), 
              SE = SD/sqrt(length(OD)),
              ymin = Mean - SD,
              ymax = Mean + SD) %>%
    ungroup 




# plot different drug concentrations and groups
tsum %>%
    mutate(FU_uM = as.factor(FU_uM)) %>%
    ggplot(aes(x = Time_h, y = Mean, fill = Strain, color = Strain)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
    geom_line() +
    scale_x_continuous(breaks = seq(0, 24, by = 6)) +
    ylab("OD") +
    xlab("Time, h") +
    labs(fill = "Strain") +
    #facet_wrap(~MetaboliteU,ncol = 2)+
    # facet_wrap(vars(Media, FU_uM), ncol = 4) +
    # facet_grid(FU_uM ~ Media) +
    facet_grid(Media ~ FU_uM) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50")) 

quartz.save(type = 'pdf', 
    file = here('Plots', 'BGA_sidebyside2.pdf'), 
    width = 18, height = 8, family = 'Arial')



# helping functions

# function for plotting data
bga_curves = function(data, g = Strain){
    data %>%
    ggplot(aes_string(x = 'Time_h', y = 'Mean', fill = 'Strain', color = 'Strain')) +
    geom_ribbon(aes_string(ymin = 'ymin', ymax = 'ymax'), color = NA, alpha = 0.2) +
    geom_line() +
    scale_x_continuous(breaks = seq(0, 24, by = 6)) +
    ylab("OD") +
    xlab("Time, h") +
    labs(fill = "Strain") +
    #facet_wrap(~MetaboliteU,ncol = 2)+
    facet_wrap(vars(Media)) +
    # scale_colour_manual(name = legend_lab, values = colours) +
    # scale_fill_manual(name = legend_lab, values = colours) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50")) 
}

# test if the function works
tsum %>%
    filter(FU_uM == 0) %>%
    bga_curves(.)


nested = tsum %>%
    group_by(FU_uM) %>%
    nest()


plots = nested %>% 
    mutate(plot = map(.f = bga_curves, .x = data),
           filename = paste0('BGA_',as.character(FU_uM), 'uM_Bacto_Soy.pdf')) %>%
    select(filename, plot)

pwalk(plots, ggsave, path = here('Plots/'))




# plot different drug concentrations and groups
tsum %>%
    filter(FU_uM <= 16) %>%
    mutate(FU_uM = as.factor(FU_uM)) %>%
    ggplot(aes(x = Time_h, y = Mean, fill = FU_uM, color = FU_uM)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
    geom_line() +
    scale_x_continuous(breaks = seq(0, 24, by = 6)) +
    ylab("OD") +
    xlab("Time, h") +
    labs(fill = "FU_uM") +
    scale_colour_viridis_d() +
    scale_fill_viridis_d() +
    #facet_wrap(~MetaboliteU,ncol = 2)+
    # facet_wrap(vars(Media, FU_uM), ncol = 4) +
    # facet_grid(FU_uM ~ Media) +
    facet_grid(Media ~ Strain, scales = "free") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50")) 

quartz.save(type = 'pdf', 
    file = here('Plots', 'BGA_FU_uM.pdf'), 
    width = 18, height = 10, family = 'Arial')







####################################
## new bacterial growth assay from 23/05/2019
## Tesgint IPTG different concentration and prpR overexpression
####################################

# working directory: /Users/dmarti14/Documents/MRC_Postdoc/Projects/Leo/BGA/04_18_19

# BGA analysis

library(tidyverse)
library(readxl)
library(ggrepel)
library(here)


options(width = 220)

# my own library of functions
source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')


dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")





# Get timeseries data
# It will take a while
time.data = read_csv('Output/Timeseries.csv', quote = "\"") %>%
    filter(Data == '595nm_f') %>%
    gather(Time_s, OD, matches('\\d')) %>%
    select(c(-File, -Pattern, -Replicate_x, -Reader)) %>%
    filter(!is.na(Media)) %>% # Remove empty values if there are missmatches
    rename(Strain = Media,
           IPTG = var,
           Plasmid = Metformin_mM) %>%
    mutate(Time_s = as.numeric(Time_s),
           Time_h = Time_s/3600,
           Plasmid = ifelse(Plasmid == 'prpR', 'prpR', 'Control'),
           Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
           Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
           Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
           Col = factor(Col, levels = LETTERS[1:8])) %>%
    rename(Replicate = Replicate_y) %>%
    mutate_at(c('Well', 'Strain', 'Plasmid', 'IPTG'), as.factor) 


tsum = time.data %>%
    group_by(Strain, IPTG, Plasmid, Time_h) %>%
    summarise(Mean = mean(OD), 
              SD = sd(OD), 
              SE = SD/sqrt(length(OD)),
              ymin = Mean - SD,
              ymax = Mean + SD) %>%
    ungroup 



# plot different drug concentrations and groups
tsum %>%
    ggplot(aes(x = Time_h, y = Mean, fill = IPTG, color = IPTG)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, linetype = Plasmid), color = NA, alpha = 0.2) +
    geom_line(aes(linetype = Plasmid)) +
    scale_x_continuous(breaks = seq(0, 24, by = 6)) +
    ylab("OD") +
    xlab("Time, h") +
    labs(fill = "Strain") +
    facet_wrap(~Strain,ncol = 2)+
    # facet_wrap(vars(Media, FU_uM), ncol = 4) +
    # facet_grid(FU_uM ~ Media) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "grey50")) 



quartz.save(type = 'pdf', 
    file = here('Plots', 'BGA_IPTG.pdf'), 
    width = 14, height = 18, family = 'Arial')




















