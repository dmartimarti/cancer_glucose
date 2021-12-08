## Script to analyse Seahorse data from Tanara

# The experiment tested HCT116 cells (WT and -/-p53) with different
# concentrations of micit

# libraries

library(tidyverse)
library(readxl)
library(cowplot)
library(viridis)


theme_set(theme_cowplot(14) + background_grid())

# data load ---------------------------------------------------------------

sh_rate = read_excel("Seahorse_Exp1_021221/Files from machine/rates_for_R.xlsx", 
                                                 sheet = "Sheet1")


# Shaping the data 
sh_rate = sh_rate %>% 
  filter(Group != 'Background') %>% 
        # get p53 gropus
  mutate(p53 = case_when(str_sub(Group, start = -3) == 'p53' ~ 'p53',
                    TRUE ~ 'WT'),
         # remove p53 from names
         Group = case_when(str_sub(Group, start = -3) == 'p53' ~ str_sub(Group, end = -5),
                           TRUE ~ Group),
         # get replicates
         Replicate = factor(
           str_sub(Group, start = -1)),
         
         Group = str_sub(Group, end = -3)
         # get micit concentration
  ) %>% 
  separate(Group , into = c('Group', 'Micit'),sep = ' ') %>% 
  mutate(Micit = str_sub(Micit, end = -3)) %>% 
  replace_na(list(Micit = 0)) %>% 
  mutate(Micit = factor(Micit, levels = c(0,1,5,10)),
         p53 = factor(p53, levels = c('WT', 'p53')),
         Group = factor(Group, levels = c('Control', 'Micit')),
         Well = as_factor(Well)) %>% 
  select(-PER) %>% 
  pivot_longer(cols = OCR:ECAR, names_to = 'Measure', values_to = 'Value') %>% 
  drop_na(Group)
  

# summarise data

sh_sum = sh_rate %>% 
  group_by(Measurement, Group, Micit, Time, p53, Measure) %>% 
  summarise(Mean = mean(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE)) %>% 
  ungroup


sh_sum %>% 
  filter(Measure == 'OCR') %>% 
  ggplot(aes(x = Time, y = Mean, color = Micit)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD)) +
  geom_point(size = 3) +
  geom_line(size = 1.5) +
  facet_wrap(~p53) +
  scale_color_viridis(discrete = T)













