## Script to analyse Seahorse data from Tanara

# The experiment tested HCT116 cells (WT and -/-p53) with different
# concentrations of micit

# libraries

library(tidyverse)
library(readxl)
library(cowplot)
library(viridis)
library(broom)
library(here)
library(openxlsx)

theme_set(theme_cowplot(14) + background_grid())

# data load ---------------------------------------------------------------

sh_rate = read_excel("Seahorse_Exp1_021221/Files from machine/rates_for_R.xlsx", 
                                                 sheet = "Sheet1")


# Shaping the data 
sh_rate = sh_rate %>% 
  filter(Group != 'Background') %>% 
        # get p53 gropus
  mutate(Genotype = case_when(str_sub(Group, start = -3) == 'p53' ~ 'p53',
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
         Genotype = factor(Genotype, levels = c('WT', 'p53')),
         Group = factor(Group, levels = c('Control', 'Micit')),
         Well = as_factor(Well)) %>% 
  select(-PER) %>% 
  pivot_longer(cols = OCR:ECAR, names_to = 'Measure', values_to = 'Value') %>% 
  drop_na(Group) %>% 
  # filter out bad wells
  filter(!(Well %in% c('A02','A03','A04','A11','E07'))) 

sh_rate = sh_rate %>% 
  # name the different steps %>% 
  mutate(
    stage = case_when(Measurement %in% c(1,2,3) ~ 'stage_1',
                      Measurement %in% c(4,5,6) ~ 'stage_2',
                      Measurement %in% c(7,8,9) ~ 'stage_3',
                      Measurement %in% c(10,11,12) ~ 'stage_4')
  )


# Main plots --------------------------------------------------------------



# summarise data
sh_rate %>% 
  filter(Time < 2, Measure == 'OCR', 
         Genotype == 'WT',
         Group == 'Control')


sh_sum = sh_rate %>% 
  group_by(Measurement, Group, Micit, Time, Genotype, Measure) %>% 
  summarise(Mean = mean(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE),
            SEM = SD/sqrt(n())) %>% 
  ungroup




sh_sum %>%
  # filter(Measure == 'OCR') %>% 
  ggplot(aes(x = Time, y = Mean, color = Micit)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 1.4) +
  geom_point(size = 3) +
  geom_line(size = 1.5) +
  # facet_wrap(~Genotype*Measure) +
  facet_grid(vars(Measure),vars(Genotype)) +
  scale_color_viridis(discrete = T) + 
  labs(y = 'pmol/min (mean +/− SEM)',
       x = 'Time (in mins)')


ggsave(here('summary', 'Main_plot.pdf'),
       height = 7, width = 9)



# summary stats -------------------------------------------------------------------

# according to the Seahorse manual, there are several data that can be 
# calculated from the values obtained in the different stages

# min vals
### non-mitochondrial oxygen consumption ####
nmoc = sh_rate %>% 
  filter(Measure == 'OCR') %>% 
  filter(stage == 'stage_4') %>% 
  group_by(Group, Micit, Genotype, Replicate) %>% 
  summarise(min_val = min(Value)) %>% 
  ungroup %>% 
  mutate(description = 'Non-mitochondrial Oxygen Consumption')

### basal respiration ####
bs_resp = sh_rate %>% 
  filter(Measure == 'OCR', Measurement == 3) %>%
  left_join(nmoc %>% select(-description)) %>%
  mutate(basal_resp = Value - min_val) %>% 
  mutate(description = 'Basal Respiration')

### maximal respiration ####
mresp_temp = sh_rate %>% 
  filter(Measure == 'OCR') %>% 
  filter(stage == 'stage_3') %>% 
  group_by(Group, Micit, Genotype, Replicate) %>% 
  summarise(max_val = max(Value)) %>% 
  ungroup 


mresp = mresp_temp %>% 
  left_join(nmoc) %>% 
  mutate(max_resp = max_val - min_val) %>% 
  mutate(description = 'Maximal Respiration')


### proton leak ####

# this is how I think it should be done
# min_oligo = sh_rate %>% 
#   filter(Measure == 'OCR') %>% 
#   filter(stage == 'stage_2') %>% 
#   group_by(Group, Micit, Genotype, Replicate) %>% 
#   summarise(min_oligo = min(Value)) %>% 
#   ungroup

# this is how it's in the instructions form Agilent
min_oligo = sh_rate %>% 
  filter(Measure == 'OCR') %>% 
  filter(Measurement == 6) %>% 
  rename(min_oligo = Value) %>% 
  select(-stage,-Time,-Measure,-Measurement)

proton_leak =
  min_oligo %>% 
  left_join(nmoc) %>% 
  mutate(prot_leak = min_oligo - min_val) %>% 
  mutate(description = 'Proton Leak')

### ATP production ####

last_oligo = sh_rate %>% 
  filter(Measure == 'OCR') %>% 
  filter(Measurement == 3) %>% 
  select(-stage,-Time,-Measure,-Measurement) %>% 
  rename(last_oligo = Value)

atp_prod = min_oligo %>% 
  left_join(last_oligo) %>% 
  mutate(atp_prod = last_oligo - min_oligo) %>% 
  mutate(description = 'ATP production')



### spare respiratory capacity ####
spare_resp = mresp %>% 
  select(-description) %>% 
  left_join(bs_resp %>% select(-description)) %>% 
  mutate(spare_resp = max_resp - basal_resp)

### spare respiratory capacity as % ####
spare_resp_per = mresp %>% 
  select(-description) %>% 
  left_join(bs_resp %>% select(-description)) %>% 
  mutate(spare_resp = (max_resp/basal_resp)*100)

### coupling efficiency ####
## CHECK THIS, COULD BE WRONG ##
coup_eff = atp_prod %>% 
  select(-description,-Well) %>% 
  left_join(bs_resp %>% select(-description,-Well,
                               -Time,-min_val,-Value,
                               -stage,-Measurement)) %>% 
  mutate(coup_eff = (atp_prod/basal_resp)*100)

atp_prod %>% 
  select(-description,-Well) %>% 
  left_join(bs_resp %>% select(-description,-Well,
                               -Time,-min_val,-Value,
                               -stage,-Measurement)) %>% 
  mutate(coup_eff = (atp_prod/basal_resp)*100)




### ploting the different values ####

nmoc %>% 
  ggplot(aes(x = Micit, y = min_val, color = Micit)) +
  stat_summary(fun.data = "mean_cl_boot", size = 1.5) +
  geom_point(color = 'black') +
  facet_wrap(~Genotype) +
  scale_color_viridis(discrete = T) +
  labs(title = 'Non-mitochondrial Oxygen Consumption',
       y = 'OCR (pmol/min)',
       x = 'Micit (mM)')

ggsave(here('summary', 'non_mitochondrial_ox_consumption.pdf'),
       height = 7, width = 9)


bs_resp %>% 
  ggplot(aes(x = Micit, y = basal_resp, fill = Micit,
             color = Micit)) +
  stat_summary(fun.data = "mean_cl_boot", size = 2) +
  geom_point() +
  facet_wrap(~Genotype) +
  ylim(0,250) +
  scale_color_viridis(discrete = T) +
  labs(title = 'Basal Respiration',
       y = 'OCR (pmol/min)',
       x = 'Micit (mM)')

ggsave(here('summary', 'basal_respiration.pdf'),
       height = 7, width = 9)


  
mresp %>% 
  ggplot(aes(x = Micit, y = max_resp, fill = Micit,
             color = Micit)) +
  stat_summary(fun.data = "mean_cl_boot", size = 1.5) +
  geom_point(size = 3, alpha = .7) +
  facet_wrap(~Genotype) +
  scale_color_viridis(discrete = T) +
  labs(title = 'Maximal Respiration',
       y = 'OCR (pmol/min)',
       x = 'Micit (mM)')
  
ggsave(here('summary', 'maximal_respiration.pdf'),
       height = 7, width = 9)



proton_leak %>% 
  ggplot(aes(x = Micit, y = prot_leak, fill = Micit,
             color = Micit)) +
  stat_summary(fun.data = "mean_cl_boot", size = 1.5) +
  geom_point(size = 3, alpha = .7) +
  facet_wrap(~Genotype) +
  scale_color_viridis(discrete = T) +
  labs(title = 'Proton Leak',
       y = 'OCR (pmol/min)',
       x = 'Micit (mM)')

ggsave(here('summary', 'proton_leak.pdf'),
       height = 7, width = 9)


atp_prod %>% 
  ggplot(aes(x = Micit, y = atp_prod, fill = Micit,
             color = Micit)) +
  stat_summary(fun.data = "mean_cl_boot", size = 1.5) +
  geom_point(size = 3, alpha = .7) +
  facet_wrap(~Genotype) +
  scale_color_viridis(discrete = T) +
  labs(title = 'ATP Production',
       y = 'OCR (pmol/min)',
       x = 'Micit (mM)')

ggsave(here('summary', 'atp_production.pdf'),
       height = 7, width = 9)



spare_resp %>% 
  ggplot(aes(x = Micit, y = spare_resp, fill = Micit,
             color = Micit)) +
  stat_summary(fun.data = "mean_cl_boot", size = 1.5) +
  geom_point(size = 3, alpha = .7) +
  facet_wrap(~Genotype) +
  scale_color_viridis(discrete = T) + 
  labs(title = 'Spare Respiration',
       y = 'OCR (pmol/min)',
       x = 'Micit (mM)')
  
ggsave(here('summary', 'spare_respiration.pdf'),
       height = 7, width = 9)


# spare_resp_per %>%
#   ggplot(aes(x = Micit, y = spare_resp, fill = Micit,
#              color = Micit)) +
#   stat_summary(fun.data = "mean_cl_boot", size = 1.5) +
#   geom_point(size = 3, alpha = .7) +
#   facet_wrap(~Genotype) +
#   scale_color_viridis(discrete = T) +
#   labs(title = 'Spare respiration as %')

# 
# coup_eff %>% 
#   ggplot(aes(x = Micit, y = coup_eff, fill = Micit,
#              color = Micit)) +
#   stat_summary(fun.data = "mean_cl_boot", size = 1.5) +
#   geom_point(size = 3, alpha = .7) +
#   facet_wrap(~Genotype) +
#   scale_color_viridis(discrete = T) + 
#   labs(title = 'Coupling efficiency as %')







# t-test ------------------------------------------------------------------

library(rstatix)
library(ggpubr)

nmoc_stats = nmoc %>% 
  group_by(Genotype) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~pairwise_t_test(min_val ~ Micit, 
                                       data = ., 
                                       p.adjust.method = 'none')) 
  ) %>% 
  select(Genotype, model) %>% 
  unnest(model)


bs_resp_stats = bs_resp  %>% 
  group_by(Genotype) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~pairwise_t_test(basal_resp ~ Micit, 
                                       data = ., 
                                       p.adjust.method = 'none')) 
  ) %>% 
  select(Genotype, model) %>% 
  unnest(model)



mresp_stats = mresp %>% 
  group_by(Genotype) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~pairwise_t_test(max_resp ~ Micit, 
                                       data = ., 
                                       p.adjust.method = 'none')) 
  ) %>% 
  select(Genotype, model) %>% 
  unnest(model)


prot_leak_stats = proton_leak %>% 
  group_by(Genotype) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~pairwise_t_test(prot_leak ~ Micit, 
                                       data = ., 
                                       p.adjust.method = 'none')) 
  ) %>% 
  select(Genotype, model) %>% 
  unnest(model)


atp_prod_stats = atp_prod %>% 
  group_by(Genotype) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~pairwise_t_test(atp_prod ~ Micit, 
                                       data = ., 
                                       p.adjust.method = 'none')) 
  ) %>% 
  select(Genotype, model) %>% 
  unnest(model)

spare_resp_stats = spare_resp %>% 
  group_by(Genotype) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~pairwise_t_test(spare_resp ~ Micit, 
                                       data = ., 
                                       p.adjust.method = 'none')) 
  ) %>% 
  select(Genotype, model) %>% 
  unnest(model)

#### save stats ####

list_of_datasets = list(
  'non_mito ox consumption' = nmoc_stats,
  'basal respiration' = bs_resp_stats,
  'maximal respiration' = mresp_stats,
  'proton leak' = prot_leak_stats,
  'ATP production' = atp_prod_stats,
  'spare respiration' = spare_resp_stats
)

write.xlsx(list_of_datasets,
           here('summary','parameter_stats.xlsx'),
           overwrite = TRUE)
