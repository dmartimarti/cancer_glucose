## Script to analyse Seahorse data from Tanara

# The experiment tested HCT116 cells with different
# concentrations of micit in an acute exposure experiment

# libraries ####

library(tidyverse)
library(readxl)
library(cowplot)
library(viridis)
library(broom)
library(here)
library(openxlsx)
library(ggpubr)
library(rstatix)
library(ggpubr)
library(ggprism)


theme_set(theme_cowplot(14) + background_grid())



# dir.create('summary')

# data load ---------------------------------------------------------------

sh_rate = read_excel("Files from machine/rates_acute_1.xlsx", 
                     sheet = "Sheet1")


# Shaping the data 
sh_rate = sh_rate %>% 
  filter(!(Group %in% c('Background', 'Unassigned'))) %>% 
  separate(Group , into = c('Group', 'Micit'),sep = ' ') %>% 
  mutate(Micit = str_sub(Micit, end = -3)) %>% 
  replace_na(list(Micit = 0)) %>% 
  mutate(Micit = factor(Micit, levels = c(0,1,5,10)),
         # Genotype = factor(Genotype, levels = c('WT', 'p53')),
         Group = factor(Group, levels = c('Control', 'Micit')),
         Well = as_factor(Well),
         Digitonin = factor(Digitonin, 
                            levels = c(0,10,25,75))) %>% 
  select(-PER) %>% 
  pivot_longer(cols = OCR:ECAR, names_to = 'Measure', values_to = 'Value') %>% 
  drop_na(Group) 

sh_rate = sh_rate %>% 
  # name the different steps %>% 
  mutate(
    stage = case_when(Measurement %in% c(1,2,3) ~ 'stage_1',
                      Measurement %in% c(4,5,6) ~ 'stage_2',
                      Measurement %in% c(7,8,9) ~ 'stage_3',
                      Measurement %in% c(10,11,12) ~ 'stage_4')
  )



## this overly complicated and dirty piece of code works to add
# replicates to the dataset, as Tanara didn't include that info
# in the files

reps = function(df, var){
    replicates = seq(1,  dim(df)[1], 1)
    df[var] = replicates
    return(df)
  }




groups = unique(sh_rate$Group)
micit = unique(sh_rate$Micit)
digi = unique(sh_rate$Digitonin)
times = unique(sh_rate$Time)
measures = unique(sh_rate$Measure)


sh_rate2 = sh_rate

final_sh = tibble()
for (group in groups){
  for (mic in micit){
    for (dig in digi){
      for (time in times) {
        for (measure in measures) {
          target = sh_rate2 %>% 
            filter(Group == group, Micit == mic,
                   Digitonin == dig, Time == time,
                   Measure == measure)
          
          if(dim(target)[1] > 0){
            target = reps(target, 'Replicate')
            final_sh = final_sh %>% bind_rows(target)
          } else {
            next
          }
          
        }
      }
    }
  }
}


sh_rate = final_sh


# Main plots --------------------------------------------------------------



# summarise data
sh_rate %>% 
  filter(Time < 2, Measure == 'OCR', 
         Digitonin == 0,
         # Genotype == 'WT',
         Group == 'Control')


sh_sum = sh_rate %>% 
  group_by(Measurement, Group, Micit, Time, Digitonin, Measure) %>% 
  summarise(Mean = mean(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE),
            SEM = SD/sqrt(n())) %>% 
  ungroup




sh_sum %>%
  # filter(Measure == 'OCR') %>% 
  ggplot(aes(x = Time, y = Mean, color = Micit)) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                width = 1.4) +
  geom_point(size = 3) +
  geom_line(size = 1.5) +
  # facet_wrap(~Genotype*Measure) +
  facet_grid(vars(Measure),vars(Digitonin)) +
  scale_color_viridis(discrete = T) + 
  labs(y = 'pmol/min (mean +/− SEM)',
       x = 'Time (in mins)')


ggsave(here('summary', 'Main_plot.pdf'),
       height = 7, width = 9)




# let's keep only Digitonin %in% c(0, 10, 25)

sh_rate = sh_rate %>% 
  filter(Digitonin != 75)


# summary stats -------------------------------------------------------------------

# according to the Seahorse manual, there are several data that can be 
# calculated from the values obtained in the different stages

# min vals
### non-mitochondrial oxygen consumption ####
nmoc = sh_rate %>% 
  filter(Measure == 'OCR') %>% 
  filter(stage == 'stage_4') %>% 
  group_by(Group, Micit, Digitonin, Replicate) %>% 
  summarise(min_val = min(Value)) %>% 
  ungroup %>% 
  mutate(description = 'Non-mitochondrial Oxygen Consumption')

# TODO: work from here!!!!

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
  group_by(Group, Micit, Digitonin, Replicate) %>% 
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




# ploting the different values ####

nmoc %>% 
  ggplot(aes(x = Micit, y = min_val, color = Micit)) +
  stat_summary(fun.data = "mean_cl_boot", size = 1.5) +
  geom_point(color = 'black') +
  facet_wrap(~Digitonin) +
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
  facet_wrap(~Digitonin) +
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
  facet_wrap(~Digitonin) +
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

nmoc_stats = nmoc %>% 
  group_by(Genotype) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~pairwise_t_test(min_val ~ Micit, 
                                       data = ., 
                                       p.adjust.method = 'fdr')) 
  ) %>% 
  select(Genotype, model) %>% 
  unnest(model)


bs_resp_stats = bs_resp  %>% 
  group_by(Genotype) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~pairwise_t_test(basal_resp ~ Micit, 
                                       data = ., 
                                       p.adjust.method = 'fdr')) 
  ) %>% 
  select(Genotype, model) %>% 
  unnest(model)



mresp_stats = mresp %>% 
  group_by(Genotype) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~pairwise_t_test(max_resp ~ Micit, 
                                       data = ., 
                                       p.adjust.method = 'fdr')) 
  ) %>% 
  select(Genotype, model) %>% 
  unnest(model)


prot_leak_stats = proton_leak %>% 
  group_by(Genotype) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~pairwise_t_test(prot_leak ~ Micit, 
                                       data = ., 
                                       p.adjust.method = 'fdr')) 
  ) %>% 
  select(Genotype, model) %>% 
  unnest(model)


atp_prod_stats = atp_prod %>% 
  group_by(Genotype) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~pairwise_t_test(atp_prod ~ Micit, 
                                       data = ., 
                                       p.adjust.method = 'fdr')) 
  ) %>% 
  select(Genotype, model) %>% 
  unnest(model)

spare_resp_stats = spare_resp %>% 
  group_by(Genotype) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~pairwise_t_test(spare_resp ~ Micit, 
                                       data = ., 
                                       p.adjust.method = 'fdr')) 
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





# plots with pvals --------------------------------------------------------



# version with p values


###


df_p_val = bs_resp %>% 
  rstatix::group_by(Genotype) %>% 
  rstatix::pairwise_t_test(basal_resp ~ Micit, p.adjust.method = 'fdr') %>% 
  rstatix::add_y_position()


bs_resp %>% 
  ggplot(aes(x = Micit, y = basal_resp)) +
  stat_summary(aes(fill = Micit,
                   color = Micit),
               fun.data = "mean_cl_boot", size = 2) +
  geom_point(aes(fill = Micit)) +
  facet_wrap(~Genotype) +
  ylim(0,250) +
  scale_color_viridis(discrete = T) +
  labs(title = 'Basal Respiration',
       y = 'OCR (pmol/min)',
       x = 'Micit (mM)') +
  add_pvalue(df_p_val, fontface = "bold",tip.length = 0.01)


ggsave(here('summary', 'pval_plots/basal_respiration.pdf'),
       height = 7, width = 9)


###


df_p_val = nmoc %>% 
  rstatix::group_by(Genotype) %>% 
  rstatix::pairwise_t_test(min_val ~ Micit, p.adjust.method = 'fdr') %>% 
  rstatix::add_y_position()



nmoc %>% 
  ggplot(aes(x = Micit, y = min_val)) +
  stat_summary(aes(color = Micit),fun.data = "mean_cl_boot", size = 1.5) +
  geom_point(color = 'black') +
  facet_wrap(~Genotype) +
  scale_color_viridis(discrete = T) +
  labs(title = 'Non-mitochondrial Oxygen Consumption',
       y = 'OCR (pmol/min)',
       x = 'Micit (mM)') +
  add_pvalue(df_p_val, fontface = "bold", tip.length = 0.01)

ggsave(here('summary', 'pval_plots/non_mitochondrial_ox_consumption.pdf'),
       height = 7, width = 9)



###


df_p_val = mresp %>% 
  rstatix::group_by(Genotype) %>% 
  rstatix::pairwise_t_test(max_resp ~ Micit, p.adjust.method = 'fdr') %>% 
  rstatix::add_y_position()


mresp %>% 
  ggplot(aes(x = Micit, y = max_resp)) +
  stat_summary(aes( fill = Micit,
                    color = Micit), fun.data = "mean_cl_boot", size = 1.5) +
  geom_point(size = 3, alpha = .7) +
  facet_wrap(~Genotype) +
  scale_color_viridis(discrete = T) +
  labs(title = 'Maximal Respiration',
       y = 'OCR (pmol/min)',
       x = 'Micit (mM)') +
  add_pvalue(df_p_val, fontface = "bold", tip.length = 0.01)


ggsave(here('summary', 'pval_plots/maximal_respiration.pdf'),
       height = 7, width = 9)



###

df_p_val = proton_leak %>% 
  rstatix::group_by(Genotype) %>% 
  rstatix::pairwise_t_test(prot_leak ~ Micit, p.adjust.method = 'fdr') %>% 
  rstatix::add_y_position()


proton_leak %>% 
  ggplot(aes(x = Micit, y = prot_leak)) +
  stat_summary(aes(fill = Micit,
                   color = Micit),fun.data = "mean_cl_boot", size = 1.5) +
  geom_point(size = 3, alpha = .7) +
  facet_wrap(~Genotype) +
  scale_color_viridis(discrete = T) +
  labs(title = 'Proton Leak',
       y = 'OCR (pmol/min)',
       x = 'Micit (mM)') +
  add_pvalue(df_p_val, fontface = "bold", tip.length = 0.01)

ggsave(here('summary', 'pval_plots/proton_leak.pdf'),
       height = 7, width = 9)



###



df_p_val = atp_prod %>% 
  rstatix::group_by(Genotype) %>% 
  rstatix::pairwise_t_test(atp_prod ~ Micit, p.adjust.method = 'fdr') %>% 
  rstatix::add_y_position()


atp_prod %>% 
  ggplot(aes(x = Micit, y = atp_prod)) +
  stat_summary(aes(fill = Micit,
                   color = Micit), fun.data = "mean_cl_boot", size = 1.5) +
  geom_point(size = 3, alpha = .7) +
  facet_wrap(~Genotype) +
  scale_color_viridis(discrete = T) +
  labs(title = 'ATP Production',
       y = 'OCR (pmol/min)',
       x = 'Micit (mM)') +
  add_pvalue(df_p_val, fontface = "bold", tip.length = 0.01)

ggsave(here('summary', 'pval_plots/atp_production.pdf'),
       height = 7, width = 9)



###


df_p_val = spare_resp %>% 
  rstatix::group_by(Genotype) %>% 
  rstatix::pairwise_t_test(spare_resp ~ Micit, p.adjust.method = 'fdr') %>% 
  rstatix::add_y_position()


spare_resp %>% 
  ggplot(aes(x = Micit, y = spare_resp)) +
  stat_summary(aes(fill = Micit,
                   color = Micit), fun.data = "mean_cl_boot", size = 1.5) +
  geom_point(size = 3, alpha = .7) +
  facet_wrap(~Genotype) +
  scale_color_viridis(discrete = T) + 
  labs(title = 'Spare Respiration',
       y = 'OCR (pmol/min)',
       x = 'Micit (mM)')+
  add_pvalue(df_p_val, fontface = "bold", tip.length = 0.01, bracket.size = 0.9)

ggsave(here('summary', 'pval_plots/spare_respiration.pdf'),
       height = 7, width = 9)





# ref pval version --------------------------------------------------------


###

df_p_val = spare_resp %>% 
  rstatix::group_by(Genotype) %>% 
  rstatix::pairwise_t_test(spare_resp ~ Micit, p.adjust.method = 'fdr', ref.group = "0") %>% 
  rstatix::add_y_position()


spare_resp %>% 
  ggplot(aes(x = Micit, y = spare_resp)) +
  stat_summary(aes(fill = Micit,
                   color = Micit), fun.data = "mean_cl_boot", size = 1.5) +
  geom_point(size = 3, alpha = .7) +
  facet_wrap(~Genotype) +
  scale_color_viridis(discrete = T) + 
  labs(title = 'Spare Respiration',
       y = 'OCR (pmol/min)',
       x = 'Micit (mM)')+
  add_pvalue(df_p_val, label = "p.adj.signif", remove.bracket = TRUE,
             y.position = 310, label.size = 5) +
  theme_prism()

ggsave(here('summary', 'pval_plots/spare_respiration_simple.pdf'),
       height = 7, width = 9)


###


df_p_val = mresp %>% 
  rstatix::group_by(Genotype) %>% 
  rstatix::pairwise_t_test(max_resp ~ Micit, p.adjust.method = 'fdr', ref.group = "0") %>% 
  rstatix::add_y_position()


mresp %>% 
  ggplot(aes(x = Micit, y = max_resp)) +
  stat_summary(aes( fill = Micit,
                    color = Micit), fun.data = "mean_cl_boot", size = 1.5) +
  geom_point(size = 3, alpha = .7) +
  facet_wrap(~Genotype) +
  scale_color_viridis(discrete = T) +
  labs(title = 'Maximal Respiration',
       y = 'OCR (pmol/min)',
       x = 'Micit (mM)') +
  add_pvalue(df_p_val, label = "p.adj.signif", remove.bracket = TRUE,
             y.position = 310, label.size = 5) +
  theme_prism()



ggsave(here('summary', 'pval_plots/maximal_respiration_simple.pdf'),
       height = 7, width = 9)




# interaction terms ####


###
bs_resp_int = bs_resp %>% 
  rstatix::anova_test(basal_resp ~ Genotype * Micit) 

bs_resp_micit = bs_resp %>% 
  rstatix::group_by(Micit) %>% 
  rstatix::anova_test(basal_resp ~ Genotype) 


###
nmoc_int = nmoc %>% 
  rstatix::anova_test(min_val ~ Genotype * Micit)

nmoc_micit = nmoc %>% 
  rstatix::group_by(Micit) %>% 
  rstatix::anova_test(min_val ~ Genotype)


###
mresp_int = mresp %>% 
  rstatix::anova_test(max_resp ~ Genotype * Micit)

mresp_micit = mresp %>% 
  rstatix::group_by(Micit) %>% 
  rstatix::anova_test(max_resp ~ Genotype)

###

proton_leak_int = proton_leak %>% 
  rstatix::anova_test(prot_leak ~ Genotype * Micit)

proton_leak_micit = proton_leak %>% 
  rstatix::group_by(Micit) %>% 
  rstatix::anova_test(prot_leak ~ Genotype)


###
atp_prod_int = atp_prod %>% 
  rstatix::anova_test(atp_prod ~ Genotype * Micit)

atp_prod_micit = atp_prod %>% 
  rstatix::group_by(Micit) %>% 
  rstatix::anova_test(atp_prod ~ Genotype)

###


spare_resp_int = spare_resp %>% 
  rstatix::anova_test(spare_resp ~ Genotype * Micit)


spare_resp_micit = spare_resp %>% 
  rstatix::group_by(Micit) %>% 
  rstatix::anova_test(spare_resp ~ Genotype)


# save interaction terms

list_of_datasets = list(
  'non_mito ox consumption' = nmoc_int,
  'basal respiration' = bs_resp_int,
  'maximal respiration' = mresp_int,
  'proton leak' = proton_leak_int,
  'ATP production' = atp_prod_int,
  'spare respiration' = spare_resp_int
)

write.xlsx(list_of_datasets,
           here('summary','parameter_interaction_terms.xlsx'),
           overwrite = TRUE)


# save comparisons of genotype by micit



list_of_datasets = list(
  'non_mito ox consumption' = nmoc_micit,
  'basal respiration' = bs_resp_micit,
  'maximal respiration' = mresp_micit,
  'proton leak' = proton_leak_micit,
  'ATP production' = atp_prod_micit,
  'spare respiration' = spare_resp_micit
)

write.xlsx(list_of_datasets,
           here('summary','parameter_genotype.xlsx'),
           overwrite = TRUE)

