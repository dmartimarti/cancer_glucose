# libraries

library(tidyverse)
library(readxl)
library(here)
library(openxlsx)

theme_set(theme_classic())



# load data ---------------------------------------------------------------

rep1 = read_excel("rep1/Biolog 241120.xlsx", sheet = "Data") %>% 
  mutate(Micit = as.factor(Micit),
         Plate = as.factor(Plate))

drugs = read_excel("biolog_metabolites_cancer.xlsx") %>% 
  mutate(Plate = as.factor(Plate))

data = left_join(rep1, drugs) %>% 
  mutate(Well = as.factor(Well),
         Drug = as.factor(Drug))



# exploration -------------------------------------------------------------

# how are the Neg Controls behaving?
# 
# data %>% 
#   filter(Drug %in% c('Negative Control 1','Negative Control 2',
#                      'Negative Control 3','Negative Control 4')) %>% 
#   ggplot(aes(y = Value, x = Micit, fill = Micit)) +
#     geom_boxplot()+
#   geom_point(position = position_jitterdodge())


# summarise controls
controls = data %>% 
  filter(DrugU %in% c('Negative Control|1','Negative Control|2',
                     'Negative Control|3','Negative Control|4')) %>% 
  group_by(Plate, Micit) %>% 
  summarise(Mean_control = mean(Value, na.rm = TRUE),
            SD_control = sd(Value, na.rm = TRUE))

controls

# control dispersion between neg controls
controls %>% 
  ggplot(aes(x = Micit, y = Mean_control, fill = Micit)) +
  geom_histogram(stat = 'identity') +
  geom_errorbar(aes(ymin = Mean_control - SD_control, ymax = Mean_control + SD_control), width = 0.2) +
  facet_wrap(~Plate)

### only use as a control the one of micit = 0, as we want a global control for all the treatments. 
data = data %>% 
  left_join(controls %>% filter(Micit == 0) %>% select(Plate, Mean_control)) %>% 
  mutate(Viability = (Value / Mean_control) * 100,
         Drug_conc = str_sub(DrugU, -1),
         Drug_conc = as.factor(Drug_conc)) %>% 
  select(Plate, Well, Micit, Drug, Drug_conc, Value, Viability, Mean_control) 


drug_list = unique(as.character(data$Drug))



drug = 'Acriflavinium Hydrochloride'
data %>% filter(Drug == drug) %>% 
  ggplot(aes(x = Micit, y = Drug_conc, fill = Viability)) +
  geom_tile() +
  scale_fill_gradient(name = "Viability",
                      low = "#FFFFFF",
                      high = "#012345") +
  ggtitle(label = drug) +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5))

for (drug in drug_list){
  p = data %>% filter(Drug == drug) %>% 
    ggplot(aes(x = Micit, y = Drug_conc, fill = Viability)) +
    geom_tile() +
    # scale_fill_gradient(name = "Viability",
    #                     low = "#FFFFFF",
    #                     high = "#012345") +
    scale_fill_gradientn(name = 'Viability',
                         colours = c('#FFFFFF', '#012345'), limits = c(13,130)) +
    ggtitle(label = drug) +
    theme(strip.placement = "outside",
          plot.title = element_text(hjust = 0.5))
  ggsave(plot = p, here('exploration/heatmaps',filename = paste0(drug,'_heatmap.pdf') ),device = 'pdf',
         width = 12, height = 11, units = 'cm')
}


data %>% 
  ggplot(aes(x = Micit, y = Drug_conc, fill = Viability)) +
  geom_tile() +
  scale_fill_gradient(name = "Viability",
                      low = "#FFFFFF",
                      high = "#012345") +
  facet_wrap(~Drug, ncol = 9) +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5))


ggsave(here('exploration',filename = 'Complete_heatmap_rawValues.pdf') ,device = 'pdf',
       width = 40, height = 26, units = 'cm')

# Synnergy finder ---------------------------------------------------------


sfinder = data %>% 
  rename(Response = Viability,
         Drug2 = Drug,
         Conc1 = Micit,
         Conc2 = Drug_conc) %>% 
  mutate(Drug1 = 'Micit',
         ConcUnit = 'A.U.',
         PairIndex = 0) %>% 
  select(-Well, -Value, -Mean_control, -Plate) %>% 
  arrange(Drug2, Conc2) 


viab_ctr = sfinder %>% filter(Drug2 == 'Negative Control') %>%
  group_by(Conc1) %>% 
  summarise(Mean = mean(Response))



drug_list = unique(as.character(sfinder$Drug2))

sfinder_exp = expand_grid(Conc1 = c(0,1,5,10),
                          Conc2 = c(0,1,2,3,4),
                          Drug1 = 'Micit',
                          Drug2 = drug_list,
                          ConcUnit = 'A.U.')

sfinder = sfinder %>% 
  full_join(sfinder_exp %>% 
              mutate(Conc1 = as.factor(Conc1),
                     Conc2 = as.factor(Conc2))) %>%
  filter(Drug2 != 'Negative Control')

for (drug in drug_list) {
  sfinder[sfinder$Drug2 == drug & sfinder$Conc2 == 0 & sfinder$Conc1 == 0,]$Response = viab_ctr[1,2]$Mean
  sfinder[sfinder$Drug2 == drug & sfinder$Conc2 == 0 & sfinder$Conc1 == 1,]$Response = viab_ctr[2,2]$Mean
  sfinder[sfinder$Drug2 == drug & sfinder$Conc2 == 0 & sfinder$Conc1 == 5,]$Response = viab_ctr[3,2]$Mean
  sfinder[sfinder$Drug2 == drug & sfinder$Conc2 == 0 & sfinder$Conc1 == 10,]$Response = viab_ctr[4,2]$Mean

}


for (i in 1:length(drug_list)) {
  sfinder[sfinder$Drug2 == drug_list[i],]$PairIndex = i
}

sfinder = sfinder  %>% arrange(Drug2, Conc1, Conc2)

# Save statistical analysis results

list_of_datasets = list('PM-M1' = sfinder)

write.xlsx(list_of_datasets, here('exploration', 'SynergyFinder_grid.xlsx'), colNames = T, rowNames = F) 



# zip scores --------------------------------------------------------------


library(readxl)
result_ZIP = read_excel("exploration/zip_results/result_ZIP_2020-12-03.xlsx") %>% 
  rename(Drug = Drug.combination) %>% 
  mutate(Drug = str_sub(Drug, 1,-9))


data.zip = data %>% left_join(result_ZIP)  %>% 
  mutate(Syn.direction = case_when(Synergy.score < 0 ~ 'Antagonic',
                                   Synergy.score >= 0 ~ 'Synergic'))

drug = 'Acriflavinium Hydrochloride'
data.zip %>% filter(Drug == drug) %>% 
  ggplot(aes(x = Micit, y = Drug_conc, fill = Viability)) +
  geom_tile() +
  scale_fill_gradient(name = "Viability",
                      low = "#FFFFFF",
                      high = "#012345") +
  geom_text(aes(2,2, label=Synergy.score), size = 13) +
  ggtitle(label = drug) +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5))

for (drug in drug_list){
  p = data.zip %>% filter(Drug == drug) %>% 
    ggplot(aes(x = Micit, y = Drug_conc, fill = Viability)) +
    geom_tile() +
    # scale_fill_gradient(name = "Viability",
    #                     low = "#FFFFFF",
    #                     high = "#012345") +
    scale_fill_gradientn(name = 'Viability',
                         colours = c('#FFFFFF', '#012345'), limits = c(13,130)) +
    geom_text(aes(2,2, label=Synergy.score), size = 13) +
    ggtitle(label = drug) +
    theme(strip.placement = "outside",
          plot.title = element_text(hjust = 0.5))
  ggsave(plot = p, here('exploration/heatmaps_zip',filename = paste0(drug,'_heatmap.pdf') ),device = 'pdf',
         width = 12, height = 11, units = 'cm')
}


data.zip %>% 
  ggplot(aes(x = Micit, y = Drug_conc, fill = Viability)) +
  geom_tile() +
  scale_fill_gradient(name = "Viability",
                      low = "#FFFFFF",
                      high = "#012345") +
  facet_wrap(~Drug, ncol = 9) +
  geom_text(aes(2,2, label=Synergy.score, color = Syn.direction), size = 5) +
  scale_color_manual(name = 'Synergy\n direction', values = c('#00A325', '#F0011A')) +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5))


ggsave(here('exploration',filename = 'Complete_heatmap_ZIP.pdf') ,device = 'pdf',
       width = 40, height = 26, units = 'cm')



# SynergyFinder R ---------------------------------------------------------

library(synergyfinder)

sfinderR = sfinder %>% rename(
  conc_r = Conc1,
  conc_c = Conc2,
  drug_row = Drug1,
  drug_col = Drug2,
  block_id = PairIndex,
  response = Response,
  conc_r_unit = ConcUnit
) %>% mutate(conc_c_unit = 'A.U.')


dose.response.mat <- ReshapeData(sfinderR,
                                 data.type = "viability",
                                 impute = TRUE,
                                 noise = TRUE,
                                 correction = "non")

PlotDoseResponse(dose.response.mat, save.file = TRUE)



synergy.score <- CalculateSynergy(data = dose.response.mat,
                                  method = "ZIP")

str(synergy.score)

PlotSynergy(synergy.score, type = "all", save.file = TRUE)




