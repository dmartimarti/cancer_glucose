# libraries

library(tidyverse)
library(readxl)
library(here)
library(openxlsx)
library(cowplot)

theme_set(theme_cowplot(14))



# load data ---------------------------------------------------------------

rep1 = read_excel("rep1/Biolog 241120.xlsx", sheet = "Data") %>% 
  mutate(Micit = as.factor(Micit),
         Plate = as.factor(Plate))


rep2 = read_excel("rep2/Biolog 170721.xlsx", sheet = "Data") %>% 
  mutate(Micit = as.factor(Micit),
         Plate = as.factor(Plate))

rep3 = read_excel("rep3/Biolog_220721.xlsx", sheet = "Data") %>% 
  mutate(Micit = as.factor(Micit),
         Plate = as.factor(Plate))


drugs = read_excel("biolog_metabolites_cancer.xlsx") %>% 
  mutate(Plate = as.factor(Plate))


# load cancer drugs db from: https://www.anticancerfund.org/en/cancerdrugs-db
cancerdrugsdb = read_delim("cancerdrugsdb.txt", 
           delim = "\t", escape_double = FALSE, 
           trim_ws = TRUE)

# read file with drugbank compounds with FP 
drugbank_fp = read_csv("DrugBank_drugs_FP.csv")




rep1 = left_join(rep1, drugs) %>% 
  mutate(Well = as.factor(Well),
         Drug = as.factor(Drug), 
         Replicate = 1, .before = Micit)

rep2 = left_join(rep2, drugs) %>% 
  mutate(Well = as.factor(Well),
         Drug = as.factor(Drug),
         Replicate = 2, .before = Micit)

rep3 = left_join(rep3, drugs) %>% 
  mutate(Well = as.factor(Well),
         Drug = as.factor(Drug),
         Replicate = 3, .before = Micit)


# rep3 fix
# change Micit 1 to 0
# rep3_micit_0 = rep3 %>% filter(Micit == 1) %>% 
#   mutate(Micit = 0, Micit = factor(Micit, levels = c(0,1,5,10)))
# # change Micit 0 to 1
# rep3_micit_1 = rep3 %>% filter(Micit == 0) %>% 
#   mutate(Micit = 1, Micit = factor(Micit, levels = c(0,1,5,10)))
# 
# rep3 = rep3 %>% 
#   filter(Micit %in% c(5,10)) %>% 
#   bind_rows(rep3_micit_0, rep3_micit_1)


# exploration -------------------------------------------------------------

data = bind_rows(rep1,rep2,rep3)


# how are the Neg Controls behaving?
data %>%
  filter(DrugU %in% c('Negative Control|1','Negative Control|2',
                     'Negative Control|3','Negative Control|4')) %>%
  ggplot(aes(y = Value, x = Micit, fill = Micit)) +
    geom_boxplot()+
  geom_point(position = position_jitterdodge())

ggsave(here('exploration', 'controls_boxplots.pdf'), height = 10, width = 11)

# summarise controls
controls = data %>% 
  filter(DrugU %in% c('Negative Control|1','Negative Control|2',
                     'Negative Control|3','Negative Control|4')) %>% 
  group_by(Plate, Replicate, Micit) %>% 
  summarise(Mean_control = mean(Value, na.rm = TRUE),
            SD_control = sd(Value, na.rm = TRUE)) %>% 
  ungroup

controls

# control dispersion between neg controls
controls %>% 
  ggplot(aes(x = Micit, y = Mean_control, fill = Micit)) +
  geom_histogram(stat = 'identity') +
  geom_point() +
  geom_errorbar(aes(ymin = Mean_control - SD_control, ymax = Mean_control + SD_control), width = 0.2) +
  facet_grid(vars(Plate), vars(Replicate))

ggsave(here('exploration', 'controls_barplots_replicates.pdf'), height = 10, width = 11)



data %>% 
  filter(DrugU %in% c('Negative Control|1','Negative Control|2',
                      'Negative Control|3','Negative Control|4')) %>% 
  group_by(Plate, Micit) %>% 
  summarise(Mean_control = mean(Value, na.rm = TRUE),
            SD_control = sd(Value, na.rm = TRUE)) %>% 
  ggplot(aes(x = Micit, y = Mean_control, fill = Micit)) +
  geom_histogram(stat = 'identity') +
  geom_point() +
  geom_errorbar(aes(ymin = Mean_control - SD_control, ymax = Mean_control + SD_control), width = 0.2) +
  facet_wrap(vars(Plate))

ggsave(here('exploration', 'controls_barplots.pdf'), height = 10, width = 11)



#### calculate viability ####
### only use as a control the one of micit = 0, as we want a global control for all the treatments. 
data = data %>% 
  left_join(controls %>% filter(Micit == 0) %>% select(Plate, Replicate, Mean_control)) %>% 
  mutate(Viability = (Value / Mean_control) * 100,
         Drug_conc = str_sub(DrugU, -1),
         Drug_conc = as.factor(Drug_conc)) %>% 
  select(Plate, Well, Replicate, Micit, Drug, Drug_conc, Value, Viability, Mean_control) 


drug_list = unique(as.character(data$Drug))



#### data summary ####

data_sum = data %>%  
  group_by(Drug, Drug_conc, Micit) %>% 
  summarise(Mean = mean(Viability),
            SD = sd(Viability))


drug_list = unique(as.character(data$Drug))


# plot mean of drug viability
drug = 'Fluorouracil'
data_sum %>% 
  filter(Drug == drug) %>% 
  ggplot(aes(x = Micit, y = Drug_conc, fill = Mean)) +
  geom_tile() +
  scale_fill_gradient(name = "Viability",
                      low = "#FFFFFF",
                      high = "#F57A20") +
  ggtitle(label = drug) +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5))





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
  arrange(Conc1, Drug2, Conc2) 


viab_ctr = sfinder %>% filter(Drug2 == 'Negative Control') %>%
  group_by(Conc1, Replicate) %>% 
  summarise(Mean = mean(Response))



drug_list = unique(as.character(sfinder$Drug2))

# sfinder_exp = expand_grid(Conc1 = c(0,1,5,10),
#                           Conc2 = c(0,1,2,3,4),
#                           Drug1 = 'Micit',
#                           Drug2 = drug_list,
#                           Replicate = c(1,2,3),
#                           ConcUnit = 'A.U.')

# put the controls in an expand grid with all drugs and concs
sfinder_ctrl_exp = expand_grid(Conc1 = c(0,1,5,10),
                          Conc2 = c(0),
                          Drug1 = 'Micit',
                          Drug2 = drug_list,
                          Replicate = c(1,2,3),
                          ConcUnit = 'A.U.')

sfinder_ctrl_exp = sfinder_ctrl_exp %>% 
  mutate(Conc1 = as.factor(Conc1),
         Conc2 = as.factor(Conc2)) %>% 
  left_join(viab_ctr %>% 
              rename(Response = Mean)) 

sfinder = sfinder %>% 
  bind_rows(sfinder_ctrl_exp) %>%
  filter(Drug2 != 'Negative Control')
  
  
# 
# 
# for (drug in drug_list) {
#   sfinder[sfinder$Drug2 == 'Negative Control' & sfinder$Conc2 == 0 & sfinder$Conc1 == 0,]$Response = viab_ctr[1,3]$Mean
#   sfinder[sfinder$Drug2 == 'Negative Control' & sfinder$Conc2 == 0 & sfinder$Conc1 == 1,]$Response = viab_ctr[2,3]$Mean
#   sfinder[sfinder$Drug2 == 'Negative Control' & sfinder$Conc2 == 0 & sfinder$Conc1 == 5,]$Response = viab_ctr[3,3]$Mean
#   sfinder[sfinder$Drug2 == 'Negative Control' & sfinder$Conc2 == 0 & sfinder$Conc1 == 10,]$Response = viab_ctr[4,3]$Mean
# 
# }


for (i in 1:length(drug_list)) {
  sfinder[sfinder$Drug2 == drug_list[i],]$PairIndex = i
}

sfinder = sfinder  %>% arrange(Drug2, Conc1, Conc2) %>% 
  select(-Replicate)

# Save statistical analysis results

list_of_datasets = list('PM-M1'= sfinder)

write.xlsx(list_of_datasets, here('exploration', 'SynergyFinder_grid_3reps_No_NC.xlsx'), 
           colNames = T, rowNames = F, overwrite = TRUE) 







# test how it looks in a heatmap


sfinder_sum = sfinder %>%  
  group_by(Drug2, Drug1, Conc2, Conc1) %>% 
  summarise(Mean = mean(Response),
            SD = sd(Response)) %>% 
  ungroup %>% 
  mutate(Conc2 = factor(Conc2, levels = c(0,1,2,3,4)))


drug_list = unique(as.character(sfinder_sum$Drug2))


# plot mean of drug viability
drug = 'Rapamycin'
sfinder_sum %>% 
  filter(Drug2 == drug) %>% 
  ggplot(aes(x = Conc1, y = Conc2, fill = Mean)) +
  geom_tile() +
  scale_fill_gradient(name = "Viability",
                      low = "#FFFFFF",
                      high = "#F57A20") +
  ggtitle(label = drug) +
  labs(x = 'Micit (mM)',
       y = 'Query drug (A.U.)') +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5))


  
  
# plot a heatmap per drug
# this represents the MEAN Viability
for (drug in drug_list){
    p = sfinder_sum %>% 
      filter(Drug2 == drug) %>% 
      ggplot(aes(x = Conc1, y = Conc2, fill = Mean)) +
      geom_tile() +
      scale_fill_gradient(name = "Viability",
                          low = "#FFFFFF",
                          high = "#F57A20") +
      ggtitle(label = drug) +
      labs(x = 'Micit (mM)',
           y = 'Query drug (A.U.)') +
    theme(strip.placement = "outside",
          plot.title = element_text(hjust = 0.5))
    
    ggsave(plot = p, here('exploration/heatmaps',paste0(drug,'_heatmap.pdf') ),device = 'pdf',
           width = 12, height = 11, units = 'cm')
  }
  
  
  sfinder_sum %>% 
    ggplot(aes(x = Conc1, y = Conc2, fill = Mean)) +
    geom_tile()  +
    scale_fill_gradient(name = "Viability",
                        low = "#FFFFFF",
                        high = "#F57A20",
                        limits = c(10,150)) +
    labs(x = 'Micit (mM)',
         y = 'Query drug (A.U.)') +
    facet_wrap(~Drug2, ncol = 9) +
    theme(strip.placement = "outside",
          plot.title = element_text(hjust = 0.5))
  
  
  ggsave(here('exploration','Complete_heatmap_MeanViability.pdf') ,device = 'pdf',
         width = 40, height = 26, units = 'cm')
  
  
  
  
# zip scores --------------------------------------------------------------


library(readxl)
# 
# result_ZIP = read_excel("exploration/zip_results_rep1/result_ZIP_2020-12-03.xlsx") %>% 
#   rename(Drug = Drug.combination) %>% 
#   mutate(Drug = str_sub(Drug, 1,-9))

result_ZIP = read_excel("exploration/zip_results/result_ZIP_2021-07-30.xlsx") %>% 
  rename(Drug = Drug.combination) %>% 
  mutate(Drug = str_sub(Drug, 1,-9))
  
# THE GOOD ONE
# result_ZIP = read_excel("exploration/zip_results/result_ZIP_2021-07-30.xlsx") %>% 
#   rename(Drug = Drug.combination) %>% 
#   mutate(Drug = str_sub(Drug, 1,-9))


data.zip = data %>% left_join(result_ZIP)  %>% 
  mutate(Syn.direction = case_when(Synergy.score < 0 ~ 'Antagonic',
                                   Synergy.score >= 0 ~ 'Synergic'))

drug = 'Fluorouracil'
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



# nucleotide stats --------------------------------------------------------


result_ZIP = read_excel("exploration/zip_results/result_ZIP_2021-07-30.xlsx") %>% 
  rename(Drug = Drug.combination) %>% 
  mutate(Drug = str_sub(Drug, 1,-9))



nucl = c('Azathioprine', 'Fluorouracil',
         'Mercaptopurine', 'Thioguanine',
         "5-Fluoro-5'- DeoxyurIdine", 'Zidovudine',
         'Azacytidine', 'Carmofur',
         'Floxuridine', "Cytosine-Beta-DArabinofuranoside",
         'Methotrexate')

res_ZIP_cats = result_ZIP %>% 
  mutate(Category = case_when(Drug %in% nucl ~ 'Nucleotides',
                              TRUE ~ 'Other'))


res_ZIP_cats %>% filter(Category == 'Nucleotides')


#### stats ####
library(rstatix)

# homoscedasticity test
CG = subset(res_ZIP_cats, Category == "Nucleotides")$Most.synergistic.area.score
TG = subset(res_ZIP_cats, Category != "Nucleotides")$Most.synergistic.area.score
tidy(var.test(CG, TG))

# as they suffer from heterocedasticity, var.equal to FALSE
nucl_stats = t.test(CG, TG, var.equal = T)
library(broom)
pval = tidy(nucl_stats)$p.value

res_ZIP.sum = res_ZIP_cats %>%
  group_by(Category) %>%
  summarise(Mean = mean(Most.synergistic.area.score),
            SD = sd(Most.synergistic.area.score),
            SEM = SD/sqrt(n()))

## MAIN PLOT
res_ZIP_cats %>%
  ggplot(aes(x = Category, y = Most.synergistic.area.score, fill = Category)) +
  geom_violin() +
  geom_jitter(size = 0.5,height = 0, width = 0.1) +
  geom_point(data = res_ZIP.sum,
             size = 3.5,
             color = 'black',
             aes(x = Category, y = Mean)) +
  geom_errorbar(data = res_ZIP.sum,
                aes(x = Category, y = Mean,
                    ymin = Mean - SEM, ymax = Mean + SEM),
                width = 0.03) +
  geom_segment(aes(x = 0, y = res_ZIP.sum$Mean[1], xend = 1, yend = res_ZIP.sum$Mean[1]),
               linetype="dashed", colour = 'grey50') +
  geom_segment(aes(x = 0, y = res_ZIP.sum$Mean[2], xend = 2, yend = res_ZIP.sum$Mean[2]),
               linetype="dashed", colour = 'grey50') +
  # pval annotation
  geom_segment(aes(x = 1, xend = 2, y = 16.5, yend = 16.5)) +
  annotate("text", x = 1.5, y = 18, label = paste('P-value: ',round(pval, 5))) +
  # annotate("text", x = 1.5, y = 18, label = paste('P-value < ','0.0001')) +
  scale_fill_manual(values = c('#F5EC49',
                                '#3D9CE6'),
                    labels = c('DNA Damage Drugs',
                               'Other')) + 
  labs(x = 'Drug category',
       y = 'Most synergistic ZIP score') +
  scale_x_discrete(name = '',
                   labels = c("Nucleotides" = "DNA Damage Drugs","Other" = "Other")) +
  theme(axis.text.x = element_text(face = "bold", size = 13, color = 'black'),
        axis.text.y = element_text(face = "bold", size = 13, color = 'black'))

ggsave(here('exploration', 'violin_nucleotides_biolog_MostSynergyScore.pdf'), height = 8, width = 9)


res_ZIP.sum %>%
  ggplot(aes(x = Category, y = Mean)) +
  geom_point(size = 5, aes(color = Category)) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.1)



#### stats ####
library(rstatix)

# homoscedasticity test
CG = subset(res_ZIP_cats, Category == "Nucleotides")$Synergy.score
TG = subset(res_ZIP_cats, Category != "Nucleotides")$Synergy.score
tidy(var.test(CG, TG))

# as they suffer from heterocedasticity, var.equal to FALSE
nucl_stats = t.test(CG, TG, var.equal = T)
library(broom)
pval = tidy(nucl_stats)$p.value

res_ZIP.sum = res_ZIP_cats %>%
  group_by(Category) %>%
  summarise(Mean = mean(Synergy.score),
            SD = sd(Synergy.score),
            SEM = SD/sqrt(n()))

## MAIN PLOT
res_ZIP_cats %>%
  ggplot(aes(x = Category, y = Synergy.score, fill = Category)) +
  geom_violin() +
  geom_jitter(size = 0.5,height = 0, width = 0.1) +
  geom_point(data = res_ZIP.sum,
             size = 3.5,
             color = 'black',
             aes(x = Category, y = Mean)) +
  geom_errorbar(data = res_ZIP.sum,
                aes(x = Category, y = Mean,
                    ymin = Mean - SEM, ymax = Mean + SEM),
                width = 0.03) +
  geom_segment(aes(x = 0, y = res_ZIP.sum$Mean[1], xend = 1, yend = res_ZIP.sum$Mean[1]),
               linetype="dashed", colour = 'grey50') +
  geom_segment(aes(x = 0, y = res_ZIP.sum$Mean[2], xend = 2, yend = res_ZIP.sum$Mean[2]),
               linetype="dashed", colour = 'grey50') +
  # pval annotation
  geom_segment(aes(x = 1, xend = 2, y = 16.5, yend = 16.5)) +
  annotate("text", x = 1.5, y = 18, label = paste('P-value: ',round(pval, 5))) +
  # annotate("text", x = 1.5, y = 18, label = paste('P-value < ','0.0001')) +
  scale_fill_manual(values = c('#F5EC49',
                               '#3D9CE6'),
                    labels = c('DNA Damage Drugs',
                               'Other')) + 
  labs(x = 'Drug category',
       y = 'ZIP score') +
  scale_x_discrete(name = '',
                   labels = c("Nucleotides" = "DNA Damage Drugs","Other" = "Other")) +
  theme(axis.text.x = element_text(face = "bold", size = 13, color = 'black'),
        axis.text.y = element_text(face = "bold", size = 13, color = 'black'))

ggsave(here('exploration', 'violin_nucleotides_biolog.pdf'), height = 8, width = 9)




# Weighted means and standard deviations
# see my paper from 2017 for the equations to calculate these values
weigted_stuff = res_ZIP_cats %>% 
  group_by(Drug) %>% 
  mutate(SD = ( sqrt(3) * ((abs(Synergy.score) + `95% CI`) - (abs(Synergy.score) - `95% CI`))/3.92),
         Var = SD**2) %>% 
  group_by(Category) %>% 
  mutate(w = (1/Var)/(sum(1/Var)),
         Wmean = sum(Synergy.score * w),
         Wsd = sqrt((1/(1-sum(w**2))) * sum(w * ((Synergy.score - Wmean)**2)) )) 

weigted_stuff.sum = weigted_stuff %>% distinct(Category, .keep_all = T) %>% 
  select(Wmean, Wsd, Category)

res_ZIP_cats %>% 
  group_by(Category) %>% 
  count()



x = weigted_stuff %>% filter(Category == 'Other') %>% pull(Synergy.score)
y = weigted_stuff %>% filter(Category == 'Nucleotides') %>% pull(Synergy.score)

wx = weigted_stuff %>% filter(Category == 'Other') %>% pull(w)
wy = weigted_stuff %>% filter(Category == 'Nucleotides') %>% pull(w)

weights::wtd.t.test(x=x, y=y, weight=wx, weighty=wy, samedata=F)



# Heatmap for Tanara ------------------------------------------------------





drug = 'Fluorouracil'
data %>% filter(Drug == drug) %>% 
  ggplot(aes(x = Micit, y = Drug_conc, fill = Viability)) +
  geom_tile() +
  scale_fill_gradient(name = "Viability",
                      low = "#FFFFFF",
                      high = "#FA7235") +
  # ggtitle(label = drug) +
  labs(x = 'Metabolite (mM)',
       y = '5-FU (A.U.)') +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(face = "bold", size = 13, color = 'black'),
        axis.title.y = element_text(face = "bold", size = 13, color = 'black'))

ggsave(here('exploration', 'micit_heatmap_conference.pdf'), height = 8, width = 9)





# Fingerprints ---------------------------------------------------

dir.create(here('exploration', 'Fingerprints_results'))


### TEST ZONE ####



library(broom)
library(plotly)
library(cluster)
library(factoextra)

## load the Morgan fingerprints produced by the script: smiles2morgan.ipynb 

morgan = read_csv("drug_cancer_biolog_morganFP.csv")


### PCA ####

pca_fit = morgan %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(scale = F) # do PCA on scaled data


pca_fit %>%
  augment(morgan) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, .fittedPC2)) + 
  geom_point(size = 1.5) +
  geom_label(aes(label = Drug), alpha = 0.6) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()

ggsave(here('exploration/Fingerprints_results', 'PCA_fgroups.pdf'), 
       height = 13, width = 14)


pca_fit %>%
  augment(morgan) %>% 
  plot_ly(x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3,
          text = ~Drug) %>% 
  add_markers()


pca_fit %>%
  tidy(matrix = "eigenvalues")

# barplot of percent explained by PC
pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

ggsave(here('exploration/Fingerprints_results', 'percent_explained.pdf'), 
       height = 9, width = 10)

# cumulative expl by PC 
pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, cumulative)) +
  geom_line(color = "#56B4E9", alpha = 0.8) +
  geom_point(color = '#56B4E9', alpha = 0.6) +
  geom_text(aes(label = round(PC,2)), nudge_x = 0, nudge_y = 0.05) +
  # scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

ggsave(here('exploration/Fingerprints_results', 'percent_explained_cumulative.pdf'), 
       height = 9, width = 10)

### distances ####

drug_dists = morgan %>% select(where(is.numeric)) %>%  data.frame %>% get_dist
fviz_dist(drug_dists, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

ggsave(here('exploration/Fingerprints_results', 'distances.pdf'), 
       height = 9, width = 10)

### K-means with morgan FP ####

morgan_mat = morgan %>% 
  select(where(is.numeric)) 



kclust = morgan_mat %>% kmeans(centers = 30)

pca_fit %>%
  augment(morgan) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  ggplot(aes(.fittedPC1, .fittedPC2, color = cluster)) + 
  geom_point(size = 3.5) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()


pca_fit %>%
  augment(morgan) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  plot_ly(x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3,
          text = ~Drug, color=~cluster) %>% 
  add_markers()


pca_fit %>%
  augment(morgan) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster)) %>% 
  select(cluster, Drug) %>% 
  left_join(drugs) %>% 
  distinct(Drug, .keep_all=T) %>% 
  filter(Drug %in% nucl) %>% view


kclusts <- 
  tibble(k = 1:80) %>%
  mutate(
    kclust = map(k, ~kmeans(morgan_mat, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, morgan_mat)
  )

clusterings <- 
  kclusts %>%
  unnest(cols = c(glanced))

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  geom_text(aes(label = k), nudge_y = 90, nudge_x = 0.3)

ggsave(here('exploration/Fingerprints_results', 'silhouette_plot_kmeans.pdf'), 
       height = 9, width = 10)


### K-means with PC from PCA ####

pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, cumulative)) +
  geom_line(color = "#56B4E9", alpha = 0.8) +
  geom_point(color = '#56B4E9', alpha = 0.6) +
  # scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)


pca_PC = pca_fit %>%
  augment(morgan) %>% 
  # select(Drug, .fittedPC1:.fittedPC89)
  select(Drug, .fittedPC1:.fittedPC30)


pca_pc_mat = pca_PC %>% 
  select(where(is.numeric)) 

kclust = pca_pc_mat %>% 
  kmeans(centers = 10)


pca_fit %>%
  augment(morgan) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  ggplot(aes(.fittedPC1, .fittedPC2, color = cluster)) + 
  geom_point(size = 3.5) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()



pca_fit %>%
  augment(morgan) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster)) %>% 
  select(cluster, Drug) %>% 
  left_join(drugs) %>% 
  distinct(Drug, .keep_all=T) %>% 
  filter(Drug %in% nucl) %>% view


kclusts <- 
  tibble(k = 1:50) %>%
  mutate(
    kclust = map(k, ~kmeans(pca_pc_mat, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, morgan_mat)
  )

clusterings <- 
  kclusts %>%
  unnest(cols = c(glanced))

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point()






# Functional groups -------------------------------------------------------

dir.create(here('exploration', 'Functional_groups_results'))

## load the Morgan fingerprints produced by the script: pyMol2FuncGroup.py 

fgroups = read_csv("func_groups_biolog.csv")

# remove dactinomycin as it's very different from all the others
fgroups = fgroups %>% 
  filter(Drug != 'Dactinomycin')

### PCA ####

pca_fit = fgroups %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(scale = F) # do PCA on scaled data


pca_fit %>%
  augment(fgroups) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, .fittedPC2)) + 
  geom_point(size = 1.5) +
  geom_label(aes(label = Drug)) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()

ggsave(here('exploration/Functional_groups_results', 'PCA_fgroups.pdf'), 
       height = 13, width = 14)


pca_fit %>%
  augment(fgroups) %>% 
  plot_ly(x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3,
          text = ~Drug) %>% 
  add_markers()



pca_fit %>%
  tidy(matrix = "eigenvalues")

# barplot of percent explained by PC
pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

ggsave(here('exploration/Functional_groups_results', 'percent_explained.pdf'), 
       height = 9, width = 10)

# cumulative expl by PC 
pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, cumulative)) +
  geom_line(color = "#56B4E9", alpha = 0.8) +
  geom_point(color = 'black', alpha = 0.7) +
  geom_text(aes(label = round(PC,2)), nudge_x = 0, nudge_y = 0.05) +
  # scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

ggsave(here('exploration/Functional_groups_results', 'percent_explained_cumulative.pdf'), 
       height = 9, width = 10)

### distances ####

drug_dists = fgroups %>% select(where(is.numeric)) %>%  data.frame %>% get_dist
fviz_dist(drug_dists, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

ggsave(here('exploration/Functional_groups_results', 'distances.pdf'), 
       height = 9, width = 10)


### K-means with morgan FP ####

morgan_mat = fgroups %>% 
  select(where(is.numeric)) 



kclust = morgan_mat %>% kmeans(centers = 10)

pca_fit %>%
  augment(fgroups) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  ggplot(aes(.fittedPC1, .fittedPC2, color = cluster)) + 
  geom_point(size = 3.5) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()


# 3D plot 
pca_fit %>%
  augment(fgroups) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  plot_ly(x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3,
          text = ~Drug, color=~cluster) %>% 
  add_markers()



pca_fit %>%
  augment(fgroups) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster)) %>% 
  select(cluster, Drug) %>% 
  left_join(drugs) %>% 
  distinct(Drug, .keep_all=T) %>% 
  filter(Drug %in% nucl) %>% view


kclusts <- 
  tibble(k = 1:50) %>%
  mutate(
    kclust = map(k, ~kmeans(morgan_mat, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, morgan_mat)
  )

clusterings <- 
  kclusts %>%
  unnest(cols = c(glanced))

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  geom_text(aes(label = k), nudge_y = 90, nudge_x = 0.3)

ggsave(here('exploration/Functional_groups_results', 'silhouette_plot_kmeans.pdf'), 
       height = 9, width = 10)






# BIG PCA -----------------------------------------------------------------



dir.create(here('exploration', 'DrugBank_biolog_drugs'))

## load the Morgan fingerprints produced by the script: pyMol2FuncGroup.py 

all_compounds = fgroups %>% 
  mutate(DB = 'Biolog', .before = 'Al_COO') %>% 
  bind_rows(drugbank_fp %>% 
              mutate(DB = 'DrugBank', .before = 'Al_COO'))

# remove dactinomycin as it's very different from all the others
all_compounds = all_compounds %>% 
  filter(Drug != 'Dactinomycin')

### PCA ####

pca_fit = all_compounds %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(scale = F) # do PCA on scaled data


pca_fit %>%
  augment(all_compounds) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, .fittedPC2)) + 
  geom_point(size = 1.5, alpha = 0.3) +
  # geom_label(aes(label = Drug)) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()

ggsave(here('exploration/DrugBank_biolog_drugs', 'PCA_fgroups.pdf'), 
       height = 13, width = 14)


pca_fit %>%
  augment(all_compounds) %>% 
  plot_ly(x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3,
          text = ~Drug) %>% 
  add_markers()



pca_fit %>%
  tidy(matrix = "eigenvalues")

# barplot of percent explained by PC
pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

ggsave(here('exploration/DrugBank_biolog_drugs', 'percent_explained.pdf'), 
       height = 9, width = 10)

# cumulative expl by PC 
pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, cumulative)) +
  geom_line(color = "#56B4E9", alpha = 0.8) +
  geom_point(color = 'black', alpha = 0.7) +
  geom_text(aes(label = round(PC,2)), nudge_x = 0, nudge_y = 0.05) +
  # scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

ggsave(here('exploration/DrugBank_biolog_drugs', 'percent_explained_cumulative.pdf'), 
       height = 9, width = 10)

### distances ####

drug_dists = all_compounds %>% select(where(is.numeric)) %>%  data.frame %>% get_dist
fviz_dist(drug_dists, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

ggsave(here('exploration/DrugBank_biolog_drugs', 'distances.pdf'), 
       height = 19, width = 21)


### K-means with morgan FP ####

morgan_mat = all_compounds %>% 
  select(where(is.numeric)) 



kclust = morgan_mat %>% kmeans(centers = 10)

pca_fit %>%
  augment(all_compounds) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  ggplot(aes(.fittedPC1, .fittedPC2, color = cluster)) + 
  geom_point(size = 3.5) +
  labs(
    x = 'PC1',
    y = 'PC2'
  ) +
  theme_half_open(12) + 
  background_grid()


# 3D plot 
pca_fit %>%
  augment(all_compounds) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  plot_ly(x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3,
          text = ~Drug, color=~cluster) %>% 
  add_markers()


# 3D plot with Biolog/DrugBank info
pca_fit %>%
  augment(all_compounds) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster))  %>% 
  plot_ly(x = ~.fittedPC1, y = ~.fittedPC2, z = ~.fittedPC3,
          text = ~Drug, color=~DB) %>% 
  add_markers()



pca_fit %>%
  augment(fgroups) %>% # add original dataset back in
  mutate(cluster = kclust$cluster, .before = Drug,
         cluster = as.factor(cluster)) %>% 
  select(cluster, Drug) %>% 
  left_join(drugs) %>% 
  distinct(Drug, .keep_all=T) %>% 
  filter(Drug %in% nucl) %>% view


kclusts <- 
  tibble(k = 1:50) %>%
  mutate(
    kclust = map(k, ~kmeans(morgan_mat, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, morgan_mat)
  )

clusterings <- 
  kclusts %>%
  unnest(cols = c(glanced))

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  geom_text(aes(label = k), nudge_y = 90, nudge_x = 0.3)

ggsave(here('exploration/Functional_groups_results', 'silhouette_plot_kmeans.pdf'), 
       height = 9, width = 10)























