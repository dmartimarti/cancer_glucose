
# libraries ---------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(here)
library(readxl)
library(DescTools)
library(viridisLite)


theme_set(theme_cowplot(15))

# load data ---------------------------------------------------------------


clean_data = read_excel("clean_data.xlsx", 
                         sheet = "Sheet1")


# make it longer
clean_data = clean_data %>% 
  pivot_longer(cols = `0`:`10`, 
               names_to = 'FU', values_to = 'confluence') %>% 
  mutate(Micit = factor(Micit, levels = c(10, 5, 1, 0)),
         FU = factor(FU, levels = c(0, 0.04, 0.08, 0.16, 0.31, 0.63, 1.25, 2.5, 5, 10)))


# plot 
clean_data %>% 
  mutate(Elapsed = round(Elapsed/60, 1)) %>% 
  ggplot(aes(x = Elapsed, y = confluence)) +
  geom_line(aes(color = Micit)) +
  # facet_wrap(vars(Micit, FU), ncol = 10)
  facet_grid(rows = vars(Micit), 
             cols = vars(FU)) +
  labs(
    x = 'Time (hours)',
    y = 'Confluence (%)'
  )



ggsave(here('summary', 'growth_curves_rep1.pdf'),
       height = 6, width = 8)


# calculate AUCs ----------------------------------------------------------


conf = clean_data %>% 
  filter(FU == 0, Micit == 0) %>% 
  pull(confluence)

times = clean_data %>% 
  filter(FU == 0, Micit == 0) %>% 
  pull(Elapsed)



data_auc = clean_data %>%
  group_by(Micit, FU, Replicate) %>% 
  summarise(auc = AUC(Elapsed, confluence))


control = data_auc %>% 
  filter(FU == 0, Micit == 0) %>% pull(auc)



# calculate viability
data_viab = data_auc %>% 
  mutate(control = data_auc %>% 
           filter(FU == 0, Micit == 0) %>% 
           pull(auc)) %>% 
  mutate(viability = ((auc/control) * 100))



data_viab %>% 
  ggplot(aes(y = Micit, x = FU, fill = viability)) +
  geom_raster() +
  scale_fill_viridis_c() +
  labs(
    x = '5-FU',
    y = 'Micit'
  )

ggsave(here('summary', 'auc_heatmap_rep1.pdf'),
       height = 6, width = 8)



synfinder = data_viab %>% 
  rename(conc1 = Micit, 
         conc2 = FU,
         response = viability) %>% 
  mutate(drug1 = 'Micit',
         drug2 = '5-FU',
         conc_unit = 'A.U.',
         block_id = 1) %>% 
  select(drug1, conc1, drug2, conc2, response, conc_unit, block_id)


synfinder %>% 
  write_csv(here('summary', 'SynergyFinder_incucyte.csv'))













