# libraries --------------------------------------------------------------

library(tidyverse)
library(readr)
library(readxl)
theme_set(theme_light())
library(openxlsx)
library(here)
library(ggrepel)
# library(tidymodels)

# load data ---------------------------------------------------------------

keio.rat = read_excel("D:/MRC_Postdoc/Pangenomic/biolog/KEIO/fixed_data/KEIO_output/summary/AUC_KEIO_screen.xlsx", 
                              sheet = "ratios") %>% select(-`...1`)

aska.rat = read_excel("D:/MRC_Postdoc/Pangenomic/biolog/ASKA/fixed_data/Output/summary/AUC_raw_ASKAscreen.xlsx", 
                                 sheet = "ratios_by_col") %>% select(-`...1`)

# select data
keio.rat = keio.rat %>% select(Strain, keio_mean = Mean_correct)
aska.rat = aska.rat %>% select(Strain, aska_mean = Mean_correct)

merg_data = keio.rat %>% left_join(aska.rat) %>% drop_na(keio_mean, aska_mean)


# mean as the center of the circle
x.1 = keio.rat %>% summarise(Mean = mean(keio_mean))
y.1 = aska.rat %>% summarise(Mean = mean(aska_mean))

# data transformations to separate hits
merg_data = merg_data %>% 
  mutate(x = x.1$Mean,
         y = y.1$Mean,
         distance = sqrt((x - keio_mean)**2 + (y - aska_mean)**2),
    ratios = keio_mean/aska_mean,
         lograt = log10(ratios),
         rest = abs(keio_mean - aska_mean)) 

# get top resistant/sensitive hits
n = 40
bot_keio = merg_data %>% arrange(keio_mean) %>% head(n) %>% select(Strain) %>% t %>% as.character()
top_keio = merg_data %>% arrange(desc(keio_mean)) %>% head(n) %>% select(Strain) %>% t %>% as.character()
bot_aska = merg_data %>% arrange(aska_mean) %>% head(n) %>% select(Strain) %>% t %>% as.character()
top_aska = merg_data %>% arrange(desc(aska_mean)) %>% head(n) %>% select(Strain) %>% t %>% as.character()

lab = unique(c(bot_aska, bot_keio, top_aska, top_keio))

# label the 
merg_data = merg_data %>% 
  mutate(labels = case_when(rest > 0.38 ~ Strain,
                            Strain %in% lab ~ Strain))


merg_data %>% 
  # filter(!(Strain %in% c('yqaA', 'yajC'))) %>% 
  ggplot(aes(x = aska_mean, y = keio_mean, size = rest, colour = distance)) +
  geom_point() +
  geom_text_repel(aes(label = labels))



