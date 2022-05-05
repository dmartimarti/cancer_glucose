
# libraries ---------------------------------------------------------------

library(tidyverse)
library(openxlsx)
library(readr)
library(plotly)
library(cowplot)
library(glue)

theme_set(theme_cowplot(15))


# read teh data ------------------------------------------------------------

# B1. Mechanism of action
# C4. Biological processes
# D1. Gene expression
# D2. Cancer cell lines

B1 = read_csv("drugbank_B1_CC.csv")
C4 = read_csv("drugbank_C4_CC.csv")
D1 = read_csv("drugbank_D1_CC.csv")
D2 = read_csv("drugbank_D2_CC.csv")
glob_cc = read_csv("drugbank_GLOBAL_CC.csv")


# plot some samples -------------------------------------------------------

B1 %>%
  plot_ly(x = ~tsne_1, y = ~tsne_2, z = ~tsne_3,
          text = ~Drug, color = ~origin,
          colors = c("#E6B324", "#2550E6"),
          alpha = 0.8
          ) %>% 
  add_markers()


C4 %>%
  plot_ly(x = ~tsne_1, y = ~tsne_2, z = ~tsne_3,
          text = ~Drug, color = ~origin,
          colors = c("#E6B324", "#2550E6"),
          alpha = 0.8
  ) %>% 
  add_markers()


D1 %>%
  plot_ly(x = ~tsne_1, y = ~tsne_2, z = ~tsne_3,
          text = ~Drug, color = ~origin,
          colors = c("#E6B324", "#2550E6"),
          alpha = 0.8
  ) %>% 
  add_markers()



D2 %>%
  plot_ly(x = ~tsne_1, y = ~tsne_2, z = ~tsne_3,
          text = ~Drug, color = ~origin,
          colors = c("#E6B324", "#2550E6"),
          alpha = 0.8
  ) %>% 
  add_markers()



glob_cc %>%
  plot_ly(x = ~tsne_1, y = ~tsne_2, z = ~tsne_3,
          text = ~Drug, color = ~origin,
          colors = c("#E6B324", "#2550E6"),
          alpha = 0.8
  ) %>% 
  add_markers()



# euclidean distances --------------------------------------

# calculate euclidean distances from Micit

euclidean <- function(a, b) sqrt(sum((a - b)^2))


## test with two compounds ####

micit = B1 %>% 
  filter(Drug == 'Micit') %>% 
  select(tsne_1:tsne_3) %>% 
  pivot_longer(cols = tsne_1:tsne_3,
               names_to = 'dims', values_to = 'vals') %>% 
  pull(vals)

# far compound
hist = B1 %>% 
  filter(Drug == 'Histidine') %>% 
  select(tsne_1:tsne_3) %>% 
  pivot_longer(cols = tsne_1:tsne_3,
               names_to = 'dims', values_to = 'vals') %>% 
  pull(vals)

# close compound
acal = B1 %>% 
  filter(Drug == 'Acalabrutinib') %>% 
  select(tsne_1:tsne_3) %>% 
  pivot_longer(cols = tsne_1:tsne_3,
               names_to = 'dims', values_to = 'vals') %>% 
  pull(vals)

euclidean(micit, hist)
euclidean(micit, acal)



## B1 distances  ####
# computing the distances for all compounds in the table

micit_B1 = B1 %>% 
  filter(Drug == 'Micit') %>% 
  select(tsne_1:tsne_3) %>% 
  pivot_longer(cols = tsne_1:tsne_3,
               names_to = 'dims', values_to = 'vals') %>% 
  pull(vals)

B1_dist = B1 %>% 
  select(-tsne_2d_1, -tsne_2d_2) %>% 
  group_by(Drug, smiles, origin) %>% 
  nest() %>% 
  mutate(dist = map(.x= data, .f = euclidean, b = micit_B1)) %>% 
  unnest(cols = c(data,dist))

B1_dist %>% filter(Drug == 'Acalabrutinib')

closest_B1 = B1_dist %>% 
  arrange(dist) %>% 
  head(30)

B1_dist %>% 
  filter(Drug != 'Micit') %>% 
  ggplot(aes(x = fct_reorder(Drug, dist), y = dist)) + 
  geom_point(alpha = 0.2) +
  theme(axis.text.x = element_blank(),
        axis.ticks.length.x = NULL)



## C4 distances  ####
# computing the distances for all compounds in the table

micit_C4 = C4 %>% 
  filter(Drug == 'Micit') %>% 
  select(tsne_1:tsne_3) %>% 
  pivot_longer(cols = tsne_1:tsne_3,
               names_to = 'dims', values_to = 'vals') %>% 
  pull(vals)

C4_dist = C4 %>% 
  select(-tsne_2d_1, -tsne_2d_2) %>% 
  group_by(Drug, smiles, origin) %>% 
  nest() %>% 
  mutate(dist = map(.x= data, .f = euclidean, b = micit_C4)) %>% 
  unnest(cols = c(data,dist))

closest_C4 = C4_dist %>% 
  arrange(dist) %>% 
  head(30)

C4_dist %>% 
  filter(Drug != 'Micit') %>% 
  ggplot(aes(x = fct_reorder(Drug, dist), y = dist)) + 
  geom_point(alpha = 0.2) +
  theme(axis.text.x = element_blank(),
        axis.ticks.length.x = NULL)



## D1 distances  ####
# computing the distances for all compounds in the table

micit_D1 = D1 %>% 
  filter(Drug == 'Micit') %>% 
  select(tsne_1:tsne_3) %>% 
  pivot_longer(cols = tsne_1:tsne_3,
               names_to = 'dims', values_to = 'vals') %>% 
  pull(vals)

D1_dist = D1 %>% 
  select(-tsne_2d_1, -tsne_2d_2) %>% 
  group_by(Drug, smiles, origin) %>% 
  nest() %>% 
  mutate(dist = map(.x= data, .f = euclidean, b = micit_C4)) %>% 
  unnest(cols = c(data,dist))

closest_D1 = D1_dist %>% 
  arrange(dist) %>% 
  head(30)

D1_dist %>% 
  filter(Drug != 'Micit') %>% 
  ggplot(aes(x = fct_reorder(Drug, dist), y = dist)) + 
  geom_point(alpha = 0.2) +
  theme(axis.text.x = element_blank(),
        axis.ticks.length.x = NULL)


## D2 distances  ####
# computing the distances for all compounds in the table

micit_D2 = D2 %>% 
  filter(Drug == 'Micit') %>% 
  select(tsne_1:tsne_3) %>% 
  pivot_longer(cols = tsne_1:tsne_3,
               names_to = 'dims', values_to = 'vals') %>% 
  pull(vals)

D2_dist = D2 %>% 
  select(-tsne_2d_1, -tsne_2d_2) %>% 
  group_by(Drug, smiles, origin) %>% 
  nest() %>% 
  mutate(dist = map(.x= data, .f = euclidean, b = micit_C4)) %>% 
  unnest(cols = c(data,dist))

closest_D2 = D2_dist %>% 
  arrange(dist) %>% 
  head(30)

D2_dist %>% 
  filter(Drug != 'Micit') %>% 
  ggplot(aes(x = fct_reorder(Drug, dist), y = dist)) + 
  geom_point(alpha = 0.2) +
  theme(axis.text.x = element_blank(),
        axis.ticks.length.x = NULL)

