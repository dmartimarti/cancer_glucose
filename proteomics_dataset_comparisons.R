
# libraries ---------------------------------------------------------------

library(tidyverse)
library(readr)
library(tidymodels)

theme_set(theme_light())

# read data ---------------------------------------------------------------

# data from https://gygi.hms.harvard.edu/publications/ccle.html
protein_quant = read_csv("summed_sn_non_normalized.csv")
protein_quant_norm = read_csv("protein_quant_current_normalized.csv")

# our data

data = read_csv("means_samples.csv")


## some data cleanup

protein_quant_norm = protein_quant_norm %>%
  select(-c(TenPx01_Peptides:TenPx24_Peptides),
         -Group_ID) %>%
  rename(Gene_names = Gene_Symbol)


protein_quant = protein_quant %>%
  select(-(starts_with('bridge'))) %>%
  rename(Gene_names = Gene.Symbol)





## check the gene names, do my gene names match all present in the dataset?




# there are names that are collapsed as "name1;name2", need to separate rows

data = data %>%
  separate_rows(Gene_names, sep = ';') 


my_names = unique(data$Gene_names)
other_names = unique(protein_quant$Gene_names)

setdiff(my_names, other_names)

# now is better as we gain some info. 





# PCA from big df -----------------------------------------------------------



# names(protein_quant)

prot_mat_full = protein_quant %>% 
  select(MDAMB468_BREAST:THP1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE) %>%
  replace(is.na(.), 0) 

# merge columns
protein_quant_full = bind_cols(protein_quant %>% 
            select(Protein.Id:Gene_names),
          prot_mat_full)



# select numeric columns, make them a numeric matrix
prot_mat = protein_quant_full %>%
  select(MDAMB468_BREAST:THP1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE) %>%
  as.matrix

# transpose the matrix
prot_mat = t(prot_mat)

# get the name of the cell lines and genes
cells = names(prot_mat_full)
genes = protein_quant_full$Protein.Id

# check dimensions
dim(prot_mat)
length(cells)
length(genes)

prot_mat = as_tibble(prot_mat) %>%
  mutate(Cells = cells, .before = V1)

colnames(prot_mat) = c('Cells',genes)

# FINAL VERSION OF THE MATRIX
prot_mat




pca_rec = recipe(~., data = prot_mat) %>%
  update_role(Cells, new_role = "id") %>%
  step_normalize(all_numeric()) %>%
  step_pca(all_numeric())

pca_prep = prep(pca_rec)

pca_prep


# tidied_pca = tidy(pca_prep, 2)
# 
# tidied_pca %>%
#   filter(component %in% paste0("PC", 1:5)) %>%
#   mutate(component = fct_inorder(component)) %>%
#   ggplot(aes(value, terms, fill = terms)) +
#   geom_col(show.legend = FALSE) +
#   facet_wrap(~component, nrow = 1) +
#   labs(y = NULL)



library(tidytext)

tidied_pca %>%
  filter(component %in% paste0("PC", 1:4)) %>%
  group_by(component) %>%
  top_n(8, abs(value)) %>%
  ungroup() %>%
  mutate(terms = reorder_within(terms, abs(value), component)) %>%
  ggplot(aes(abs(value), terms, fill = value > 0)) +
  geom_col() +
  facet_wrap(~component, scales = "free_y") +
  scale_y_reordered() +
  labs(
    x = "Absolute value of contribution",
    y = NULL, fill = "Positive?"
  )


# plot the PCA
juice(pca_prep) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(alpha = 0.7, size = 2)









#  join datasets ----------------------------------------------------------

### I need to scale the datasets to the same extent, our dataset is in log2
# format, whereas the others are in a non-normalised format. First, I need to 
# get the min value (not 0), substitute all 0s by a fraction of that value, 
# and then divide everything by that fraction so the 0s will be 1s, and the log2
# will be then 0, and not a negative value


# data %>% full_join(protein_quant)
removals = setdiff(my_names, other_names)



prot_mat_full = protein_quant %>% 
  filter(Gene_names %in% data$Gene_names) %>%
  select(MDAMB468_BREAST:THP1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE) 

# ugly code to get the min value
min_val = sort(unlist(unname(prot_mat_full)))[sort(unlist(unname(prot_mat_full))) > 0][1]

min_val = min_val / 1

prot_mat_full = prot_mat_full %>%
  replace(is.na(.), min_val) 

# replace 0
prot_mat_full[prot_mat_full == 0] = min_val
# log2 transform
prot_mat_full = log2(prot_mat_full / min_val)

# Normalized Data
prot_mat_full = (prot_mat_full-min(prot_mat_full))/(max(prot_mat_full)-min(prot_mat_full))

# set up as tibble
prot_mat_full = as_tibble(prot_mat_full)


# merge columns



protein_quant_full = bind_cols(protein_quant %>% 
                                 select(Protein.Id:Gene_names) %>%
                                 filter(Gene_names %in% data$Gene_names),
                               prot_mat_full)

data_mat = data %>%select(Micit_10:Control)

# Normalized Data
data_mat = (data_mat-min(data_mat))/(max(data_mat)-min(data_mat))


sort(unlist(unname(data_mat)))


data_norm = bind_cols(data$Gene_names,as_tibble(data_mat)) %>% 
  rename(Gene_names=`...1`)

 

protein_quant_full = data_norm %>% full_join(protein_quant_full) %>%
  filter(!(Gene_names %in% removals)) %>%
  select(Gene_names,Protein.Id,
    Micit_10:Control,
    MDAMB468_BREAST:THP1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE)







# select numeric columns, make them a numeric matrix
prot_mat = protein_quant_full %>%
  select(Micit_10:THP1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE) %>%
  as.matrix

# transpose the matrix
prot_mat = t(prot_mat)

# get the name of the cell lines and genes
cells = names(protein_quant_full[,3:386])
genes = protein_quant_full$Protein.Id

# check dimensions
dim(prot_mat)
length(cells)
length(genes)

prot_mat = as_tibble(prot_mat) %>%
  mutate(Cells = cells, .before = V1)

colnames(prot_mat) = c('Cells',genes)

# FINAL VERSION OF THE MATRIX
prot_mat


prot_mat = prot_mat %>% 
  replace(is.na(.), 0) 


prot_mat %>% select(Cells, `sp|Q2M2I8|AAK1_HUMAN`) %>%
  # filter(!(Cells %in% lines)) %>%
  ggplot(aes(x = `sp|Q2M2I8|AAK1_HUMAN`)) +
  geom_histogram()


### CALCULATE PCA

pca_rec = recipe(~., data = prot_mat) %>%
  update_role(Cells, new_role = "id") %>%
  step_normalize(all_numeric()) %>%
  step_pca(all_numeric())

pca_prep = prep(pca_rec)

pca_prep



library(tidytext)

tidied_pca %>%
  filter(component %in% paste0("PC", 1:4)) %>%
  group_by(component) %>%
  top_n(8, abs(value)) %>%
  ungroup() %>%
  mutate(terms = reorder_within(terms, abs(value), component)) %>%
  ggplot(aes(abs(value), terms, fill = value > 0)) +
  geom_col() +
  facet_wrap(~component, scales = "free_y") +
  scale_y_reordered() +
  labs(
    x = "Absolute value of contribution",
    y = NULL, fill = "Positive?"
  )


# plot the PCA

lines = c('Micit_10','Micit_1',
  'FU','FU_10mM_Micit',
  'FU_1mM_Micit',
  'Control',
  'HCT116_LARGE_INTESTINE_TenPx04')

juice(pca_prep) %>%
  mutate(labels = case_when(Cells %in% lines ~ Cells),
         Class = case_when(Cells %in% lines ~ 'Hightligh',
                           TRUE ~ 'Others')) %>%
  ggplot(aes(PC1, PC2, label = labels)) +
  geom_text(check_overlap = F) +
  geom_point(aes(color = Class), size = 2)




# scaling data -------------------------------------------------------------


protein_quant %>%
  select(Gene_names, HCT116_LARGE_INTESTINE)%>%
  drop_na() %>%
  filter(Gene_names %in% data$Gene_names) %>%
  mutate(lg2 = log2(HCT116_LARGE_INTESTINE)) %>%
  ggplot(aes(lg2)) +
  geom_histogram()






data %>%
  ggplot(aes(Control)) +
  geom_histogram()




















