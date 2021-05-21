library(readr)
library(multcomp)
library(tidyverse)
library(here)
library(ComplexHeatmap)
library(openxlsx)
# library(PFun)
library(broom)
# library(multidplyr)
library(tau)
library(fmsb)
library(readxl)



theme_set(theme_classic())

# load the data -----------------------------------------------------------

data = read_delim("stats_perseus_full.tsv", 
                  "\t", escape_double = FALSE, trim_ws = TRUE)

names(data)

names(data) = gsub(' ', '_', names(data))

data




### TABLE TIDYING
# create groups depending on their DE profiles

data = data %>% 
  mutate(FU_C = case_when(`Student's_T-test_Significant_5FU_Control` == '+' & `Student's_T-test_Difference_5FU_Control` < 0 ~ -1,
                          `Student's_T-test_Significant_5FU_Control` == '+' & `Student's_T-test_Difference_5FU_Control` > 0 ~ 1,
                          TRUE ~ 0),
         micit1mM_C = case_when(`Student's_T-test_Significant_1mM_Micit_Control` == '+' & `Student's_T-test_Difference_1mM_Micit_Control` < 0 ~ -1,
                                `Student's_T-test_Significant_1mM_Micit_Control` == '+' & `Student's_T-test_Difference_1mM_Micit_Control` > 0 ~ 1,
                                 TRUE ~ 0),
         micit10mM_C = case_when(`Student's_T-test_Significant_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_10mM_Micit_Control` < 0 ~ -1,
                                 `Student's_T-test_Significant_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_10mM_Micit_Control` > 0 ~ 1,
                                  TRUE ~ 0),
         micit1mM_5FU_C = case_when(`Student's_T-test_Significant_5FU_1mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_1mM_Micit_Control` < 0 ~ -1,
                                    `Student's_T-test_Significant_5FU_1mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_1mM_Micit_Control` > 0 ~ 1,
                                     TRUE ~ 0),
         micit10mM_5FU_C = case_when(`Student's_T-test_Significant_5FU_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_10mM_Micit_Control` < 0 ~ -1,
                                     `Student's_T-test_Significant_5FU_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_10mM_Micit_Control` > 0 ~ 1,
                                     TRUE ~ 0),
         .after = "KEGG_name") 


# let's rename columns and remove the ones that are useless

names(data)

data_tidy = data %>% 
  select(Protein_IDs:Gene_names, Control:KEGG_name, -`Q-value`,
         mol_weight = `Mol._weight_[kDa]`,
         # 5FU_control
         logPval_FU_C = `-Log_Student's_T-test_p-value_5FU_Control`,
         Qval_FU_C = `Student's_T-test_q-value_5FU_Control`,
         difference_FU_C = `Student's_T-test_Difference_5FU_Control`,
         # 1 mM micit
         logPval_micit1mM_C = `-Log_Student's_T-test_p-value_1mM_Micit_Control`,
         Qval_micit1mM_C = `Student's_T-test_q-value_1mM_Micit_Control`,
         difference_micit1mM_C = `Student's_T-test_Difference_1mM_Micit_Control`,
         # 10 mM micit
         logPval_micit10mM_C = `-Log_Student's_T-test_p-value_10mM_Micit_Control`,
         Qval_micit10mM_C = `Student's_T-test_q-value_10mM_Micit_Control`,
         difference_micit10mM_C = `Student's_T-test_Difference_10mM_Micit_Control`,
         # 5FU + 1 mM micit
         logPval_micit1mM_5FU_C = `-Log_Student's_T-test_p-value_5FU_1mM_Micit_Control`,
         Qval_micit1mM_5FU_C = `Student's_T-test_Test_statistic_5FU_1mM_Micit_Control`,
         difference_micit1mM_5FU_C = `Student's_T-test_Difference_5FU_1mM_Micit_Control`,
         # 5FU + 10 mM micit
         logPval_micit10mM_5FU_C = `-Log_Student's_T-test_p-value_5FU_10mM_Micit_Control`,
         Qval_micit10mM_5FU_C = `Student's_T-test_q-value_5FU_10mM_Micit_Control`,
         difference_micit10mM_5FU_C = `Student's_T-test_Difference_5FU_10mM_Micit_Control`
  )


write_csv(data_tidy, here('summary','Stats_with_differences_full.csv'))




# convert data to long format

names(data)

data_long = data_tidy %>% 
  pivot_longer(Control:`5FU_10mM_Micit_4`, names_to = 'Sample', values_to = 'Intensity') 

# fix variable names
data_long = data_long %>% select(Gene_names, Sample, Intensity, everything()) %>% 
  mutate(Sample = case_when(Sample == 'Control_1' | Sample == 'Control_2' | Sample == 'Control_3' | Sample == 'Control_4' ~ 'Control',
                            Sample == '5FU_1' | Sample == '5FU_2' | Sample == '5FU_3' | Sample == '5FU_4' ~ '5FU',
                            Sample == '1mM_Micit_1' | Sample == '1mM_Micit_2' | Sample == '1mM_Micit_3' | Sample == '1mM_Micit_4' ~ '1mM_Micit',
                            Sample == '5FU_1mM_Micit_1' | Sample == '5FU_1mM_Micit_2' | Sample == '5FU_1mM_Micit_3' | Sample == '5FU_1mM_Micit_4' ~ '5FU_1mM_Micit',
                            Sample == '10mM_Micit_1' | Sample == '10mM_Micit_2' | Sample == '10mM_Micit_3' | Sample == '10mM_Micit_4' ~ '10mM_Micit',
                            Sample == '5FU_10mM_Micit_1' | Sample == '5FU_10mM_Micit_2' | Sample == '5FU_10mM_Micit_3' | Sample == '5FU_10mM_Micit_4' ~ '5FU_10mM_Micit',
                            TRUE ~ Sample),
         Sample = factor(Sample)) %>% 
  arrange(desc(Sample))

data_long

data_long %>% 
  select(Gene_names, Sample, Intensity) %>% 
  write_csv(here('summary','data_long.csv'))


# UpSet general -----------------------------------------------------------

library(UpSetR)


### global differences ####

data

FU_genes = data %>% 
  filter(`Student's_T-test_Significant_5FU_Control` == '+') %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character

micit_10mM_genes = data %>% 
  filter(`Student's_T-test_Significant_10mM_Micit_Control` == '+') %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character

micit_1mM_5FU_genes = data %>% 
  filter(`Student's_T-test_Significant_5FU_1mM_Micit_Control` == '+') %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character

micit_10mM_5FU_genes = data %>% 
  filter(`Student's_T-test_Significant_5FU_10mM_Micit_Control` == '+') %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character




list_full = list(
  '5FU' = FU_genes,
  'Micit 10mM' = micit_10mM_genes,
  '5FU+Micit 1mM' = micit_1mM_5FU_genes,
  '5FU+Micit 10mM' = micit_10mM_5FU_genes
)


# generate a combination matrix
m1 = make_comb_mat(list_full)

m2 = make_comb_mat(list_full, mode = "intersect")

# size of the groups
comb_size(m1)
comb_size(m2)


UpSet(m1)


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'upset_full_distinct.pdf'),
             width = 8, height = 6, useDingbats = FALSE)


UpSet(m2)


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'upset_full_intersect.pdf'),
             width = 8, height = 6, useDingbats = FALSE)







# upset up/down -----------------------------------------------------------


# 5FU
FU_genes_down = data %>% 
  filter(`Student's_T-test_Significant_5FU_Control` == '+' & `Student's_T-test_Difference_5FU_Control` < 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character

FU_genes_up = data %>% 
  filter(`Student's_T-test_Significant_5FU_Control` == '+' & `Student's_T-test_Difference_5FU_Control` > 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character


# micit 10 mM
micit_10mM_genes_down = data %>% 
  filter(`Student's_T-test_Significant_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_10mM_Micit_Control` < 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character

micit_10mM_genes_up = data %>% 
  filter(`Student's_T-test_Significant_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_10mM_Micit_Control` > 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character


# 5FU + 1mM micit
micit_1mM_5FU_genes_down = data %>% 
  filter(`Student's_T-test_Significant_5FU_1mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_1mM_Micit_Control` < 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character

micit_1mM_5FU_genes_up = data %>% 
  filter(`Student's_T-test_Significant_5FU_1mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_1mM_Micit_Control` > 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character


# 5FU + 10mM micit
micit_10mM_5FU_genes_down = data %>% 
  filter(`Student's_T-test_Significant_5FU_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_10mM_Micit_Control` < 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character

micit_10mM_5FU_genes_up = data %>% 
  filter(`Student's_T-test_Significant_5FU_10mM_Micit_Control` == '+' & `Student's_T-test_Difference_5FU_10mM_Micit_Control` > 0) %>% 
  select(Gene_names)  %>% drop_na(Gene_names) %>% t %>% as.character




list_directions = list(
  '5FU down' = FU_genes_down,
  '5FU up' = FU_genes_up,
  'Micit 10 mM down' = micit_10mM_genes_down,
  'Micit 10 mM up' = micit_10mM_genes_up,
  '5FU + Micit 1 mM down' = micit_1mM_5FU_genes_down,
  '5FU + Micit 1 mM up' = micit_1mM_5FU_genes_up,
  '5FU + Micit 10 mM down' = micit_10mM_5FU_genes_down,
  '5FU + Micit 10 mM up' = micit_10mM_5FU_genes_up
)


# generate a combination matrix
m1 = make_comb_mat(list_directions)

m2 = make_comb_mat(list_directions, mode = "intersect")




UpSet(m1)


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'upset_directions_distinct.pdf'),
             width = 8, height = 6, useDingbats = FALSE)


UpSet(m2)


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'upset_directions_intersect.pdf'),
             width = 8, height = 6, useDingbats = FALSE)
  




## extracting gene sets ####-----------------------------------------------

# gene order:
# 5FU || 10mM Micit || 5FU+1mM Micit || 5FU+10mM Micit
# down/up

comb_name(m2, readable = T)

# 5FU ∩ 10 mM Micit
FU_Micit1_down = extract_comb(m2, "10100000")
FU_Micit1_up = extract_comb(m2, "01010000")

# 5FU ∩ 5FU+1mM Micit
FU_Micit1FU_down = extract_comb(m2, "10001000")
FU_Micit1FU_up = extract_comb(m2, "01000100")

# 5FU ∩ 5FU+10mM Micit
FU_Micit10FU_down = extract_comb(m2, "10000010")
FU_Micit10FU_up = extract_comb(m2, "01000001")

# 10mM micit ∩ 5FU+10mM Micit
Micit10_Micit10FU_down = extract_comb(m2, "00100010")
Micit10_Micit10FU_up = extract_comb(m2, "00010001")

## Specific combinations
# 5FU+micit 10mM ∩ 5FU+micit 1mM
Micit10FU_Micit1FU_down = extract_comb(m2, "00001010")
Micit10FU_Micit1FU_up = extract_comb(m2, "00000101")

# 5FU+10mM micit ∩ 5FU+1mM micit U 10mM Micit
trio_down = extract_comb(m2, "00101010")
trio_up = extract_comb(m2, "00010101")


# 5FU + 1mM ! 5FU


setdiff(micit_1mM_5FU_genes_up, FU_genes_up)
length(setdiff(micit_1mM_5FU_genes_up, FU_genes_up))
length(micit_1mM_5FU_genes_up)
length(FU_genes_up)

setdiff(micit_1mM_5FU_genes_down, FU_genes_down)
length(setdiff(micit_1mM_5FU_genes_down, FU_genes_down))
length(micit_1mM_5FU_genes_down)
length(FU_genes_down)



FU_Micit1_diff_down = setdiff(micit_1mM_5FU_genes_down, FU_genes_down)
FU_Micit1_diff_up = setdiff(micit_1mM_5FU_genes_up, FU_genes_up)




# make a file of genes in each contrast
explanations = c('5FU ∩ 10 mM Micit', '5FU ∩ 5FU+1mM Micit',
                 ' 5FU ∩ 5FU+10mM Micit', '10mM micit ∩ 5FU+10mM Micit',
                 '5FU+micit 10mM ∩ 5FU+micit 1mM', 
                 '5FU+10mM micit ∩ 5FU+1mM micit ∩ 10mM Micit',
                 '5FU + 1mM micit ! 5FU')

table_names = c('FU_Micit1', 'FU_Micit1FU', 'FU_Micit10FU', 'Micit10_Micit10FU',
                'Micit10FU_Micit1FU', 'trio', 'FU_Micit1_diff')

exp_df = tibble(explanations, table_names)

library(openxlsx)

list_of_tables = list(
  'Metadata' = exp_df,
  FU_Micit1_down = unlist(str_split(FU_Micit1_down, pattern = ';')),
  FU_Micit1_up = unlist(str_split(FU_Micit1_up, pattern = ';')),
  FU_Micit1FU_down = unlist(str_split(FU_Micit1FU_down, pattern = ';')),
  FU_Micit1FU_up = unlist(str_split(FU_Micit1FU_up, pattern = ';')),
  FU_Micit10FU_down = unlist(str_split(FU_Micit10FU_down, pattern = ';')),
  FU_Micit10FU_up = unlist(str_split(FU_Micit10FU_up, pattern = ';')),
  Micit10_Micit10FU_down = unlist(str_split(Micit10_Micit10FU_down, pattern = ';')),
  Micit10_Micit10FU_up = unlist(str_split(Micit10_Micit10FU_up, pattern = ';')),
  Micit10FU_Micit1FU_down = unlist(str_split(Micit10FU_Micit1FU_down, pattern = ';')),
  Micit10FU_Micit1FU_up =unlist(str_split(Micit10FU_Micit1FU_up, pattern = ';')),
  trio_down = unlist(str_split(trio_down, pattern = ';')),
  trio_up = unlist(str_split(trio_up, pattern = ';')),
  FU_Micit1_diff_down = unlist(str_split(FU_Micit1_diff_down, pattern = ';')),
  FU_Micit1_diff_up = unlist(str_split(FU_Micit1_diff_up, pattern = ';'))
)

write.xlsx(list_of_tables, here('summary','Protein_sets.xlsx'))


for (i in 2:length(list_of_tables)) {
  file_name = names(list_of_tables)[i]
  write.table(list_of_tables[[i]],here('summary/protein_groups_intersection', 
                                     paste0(file_name,'.txt')),
              quote=F,col.names = F,row.names = F)
}



# saving the original protein differences
list_of_tables = list(
  FU_DOWN = unlist(str_split(FU_genes_down, pattern = ';')),
  FU_UP = unlist(str_split(FU_genes_up, pattern = ';')),
  micit10_DOWN = unlist(str_split(micit_10mM_genes_down, pattern = ';')),
  micit10_UP = unlist(str_split(micit_10mM_genes_up, pattern = ';')),
  micit1_5FU_DOWN = unlist(str_split(micit_1mM_5FU_genes_down, pattern = ';')),
  micit1_5FU_UP = unlist(str_split(micit_1mM_5FU_genes_up, pattern = ';')),
  micit10_5FU_DOWN = unlist(str_split(micit_10mM_5FU_genes_down, pattern = ';')),
  micit10_5FU_UP = unlist(str_split(micit_10mM_5FU_genes_up, pattern = ';'))
)

write.xlsx(list_of_tables, here('summary','Proteins_directions_MainGroups.xlsx'))

# loop to save everything in a folder
for (i in 1:length(list_of_tables)) {
  file_name = names(list_of_tables)[i]
  write.table(list_of_tables[[i]],here('summary/protein_groups', 
                                       paste0(file_name,'.txt')),
              quote=F,col.names = F,row.names = F)
}



# plot examples of genes
gene = c('RPL5','EIF5A;EIF5AL1')
data_long %>% filter(Gene_names %in%  gene) %>% 
  # filter(Sample != 'Control') %>% 
  ggplot(aes(x = Sample, y = Intensity, fill = Sample)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  labs(title = gene) +
  facet_wrap(~Gene_names)





# radarplot ---------------------------------------------------------------


### LOAD DATA
# 1mM micit + 5FU
micit_1mM_5FU_enrich_down = read_excel("summary/protein_groups/output/micit_1mM_5FU_genes_down_output.xlsx", 
                                              sheet = "Process") %>% 
  select(-`...1`)

micit_1mM_5FU_enrich_up = read_excel("summary/protein_groups/output/micit_1mM_5FU_genes_up_output.xlsx", 
                                       sheet = "Process") %>% 
  select(-`...1`)

# 10mM micit + 5FU
micit_10mM_5FU_enrich_down = read_excel("summary/protein_groups/output/micit_10mM_5FU_genes_down_output.xlsx", 
                                               sheet = "Process")%>% 
  select(-`...1`)

micit_10mM_5FU_enrich_up = read_excel("summary/protein_groups/output/micit_10mM_5FU_genes_up_output.xlsx", 
                                        sheet = "Process")%>% 
  select(-`...1`)

# 10mM micit
micit_10mM_enrich_down = read_excel("summary/protein_groups/output/micit_10mM_genes_down_output.xlsx", 
                                        sheet = "Process")%>% 
  select(-`...1`)

micit_10mM_enrich_up = read_excel("summary/protein_groups/output/micit_10mM_genes_up_output.xlsx", 
                                      sheet = "Process")%>% 
  select(-`...1`)

# 5FU
FU_enrich_down = read_excel("summary/protein_groups/output/FU_genes_down_output.xlsx", 
                                   sheet = "Process") %>% 
  select(-`...1`)

FU_enrich_up = read_excel("summary/protein_groups/output/FU_genes_up_output.xlsx", 
                             sheet = "Process") %>% 
  select(-`...1`)


### FUNCTIONS 
terms_in_words = function(terms, elms=10){
  word_count = textcnt(terms,  split = ' ',method = "string", 
                       n = 1L)
  
  df = data.frame(matrix(word_count)) %>%
    mutate(words = names(word_count)) 
  
  names(df) = c('count', 'words')
  
  df = tibble(df)
  
  df = df %>% arrange(desc(count)) %>% 
    filter(!words %in% c('to', 'of', 'activity','process', 'metabolic', 
                         'regulation', 'cellular', 'cell', 'compound',
                         'in','to','from','g1/s','in.')) %>% 
    head(elms) %>% 
    mutate(count = count/sum(count))
  
  return(df)
}


down = terms_in_words(micit_1mM_5FU_enrich_down$description) %>% mutate(direction = 'down')
up = terms_in_words(micit_1mM_5FU_enrich_up$description) %>% mutate(direction = 'up')


radar_plot = function(down,up){
  test_radar = down %>% bind_rows(up)
  test_radar = test_radar %>% 
    pivot_wider(names_from = words, values_from = count) %>% 
    replace(is.na(.), 0) %>% 
    data.frame
  
  rownames(test_radar) = test_radar[,1]
  test_radar[,1] = NULL
  
  test_radar['Max',] = max(test_radar)
  test_radar['Min',] = 0
  
  test_radar = test_radar[c(3,4,1,2),]
  
  radarchart(test_radar)
  op <- par(mar = c(1, 2, 2, 1))
  create_beautiful_radarchart(test_radar, caxislabels = c(0, 0.25, 0.5, 0.75, 1),
                              color = c("#00AFBB", "#E7B800"))
  legend(
    x = "bottom", legend = rownames(test_radar[-c(1,2),]), horiz = TRUE,
    bty = "n", pch = 20 , col = c("#00AFBB", "#E7B800"),
    text.col = "black", cex = 1, pt.cex = 1.5
  )
  par(op)
  
}


create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}


### PLOT RESULTS
# 1mM micit + 5FU
down = terms_in_words(micit_1mM_5FU_enrich_down$description, elms = 15) %>% mutate(direction = 'down')
up = terms_in_words(micit_1mM_5FU_enrich_up$description, elms = 15) %>% mutate(direction = 'up')

radar_plot(down,up)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'micit1mM_5FU_radarplot_process.pdf'),
             width = 8, height = 8, useDingbats = FALSE)

# 10mM micit + 5FU
down = terms_in_words(micit_10mM_5FU_enrich_down$description, elms = 15) %>% mutate(direction = 'down')
up = terms_in_words(micit_10mM_5FU_enrich_up$description, elms = 15) %>% mutate(direction = 'up')

radar_plot(down,up)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'micit10mM_5FU_radarplot_process.pdf'),
             width = 8, height = 8, useDingbats = FALSE)

# 10mM micit
down = terms_in_words(micit_10mM_enrich_down$description, elms = 15) %>% mutate(direction = 'down')
up = terms_in_words(micit_10mM_enrich_up$description, elms = 15) %>% mutate(direction = 'up')

radar_plot(down,up)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'micit10mM_radarplot_process.pdf'),
             width = 8, height = 8, useDingbats = FALSE)


# 5FU
down = terms_in_words(FU_enrich_down$description, elms = 15) %>% mutate(direction = 'down')
up = terms_in_words(FU_enrich_up$description, elms = 15) %>% mutate(direction = 'up')

radar_plot(down,up)

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', '5FU_radarplot_process.pdf'),
             width = 8, height = 8, useDingbats = FALSE)





# t-test with interactions ------------------------------------------------


# first do a test

test = data_long %>% 
  filter(Gene_names == 'RPL5') %>% 
  mutate(Sample = factor(Sample, levels = c('Control', 
                                            '5FU',
                                            '1mM_Micit',
                                            '10mM_Micit',
                                            '5FU_1mM_Micit', 
                                            '5FU_10mM_Micit')))

model1 = lm(Intensity ~ 0 +Sample, data = test)
summary(model1)

test %>% 
  group_by(Sample) %>% 
  summarise(Mean = mean(Intensity))

means = test %>% 
  group_by(Sample) %>% 
  summarise(Mean = mean(Intensity))

K = matrix(c(-1,   # Control
             1,  # 5FU
             0,  # 1mM_Micit
             0,  # 10mM_Micit
             0,  # 5FU_1mM_Micit
             0),  # 5FU_10mM_Micit
           1)


t <- glht(model1, linfct = mcp(Sample = K), test = adjusted('none'))
tidy(summary(t, test = adjusted('none')))


# plot examples of genes
gene = c('RRM2')
data_long %>% filter(Gene_names %in%  gene) %>% 
  mutate(Sample = factor(Sample, levels = c('Control', 
                                            '5FU',
                                            '1mM_Micit',
                                            '10mM_Micit',
                                            '5FU_1mM_Micit', 
                                            '5FU_10mM_Micit'))) %>% 
  # filter(Sample != 'Control') %>% 
  ggplot(aes(x = Sample, y = Intensity, fill = Sample)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  labs(title = gene) +
  facet_wrap(~Gene_names)




### start building our test
# first let's remove things that have n <= 2

data_long %>% 
  group_by(Gene_names) %>% 
  drop_na(Intensity, Gene_names) %>% 
  count(Sample) %>% 
  arrange(desc(n)) %>% 
  ggplot(aes(n)) +
  geom_histogram()



## FILTER LIST
# get the list of genes with only 1 datapoint or less in Control 

removals = data_long %>% 
  group_by(Gene_names) %>% 
  drop_na(Intensity, Gene_names) %>% 
  count(Sample) %>% 
  filter(n < 2) %>% 
  distinct(Gene_names) %>% t %>%  as.character

# removing NA groups
removals_na = data_long %>% group_by(Gene_names, Sample) %>% 
  summarise(na_count = sum(is.na(Intensity))) %>% 
  arrange(desc(na_count)) %>% 
  filter(na_count > 3) %>% distinct(Gene_names) %>% t %>% as.character

# final list of removals
removals = unique(c(removals,removals_na))

## read contrasts
library(readxl)
Contrasts = read_excel("Contrasts.xlsx")

contr_mat = as.matrix(Contrasts[,5:10])

contr_factors = Contrasts[,1:4]


row.names(contr_mat) = Contrasts$contrast


statsR_raw = data_long %>% 
  mutate(Sample = factor(Sample, levels = c('Control', 
                                            '5FU',
                                            '1mM_Micit',
                                            '10mM_Micit',
                                            '5FU_1mM_Micit', 
                                            '5FU_10mM_Micit'))) %>% 
  filter(!(Gene_names %in% removals)) %>% 
  drop_na(Intensity) %>% 
  group_by(Gene_names) %>% 
  nest() %>% 
  mutate(model = map(data, lm, formula = 'Intensity ~ 0 + Sample'),
         inter = map(model, glht, linfct = contr_mat, test = adjusted('none')),
         inter_sum = map(inter, summary, test = adjusted('none')),       ## this way I can extract raw p-values
         inter_tidy = map(inter_sum, tidy))



# summary(statsR_raw$inter[[1]])
# statsR_raw$inter_sum[[1]]



statsR = statsR_raw %>% 
  select(Gene_names, inter_tidy) %>% 
  unnest(cols = c(inter_tidy)) %>% 
  mutate(p.stars = gtools::stars.pval(p.value)) %>% 
  group_by(contrast) %>% 
  mutate(FDR = p.adjust(p.value, method = 'fdr')) %>% 
  ungroup %>% 
  mutate(FDR_stars = gtools::stars.pval(FDR))


statsR = statsR %>% left_join(contr_factors) %>% 
  select(Gene_names, Contrast_type:Target,contrast,everything())

write.xlsx(statsR, here('summary','Stats_R_interactions.xlsx'))





### protein groups from R stats ####

# this xlsx will be fed to the python script to extract categories and so on
# save protein groups
FU_UP = statsR %>% filter(contrast == '5FU-Ctrl', FDR < 0.05, estimate > 0) %>% 
  select(genes = Gene_names) %>% separate_rows(genes, sep = ';')
FU_DOWN = statsR %>% filter(contrast == '5FU-Ctrl', FDR < 0.05, estimate < 0)%>% 
  select(genes = Gene_names) %>% separate_rows(genes, sep = ';')

micit10_UP = statsR %>% filter(contrast == '10mMmicit - Ctrl', FDR < 0.05, estimate > 0)%>% 
  select(genes = Gene_names) %>% separate_rows(genes, sep = ';')
micit10_DOWN = statsR %>% filter(contrast == '10mMmicit - Ctrl', FDR < 0.05, estimate < 0)%>% 
  select(genes = Gene_names)%>% separate_rows(genes, sep = ';')

micit1_5FU_UP = statsR %>% filter(contrast == '5FU_1mM_Micit - Ctrl', FDR < 0.05, estimate > 0)%>% 
  select(genes = Gene_names) %>% separate_rows(genes, sep = ';')
micit1_5FU_DOWN = statsR %>% filter(contrast == '5FU_1mM_Micit - Ctrl', FDR < 0.05, estimate < 0)%>% 
  select(genes = Gene_names)%>% separate_rows(genes, sep = ';')

micit10_5FU_UP = statsR %>% filter(contrast == '5FU_10mM_Micit - Ctrl', FDR < 0.05, estimate > 0)%>% 
  select(genes = Gene_names) %>% separate_rows(genes, sep = ';')
micit10_5FU_DOWN = statsR %>% filter(contrast == '5FU_10mM_Micit - Ctrl', FDR < 0.05, estimate < 0)%>% 
  select(genes = Gene_names) %>% separate_rows(genes, sep = ';')

list_of_tables = list(
  FU_UP = FU_UP,
  FU_DOWN = FU_DOWN,
  micit10_UP = micit10_UP,
  micit10_DOWN = micit10_DOWN,
  micit1_5FU_UP = micit1_5FU_UP,
  micit1_5FU_DOWN = micit1_5FU_DOWN,
  micit10_5FU_UP = micit10_5FU_UP,
  micit10_5FU_DOWN = micit10_5FU_DOWN
)

write.xlsx(list_of_tables, here('summary','Protein_sets.xlsx'))



# exploratory plots -------------------------------------------------------



# plot examples of genes
gene = c('MYC')
data_long %>% filter(Gene_names %in%  gene) %>% 
  mutate(Sample = factor(Sample, levels = c('Control', 
                                            '5FU',
                                            '1mM_Micit',
                                            '10mM_Micit',
                                            '5FU_1mM_Micit', 
                                            '5FU_10mM_Micit'))) %>% 
  # filter(Sample != 'Control') %>% 
  ggplot(aes(x = Sample, y = Intensity, fill = Sample)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) +
  # labs(title = gene) +
  facet_wrap(~Gene_names) +
  theme(axis.text.x = element_text(hjust = 1, angle=45))




gene = c('PTPRF','PDCD4','TYMS','RRM2','KIAA0101','CDKN1A')
for (gen in gene){
  data_long %>% filter(Gene_names %in%  gen) %>% 
    mutate(Sample = factor(Sample, levels = c('Control', 
                                              '5FU',
                                              '1mM_Micit',
                                              '10mM_Micit',
                                              '5FU_1mM_Micit', 
                                              '5FU_10mM_Micit'))) %>% 
    # filter(Sample != 'Control') %>% 
    ggplot(aes(x = Sample, y = Intensity, fill = Sample)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge()) +
    # labs(title = gene) +
    facet_wrap(~Gene_names) +
    theme(axis.text.x = element_text(hjust = 1, angle=45))
  ggsave(here('summary',paste0(gen,'_boxplot.pdf')), height = 6, width = 8)
}




# enrichment exploration --------------------------------------------------

library(readxl)
FU_enrich = read_excel("summary/Enrichment_RStats/FU/FU_output.xlsx", 
                        sheet = "Process") %>% 
  select(-`...1`)


FU_enrich %>%  
  mutate(strength = -log10(number_of_genes/number_of_genes_in_background)) %>% 
  filter(fdr <= 0.01)










# Heatmaps ----------------------------------------------------------------

# this specifies the groups to filter in the Rstats
contr_groups = c('5FU - 1mM_Micit', '5FU-Ctrl','10mMmicit - Ctrl','1mMmicit - Ctrl',
                 '5FU_1mM_Micit - Ctrl', '5FU_10mM_Micit - Ctrl')

### top 100 ####


large_effect_prots = statsR %>% 
  mutate(estimate = abs(estimate)) %>% 
  arrange(desc(estimate)) %>% 
  filter(FDR < 0.05) %>% 
  distinct(Gene_names, .keep_all = TRUE) %>% 
  head(100) %>% 
  separate_rows(Gene_names, sep = ';') %>% 
  select(Gene_names) %>% t %>% as.character()

write.table(large_effect_prots, 'large_effect_prots.txt',
            quote = F, row.names = F, col.names = F)

# KEGG pathways
data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% large_effect_prots) %>% 
  distinct(Gene_names,.keep_all = T) %>% 
  select(Gene_names, KEGG_name)  %>% 
  # separate(KEGG_name,sep=';',into=)
  separate_rows(KEGG_name, sep = ';') %>% 
  filter(!(KEGG_name %in% c('Phototransduction - fly','Cell cycle - yeast',
                          'Meiosis - yeast','Amoebiasis','Toxoplasmosis',
                          'Vascular smooth muscle contraction',
                          'Proximal tubule bicarbonate reclamation','Malaria'))) %>% 
  count(Gene_names) %>% arrange(desc(n))


heat_zscore_wide = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% large_effect_prots) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)

pval_top = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% large_effect_prots, contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)

# generate matrices
heat_mat = as.matrix(heat_zscore_wide[,2:7])
row.names(heat_mat) = heat_zscore_wide$Gene_names

pval_mat = as.matrix(pval_top[,2:7])
row.names(pval_mat) = pval_top$Gene_names


Heatmap(heat_mat,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        name = "z-score",
        column_km = 3,
        # row_km = 3,
        row_names_gp = gpar(fontsize = 5),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 4))


dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_top100.pdf'),
             width = 5, height = 9, useDingbats = FALSE)



### global heatmap ####


heat_zscore_wide_all = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  # filter(Gene_names %in% large_effect_prots) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)




heat_mat_all = as.matrix(heat_zscore_wide_all[,2:7])

row.names(heat_mat_all) = heat_zscore_wide_all$Gene_names


Heatmap(heat_mat_all,
        name = "z-score",
        column_km = 3,
        # row_km = 3,
        na_col = "black",
        row_names_gp = gpar(fontsize = 4),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 10),
        show_row_names = FALSE)


dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_total.pdf'),
             width = 5, height = 6, useDingbats = FALSE)




### save dataset for comparison ####


comparison_df = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  pivot_wider(names_from = Sample, values_from = Mean)


colnames(comparison_df) = c("Gene_names","Micit_10",
                                   "Micit_1","FU","FU_10mM_Micit",
                                   "FU_1mM_Micit","Control")


write_csv(comparison_df,here('summary', 'means_samples.csv'))
write_csv(comparison_df,here('dataset_comparison', 'means_samples.csv'))




### TCA genes ####
tca = data %>% 
  filter(str_detect(KEGG_name,'TCA') ) %>% 
  select(Gene_names) %>% 
  t %>% as.character()


# calculate means and zscore and put it wider
heat_tca = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% tca) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)

pval_top = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% tca, 
         contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)


heat_mat = as.matrix(heat_tca[,2:7])
row.names(heat_mat) = heat_tca$Gene_names

pval_mat = as.matrix(pval_top[,2:7])
row.names(pval_mat) = pval_top$Gene_names


Heatmap(heat_mat,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        name = "z-score",
        column_km = 3,
        # row_km = 3,
        row_names_gp = gpar(fontsize = 4),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 8))

dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_TCA.pdf'),
             width = 5, height = 7, useDingbats = FALSE)



### mTOR genes ####

mtor = data %>% 
  filter(str_detect(KEGG_name,'mTOR') | str_detect(GOBP_name, 'mTOR') |  
           str_detect(GOMF_name, 'mTOR')) %>% 
  select(Gene_names) %>% 
  t %>% as.character()



# calculate means and zscore and put it wider
heat_mtor = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% mtor) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)

pval_mtor = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% mtor, contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)
  
# generate matrices
heat_mat = as.matrix(heat_mtor[,2:7])
row.names(heat_mat) = heat_mtor$Gene_names

pval_mat = as.matrix(pval_mtor[,2:7])
row.names(pval_mat) = pval_mtor$Gene_names


Heatmap(heat_mat, 
        name = "z-score",
        column_km = 3,
        # row_km = 3,
        row_names_gp = gpar(fontsize = 4),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 7),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })

dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_mTOR.pdf'),
             width = 5, height = 6, useDingbats = FALSE)







### p53 genes ####

p53 = data %>% 
  filter(str_detect(KEGG_name,'p53') | str_detect(GOBP_name, 'p53') |  
           str_detect(GOMF_name, 'p53')) %>% 
  select(Gene_names) %>% 
  t %>% as.character()



# calculate means and zscore and put it wider
heat = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% p53) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)

pval_mtor = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% p53, contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)

# generate matrices
heat_mat = as.matrix(heat[,2:7])
row.names(heat_mat) = heat$Gene_names

pval_mat = as.matrix(pval_mtor[,2:7])
row.names(pval_mat) = pval_mtor$Gene_names


Heatmap(heat_mat, 
        name = "z-score",
        column_km = 3,
        # row_km = 3,
        row_names_gp = gpar(fontsize = 4),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 6),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 7))
        })

dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_p53.pdf'),
             width = 5, height = 9, useDingbats = FALSE)





### cell cycle genes ####

prots = data %>% 
  filter(str_detect(KEGG_name,'cell cycle') | str_detect(GOBP_name, 'cell cycle') |  
           str_detect(GOMF_name, 'cell cycle')) %>% 
  select(Gene_names) %>% 
  t %>% as.character()



# calculate means and zscore and put it wider
heat = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)

pval_mtor = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots, contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)

# generate matrices
heat_mat = as.matrix(heat[,2:7])
row.names(heat_mat) = heat$Gene_names

pval_mat = as.matrix(pval_mtor[,2:7])
row.names(pval_mat) = pval_mtor$Gene_names


Heatmap(heat_mat, 
        name = "z-score",
        column_km = 3,
        # row_km = 3,
        row_names_gp = gpar(fontsize = 4),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 3),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 7))
        })

dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_cell_cycle.pdf'),
             width = 5, height = 20, useDingbats = FALSE)





### mitochondrion genes ####

prots = data %>% 
  filter(str_detect(KEGG_name,'mitochondrion') | str_detect(GOBP_name, 'mitochondrion') |  
           str_detect(GOMF_name, 'mitochondrion')) %>% 
  select(Gene_names) %>% 
  t %>% as.character()



# calculate means and zscore and put it wider
heat = data_long %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots) %>% 
  group_by(Gene_names, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>% 
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score)

pval_mtor = statsR %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots, contrast %in% contr_groups) %>% 
  # create dummy variable for Control
  mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
                               TRUE ~ FDR_stars),
         Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
                            TRUE ~ Target)) %>% 
  select(Gene_names,Target,FDR_stars) %>% 
  pivot_wider(names_from = Target, values_from = FDR_stars)

# generate matrices
heat_mat = as.matrix(heat[,2:7])
row.names(heat_mat) = heat$Gene_names

pval_mat = as.matrix(pval_mtor[,2:7])
row.names(pval_mat) = pval_mtor$Gene_names


Heatmap(heat_mat, 
        name = "z-score",
        # column_km = 3,
        # row_km = 3,
        row_names_gp = gpar(fontsize = 4),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 7),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", pval_mat[i, j]), x, y, gp = gpar(fontsize = 7))
        })

dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_mitochondrion.pdf'),
             width = 5, height = 20, useDingbats = FALSE)









### category selection genes ####

kegg_select = c('Cell cycle','p53 signaling pathway', 'Citrate cycle (TCA cycle)',
                'mTOR signaling pathway', 'RNA transport', 'Purine metabolism',
                'Pyrimidine metabolism', 'Ribosome', 'Fatty acid metabolism',
                'Focal adhesion', 'Glycolysis','Gap junction','ErbB signaling pathway',
                'MAPK signaling pathway')


prots = data %>% 
  separate_rows(KEGG_name, sep = ';') %>% 
  filter(KEGG_name %in% kegg_select) %>% 
  # select(Gene_names) %>% 
  distinct(Gene_names) %>%
  t %>% as.character()


# cluster = new_cluster(8)
# calculate means and zscore and put it wider
heat = data_long %>% 
  separate_rows(KEGG_name, sep=';') %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots) %>% 
  arrange(Gene_names) %>% 
  filter(KEGG_name %in%  kegg_select) %>% 
  unite(ID, Gene_names, KEGG_name, remove = FALSE) %>% 
  group_by(ID, Sample) %>% 
  # partition(cluster) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  # collect() %>%
  mutate(z_score = scale(Mean)) %>%
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score) %>% 
  separate(ID, into = c('Gene_names', 'Pathway'), sep = '_') %>% 
  arrange(Pathway)

# pval_mtor = statsR %>% 
#   filter(!(Gene_names %in% removals)) %>% 
#   filter(Gene_names %in% prots, contrast %in% contr_groups) %>% 
#   # create dummy variable for Control
#   mutate(FDR_stars = case_when(contrast == '5FU - 1mM_Micit' ~ '',
#                                TRUE ~ FDR_stars),
#          Target = case_when(contrast == '5FU - 1mM_Micit' ~ 'Control',
#                             TRUE ~ Target)) %>% 
#   select(Gene_names,Target,FDR_stars) %>% 
#   pivot_wider(names_from = Target, values_from = FDR_stars)

# generate matrices
heat_mat = as.matrix(heat[,3:8])
row.names(heat_mat) = heat$Gene_names

# pval_mat = as.matrix(pval_mtor[,2:7])
# row.names(pval_mat) = pval_mtor$Gene_names

library(circlize)
col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

ha = rowAnnotation(Pathway = heat$Pathway)

h1 = Heatmap(heat_mat, 
        name = "z-score",
        # column_km = 3,
        # row_km = 3,
        cluster_rows = FALSE, # turns off row clustering
        show_row_names = FALSE,
        row_names_gp = gpar(fontsize = 4),
        right_annotation = ha,
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 7))

heat$Pathway

# plot every pathway as a subheatmap
# LONG AND UGLY AND LONG AND UGLY CODE
# BUT IT WORKS, OK? 

cell_ht = Heatmap(heat_mat[1:76,], 
        name = "Z-score",
        show_row_names = FALSE,
        row_names_gp = gpar(fontsize = 4),
        right_annotation = rowAnnotation(
          Pathway = heat$Pathway[1:76],
          show_annotation_name = F),
        column_names_rot =30, 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 7))

tca_ht = Heatmap(heat_mat[77:103,], 
                  name = "Z-score",
                  show_row_names = FALSE,
                  row_names_gp = gpar(fontsize = 4),
                  right_annotation = rowAnnotation(
                    Pathway = heat$Pathway[77:103],
                    show_annotation_name = F),
                  column_names_rot =30, 
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 7))

erbB_ht = Heatmap(heat_mat[104:139,], 
                 name = "Z-score",
                 show_row_names = FALSE,
                 row_names_gp = gpar(fontsize = 4),
                 right_annotation = rowAnnotation(
                   Pathway = heat$Pathway[104:139],
                   show_annotation_name = F),
                 column_names_rot =30, 
                 column_names_side = "top",
                 column_names_gp = gpar(fontsize = 7))

fa_ht = Heatmap(heat_mat[140:168,], 
                  name = "Z-score",
                  show_row_names = FALSE,
                  row_names_gp = gpar(fontsize = 4),
                  right_annotation = rowAnnotation(
                    Pathway = heat$Pathway[140:168],
                    show_annotation_name = F),
                  column_names_rot =30, 
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 7))

focal_ht = Heatmap(heat_mat[169:243,], 
                name = "Z-score",
                show_row_names = FALSE,
                row_names_gp = gpar(fontsize = 4),
                right_annotation = rowAnnotation(
                  Pathway = heat$Pathway[169:243],
                  show_annotation_name = F),
                column_names_rot =30, 
                column_names_side = "top",
                column_names_gp = gpar(fontsize = 7))

gap_ht = Heatmap(heat_mat[244:280,], 
                   name = "Z-score",
                   show_row_names = FALSE,
                   row_names_gp = gpar(fontsize = 4),
                   right_annotation = rowAnnotation(
                     Pathway = heat$Pathway[244:280],
                     show_annotation_name = F),
                   column_names_rot =30, 
                   column_names_side = "top",
                   column_names_gp = gpar(fontsize = 7))

mapk_ht = Heatmap(heat_mat[281:361,], 
                 name = "Z-score",
                 show_row_names = FALSE,
                 row_names_gp = gpar(fontsize = 4),
                 right_annotation = rowAnnotation(
                   Pathway = heat$Pathway[281:361],
                   show_annotation_name = F),
                 column_names_rot =30, 
                 column_names_side = "top",
                 column_names_gp = gpar(fontsize = 7))

mtor_ht = Heatmap(heat_mat[362:379,], 
                  name = "Z-score",
                  show_row_names = FALSE,
                  row_names_gp = gpar(fontsize = 4),
                  right_annotation = rowAnnotation(
                    Pathway = heat$Pathway[362:379],
                    show_annotation_name = F),
                  column_names_rot =30, 
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 7))

p53_ht = Heatmap(heat_mat[380:405,], 
                  name = "Z-score",
                  show_row_names = FALSE,
                  row_names_gp = gpar(fontsize = 4),
                  right_annotation = rowAnnotation(
                    Pathway = heat$Pathway[380:405],
                    show_annotation_name = F),
                  column_names_rot =30, 
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 7))

purine_ht = Heatmap(heat_mat[406:487,], 
                 name = "Z-score",
                 show_row_names = FALSE,
                 row_names_gp = gpar(fontsize = 4),
                 right_annotation = rowAnnotation(
                   Pathway = heat$Pathway[406:487],
                   show_annotation_name = F),
                 column_names_rot =30, 
                 column_names_side = "top",
                 column_names_gp = gpar(fontsize = 7))

pyr_ht = Heatmap(heat_mat[488:555,], 
                    name = "Z-score",
                    show_row_names = FALSE,
                    row_names_gp = gpar(fontsize = 4),
                    right_annotation = rowAnnotation(
                      Pathway = heat$Pathway[488:555],
                      show_annotation_name = F),
                    column_names_rot =30, 
                    column_names_side = "top",
                    column_names_gp = gpar(fontsize = 7))

ribosome_ht  = Heatmap(heat_mat[556:635,], 
                 name = "Z-score",
                 show_row_names = FALSE,
                 row_names_gp = gpar(fontsize = 4),
                 right_annotation = rowAnnotation(
                   Pathway = heat$Pathway[556:635],
                   show_annotation_name = F),
                 column_names_rot =30, 
                 column_names_side = "top",
                 column_names_gp = gpar(fontsize = 7))

rnatransport_ht = Heatmap(heat_mat[636:757,], 
                      name = "Z-score",
                      show_row_names = FALSE,
                      row_names_gp = gpar(fontsize = 4),
                      right_annotation = rowAnnotation(
                        Pathway = heat$Pathway[636:757],
                        show_annotation_name = F),
                      column_names_rot =30, 
                      column_names_side = "top",
                      column_names_gp = gpar(fontsize = 7))



ht_list = cell_ht %v% tca_ht %v% erbB_ht %v% 
  fa_ht %v% focal_ht %v% gap_ht %v% mapk_ht %v% 
  mtor_ht %v% p53_ht %v% purine_ht %v% pyr_ht %v%  
  ribosome_ht %v% rnatransport_ht

lgd = Legend(labels = unique(heat$Pathway), title = "", 
             legend_gp = gpar(fill = 1:13),
             grid_height = unit(0, "cm"), grid_width = unit(0, "mm"),
             labels_gp = gpar(col = "white", fontsize = 0))

draw(ht_list, annotation_legend_list = lgd) 



dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_KEGG.pdf'),
             width = 8, height = 14, useDingbats = FALSE)






###  sig category selection genes ####

prots = data %>% 
  separate_rows(KEGG_name, sep = ';') %>% 
  filter(KEGG_name %in% kegg_select) %>% 
  # select(Gene_names) %>% 
  distinct(Gene_names) %>%
  t %>% as.character()

# take the proteins and see which ones are at least significant
# in one of the groups in respect to the control
prots_sig = statsR %>%
  filter(!(Gene_names %in% removals)) %>%
  filter(Gene_names %in% prots, contrast %in% contr_groups) %>% 
  filter(FDR<=0.05) %>% 
  distinct(Gene_names) %>% t %>% as.character
  


# calculate means and zscore and put it wider
heat = data_long %>% 
  separate_rows(KEGG_name, sep=';') %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots_sig) %>% 
  arrange(Gene_names) %>% 
  filter(KEGG_name %in%  kegg_select) %>% 
  unite(ID, Gene_names, KEGG_name, remove = FALSE) %>% 
  group_by(ID, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>%
  select(-Mean) %>% 
  pivot_wider(names_from = Sample, values_from = z_score) %>% 
  separate(ID, into = c('Gene_names', 'Pathway'), sep = '_') %>% 
  arrange(Pathway)



# generate matrices
heat_mat = as.matrix(heat[,3:8])
row.names(heat_mat) = heat$Gene_names

# pval_mat = as.matrix(pval_mtor[,2:7])
# row.names(pval_mat) = pval_mtor$Gene_names


ha = rowAnnotation(Pathway = heat$Pathway)

h1 = Heatmap(heat_mat, 
             name = "z-score",
             # column_km = 3,
             # row_km = 3,
             cluster_rows = FALSE, # turns off row clustering
             show_row_names = FALSE,
             row_names_gp = gpar(fontsize = 4),
             right_annotation = ha,
             column_names_rot =30, 
             column_names_side = "top",
             column_names_gp = gpar(fontsize = 7))

h1


# get indexes to subset 
require(rlist)
cat_list = list()
for (elm in unique(heat$Pathway)){
  min_index = min(str_which(heat$Pathway,as.character(elm)))
  max_index = max(str_which(heat$Pathway,as.character(elm)))
  cat_list = list.append(cat_list, elm = c(min_index,max_index))
}
# fix a couple things
names(cat_list) = unique(heat$Pathway)
cat_list$`Citrate cycle (TCA cycle)` = c(39,64)


# this loop will create many variables for each pathway
# each one of them will be a heatmap
paths = unique(heat$Pathway)
list_of_vars = c()
for (i in 1:length(paths)){
  var_name = paste0(str_replace_all(paths[i],' ','_'),'_ht')
  print(var_name)
  list_of_vars = c(var_name,list_of_vars)
  
  min_index = cat_list[[i]][1]
  max_index = cat_list[[i]][2]
  print(min_index)
  print(max_index)
  
  var_name_heat = Heatmap(heat_mat[min_index:max_index,], 
                    name = "Z-score",
                    show_row_names = FALSE,
                    row_names_gp = gpar(fontsize = 4),
                    right_annotation = rowAnnotation(
                      Pathway = heat$Pathway[min_index:max_index],
                      show_annotation_name = F),
                    column_names_rot =30, 
                    column_names_side = "top",
                    column_names_gp = gpar(fontsize = 7))
  
  assign(var_name, var_name_heat)
  
}

dummy = c()
for (i in list_of_vars){
  cosa = paste(i,'%v% ')
  dummy = c(cosa,dummy)
}
str_c(dummy, collapse = '')

ht_list = Cell_cycle_ht %v% `Citrate_cycle_(TCA_cycle)_ht` %v% 
  ErbB_signaling_pathway_ht %v% Fatty_acid_metabolism_ht %v% 
  Focal_adhesion_ht %v% Gap_junction_ht %v% MAPK_signaling_pathway_ht %v%
  mTOR_signaling_pathway_ht %v% 
  p53_signaling_pathway_ht %v% Purine_metabolism_ht %v% 
  Pyrimidine_metabolism_ht %v% Ribosome_ht %v% RNA_transport_ht

lgd = Legend(labels = unique(heat$Pathway), title = "", 
             legend_gp = gpar(fill = 1:13),
             grid_height = unit(0, "cm"), grid_width = unit(0, "mm"),
             labels_gp = gpar(col = "white", fontsize = 0))

draw(ht_list, annotation_legend_list = lgd) 



dev.copy2pdf(device = cairo_pdf,
             file = here('summary/heatmaps', 'heatmap_KEGG_sigProts.pdf'),
             width = 8, height = 14, useDingbats = FALSE)






# radar plot zscores --------------------------------------------------------------

### barplots of selected pathways ####

zscores_paths = data_long %>% 
  separate_rows(KEGG_name, sep=';') %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots_sig) %>% 
  arrange(Gene_names) %>% 
  filter(KEGG_name %in%  kegg_select) %>% 
  unite(ID, Gene_names, KEGG_name, remove = FALSE) %>% 
  group_by(ID, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>%
  select(-Mean) %>% 
  separate(ID, into = c('Gene_names', 'Pathway'), sep = '_') %>% 
  arrange(Pathway) %>% 
  group_by(Pathway, Sample) %>% 
  summarise(zMean = mean(z_score,na.rm = TRUE),
            zSD = sd(z_score, na.rm = TRUE),
            zSEM = zSD/sqrt(n())) %>% 
  ungroup




zscores_paths %>% 
  ungroup %>% 
  mutate(Sample = factor(Sample, 
                         levels = c('Control', 
                                    '5FU',
                                    '1mM_Micit',
                                    '5FU_1mM_Micit',
                                    '10mM_Micit',
                                    '5FU_10mM_Micit'))) %>% 
  ggplot(aes(Pathway, zMean, fill = Sample)) +
  geom_histogram(stat='identity', position = 'dodge2') +
  geom_errorbar(aes(ymin = zMean - zSEM, ymax = zMean + zSEM), 
                position = position_dodge(0.9), width = 0) +
  scale_fill_manual(values = c('grey50',
                               '#FA2713', # 5FU
                               '#F0C600', # 1 micit
                               '#E05C10', # 5fu + 1 micit
                               '#F79E05', # 10 micit
                               '#AD4413'  # 5fu + 10 micit
  )) +
  labs(
    x = 'KEGG Pathway',
    y = 'z-score average (+- SEM)',
    fill = 'Condition'
  ) +
  facet_wrap(~Pathway, scales = 'free')


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'zScore_means_sigProts.pdf'),
             width = 10, height = 9, useDingbats = FALSE)




### calculate zscore mean for everything ####

# take the proteins and see which ones are at least significant
# in one of the groups in respect to the control
prots_sig = statsR %>%
  filter(!(Gene_names %in% removals)) %>%
  filter(contrast %in% contr_groups) %>% 
  filter(FDR<=0.05) %>% 
  distinct(Gene_names) %>% t %>% as.character


zscores_paths_all = data_long %>% 
  separate_rows(KEGG_name, sep=';') %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots_sig) %>% 
  arrange(Gene_names) %>% 
  # filter(KEGG_name %in%  kegg_select) %>% 
  unite(ID, Gene_names, KEGG_name, remove = FALSE) %>% 
  group_by(ID, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  mutate(z_score = scale(Mean)) %>%
  select(-Mean) %>% 
  separate(ID, into = c('Gene_names', 'Pathway'), sep = '_') %>% 
  arrange(Pathway) %>% 
  group_by(Pathway, Sample) %>% 
  summarise(zMean = mean(z_score,na.rm = TRUE),
            zSD = sd(z_score, na.rm = TRUE),
            zSEM = zSD/sqrt(n())) %>% 
  ungroup





zscores_paths_all %>% 
  ungroup %>% 
  mutate(Sample = factor(Sample, 
                         levels = c('Control', 
                                    '5FU', 
                                    '1mM_Micit', '10mM_Micit',
                                    '5FU_1mM_Micit', '5FU_10mM_Micit'))) %>% 
  ggplot(aes(Pathway, zMean, fill = Sample)) +
  geom_histogram(stat='identity', position = 'dodge2') +
  geom_errorbar(aes(ymin = zMean - zSEM, ymax = zMean + zSEM), 
                position = position_dodge(0.9), width = 0) +
  scale_fill_manual(values = c('grey50',
                               '#FA2713', # 5FU
                               '#F0C600', # 1 micit
                               '#F79E05', # 10 micit
                               '#E05C10', # 5fu + 1 micit
                               '#AD4413'  # 5fu + 10 micit
  )) +
  labs(
    x = 'KEGG Pathway',
    y = 'z-score average (+- SEM)',
    fill = 'Condition'
  ) +
  facet_wrap(~Pathway, scales = 'free')


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'zScore_means_sigProts_ALL.pdf'),
             width = 70, height = 60, useDingbats = FALSE)







## create the spider plot #####


mins = c(min(zscores_paths$zMean))
maxs = c(max(zscores_paths$zMean))

# define values that are symmetric
mins = -1.5
maxs = 1.5

zscores_wide = zscores_paths %>% 
  select(Pathway:zMean) %>% 
  pivot_wider(names_from = Sample, values_from = zMean) %>%
  mutate(Max = maxs, Min = mins, Limit = 0, .before = `10mM_Micit`)

zscores_mat = as.matrix(zscores_wide[,2:10])

rownames(zscores_mat) = zscores_wide$Pathway

zscores_mat = t(zscores_mat)



zscores_mat = as.data.frame(zscores_mat[1:9,])


#~~~~~~~~~~~~~~~~~

## This piece of code will plot all Cell samples
# per separate in a single plot, might be good 
# for comparison

opar = par() 
# Define settings for plotting in a 3x4 grid, with appropriate margins:
par(mar = rep(0.8,4))
par(mfrow = c(3,2))
# Produce a radar-chart for each student
for (i in 4:nrow(zscores_mat)) {
  radarchart(
    zscores_mat[c(1:3, i), ],
    pfcol = c("#99999980",NA),
    pcol= c(NA,2), plty = 1, plwd = 2,
    title = row.names(zscores_mat)[i]
  )
}
# Restore the standard par() settings
par <- par(opar) 


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'radar_ALL_samples.pdf'),
             width = 15, height = 10, useDingbats = FALSE)


#~~~~~~~~~~~~

select_zscores = zscores_mat[c(1:3,4,6,7),] 


opar = par() 
# Define settings for plotting in a 3x4 grid, with appropriate margins:
par(mar = rep(1,4))

radarchart(
  select_zscores,
  caxislabels = c(-1.5, -0.75, 0, 0.75, 1.5),
  pfcol = c("#99999980",NA,NA,NA),
  pcol= c(NA,
          '#F2AB0C',
          '#E60B1A',
          '#1700F2'), plty = 1, plwd = 2
)

legend(
  # x = "bottom",
  x = -1,y = -1.2,
  legend = rownames(select_zscores[-c(1,2,3),]), horiz = TRUE,
  bty = "n", pch = 20 ,col= c('#F2AB0C',
                          '#E60B1A',
                          '#1700F2'),
  text.col = "black", cex = 1, pt.cex = 2.5
)

par <- par(opar) 

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'radar_5FU_Micit_comparison.pdf'),
             width = 13, height = 11.5, useDingbats = FALSE)




#~~~~~~~~~~~~~~
### ggradar version ####

library(ggradar)

zscores_wide = zscores_paths %>% 
  select(Pathway:zMean) %>% 
  pivot_wider(names_from = Sample, values_from = zMean) 

zscores_mat = as.matrix(zscores_wide[,2:7])

rownames(zscores_mat) = zscores_wide$Pathway

zscores_mat = t(zscores_mat)

zscores_mat = zscores_mat %>% as_tibble() %>% 
  mutate(Sample = unique(zscores_paths$Sample), .before = `Cell cycle`)


zscores_mat %>% 
  # filter(Pathway == 'Cell cycle') %>% 
  mutate(Sample = factor(Sample, 
                         levels = c('Control', 
                                    '5FU', 
                                    '1mM_Micit', '10mM_Micit',
                                    '5FU_1mM_Micit', '5FU_10mM_Micit'))) %>% 
  ggradar(
    values.radar = c("-1.5", "0", "1.5"),
    grid.min = -1.5, grid.mid = 0, grid.max = 1.5,
    group.line.width = 1, 
    group.colours = c('grey50',
                      '#FA2713', # 5FU
                      '#F0C600', # 1 micit
                      '#F79E05', # 10 micit
                      '#E05C10', # 5fu + 1 micit
                      '#AD4413'),
    group.point.size = 3,
    background.circle.colour = "white",
    gridline.mid.colour = "grey",
    legend.position = "bottom"
  )



names(zscores_mat)

zscores_mat %>% 
  add_row(Sample = 'DUMMY', `Cell cycle` = 0, `Citrate cycle (TCA cycle)` = 0,
          `ErbB signaling pathway` = 0, `Fatty acid metabolism` = 0,
          `Focal adhesion` = 0, `Gap junction`=0,`MAPK signaling pathway`=0,
          `mTOR signaling pathway`=0,`p53 signaling pathway`=0,
          `Purine metabolism`=0,`Pyrimidine metabolism`=0,
          Ribosome = 0, `RNA transport`=0) %>% 
  filter(Sample %in% c('DUMMY','5FU','10mM_Micit','5FU_10mM_Micit')) %>% 
  mutate(Sample = factor(Sample, 
                         levels = c('DUMMY',
                                    '5FU', 
                                    '10mM_Micit', 
                                    '5FU_10mM_Micit'
                                    ))) %>% 
  ggradar(
    values.radar = c("-1.5", "0", "1"),
    grid.min = -1.5, grid.mid = 0, grid.max = 1,
    group.line.width = 1, 
    group.colours = c('grey70',
                      '#FA2713', # 5FU
                      '#F79E05', # 10 micit
                      '#AD4413'),
    group.point.size = 3,
    background.circle.colour = "white",
    gridline.mid.colour = "grey",
    legend.position = "bottom"
  ) 

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'radar_categories_v1.pdf'),
             width = 12, height = 11, useDingbats = FALSE)







# protein means comparison ------------------------------------------------

# Here I'll do the same comparison as for the radar plots, but 
# using the ratio of Treatment over control, might be more visual
# It's important to notice that this will NOT be a different analysis
# only a different way to represent things


ratios_paths = data_long %>% 
  separate_rows(KEGG_name, sep=';') %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots_sig) %>% 
  arrange(Gene_names) %>% 
  filter(KEGG_name %in%  kegg_select) %>% 
  unite(ID, Gene_names, KEGG_name, remove = FALSE) %>% 
  group_by(ID, Sample) %>% 
  summarise(Mean = mean(Intensity, na.rm = TRUE)) %>% 
  group_by(ID) %>% 
  mutate(ratio = Mean - Mean[Sample == 'Control']) %>% 
  filter(Sample != 'Control') %>% 
  select(-Mean) %>% 
  separate(ID, into = c('Gene_names', 'Pathway'), sep = '_') %>% 
  arrange(Pathway) %>% 
  group_by(Pathway, Sample) %>% 
  summarise(Mean = mean(ratio ,na.rm = TRUE),
            SD = sd(ratio, na.rm = TRUE),
            SEM = SD/sqrt(n())) %>% 
  ungroup



ratios_paths %>% 
  ungroup %>% 
  mutate(Sample = factor(Sample, 
                         levels = c('5FU',
                                    '1mM_Micit',
                                    '10mM_Micit',
                                    '5FU_1mM_Micit', '5FU_10mM_Micit'))) %>% 
  ggplot(aes(Pathway, Mean, fill = Sample)) +
  geom_histogram(stat='identity', position = 'dodge2') +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), 
                position = position_dodge(0.9), width = 0) +
  scale_fill_manual(values = c(
                               '#FA2713', # 5FU
                               '#F0C600', # 1 micit
                               '#F79E05', # 10 micit
                               '#E05C10', # 5fu + 1 micit
                               '#AD4413'  # 5fu + 10 micit
  )) +
  labs(
    x = 'KEGG Pathway',
    y = 'Ratio agains control average (+- SEM)',
    fill = 'Condition'
  ) +
  facet_wrap(~Pathway, scales = 'free')


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'ratios_means_sigProts.pdf'),
             width = 10, height = 9, useDingbats = FALSE)






## create the spider plot #####


mins = c(min(ratios_paths$Mean))
maxs = c(max(ratios_paths$Mean))

# define values that are symmetric
mins = -0.5
maxs = 0.25

ratios_wide = ratios_paths %>% 
  select(Pathway:Mean) %>% 
  pivot_wider(names_from = Sample, values_from = Mean) %>%
  mutate(Max = maxs, Min = mins, Limit = 0, .before = `10mM_Micit`)

ratios_mat = as.matrix(ratios_wide[,2:9])

rownames(ratios_mat) = ratios_wide$Pathway

ratios_mat = t(ratios_mat)



ratios_mat = as.data.frame(ratios_mat[1:8,])


#~~~~~~~~~~~~~~~~~

## This piece of code will plot all Cell samples
# per separate in a single plot, might be good 
# for comparison

opar = par() 
# Define settings for plotting in a 3x4 grid, with appropriate margins:
par(mar = rep(0.8,4))
par(mfrow = c(3,2))
# Produce a radar-chart for each student
for (i in 4:nrow(ratios_mat)) {
  radarchart(
    ratios_mat[c(1:3, i), ],
    pfcol = c("#99999980",NA),
    pcol= c(NA,2), plty = 1, plwd = 2,
    title = row.names(ratios_mat)[i]
  )
}
# Restore the standard par() settings
par <- par(opar) 


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'radar_ratios_ALL_samples.pdf'),
             width = 15, height = 10, useDingbats = FALSE)


#~~~~~~~~~~~~

ratios_mat_select = ratios_mat[c(1:3,4,6,7),] 


opar = par() 
# Define settings for plotting in a 3x4 grid, with appropriate margins:
par(mar = rep(1,4))

radarchart(
  ratios_mat_select,
  caxislabels = c(-1.5, -0.75, 0, 0.75, 1.5),
  pfcol = c("#99999980",NA,NA,NA),
  pcol= c(NA,
          '#F2AB0C',
          '#E60B1A',
          '#1700F2'), plty = 1, plwd = 2
)

legend(
  # x = "bottom",
  x = -1,y = -1.2,
  legend = rownames(ratios_mat_select[-c(1,2,3),]), horiz = TRUE,
  bty = "n", pch = 20 ,col= c('#F2AB0C',
                              '#E60B1A',
                              '#1700F2'),
  text.col = "black", cex = 1, pt.cex = 2.5
)

par <- par(opar) 

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'radar_ratios_5FU_Micit_comparison.pdf'),
             width = 13, height = 11.5, useDingbats = FALSE)




# Random Forests ----------------------------------------------------------

library(tidymodels)


rf_data = data_long %>% 
  select(Gene_names, Sample, Intensity, KEGG_name, mol_weight) %>% 
  separate_rows(KEGG_name, sep=';') %>% 
  filter(!(Gene_names %in% removals)) %>% 
  filter(Gene_names %in% prots_sig) %>% 
  arrange(Gene_names) %>% 
  filter(KEGG_name %in%  kegg_select)

data_rec = recipe(Sample ~ ., data = rf_data) %>%
  # update_role(Gene_names, new_role = "ID") %>%
  step_dummy(all_nominal(), -all_outcomes())

data_prep = prep(data_rec)
juiced = juice(data_prep)

data_ranger =  rand_forest(trees = 100, mode = "classification") %>%
  set_engine("randomForest") %>%
  fit(Sample ~ ., data = juiced)

data_ranger$lvl




