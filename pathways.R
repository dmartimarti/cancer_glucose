
### KEGG networks and paths


# setwd("/Users/dmarti14/Documents/MRC_Postdoc/Projects/Leo")

# libraries
library(tidyverse)
library(readxl)
library(ggrepel)
library(pathview)

# my library of functions
source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')

# session options
options(width = 220)

#####################################
### WARNING!!! loading pathview library will mask a lot of different functions from other libraries and from R base
### such as: 
###   dplyr & tidyr -> combine, intersect, setdiff, union, first, rename, expand, collapse, desc, slice, select
#####################################

# for more information about pathview library, visit: https://www.bioconductor.org/packages/devel/bioc/vignettes/pathview/inst/doc/pathview.pdf

# load data

# every plot will go to Summary folder
odir <- 'Summary'

data_dev <- read_xlsx('Develop_data/Summary_Keio Sublibrary_Glucose Supp_10mM_N2 worms_08-12-18.xlsx', sheet = 'Summary') 
data_dev$X__1 = NULL

data_dev = data_dev %>%
	gather(Drug, Score, `0`, `1`, `2.5`, `5`) %>% # make table long
	unite(ID, Genotype, Supplement, Replicate, remove = FALSE) %>% 
	mutate(Replicate = as.numeric(Replicate),
		Supplement_mM = ifelse(Supplement == 'Control', 0, 10),
		Supplement_mM = as.factor(Supplement_mM),
		Supplement = as.factor(Supplement),
		Genotype = as.factor(Genotype),
		Pathway = as.factor(Pathway),
		ID = as.factor(ID),
		Drug = as.factor(Drug)) %>%
	select(Genotype, Well, ID, Pathway, Replicate, Supplement, Supplement_mM, Drug, Score)


# summary of the data

data.sum <- data_dev %>%
  group_by(Supplement, Supplement_mM, Drug, Genotype) %>%
  summarise(Median_Score = median(Score, na.rm = TRUE),
  			Mean = mean(Score, na.rm = TRUE),
  			SD = sd(Score, na.rm = TRUE)) %>%
  mutate(BW_Score = Median_Score[Genotype == 'BW'],
         BW_Mean = Mean[Genotype == 'BW']) %>%
  ungroup %>%
  group_by(Genotype, Drug) %>%
  mutate(No_supplement = ifelse("0" %in% Supplement_mM, Median_Score[Supplement_mM == '0'], NA)) %>% # Get phenotype values at Supplement_mM=0 
  ungroup %>%
  mutate(BW_norm = Median_Score - BW_Score,
         BW_Mean_norm = Mean - BW_Mean,
         Suppl_norm = Median_Score - No_supplement,
         Interaction = Median_Score - BW_Score - No_supplement,
         Supplement = 'Glucose')


# generate two different test datasets 


# read ID list

ids = read_xlsx('KEIO_dev_data/Keio Library.xlsx', sheet = 'ID') 

data.sum2 = data.sum %>% mutate(Genotype = as.character(Genotype))

data.sum2[data.sum2$Genotype == 'fumD',]$Genotype = 'ydhZ'
data.sum2[data.sum2$Genotype == 'fumE',]$Genotype = 'yggD'
data.sum2[data.sum2$Genotype == 'gpmM',]$Genotype = 'gpmI'
data.sum2[data.sum2$Genotype == 'nudl',]$Genotype = 'yfaO'
data.sum2[data.sum2$Genotype == 'ppnP',]$Genotype = 'yaiE'
data.sum2[data.sum2$Genotype == 'ppsa',]$Genotype = 'pps'
data.sum2[data.sum2$Genotype == 'pyrl',]$Genotype = 'pyrL'

genes = unique(as.character(data.sum2$Genotype))
id2 = ids[ids$ecocyc %in% genes,]
colnames(id2)[2] = 'Genotype'

# there are some genes without entry
setdiff(genes, id2$Genotype)


genes2 = data.sum2 %>%
	filter(Genotype != 'BW') %>%
	left_join(id2, by = 'Genotype') %>%
	dplyr::filter(Drug == 2.5, Supplement_mM == 0) %>%
	dplyr::select(KEGG, BW_norm) %>%
	data.frame

rownames(genes2) = genes2[,1]; genes2[,1] = NULL


# to see all the pathways
data(paths.hsa)
paths.hsa
# remember to change the three letter code for the species
species = 'ecj'

pv.out = pathview(gene.data = genes2, pathway.id = 'ecj00020', species = "ecj", gene.idtype = "kegg", out.suffix = "output", 
	keys.align = "y", kegg.native = T)





























