# This script analyses different data sets of Devlopmental Assays with different compounds, 
# and with different conditions. This version uses some of my (simple) R functions, as the 
# color selection for plots with several variables
#  
# For more information, visit my github repository: https://github.com/dmartimarti/R_functions

# setwd("/Users/dmarti14/Documents/MRC_Postdoc/Projects/Leo")

# libraries
library(tidyverse)
library(readxl)
library(ggrepel)
library(caret)
library(Rtsne)

# my library of functions
source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')

# session options
options(width = 220)

# every plot will go to Summary folder
odir <- 'Summary'
keiodir = '/KEIO'
dir.create(paste0(odir,keiodir), showWarnings = TRUE, recursive = FALSE, mode = "0777")



# load data
## REMEMBER, PLATE A AND PLATE B ARE NOW 97 AND 99
data_dev = read_xlsx('KEIO_dev_data/Keio Library.xlsx', sheet = 'ALL_final') 

data_dev = data_dev %>%
	filter(!is.na(Gene), !(Gene %in% c('WT', 'present', 'no bact' )))
## running some checks

# is there any gene repeated? 
# yes we have
subset = data_dev %>% filter(Gene != 'Control')
length(unique(subset$Gene)) == length(subset$Gene)

# dup = data_dev$Gene[duplicated(data_dev$Gene)]

# for (name in dup){
# 	print(paste0('changing gene ', name))
# 	for (i in (1:dim(filter(data_dev, Gene == name))[1])) {
# 		print(as.character(dim(filter(data_dev, Gene == name))[1]))
# 		print(paste0('looping over ', name, ' ,iteration n. ', i))
# 		data_dev[data_dev$Gene == name, ]$Gene[i] = paste0(name, as.character(i))
# 	}
# }


# transform data into 1 column data frame
data_dev = data_dev %>% 
	gather(Sample, Score, Control_1:GluFU_3) %>%
	separate(Sample, c('Sample', 'Replicate'), sep = '_') %>%
	mutate_at(c('Plate', 'Column', 'Row', 'Gene', 'Sample'), as.factor) %>%
	mutate(Score = as.numeric(Score))


# Controls median
controls = data_dev %>%
	filter(Gene == 'Control') %>%
	mutate(Score = as.numeric(Score)) %>%
	group_by(Sample) %>%
	summarise(Median = median(Score))

# we can use the values of the first row of the Controls as the median, it has the same values
controls = filter(data_dev, Gene == 'Control', Sample == 'Control')[1,]
controls = rbind(controls, filter(data_dev, Gene == 'Control', Sample == 'FU')[1,])
controls = rbind(controls, filter(data_dev, Gene == 'Control', Sample == 'Glu')[1,])
controls = rbind(controls, filter(data_dev, Gene == 'Control', Sample == 'GluFU')[1,])


## summarise the data
data.sum = data_dev %>%
	filter(Gene != 'Control') %>%
	rbind(., controls) %>%
	group_by(Sample, Gene) %>%
	summarise(Median_Score = median(Score, na.rm = TRUE),
			  Mean = mean(Score, na.rm = TRUE),
			  SD = sd(Score, na.rm = TRUE)) %>%
	mutate(BW_Score = Median_Score[Gene == 'Control'],
	       BW_Mean = Mean[Gene == 'Control']) %>%
	ungroup %>%
	mutate(BW_norm = Median_Score - BW_Score,
	       BW_Mean_norm = Mean - BW_Mean,
	       Supplement = 'Glucose')

# adding a helper column
data.sum['Supplement_mM'] = 0
data.sum[data.sum$Sample == 'Glu',]$Supplement_mM = 10
data.sum[data.sum$Sample == 'GluFU',]$Supplement_mM = 10

###
### scatter plot with everything
###

pos = position_jitter(width = 0.15, height = 0.15, seed = 1) # to plot names in jitter positions
data.sum %>% 
	filter(Sample %in% c('FU', 'GluFU')) %>%
	select(Supplement, Supplement_mM, Gene, BW_norm) %>%
	unite(Supp, Supplement, Supplement_mM) %>%
	spread(Supp, BW_norm) %>%
	ggplot(aes(x = Glucose_0, y = Glucose_10)) +
	geom_hline(yintercept = 0, colour = 'grey30') +
	geom_vline(xintercept = 0, colour = 'grey30') +
	geom_point(position = pos, size = 2, alpha = .2) +
	# geom_text_repel(aes(label = ifelse(Glucose_0 > 1.5 & abs(Glucose_10) > 0.5, as.character(Gene), '')), position = pos) +
	# geom_text_repel(aes(label = ifelse((abs(Glucose_0) > 1.5 | abs(Glucose_10) > 0.5), as.character(Gene), '')), position = pos) +
	labs(title = "5FU + Supplement (Glucose) effect on different BW mutants",
		 x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype', sep = ' ')),
		 y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype ' , bold('(Glucose)'), sep = ' '))) +
	theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
		panel.grid.major = element_line(colour = "grey90"),
		panel.background = element_rect(fill = "white", colour = "grey50")) 

dev.copy2pdf(device = cairo_pdf,
             file = paste(odir, keiodir,"/ALL_scatterplot.pdf", sep=''),
             width = 14, height = 10, useDingbats = FALSE)


# how many points do we have in the big group?
big_group = data.sum %>% 
	filter(Sample %in% c('FU', 'GluFU')) %>%
	select(Supplement, Supplement_mM, Gene, BW_norm) %>%
	unite(Supp, Supplement, Supplement_mM) %>%
	spread(Supp, BW_norm) %>%
	filter(Glucose_0 == 0, Glucose_10 == -1)



########
# plot with sublibrary names

data_dev_sub = read_xlsx('Develop_data/Summary_Keio Sublibrary_Glucose Supp_10mM_N2 worms_08-12-18.xlsx', sheet = 'Summary') 
data_dev_sub$X__1 = NULL

sub_names = unique(data_dev_sub$Genotype)
sub_names = data.frame(sub_names)
colnames(sub_names) = 'Gene'
sub_names['Genotype'] = unique(data_dev_sub$Genotype)

sublibrary = data.sum %>% 
	filter(Sample %in% c('FU', 'GluFU')) %>%
	select(Supplement, Supplement_mM, Gene, BW_norm) %>%
	unite(Supp, Supplement, Supplement_mM) %>%
	spread(Supp, BW_norm) %>%
	left_join(sub_names)



pos = position_jitter(width = 0.15, height = 0.15, seed = 1) # to plot names in jitter positions
sublibrary %>% 
	ggplot(aes(x = Glucose_0, y = Glucose_10)) +
	geom_hline(yintercept = 0, colour = 'grey30') +
	geom_vline(xintercept = 0, colour = 'grey30') +
	geom_point(position = pos, size = 2, alpha = .2) +
	# geom_text_repel(aes(label = ifelse(Glucose_0 > 1.5 & abs(Glucose_10) > 0.5, as.character(Gene), '')), position = pos) +
	# geom_text_repel(aes(label = ifelse((abs(Glucose_0) > 1.5 | abs(Glucose_10) > 0.5), as.character(Gene), '')), position = pos) +
	labs(title = "5FU + Supplement (Glucose) effect on different BW mutants",
		 x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype', sep = ' ')),
		 y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype ' , bold('(Glucose)'), sep = ' '))) +
	theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
		panel.grid.major = element_line(colour = "grey90"),
		panel.background = element_rect(fill = "white", colour = "grey50")) +
	geom_text_repel(aes(label = Genotype))







##################### 
### batch effect? ###
#####################


# load data
## REMEMBER, PLATE A AND PLATE B ARE NOW 97 AND 99
data_dev2 = read_xlsx('Keio Library.xlsx', sheet = 'ALL_final_batch') 

data_dev2 = data_dev2 %>%
	filter(!is.na(Gene), !(Gene %in% c('WT', 'present', 'no bact' ))) # remove empty wells, WT, and other things

## running some checks

# is there any gene repeated? 
# yes we have
subset = data_dev2 %>% filter(Gene != 'Control')
data_dev2 = data_dev2 %>% filter(Gene != 'Control')
length(unique(subset$Gene)) == length(subset$Gene)

# dup = data_dev$Gene[duplicated(data_dev$Gene)]

# for (name in dup){
# 	print(paste0('changing gene ', name))
# 	for (i in (1:dim(filter(data_dev, Gene == name))[1])) {
# 		print(as.character(dim(filter(data_dev, Gene == name))[1]))
# 		print(paste0('looping over ', name, ' ,iteration n. ', i))
# 		data_dev[data_dev$Gene == name, ]$Gene[i] = paste0(name, as.character(i))
# 	}
# }


data_dev2[data_dev2$Batch == 1,]$Batch = 'Batch_1'
data_dev2[data_dev2$Batch == 2,]$Batch = 'Batch_2'
data_dev2[data_dev2$Batch == 3,]$Batch = 'Batch_3'
data_dev2[data_dev2$Batch == 4,]$Batch = 'Batch_4'
data_dev2[data_dev2$Batch == 5,]$Batch = 'Batch_5'
data_dev2[data_dev2$Batch == 6,]$Batch = 'Batch_6'


batch = data_dev2 %>%
	filter(Gene != 'Control') %>%
	mutate(Batch = as.factor(Batch), 
		Plate = as.factor(Plate)) %>%
	select(Gene, Batch) %>%
	unique(.)


# transform data into 1 column data frame
data_dev2 = data_dev2 %>% 
	select(-Batch) %>%
	gather(Sample, Score, Control_1:GluFU_3) %>%
	separate(Sample, c('Sample', 'Replicate'), sep = '_') %>% 
	mutate_at(c('Plate', 'Column', 'Row', 'Gene', 'Sample'), as.factor) 


# # Controls median
# controls = data_dev2 %>%
# 	filter(Gene == 'Control') %>%
# 	group_by(Sample) %>%
# 	summarise(Median = median(Score))

# # we can use the values of the first row of the Controls as the median, it has the same values
# controls = filter(data_dev2, Gene == 'Control', Sample == 'Control')[1,]
# controls = rbind(controls, filter(data_dev2, Gene == 'Control', Sample == 'FU')[1,])
# controls = rbind(controls, filter(data_dev2, Gene == 'Control', Sample == 'Glu')[1,])
# controls = rbind(controls, filter(data_dev2, Gene == 'Control', Sample == 'GluFU')[1,])


## summarise the data
data.sum2 = data_dev2 %>%
	filter(Gene != 'Control') %>%
	rbind(., controls) %>%
	group_by(Sample, Gene) %>%
	summarise(Median_Score = median(Score, na.rm = TRUE),
			  Mean = mean(Score, na.rm = TRUE),
			  SD = sd(Score, na.rm = TRUE)) %>%
	mutate(BW_Score = Median_Score[Gene == 'Control'],
	       BW_Mean = Mean[Gene == 'Control']) %>%
	ungroup %>%
	mutate(BW_norm = Median_Score - BW_Score,
	       BW_Mean_norm = Mean - BW_Mean,
	       Supplement = 'Glucose') %>%
	filter(Gene != 'Control') %>%
	left_join(batch)

# adding a helper column
data.sum2['Supplement_mM'] = 0
data.sum2[data.sum2$Sample == 'Glu',]$Supplement_mM = 10
data.sum2[data.sum2$Sample == 'GluFU',]$Supplement_mM = 10

###
### scatter plot with everything
###

pos = position_jitter(width = 0.15, height = 0.15, seed = 1) # to plot names in jitter positions
data.sum2 %>% 
	filter(Sample %in% c('FU', 'GluFU')) %>%
	select(Supplement, Supplement_mM, Gene, BW_norm, Batch) %>%
	unite(Supp, Supplement, Supplement_mM) %>%
	spread(Supp, BW_norm) %>%
	ggplot(aes(x = Glucose_0, y = Glucose_10, color = Batch)) +
	geom_hline(yintercept = 0, colour = 'grey30') +
	geom_vline(xintercept = 0, colour = 'grey30') +
	geom_point(position = pos, size = 2, alpha = .5) +
	# geom_text_repel(aes(label = ifelse(Glucose_0 > 1.5 & abs(Glucose_10) > 0.5, as.character(Gene), '')), position = pos) +
	# geom_text_repel(aes(label = ifelse((abs(Glucose_0) > 1.5 | abs(Glucose_10) > 0.5), as.character(Gene), '')), position = pos) +
	labs(title = "5FU + Supplement (Glucose) effect on different BW mutants",
		 x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype', sep = ' ')),
		 y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype ' , bold('(Glucose)'), sep = ' '))) +
	theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
		panel.grid.major = element_line(colour = "grey90"),
		panel.background = element_rect(fill = "white", colour = "grey50")) 

dev.copy2pdf(device = cairo_pdf,
             file = paste(odir,"/ALL_scatterplot_batch.pdf", sep=''),
             width = 14, height = 10, useDingbats = FALSE)


big_group = data.sum2 %>% 
	filter(BW_norm == -1, Sample == 'GluFU') 

big_group %>%
	group_by(Batch) %>%
	summarise(Count = n()) %>%
	ggplot(aes(y = Count, x = Batch, fill = Batch)) + geom_bar(stat = 'identity')











############################
### t-SNE transformation ###
############################


# transform the data, and then, plot
# change the BW_norm if you want to work with median instead of means
tsne_data = data.sum %>%
	filter(!is.na(BW_norm)) %>% 
	select(Sample, Gene, BW_norm) %>%
	spread(Sample, BW_norm) 
tsne_mat = as.matrix(tsne_data[,2:length(tsne_data)])
rownames(tsne_mat) = tsne_data$Gene

# t-SNE model
tsne_model_1 = Rtsne(as.matrix(tsne_mat), check_duplicates = FALSE, pca = TRUE, perplexity = 50, theta = 0.5, dims = 2, max_iter = 5000, num_threads = 6)

## getting the two dimension matrix
d_tsne_1 = as.data.frame(tsne_model_1$Y)
d_tsne_1['Gene'] = tsne_data$Gene

## plotting the results without clustering
d_tsne_1 %>%
    ggplot(aes(x = V1, y = V2)) +
    geom_point(size = 1, alpha = 0.4) +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    xlab("") + ylab("") +
    ggtitle("t-SNE") +
    theme_light(base_size = 20) # +
    # geom_text_repel((aes(label = Gene)))

dev.copy2pdf(device = cairo_pdf,
             file = paste(odir,"/tSNE_worm_glucose.pdf", sep=''),
             width = 14, height = 10, useDingbats = FALSE)


# filter data from the islands (this will change every time you run the tsne model)

g1 = d_tsne_1 %>% filter(V2 > 0, V1 < -20) 
g2 = d_tsne_1 %>% filter(V2 > 20) 
g3 = d_tsne_1 %>% filter(V1 > 20) 
g4 = d_tsne_1 %>% filter(V2 < -20)




#############################################
### Pathways and gene enrichment analysis ###
#############################################



# read ID list

ids = read_xlsx('Keio Library.xlsx', sheet = 'ID') 

genes = unique(as.character(data.sum$Gene))
id2 = ids[ids$ecocyc %in% genes,]
colnames(id2)[2] = 'Gene'

# there are some genes without entry
setdiff(genes, id2$Gene)


genes = data.sum %>%
	filter(Gene != 'Control') %>%
	left_join(id2, by = 'Gene') %>%
	select(Sample, Gene, KEGG, everything())



test = genes %>%
	filter(Sample == 'FU') %>%
	select(KEGG, BW_norm) %>%
	data.frame

rownames(test) = test[,1]; test[,1] = NULL


# to see all the pathways
pathview::data(paths.ecj)
# remember to change the three letter code for the species
species = 'ecj'

# plot ALL KEGG map, not very useful
data("bods", package = "pathview")
pv.out = pathview::pathview(gene.data = test, pathway.id = 'ecj01100', species = "ecj", gene.idtype = "kegg", out.suffix = "output", 
	keys.align = "y", kegg.native = T)

# plot KEGG map
data("bods", package = "pathview")
pv.out = pathview::pathview(gene.data = test, pathway.id = 'ecj01200', species = "ecj", gene.idtype = "kegg", out.suffix = "output", 
	keys.align = "y", kegg.native = T)






























