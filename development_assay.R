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

# write a backup of the data to use
data_dev %>%
	write_csv(paste0(odir,'/KeioSublibrary_08-12-18_data_raw.csv'))


# table with genotypes and pathways
pathways = data_dev %>% select(Genotype, Pathway) %>% distinct(Genotype, .keep_all = TRUE)
pathways %>%
	write_csv(paste0(odir,'/EcoCyc_pathways_KeioSublibrary_08-12-18.csv'))

# after modifying the pathways file, we have reduced the categories
pathways = read_xlsx('Summary/EcoCyc_pathways_KeioSublibrary_08-12-18.xlsx', sheet = 'paths')

pathways = pathways %>% 
	mutate(Genotype = as.factor(Genotype),
		   Pathway = as.factor(Pathway),
		   KEGG1 = as.factor(KEGG1),
		   KEGG2 = as.factor(KEGG2),
		   KEGG3 = as.factor(KEGG3))

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


# test data
data.sum %>% filter(Genotype == 'aceE')


# plot heatmaps
gradcolours<-c('black', 'yellow', 'orange', 'red')

diffcolours<-c('purple', 'blue', 'steelblue', 'black', 'yellow', 'orange', 'red')
breaks<-c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
limits<-c(-4 , 4)

# function to plot heatmap, from Pov
plotHeatmap <- function(data, x, y, fill, grad = gradcolours, breaks = c(0,1,2,3,4), limits = c(0,4)){
  ggplot(data,aes_string(x = x, y = y, fill = fill))+
    geom_tile()+
    xlab('5FU, uM')+
    theme(axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(colour = NA),
          axis.line.x = element_line(colour = NA),
          axis.line.y = element_line(colour = NA),
          strip.text = element_text(colour = 'black', face='bold',size = 10),
          axis.text.x= element_text(face = 'bold', colour = 'black', size = 10, angle = 90, hjust = 1),
          axis.text.y= element_text(face = 'bold', colour = 'black', size = 10))+
    scale_fill_gradientn(colours = grad, breaks = breaks, limits = limits, guide = "legend")
}


#### heatmaps


data.sum %>%
  plotHeatmap(x = 'Drug', y = 'Genotype', fill = 'Median_Score')+
  ylab('Gene')+
  labs(fill = 'C. elegans phenotype')+
  facet_wrap(~ Supplement + Supplement_mM, ncol = 2)


dev.copy2pdf(device = cairo_pdf,
             file = paste(odir,"/Worm_phenotype_comparison_Glucose.pdf", sep=''),
             width = 8, height = 20, useDingbats = FALSE)



data.sum %>%
  filter(Genotype != 'BW') %>%
  plotHeatmap(x = 'Drug', y = 'Supplement_mM', fill = 'BW_norm', grad = diffcolours, breaks = breaks, limits = limits)+
  ylab('Supplement, mM')+
  labs(fill='C. elegans phenotype')+
  facet_wrap(~Genotype, ncol = 8)

dev.copy2pdf(device = cairo_pdf,
             file = paste0(odir,"/Worm_phenotype_comparison_N2_vs_BW_Glucose.pdf"),
             width = 10, height = 10, useDingbats = FALSE)


#### Scatterplots 

# Score vs BW


KOscores = data.sum %>%
  filter(Genotype != 'BW') %>%
  mutate(Suppl_mM = as.numeric(as.character(Supplement_mM)),
         Drug = as.numeric(as.character(Drug))) %>%
  group_by(Supplement, Supplement_mM, Genotype) %>%
  summarise(Sum = sum(BW_norm),
            Mean = mean(BW_norm),
            Median = median(BW_norm),
            WMean = mean(Drug * BW_norm)/sum(Drug),
            WSum = sum(Drug * BW_norm)/sum(Drug)) %>%
  gather(Stat, Value, Sum, Mean, Median, WMean, WSum) %>%
  group_by(Supplement, Supplement_mM, Stat) %>%
  mutate(Quantile = ecdf(Value)(Value) * 100,
         Q = Quantile - 50,
         Z = scale(Value)) 

KOscores.sum = KOscores %>%
  filter(Stat %in% c('Sum','WSum','WMean')  ) %>%
  gather(Measure,Value,Value,Quantile,Z,Q) %>%
  unite(MS, Measure, Supplement_mM) %>%
  spread(MS,Value) %>%
  mutate(IQ=sqrt(Q_0^2+Q_10^2),
         IQQ=ecdf(IQ)(IQ)*100)


PlotScatter <- function(data, x = 'Value_0', y = 'Value_10'){
    data %>%
    ggplot(aes_string(x = x, y = y))+
    geom_hline(aes(yintercept = 0), color = 'red', alpha = 0.5)+
    geom_vline(aes(xintercept = 0), color = 'red', alpha = 0.5)+
    geom_jitter(aes(color = IQQ))+
    # ggtitle(paste0('Consistency of KO effects on NGM and +10mM Supplement (BW normalised): ',as.character(unique(data$Stat))),,
    #         subtitle='KOs are marked by their interaction quantile, top: 5% - dark red, 10% - red, 20% - cyan')+
    ggtitle('KOs are marked by their interaction quantile, top: 5% - dark red, 10% - red, 20% - light red')+
    scale_x_continuous(breaks = seq(-100, 100, by = 0.5))+
    scale_y_continuous(breaks = seq(-100, 100, by = 0.5))+
    xlab('C. elegans phenotype vs BW on NGM')+
    ylab('C. elegans phenotype vs BW on NGM + 10mM Supplement')+
    labs(color='Interaction\nquantile')+
    scale_color_gradientn(colours = c('gray80','red4'), breaks = seq(0, 100, by = 25), guide = "legend")+
    geom_text_repel(aes(label = ifelse(IQQ >= 80 & IQQ < 90, as.character(Genotype), '')), color = 'salmon', size = 4)+
    geom_text_repel(aes(label = ifelse(IQQ >= 90 & IQQ < 95, as.character(Genotype), '')), color = 'red2', size = 4)+
    geom_text_repel(aes(label = ifelse(IQQ >= 95, as.character(Genotype),'')), color = 'red4', size = 4)
}


#KOGLucose
KOGlucoseplots <- KOscores.sum %>%
  filter(Supplement == "Glucose") %>%
  group_by(Stat) %>%
  do(plot = PlotScatter(.)) 

map2(paste0(odir,"/Scatter_KO_AllWorms_Glucose_",as.character(KOGlucoseplots$Stat),".pdf"),KOGlucoseplots$plot,
   width = 6, height = 6, useDingbats = FALSE, ggsave)

# code to colour different pathways

'#000000' # control
'#0A2CFF' # dark blue
'#008BF6' # azure for de novo ribonucleotide
'#FF1500' # red for salvage ribonucleotide
'#5EAA44' # salvage deoxy
'#8C158C' # denovo deoxy

# range of colours
colours17 = c('#000000','#8C158C','#8C158C','#0A2CFF','#006006',
			'#FF008C','#9D4300','#7B008C','#008BF6','#5EAA44',
			'#FF1500','#A42900','#006A54','#00CDFF','#F0BF8F',
			'#B5B0B8','#FAE603')
colours_kegg2 = c('#F0BF8F','#000000','#8C158C','#0A2CFF','#006006',
				  '#008BF6','#5EAA44','#FF1500','#CC80E6','#BAE303')
colours_kegg3 = c('#8A7090','#000000','#53BAE2','#EF6713','#58A32D','#EF1613','#F45A95','#8C5A3A')


pos = position_jitter(width = 0.05, height = 0.05, seed = 1) # to plot names in jitter positions

# function to plot. Remember to change the color coding!
PlotScatter2 = function(data, x = 'Glucose_0', y = 'Glucose_10', grad = colours_kegg3){
    data %>%
    ggplot(aes_string(x = x, y = y)) + 
	geom_hline(yintercept = 0, colour = 'grey30') +
	geom_vline(xintercept = 0, colour = 'grey30') +
	geom_point(aes(colour = KEGG3), position = pos, size = 2) + 
	scale_color_manual(values = grad) + 
    scale_fill_manual(values = grad) +
	# geom_text_repel(aes(label = ifelse(Pathway == 'ribonucleotides de novo biosynthesis', as.character(Genotype), '')), position = pos)
	geom_text_repel(aes(label = Genotype, colour = KEGG3), position = pos) +
	labs(title = "5FU + Supplement (Glucose) effect on different BW mutants",
		 x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype', sep = ' ')),
		 y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype ' , bold('(Glucose)'), sep = ' '))) +
	theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
		panel.grid.major = element_line(colour = "grey90"),
		panel.background = element_rect(fill = "white", colour = "grey50")) + 
	guides(colour = guide_legend(override.aes = list(size = 4))) # make lengend points larger
}

# transform the data, and then, plot
data.sum %>% 
	left_join(pathways) %>%
	filter(Drug == 5) %>% 
	select(Supplement, Supplement_mM, Genotype, BW_norm, Pathway:KEGG3) %>%
	unite(Supp, Supplement, Supplement_mM) %>%
	spread(Supp, BW_norm) %>%
	PlotScatter2(.)


dev.copy2pdf(device = cairo_pdf,
             file = paste(odir,"/Scatter_worm_glucose_5uM_5FU_8cat.pdf", sep=''),
             width = 14, height = 10, useDingbats = FALSE)


### spectral plot


# transform the data, and then, plot
data.sum %>% 
    left_join(pathways) %>%
    select(Supplement, Supplement_mM, Genotype, Drug, BW_norm, Pathway:KEGG3) %>%
    unite(Supp, Supplement, Supplement_mM, Drug) %>%
    spread(Supp, BW_norm) %>% 
    ggplot(aes(x = Glucose_0_0, y = Glucose_10_0)) + 
    geom_hline(yintercept = 0, colour = 'grey30') +
    geom_vline(xintercept = 0, colour = 'grey30') +
    # geom_point(position = pos, size = 2, colour = 'black') +
    # geom_point(aes(x = Glucose_0_1, y = Glucose_10_1), position = pos, size = 2, colour = 'red') +
    # geom_point(aes(x = Glucose_0_2.5, y = Glucose_10_2.5), position = pos, size = 2, colour = 'blue') +
    geom_point(aes(x = Glucose_0_5, y = Glucose_10_5), position = pos2, size = 2, colour = '#8D3B72') +
    geom_segment(aes(x = Glucose_0_0, xend = Glucose_0_1, y = Glucose_10_0, yend = Glucose_10_1), position = pos, size = 1, colour = 'black', alpha = 0.2, arrow = arrow(length = unit(0.01, "npc"))) +
    geom_segment(aes(x = Glucose_0_1, xend = Glucose_0_2.5, y = Glucose_10_1, yend = Glucose_10_2.5), position = pos2, size = 1, colour = 'red', alpha = 0.2, arrow = arrow(length = unit(0.01, "npc"))) +
    geom_segment(aes(x = Glucose_0_2.5, xend = Glucose_0_5, y = Glucose_10_2.5, yend = Glucose_10_5), position = pos2, size = 1, colour = 'blue', alpha = 0.2, arrow = arrow(length = unit(0.01, "npc"))) +
    labs(title = expression(paste("5FU + Glucose + gltA mutants effect on ", italic('C. elegans'), " phenotype", sep = '')),
         x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype', sep = ' ')),
         y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype ' , bold('(Glucose)'), sep = ' '))) +
    geom_text_repel(aes(x = Glucose_0_5, y = Glucose_10_5, label = Genotype), position = pos2) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_rect(fill = "white", colour = "grey50")) + 
    guides(colour = guide_legend(override.aes = list(size = 4))) # make lengend points larger




### t-SNE plot (test) ###



# transform the data, and then, plot
# change the BW_norm if you want to work with median instead of means
tsne_data = data.sum %>% 
    left_join(pathways) %>%
    select(Supplement, Supplement_mM, Genotype, Drug, BW_norm) %>%
    unite(Supp, Supplement, Supplement_mM, Drug) %>%
    spread(Supp, BW_norm)

tsne_mat = as.matrix(tsne_data[,2:9])
rownames(tsne_mat) = tsne_data$Genotype

# t-SNE model
tsne_model_1 = Rtsne(as.matrix(tsne_mat), check_duplicates = FALSE, pca = TRUE, perplexity = 15, theta = 0.1, dims = 2, max_iter = 10000)

## getting the two dimension matrix
d_tsne_1 = as.data.frame(tsne_model_1$Y)
d_tsne_1['Genotype'] = tsne_data$Genotype

## plotting the results without clustering
d_tsne_1 %>%
    left_join(pathways) %>%
    ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(colour = KEGG3), size = 2, alpha = 0.8) +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    xlab("") + ylab("") +
    labs(title = expression(paste("t-SNE of KEIO sub-library from ", italic('C. elegans'), " phenotype", sep = ''))) +
    theme_light(base_size = 20) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
    geom_text_repel((aes(label = Genotype, colour = KEGG3))) +
    scale_color_manual(values = colours_kegg3) + 
    scale_fill_manual(values = colours_kegg3)

dev.copy2pdf(device = cairo_pdf,
             file = paste(odir,"/tSNE_worm_glucose_8cat.pdf", sep=''),
             width = 14, height = 10, useDingbats = FALSE)


### adding some clustering to our data, but not really sure if this works
## keeping original data
d_tsne_1_original=d_tsne_1

## Creating k-means clustering model, and assigning the result to the data used to create the tsne
fit_cluster_kmeans=kmeans(scale(d_tsne_1[,1:2]), 5)
d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)

## Creating hierarchical cluster model, and assigning the result to the data used to create the tsne
fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1[,1:2])))

## setting 3 clusters as output
d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=5))

## plotting the results without clustering
d_tsne_1_original %>%
    left_join(pathways) %>%
    ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(colour = cl_hierarchical), size = 1) +
    # scale_color_manual(values = colsel(8)) + 
    # scale_fill_manual(values = colsel(8)) +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    xlab("") + ylab("") +
    ggtitle("t-SNE") +
    theme_light(base_size = 20) +
    geom_text_repel((aes(label = Genotype, colour = KEGG3)))








####################################
#### Transcription factor stuff ####
####################################

# load data

tf_data = read_xlsx('Develop_data/Summary_TF Sublibrary_Glucose Supp_10mM_N2 worms.xlsx', sheet = 'Summary') 

tf_data = tf_data %>%
	gather(Drug, Score, `0`, `1`, `2.5`, `5`) %>% # make table long
	unite(ID, Gene, Supplement, Replicate, remove = FALSE) %>% 
	mutate(Replicate = as.numeric(Replicate),
		Supplement_mM = ifelse(Supplement == 'Control', 0, 10),
		Supplement_mM = as.factor(Supplement_mM),
		Supplement = as.factor(Supplement),
		Gene = as.factor(Gene),
		ID = as.factor(ID),
		Drug = as.factor(Drug)) %>%
	select(Gene, Well, ID, Replicate, Supplement, Supplement_mM, Drug, Score)

# write a backup of the data to use
tf_data %>%
	write_csv(paste0(odir,'/TFSublibrary_12-12-18_data_raw.csv'))



# summary of the data

tf.sum <- tf_data %>%
  group_by(Supplement, Supplement_mM, Drug, Gene) %>%
  summarise(Median_Score = median(Score, na.rm = TRUE),
  			Mean = mean(Score, na.rm = TRUE),
  			SD = sd(Score, na.rm = TRUE)) %>%
  mutate(BW_Score = Median_Score[Gene == 'BW']) %>%
  ungroup %>%
  group_by(Gene, Drug) %>%
  mutate(No_supplement = ifelse("0" %in% Supplement_mM, Median_Score[Supplement_mM == '0'], NA)) %>% # Get phenotype values at Supplement_mM=0 
  ungroup %>%
  mutate(BW_norm = Median_Score - BW_Score,
         Suppl_norm = Median_Score - No_supplement,
         Interaction = Median_Score - BW_Score - No_supplement,
         Supplement = 'Glucose')


# test data
tf.sum %>% filter(Gene == 'AgaR')

# get scores
TFscores = tf.sum %>%
  filter(Gene != 'BW') %>%
  mutate(Suppl_mM = as.numeric(as.character(Supplement_mM)),
         Drug = as.numeric(as.character(Drug))) %>%
  group_by(Supplement, Supplement_mM, Gene) %>%
  summarise(Sum = sum(BW_norm),
            Mean = mean(BW_norm),
            Median = median(BW_norm),
            WMean = mean(Drug * BW_norm)/sum(Drug),
            WSum = sum(Drug * BW_norm)/sum(Drug)) %>%
  gather(Stat, Value, Sum, Mean, Median, WMean, WSum) %>%
  group_by(Supplement, Supplement_mM, Stat) %>%
  mutate(Quantile = ecdf(Value)(Value) * 100,
         Q = Quantile - 50,
         Z = scale(Value)) 

# summary with weird quantile stuff
TFscores.sum = TFscores %>%
  filter(Stat %in% c('Sum','WSum','WMean')  ) %>%
  gather(Measure,Value,Value,Quantile,Z,Q) %>%
  unite(MS, Measure, Supplement_mM) %>%
  spread(MS,Value) %>%
  mutate(IQ = sqrt(Q_0^2+Q_10^2),
         IQQ = ecdf(IQ)(IQ)*100)


# transform the data, and then, plot
pos = position_jitter(width = 0.05, height = 0.05, seed = 1) # to plot names in jitter positions
tf.sum %>% 
	filter(Drug == 5) %>% 
	select(Supplement, Supplement_mM, Gene, BW_norm) %>%
	unite(Supp, Supplement, Supplement_mM) %>%
	spread(Supp, BW_norm) %>%
	ggplot(aes(x = Glucose_0, y = Glucose_10)) + 
	geom_hline(yintercept = 0, colour = 'grey30') +
	geom_vline(xintercept = 0, colour = 'grey30') +
	geom_point(position = pos, size = 2) + 
	# scale_color_manual(values = grad) + 
    # scale_fill_manual(values = grad) +
	geom_text_repel(aes(label = ifelse((abs(Glucose_0) > 0.5 | abs(Glucose_10) > 0.5) & Gene != 'CRP', as.character(Gene), '')), position = pos) +
	geom_text_repel(aes(label = ifelse(Gene == 'CRP', 'CRP', '')), size = 6, color = 'red', position = pos) +
	# geom_text_repel(aes(label = Gene), position = pos) +
	labs(title = expression(paste("5FU + Glucose + TF mutants effect on ", italic('C. elegans'), " phenotype", sep = '')),
		 x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype', sep = ' ')),
		 y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype ' , bold('(Glucose)'), sep = ' '))) +
	theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
		panel.grid.major = element_line(colour = "grey90"),
		panel.background = element_rect(fill = "white", colour = "grey50")) + 
	guides(colour = guide_legend(override.aes = list(size = 4))) # make lengend points larger

dev.copy2pdf(device = cairo_pdf,
             file = paste(odir,"/Scatter_worm_TF_glucose_5uM_5FU.pdf", sep=''),
             width = 14, height = 10, useDingbats = FALSE)



### t-SNE transformation


# transform the data, and then, plot
# change the BW_norm if you want to work with median instead of means
tsne_data = tf.sum %>% 
    select(Supplement, Supplement_mM, Gene, Drug, BW_norm) %>%
    unite(Supp, Supplement, Supplement_mM, Drug) %>%
    spread(Supp, BW_norm)

tsne_mat = as.matrix(tsne_data[,2:9])
rownames(tsne_mat) = tsne_data$Gene

# t-SNE model
tsne_model_1 = Rtsne(as.matrix(tsne_mat), check_duplicates = FALSE, pca = TRUE, perplexity = 25, theta = 0.3, dims = 2, max_iter = 5000)

## getting the two dimension matrix
d_tsne_1 = as.data.frame(tsne_model_1$Y)
d_tsne_1['Gene'] = tsne_data$Gene

## plotting the results without clustering
d_tsne_1 %>%
    ggplot(aes(x = V1, y = V2)) +
    geom_point(, size = 1) +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    xlab("") + ylab("") +
    ggtitle("t-SNE") +
    theme_light(base_size = 20) +
    geom_text_repel((aes(label = Gene)))

dev.copy2pdf(device = cairo_pdf,
             file = paste(odir,"/tSNE_worm_glucose_8cat.pdf", sep=''),
             width = 14, height = 10, useDingbats = FALSE)



###
### binding both datasets, development and TF 

tf_plot = tf.sum %>% 
	select(Supplement, Supplement_mM, Drug, Gene, BW_norm) %>%
	unite(Supp, Supplement, Supplement_mM, Drug) %>%
	spread(Supp, BW_norm) %>%
    mutate(Category = 'TFactor') %>%
    rename(Genotype = Gene)
dev_plot = data.sum %>% 
	left_join(pathways) %>%
	select(Supplement, Supplement_mM, Genotype, Drug, BW_norm, KEGG3) %>%
	unite(Supp, Supplement, Supplement_mM, Drug) %>%
	spread(Supp, BW_norm) %>%
    select(Genotype, Glucose_0_0:Glucose_10_5, KEGG3) %>%
    rename(Category = KEGG3)

dev.tf_data = rbind(tf_plot, dev_plot)

dev.tf_data = dev.tf_data %>%
    mutate(Group = ifelse((Glucose_0_5 > 1.01 | abs(Glucose_10_5) > 0.5), as.character(Genotype), ''))


colours = c('#8A7090','#000000','#53BAE2','#EF6713','#58A32D','#EF1613','#F45A95','#52489C','#8C5A3A') 
# a simple scatterplot with 5 uM of 5FU
pos = position_jitter(width = 0.05, height = 0.05, seed = 1)
dev.tf_data %>% ggplot(aes(x = Glucose_0_5, y = Glucose_10_5, colour = Category)) + 
    geom_hline(yintercept = 0, colour = 'grey30') +
    geom_vline(xintercept = 0, colour = 'grey30') +
	geom_point(position = pos, size = 2, alpha = 0.8) +
    # text repel with conditionals to filter most of the rubish
	geom_text_repel(aes(label = ifelse(Genotype == 'CRP', 'CRP', '')), size = 3, color = 'red', position = pos) +
    geom_text_repel(aes(label = ifelse((abs(Glucose_0_5) > 0.5 | abs(Glucose_10_5) > 0.5) & Genotype != 'CRP' & 
        (Glucose_0_5 > 1.01 | abs(Glucose_10_5) > 0.5), as.character(Genotype), '')), position = pos) +
    scale_color_manual(values = colours) + 
    scale_fill_manual(values = colours) +
    labs(title = expression(paste("5FU + Glucose + TF/Sublibrary mutants effect on ", italic('C. elegans'), " phenotype", sep = '')),
         x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype', sep = ' ')),
         y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype ' , bold('(Glucose)'), sep = ' ')))

dev.copy2pdf(device = cairo_pdf,
             file = paste(odir,"/Scatter_worm_Merge_glucose.pdf", sep=''),
             width = 14, height = 10, useDingbats = FALSE)



#### t-SNE of the merged dataset

# transform the data, and then, plot
# change the BW_norm if you want to work with median instead of means
tsne_data = dev.tf_data

tsne_mat = as.matrix(tsne_data[,2:9])
rownames(tsne_mat) = tsne_data$Genotype

# t-SNE model
tsne_model_1 = Rtsne(as.matrix(tsne_mat), check_duplicates = FALSE, pca = TRUE, perplexity = 25, theta = 0.3, dims = 2, max_iter = 5000)

## getting the two dimension matrix
d_tsne_1 = as.data.frame(tsne_model_1$Y)
d_tsne_1['Genotype'] = tsne_data$Genotype
d_tsne_1['Category'] = tsne_data$Category
d_tsne_1['Group'] = dev.tf_data$Group

colours_backup = colours
colours = c('#8A7090','#000000','#53BAE2','#EF6713','#58A32D','#EF1613','#F45A95','#52489C','#8C5A3A') 
#colours = colsel(length(unique(d_tsne_1$Category)), palette = 'total', mode = 'random')

## plotting the results without clustering
d_tsne_1 %>%
    ggplot(aes(x = V1, y = V2, colour = Category)) +
    geom_point(size = 2, alpha = 0.8) +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    xlab("") + ylab("") +
    geom_text_repel((aes(label = Group))) +
    scale_color_manual(values = colours) +
    labs(title = expression(paste("t-SNE of KEIO sub-library and TF from ", italic('C. elegans'), " phenotype", sep = ''))) +
    theme_light(base_size = 20) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))

dev.copy2pdf(device = cairo_pdf,
             file = paste(odir,"/tSNE_worm_glucose_Merge.pdf", sep=''),
             width = 14, height = 10, useDingbats = FALSE)



###########################################
#### Methylcitrate stuff (Pyr version) ####
###########################################

# every plot will go to Summary folder

pyr_odir <- '/Methyl_Pyr'
dir.create(paste0(odir, pyr_odir), showWarnings = TRUE, recursive = FALSE, mode = "0777")

# load data

pyr_data = read_xlsx('Develop_data/10. DevelopAssay_glta_2-methylcitrate cycle_sublibrary_5reps_10mM_Pyruvate_13-11-18.xlsx', sheet = 'Summary') 

pyr_data = pyr_data %>%
	gather(Drug, Score, `0`, `1`, `2.5`, `5`) %>% # make table long
	unite(ID, Genotype, Supplement, Replicate, remove = FALSE) %>% 
	mutate(Replicate = as.numeric(Replicate),
		Supplement_mM = ifelse(Supplement == 'Control', 0, 10),
		Supplement_mM = as.factor(Supplement_mM),
		Supplement = as.factor(Supplement),
		Genotype = as.factor(Genotype),
		ID = as.factor(ID),
		Drug = as.factor(Drug)) %>%
	select(Genotype, Well, ID, Replicate, Supplement, Supplement_mM, Drug, Score)



# summary of the data

pyr.sum = pyr_data %>%
  group_by(Supplement, Supplement_mM, Drug, Genotype) %>%
  summarise(Median_Score = median(Score, na.rm = TRUE),
            Mean = mean(Score, na.rm = TRUE),
            SD = sd(Score, na.rm = TRUE)) %>%
  mutate(BW_Score = Median_Score[Genotype == 'BW']) %>%
  ungroup %>%
  group_by(Genotype, Drug) %>%
  mutate(No_supplement = ifelse("0" %in% Supplement_mM, Median_Score[Supplement_mM == '0'], NA)) %>% # Get phenotype values at Supplement_mM=0 
  ungroup %>%
  mutate(BW_norm = Median_Score - BW_Score,
         Suppl_norm = Median_Score - No_supplement,
         Interaction = Median_Score - BW_Score - No_supplement,
         Supplement = 'Pyruvate')


# test data
pyr.sum %>% filter(Genotype == 'ΔgltA:ΔacnB')



# transform the data, and then, plot
pos = position_jitter(width = 0.05, height = 0.05, seed = 1) # to plot names in jitter positions
pyr.sum %>% 
    filter(Drug == 5) %>% 
    select(Supplement, Supplement_mM, Genotype, BW_norm) %>%
    unite(Supp, Supplement, Supplement_mM) %>%
    spread(Supp, BW_norm) %>%
    # filter(!Genotype == 'Aglta AydbK::K') %>%
    ggplot(aes(x = Pyruvate_0, y = Pyruvate_10)) + 
    geom_hline(yintercept = 0, colour = 'grey30') +
    geom_vline(xintercept = 0, colour = 'grey30') +
    geom_point(position = pos, size = 2) + 
    # scale_color_manual(values = grad) + 
    # scale_fill_manual(values = grad) +
    labs(title = expression(paste("5FU + Pyruvate + gltA mutants effect on ", italic('C. elegans'), " phenotype", sep = '')),
         x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype', sep = ' ')),
         y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype ' , bold('(Pyruvate)'), sep = ' '))) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_rect(fill = "white", colour = "grey50")) + 
    guides(colour = guide_legend(override.aes = list(size = 4))) + # make lengend points larger
    geom_text_repel(aes(label = Genotype), position = pos) 

ggsave(file = paste(odir,pyr_odir,"/Scatter_worm_gltA_pyruvate_5uM_5FU.pdf", sep = ''),
       width = 140, height = 120, units = 'mm', scale = 2, device = cairo_pdf, family = "Arial")






####################
### gltA mutants ###
####################



# read data
gltA_data = read_xlsx('Develop_data/6. DevelopAssay_gltA_DMs_sublibrary_4reps_10mM_Glu_06-08-18.xlsx', sheet = 'Summary_deltas') 


# manipulate data
gltA_data = gltA_data %>% 
    gather(Drug, Score, `0`, `1`, `2.5`, `5`) %>% # make table long
    unite(ID, Genotype, Supplement, Replicate, remove = FALSE) %>% 
    mutate(Replicate = as.numeric(Replicate),
        Supplement_mM = ifelse(Supplement == 'Control', 0, 10),
        Supplement_mM = as.factor(Supplement_mM),
        Supplement = as.factor(Supplement),
        Genotype = as.factor(Genotype),
        ID = as.factor(ID),
        Drug = as.factor(Drug)) %>%
    select(Genotype, Well, ID, Replicate, Supplement, Supplement_mM, Drug, Score)





# summary of the data

gltA.sum <- gltA_data %>%
  group_by(Supplement, Supplement_mM, Drug, Genotype) %>%
  summarise(Median_Score = median(Score, na.rm = TRUE),
            Mean = mean(Score, na.rm = TRUE),
            SD = sd(Score, na.rm = TRUE)) %>%
  mutate(BW_Score = Median_Score[Genotype == 'BW']) %>%
  ungroup %>%
  group_by(Genotype, Drug) %>%
  mutate(No_supplement = ifelse("0" %in% Supplement_mM, Median_Score[Supplement_mM == '0'], NA)) %>% # Get phenotype values at Supplement_mM=0 
  ungroup %>%
  mutate(BW_norm = Median_Score - BW_Score,
         Suppl_norm = Median_Score - No_supplement,
         Interaction = Median_Score - BW_Score - No_supplement,
         Supplement = 'Glucose')


# test data
gltA.sum %>% filter(Genotype == 'Aglta AalaA::K')



# transform the data, and then, plot
pos = position_jitter(width = 0.05, height = 0.05, seed = 1) # to plot names in jitter positions
gltA.sum %>% 
    filter(Drug == 5) %>% 
    select(Supplement, Supplement_mM, Genotype, BW_norm) %>%
    unite(Supp, Supplement, Supplement_mM) %>%
    spread(Supp, BW_norm) %>%
    # filter(!Genotype == 'Aglta AydbK::K') %>%
    ggplot(aes(x = Glucose_0, y = Glucose_10)) + 
    geom_hline(yintercept = 0, colour = 'grey30') +
    geom_vline(xintercept = 0, colour = 'grey30') +
    geom_point(position = pos, size = 2) + 
    # scale_color_manual(values = grad) + 
    # scale_fill_manual(values = grad) +
    labs(title = expression(paste("5FU + Glucose + gltA mutants effect on ", italic('C. elegans'), " phenotype", sep = '')),
         x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype', sep = ' ')),
         y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype ' , bold('(Glucose)'), sep = ' '))) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_rect(fill = "white", colour = "grey50")) + 
    guides(colour = guide_legend(override.aes = list(size = 4))) + # make lengend points larger
    geom_text_repel(aes(label = ifelse(Glucose_0 == 2, as.character(Genotype), '')), position = pos, colour = 'blue') +
    geom_text_repel(aes(label = ifelse(Glucose_0 < 2, as.character(Genotype), '')), position = pos, colour = 'red') +
    geom_text_repel(aes(label = ifelse(Genotype == 'Δglta:ΔprpB', as.character(Genotype), '')), position = pos, colour = 'blue')

ggsave(file = paste(odir,"/Scatter_worm_gltA_glucose_5uM_5FU_2_colours.pdf", sep = ''),
       width = 160, height = 160, units = 'mm', scale = 2, device = cairo_pdf, family = "Arial")





# spectral plot

pos = position_jitter(width = 0.05, height = 0.05, seed = 1) # to plot names in jitter positions
pos2 = position_jitter(width = 0.05, height = 0.05, seed = 1)
gltA.sum %>% 
    select(Supplement, Supplement_mM, Drug, Genotype, BW_norm) %>%
    unite(Supp, Supplement, Supplement_mM, Drug) %>%
    spread(Supp, BW_norm) %>%
    filter(!Genotype == 'Aglta AydbK::K') %>%
    ggplot(aes(x = Glucose_0_0, y = Glucose_10_0)) + 
    geom_hline(yintercept = 0, colour = 'grey30') +
    geom_vline(xintercept = 0, colour = 'grey30') +
    # geom_point(position = pos, size = 2, colour = 'black') +
    # geom_point(aes(x = Glucose_0_1, y = Glucose_10_1), position = pos, size = 2, colour = 'red') +
    # geom_point(aes(x = Glucose_0_2.5, y = Glucose_10_2.5), position = pos, size = 2, colour = 'blue') +
    geom_point(aes(x = Glucose_0_5, y = Glucose_10_5), position = pos2, size = 2, colour = '#8D3B72') +
    geom_segment(aes(x = Glucose_0_0, xend = Glucose_0_1, y = Glucose_10_0, yend = Glucose_10_1), position = pos, size = 1, colour = 'black', alpha = 0.2, arrow = arrow(length = unit(0.01, "npc"))) +
    geom_segment(aes(x = Glucose_0_1, xend = Glucose_0_2.5, y = Glucose_10_1, yend = Glucose_10_2.5), position = pos2, size = 1, colour = 'red', alpha = 0.2, arrow = arrow(length = unit(0.01, "npc"))) +
    geom_segment(aes(x = Glucose_0_2.5, xend = Glucose_0_5, y = Glucose_10_2.5, yend = Glucose_10_5), position = pos2, size = 1, colour = 'blue', alpha = 0.2, arrow = arrow(length = unit(0.01, "npc"))) +
    labs(title = expression(paste("5FU + Glucose + gltA mutants effect on ", italic('C. elegans'), " phenotype", sep = '')),
         x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype', sep = ' ')),
         y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype ' , bold('(Glucose)'), sep = ' '))) +
    geom_text_repel(aes(x = Glucose_0_5, y = Glucose_10_5, label = Genotype), position = pos2) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_rect(fill = "white", colour = "grey50")) + 
    guides(colour = guide_legend(override.aes = list(size = 4))) # make lengend points larger






########################
#### Pyruvate stuff ####
########################

# load data

mc_data = read_xlsx('Develop_data/9. DevelopAssay_glta_2-methylcitrate cycle_ASKA_sublibrary_5reps_10mM_Glu_13-11-18.xlsx', sheet = 'Summary') 

mc_data = mc_data %>%
  gather(Drug, Score, `0`, `1`, `2.5`, `5`) %>% # make table long
  unite(ID, Genotype, Supplement, Replicate, remove = FALSE) %>% 
  mutate(Replicate = as.numeric(Replicate),
    Supplement_mM = ifelse(Supplement == 'Control', 0, 10),
    Supplement_mM = as.factor(Supplement_mM),
    Supplement = as.factor(Supplement),
    Genotype = as.factor(Genotype),
    ID = as.factor(ID),
    Drug = as.factor(Drug)) %>%
  select(Genotype, Plasmid, Well, ID, Replicate, Supplement, Supplement_mM, Drug, Score)

# write a backup of the data to use
mc_data %>%
  write_csv(paste0(odir,'/Methylcitrate_ASKAsublibrary_12-12-18_data_raw.csv'))

# boxplots!

# set colours
colours = wesanderson::wes_palette(n = 2, name = "FantasticFox1")
mc_data %>%
  filter(Plasmid %in% c('none', 'prpR') & Drug %in% c(5)) %>%
  ggplot(aes(x = Genotype, y = Score, colour = Plasmid)) +
  geom_boxplot(position = position_dodge(0.9), na.rm = TRUE, outlier.shape = NA) +
  geom_dotplot(aes(fill = Plasmid), colour = 'black', binaxis = 'y', 
    position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1, dodge.width = 0.9), stackdir = 'center', binwidth = 0.03) +
  facet_wrap(vars(Supplement), ncol = 2) +
  scale_color_manual(values = colours) + 
  scale_fill_manual(values = colours) +
  theme(panel.background = element_rect(fill = "white"))


dev.copy2pdf(device = cairo_pdf,
             file = paste(odir,"/Boxplot_gltA_pprP_5uM_5FU.pdf", sep=''),
             width = 8, height = 5, useDingbats = FALSE)

cairo_pdf(file = "ggplot-greek.pdf", width = 8, height = 5)
mc_data %>%
    filter(Plasmid %in% c('none', 'prpR') & Drug %in% c(5)) %>%
    ggplot(aes(x = Genotype, y = Score, colour = Plasmid)) +
    geom_boxplot(position = position_dodge(0.9), na.rm = TRUE, outlier.shape = NA) +
    geom_dotplot(aes(fill = Plasmid), colour = 'black', binaxis = 'y', 
        position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1, dodge.width = 0.9), stackdir = 'center', binwidth = 0.03) +
    facet_wrap(vars(Supplement), ncol = 2) +
    scale_color_manual(values = colours) + 
    scale_fill_manual(values = colours) +
    theme(panel.background = element_rect(fill = "white"))
dev.off()





###########################################
#### Bliss model for drug* interaction ####
###########################################


# let's play with the sublibrary

unique(data.sum$Genotype)


# we will only do things with Drug == 0 & 5

bliss = data.sum %>% filter(Drug %in% c(0,5)) %>% mutate(Median_Score = Median_Score / 4)

# create sub-dataframes with all the info required

yglu = bliss[bliss$Drug == 0,]
yfu = bliss[bliss$Supplement_mM == 0,]
yglufu = bliss[bliss$Drug == 5 & bliss$Supplement_mM == 10,]

# quick check
yglu$Genotype == yfu$Genotype


x = yfu[yfu$Drug == 0,]$Median_Score - yfu[yfu$Drug == 5,]$Median_Score 
y = yglu[yglu$Supplement_mM == 0,]$Median_Score - yglu[yglu$Supplement_mM == 10,]$Median_Score

# calculate the prediction of Bliss model
ypredict = x + y - (x * y)


yglufu['yfu'] = x
yglufu['yglu'] = y
yglufu['ypred'] = ypredict

yglufu = yglufu %>% mutate(eps = ypred - Median_Score)

# function to add names into one more column to see the effect of the interaction
bliss_names = function(x){
    out = c()
    for (i in 1:length(x)){
        if (x[i] >= -0.1 & x[i] <= 0.1) {
            out = c(out, 'independent')
        } else if (x[i] > 0.1){
            out = c(out, 'synergy')
        } else if (x[i] < -0.1) {
            out = c(out, 'antagonistic')
        }
    }
    return(out)
}
# function to plot. Remember to change the color coding!
PlotScatter3 = function(data, x = 'Glucose_0', y = 'Glucose_10', grad = colours_kegg3){
    data %>%
    ggplot(aes_string(x = x, y = y)) + 
    geom_hline(yintercept = 0, colour = 'grey30') +
    geom_vline(xintercept = 0, colour = 'grey30') +
    geom_point(aes(colour = Bliss), position = pos, size = 2) + 
    # scale_color_manual(values = grad) + 
    # scale_fill_manual(values = grad) +
    # geom_text_repel(aes(label = ifelse(Pathway == 'ribonucleotides de novo biosynthesis', as.character(Genotype), '')), position = pos)
    geom_text_repel(aes(label = Genotype, colour = Bliss), position = pos) +
    labs(title = "5FU + Supplement (Glucose), Bliss drug interaction values",
         x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype', sep = ' ')),
         y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype ' , bold('(Glucose)'), sep = ' '))) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_rect(fill = "white", colour = "grey50")) + 
    guides(colour = guide_legend(override.aes = list(size = 4))) # make lengend points larger
}

# final output from a simple Bliss model
yglufu = yglufu %>% group_by(Genotype) %>% mutate(Bliss = bliss_names(eps))


# let's join data.sum and this toghether for plotting purposes

data.sum.bliss = data.sum %>%
    left_join(yglufu %>% select(Genotype, Bliss))

data.sum.bliss %>% 
    left_join(pathways) %>%
    filter(Drug == 5) %>% 
    select(Supplement, Supplement_mM, Genotype, BW_norm, KEGG3, Bliss) %>%
    unite(Supp, Supplement, Supplement_mM) %>%
    spread(Supp, BW_norm) %>%
    PlotScatter3(.)


dev.copy2pdf(device = pdf,
             file = paste(odir,"/Bliss_example.pdf", sep=''),
             width = 14, height = 10, useDingbats = FALSE)










