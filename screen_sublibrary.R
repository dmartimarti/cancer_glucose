# library screening
# This code is meant to analyse the sub-library that Leo created with every
# gene related to glucose in E.coli, and also to use part of our big screen

# MAIN DIRECTORY: "/Users/dmarti14/Documents/MRC_Postdoc/Projects/Leo/screenings"

# libraries
library(tidyverse)
library(readxl)
library(ggrepel)
library(colorspace)
library(Rtsne)
library(here)

# my library of functions
source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')

# session options
options(width = 220)


# creation of directories in which all plots will go
# DO ONLY ONCE

dirs = c('exploration', 'analysis', 'summary')
for (dir in dirs){
	dir.create(dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
}





### load the data

sublib = read_xlsx(here('data','Summary_Keio Sublibrary_Glucose Supp_10mM_N2 worms_08-12-18.xlsx'), 
	sheet = 'Summary') 
sublib[,10:11] = NULL


# transform original data into something usable
sublib = sublib %>%
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
sublib %>%
	write_csv(here('summary','/KeioSublibrary_08-12-18_data_raw.csv'))


# table with genotypes and pathways
pathways = sublib %>% select(Genotype, Pathway) %>% distinct(Genotype, .keep_all = TRUE)
pathways %>%
	write_csv(here('summary','EcoCyc_pathways_KeioSublibrary_08-12-18.csv'))




# load the file with the pathways. Remember there are 3 version, the last one (kegg3)
pathways = read_xlsx(here('data','EcoCyc_pathways_KeioSublibrary_08-12-18.xlsx'), 
	sheet = 'paths')

pathways = pathways %>% 
	mutate(Genotype = as.factor(Genotype),
		   Pathway = as.factor(Pathway),
		   KEGG1 = as.factor(KEGG1),
		   KEGG2 = as.factor(KEGG2),
		   KEGG3 = as.factor(KEGG3))

# summary of the data

sublib.sum <- sublib %>%
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


# range of colours
colours17 = c('#000000','#8C158C','#8C158C','#0A2CFF','#006006',
			'#FF008C','#9D4300','#7B008C','#008BF6','#5EAA44',
			'#FF1500','#A42900','#006A54','#00CDFF','#F0BF8F',
			'#B5B0B8','#FAE603')
colours_kegg2 = c('#F0BF8F','#000000','#8C158C','#0A2CFF','#006006',
				  '#008BF6','#5EAA44','#FF1500','#CC80E6','#BAE303')
colours_kegg3 = c('#8A7090','#000000','#53BAE2','#EF6713','#58A32D','#EF1613','#F45A95','#8C5A3A')


pos = position_jitter(width = 0.05, height = 0.05, seed = 1) # to plot names in jitter positions

# test if the data is OK

sublib.sum %>% 
	left_join(pathways) %>%
	filter(Drug == 5) %>% 
	select(Supplement, Supplement_mM, Genotype, BW_norm, Pathway:KEGG3) %>%
	unite(Supp, Supplement, Supplement_mM) %>%
	spread(Supp, BW_norm) %>%
    ggplot(aes(x = Glucose_0, y = Glucose_10)) + 
	geom_hline(yintercept = 0, colour = 'grey30') +
	geom_vline(xintercept = 0, colour = 'grey30') +
	geom_point(aes(colour = KEGG3), position = pos, size = 2) + 
	scale_color_manual(values = colours_kegg3) + 
    scale_fill_manual(values = colours_kegg3) +
	# geom_text_repel(aes(label = ifelse(Pathway == 'ribonucleotides de novo biosynthesis', as.character(Genotype), '')), position = pos)
	geom_text_repel(aes(label = Genotype, colour = KEGG3), position = pos) +
	labs(title = "5FU + Supplement (Glucose) effect on different BW mutants",
		 x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype', sep = ' ')),
		 y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype ' , bold('(Glucose)'), sep = ' '))) +
	theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
		panel.grid.major = element_line(colour = "grey90"),
		panel.background = element_rect(fill = "white", colour = "grey50")) + 
	guides(colour = guide_legend(override.aes = list(size = 4))) # make lengend points larger


quartz.save(file = here('exploration', 'sublibrary_8categories.pdf'),
	type = 'pdf', dpi = 300, height = 8, width = 10)




###
### load the KEIO library data



# load data
## REMEMBER, PLATE A AND PLATE B ARE NOW 97 AND 99
keiolib = read_xlsx(here('data','Keio Library.xlsx'), sheet = 'ALL_final') 

keiolib = keiolib %>%
	filter(!is.na(Gene), !(Gene %in% c('WT', 'present', 'no bact' )))


## running some checks

# is there any gene repeated? 

subset = keiolib %>% filter(Gene != 'Control')
length(unique(subset$Gene)) == length(subset$Gene)
# if TRUE, everything is fine



# transform data into 1 column data frame
keiolib = keiolib %>% 
	gather(Sample, Score, Control_1:GluFU_3) %>%
	separate(Sample, c('Sample', 'Replicate'), sep = '_') %>%
	mutate_at(c('Plate', 'Column', 'Row', 'Gene', 'Sample'), as.factor) %>%
	mutate(Score = as.numeric(Score))


# Controls median
controls = keiolib %>%
	filter(Gene == 'Control') %>%
	mutate(Score = as.numeric(Score)) %>%
	group_by(Sample) %>%
	summarise(Median = median(Score))

# we can use the values of the first row of the Controls as the median, it has the same values
controls = filter(keiolib, Gene == 'Control', Sample == 'Control')[1,]
controls = rbind(controls, filter(keiolib, Gene == 'Control', Sample == 'FU')[1,])
controls = rbind(controls, filter(keiolib, Gene == 'Control', Sample == 'Glu')[1,])
controls = rbind(controls, filter(keiolib, Gene == 'Control', Sample == 'GluFU')[1,])


## summarise the data
keio.sum = keiolib %>%
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
keio.sum['Supplement_mM'] = 0
keio.sum[keio.sum$Sample == 'Glu',]$Supplement_mM = 10
keio.sum[keio.sum$Sample == 'GluFU',]$Supplement_mM = 10



### test data if everything is OK
pos = position_jitter(width = 0.15, height = 0.15, seed = 1) # to plot names in jitter positions
keio.sum %>% 
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



####
## now load the info from the supplementary material from Leo's previous work


# load data
supp = read_xlsx(here('data','mmc2.xlsx'), sheet = 'Enrichment_with_genes') 

ox_res = supp %>% 
	select(Category, Term_ID, Term, Gene, Plate) %>%
	filter(Term == 'Oxidative phosphorylation')


pur_met = supp %>% 
	select(Category, Term_ID, Term, Gene, Plate) %>%
	filter(Term == 'Purine metabolism')




# is any of these genes in the original sublibrary?
intersect(unique(sublib$Genotype), ox_res$Gene)
# there are four genes found to be repeated in the sublibrary: "sdhA", "sdhB", "sdhC", "sdhD"

# subset the keiolib with the genes from the oxidative phosphorilation

keio.genes = keio.sum %>% 
	mutate(Gene = as.character(Gene)) %>%
	filter(Gene %in% ox_res$Gene) %>% 
	mutate(Gene = recode(Gene, sdhA = 'sdhA_K',
		   					   sdhB = 'sdhB_K',
		   					   sdhC = 'sdhC_K',
		   					   sdhD = 'sdhD_K'), 
		   Drug = 0,
		   Drug = ifelse(Sample %in% c('FU', 'GluFU'), 5, 0),
		   Pathway = 'Oxidative phosphorylation') %>%
	rename(Genotype = Gene) %>%
	select(-Sample)




# intersection between purine metabolism and sublibrary
intersect(unique(sublib$Genotype), pur_met$Gene)

met = intersect(unique(sublib$Genotype), pur_met$Gene)

keio.pur.genes = keio.sum %>% 
	mutate(Gene = as.character(Gene)) %>%
	filter(Gene %in% pur_met$Gene) %>% 
	mutate(Gene = recode(Gene, pykA = 'pykA_K',
		   					   pgm = 'pgm_K',
		   					   nrdF = 'nrdF_K',
		   					   nrdE = 'nrdE_K',
		   					   nrdD = 'nrdD_K',
		   					   mazG = 'mazG_K',
		   					   ndk = 'ndk_K',
		   					   yjjG = 'yjjG_K',
		   					   pykF = 'pykF_K'),
		   Drug = 0,
		   Drug = ifelse(Sample %in% c('FU', 'GluFU'), 5, 0),
		   Pathway = 'Purine metabolism') %>%
	rename(Genotype = Gene) %>%
	select(-Sample)




# join both datasets
join = sublib.sum %>% 
	select(Genotype, Median_Score, Mean, SD, BW_Score, BW_Mean, BW_norm, BW_Mean_norm, Supplement, Supplement_mM, Drug) %>%
	left_join(pathways[,c(1,5)]) %>% 
	rename(Pathway = KEGG3) %>%
	rbind(keio.genes, keio.pur.genes)



joined_colours = c('#8A7090','#000000','#53BAE2','#EF6713','#58A32D','#EF1613','#F45A95','#8C5A3A', 'grey50', '#FF00DC')

join %>% 
	filter(Drug == 5) %>% 
	select(Supplement, Supplement_mM, Genotype, BW_norm, Pathway:Pathway) %>%
	unite(Supp, Supplement, Supplement_mM) %>%
	spread(Supp, BW_norm) %>%
    ggplot(aes(x = Glucose_0, y = Glucose_10)) + 
	geom_hline(yintercept = 0, colour = 'grey30') +
	geom_vline(xintercept = 0, colour = 'grey30') +
	geom_point(aes(colour = Pathway), position = pos, size = 2) + 
	scale_color_manual(values = joined_colours) + 
    scale_fill_manual(values = joined_colours) +
	# geom_text_repel(aes(label = ifelse(Pathway == 'ribonucleotides de novo biosynthesis', as.character(Genotype), '')), position = pos)
	geom_text_repel(aes(label = Genotype, colour = Pathway), position = pos) +
	labs(title = "5FU + Supplement (Glucose) effect on different BW mutants",
		 x = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype', sep = ' ')),
		 y = expression(paste('Normalised median scores of ', italic('C. elegans'), ' phenotype ' , bold('(Glucose)'), sep = ' '))) +
	theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
		panel.grid.major = element_line(colour = "grey90"),
		panel.background = element_rect(fill = "white", colour = "grey50")) + 
	guides(colour = guide_legend(override.aes = list(size = 4))) # make lengend points larger

quartz.save(file = here('exploration', 'sublibrary_extra_pathways.pdf'),
	type = 'pdf', dpi = 300, height = 13, width = 15)






















