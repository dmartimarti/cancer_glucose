#### analysis of data from Leo

#### Index
#### BACTERIAL DATA
#### Load of data
###
#### Summary statistics
###
#### Time series

#### WORM DATA
#### Load of data

#### PCA
###
#### Heatmap
#### Nutritional recovery
###
###



# Load necessary packages
# if something is not loading, try installing package using as in example:
# install.packages('package_name')

library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(PFun)
library(forcats)
library(FactoMineR) # for PCA
library(factoextra) # for PCA


# Set to custom clean theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau

# set working directory to selected directory
# setwd("~/Dropbox/Projects/2015-Metformin/5FU_nutrients/Biolog/")
setwd("/Users/danielmartinez/Google Drive/MRC_Postdoc/scripts/Pov/Leo")

# Set output directory, try creating it
odir <- 'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

# for test plots
tdir = 'Test_plots'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


# If necessary load workspace from previous analysis (had to be saved first)
# load('5FU_Biolog.RData')
# save.image('5FU_Biolog.RData')

########################
#### BACTERIAL DATA ####
########################

#### Load of data

# Info open
# Provide a link to biolog metabolites file. It can be in a different location on your computer
# I have the same file from growth analysis, so if you have all the packages installed and loaded - it will work
# Otherwise, check you assumptions about the directory you are in
# Where the files are, etc.
# I have found that at the beginning it's one of the most common causes of frustration
# Always check you assumptions!
info <- read_csv('Biolog_metabolites.csv',quote = "\"")

info %>%
  filter(Plate == 'PM5')

# Here ~ stands for users home directory
# In each operating system it will be in a different location

# Get bacterial growth summary data
# It should work if you have set your working directory correctly!
# Location relative to the working directory set at the beginning
data.b <- read_csv('Bacteria data/All_files/Data_Leo5/Summary.csv', quote = "\"") %>%
  rename(logAUC_raw = `750nm_f_logAUC`) %>% #`750nm_f_logAUC` data column is what we need for logAUC values
  mutate(Strain = as.character(Strain),
         Strain = recode(Strain, 'BW25113' = 'BW'), #Change strain namings
         ReplicateB = paste0('B', as.character(Replicate)), #Generate replicate ID which indicate that it's data for bacteria
         Type = ifelse(Type == 'Control','C','T'),
         Sample = paste(Strain, Type, sep = '_'),
         Type = factor(Type,
                     levels = c('C', 'T'),
                     labels = c('Control', 'Treatment')),
         Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
         Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
         Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
         Col = factor(Col, levels = LETTERS[1:8]),
         SampleID = paste(Sample, Replicate, sep = '_'),
         Sample = as.factor(Sample),
         Strain = as.factor(Strain)) %>% #Change Type column coding
  group_by(Strain, Type, Plate, Replicate, Group) %>%
  mutate(logAUC = logAUC_raw - logAUC_raw[Metabolite == 'Negative Control']) %>% # Normalisation agains negative control in each plate
  select(SampleID, Sample, Strain, Type, Replicate, ReplicateB, Index, Plate, 
    Well, Row, Col, Metabolite, MetaboliteU, EcoCycID, KEGG_ID, Group, Class, logAUC_raw, logAUC) %>%
  ungroup 
  

# deleting replicate 3 of pyrE treatment (pyrE_T_3), Biospa failed to deliver data correctly 
data.b = data.b %>%
  filter(SampleID != "pyrE_T_3")



head(data.b)


# Error checking
# Plates used
unique(data.b$Plate)

# Strains
unique(data.b$Strain)

# Experiment type
unique(data.b$Type)

# Replicates
unique(data.b$Replicate)

# metabolites
unique(data.b$MetaboliteU)

# Replicates per plate
data.b %>%
  group_by(Strain, Plate, Type) %>% #For each unique combination of values in selected columns
  summarise(Count = n()/96) #Calculate number of rows in a subgroup (corresponds to the 96 wells * number of replicates), to get the number of replicates


#### Summary statistics

# Bacterial data comparison side by side
# Put all bacterial data into one really wide table which would show well values
colnames(data.b)

bac.byrep <- data.b %>%
  unite(STR, Strain, Type, ReplicateB) %>% #Generete unique IDs by concatenating Strain, Type and Replicate values
  # Select only the columns which are necessary: descriptions, values, exclude Replicate (Now it's ReplicateB), Comment, Reader
  select(STR, Plate, Well, Row, Col, Index:Class,logAUC_raw, -Replicate) %>% 
  spread(STR, logAUC_raw) #Spread values into a wide table

# Check colnames of the original table
colnames(data.b)

# Check header of the generated table
head(bac.byrep)
colnames(bac.byrep)



#Replicates per plate
data.breps <- data.b %>%
  group_by(Strain, Plate, Type) %>%
  summarise(Replicates = n()/96)

data.breps
data.breps %>%
  write_csv(paste(odir,'/Bacterial_replicates_total.csv',sep = ''))



head(data.b)



data.b %>%
  group_by(Strain, Plate, Type, Replicate) %>%
  summarise #%>%
#  View

glimpse(data.b)

unique(data.b$SampleID)


# Find duplicated metabolites if any
data.b %>%
  group_by(SampleID, MetaboliteU) %>%
  summarise(Count = n()) %>%
  filter(Count > 1)




# Generate summary statistics for each group
bac.sum <- data.b %>%
  group_by(Strain, Plate, Replicate, Group) %>%
  gather(Measure, Value, logAUC_raw, logAUC) %>%
  group_by(Measure, Strain, Type, Plate, Well, Index, Metabolite, MetaboliteU) %>% #Work on subsets based on selected grouping variables
  summarise(Mean = mean(Value, na.rm = TRUE), # Do selected summary over group items - replicates, ignore missing values
            SD = sd(Value, na.rm = TRUE),
            SE = SD/sqrt(n())) %>% # standard error of the mean
  ungroup %>%
  mutate(PE = Mean+SE,
         NE = Mean-SE)




# With additional info
bac.sum2 <- data.b %>%
  filter(Strain == 'BW' & Plate %in% c('PM1','PM2A', 'PM5')) %>%
  group_by(Strain, Plate, Replicate, Group) %>%
  mutate(logAUC_C = logAUC_raw - logAUC_raw[Metabolite == 'Negative Control' & Type == 'Control']) %>% # Normalisation agains negative control in control
  ungroup %>%
  gather(Measure, Value, logAUC_C, logAUC_raw, logAUC) %>%
  group_by(Measure, Strain, Type, Plate, Well, Index, Metabolite, MetaboliteU) %>% #Work on subsets based on selected grouping variables
  summarise(Mean = mean(Value, na.rm = TRUE), #Do selected summary over group items - replicates, ignore missing values
            SD = sd(Value, na.rm = TRUE),
            SE = SD/sqrt(n())) %>%
  ungroup %>%
  mutate(PE = Mean + SE,
         NE = Mean - SE)




bac.sum.sbs <- bac.sum %>% #Remove table grouping. Can be checked via glimpse(bac.sum)
  filter(Measure == 'logAUC_raw') %>%
  select(Strain, Type, Plate, Well, Index, Mean, SD) %>%
  #Transform the table wide to long, Mean and SD become values of a column named Value, while Mean and SD coding becomes Stat column values
  gather(Stat, Value, Mean, SD) %>% 
  unite(STS, Strain, Type, Stat) %>% # Generate unique IDs for summary columns based on Strain, Type and Stat - statistic values
  spread(STS, Value) # Spread the table to wide, showing statistics in each group


bacall <- bac.byrep %>% 
  left_join(bac.sum.sbs) %>% 
  select(-c(Index, Metabolite, MetaboliteU)) %>% # Join with raw replicate data
  arrange(Plate, Col, Row) %>%
  select(-c('Row', 'Col'))


# Open table to visually inspect
# View(bacall)

# Check column names
# glimpse(bacall)

#Write output
bacall %>%
  write_csv(paste0(odir,'/Biolog_Bacteria_raw_summary.csv'))

### Bacterial data now is ready



# Get timeseries data
# It will take a while
data.b.ts <- read_csv('Bacteria data/All_files/Data_Leo5/Timeseries.csv',quote = "\"") %>%
  filter(Data == '750nm_f') %>% #Select only the relevant data to speed it up
  gather(Time_s, OD, matches('\\d')) %>%
  filter(!is.na(OD)) %>% #Remove empty values if there are missmatches
  mutate(Strain = as.character(Strain),
         Strain = recode(Strain,'BW25113'='BW'), #Change strain namings
         ReplicateB = paste0('B', as.character(Replicate)), #Generate replicate ID which indicate that it's data for bacteria
         Type = ifelse(Type == 'Control','C','T'),
         Type = factor(Type,
                     levels = c('C','T'),
                     labels = c('Control','Treatment')),
         Time_s = as.numeric(Time_s),
         Time_h = Time_s/3600,
         Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
         Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
         Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
         Col = factor(Col, levels = LETTERS[1:8])) #Make them categorical variable with set order them

# [:digit:]{1,} - one or more digits
# This is a 'regular expression'
# More info here https://stringr.tidyverse.org/articles/regular-expressions.html

data.b.ts = data.b.ts %>% mutate(OD = as.numeric(OD))


# Error checking
# Plates used?
unique(data.b.ts$Plate)

# Strains?
unique(data.b.ts$Strain)

# Experiment type?
unique(data.b.ts$Type)

# Replicates
unique(data.b.ts$Replicate)


# Observations per well
data.b.ts %>%
  group_by(Strain,Plate,Type,Replicate) %>%
  summarise(Count = n()/96) %>%
  data.frame
  #View # You pipe into table view function directly!
# Should be interesting to check how many observations you have collected with BioSpa ;)





# Growth curves for summary

tsum <- data.b.ts %>%
  group_by(Strain, Type, Index, Plate, Well, Metabolite, MetaboliteU, Time_h) %>%
  summarise(Mean = mean(OD), SD = sd(OD),SE = SD/sqrt(length(OD))) %>%
  ungroup 

unique(tsum$Type)



Metcols <- c("black", "red")
names(Metcols) <- c("Control","Treatment")
Metlab <- 'Type'



mlevels <- c('Negative Control_C','Uridine','L-Arabinose','Adenosine')

mlabels <- c('NGM','Uridine','L-Arabinose','Adenosine')

# This may take a while
# plot growth curves by different conditions


# tsum %>%
#   filter(Strain == 'BW') %>%
#   filter(MetaboliteU %in% mlevels) %>%
#   mutate(MetaboliteU = factor(MetaboliteU, levels = mlevels, labels = mlabels)) %>%
#   ggplot(aes(x = Time_h, y = Mean, fill = Type, color = Type)) +
#   geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
#   geom_line() +
#   scale_x_continuous(breaks = seq(0, 24, by = 12)) +
#   ylab("OD") +
#   xlab("Time, h") +
#   labs(fill="Type") +
#   #facet_wrap(~MetaboliteU,ncol = 2)+
#   facet_wrap(~MetaboliteU,ncol = 4) +
#   scale_colour_manual(name = Metlab, values = Metcols) +
#   scale_fill_manual(name = Metlab, values = Metcols) +
#   theme(legend.position = "top")

## 

tsum %>%
  filter(MetaboliteU %in% mlevels) %>%
  mutate(MetaboliteU = factor(MetaboliteU, levels = mlevels, labels = mlabels)) %>%
  mutate(Strain = factor(Strain)) %>%
  ggplot( aes(x = Time_h, y = Mean, fill = Type, color = Type)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 12)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Type") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  facet_grid((MetaboliteU ~ Strain)) +
  scale_colour_manual(name = Metlab, values = Metcols) +
  scale_fill_manual(name = Metlab, values = Metcols) +
  theme(legend.position = "top", strip.text = element_text(size = 13))



ggsave(file = paste0(odir,"/Summary_Growth_curves.pdf"),
       width = 110, height = 41, units = 'mm', scale = 2, device = cairo_pdf, family = "Arial")



dwidth <- 1
BWbacsum <- bac.sum2 %>%
  filter(Strain == 'BW' & Measure == 'logAUC_C')

BWbacsum %>%
  filter(MetaboliteU %in% mlevels) %>%
  mutate(MetaboliteU = factor(MetaboliteU, levels = mlevels, labels = mlabels)) %>%
  ggplot(aes(x = MetaboliteU, y = Mean, color = Type)) +
  geom_hline(aes(yintercept = as.double(BWbacsum[BWbacsum$Index == 'PM1-A1' & BWbacsum$Type == 'Control','Mean'])),
    color = 'red4',linetype = 'longdash', alpha = 0.5) +
  geom_hline(aes(yintercept = as.double(BWbacsum[BWbacsum$Index == 'PM1-A1' & BWbacsum$Type == 'Treatment','Mean'])),
    color = 'blue4',linetype = 'longdash', alpha = 0.5) +
  #scale_y_continuous(breaks=seq(-20,20,by=1),limits=c(-4,0))+
  geom_errorbar(aes(ymin = NE, ymax = PE), position = position_dodge(width = dwidth), width = 0.25) +
  geom_point(position = position_dodge(width = dwidth),size = 2) +
  scale_colour_manual(name = Metlab,values = Metcols)+
  labs(color = "Metformin, mM")+
  ylab('Growth AUC vs NGM control, logFC')+
  xlab('Metabolite')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


### FUNCTION TO PLOT YOUR OWN DATA
# This function allow you to plot the strain you want, just pass a variable into the function as a string,
# for example, if you want to plot BW, you should type plot1_fun("BW", c('PM1','PM2A','PM5'))

plot1_fun = function(strain, plate = 'PM1') {

  bac.sum2 <- data.b %>%
  filter(Strain == as.character(strain) & Plate %in% plate) %>%
  group_by(Strain, Plate, Replicate, Group) %>%
  mutate(logAUC_C = logAUC_raw - logAUC_raw[Metabolite == 'Negative Control' & Type == 'Control']) %>% # Normalisation agains negative control in control
  ungroup %>%
  gather(Measure, Value, logAUC_C, logAUC_raw, logAUC) %>%
  group_by(Measure, Strain, Type, Plate, Well, Index, Metabolite, MetaboliteU) %>% #Work on subsets based on selected grouping variables
  summarise(Mean = mean(Value, na.rm = TRUE), #Do selected summary over group items - replicates, ignore missing values
            SD = sd(Value, na.rm = TRUE),
            SE = SD/sqrt(n())) %>%
  ungroup %>%
  mutate(PE = Mean + SE,
         NE = Mean - SE)
  BWbacsum <- bac.sum2 %>%
  filter(Strain == as.character(strain) & Measure == 'logAUC_C')

  BWbacsum %>%
  filter(MetaboliteU %in% mlevels) %>%
  mutate(MetaboliteU = factor(MetaboliteU, levels = mlevels, labels = mlabels)) %>%
  ggplot(aes(x = MetaboliteU, y = Mean, color = Type)) +
  geom_hline(aes(yintercept = as.double(BWbacsum[BWbacsum$Index == 'PM1-A1' & BWbacsum$Type == 'Control','Mean'])),
    color = 'red4',linetype = 'longdash', alpha = 0.5) +
  geom_hline(aes(yintercept = as.double(BWbacsum[BWbacsum$Index == 'PM1-A1' & BWbacsum$Type == 'Treatment','Mean'])),
    color = 'blue4',linetype = 'longdash', alpha = 0.5) +
  #scale_y_continuous(breaks=seq(-20,20,by=1),limits=c(-4,0))+
  geom_errorbar(aes(ymin = NE, ymax = PE), position = position_dodge(width = dwidth), width = 0.25) +
  geom_point(position = position_dodge(width = dwidth),size = 2) +
  scale_colour_manual(name = Metlab,values = Metcols)+
  labs(color = "Metformin, mM")+
  ylab('Growth AUC vs NGM control, logFC')+
  xlab('Metabolite')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
}


# this is for saving the previous plot

ggsave(file = paste0(odir,"/Summary_Ecoli_growth_comparison_large.pdf"),
       width = 110, height = 82, units = 'mm', scale = 2, device = cairo_pdf, family = "Arial")






###################
#### WORM DATA ####
###################


# Let's organise worm data
# Worm data
data.bw <- read_xlsx('Worm_data/Biolog Worm nutrient screen - BW.xlsx', sheet = 'All')
data.crp <- read_xlsx('Worm_data/Biolog Worm nutrient screen - CRP.xlsx',sheet= 'All')
data.pyrE <- read_xlsx('Worm_data/Biolog Worm nutrient screen - pyrE.xlsx',sheet = 'All')
data.TM <- read_xlsx('Worm_data/Biolog Worm nutrient screen - TM.xlsx',sheet = 'All')

# View(data.UV)
# View(data.TM)
# View(data.bw)

### Join all data, create scores, plates, experiments and worm replicate names

data.wa <- rbind(data.bw, data.crp, data.pyrE, data.TM) %>%
  gather(Column, Worm_phenotype,`1`:`12`) %>%
  mutate(Score = as.numeric(ifelse(grepl('[0-9]', Worm_phenotype), as.numeric(Worm_phenotype), NA)),
         Plate = recode(Plate, 'PM2'='PM2A', 'PM3'='PM3B', 'PM4'='PM4A'),
         Replicate = as.character(Replicate),
         ReplicateW = paste0('W',Replicate),
         Experiment = paste0(ifelse(Type == 'Control', 'C','T'),`5FU_uM`)) %>%
  unite(Well, Row, Column, sep = '') %>%
  unite(SER, Strain, Experiment, ReplicateW, remove = FALSE) %>%
  left_join(info) %>%
  select(Strain, Type,`5FU_uM`, Experiment, Replicate, ReplicateW, SER, Plate, Well, Index:Class, Worm_phenotype, Score)

head(data.wa)

# Plates used?
unique(data.wa$Plate)

# Strains?
unique(data.wa$Strain)

# Experiment type?
unique(data.wa$Type)

# 5FU concentration?
unique(data.wa$`5FU_uM`)

# Experiment type?
unique(data.wa$Type)

# Experiment ID?
unique(data.wa$Experiment)

# Are there duplicates?
data.wa[duplicated(data.wa[,c('Strain', 'Type', 'Plate', 'Well', '5FU_uM', 'Replicate')]),]


# All is treatemnt in new dataset
# data.wa<-subset(data.wa,Type=='T')


# Worm raw data side by side comparison
worm.byrep <- data.wa %>%
  #left_join(info[,c('Plate','Well','MetaboliteU')]) %>%
  select(Plate, Well, Index, Metabolite, MetaboliteU, EcoCycID, KEGG_ID, SER, Score) %>%
  spread(SER, Score)
  

head(worm.byrep)




worm.sum <- data.wa %>%
  group_by(Strain, Plate, Well, Experiment) %>%
  summarise(Median = median(Score, na.rm = TRUE),
            Mean = mean(Score, na.rm = TRUE),
            N = sum(!is.na(Score)),
            SD = sd(Score, na.rm = TRUE))
  

head(worm.sum)

# Unique experiments
worm.sum %>%
  group_by(Strain,Experiment) %>%
  summarise(Plate_count = n()/96)



worm.sum %>%
  group_by(Strain,Experiment,Plate) %>%
  summarise(Missing_Median = all(is.na(Median)))




worm.sum %>%
  gather(Stat,Value,Median:SD) %>%
  group_by(Strain,Experiment,Plate,Stat) %>%
  summarise(Missing = all(is.na(Value) | Value == 0) )
  
  



worm.sum.clean <- worm.sum %>%
  gather(Stat, Value, Median:SD) %>%
  unite(SES, Strain, Experiment, Stat) %>%
  filter(!( SES == 'CRP_C0_SD' |
            SES == 'pyrE_C0_SD' 
            ) ) %>%
  select(Plate, Well, SES, Value) %>%
  spread(SES, Value)
         

# worm.cast<-dcast(worm.sum.mc,Plate+Well~Strain+Experiment+Stat,
#                  fun.aggregate = NULL,value.var = 'Value',fill = as.numeric(NA),drop = TRUE)



## this creates a HUGE table with 64 variables 
wormall <- worm.byrep %>%
  left_join(worm.sum.clean)

summary(wormall)

wormall %>%
  write_csv(paste0(odir,'/Biolog_Worms_raw_summary.csv'))








# Data collection and prep complete
# Let's move on to preprocessing
# Take all plates


#& !SampleID=='TM_C_3'
data.clean <- data.b

# Here some samples which are considered deviant (based on PCA results) are removed
data.clean <- data.b %>%
  filter(! (Plate == 'PM3B' & SampleID == 'BW_C_1')) %>%
  filter(! (Plate == 'PM3B' & SampleID == 'TM_C_1')) %>%
  filter(! (Plate == 'PM4A' & SampleID == 'BW_C_1')) %>%
  filter(! (Plate == 'PM4A' & SampleID == 'BW_C_2')) #%>%
  #filter(! (Plate=='PM4A' & SampleID=='BW_T_3') ) %>%
  # filter(! (Plate=='PM4A' & SampleID=='TM_C_3') )



data.brepsc <- data.clean %>%
  group_by(Strain, Plate, Type) %>%
  summarise(Replicates = n()/96)

data.brepsc
data.brepsc %>%
  write_csv(paste(odir,'/Bacterial_replicates_total_clean.csv', sep = ''))



#### PCA plots

#pltsel<-c('All','PM12','PM1','PM2A','PM3B','PM4A','PM5','PM125')

#pltsel<-c('PM12')

# plt<-'PM12'
# 
# plt<-'PM3B'
# 
# plt<-'PM4A'

#Generate infor tables which describes samples
bioinfo <- data.b %>% 
  group_by(SampleID, Sample, Strain, Type) %>%
  summarise %>%
  data.frame

rownames(bioinfo) <- bioinfo$SampleID



fillcols <- c("#FF0000",NA)
names(fillcols) <- c("Control","Treatment")
filllab <- 'Type'

# Definition of PCAplot function
PCAplot <- function(PCAres) {
  ellipses <- PCAres$Ellipses
  pcadata <- PCAres$pcadata
  ggplot(pcadata, aes(x = PC1, y = PC2, colour = Type))+
    xlab(paste('PC1 - ',PCAres$PC1prc,'% of variance',sep=''))+
    ylab(paste('PC2 - ',PCAres$PC2prc,'% of variance',sep=''))+
    geom_path(data = ellipses, aes(x = x, y = y, group = interaction(Sample), linetype = Type), size = 1)+
    geom_point(aes(fill = Type), size = 2, stroke = 1, shape = 21)+
    scale_linetype_manual("Type", values = c("Control" = 1,"Treatment" = 3))+
    scale_fill_manual(name = filllab,values = fillcols)+
    scale_colour_manual(name = Metlab,values = Metcols)+
    guides(linetype = guide_legend(override.aes = list(shape = c(21,21), size = 1,linetype = c(1,3), fill = c(2,NA))))+
    labs(colour = 'Type')+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}


# Plot PCA for selected data
data.b %>%
  filter(Plate %in% c('PM1','PM2A','PM4A', 'PM3B') & Strain == 'BW' & Metabolite != 'Negative Control') %>%
  PCAprep('SampleID', 'MetaboliteU', 'logAUC', bioinfo %>% filter(Strain == 'BW') ) %>% 
  # That's custom function that I wrote. Nothing fancy here, it just facilitates integration with my plotting approaches
  PCAplot


ggsave(file = paste(odir,'/PCA_BW_PM125.pdf',sep = ''),
       width = 55, height = 41, units = 'mm', scale = 2, device = cairo_pdf, family = "Arial")



### this is a more extensive study of PCA analysis
# for this part it is needed:

require(FactoMineR)
library(factoextra)

# pca_b_data = data.b %>%
#   filter(Plate %in% c('PM1','PM2A') & Strain %in% c('BW', 'TM') & Metabolite != 'Negative Control')

# pca_b_data = pca_b_data %>%
#   select(SampleID, MetaboliteU, logAUC) %>%
#   spread(MetaboliteU, logAUC) %>%
#   data.frame(check.names = F)


# rownames(pca_b_data) = pca_b_data[,1]; pca_b_data = pca_b_data[,-1]

# res.pca = PCA(pca_b_data, scale.unit = TRUE, ncp = 5, graph = TRUE)

# meta_var = bioinfo %>% filter(Strain %in% c('BW', 'TM'))

# corrplot(var$cos2, is.corr = FALSE, tl.cex = 0.2)

# # shitty plot
# dev.copy2pdf(device = cairo_pdf,
#              file = paste0(odir,'/PCA_var.pdf'),
#              width = 5, height = 40, useDingbats = FALSE)

# fviz_cos2(res.pca, choice = "var", axes = 1:2)
# dev.copy2pdf(device = cairo_pdf,
#              file = paste0(odir,'/cos2_var.pdf'),
#              width = 15, height = 5, useDingbats = FALSE)


# corrplot(var$contrib, is.corr=FALSE, tl.cex = 0.2) 


# fviz_pca_ind(res.pca,
#              geom.ind = c("text","point"), # show points only (nbut not "text")
#              col.ind = c(meta_var$Strain, meta_var$Type), # color by groups
#              palette = c("#00AFBB", "#E7B800"),
#              addEllipses = TRUE, ellipse.type = "confidence", # Concentration ellipses
#              legend.title = "Groups"
#              )


# PCA Analysis start here!

# filter by strain and Plate if you want
pca_b_data = data.b %>%
  filter(Strain %in% c('pyrE', 'BW') & Metabolite != 'Negative Control')

# some data frame transformations
pca_b_data = pca_b_data %>%
  select(SampleID, MetaboliteU, logAUC) %>%
  spread(MetaboliteU, logAUC) %>%
  data.frame(check.names = F)

rownames(pca_b_data) = pca_b_data[,1]; pca_b_data = pca_b_data[,-1]

# lets compute the PCA
res.pca = PCA(pca_b_data, scale.unit = TRUE, ncp = 5, graph = F)

# metadata 
meta_var = bioinfo %>% filter(Strain %in% c('BW','pyrE'))

# extract info about the individuals
ind = get_pca_ind(res.pca)
ind_df = data.frame(ind$coord[,1], ind$coord[,2], ind$coord[,3], meta_var$Type, 
  meta_var$Strain)

colnames(ind_df) = c('Dim1', 'Dim2', 'Dim3', 'Type', 'Strain')

# get ellipses based on the correlation
getellipse = function(x, y, sc = 1) {
  as.data.frame(ellipse::ellipse(cor(x, y),
                                  scale = c(sd(x) * sc, sd(y) * sc),
                                  centre = c(mean(x), mean(y))))
}

# make a data frame from ellipses
ell = ind_df %>% group_by(Type, Strain) %>% do(getellipse(.$Dim1, .$Dim2, 1)) %>% data.frame

# plot!
ggplot(ind_df, aes(x = Dim1, y = Dim2, color = Strain, group = interaction(Type, Strain))) + 
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell, aes(x = x, y = y, group = interaction(Type, Strain), linetype = Type), size = 1) +
  geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Type, Strain), 
    linetype = Type, fill = Strain), size = 1, alpha = 0.3) +
  xlab(paste("PC1 - ", round(res.pca$eig[1,2], 1), " % of variance", sep = "")) + 
  ylab(paste("PC2 - ", round(res.pca$eig[2,2], 1), " % of variance", sep = "")) +
  ggtitle(paste("PCA of strain:", paste(unique(ell$Strain), collapse = '_'))) + # paste inside a paste to collapse strain names
  theme(plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())


ggsave(file = paste(odir,'/PCA_BW_pyrE.pdf', sep = ''),
       width = 100, height = 80, units = 'mm', scale = 2, device = 'pdf')



# this big and dummy function will plot the PCA of the selected strains and save them in a PDF file
pca_plot = function(strain = "BW"){
  # filter by strain and Plate if you want
  pca_b_data = data.b %>%
    filter(Strain %in% c(strain) & Metabolite != 'Negative Control')

  # some data frame transformations
  pca_b_data = pca_b_data %>%
    select(SampleID, MetaboliteU, logAUC) %>%
    spread(MetaboliteU, logAUC) %>%
    data.frame(check.names = F)

  rownames(pca_b_data) = pca_b_data[,1]; pca_b_data = pca_b_data[,-1]

  # lets compute the PCA
  res.pca = PCA(pca_b_data, scale.unit = TRUE, ncp = 5, graph = F)

  # metadata 
  meta_var = bioinfo %>% filter(Strain %in% c(strain))

  # extract info about the individuals
  ind = get_pca_ind(res.pca)
  ind_df = data.frame(ind$coord[,1], ind$coord[,2], ind$coord[,3], meta_var$Type, 
    meta_var$Strain)

  colnames(ind_df) = c('Dim1', 'Dim2', 'Dim3', 'Type', 'Strain')

  # get ellipses based on the correlation
  getellipse = function(x, y, sc = 1) {
    as.data.frame(ellipse::ellipse(cor(x, y),
                                    scale = c(sd(x) * sc, sd(y) * sc),
                                    centre = c(mean(x), mean(y))))
  }

  # make a data frame from ellipses
  ell = ind_df %>% group_by(Type, Strain) %>% do(getellipse(.$Dim1, .$Dim2, 1)) %>% data.frame

  # plot!
  ggplot(ind_df, aes(x = Dim1, y = Dim2, color = Strain, group = interaction(Type, Strain))) + 
    geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
    geom_path(data = ell, aes(x = x, y = y, group = interaction(Type, Strain), linetype = Type), size = 1) +
    geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Type, Strain), 
      linetype = Type, fill = Strain), size = 1, alpha = 0.3) +
    xlab(paste("PC1 - ", round(res.pca$eig[1,2], 1), " % of variance", sep = "")) + 
    ylab(paste("PC2 - ", round(res.pca$eig[2,2], 1), " % of variance", sep = "")) +
    ggtitle(paste("PCA of strain:", unique(ell$Strain))) +
    theme(plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())

  # save the file
  ggsave(file = here('Summary',paste('PCA_', strain , '.pdf',sep = '')),
  width = 100, height = 80, units = 'mm', scale = 2, device = cairo_pdf, family = "Arial")
}

pca_plot("BW")
pca_plot("TM")
pca_plot("pyrE")
pca_plot(c("crp", "BW", "TM", "pyrE"))
# change name of generated file before next line
pca_plot("crp")


##################
#### Heatmaps ####
##################


hinfo <- data.b %>%
  group_by(SampleID, Strain, Type) %>%
  summarise %>%
  data.frame

rownames(hinfo) <- hinfo$SampleID
hinfo$SampleID <- NULL



cols <- list(Strain = c('BW' = "blue", 'TM' = "brown"),
           Type = c('Treatment' = 'red', 'Control' = 'black'))


oldinfo <- hinfo %>%
  select(-Strain)

# # Complete
# data.b %>%
#   filter(Plate %in% c('PM1', 'PM2A') & Strain == 'TM' & Metabolite != 'Negative Control') %>%
#   HMap('SampleID','MetaboliteU','logAUC', oldinfo, cols = cols) 
#   #Generate heatmap using ggplot. Again, nothing fancy. Code available in PFun R package



# dev.copy2pdf(device = cairo_pdf,
#              file = paste0(odir,'/Heatmap_BW_PM125.pdf'),
#              width = 5, height = 4, useDingbats = FALSE)



#### Heatmap ####

# metadata for columns in Heatmap
hinfo <- data.b %>%
  group_by(SampleID, Strain, Type) %>%
  summarise %>%
  data.frame

rownames(hinfo) <- hinfo$SampleID
hinfo$SampleID <- NULL

heatshape <- data.b %>%
  arrange(Plate) %>%
  select(SampleID, MetaboliteU, logAUC) %>%
  spread(SampleID, logAUC) %>%
  data.frame(check.names = FALSE)
rownames(heatshape) = heatshape[,1]
heatshape = heatshape[,-1] # delete names
colnames(heatshape) = unique(names(heatshape))


hinfo_rows <- data.b %>%
  group_by(MetaboliteU, Plate) %>%
  summarise %>%
  arrange((Plate)) %>%
  data.frame

rownames(hinfo_rows) <- hinfo_rows$MetaboliteU
hinfo_rows$MetaboliteU <- NULL

# ordering rows by plate 
heatshape = heatshape[rownames(hinfo_rows),]

# oldinfo <- hinfo %>%
#   select(-Strain)


cols <- list(Strain = c('BW' = "blue", 'TM' = "brown", 'crp' = "yellow", 'pyrE' = "green"),
           Type = c('Treatment' = 'red', 'Control' = 'black'))

cols2 = list(Plate = c('PM1' = '#003049', 'PM2A' = '#D62828', 'PM3B' = '#F77F00', 'PM4A' = '#FCC459'))


ha = HeatmapAnnotation(df = hinfo, col = cols)
ha2 = HeatmapAnnotation(df = hinfo_rows, col = cols2, which = "row")

# labels = c("BW", paste("BW", italic("gcvP"), sep = ""))

h1 = Heatmap(heatshape, 
 #   col = colorRamp2(c(-2, -1, 0, 2), c("red", "orange","yellow", "green4")), # specifies the color range
    column_title = 'Heatmap', # plot title
    cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", heatshape[i, j]), x, y, gp = gpar(fontsize = 10))}, # plots numbers inside the cell
    column_title_gp = gpar(fontsize = 26, fontface = "bold"), # title options
    name = "log2FC", # name of the legend
    clustering_distance_rows = "euclidean", # column clustering distance
    clustering_distance_columns = "euclidean", # row clustering distance
    clustering_method_rows = "ward.D2", # clustering method
    clustering_method_columns = "ward.D2", 
    # cluster_columns = FALSE, # this disable column order
    # cluster_rows = FALSE,  # this disable row order
    row_dend_reorder = FALSE, 
    top_annotation = ha, # top annotation with two classes 
    top_annotation_height = unit(1, "cm"),
    row_names_side = "left", 
    show_column_names = T,
    column_names_side = "top",
    row_names_gp = gpar(fontsize = 14) # row text size, change this!
    )

# creates a heatmap with row annotations
# HUGE and HEAVY heatmap
h1 + ha2

dev.copy2pdf(device = cairo_pdf,
  file = paste0(odir,"/Heatmap_all.pdf"),
  width = 25, height = 55, useDingbats = FALSE)



# function to plot heatmaps per strain and per plate

heat_plot = function(strain = c('BW','pyrE','crp','TM'), plate = 'PM1', num = FALSE) {

  ## data preparation
  hinfo <- data.b %>%
  filter(Strain %in% strain, Plate %in% plate) %>%
  group_by(SampleID, Strain, Type) %>%
  summarise %>%
  data.frame

  rownames(hinfo) <- hinfo$SampleID
  hinfo$SampleID <- NULL
  
  heatshape <- data.b %>%
    filter(Strain %in% strain, Plate %in% plate) %>%
    arrange(Plate) %>%
    select(SampleID, MetaboliteU, logAUC) %>%
    spread(SampleID, logAUC) %>%
    data.frame(check.names = FALSE)

  rownames(heatshape) = heatshape[,1]
  heatshape = heatshape[,-1] # delete names
  colnames(heatshape) = unique(names(heatshape))
  
  
  hinfo_rows <- data.b %>%
    filter(Strain %in% strain, Plate %in% plate) %>%
    group_by(MetaboliteU, Plate) %>%
    summarise %>%
    arrange((Plate)) %>%
    data.frame
  
  rownames(hinfo_rows) <- hinfo_rows$MetaboliteU
  hinfo_rows$MetaboliteU <- NULL
  
  # ordering rows by plate 
  heatshape = heatshape[rownames(hinfo_rows),]
  
  # oldinfo <- hinfo %>%
  #   select(-Strain)
  
  
  cols <- list(Strain = c('BW' = "blue", 'TM' = "brown", 'crp' = "yellow", 'pyrE' = "green"),
             Type = c('Treatment' = 'red', 'Control' = 'black'))
  
  cols2 = list(Plate = c('PM1' = '#003049', 'PM2A' = '#D62828', 'PM3B' = '#F77F00', 'PM4A' = '#FCC459'))
  
  
  ha = HeatmapAnnotation(df = hinfo, col = cols)
  ha2 = HeatmapAnnotation(df = hinfo_rows, col = cols2, which = "row")
  
  # labels = c("BW", paste("BW", italic("gcvP"), sep = ""))
  
  if (num == FALSE) {
  h1 = Heatmap(heatshape, 
   #   col = colorRamp2(c(-2, -1, 0, 2), c("red", "orange","yellow", "green4")), # specifies the color range
      column_title = paste("Heatmap of plates: ", paste(plate, collapse = ', '), '; and strains: ', paste(strain, collapse = ', ')), # plot title
      # cell_fun = function(j, i, x, y, width, height, fill) {
      # grid.text(sprintf("%.1f", heatshape[i, j]), x, y, gp = gpar(fontsize = 10))}, # plots numbers inside the cell
      column_title_gp = gpar(fontsize = 26, fontface = "bold"), # title options
      name = "log2FC", # name of the legend
      clustering_distance_rows = "euclidean", # column clustering distance
      clustering_distance_columns = "euclidean", # row clustering distance
      clustering_method_rows = "ward.D2", # clustering method
      clustering_method_columns = "ward.D2", 
      # cluster_columns = FALSE, # this disable column order
      # cluster_rows = FALSE,  # this disable row order
      row_dend_reorder = FALSE, 
      top_annotation = ha, # top annotation with two classes 
      top_annotation_height = unit(1, "cm"),
      row_names_side = "left", 
      show_column_names = T,
      column_names_side = "top",
      row_names_gp = gpar(fontsize = 14) # row text size, change this!
      )
  } else {
    h1 = Heatmap(heatshape, 
   #   col = colorRamp2(c(-2, -1, 0, 2), c("red", "orange","yellow", "green4")), # specifies the color range
      column_title = paste("Heatmap of plates: ", paste(plate, collapse = ', '), '; and strains: ', paste(strain, collapse = ', ')), # plot title
      cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.1f", heatshape[i, j]), x, y, gp = gpar(fontsize = 10))}, # plots numbers inside the cell
      column_title_gp = gpar(fontsize = 26, fontface = "bold"), # title options
      name = "log2FC", # name of the legend
      clustering_distance_rows = "euclidean", # column clustering distance
      clustering_distance_columns = "euclidean", # row clustering distance
      clustering_method_rows = "ward.D2", # clustering method
      clustering_method_columns = "ward.D2", 
      # cluster_columns = FALSE, # this disable column order
      # cluster_rows = FALSE,  # this disable row order
      row_dend_reorder = FALSE, 
      top_annotation = ha, # top annotation with two classes 
      top_annotation_height = unit(1, "cm"),
      row_names_side = "left", 
      show_column_names = T,
      column_names_side = "top",
      row_names_gp = gpar(fontsize = 14) # row text size, change this!
      )
  }
  # creates a heatmap with row annotations
  # HUGE and HEAVY heatmap
  h1 + ha2
}

heat_plot(plate = 'PM1')
dev.copy2pdf(device = cairo_pdf, file = paste0(odir,"/Heatmap_PM1.pdf"), width = 25, height = 45, useDingbats = FALSE)
heat_plot(plate = 'PM2A')
dev.copy2pdf(device = cairo_pdf, file = paste0(odir,"/Heatmap_PM2A.pdf"), width = 25, height = 45, useDingbats = FALSE)

heat_plot(plate = 'PM3B')
dev.copy2pdf(device = cairo_pdf, file = paste0(odir,"/Heatmap_PM3B.pdf"), width = 25, height = 45, useDingbats = FALSE)

heat_plot(plate = 'PM4A')
dev.copy2pdf(device = cairo_pdf, file = paste0(odir,"/Heatmap_PM4A.pdf"), width = 25, height = 45, useDingbats = FALSE)






##############################
#### Nutritional recovery ####
##############################


# see what is the mean of logAUC_raw of negative control per plate
# this is important because we use these values to plot the line that
# separates the nutrients that have the higher recovery in mutants

neg_c_sum = data.b %>%
    filter(Metabolite == "Negative Control") %>%
    group_by(Strain, Type, Plate) %>%
    summarise(Mean = mean(logAUC_raw),
        Sd = sd(logAUC_raw)) %>%
    arrange(Strain) %>% data.frame

# dummy function to substract one value to other, only works here
dumm = function(x) {
    return(x[1] - x[2])
}

# creates the intercepts of lines that mark the nutrients
intercepts = data.b %>%
    filter(Metabolite == "Negative Control") %>%
    group_by(Strain, Type, Plate) %>%
    summarise(Mean = mean(logAUC_raw),
        Sd = sd(logAUC_raw)) %>%
    arrange(Plate) %>%
    ungroup %>%
    group_by(Strain, Plate) %>%
    summarise(intercept = dumm(Mean)) %>%
    data.frame

## NegCtrl_BW - NegCtrl_gcvP 


# creating the table with values to plot
# this is a hell of a table
growth = data.b %>%
    unite(STR, Strain, Type, ReplicateB) %>% # Generete unique IDs by concatenating Strain, Type and Replicate values
    select(STR, Plate, Well, Row, Col, Index:Class,logAUC, -Replicate) %>% 
    spread(STR, logAUC)  %>%
    group_by(Plate, Well, Index, Metabolite, MetaboliteU) %>%
    mutate(BW_C_mean = mean(c(BW_Control_B1, BW_Control_B2, BW_Control_B3, BW_Control_B4)),
           BW_C_sd = sd(c(BW_Control_B1, BW_Control_B2, BW_Control_B3, BW_Control_B4)),
           BW_T_mean = mean(c(BW_Treatment_B1, BW_Treatment_B2, BW_Treatment_B3, BW_Treatment_B4)),
           BW_T_sd = sd(c(BW_Treatment_B1, BW_Treatment_B2, BW_Treatment_B3, BW_Treatment_B4)),
           crp_C_mean = mean(c(crp_Control_B1, crp_Control_B2, crp_Control_B3, crp_Control_B4)),
           crp_C_sd = sd(c(crp_Control_B1, crp_Control_B2, crp_Control_B3, crp_Control_B4)),
           crp_T_mean = mean(c(crp_Treatment_B1, crp_Treatment_B2, crp_Treatment_B3, crp_Treatment_B4)),
           crp_T_sd = sd(c(crp_Treatment_B1, crp_Treatment_B2, crp_Treatment_B3, crp_Treatment_B4)),
           pyrE_C_mean = mean(c(pyrE_Control_B1, pyrE_Control_B2, pyrE_Control_B3, pyrE_Control_B4)),
           pyrE_C_sd = sd(c(pyrE_Control_B1, pyrE_Control_B2, pyrE_Control_B3, pyrE_Control_B4)),
           pyrE_T_mean = mean(c(pyrE_Treatment_B1, pyrE_Treatment_B2, pyrE_Treatment_B4)),
           pyrE_T_sd = sd(c(pyrE_Treatment_B1, pyrE_Treatment_B2, pyrE_Treatment_B4)),
           TM_C_mean = mean(c(TM_Control_B1, TM_Control_B2, TM_Control_B3, TM_Control_B4)),
           TM_C_sd = sd(c(TM_Control_B1, TM_Control_B2, TM_Control_B3, TM_Control_B4)),
           TM_T_mean = mean(c(TM_Treatment_B1, TM_Treatment_B2, TM_Treatment_B3, TM_Treatment_B4)),
           TM_T_sd = sd(c(TM_Treatment_B1, TM_Treatment_B2, TM_Treatment_B3, TM_Treatment_B4))) %>% 
    ungroup 



# these are the best nutrients, remember to change plate and value!!!

best1  = growth %>% filter(BW_T_mean   >= (BW_C_mean + 2.068875),   Plate == "PM1") 
best2  = growth %>% filter(pyrE_T_mean >= (pyrE_C_mean + 1.835292), Plate == "PM1")
best3  = growth %>% filter(crp_T_mean  >= (crp_C_mean + 0.456400),  Plate == "PM1")
best4  = growth %>% filter(TM_T_mean   >= (TM_C_mean + 0.182975),   Plate == "PM1")
best5  = growth %>% filter(BW_T_mean   >= (BW_C_mean + 2.423125),   Plate == "PM2A")
best6  = growth %>% filter(pyrE_T_mean >= (pyrE_C_mean + 1.972075), Plate == "PM2A")
best7  = growth %>% filter(crp_T_mean  >= (crp_C_mean + 0.529125),  Plate == "PM2A")
best8  = growth %>% filter(TM_T_mean   >= (TM_C_mean + 0.122100),   Plate == "PM2A")
best9  = growth %>% filter(BW_T_mean   >= (BW_C_mean + 2.713150),   Plate == "PM3B")
best10 = growth %>% filter(pyrE_T_mean >= (pyrE_C_mean + 2.196950), Plate == "PM3B")
best11 = growth %>% filter(crp_T_mean  >= (crp_C_mean + 0.864425),  Plate == "PM3B")
best12 = growth %>% filter(TM_T_mean   >= (TM_C_mean + 0.291775),   Plate == "PM3B")
best13 = growth %>% filter(BW_T_mean   >= (BW_C_mean + 2.708600),   Plate == "PM4A")
best14 = growth %>% filter(pyrE_T_mean >= (pyrE_C_mean + 2.108979), Plate == "PM4A")
best15 = growth %>% filter(crp_T_mean  >= (crp_C_mean + 0.537950),  Plate == "PM4A")
best16 = growth %>% filter(TM_T_mean   >= (TM_C_mean + 0.303275),   Plate == "PM4A")

best_BW = rbind(best1, best5, best9, best13)
best_pyrE = rbind(best2, best6, best10, best14)
best_crp = rbind(best3, best7, best11, best15)
best_TM = rbind(best4, best8, best12, best16)

best_str = c(rep("BW", dim(best_BW)[1]), rep("pyrE", dim(best_pyrE)[1]), rep("crp", dim(best_crp)[1]), rep("TM", dim(best_TM)[1]))

# build a table with the best metabolites, add the strain variable at the end
best = rbind(best_BW, best_pyrE, best_crp, best_TM)
best = cbind(best, best_str); colnames(best)[59] = "Strain"
best = best %>% tbl_df() %>% select(Plate, Strain, Index, MetaboliteU, Metabolite, EcoCycID, KEGG_ID, Group, Class)

#### Scatter plot

# this function takes three arguments: the dataframe with means and std 
# of samples; the strain; the plate. It NEEDS the 'best' variable with the
# best metabolites to write them

scaplot = function(data, strain = "BW", plate = "PM1", lvls = TRUE){

  data = data %>% filter(Plate %in% plate) %>% data.frame

  # selecting index of matrix
  x  =  grep(paste(strain, "_C_mean", sep = ''), colnames(data)) 
  y  =  grep(paste(strain, "_T_mean", sep = ''), colnames(data))
  ysd = grep(paste(strain, "_T_sd", sep = ''), colnames(data))
  xsd = grep(paste(strain, "_C_sd", sep = ''), colnames(data))
  

  # build a data frame
  df = data.frame(data[,x], data[,y], data[,xsd], data[,ysd])
  colnames(df) = c('x', 'y', 'xsd', 'ysd')

  # intercept of blue line
  inter = filter(intercepts, Plate %in% plate, Strain %in% strain)[3]

  # best metabolites in vector mode
  best_m = best %>% filter(Strain %in% strain, Plate %in%plate) %>% 
    select(MetaboliteU) %>% as.matrix %>% as.vector


  # creates an extra column with the metabolite labels to plot
  df$metb = ""
  ix_label = match(best_m, filter(data, Plate %in%  plate)$MetaboliteU)
  df$metb[ix_label] <- best_m
  df$metb = as.factor(df$metb)

  # initiate the plot
  p = ggplot(df, aes(x = x, y = y)) +
    geom_errorbar(aes(ymin = (y - ysd), ymax = (y + ysd)), 
      colour = "grey50", alpha = 0.5) +
    geom_errorbarh(aes(xmin = x - xsd, xmax = x + xsd), 
      colour = "grey50", alpha = 0.5) +
    geom_point(alpha = 0.8, size = 1.5) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "grey50")) +
    labs(
    title = paste(strain, "growth,", plate, "plate", sep = ' '),
    x = paste(strain, "growth", sep = ' '),
    y = paste(strain, "growth + 5FU", sep = ' ')) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", alpha = 0.8) + 
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50", alpha = 0.8) +
          geom_text_repel(aes(label = metb))
    
    if (lvls == TRUE) {
        p + geom_abline(intercept = inter, slope = 1, linetype = "dashed", colour = "blue", alpha = 0.7) 
    } else if (lvls == FALSE){
       p
    } else {
      print('Bad argument to define the lines. Please type TRUE or FALSE')
    }
}

# example of use
p = scaplot(growth, "crp", "PM1") 
p + xlim(-1, 1) + ylim(-0.5, 1.3)



# loop to generate all files
for (i in unique(data.b$Strain)){
  for (j in unique(data.b$Plate)) {
    scaplot(growth, i, j)
    ggsave(file = paste(odir, '/Scatter_', i, '_', j, '.pdf', sep = ''),
    width = 100, height = 80, units = 'mm', scale = 2, device = 'pdf', family = "Arial")
  }
}


# this function draws the average line AFTER the scaplot function has been stored in a 
# variable. 
# WARNING: data here should be the ORIGINAL data, in this case, data.b

avlines = function(data, strain = 'BW', plate = c("PM1", "PM2A", "PM3B", "PM4A"), plot){
  # filter data
  data = data %>% filter(Strain %in% strain, Plate %in% plate, Type == 'Control') %>% 
  select(logAUC) %>% summarise(Max = max(logAUC) + (max(logAUC) * 0.15), Min = min(logAUC) + (min(logAUC) * 0.15))
  
  inter = intercepts %>% filter(Strain == strain) %>% summarise(Mean = mean(intercept), SD = sd(intercept))

  x = seq(data$Min, data$Max, 0.1)
  y = x + inter$Mean
  df = data.frame(x,y)

  plot = plot + geom_ribbon(data = df, aes(ymin = y - inter$SD, ymax = y + inter$SD), fill = 'gray', alpha = 0.2) + 
  geom_line(data = df, aes(x , y), alpha = 0.2)

  return(plot)
}

# example of use

strain = 'pyrE'

strains = c('BW', 'TM', 'crp', 'pyrE')

for (strain in strains){

p = scaplot(growth, strain = strain, plate = c("PM1", "PM2A", "PM3B", "PM4A"), lvls = FALSE)
p = avlines(data.b, strain = strain, plot = p)
p = p +labs(title = paste(strain, ' - all plates', sep = ''))
ggsave(file = paste(odir, '/Scatter_', strain, '_all_plates.pdf', sep = ''),
    width = 100, height = 80, units = 'mm', scale = 2, device = 'pdf')
}




####

# correlation plots
require(corrplot)
require(RColorBrewer)
# generate Pearson correlation coefficients
cor_data = growth %>% 
  select(BW_C_mean, BW_T_mean, crp_C_mean, crp_T_mean, pyrE_C_mean, pyrE_T_mean, TM_C_mean, TM_T_mean) %>%
  cor()

corrplot(cor_data, type = "upper", method = 'square', col = brewer.pal(n = 8, name = "RdBu"), tl.col = "black")

dev.copy2pdf(device = cairo_pdf, file = paste0(odir,"/cor_strains.pdf"), width = 10, height = 10, useDingbats = FALSE)


# statistical test for significance
res1 <- cor.mtest(cor_data, conf.level = .95)

# corrplot with significance
corrplot(cor_data, p.mat = res1$p, insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, 
  type = "upper", method = 'square', col = brewer.pal(n = 8, name = "RdBu"), tl.col = "black")

dev.copy2pdf(device = cairo_pdf, file = paste0(odir,"/cor_strains_sig.pdf"), width = 10, height = 10, useDingbats = FALSE)






####
# linear modeling



contrasts <- read.contrasts('!Contrasts.xlsx') # This table describes comparisons made in data

# loads the metadata table
contrasts.table <- contrasts$Contrasts.table
# loads the contrasts matrix (1s and 0s)
contr.matrix <- contrasts$Contrasts.matrix

# shortens the metadata data frame
contr.desc <- contrasts.table %>%
  select(Description:Strain)



contr.matrix

# Initialise adjustments
group <- c('All','PM4A','PM3B','PM2A','PM1')
strain <- c('BW','TM','crp','pyrE')


adjustments <- expand.grid('Group' = group, 'Strain' = strain) %>%
  mutate(a = -1, b = 0,
         Contrast = paste0(Strain, '_5FU_', Group, 'n')) %>%
  select(Contrast, Group:b) 



# transform the contr.matrix into a real matrix (not data frame)
contr.matrix <- contr.matrix %>%
  as.matrix()


# this piece of code performs a lot of different things, be careful

allresults <- data.b %>%
  #filter(Strain == 'BW') %>%
  filter(!Metabolite %in% c('Negative Control','Positive Control') ) %>%
  # filter(! (Plate == 'PM3B' & SampleID == 'BW_C_1') ) %>%
  # filter(! (Plate == 'PM4A' & SampleID == 'TM_C_3') ) %>%
  group_by(Plate, Well, Index, Metabolite, MetaboliteU, EcoCycID, KEGG_ID) %>%
  # Another creation of mine. Nothing fancy, just a wrapper for the use of contrasts in LM. 
  # But integrates nicely with tidyverse workflow. Code in PFun
  do(hypothesise(.,"logAUC~0+Sample", contr.matrix)) %>% 
  group_by(Plate, Well, Metabolite, MetaboliteU, EcoCycID, KEGG_ID) %>%
  getresults(contr.desc) # Calculate some additional parameters. Code in PFun

#get results in different shapes
results <- allresults$results
results.cast <- allresults$cast
results.castfull <- allresults$castfull

results$Contrast %>% unique %>% as.character



# metrem <- c('2-Hydroxy Benzoic Acid','L-Leucine')

# Multiplex is another handy tool I have written. It allows stacking table with itself side_by_side as many times as you wish.
# The generated table can map any data entry to ggplot
# In such framework, you can analyse/plot any relationship between data entries in the table
# Integrates very well with tidyverse framework


results.cor <- multiplex(results,c("Strain", "Plate", "Well", "Index", "Metabolite", "MetaboliteU", "EcoCycID", "KEGG_ID"), 2) %>% 
  filter(str_detect(x_Contrast,'_C') &
           str_detect(y_Contrast,'_T'))
            # &
           # !Metabolite %in% metrem )



results.cor$x_Contrast %>% unique
results.cor$y_Contrast %>% unique

glimpse(results.cor)

# results.cor %>%
#   filter(Metabolite %in% metrem)



# Some testing

results.cor %>%
  filter(Strain == 'BW') %>%
  filter(Plate %in% c('PM1','PM2A', 'PM3B', 'PM4A')) %>%
  lm(y_logFC~x_logFC,data=.) %>%
  summary

results.cor %>%
  filter(Strain == 'TM') %>%
  filter(Plate %in% c('PM1','PM2A', 'PM3B', 'PM4A')) %>%
  lm(y_logFC~x_logFC,data=.) %>%
  summary

results.cor %>%
  filter(Strain == 'pyrE') %>%
  filter(Plate %in% c('PM1','PM2A', 'PM3B', 'PM4A')) %>%
  lm(y_logFC~x_logFC,data=.) %>%
  summary

results.cor %>%
  filter(Strain == 'crp') %>%
  filter(Plate %in% c('PM1','PM2A', 'PM3B', 'PM4A')) %>%
  lm(y_logFC~x_logFC,data=.) %>%
  summary

results.cor %>%
  filter(Strain %in% c('BW', 'TM')) %>%
  filter(Plate %in% c('PM1','PM2A', 'PM3B', 'PM4A')) %>%
  lm(y_logFC~x_logFC,data=.) %>%
  summary

adjustments

dim(results.cor)

# Step No 2
# Find contrast adjustments
confounding <- function(adjustments, data) {
  
  grp <- as.character(adjustments$Group)
  strn <<- as.character(adjustments$Strain)


  print(paste0("Strain: ",strn) )
  print(paste0("Group: ",grp) )
  
  seldata<-data %>%
    filter(Strain==strn)
  
  if (grp=='All'){
    seldata<-seldata
  } else if (grp=='PM1'){
    seldata<-seldata %>%
      filter(Plate %in% c('PM1'))
  } else if (grp=='PM3B') {
    seldata<-seldata %>%
      filter(Plate %in% c('PM3B'))
  } else if (grp=='PM4A') {
    seldata<-seldata %>%
      filter(Plate %in% c('PM4A'))
  } else if (grp=='PM2A') {
    seldata<-seldata %>%
      filter(Plate %in% c('PM2A'))
  } 

  print(paste0('Plates: ', paste0(seldata$Plate %>% unique %>% print,collapse=';') ) )
  print(paste0('Rows: ',dim(seldata)[1]) )
  
  
  if (dim(seldata)[1] == 0){
    result<-data.frame('a' = NA,
                       'neg_a' = NA,
                       'b' = NA,
                       'a_p' = NA,
                       'b_p' = NA,
                       'r2' = NA)
    return(result)
  }

  print('First round')
  
  # print(unique(seldata$x_Contrast))
  # print(unique(seldata$y_Contrast))
  
  # frml<-paste0(strn,"_T_logFC~",strn,"_C_logFC")
  # print(paste0("Formula: ",frml) )
  
  #genfit<-lm(seldata[,paste0(strn,"_T_logFC")]~seldata[,paste0(strn,"_C_logFC")])
  
  #attach(seldata)
  attach(seldata)
  genfit <- lm("y_logFC~x_logFC")
  
  genres <- summary(genfit)
  
  car::outlierTest(genfit) # Automatically remove outliers, which could skew the correlation. Especially valuable in the case of 5-FU
  ot <- car::outlierTest(genfit)
  car::qqPlot(genfit, main = "QQ Plot")
  outlist<-c(names(ot$rstudent))
  #Outliers:
  print(paste('Outliers in ', grp, sep = ''))
  print(seldata[outlist,])
  #
  #
  filtseldata <- seldata[setdiff(rownames(seldata),outlist),]
  
  detach(seldata)

  print('Second round')
  # 
  #Get final fit
  attach(filtseldata)
  genfit2 <- lm("y_logFC~x_logFC")
  genres <- summary(genfit2)
  
  result <- data.frame('a' = genres$coefficients[[2]],
                     'neg_a' = -genres$coefficients[[2]],
                     'b' = genres$coefficients[[1]],
                     'a_p' = genres$coefficients[[8]],
                     'b_p' = genres$coefficients[[7]],
                     'r2' = genres$r.squared)
  

  return(result)

}

# Update adjustments based on identified correlations. Go back to line 806 and rerun linear modeling with adjustments
adjustments <- adjustments %>%
  group_by(Contrast,Group,Strain) %>%
  do(confounding(adjustments = ., data = results.cor)) #



adjustments


# Save adjustments for future reference
adjustments %>%
  write_csv(paste0(odir,'/Contrast_adjustments.csv'))




# Update this based on identified adjustment values
# Only run this after the first round of analysis
#######################################


contr.matrix <- contr.matrix %>%
  data.frame()


contr.matrix['BW_5FU_Alln',c('BW_C','m')] <- adjustments[adjustments$Contrast ==  'BW_5FU_Alln', c('neg_a','b')]
contr.matrix['TM_5FU_Alln',c('TM_C','m')] <- adjustments[adjustments$Contrast ==  'TM_5FU_Alln', c('neg_a','b')]
contr.matrix['crp_5FU_Alln',c('crp_C','m')] <- adjustments[adjustments$Contrast == 'crp_5FU_Alln', c('neg_a','b')]
contr.matrix['pyrE_5FU_Alln',c('pyrE_C','m')]  <- adjustments[adjustments$Contrast == 'pyrE_5FU_Alln', c('neg_a','b')]

adjustments

contr.matrix

#######################################


#### new round of statistical analyses
contr.matrix <- contr.matrix %>%
  as.matrix()
# this piece of code performs a lot of different things, be careful

allresults <- data.b %>%
  #filter(Strain == 'BW') %>%
  filter(!Metabolite %in% c('Negative Control','Positive Control') ) %>%
  # filter(! (Plate == 'PM3B' & SampleID == 'BW_C_1') ) %>%
  # filter(! (Plate == 'PM4A' & SampleID == 'TM_C_3') ) %>%
  group_by(Plate, Well, Index, Metabolite, MetaboliteU, EcoCycID, KEGG_ID) %>%
  # Another creation of mine. Nothing fancy, just a wrapper for the use of contrasts in LM. 
  # But integrates nicely with tidyverse workflow. Code in PFun
  do(hypothesise(.,"logAUC~0+Sample", contr.matrix)) %>% 
  group_by(Plate, Well, Metabolite, MetaboliteU, EcoCycID, KEGG_ID) %>%
  getresults(contr.desc) # Calculate some additional parameters. Code in PFun



#get results in different shapes
results <- allresults$results
results.cast <- allresults$cast
results.castfull <- allresults$castfull

results$Contrast %>% unique %>% as.character



## some tests with results2 resulting from analysis with adjustments

cosa1 = results %>%
  filter(FDR < 0.05 & Contrast == "TM_5FU_Alln") %>%
  arrange(FDR)
cosa2 = results2 %>%
  filter(FDR < 0.05 & Contrast == "TM_5FU_Alln") %>%
  arrange(FDR) 




results %>%
  filter(Contrast %in% c('BW_C', 'BW_T', 'BW_5FU') & Metabolite %in% c('Uridine'))
results2 %>%
  filter(Contrast %in% c('BW_C', 'BW_T', 'BW_5FU') & Metabolite %in% c('Uridine'))


# Save statistical analysis results
write.csv(results,paste0(odir,'/Ecoli_results.csv'),row.names = FALSE)
write.csv(results.cast,paste(odir,'/Ecoli_results_sidebyside.csv', sep = ''), row.names = FALSE)
write.csv(results.castfull,paste(odir,'/Ecoli_results_sidebyside_full.csv',sep=''), row.names = FALSE)



#### Boxplots

# joining data

box.data = growth %>% 
  filter(Metabolite != 'Negative Control') %>%
  left_join(results.cast)

# look how many variables, it's a mess
names(box.data)

data.b %>%
  filter(MetaboliteU == 'Methylamine', Strain %in% c('BW', 'pyrE')) %>%
  ggplot(aes(x = Strain, y = logAUC)) + 
  geom_boxplot(aes(fill = Type), position = position_dodge(0.9), outlier.color = 'red') +
  scale_fill_manual(values = c("#999999", "#E69F00"))




# box plots of all Glucose nutrients 

metabolites = growth$MetaboliteU
glu_nut = metabolites[grep('Glucose', metabolites)] 

data.b %>%
  filter(MetaboliteU %in% glu_nut, Strain %in% c('BW', 'pyrE', 'crp')) %>%
  ggplot(aes(x = Sample, y = logAUC)) + 
  geom_boxplot(aes(fill = Type), position = position_dodge(0.9), outlier.color = 'red') +
  labs(title = "Boxplots of logFC",
  x = "Strain", y = 'logAUC' ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 7)) +
  facet_wrap(~MetaboliteU, ncol = 4) +
  scale_fill_manual(values = c("#CE372F", "#00A8E8"))


ggsave(file = paste(odir, '/boxplot_glucose_BWvspyrE_.pdf', sep = ''),
    width = 100, height = 100, units = 'mm', scale = 2, device = 'pdf')


#### time series curves for glucose compounds



# Growth curves for summary

tsum <- data.b.ts %>%
  group_by(Strain, Type, Index, Plate, Well, Metabolite, MetaboliteU, Time_h) %>%
  summarise(Mean = mean(OD), SD = sd(OD),SE = SD/sqrt(length(OD))) %>%
  ungroup 

Metcols <- c("black", "red")
names(Metcols) <- c("Control","Treatment")
Metlab <- 'Type'

mlevels <- glu_nut
mlabels <- glu_nut

tsum %>%
  filter(MetaboliteU %in% glu_nut, Strain %in% c('BW', 'pyrE')) %>%
  mutate(MetaboliteU = factor(MetaboliteU, levels = glu_nut, labels = mlabels)) %>%
  mutate(Strain = factor(Strain)) %>%
  ggplot( aes(x = Time_h, y = Mean, fill = Type, color = Type)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 12)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Type") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  facet_grid((MetaboliteU ~ Strain)) +
  scale_colour_manual(name = Metlab, values = Metcols) +
  scale_fill_manual(name = Metlab, values = Metcols) +
  theme(legend.position = "top", strip.text = element_text(size = 8))



ggsave(file = paste0(odir,"/glucose_Growth_curves.pdf"),
       width = 80, height = 200, units = 'mm', scale = 2, device = 'pdf')












##################
### 4-way plot ###
##################


worm.sum[worm.sum$Strain == 'CRP', ]$Strain = 'crp'
worm.sum$Strain = as.factor(worm.sum$Strain)

worm.sum = worm.sum %>% mutate(
              Median = ifelse(Median == 0, 1, Median))

unique(worm.sum$Strain)
unique(worm.sum$Experiment)

for (name in unique(worm.sum$Strain)){
  print(name)
  print(unique((worm.sum %>% filter(Strain == name))$Experiment))
}
# select only one experiment and strain
# contrasts and other main variables of this part

strain = 'crp'
experiment = "T86"
str_C = paste(strain, '_C', sep = '')
str_T = paste(strain, '_T', sep = '')
str_CT = paste(strain, '_5FU', sep = '')
alln = paste(strain, '_5FU_Alln', sep = '')
adjs = paste(strain, '_5FU_adj', sep = '')


# for example, to mimic Pov's code
worm.sum %>% filter(Experiment == experiment, Strain == strain)

worm.old = worm.sum %>%
  filter(Experiment == experiment & Strain == strain) %>% 
  ungroup %>%
  mutate(Description = 'C. elegans development at 5uM 5-FU',
         Contrast = 'Ce_Dev5',   # C. elegans developement with 5FU
         Contrast_type = 'Treatment',
         FDR = NA) %>%
  select(Description, Contrast, Contrast_type, Plate, Well, logFC = Median, SE = SD, FDR)

# check that we filtered stuff properly
worm.old %>%
  group_by(Plate) %>%
  summarise(N = n())


# Join bacterial and worm results (change contrasts)
jointresults = allresults$results %>%
  filter((Contrast %in% c(str_C, str_T, str_CT, alln))) %>%
  mutate(Contrast = fct_recode(Contrast, adjs = alln)) %>% 
  select(Description, Contrast, Contrast_type, Plate, Well, logFC, SE, FDR) %>%
  bind_rows(worm.old) %>%
  left_join(info) %>%
  select(Description:Well,Index:KEGG_ID,logFC:FDR)

### some tests 
jointresults %>%
  filter(FDR < 0.05 & Contrast == adjs) %>%
  arrange(FDR)
  #  %>%
  # View

jointresults %>%
  mutate(Direction = ifelse(logFC > 0,'Up','Down'))%>%
  filter(FDR < 0.05 & Contrast == adjs) %>%
  group_by(Direction) %>%
  summarise(N = n()) 
  
  
# should be 384
jointresults %>%
  filter(Contrast == 'Ce_Dev5')



jointcast = jointresults %>%
  select(Contrast, Plate, Well, Index, Metabolite, MetaboliteU, EcoCycID, KEGG_ID, logFC, SE, FDR) %>%
  gather(Stat, Value, logFC, SE, FDR) %>%
  unite(CS, Contrast, Stat) %>%
  spread(CS, Value) %>%
  rename(Ce_Dev5_Median = Ce_Dev5_logFC, Ce_Dev5_SD = Ce_Dev5_SE) %>%
  select(-Ce_Dev5_FDR)


# write.csv(jointresults,paste(odir,'/Ecoli_results_All_As_old_screen.csv', sep = ''),row.names = FALSE)
# write.csv(jointcast,paste(odir,'/Ecoli_results_sidebyside_As_old_screen.csv', sep = ''),row.names = FALSE)

 
# Multiplex 
jointresults.multi <- multiplex(jointresults,c("Plate","Well","Index","Metabolite","MetaboliteU","EcoCycID","KEGG_ID"))
jointresults.multi2 <- multiplex(jointresults,c("Plate","Well","Index","Metabolite","MetaboliteU","EcoCycID","KEGG_ID"),2)



# more tests
jointresults.multi2 %>%
  filter(x_Contrast == adjs & y_Contrast == 'Ce_Dev5') %>%
  filter(x_FDR < 0.05 & x_logFC > 0) %>%
  filter(y_logFC > 1)


#Up 13, with 8 Ce, 5 unique
#Down 28 with 6 Ce, 22 unique


#Ce 44-8-6=30 unique

#44


jointresults %>%
  filter(Metabolite %in% mrknames & Contrast %in% c(adjs))





########################
# Old screen data PM1,PM2A,PM5 as it was planned for Scott et al paper
# Generated for thesis
NGMa <- adjustments[adjustments$Contrast == alln, ]$a
NGMb <- adjustments[adjustments$Contrast == alln, ]$b

NGMa
NGMb

# adjustments


adjustments[adjustments$Contrast == alln,]

# gradcolours<-c('black','yellow','orange','red')
# gradcolours2<-c('blue','cyan','black','yellow','red')
#gradcolours = c('#000000', '#4E4972', '#63D471', '#B55F00', '#C21C22')


# some tests

# jointresults %>%
#   filter(Contrast=='BW_5FU_adj' & Metabolite %in% mrknames) %>%
#   arrange(FDR)


# jointresults %>%
#   filter(Contrast=='BW_5FU_adj' & logFC>0 & FDR<0.05) %>%
#   arrange(FDR)



# jointcast %>%
#   filter(BW_C_logFC>1.25 & BW_T_logFC<0 & Ce_Dev5_Median==1)

gradcolours = c('#71B83B','yellow','orange','red')


jointresults.multi %>%
  filter(x_Contrast == str_C & y_Contrast == str_T & z_Contrast == 'Ce_Dev5') %>%
  filter(z_logFC != is.na(z_logFC)) %>%
  # mutate(z_logFC = 4) %>%
  ggplot(aes(x = x_logFC, y = y_logFC)) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.3, color = 'grey', linetype = 'longdash') +
  geom_abline(aes(intercept = NGMb, slope = NGMa), alpha = 0.6, color = 'red') +
  geom_vline(aes(xintercept = 0), alpha = 0.9, color = 'grey')+
  geom_hline(aes(yintercept = 0), alpha = 0.9, color = 'grey')+
  geom_errorbarh(aes(xmax = x_logFC + x_SE, xmin = x_logFC - x_SE), height = 0, alpha = 0.3, color = 'grey50') +
  geom_errorbar(aes(ymax = y_logFC + y_SE, ymin = y_logFC - y_SE), width = 0, alpha = 0.3, color = 'grey50') +
  geom_point(aes(colour = z_logFC), alpha = 0.9, size = 2) +
  scale_size(range = c(1,5), name = 'C. elegans\nphenotype') +
  scale_colour_gradientn(colours = gradcolours,
                         breaks = c(1,2,3,4), limits = c(1,4), guide = "legend", name = 'C. elegans\nphenotype') +
  scale_y_continuous(breaks = -5:5)+
  coord_cartesian(xlim = c(-4, 3), ylim = c(-3, 4)) +
  labs(title = paste(strain, " growth with 5FU treatment at 86 uM", sep = ''),
  x = "E. coli growth vs NGM - Control, logFC", 
  y = 'E. coli growth vs NGM - 5-FU Treatment, logFC' ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # geom_text_repel(aes(label = ifelse(z_logFC >= 4, MetaboliteU, '')), box.padding = unit(0.6, "lines"), segment.alpha = 0.4) + # only name those nutr with C. elegans phenotype 3 or more
  geom_text_repel(aes(label = ifelse(z_logFC >= 3 | z_logFC == 0, MetaboliteU, '')), box.padding = unit(0.6, "lines"), segment.alpha = 0.4) + # only name those nutr with C. elegans phenotype 3 or more
  guides(color = guide_legend()) +
  theme(panel.grid.major = element_line(size = 0.06, colour = "grey50"),
      # panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.background = element_rect(fill = "white", colour = "grey50"),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  annotate(geom = 'text', x = -2, y = 1.5, label = paste('y = ', round(NGMa,3), 'x +', round(NGMb,3), sep = ''), size = 5)

ggsave(file = paste0(odir,"/Scatter_", strain, "_86uM_screen.pdf"),
       width = 120, height = 100, units = 'mm', scale = 2, device = 'pdf')






###########################
### Profile plots of BW ###
###########################


# let's gather the data. we should have a table with nutrients as rows, and 
# the median/mean of the worm phenotype as columns

# select the data from BW
prof = worm.sum %>%
  filter(Strain == 'BW')


# join data with metabolite info

joinprof = data.b %>%
  filter(Strain == 'BW', Type == 'Control', Replicate == 1) %>%
  select(Well, Plate, MetaboliteU) %>%
  left_join(prof)

# fix one point which Median is 3.25
joinprof[joinprof$MetaboliteU == 'Ala-Leu' & joinprof$Experiment == 'T1',]$Median = 3

# check that everything is ok
unique(joinprof$MetaboliteU)

# plot!
joinprof %>% 
  ggplot(aes(x = Experiment, y = Median, group = MetaboliteU)) +
  geom_line(alpha = 0.1, size = 2, colour = "black") + 
  geom_jitter(aes(colour = Experiment), width = 0.1, height = 0.1, alpha = 0.8) +
  labs(
  title = 'C. elegans phenotype transition',
  x = expression(paste("Concentration of 5FU, in ", mu, "M", sep = '')), 
  y = 'Worm phenotype (median)' ) +
  theme(
  plot.title = element_text(hjust = 0.5, face = "bold"),
  panel.grid.major = element_line(size = 0.06, colour = "grey50"),
  panel.grid.minor.y = element_blank(),
  panel.background = element_rect(fill = "white", colour = "grey50")) +
  scale_x_discrete(labels = c("T1" = '1', "T1.5" = "1.5", "T5" = "5"))

ggsave(file = paste0(odir,"/worm_phenotype_profile_BW.pdf"),
       width = 160, height = 100, units = 'mm', scale = 2, device = 'pdf')



##########################
### Development assay  ###
##########################

# read table
dev = read_xlsx('Develop_data/test_table.xlsx', sheet = 'Summary')

# data table creation
dev = dev %>% 
  mutate(Supplement = as.factor(Supplement),
         Genotype = as.factor(Genotype)) %>% 
  gather(Drug, Score, `0`, `1`, `2.5`, `5`) %>%
  mutate(Drug = as.factor(Drug))

# compute some summary statistics
dev.sum = dev %>% 
  group_by(Supplement, Genotype) %>%
  summarise(Mean = mean(Score, na.rm = TRUE), SD = sd(Score, na.rm = TRUE), Median = median(Score, na.rm = TRUE))

dev %>% 
  filter(Genotype %in% c('BW25113', 'ΔacnB::K')) %>%
  ggplot(aes(x = Drug, y = Score, fill = Genotype)) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.height = 0.1)) +
  geom_boxplot(position = position_dodge(0.9), na.rm = TRUE) +
  facet_wrap(~Supplement) +
  scale_color_manual(values = c("#CE372F", "#00A8E8")) + 
  scale_fill_manual(values = c("#CE372F", "#00A8E8")) +
  labs(title = "Worm phenotype",
  x = "Concentration of 5FU", y = 'Worm Score' ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 7))


# multiplot with single genotypes

# lets reorder the data
genotypes = unique(dev$Genotype)

gen_order = c('BW25113', 'ΔgltA::K', 'ΔacnB::K', 'ΔgltA ΔacnB::K', 'Δacs::K', 'ΔgltA Δacs::K', 'ΔprpB::K', 'Δglta ΔprpB::K',
  'ΔprpC::K', 'ΔgltA ΔprpC::K', 'ΔprpD::K', 'ΔgltA ΔprpD::K', 'ΔprpE::K', 'ΔgltA ΔprpE::K')

dev$Genotype = factor(dev$Genotype, levels = gen_order)

# set colours
colours = terrain.colors(14)
colours = c('#403F4C', '#3185FC', '#F3A712', '#F2D921', '#F0CEA0', '#965A6D', '#245428', '#243228',
 '#276FBF', '#00489B', '#420039', '#911580', '#635255', '#362C2E')

dev %>% 
  ggplot(aes(x = Drug, y = Score, fill = Genotype)) +
  geom_boxplot(position = position_dodge(0.9), na.rm = TRUE) +
  geom_point(pch = 21, size = 2, position = position_jitterdodge(jitter.height = 0.1), na.rm = TRUE) +
  facet_wrap(vars(Genotype, Supplement), ncol = 4, nrow = 7,  labeller = labeller(.multi_line = FALSE)) +
  scale_color_manual(values = colours) + 
  scale_fill_manual(values = colours) +
  labs(title = "Worm phenotype",
  x = "Concentration of 5FU", y = 'Worm Score' ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 7),
    strip.background.x = element_rect(colour = "black", fill = "white"))

ggsave(file = paste0(odir,"/Worm_phenotype_gltA.pdf"),
       width = 170, height = 150, units = 'mm', scale = 2, device = cairo_pdf)



##### Balloon plot

# shape data

dev %>% spread(Supplement, Drug, Score)











#############
### Worms ###
#############

worm.hdes <- data.wa %>%
  group_by(SER, Strain, Type, `5FU_uM`) %>%
  summarise %>%
  data.frame


wormshape <- worm.byrep %>%
  filter(Metabolite != 'Negative Control') %>%
  select(MetaboliteU, contains('BW'), contains('TM'), contains('crp'), contains('pyrE')) %>%
  data.frame(check.names = FALSE)

rownames(wormshape) <- wormshape$MetaboliteU
wormshape$MetaboliteU <- NULL

rownames(worm.hdes) <- worm.hdes$SER
worm.hdes$SER <- NULL


whanot <- worm.hdes[colnames(wormshape),]
wha <- HeatmapAnnotation(df = whanot, col = list(Strain = c('BW' = "green", 'TM' = "red",'crp' = 'blue', 'pyrE' = 'purple'),
      Type = c('Control'='black','Treatment'='white'),
      X5FU_uM = colorRamp2(c(0,20), c("white","red4")) ))
head(wormshape)
# make annotated heatmap
# This owul not work with ggplot, so I'm using Heatmap from ComplexHeatmap package
Heatmap(wormshape, name = "Phenotype",
            col = colorRamp2(c(0,2,5), c("black",'orange',"red")),
            column_names_side = 'top',
            clustering_method_rows = 'ward.D2',
            #clustering_method_columns ='ward.D2',
            top_annotation = wha,
            row_dend_reorder = TRUE,
            column_dend_reorder = TRUE,
            row_names_max_width = unit(10, "cm"))


dev.copy2pdf(device = 'pdf',
             file = paste0(odir,"/Heatmap_Worms.pdf"),
             width = 15, height = 60, useDingbats = FALSE)






#Join bacterial and worm results to one massive table
jointresults<-allresults$results %>%
  filter( (Contrast %in% c('BW_C','BW_T','BW_5FU','BW_5FU_PM125n') & 
    Plate %in% c('PM1','PM2A','PM5')) |  ( (Contrast=="BW_5FU_PM12n" & Plate %in% c('PM1','PM2A')) | 
    (Contrast=="BW_5FU_PM5n" & Plate %in% c('PM5')) ) ) %>%
  # In this case plate specific adjustments were used due to distinct data stracture. In the future, with more data, global adjustment should be used.
  mutate(Contrast=fct_recode(Contrast,'BW_5FU_adj'='BW_5FU_PM12n','BW_5FU_adj'='BW_5FU_PM5n')) %>% 
  select(Description,Contrast,Contrast_type,Plate,Well,logFC,SE,FDR) %>%
  bind_rows(worm.old) %>%
  left_join(info) %>%
  select(Description:Well,Index:KEGG_ID,logFC:FDR)





# Summaries and comparisons

head(data.b)

#### this is not working, undefined MInMeanSDMax variable ####

data.b %>%
  ggplot(aes(x = Strain, y = logAUC, color = Type))+
  stat_summary(fun.data = MinMeanSDMMax, geom = "boxplot", position = "identity") +
  geom_point() +
  #scale_y_continuous(breaks=seq(-20,20,by=1))+
  geom_text(aes(label = Replicate), color = 'black', size = 2)+
  ylab('Normalised log2 AUC, OD*h') +
  xlab('Strain') +
  labs(color = '5-FU') +
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~MetaboliteU, ncol = 8, scales = 'free_y')

dev.copy2pdf(device = cairo_pdf,
             file = paste(odir,"/Raw_logAUC_by_metabolite.pdf",sep=''),
             width = 15,height=80, useDingbats=FALSE)




# careful, heavy plot here!
# Histograms
data.b %>%
  ggplot(aes(x = logAUC))+
  ggtitle('Distribution of normalised log2 AUC')+
  geom_density(aes(y = ..scaled..),fill = 'red', alpha = 0.5)+
  geom_rug(aes(color = Strain, linetype = Type))+
  scale_x_continuous(breaks = seq(-20 , 20, by = 1))+
  ylab('Scaled density') +
  xlab('Normalised log2 AUC, OD*h') +
  labs(color = 'Strain',
       linetype = '5-FU') + 
  facet_wrap(~Metabolite, ncol = 4, scales = 'free_x')


dev.copy2pdf(device = cairo_pdf,
             file = paste(odir,"/Density_logConc_by_metabolite.pdf",sep = ''),
             width = 15, height = 160, useDingbats = FALSE)




# Summarise variance
# Find experimental groups which are most variable

sum.c <- data.b %>%
  group_by(MetaboliteU, Strain, Type) %>%
  summarise(Mean = mean(logAUC),
            SD = sd(logAUC)) %>%
  mutate(Index = paste(Strain, ifelse(Type == 'Control','C','T'), MetaboliteU, sep = '_'),
         VarPrc = ifelse(is.na(SD), Inf, (2^(SD) - 1) * 100 )) %>%
  # Reorder, then preserve order with factor levels
  arrange(VarPrc) %>%
  data.frame %>%
  mutate(Index = factor(Index, levels = Index, labels = Index))


sum.c %>%
  filter(VarPrc > 300)


ggplot(sum.c,aes(x = Index, y = VarPrc))+
  geom_point() +
  scale_y_continuous(breaks = seq(0, 1400, by = 50), limits = c(0, 1400))+
  ylab('In-group variation, %')+
  coord_flip()

dev.copy2pdf(device = cairo_pdf,
             file = paste(odir,"/Ingroup_variation_Percentage.pdf", sep = ''),
             width = 10,height = 200, useDingbats = FALSE)















###########################
### pyrE vs BW analysis ###
###########################


# helping function to plot linear equations 
stat_smooth_func <- function(mapping = NULL, data = NULL,
                        geom = "smooth", position = "identity",
                        ...,
                        method = "auto",
                        formula = y ~ x,
                        se = TRUE,
                        n = 80,
                        span = 0.75,
                        fullrange = FALSE,
                        level = 0.95,
                        method.args = list(),
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE,
                        xpos = NULL,
                        ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}

StatSmoothFunc <- ggproto("StatSmooth", Stat,
                      
                      setup_params = function(data, params) {
                        # Figure out what type of smoothing to do: loess for small datasets,
                        # gam with a cubic regression basis for large data
                        # This is based on the size of the _largest_ group.
                        if (identical(params$method, "auto")) {
                          max_group <- max(table(data$group))
                          
                          if (max_group < 1000) {
                            params$method <- "loess"
                          } else {
                            params$method <- "gam"
                            params$formula <- y ~ s(x, bs = "cs")
                          }
                        }
                        if (identical(params$method, "gam")) {
                          params$method <- mgcv::gam
                        }
                        
                        params
                      },
                      
                      compute_group = function(data, scales, method = "auto", formula = y~x,
                                               se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                               xseq = NULL, level = 0.95, method.args = list(),
                                               na.rm = FALSE, xpos=NULL, ypos=NULL) {
                        if (length(unique(data$x)) < 2) {
                          # Not enough data to perform fit
                          return(data.frame())
                        }
                        
                        if (is.null(data$weight)) data$weight <- 1
                        
                        if (is.null(xseq)) {
                          if (is.integer(data$x)) {
                            if (fullrange) {
                              xseq <- scales$x$dimension()
                            } else {
                              xseq <- sort(unique(data$x))
                            }
                          } else {
                            if (fullrange) {
                              range <- scales$x$dimension()
                            } else {
                              range <- range(data$x, na.rm = TRUE)
                            }
                            xseq <- seq(range[1], range[2], length.out = n)
                          }
                        }
                        # Special case span because it's the most commonly used model argument
                        if (identical(method, "loess")) {
                          method.args$span <- span
                        }
                        
                        if (is.character(method)) method <- match.fun(method)
                        
                        base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                        model <- do.call(method, c(base.args, method.args))
                        
                        m = model
                        eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                         list(a = format(coef(m)[[1]], digits = 3), 
                                              b = format(coef(m)[[2]], digits = 3), 
                                              r2 = format(summary(m)$r.squared, digits = 3)))
                        func_string = as.character(as.expression(eq))
                        
                        if(is.null(xpos)) xpos = min(data$x)*0.9
                        if(is.null(ypos)) ypos = max(data$y)*0.9
                        data.frame(x=xpos, y=ypos, label=func_string)
                        
                      },
                      
                      required_aes = c("x", "y")
)



# creates the intercepts of lines that mark the nutrients
intercepts2 = data.b %>%
    filter(Metabolite == "Negative Control") %>%
    group_by(Strain, Type, Plate) %>%
    summarise(Mean = mean(logAUC_raw),
        Sd = sd(logAUC_raw)) %>%
    arrange(Plate) %>%
    ungroup %>%
    filter(Strain %in% c('BW', 'pyrE'), Type == 'Control') %>%
    group_by(Plate) %>%
    summarise(intercept = dumm(Mean)) %>%
    data.frame

# these are the best nutrients, remember to change plate and value!!!
best_PM1 = growth %>% filter(pyrE_C_mean >= (BW_C_mean + 0.8170250),   Plate == "PM1") 
best_PM2 = growth %>% filter(pyrE_C_mean >= (BW_C_mean + 0.7615000),   Plate == "PM2A")
best_PM3 = growth %>% filter(pyrE_C_mean >= (BW_C_mean + 0.7558500),   Plate == "PM3B")
best_PM4 = growth %>% filter(pyrE_C_mean >= (BW_C_mean + 0.6805375),   Plate == "PM4A")

# compute linear models of each plate
lm_PM1 = growth %>% filter(Plate == 'PM1') %>% lm(pyrE_C_mean ~ BW_C_mean, data = .) 
lm_PM2 = growth %>% filter(Plate == 'PM2A') %>% lm(pyrE_C_mean ~ BW_C_mean, data = .) 
lm_PM3 = growth %>% filter(Plate == 'PM3B') %>% lm(pyrE_C_mean ~ BW_C_mean, data = .) 
lm_PM4 = growth %>% filter(Plate == 'PM4A') %>% lm(pyrE_C_mean ~ BW_C_mean, data = .) 



# plot, dont forget to change the text variable in geom_repel

plate = 'PM4A'
growth %>% filter(Plate == plate) %>% 
ggplot(aes(x = BW_C_mean, y = pyrE_C_mean)) + 
    geom_errorbar(aes(ymin = (pyrE_C_mean - pyrE_C_sd), ymax = (pyrE_C_mean + pyrE_C_sd)), 
      colour = "grey50", alpha = 0.5) +
    geom_errorbarh(aes(xmin = BW_C_mean - BW_C_sd, xmax = BW_C_mean + BW_C_sd), 
      colour = "grey50", alpha = 0.5) +
    geom_point(alpha = 0.8, size = 1.5) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "grey50")) +
    labs(
    title = 'test',
    x = paste('BW', "growth", sep = ' '),
    y = paste('pyrE', "growth", sep = ' ')) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", alpha = 0.8) + 
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50", alpha = 0.8) +
    geom_abline(intercept = intercepts2[4,2], slope = 1, linetype = "dashed", colour = "blue", alpha = 0.7) +
    stat_smooth(method = "lm", col = "red", size = 0.5, alpha = 0.1) +
    stat_smooth_func(geom = "text", method = "lm", hjust = 0, parse = TRUE) +
    geom_text_repel(data = best_PM4, aes(label = MetaboliteU))

  

ggsave(file = paste0(odir,"/Scatter_BW_pyrE_", plate, ".pdf"),
       width = 160, height = 120, units = 'mm', scale = 2, device = 'pdf')



## subplots with labels  in other quadrants


# plot, dont forget to change the text variable in geom_repel

plate = 'PM1'
growth %>% filter(Plate == plate, Metabolite != 'Negative Control') %>% 
ggplot(aes(x = BW_C_mean, y = pyrE_C_mean)) + 
    geom_errorbar(aes(ymin = (pyrE_C_mean - pyrE_C_sd), ymax = (pyrE_C_mean + pyrE_C_sd)), 
      colour = "grey50", alpha = 0.5) +
    geom_errorbarh(aes(xmin = BW_C_mean - BW_C_sd, xmax = BW_C_mean + BW_C_sd), 
      colour = "grey50", alpha = 0.5) +
    geom_point(alpha = 0.8, size = 1.5) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "grey50")) +
    labs(
    title = 'test',
    x = paste('BW', "growth", sep = ' '),
    y = paste('pyrE', "growth", sep = ' ')) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", alpha = 0.8) + 
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50", alpha = 0.8) +
    geom_abline(intercept = intercepts2[4,2], slope = 1, linetype = "dashed", colour = "blue", alpha = 0.7) +
    stat_smooth(method = "lm", col = "red", size = 0.5, alpha = 0.1) +
    stat_smooth_func_with_pval(geom = "text", method = "lm", hjust = 0, parse = TRUE) +
    geom_text_repel(aes(label = ifelse(BW_C_mean >= 0 & pyrE_C_mean <= 0, MetaboliteU, '')))

ggsave(file = paste0(odir,"/Scatter_BW_pyrE_", plate, "_subplot.pdf"),
       width = 160, height = 120, units = 'mm', scale = 2, device = 'pdf')












#################################
### Hightly experimental zone ###
#################################

# as we lack the form to put several plates within the same scatter plot, 
# I've tried to normalise using a quantile normalisation. Didn't work, but 
# here is the code
 
#### quantile normalisation part

# select only the numeric data from all samples

sub1 = data.b %>%
    unite(STR, Strain, Type, ReplicateB) %>% # Generete unique IDs by concatenating Strain, Type and Replicate values
    select(STR, Plate, Well, Row, Col, Index:Class,logAUC_raw, -Replicate) %>% 
    spread(STR, logAUC_raw) %>%
    select(BW_Control_B1:TM_Treatment_B4) %>%
    data.frame


# this plots the densities of the different samples
colramp = colorRampPalette(c(3,"white",2))(31)
plot(density(sub1[,1]),col = colramp[1], lwd = 3)
for(i in 2:31){lines(density(sub1[,i]),lwd=3,col=colramp[i])}

# quantile normalisation
qnor = preprocessCore::normalize.quantiles(as.matrix(sub1))
colnames(qnor) = colnames(sub1)

# this plots the densities of normalised data sets
colramp = colorRampPalette(c(3,"white",2))(31)
plot(density(qnor[,1]),col = colramp[1], lwd = 3)
for(i in 2:31){lines(density(qnor[,i]),lwd=3,col=colramp[i])}


temp1 = data.b %>%
    unite(STR, Strain, Type, ReplicateB) %>% # Generete unique IDs by concatenating Strain, Type and Replicate values
    select(STR, Plate, Well, Row, Col, Index:Class,logAUC, -Replicate) %>% 
    spread(STR, logAUC)

temp1 = temp1[,-12:-42]
growth2 = (cbind(temp1, qnor))

# normalise against negative control
for (i in 12:42){
  growth2[,i] = growth2[,i] - growth2[1,i]
}
growth2 = tbl_df(growth2)


# check
growth2[1,12:42]



growth1 = growth2 %>%
    group_by(Plate, Well, Index, Metabolite, MetaboliteU) %>%
    mutate(BW_C_mean = mean(c(BW_Control_B1, BW_Control_B2, BW_Control_B3, BW_Control_B4)),
           BW_C_sd = sd(c(BW_Control_B1, BW_Control_B2, BW_Control_B3, BW_Control_B4)),
           BW_T_mean = mean(c(BW_Treatment_B1, BW_Treatment_B2, BW_Treatment_B3, BW_Treatment_B4)),
           BW_T_sd = sd(c(BW_Treatment_B1, BW_Treatment_B2, BW_Treatment_B3, BW_Treatment_B4)),
           crp_C_mean = mean(c(crp_Control_B1, crp_Control_B2, crp_Control_B3, crp_Control_B4)),
           crp_C_sd = sd(c(crp_Control_B1, crp_Control_B2, crp_Control_B3, crp_Control_B4)),
           crp_T_mean = mean(c(crp_Treatment_B1, crp_Treatment_B2, crp_Treatment_B3, crp_Treatment_B4)),
           crp_T_sd = sd(c(crp_Treatment_B1, crp_Treatment_B2, crp_Treatment_B3, crp_Treatment_B4)),
           pyrE_C_mean = mean(c(pyrE_Control_B1, pyrE_Control_B2, pyrE_Control_B3, pyrE_Control_B4)),
           pyrE_C_sd = sd(c(pyrE_Control_B1, pyrE_Control_B2, pyrE_Control_B3, pyrE_Control_B4)),
           pyrE_T_mean = mean(c(pyrE_Treatment_B1, pyrE_Treatment_B2, pyrE_Treatment_B4)),
           pyrE_T_sd = sd(c(pyrE_Treatment_B1, pyrE_Treatment_B2, pyrE_Treatment_B4)),
           TM_C_mean = mean(c(TM_Control_B1, TM_Control_B2, TM_Control_B3, TM_Control_B4)),
           TM_C_sd = sd(c(TM_Control_B1, TM_Control_B2, TM_Control_B3, TM_Control_B4)),
           TM_T_mean = mean(c(TM_Treatment_B1, TM_Treatment_B2, TM_Treatment_B3, TM_Treatment_B4)),
           TM_T_sd = sd(c(TM_Treatment_B1, TM_Treatment_B2, TM_Treatment_B3, TM_Treatment_B4))) %>% 
    ungroup 



# infernal plot with too many variables and random stuff
p1 = growth1 %>%
    filter(Plate == "PM1", Metabolite != "Negative Control") %>%
    ggplot(aes(x = BW_C_mean, y = BW_T_mean)) +
    geom_errorbar(aes(ymin = BW_T_mean - BW_T_sd, ymax = BW_T_mean + BW_T_sd), 
      colour = "grey50", alpha = 0.5) +
    geom_errorbarh(aes(xmin = BW_C_mean - BW_C_sd, xmax = BW_C_mean + BW_C_sd), 
      colour = "grey50", alpha = 0.5) +
    geom_point(alpha = 0.8, size = 1.5) +
    # geom_point(aes(colour = Diff_mean) , alpha = 0.8) + 
    # xlim(-1.5, 1) + ylim(-1.5, 2.5) +  # delimit plot 
    # xlim(0, 0.8) + ylim(0.8, 2.2) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "grey50")) +
    labs(
    title = "Bacterial growth, PM1 plate",
    x = "BW growth",
    y = expression(paste(Delta, italic("gcvP"), " growth"), sep = "")) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", alpha = 0.8) + 
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50", alpha = 0.8)








