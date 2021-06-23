### RNA seq analysis of human cell lines from Tanara

# In this script we will analyse the RNA seq with the following cell lines:
# 	- HCT116
# 	- LoVo
# 	- SW948 
# 	- DLD-1

# useful links:
# 	- https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#10_getting_or_building_ensdb_databasespackages
# 	- https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# 	- https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts

# maybe use ulimit -s 16384 before start R console

### libraries ####
library(tximport)
library(tidyverse)
library(DESeq2)
# notice that DESeq2 library masks 'rename' function from dplyr 
# library(ensembldb) # use only if you are going to deal with db
here::set_here()
library(here)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(openxlsx)
library(viridis)
library(glue)


theme_set(theme_classic())

# the first step we need to do is to read the sample file, and parse it with 
# the data produced by salmon to get tables ready to be analysed by DESeq2



# Get sample info ---------------------------------------------------------



samples = read.delim("sampleInfo.txt") 

dir = getwd()
rownames(samples) = samples$Name

quants_dir = 'quants_103'

# load kegg tables from wormenrichr
# kegg = read.delim("KEGG_2019.txt", header = FALSE) 
# kegg = kegg[,-2]
# 
# 
# rownames(kegg) = kegg[,1] ; kegg = kegg[,-1]

# prepare a list with file names
files = file.path(dir,quants_dir, samples$Name, "quant.sf")
names(files) = samples$Name
all(file.exists(files)) # check that files exist






# Get gene labels from databases ------------------------------------------




#### manual id handling  #####

# 
# gtf_file = 'Homo_sapiens.GRCh38.104.gtf.gz'
# gtf_file = 'Homo_sapiens.GRCh38.103.gtf.gz'
# txdb = GenomicFeatures::makeTxDbFromGFF(gtf_file, organism='Homo sapiens')
# 
# 
# k = keys(txdb, keytype="TXNAME")
# tx_map = ensembldb::select(txdb, keys = k, columns="GENEID", keytype = "TXNAME")
# 
# head(tx_map)
# 
# tx2gene = tx_map

# write.csv(tx2gene,file="tx2gene.csv",row.names = FALSE,quote=FALSE)
#
# quants = read_tsv(files[1])
# 
# quants <- separate(quants, Name, c("TXNAME","Number"),remove = FALSE)
# head(quants)
# 
# quants <- left_join(quants, tx_map, by="TXNAME")
# head(quants)
# 
# tx2gene <- dplyr:::select(quants, Name, GENEID)
# head(tx2gene)
# tx2gene <- filter(tx2gene, !is.na(GENEID))




#### using ensembl genomes ####

# let's make our database from ensembldb
ah = AnnotationHub::AnnotationHub(proxy='127.0.0.1:10801')
ahDb = AnnotationHub::query(ah, pattern = c("Homo sapiens", "EnsDb", 103))
ahEdb = ahDb[[1]]
# generate the database
tx2gene.complete = transcripts(ahEdb, return.type = "DataFrame")

# fetch descriptions of genes
info = genes(ahEdb) %>%
  as_tibble() %>%
  dplyr::select(width, gene_id, gene_name, gene_biotype, description, entrezid) %>%
  unnest(cols = c(entrezid))

# join transcription info with gene ids and entrezids
info.join = tx2gene.complete %>%
  as_tibble() %>%
  dplyr::select(tx_id, tx_biotype, gene_id, tx_name) %>%
  left_join(info)

write.csv(info, here('summary','gene_ids_mapping.csv'))

# subset to have tx_id in first column, and gene_id in second
# tx2gene = data.frame(tx2gene.complete[,c(1,7)])


# subset to have tx_id in first column, and gene_id in second
tx2gene = data.frame(tx2gene.complete[,c(9,7)])

colnames(tx2gene) = c('tx_id','gene_id')

head(tx2gene)



# 
# library( "biomaRt" )
# 
# ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
# 


### TX import + DESeq ####

# import quantification data 
txi = tximport::tximport(files, type = "salmon", tx2gene = tx2gene)


### starting analysis with DESeq2
# create DESeq data type to be analysed


samples.batch = samples %>% 
  mutate(Batch = as.factor(Batch),
         Condition = as.factor(Condition),
         Sample = as.factor(Sample))

ddsTxi = DESeqDataSetFromTximport(txi, colData = samples.batch, 
                                  design = ~Batch + Sample)



#### filter by rowSums ####
# prefilter, but that might not be necessary
keep = rowSums(counts(ddsTxi)) >= 100

ddsTxi = ddsTxi[keep,]


#### 0s filter ####
# filter by presence of 0s in the samples
# STRATEGY: get sample columns per separate, calculate how many
# columns have 0s, select the genes that have 6 (or more) columns with
# 0s and change the rest of the values to 0


temp = counts(ddsTxi)

cells = unique(samples$Cell_line)

threshold = 6

# loop to cycle for every sample subset and fix weird values
for (cell in cells) {
samp_names = samples %>% 
  filter(Cell_line == cell) %>% 
  pull(Name)
flawed = rownames(temp[,samp_names][rowSums(temp[,samp_names] == 0) >= threshold,])
temp[flawed,samp_names] = 0
}

 
temp %>% view

ddsTxi.filt = DESeqDataSetFromMatrix(temp ,
                       colData = samples.batch,
                       design = ~Batch + Sample)


# 
# temp[flawed,samp_names]





ddsTxi.filt$Sample = relevel(ddsTxi$Sample, ref = "HCT116_C")


# run the Differential Expression Analysis
dds = DESeq(ddsTxi.filt)






### tidy results ####
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
gene_counts = counts(dds, normalized = TRUE)
gene_list = rownames(gene_counts)
gene_counts = gene_counts %>% cbind(gene_list,.) %>% as_tibble()

gene_counts = gene_counts %>% 
  gather(Name, counts, TVP_1:TVP_32) %>% 
  dplyr::select(gene_id = gene_list, Name, counts) %>%
  mutate(Name = as.factor(Name)) %>%
  left_join(as_tibble(samples), by = 'Name') %>%
  mutate(counts = as.double(counts),
         gene_id = as.factor(gene_id),
         Replicate = as.factor(Replicate)) %>% 
  left_join(info) 


gene = 'HSPA1B'
gene_counts%>% 
  dplyr::filter(gene_name == gene) %>% 
  dplyr::filter(Cell_line == 'HCT116') %>%
  ggplot(aes(x = Sample, y = counts, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2.5, position = position_jitterdodge(), aes(color = Replicate))+
  facet_wrap(~gene_id*Cell_line, scales = 'free')

ggsave(here('summary', glue('boxplot_{gene}.pdf')), height = 7, width = 9)

# transofrm data

vsd = vst(dds, blind = FALSE)
rld = rlog(dds, blind = FALSE)


# plot differences between different transformation data methods
df = bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("mean", "sd")  

ggplot(df, aes(x = mean, y = sd)) + 
  geom_hex(bins = 100) +
  # coord_fixed() + 
  facet_grid( . ~ transformation) +
  theme_light()

dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'transformation_comparison.pdf'),
             height = 4.5, width = 12, useDingbats = FALSE)


###
# Sample distances

sampleDists = dist(t(assay(rld)))
sampleDists = as.matrix(sampleDists)


# USE THIS TO CHANGE THE ROW/COLUMN NAMES
# names = colnames(sampleDists) %>%
#   str_split('_', simplify = T) %>%
#   data.frame %>% tbl_df() %>%
#   unite(sample, X1, X2, sep = " - ") %>%
#   dplyr::select(sample) %>%
#   t %>% as.vector

names = gene_counts %>% 
  unite(ID, Sample, Replicate) %>% 
  distinct(ID) %>% pull(ID)


colnames(sampleDists) = names; rownames(sampleDists) = names

col_fun = colorRamp2(c(0, 100), c("white", "blue"))
# col_fun(seq(0, 100))
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

Heatmap(sampleDists, name = 'Euclidean \ndistances', 
        col = colors)



dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'Euclidean_distances_samples.pdf'),
             height = 8, width = 9, useDingbats = FALSE)





# PCA ---------------------------------------------------------------------




pcaData = plotPCA(rld, intgroup = c("Cell_line", "Condition"), returnData = TRUE)
pcaData


# get ellipses based on the correlation
getellipse = function(x, y, sc = 1) {
  as.data.frame(ellipse::ellipse(cor(x, y),
                                 scale = c(sd(x) * sc, sd(y) * sc),
                                 centre = c(mean(x), mean(y))))
}



# get info for the ellipses
ell = pcaData %>% group_by(Condition, Cell_line) %>% do(getellipse(.$PC1, .$PC2, 1)) %>% data.frame

# % of variable explained by PC1 and PC2
percentVar = round(100 * attr(pcaData, "percentVar"))

# plot!
ggplot(pcaData, aes(x = PC1, y = PC2, color = Cell_line, group = interaction(Condition, Cell_line))) + 
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell, aes(x = x, y = y, group = interaction(Condition, Cell_line), linetype = Condition), size = 1) +
  geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Condition, Cell_line), 
                               linetype = Condition, fill = Cell_line), size = 1, alpha = 0.3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, colour = 'black'),
        axis.text.y = element_text(size = 13, colour = 'black')) 


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'PCA_main_rld.pdf'),
             height = 8, width = 9, useDingbats = FALSE)


# 
# 
# 
# KEGG databases
# get databases for genes and pathways from KEGG
kegg.links.entrez = limma::getGeneKEGGLinks('hsa', convert = TRUE)
kegg.links.ids = limma::getGeneKEGGLinks('hsa')
path.ids = limma::getKEGGPathwayNames('hsa', remove.qualifier = TRUE)
kegg.links = cbind(kegg.links.entrez, kegg.links.ids[,1])
colnames(kegg.links) = c('entrezid', 'PathwayID', 'KEGG_genes')

kegg.links = kegg.links %>%
  as_tibble %>%
  mutate(entrezid = as.integer(entrezid)) %>%
  left_join(path.ids) %>%
  mutate(PathwayID = str_replace(PathwayID, 'path:hsa', ''))




# General stats -----------------------------------------------------------


# get results and tidy it
res = results(dds) 


# results with different shape of contrasts, tidy
res.hct = results(dds,   contrast = c("Sample", "HCT116_M" , "HCT116_C"))  
res.hct = lfcShrink(dds, contrast = c("Sample", "HCT116_M" , "HCT116_C"), res = res.hct, type = 'ashr')

res.dld = results(dds,  contrast = c("Sample",  "DLD1_M", "DLD1_C")) 
res.dld = lfcShrink(dds, contrast = c("Sample",  "DLD1_M", "DLD1_C"), res = res.dld, type = 'ashr')

res.lovo = results(dds,  contrast = c("Sample",  "LoVo_M", "LoVo_C"))   
res.lovo = lfcShrink(dds, contrast = c("Sample", "LoVo_M", "LoVo_C"), res = res.lovo, type = 'ashr')

res.sw = results(dds, contrast = c("Sample",   "SW948_M", "SW948_C")) 
res.sw = lfcShrink(dds, contrast = c("Sample", "SW948_M", "SW948_C"), res = res.sw, type = 'ashr')

#### tidying results ####
res.hct.tidy = as_tibble(res.hct, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'HCT116',
  Contrast_description = 'Comparison of HCT116 Micit vs HCT116 Control') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.dld.tidy = as_tibble(res.dld, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'DLD-1',
  Contrast_description = 'Comparison of DLD-1 Micit vs DLD-1 Control') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.lovo.tidy = as_tibble(res.lovo, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'LoVo',
  Contrast_description = 'Comparison of LoVo Micit vs LoVo Control') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.sw.tidy = as_tibble(res.sw, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'SW948',
  Contrast_description = 'Comparison of SW948 Micit vs SW948 Control') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

results.complete = res.hct.tidy %>% rbind(res.dld.tidy, res.lovo.tidy, res.sw.tidy)

results.complete.kegg = results.complete %>% left_join(kegg.links)

# write results in excel files
list_of_datasets = list('HCT116' = res.hct.tidy, 
                        'DLD-1' = res.dld.tidy, 
                        'LoVo' = res.lovo.tidy,
                        'SW948' = res.sw.tidy)

write.xlsx(list_of_datasets, here('summary', 'complete_stats.xlsx'), 
           colNames = T, rowNames = F) 



# plot genes


gene = 'MYO15B'
gene_counts%>% 
  dplyr::filter(gene_name == gene) %>% 
  filter(Cell_line == 'HCT116') %>%
  ggplot(aes(x = Sample, y = counts, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2.5, position = position_jitterdodge(), aes(color = Replicate))+
  facet_wrap(~gene_id, scales = 'free_y')



plotCounts(res.hct, gene=1000, intgroup="Sample")



# MA plots ----------------------------------------------------------------



### MA plots for every comparison
plotMA(res.hct,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_HCT.pdf'),
             height = 8, width = 11, useDingbats = FALSE)

plotMA(res.dld,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_DLD.pdf'),
             height = 8, width = 11, useDingbats = FALSE)


plotMA(res.lovo,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_LoVo.pdf'),
             height = 8, width = 11, useDingbats = FALSE)


plotMA(res.sw,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_SW948.pdf'),
             height = 8, width = 11, useDingbats = FALSE)







# PCA (manual) ------------------------------------------------------------




###
# PCA data manual way

# Worm = Condition
# Bacteria =  Cell.line

# # pca_data = t(counts(rld, normalized = TRUE))
pca_data = t(assay(rld))

# HCT 116 samples
pca_data = pca_data[c(1,2,9,10,17,18,25,26),]
# DLD_1
pca_data = pca_data[c(3,4,11,12,19,20,27,28),]
# Lovo
pca_data = pca_data[c(3,4,11,12,19,20,27,28)+2,]
# SW948
pca_data = pca_data[c(3,4,11,12,19,20,27,28)+4,]

# # lets compute the PCA
res.pca = PCA(pca_data, scale.unit = FALSE, ncp = 5, graph = F)

# # metadata 
meta_var = samples

# HCT samples
meta_var = samples[c(1,2,9,10,17,18,25,26),]

# DLD_1
meta_var = samples[c(3,4,11,12,19,20,27,28),]

# LoVo
meta_var = samples[c(3,4,11,12,19,20,27,28)+2,]

# SW948
meta_var = samples[c(3,4,11,12,19,20,27,28)+4,]

# # extract info about the individuals
ind = get_pca_ind(res.pca)
ind_df = data.frame(ind$coord[,1], ind$coord[,2], ind$coord[,3], meta_var$Condition,
  meta_var$Cell_line)

colnames(ind_df) = c('Dim1', 'Dim2', 'Dim3', 'Condition', 'Cell_line')


# # make a data frame from ellipses
ell = ind_df %>% group_by(Condition, Cell_line) %>% do(getellipse(.$Dim1, .$Dim2, 1)) %>% data.frame

# % of variable explained by PC1 and PC2
percentVar = round(100 * attr(pcaData, "percentVar"))

# plot!
line = 'HCT116'
ggplot(ind_df, aes(x = Dim1, y = Dim2, color = Condition, group = interaction(Condition, Cell_line))) +
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell, aes(x = x, y = y, group = interaction(Condition, Cell_line), linetype = Condition), size = 1) +
  geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Condition, Cell_line), 
                               linetype = Condition, fill = Condition), size = 1, alpha = 0.3) +
  # xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  # ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, colour = 'black'),
        axis.text.y = element_text(size = 13, colour = 'black'))

ggsave(here('summary',glue('PCA_{line}.pdf')), height = 7, width = 8)



# GSEA --------------------------------------------------------------------



res.hct.tidy %>% 
  filter(str_detect(gene_name,'CPT'))

gene = 'DDIT4'
gene_counts%>% 
  dplyr::filter(gene_name == gene) %>% 
  filter(Cell_line == 'HCT116') %>%
  ggplot(aes(x = Sample, y = counts, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2.5, position = position_jitterdodge(), aes(color = Replicate))+
  facet_wrap(~gene_id, scales = 'free_y')




# DE genes (up/down) to analyse via StringDB

# helper function that extracts DE genes with upper and lower thresholds
DEgenes = function(dataset = results.complete, 
                   cell = 'HCT116', 
                   min_thrs = 0.5,
                   max_thrs = 10) {
  
  up = dataset %>%  
    filter(Contrast == cell) %>% 
    filter(log2FoldChange < max_thrs & log2FoldChange > min_thrs) %>% 
    filter(padj <= 0.05) %>%
    pull(gene_id) %>% unique
  
  
  down = dataset %>% 
    filter(Contrast == cell) %>% 
    filter(log2FoldChange > -max_thrs & log2FoldChange < -min_thrs) %>% 
    filter(padj <= 0.05) %>%
    pull(gene_id) %>% unique
  
  
  up = c('genes', up)
  down = c('genes', down)
  
  up_name = glue('{cell}_UP')
  down_name = glue('{cell}_DOWN')
  
  list_of_datasets = list(
    up_name = up,
    down_name = down
  )
  
  names(list_of_datasets) = c(up_name, down_name)

  write.xlsx(list_of_datasets, here('summary', glue('{cell}_genes_updown.xlsx')))
}

# HCT
DEgenes(results.complete, cell = 'HCT116', min_thrs = 0.5, max_thrs = 10)
# DLD-1
DEgenes(results.complete, cell = 'DLD-1', min_thrs = 0, max_thrs = 10)
# LoVo
DEgenes(results.complete, cell = 'LoVo', min_thrs = 0, max_thrs = 10)
# SW948
DEgenes(results.complete, cell = 'SW948', min_thrs = 0, max_thrs = 10)








