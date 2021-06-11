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

# analysis directory:  "/Users/dmarti14/Documents/MRC_Postdoc/Projects/Communities/RNA_seq"
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


# the first step we need to do is to read the sample file, and parse it with 
# the data produced by salmon to get tables ready to be analysed by DESeq2

samples = read.delim("sampleInfo.txt") 

dir = getwd()
rownames(samples) = samples$Name

# load kegg tables from wormenrichr
# kegg = read.delim("KEGG_2019.txt", header = FALSE) 
# kegg = kegg[,-2]
# 
# 
# rownames(kegg) = kegg[,1] ; kegg = kegg[,-1]

# prepare a list with file names
files = file.path(dir,"quants", samples$Name, "quant.sf")
names(files) = samples$Name
all(file.exists(files)) # check that files exist




gtf_file = 'Homo_sapiens.GRCh38.104.gtf.gz'

txdb = GenomicFeatures::makeTxDbFromGFF(gtf_file, organism='Homo sapiens')


k = keys(txdb, keytype="TXNAME")
tx_map = ensembldb::select(txdb, keys = k, columns="GENEID", keytype = "TXNAME")

head(tx_map)

tx2gene = tx_map

write.csv(tx2gene,file="tx2gene.csv",row.names = FALSE,quote=FALSE)

quants = read_tsv(files[1])

quants <- separate(quants, Name, c("TXNAME","Number"),remove = FALSE)
head(quants)

quants <- left_join(quants, tx_map, by="TXNAME")
head(quants)

tx2gene <- dplyr:::select(quants, Name, GENEID)
head(tx2gene)
tx2gene <- filter(tx2gene, !is.na(GENEID))



# BiocManager::install("EnsDb.Hsapiens.v86")




# # let's make our database from ensembldb 
# ah = AnnotationHub::AnnotationHub(proxy='127.0.0.1:10801')
# ahDb = AnnotationHub::query(ah, pattern = c("Homo sapiens", "EnsDb"))
# ahEdb = ahDb[[1]]
# # generate the database 
# tx2gene.complete = transcripts(ahEdb, return.type = "DataFrame")
# 
# # fetch descriptions of genes
# info = genes(ahEdb) %>% 
#   as_tibble() %>%
#   dplyr::select(width, gene_id, gene_name, gene_biotype, description, entrezid) %>%
#   unnest(cols = c(entrezid))
# 
# # join transcription info with gene ids and entrezids
# info.join = tx2gene.complete %>% 
#   tbl_df() %>%
#   dplyr::select(tx_id, tx_biotype, gene_id, tx_name) %>%
#   left_join(info)
# 
# write.csv(info, here('summary','gene_ids_mapping.csv'))
# 
# # subset to have tx_id in first column, and gene_id in second
# tx2gene = data.frame(tx2gene.complete[,c(1,7)])



# import quantification data 
txi = tximport::tximport(files, type = "salmon", tx2gene = tx2gene,
                         ignoreTxVersion = TRUE)




### starting analysis with DESeq2
# create DESeq data type to be analysed
ddsTxi = DESeqDataSetFromTximport(txi, colData = samples, design = ~ Sample)

# # prefilter, but that might not be necessary
# keep = rowSums(counts(ddsTxi)) >= 10
# ddsTxi = ddsTxi[keep,]

ddsTxi$Sample = relevel(ddsTxi$Sample, ref = "HCT116M")

# run the Differential Expression Analysis
# design(ddsTxi) <- formula(~ Bacteria + Worm)
dds = DESeq(ddsTxi)





### tidy results
# dds.tidy = tidy(ddsTxi, colData = samples$Sample)
gene_counts = counts(dds, normalized = TRUE)
gene_list = rownames(gene_counts)
gene_counts = gene_counts %>% cbind(gene_list,.) %>% as_tibble()

gene_counts = gene_counts %>% 
  gather(Name, counts, TVP_1_quant:TVP_32__quant) %>% 
  dplyr::select(gene_id = gene_list, Name, counts) %>%
  mutate(Name = as.factor(Name)) %>%
  left_join(as_tibble(samples), by = 'Name') %>%
  mutate(counts = as.double(counts),
         gene_id = as.factor(gene_id),
         Replicate = as.factor(Replicate)) %>% 
  left_join(info) 


gene = 'TYMS'
gene_counts%>% 
  dplyr::filter(gene_name == gene) %>% 
  filter(Cell.line == 'HCT116') %>%
  ggplot(aes(x = Sample, y = counts, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 5, position = position_jitterdodge(), aes(color = Replicate))


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
names = colnames(sampleDists) %>%
  str_split('_', simplify = T) %>%
  data.frame %>% tbl_df() %>%
  unite(sample, X1, X2, sep = " - ") %>%
  dplyr::select(sample) %>%
  t %>% as.vector

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










pcaData = plotPCA(rld, intgroup = c("Cell.line", "Condition"), returnData = TRUE)
pcaData


# get ellipses based on the correlation
getellipse = function(x, y, sc = 1) {
  as.data.frame(ellipse::ellipse(cor(x, y),
                                 scale = c(sd(x) * sc, sd(y) * sc),
                                 centre = c(mean(x), mean(y))))
}



# get info for the ellipses
ell = pcaData %>% group_by(Condition, Cell.line) %>% do(getellipse(.$PC1, .$PC2, 1)) %>% data.frame

# % of variable explained by PC1 and PC2
percentVar = round(100 * attr(pcaData, "percentVar"))

# plot!
ggplot(pcaData, aes(x = PC1, y = PC2, color = Cell.line, group = interaction(Condition, Cell.line))) + 
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell, aes(x = x, y = y, group = interaction(Condition, Cell.line), linetype = Condition), size = 1) +
  geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Condition, Cell.line), 
                               linetype = Condition, fill = Cell.line), size = 1, alpha = 0.3) +
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


# get results and tidy it
res = results(dds) 


# results with different shape of contrasts, tidy
res.hct = results(dds,   contrast = c("Sample", "HCT116M" , "HCT116C"))  
res.hct = lfcShrink(dds, contrast = c("Sample", "HCT116M" , "HCT116M"), res = res.hct, type = 'ashr')

res.dld = results(dds,  contrast = c("Sample",  "DLD-1M", "DLD-1C")) 
res.dld = lfcShrink(dds, contrast = c("Sample",  "DLD-1M", "DLD-1C"), res = res.dld, type = 'ashr')

res.OGRF = results(dds,  contrast = c("Sample",  "skpo_OG1RF", "N2_OG1RF"))   
res.OGRF = lfcShrink(dds, contrast = c("Sample", "skpo_OG1RF", "N2_OG1RF"), res = res.dld, type = 'ashr')

res.OP50 = results(dds, contrast = c("Sample",   "skpo_OP50", "N2_OP50")) 
res.OP50 = lfcShrink(dds, contrast = c("Sample", "skpo_OP50", "N2_OP50"), res = res.OP50, type = 'ashr')

# tidying the results
res.hct.tidy = as_tibble(res.hct, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'HCT116',
  Contrast_description = 'Comparison of HCT116 Micit vs HCT116 Control') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.skpo.tidy = as_tibble(res.skpo, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'skpo',
  Contrast_description = 'Comparison of skpo_OG1RF vs skpo_OP50') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.OGRF.tidy = as_tibble(res.OGRF, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'OG1RF',
  Contrast_description = 'Comparison of skpo_OG1RF vs N2_OG1RF') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

res.OP50.tidy = as_tibble(res.OP50, rownames = 'gene_id') %>% mutate(
  p_adj_stars = gtools::stars.pval(padj),
  Direction = ifelse(log2FoldChange > 0, 'Up', 'Down'),
  Contrast = 'OP50',
  Contrast_description = 'Comparison of skpo_OP50 vs N2_OP50') %>%
  left_join(info) %>%
  # left_join(kegg.links) %>%
  dplyr::select(Contrast_description, Contrast, gene_id, gene_name, everything())

results.complete = res.N2.tidy %>% rbind(res.skpo.tidy, res.OGRF.tidy, res.OP50.tidy)

# write results in excel files
list_of_datasets = list('N2 OG1RF_vs_OP50' = res.N2.tidy, 
                        'skpo OG1RF_vs_OP50' = res.skpo.tidy, 
                        'OP50 skpo_vs_N2' = res.OP50.tidy,
                        'OG1RF skpo_vs_N2' = res.OGRF.tidy)

write.xlsx(list_of_datasets, here('summary', 'complete_stats.xlsx'), colNames = T, rowNames = F) 



### MA plots for every comparison
plotMA(res.N2,  ylim=c(-3,3),  alpha = 0.05)
# idx <- identify(res.N2$baseMean, res.N2$log2FoldChange)
# rownames(res.N2)[idx]
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_N2.pdf'),
             height = 8, width = 11, useDingbats = FALSE)

plotMA(res.skpo,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_skpo.pdf'),
             height = 8, width = 11, useDingbats = FALSE)


plotMA(res.OGRF,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_OG1RF.pdf'),
             height = 8, width = 11, useDingbats = FALSE)


plotMA(res.OP50,  ylim=c(-3,3),  alpha = 0.05)
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'MAplot_OP50.pdf'),
             height = 8, width = 11, useDingbats = FALSE)





###
# PCA data manual way

# Worm = Condition
# Bacteria =  Cell.line

# # pca_data = t(counts(rld, normalized = TRUE))
pca_data = t(assay(rld))

# HCT 116 samples
pca_data = pca_data[c(1,2,9,16,17,24,25),]
# DLD_1
pca_data = pca_data[c(3,4,10,11,18,19,26,27),]



# # lets compute the PCA
res.pca = PCA(pca_data, scale.unit = FALSE, ncp = 5, graph = F)

# # metadata 
meta_var = samples

# HCT samples
meta_var = samples[c(1,2,9,16,17,24,25),]

# DLD_1
meta_var = samples[c(3,4,10,11,18,19,26,27),]

# # extract info about the individuals
ind = get_pca_ind(res.pca)
ind_df = data.frame(ind$coord[,1], ind$coord[,2], ind$coord[,3], meta_var$Condition,
  meta_var$Cell.line)

colnames(ind_df) = c('Dim1', 'Dim2', 'Dim3', 'Condition', 'Cell.line')


# # make a data frame from ellipses
ell = ind_df %>% group_by(Condition, Cell.line) %>% do(getellipse(.$Dim1, .$Dim2, 1)) %>% data.frame

# % of variable explained by PC1 and PC2
percentVar = round(100 * attr(pcaData, "percentVar"))

# plot!
ggplot(ind_df, aes(x = Dim1, y = Dim2, color = Cell.line, group = interaction(Condition, Cell.line))) +
  geom_point(size = 3, show.legend = NA, alpha = 0.5) + 
  geom_path(data = ell, aes(x = x, y = y, group = interaction(Condition, Cell.line), linetype = Condition), size = 1) +
  geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Condition, Cell.line), 
                               linetype = Condition, fill = Cell.line), size = 1, alpha = 0.3) +
  # xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  # ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, colour = 'black'),
        axis.text.y = element_text(size = 13, colour = 'black'))






