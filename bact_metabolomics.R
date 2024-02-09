# plotting metabolomics for the paper

library(tidyverse)
library(readxl)
library(cowplot)
library(ComplexHeatmap)

metabolomics = read_excel("Metabolomics.xlsx")


colnames(metabolomics) = c("Metabolite","BW_1","BW_2","BW_3","BW_4","BW + Glu_1",
                           "BW + Glu_2","BW + Glu_3","BW + Glu_4","ΔpyrE_1",
                           "ΔpyrE_2","ΔpyrE_3","ΔpyrE_4","ΔpyrE+ Glu_1",
                           "ΔpyrE+ Glu_2","ΔpyrE+ Glu_3","ΔpyrE+ Glu_4")


metabolomics = metabolomics %>% 
  mutate(Metabolite = case_when(Metabolite == 'missing' ~ "Orotate",
                                TRUE ~ Metabolite))

meta_df = metabolomics %>% 
  column_to_rownames("Metabolite")



ha = HeatmapAnnotation(condition = c("BW", "BW", "BW", "BW", 
                                     "BW + Glu", "BW + Glu", "BW + Glu", "BW + Glu", 
                                     "ΔpyrE", "ΔpyrE", "ΔpyrE", "ΔpyrE", 
                                     "ΔpyrE + Glu","ΔpyrE + Glu","ΔpyrE + Glu","ΔpyrE + Glu"),
                       col = list(condition = 
                                    c("BW" = "#E58419", 
                                      "BW + Glu" = "#66533E", 
                                      "ΔpyrE" = "#19AFE6",
                                      "ΔpyrE + Glu" = "#407B91")))

log2(meta_df + 1) %>% 
  Heatmap(name = "Intensity (log2)",
          show_column_names = F,
          cluster_columns = F,
          top_annotation = ha) 

quartz.save(file = here('metabolomics_log2_heatmap.pdf'),
            type = 'pdf', dpi = 300, height = 10, width = 12)

log2(meta_df + 1) %>%  t %>% 
  scale() %>% t %>% 
  Heatmap(name = "Z-score",
          show_column_names = F,
          cluster_columns = F,
          top_annotation = ha) 
  
quartz.save(file = here('metabolomics_zscore_heatmap.pdf'),
            type = 'pdf', dpi = 300, height = 10, width = 12)


  
  
