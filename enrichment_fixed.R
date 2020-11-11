
# EcoCyc  --------------------------------------------------------------


# Enrich bacteria ----------------------------------------------------------------


###
# let's do it with every set

res.BW = res %>% filter(Strain == 'BW', Contrast == 'BW_5FU')
res.pyrE = res %>% filter(Strain == 'pyrE', Contrast == 'pyrE_5FU')
res.TM = res %>% filter(Strain == 'TM', Contrast == 'TM_5FU')


# total metab
met = res.BW %>% filter(Metabolite != 'Negative Control') %>% select(EcoCycID) %>% t %>% as.vector

# significative metabolites 
sig.met.BW.up =   res.BW %>% filter(Metabolite != 'Negative Control', FDR < 0.05, logFC > 0) %>% select(EcoCycID) %>% t %>% as.vector
sig.met.BW.down = res.BW %>% filter(Metabolite != 'Negative Control', FDR < 0.05, logFC < 0) %>% select(EcoCycID) %>% t %>% as.vector

sig.met.pyrE.up =   res.pyrE %>% filter(Metabolite != 'Negative Control', FDR < 0.05, logFC > 0) %>% select(EcoCycID) %>% t %>% as.vector
sig.met.pyrE.down = res.pyrE %>% filter(Metabolite != 'Negative Control', FDR < 0.05, logFC < 0) %>% select(EcoCycID) %>% t %>% as.vector

sig.met.TM.up =   res.TM %>% filter(Metabolite != 'Negative Control', FDR < 0.05, logFC > 0) %>% select(EcoCycID) %>% t %>% as.vector
sig.met.TM.down = res.TM %>% filter(Metabolite != 'Negative Control', FDR < 0.05, logFC < 0) %>% select(EcoCycID) %>% t %>% as.vector


### initialize variables

N = length(met %in% EC_classes$EcoCycID) # total marked elements
classes = unique(EC_classes %>% filter(Plate != 'extra') %>% select(EcoCyc_Classes) %>% t %>% as.character())
### calculate enrichment values
### BW
# up 
enrich.BW.up = c()
for (class in classes){
  class.met = EC_classes %>% filter(EcoCyc_Classes == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.BW.up)
  x = length(class.met[class.met %in% sig.met.BW.up])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.BW.up = c(enrich.BW.up, fit)
}

enrich.BW.up = p.adjust(enrich.BW.up, method = 'fdr')
names(enrich.BW.up) = classes
enrich.BW.up = enrich.BW.up[enrich.BW.up < 0.05]

# down 
enrich.BW.down = c()
for (class in classes){
  class.met = EC_classes %>% filter(EcoCyc_Classes == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.BW.down)
  x = length(class.met[class.met %in% sig.met.BW.down])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.BW.down = c(enrich.BW.down, fit)
}

enrich.BW.down = p.adjust(enrich.BW.down, method = 'fdr')
names(enrich.BW.down) = classes
enrich.BW.down = enrich.BW.down[enrich.BW.down < 0.05]





### pyrE
# up
enrich.pyrE.up = c()
for (class in classes){
  class.met = EC_classes %>% filter(EcoCyc_Classes == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.pyrE.up)
  x = length(class.met[class.met %in% sig.met.pyrE.up])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.pyrE.up = c(enrich.pyrE.up, fit)
}

enrich.pyrE.up = p.adjust(enrich.pyrE.up, method = 'fdr')
names(enrich.pyrE.up) = classes
enrich.pyrE.up = enrich.pyrE.up[enrich.pyrE.up < 0.05]

# down
enrich.pyrE.down = c()
for (class in classes){
  class.met = EC_classes %>% filter(EcoCyc_Classes == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.pyrE.down)
  x = length(class.met[class.met %in% sig.met.pyrE.down])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.pyrE.down = c(enrich.pyrE.down, fit)
}

enrich.pyrE.down = p.adjust(enrich.pyrE.down, method = 'fdr')
names(enrich.pyrE.down) = classes
enrich.pyrE.down = enrich.pyrE.down[enrich.pyrE.down < 0.05]







### TM
# up
enrich.TM.up = c()
for (class in classes){
  class.met = EC_classes %>% filter(EcoCyc_Classes == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.TM.up)
  x = length(class.met[class.met %in% sig.met.TM.up])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.TM.up = c(enrich.TM.up, fit)
}

enrich.TM.up = p.adjust(enrich.TM.up, method = 'fdr')
names(enrich.TM.up) = classes
enrich.TM.up = enrich.TM.up[enrich.TM.up < 0.05]


# down
enrich.TM.down = c()
for (class in classes){
  class.met = EC_classes %>% filter(EcoCyc_Classes == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.TM.down)
  x = length(class.met[class.met %in% sig.met.TM.down])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.TM.down = c(enrich.TM.down, fit)
}

enrich.TM.down = p.adjust(enrich.TM.down, method = 'fdr')
names(enrich.TM.down) = classes
enrich.TM.down = enrich.TM.down[enrich.TM.down < 0.1]



# let's workout the results from each enrichment

enrich.BW.up      = data.frame(enrich.BW.up     ) %>% mutate(Class = rownames(.), Direction =   'up', Strain =   'BW');        
enrich.BW.down    = data.frame(enrich.BW.down   ) %>% mutate(Class = rownames(.), Direction = 'down', Strain =   'BW');       
enrich.pyrE.up    = data.frame(enrich.pyrE.up   ) %>% mutate(Class = rownames(.), Direction =   'up', Strain = 'pyrE');       
enrich.pyrE.down  = data.frame(enrich.pyrE.down ) %>% mutate(Class = rownames(.), Direction = 'down', Strain = 'pyrE');       
enrich.TM.up      = data.frame(enrich.TM.up     ) %>% mutate(Class = rownames(.), Direction =   'up', Strain =   'TM');       
enrich.TM.down    = data.frame(enrich.TM.down   ) %>% mutate(Class = rownames(.), Direction = 'down', Strain =   'TM');       


names(enrich.BW.up)[1] = 'p.value'
names(enrich.BW.down)[1] = 'p.value'
names(enrich.pyrE.up)[1] = 'p.value'
names(enrich.pyrE.down)[1] = 'p.value'
names(enrich.TM.up)[1] = 'p.value'
names(enrich.TM.down)[1] = 'p.value'


bac.enrich = rbind(enrich.BW.up, enrich.BW.down, enrich.pyrE.up, enrich.pyrE.down)

# bac.enrich %>%
#     ggplot(aes(x = Direction, y = Class)) +
#     geom_tile(aes(fill = p.value)) +
#     facet_wrap(~Strain)

# lets build a new dataframe that has all the info we need for plotting
classes = sort(unique(bac.enrich$Class))
direct = c('up', 'down')

df = expand.grid(Strain = c('BW', 'pyrE', 'TM'),
                 Class = classes,
                 Direction = c('up', 'down'))

df = left_join(df, bac.enrich) %>%
  mutate(p.value = replace_na(p.value, 1))

# enrichment procedure
enrbrks = c(0, -log(0.05, 10), 2, 3, 4, 100)
enrlbls = c('N.S.','<0.05','<0.01','<0.001','<0.0001')
enrcols = colorRampPalette(c("gray90", "steelblue1", "blue4"))(n = 6)

p.theme = theme(axis.ticks = element_blank(), panel.border = element_blank(), 
                panel.background = element_blank(), panel.grid.minor = element_blank(), 
                panel.grid.major = element_blank(), axis.line = element_line(colour = NA), 
                axis.line.x = element_line(colour = NA), axis.line.y = element_line(colour = NA), 
                strip.text = element_text(colour = "black", face = "bold", 
                                          size = 10), axis.text.x = element_text(face = "bold", 
                                                                                 colour = "black", size = 10, angle = 90, hjust = 1))


# plot enrichment p-values
df %>% 
  mutate(p.value = p.value + 0.00000001,
         logFDR = ifelse(-log10(p.value) < 0, 0, -log10(p.value)),
         logFDRbin = cut(logFDR, breaks = enrbrks, labels = enrlbls, right = FALSE),
         Class = factor(Class, levels = sort(classes, decreasing = T)),
         Direction = factor(Direction, levels = sort(direct, decreasing = T))) %>%
  ggplot(aes(x = Direction, y = Class)) +
  geom_tile(aes(fill = logFDRbin)) +
  scale_fill_manual(values = enrcols) + 
  facet_wrap(~Strain) +
  # labs(fill = 'FDR') + 
  p.theme


# quartz.save(file = here('Summary', 'Enrich_metab_class_bacteria.pdf'),
#             type = 'pdf', dpi = 300, height = 6, width = 6)

ggsave(file = here('Summary', 'Enrich_metab_class_bacteria.pdf'),
       height = 6, width = 6)











# Enrich Worm -------------------------------------------------------------








### WORM SIDE


### initialize variables

N = length(met %in% EC_classes$EcoCycID) # total marked elements
classes = unique(EC_classes %>% filter(Plate != 'extra') %>% select(EcoCyc_Classes) %>% t %>% as.character())

join.res = res %>% 
  filter(Contrast %in% c('BW_5FU', 'TM_5FU', 'pyrE_5FU')) %>%
  left_join(worm.sum) 


# calculate enrichment values
# BW
# significative metabolites 
res.BW = join.res %>% filter(Strain == 'BW', Experiment == 'T5')
sig.met.BW = res.BW %>% filter(Metabolite != 'Negative Control', Median == 4) %>% select(EcoCycID) %>% t %>% as.vector

classes = unique(EC_classes$EcoCyc_Classes)
enrich.BW = c()
for (class in classes){
  class.met = EC_classes %>% filter(EcoCyc_Classes == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.BW)
  x = length(class.met[class.met %in% sig.met.BW])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.BW = c(enrich.BW, fit)
}

enrich.BW = p.adjust(enrich.BW, method = 'fdr')
names(enrich.BW) = classes
enrich.BW = enrich.BW[enrich.BW < 0.05]




# pyrE

res.pyrE = join.res %>% filter(Strain == 'pyrE', Experiment == 'T5')
sig.met.pyrE = res.pyrE %>% filter(Metabolite != 'Negative Control', Median == 4) %>% select(EcoCycID) %>% t %>% as.vector

classes = unique(EC_classes$EcoCyc_Classes)
enrich.pyrE = c()
for (class in classes){
  class.met = EC_classes %>% filter(EcoCyc_Classes == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.pyrE)
  x = length(class.met[class.met %in% sig.met.pyrE])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.pyrE = c(enrich.pyrE, fit)
}

enrich.pyrE = p.adjust(enrich.pyrE, method = 'fdr')
names(enrich.pyrE) = classes
enrich.pyrE = enrich.pyrE[enrich.pyrE < 0.05]


## TM
res.TM = join.res %>% filter(Strain == 'TM', Experiment == 'T250')
sig.met.TM = res.TM %>% filter(Metabolite != 'Negative Control', Median == 4) %>% select(EcoCycID) %>% t %>% as.vector


classes = unique(EC_classes$EcoCyc_Classes)
enrich.TM = c()
for (class in classes){
  class.met = EC_classes %>% filter(EcoCyc_Classes == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.TM)
  x = length(class.met[class.met %in% sig.met.TM])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.TM = c(enrich.TM, fit)
}

enrich.TM = p.adjust(enrich.TM, method = 'fdr')
names(enrich.TM) = classes
enrich.TM = enrich.TM[enrich.TM < 0.05]




# let's workout the results from each enrichment

enrich.BW      = data.frame(enrich.BW     ) %>% mutate(Class = rownames(.), Strain =   'BW');        
enrich.pyrE    = data.frame(enrich.pyrE   ) %>% mutate(Class = rownames(.), Strain = 'pyrE');       
enrich.TM      = data.frame(enrich.TM     ) %>% mutate(Class = rownames(.), Strain =   'TM');       


names(enrich.BW)[1] = 'p.value'
names(enrich.pyrE)[1] = 'p.value'
names(enrich.TM)[1] = 'p.value'

worm.enrich = rbind(enrich.BW, enrich.pyrE, enrich.TM)


# lets build a new dataframe that has all the info we need for plotting
classes = sort(unique(worm.enrich$Class))

df = expand.grid(Strain = c('BW', 'pyrE', 'TM'),
                 Class = classes)

df = left_join(df, worm.enrich) %>%
  mutate(p.value = replace_na(p.value, 1))

# enrichment procedure
enrbrks = c(0, -log(0.05, 10), 2, 3, 4, 100)
enrlbls = c('N.S.','<0.05','<0.01','<0.001','<0.0001')
enrcols = colorRampPalette(c("gray90", "steelblue1", "blue4"))(n = 6)

p.theme = theme(axis.ticks = element_blank(), panel.border = element_blank(), 
                panel.background = element_blank(), panel.grid.minor = element_blank(), 
                panel.grid.major = element_blank(), axis.line = element_line(colour = NA), 
                axis.line.x = element_line(colour = NA), axis.line.y = element_line(colour = NA), 
                strip.text = element_text(colour = "black", face = "bold", 
                                          size = 10), axis.text.x = element_text(face = "bold", 
                                                                                 colour = "black", size = 10, angle = 90, hjust = 1))



# plot enrichment p-values
df %>% 
  # arrange((Class)) %>%
  mutate(p.value = p.value + 0.00000001,
         logFDR = ifelse(-log10(p.value) < 0, 0, -log10(p.value)),
         logFDRbin = cut(logFDR, breaks = enrbrks, labels = enrlbls, right = FALSE),
         Class = factor(Class, levels = sort(classes, decreasing = T))) %>%
  ggplot(aes(x = Strain, y = Class)) +
  geom_tile(aes(fill = logFDRbin)) +
  scale_fill_manual(values = enrcols) + 
  # facet_wrap(~Strain) +
  # labs(fill = 'FDR') + 
  p.theme
# 
# quartz.save(file = here('Summary', 'Enrich_metab_class_worm.pdf'),
#             type = 'pdf', dpi = 300, height = 6, width = 6)
ggsave(file = here('Summary', 'Enrich_metab_class_worm.pdf'),
       height = 6, width = 6)

















# KEGG --------------------------------------------------------------------


# Bacteria --------------------------------------------------------------------





EC_path = read_csv(here('Bacteria data',"All_results_side_by_side_withKEGGEcoCycClass.csv")) %>%
  select(Plate:EcoCyc_Classes) %>%
  separate_rows(Pathways, sep = ';')

### calculate enrichment values
### BW
# up 
classes = unique(EC_path$Pathways)
classes = classes[-16]
classes = classes[-89]




###
# let's do it with every set

res.BW = res %>% filter(Strain == 'BW', Contrast == 'BW_5FU')
res.pyrE = res %>% filter(Strain == 'pyrE', Contrast == 'pyrE_5FU')
res.TM = res %>% filter(Strain == 'TM', Contrast == 'TM_5FU')


# total metab
met = res.BW %>% filter(Metabolite != 'Negative Control') %>% select(EcoCycID) %>% t %>% as.vector

# significative metabolites 
sig.met.BW.up =   res.BW %>% filter(Metabolite != 'Negative Control', FDR < 0.05, logFC > 0) %>% select(EcoCycID) %>% t %>% as.vector
sig.met.BW.down = res.BW %>% filter(Metabolite != 'Negative Control', FDR < 0.05, logFC < 0) %>% select(EcoCycID) %>% t %>% as.vector

sig.met.pyrE.up =   res.pyrE %>% filter(Metabolite != 'Negative Control', FDR < 0.05, logFC > 0) %>% select(EcoCycID) %>% t %>% as.vector
sig.met.pyrE.down = res.pyrE %>% filter(Metabolite != 'Negative Control', FDR < 0.05, logFC < 0) %>% select(EcoCycID) %>% t %>% as.vector

sig.met.TM.up =   res.TM %>% filter(Metabolite != 'Negative Control', FDR < 0.05, logFC > 0) %>% select(EcoCycID) %>% t %>% as.vector
sig.met.TM.down = res.TM %>% filter(Metabolite != 'Negative Control', FDR < 0.05, logFC < 0) %>% select(EcoCycID) %>% t %>% as.vector


enrich.BW.up = c()
for (class in classes){
  class.met = EC_path %>% filter(Pathways == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.BW.up)
  x = length(class.met[class.met %in% sig.met.BW.up])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.BW.up = c(enrich.BW.up, fit)
}

enrich.BW.up = p.adjust(enrich.BW.up, method = 'fdr')
names(enrich.BW.up) = classes
enrich.BW.up = enrich.BW.up[enrich.BW.up < 0.05]

# down 
enrich.BW.down = c()
for (class in classes){
  class.met = EC_path %>% filter(Pathways == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.BW.down)
  x = length(class.met[class.met %in% sig.met.BW.down])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.BW.down = c(enrich.BW.down, fit)
}

enrich.BW.down = p.adjust(enrich.BW.down, method = 'fdr')
names(enrich.BW.down) = classes
enrich.BW.down = enrich.BW.down[enrich.BW.down < 0.05]





### pyrE
# up
classes = unique(EC_path$Pathways)
enrich.pyrE.up = c()
for (class in classes){
  class.met = EC_path %>% filter(Pathways == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.pyrE.up)
  x = length(class.met[class.met %in% sig.met.pyrE.up])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.pyrE.up = c(enrich.pyrE.up, fit)
}

enrich.pyrE.up = p.adjust(enrich.pyrE.up, method = 'fdr')
names(enrich.pyrE.up) = classes
enrich.pyrE.up = enrich.pyrE.up[enrich.pyrE.up < 0.05]

# down
classes = unique(EC_path$Pathways)
enrich.pyrE.down = c()
for (class in classes){
  class.met = EC_path %>% filter(Pathways == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.pyrE.down)
  x = length(class.met[class.met %in% sig.met.pyrE.down])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.pyrE.down = c(enrich.pyrE.down, fit)
}

enrich.pyrE.down = p.adjust(enrich.pyrE.down, method = 'fdr')
names(enrich.pyrE.down) = classes
enrich.pyrE.down = enrich.pyrE.down[enrich.pyrE.down < 0.05]







### TM
# up
classes = unique(EC_path$Pathways)
enrich.TM.up = c()
for (class in classes){
  class.met = EC_path %>% filter(Pathways == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.TM.up)
  x = length(class.met[class.met %in% sig.met.TM.up])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.TM.up = c(enrich.TM.up, fit)
}

enrich.TM.up = p.adjust(enrich.TM.up, method = 'fdr')
names(enrich.TM.up) = classes
enrich.TM.up = enrich.TM.up[enrich.TM.up < 0.05]


# down
classes = unique(EC_path$Pathways)
enrich.TM.down = c()
for (class in classes){
  class.met = EC_path %>% filter(Pathways == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.TM.down)
  x = length(class.met[class.met %in% sig.met.TM.down])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.TM.down = c(enrich.TM.down, fit)
}

enrich.TM.down = p.adjust(enrich.TM.down, method = 'fdr')
names(enrich.TM.down) = classes
enrich.TM.down = enrich.TM.down[enrich.TM.down < 0.05]



# let's workout the results from each enrichment

enrich.BW.up      = data.frame(enrich.BW.up     ) %>% mutate(Class = rownames(.), Direction =   'up', Strain =   'BW');        
enrich.BW.down    = data.frame(enrich.BW.down   ) %>% mutate(Class = rownames(.), Direction = 'down', Strain =   'BW');       
enrich.pyrE.up    = data.frame(enrich.pyrE.up   ) %>% mutate(Class = rownames(.), Direction =   'up', Strain = 'pyrE');       
enrich.pyrE.down  = data.frame(enrich.pyrE.down ) %>% mutate(Class = rownames(.), Direction = 'down', Strain = 'pyrE');       
enrich.TM.up      = data.frame(enrich.TM.up     ) %>% mutate(Class = rownames(.), Direction =   'up', Strain =   'TM');       
enrich.TM.down    = data.frame(enrich.TM.down   ) %>% mutate(Class = rownames(.), Direction = 'down', Strain =   'TM');       


names(enrich.BW.up)[1] = 'p.value'
names(enrich.BW.down)[1] = 'p.value'
names(enrich.pyrE.up)[1] = 'p.value'
names(enrich.pyrE.down)[1] = 'p.value'
names(enrich.TM.up)[1] = 'p.value'
names(enrich.TM.down)[1] = 'p.value'


bac.enrich = rbind(enrich.BW.up, enrich.BW.down, enrich.pyrE.up, enrich.pyrE.down)

# bac.enrich %>%
#     ggplot(aes(x = Direction, y = Class)) +
#     geom_tile(aes(fill = p.value)) +
#     facet_wrap(~Strain)

# lets build a new dataframe that has all the info we need for plotting
classes = sort(unique(bac.enrich$Class))
direct = c('up', 'down')

df = expand.grid(Strain = c('BW', 'pyrE', 'TM'),
                 Class = classes,
                 Direction = c('up', 'down'))

df = left_join(df, bac.enrich) %>%
  mutate(p.value = replace_na(p.value, 1))

# enrichment procedure
enrbrks = c(0, -log(0.05, 10), 2, 3, 4, 100)
enrlbls = c('N.S.','<0.05','<0.01','<0.001','<0.0001')
enrcols = colorRampPalette(c("gray90", "steelblue1", "blue4"))(n = 6)

p.theme = theme(axis.ticks = element_blank(), panel.border = element_blank(), 
                panel.background = element_blank(), panel.grid.minor = element_blank(), 
                panel.grid.major = element_blank(), axis.line = element_line(colour = NA), 
                axis.line.x = element_line(colour = NA), axis.line.y = element_line(colour = NA), 
                strip.text = element_text(colour = "black", face = "bold", 
                                          size = 10), axis.text.x = element_text(face = "bold", 
                                                                                 colour = "black", size = 10, angle = 90, hjust = 1))


# plot enrichment p-values
df %>% 
  mutate(p.value = p.value + 0.00000001,
         logFDR = ifelse(-log10(p.value) < 0, 0, -log10(p.value)),
         logFDRbin = cut(logFDR, breaks = enrbrks, labels = enrlbls, right = FALSE),
         Class = factor(Class, levels = sort(classes, decreasing = T)),
         Direction = factor(Direction, levels = sort(direct, decreasing = T))) %>%
  ggplot(aes(x = Direction, y = Class)) +
  geom_tile(aes(fill = logFDRbin)) +
  scale_fill_manual(values = enrcols) + 
  facet_wrap(~Strain) +
  # labs(fill = 'FDR') + 
  p.theme


# quartz.save(file = here('Summary', 'Enrich_pathways_bacteria.pdf'),
            # type = 'pdf', dpi = 300, height = 6, width = 6)

ggsave(file = here('Summary', 'Enrich_pathways_bacteria.pdf'),
       height = 6, width = 6)



# Worm --------------------------------------------------------------------




### WORM SIDE

join.res = res %>% 
  filter(Contrast %in% c('BW_5FU', 'TM_5FU', 'pyrE_5FU')) %>%
  left_join(worm.sum) 


# calculate enrichment values
# BW
# significative metabolites 
res.BW = join.res %>% filter(Strain == 'BW', Experiment == 'T5')
sig.met.BW = res.BW %>% filter(Metabolite != 'Negative Control', Median == 4) %>% select(EcoCycID) %>% t %>% as.vector

classes = unique(EC_path$Pathways)

enrich.BW = c()
for (class in classes){
  class.met = EC_path %>% filter(Pathways == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.BW)
  x = length(class.met[class.met %in% sig.met.BW])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.BW = c(enrich.BW, fit)
}

enrich.BW = p.adjust(enrich.BW, method = 'fdr')
names(enrich.BW) = classes
enrich.BW = enrich.BW[enrich.BW < 0.05]




# pyrE

res.pyrE = join.res %>% filter(Strain == 'pyrE', Experiment == 'T5')
sig.met.pyrE = res.pyrE %>% filter(Metabolite != 'Negative Control', Median == 4) %>% select(EcoCycID) %>% t %>% as.vector

enrich.pyrE = c()
for (class in classes){
  class.met = EC_path %>% filter(Pathways == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.pyrE)
  x = length(class.met[class.met %in% sig.met.pyrE])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.pyrE = c(enrich.pyrE, fit)
}

enrich.pyrE = p.adjust(enrich.pyrE, method = 'fdr')
names(enrich.pyrE) = classes
enrich.pyrE = enrich.pyrE[enrich.pyrE < 0.05]

# TM
res.TM = join.res %>% filter(Strain == 'TM', Experiment == 'T250')
sig.met.TM = res.TM %>% filter(Metabolite != 'Negative Control', Median == 4) %>% select(EcoCycID) %>% t %>% as.vector


enrich.TM = c()
for (class in classes){
  class.met = EC_path %>% filter(Pathways == class) %>% select(EcoCycID) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(sig.met.TM)
  x = length(class.met[class.met %in% sig.met.TM])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.TM = c(enrich.TM, fit)
}

enrich.TM = p.adjust(enrich.TM, method = 'fdr')
names(enrich.TM) = classes
enrich.TM = enrich.TM[enrich.TM < 0.05]




# let's workout the results from each enrichment

enrich.BW      = data.frame(enrich.BW     ) %>% mutate(Class = rownames(.), Strain =   'BW');        
enrich.pyrE    = data.frame(enrich.pyrE   ) %>% mutate(Class = rownames(.), Strain = 'pyrE');       
enrich.TM      = data.frame(enrich.TM     ) %>% mutate(Class = rownames(.), Strain =   'TM');       


names(enrich.BW)[1] = 'p.value'
names(enrich.pyrE)[1] = 'p.value'
names(enrich.TM)[1] = 'p.value'

worm.enrich = rbind(enrich.BW, enrich.pyrE, enrich.TM)


# lets build a new dataframe that has all the info we need for plotting
classes = sort(unique(worm.enrich$Class))

df = expand.grid(Strain = c('BW', 'pyrE', 'TM'),
                 Class = classes)

df = left_join(df, worm.enrich) %>%
  mutate(p.value = replace_na(p.value, 1))

# enrichment procedure
enrbrks = c(0, -log(0.05, 10), 2, 3, 4, 100)
enrlbls = c('N.S.','<0.05','<0.01','<0.001','<0.0001')
enrcols = colorRampPalette(c("gray90", "steelblue1", "blue4"))(n = 6)

p.theme = theme(axis.ticks = element_blank(), panel.border = element_blank(), 
                panel.background = element_blank(), panel.grid.minor = element_blank(), 
                panel.grid.major = element_blank(), axis.line = element_line(colour = NA), 
                axis.line.x = element_line(colour = NA), axis.line.y = element_line(colour = NA), 
                strip.text = element_text(colour = "black", face = "bold", 
                                          size = 10), axis.text.x = element_text(face = "bold", 
                                                                                 colour = "black", size = 10, angle = 90, hjust = 1))



# plot enrichment p-values
df %>% 
  # arrange((Class)) %>%
  mutate(p.value = p.value + 0.00000001,
         logFDR = ifelse(-log10(p.value) < 0, 0, -log10(p.value)),
         logFDRbin = cut(logFDR, breaks = enrbrks, labels = enrlbls, right = FALSE),
         Class = factor(Class, levels = sort(classes, decreasing = T))) %>%
  ggplot(aes(x = Strain, y = Class)) +
  geom_tile(aes(fill = logFDRbin)) +
  scale_fill_manual(values = enrcols) + 
  # facet_wrap(~Strain) +
  # labs(fill = 'FDR') + 
  p.theme

# quartz.save(file = here('Summary', 'Enrich_pathways_worm.pdf'),
#             type = 'pdf', dpi = 300, height = 6, width = 6)

ggsave(file = here('Summary', 'Enrich_pathways_worm.pdf'),
       height = 6, width = 6)












