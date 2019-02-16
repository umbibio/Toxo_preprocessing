library(openxlsx)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

abs.path <- '~/work/ToxoplasmaGondii/'
gse.file <- paste(abs.path, 'GSEA 2.12.19.xlsx', sep = '')
GeneSet.190 <- read.xlsx(gse.file)
## Make less lenghty names
colnames(GeneSet.190) <- 
  gsub("development", "devel", 
       gsub('tgo.*[[:digit:]]\\.', '', 
            gsub('Tachyzoite', 'Tachy', 
                 gsub('Bradyzoite', 'Brady', 
                      gsub('tachyzoite', 'Tachy',  
                           gsub('bradyzoite', 'Brady',
                                gsub('\\,', '', 
                                     gsub('\\(.*', '', 
                                          gsub('_', '.', 
                                               gsub('-', '.', colnames(GeneSet.190)))))))))))

GeneSet.190 <- GeneSet.190 %>% gather(key = GeneSet, value = GeneName) %>% na.omit() %>% 
  group_by(GeneSet) %>% summarise(genes = list(as.character(GeneName)), total = n())

toxo.tab.file <- paste(abs.path, 'toxo_table_batch_corrected_logCPM_expression_edgeR_DEGs_ap2_targs.xlsx', sep = '')
toxo.tab <- read.xlsx(toxo.tab.file)


#### Get FC and Qvalue in tidy format
toxo.fc   <- toxo.tab %>% dplyr::select(c(1, contains('fc')))
colnames(toxo.fc) <- gsub("\\._log_fc", "", colnames(toxo.fc))
toxo.qval <- toxo.tab %>% dplyr::select(c(1, contains('qval')))
colnames(toxo.qval) <- gsub("\\._qval", "", colnames(toxo.qval))

toxo.fc <- toxo.fc  %>% gather(key = Contrast, value = fc, -GeneName)
toxo.qval <- toxo.qval  %>% gather(key = Contrast, value = qval, -GeneName)

toxo.fc$fc <- as.numeric(toxo.fc$fc)
toxo.qval$qval <- as.numeric(toxo.qval$qval)

toxo.fc.qval <- inner_join(toxo.fc, toxo.qval, by = c('GeneName', 'Contrast'))

Passages <- strsplit(as.character(toxo.fc.qval$Contrast), split='\\.vs\\.')
Passage2 <- unlist(lapply(strsplit(unlist(lapply(Passages, `[[`, 1)), split='\\.'), `[[`,1))
Passage1 <- unlist(lapply(strsplit(unlist(lapply(Passages, `[[`, 2)), split='\\.'), `[[`,1))

Category <- gsub("RH.intra", "RH",
                   gsub("RH.extra", "RH", 
                        gsub("^P.*[[:digit:]]\\.", "",
                             gsub("\\.P.*[[:digit:]]\\.", "\\.", 
                                  gsub('P7', 'B', as.character(toxo.fc$Contrast))))))

sorted.ind1 <- sort(as.numeric(gsub("P", "", gsub('RH', 'P300', unique(Passage1)))), index.return = T)$ix
sorted.ind2 <- sort(as.numeric(gsub("P", "", gsub('RH', 'P300', unique(Passage2)))), index.return = T)$ix
Passage1 <- factor((Passage1), levels = unique(Passage1)[sorted.ind1])
Passage2 <- factor((Passage2), levels = unique(Passage2)[sorted.ind2])
Category <- factor(Category, levels = sort(unique(Category)))

toxo.fc.qval$Category <- Category
toxo.fc.qval$Passage1 <- Passage1
toxo.fc.qval$Passage2 <- Passage2
toxo.fc.qval$Contrast <- factor(toxo.fc.qval$Contrast, levels = unique(toxo.fc.qval$Contrast))

toxo.fc.qval <- toxo.fc.qval %>% na.omit()


toxo.fc.qval.sig <- toxo.fc.qval %>% dplyr::filter(abs(fc) > log2(1.5) & qval < 0.05)
toxo.fc.qval.up.reg <- toxo.fc.qval.sig %>% dplyr::filter(fc > 0)
toxo.fc.qval.down.reg <- toxo.fc.qval.sig %>% dplyr::filter(fc < 0)


getEnrichment <- function(GeneSet, GeneSet.list){
  ## This is a Hack to avoid nested for loops
  GeneSet$all <- 8637
  GeneSet.list$all <- 8637
  
  XX <- full_join(GeneSet, GeneSet.list, by = 'all')
  XX <- XX %>% rowwise() %>% mutate(overlap = length(intersect(c(genes.x), c(genes.y))),
                                    overlap.genes = list(intersect(c(genes.x), c(genes.y))))
  XX <- XX %>% rowwise() %>% 
    mutate(pvalue = fisher.test(matrix(c(overlap, total.x - overlap, total.y - overlap, 
                                         all - (total.x + total.y - overlap) ), byrow = T, ncol = 2, nrow = 2),
                                alternative = "greater")$p.value)
  #XX <- XX %>% dplyr::select(c('GeneSet', 'Contrast', 'pvalue'))
  
  contrasts <- unique(XX$Contrast)
  
  extra.vs.extra.ind <- grep('(P.*extra)(.*P.*extra)',contrasts)
  
  passage2 <- unlist(lapply(strsplit(as.character(contrasts[extra.vs.extra.ind]), split = '\\.vs\\.'), `[[`, 1))
  passage1 <- unlist(lapply(strsplit(as.character(contrasts[extra.vs.extra.ind]), split = '\\.vs\\.'), `[[`, 2))
  dd <- data.frame(p1 = as.numeric(gsub('P', '',gsub('\\.extra', '', passage1))), 
             p2 = as.numeric(gsub('P', '',gsub('\\.extra', '', passage2)))) 
  rownames(dd) <- 1:nrow(dd)
  sorted.extra.vs.extra.ind <- extra.vs.extra.ind[with(dd, order(p1, p2))]
  #sorted.extra.vs.extra.ind <- extra.vs.extra.ind[
  #  sort(as.numeric(gsub('P', '', gsub('\\.extra.*', '',contrasts[extra.vs.extra.ind]))), index.return = T)$ix]
  
  intra.vs.intra.ind <- grep('(P.*intra)(.*P.*intra)',contrasts)
  passage2 <- unlist(lapply(strsplit(as.character(contrasts[intra.vs.intra.ind]), split = '\\.vs\\.'), `[[`, 1))
  passage1 <- unlist(lapply(strsplit(as.character(contrasts[intra.vs.intra.ind]), split = '\\.vs\\.'), `[[`, 2))
  dd <- data.frame(p1 = as.numeric(gsub('P', '',gsub('\\.intra', '', passage1))), 
                   p2 = as.numeric(gsub('P', '',gsub('\\.intra', '', passage2)))) 
  rownames(dd) <- 1:nrow(dd)
  sorted.intra.vs.intra.ind <- intra.vs.intra.ind[with(dd, order(p1, p2))]
  
  #sorted.intra.vs.intra.ind <- intra.vs.intra.ind[
  #  sort(as.numeric(gsub('P', '', gsub('\\.intra.*', '',contrasts[intra.vs.intra.ind]))), index.return = T)$ix]
  
  extra.vs.intra.ind <- grep('(P.*extra)(.*P.*intra)',contrasts)
  sorted.extra.vs.intra.ind <- extra.vs.intra.ind[
    sort(as.numeric(gsub('P', '', gsub('\\.extra.*', '',contrasts[extra.vs.intra.ind]))), index.return = T)$ix]
  
  
  RH.vs.RH.ind <- grep('(RH.*)(RH.*)',contrasts)
  
  RH.vs.intra.ind <- grep('(RH\\.intra)(.*P.*intra)',contrasts)
  sorted.RH.vs.intra.ind <- RH.vs.intra.ind[
    sort(as.numeric(gsub('P', '', 
                         gsub('\\.intra.*', '',
                              gsub('RH\\.intra\\.vs\\.', '', contrasts[RH.vs.intra.ind])))), 
         index.return = T)$ix]
  
  RH.vs.extra.ind <- grep('(RH\\.extra)(.*P.*extra)',contrasts)
  sorted.RH.vs.extra.ind <- RH.vs.extra.ind[
    sort(as.numeric(gsub('P', '', 
                         gsub('\\.extra.*', '',
                              sub('RH\\.extra\\.vs\\.', '', contrasts[RH.vs.extra.ind])))), 
         index.return = T)$ix]
  
  my.inds <- c(sorted.extra.vs.extra.ind, sorted.intra.vs.intra.ind, sorted.extra.vs.intra.ind, 
               sorted.RH.vs.extra.ind, sorted.RH.vs.intra.ind, RH.vs.RH.ind )
  XX$Contrast <- factor(contrasts, levels = contrasts[my.inds])
  GeneSets <- XX$GeneSet
  
  categories <- gsub("RH.intra", "RH",gsub("RH.extra", "RH", 
                                           gsub("^P.*[[:digit:]]\\.", "", 
                                                gsub("\\.P.*[[:digit:]]\\.", "\\.", 
                                                     as.character(XX$Contrast)))))
  categories <- factor(categories, levels = unique(categories))
  XX$Category <- categories
  
  XX$GeneSet <- factor(GeneSets, levels = unique(GeneSets))
  #XX <- XX %>% dplyr::filter(pvalue < 0.01) %>% dplyr::filter(Category != 'RH.vs.RH') %>% arrange(Contrast, GeneSet)
  #XX <- XX %>% dplyr::filter(pvalue < 0.01) %>% arrange(Contrast, GeneSet)
  XX <- XX %>%  arrange(Contrast, GeneSet)
  XX$pvalue <- as.numeric(XX$pvalue)
  colnames(XX) <- c("GeneSet", "genes.in.set", "total.genes.in.set", "all.genes",
                    "Contrast", "genes.in.contrast", "total.genes.in.contrast", "overlap",
                    "overlap.genes","pvalue","Category")
  
  return(XX)
}


GeneSet.list.up <- toxo.fc.qval.up.reg %>% group_by(Contrast) %>% 
  summarise(genes = list(as.character(GeneName)), total = n())

GeneSet.list.down <- toxo.fc.qval.down.reg %>% group_by(Contrast) %>% 
  summarise(genes = list(as.character(GeneName)), total = n())

enrich.up <- getEnrichment(GeneSet.190, GeneSet.list.up)
enrich.down <- getEnrichment(GeneSet.190, GeneSet.list.down)

out.up.file <- paste(abs.path, "GeneSet_enrichment_up_genes.xlsx", sep = '')
write.xlsx(x = enrich.up, file = out.up.file , sheetName = "overlaps", row.names = F, col.names = T)
out.down.file <- paste(abs.path, "GeneSet_enrichment_down_genes.xlsx", sep = '')
write.xlsx(x = enrich.down, file = out.down.file, sheetName = "overlaps", row.names = F, col.names = T)


p1 <- ggplot(subset(enrich.up,pvalue < 0.01), aes(x = GeneSet, y = Contrast)) + 
  geom_point(aes(color = pvalue, size = -log(pvalue))) +
  theme_bw(base_size = 12) +
  scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + 
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2,size = 8)) + facet_grid(Category ~ ., scales='free') + 
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1, linetype="solid"))

plot(p1)

out.pic <- paste(abs.path, "GeneSet_enrichment_up_genes.pdf", sep = '')
ggsave(filename=out.pic, plot=p1,
       width = 12, height = 12, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)



p2 <- ggplot(subset(enrich.down, pvalue < 0.01), aes(x = GeneSet, y = Contrast)) + 
  geom_point(aes(color = pvalue, size = -log(pvalue))) +
  theme_bw(base_size = 12) +
  scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + 
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2, size = 8)) + facet_grid(Category ~ ., scales='free') + 
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1, linetype="solid"))

plot(p2)
out.pic <- paste(abs.path, "GeneSet_enrichment_down_genes.pdf", sep = '')
ggsave(filename=out.pic, plot=p2,
       width = 12, height = 12, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)



####
## Checking APIAp2 targets

AP2.targets <- data.frame(GeneName = toxo.tab$GeneName, toxo.tab[,grep('AP2.*_genes', colnames(toxo.tab))])
#AP2.targets <- AP2.targets %>% dplyr::select(-c('AP2X5_ikD_upregulated_genes')) %>%
#colnames(AP2.targets) <- gsub('XI5', 'XI-5',gsub('\\.', '-',gsub("_downregulated_genes", '', 
#                                                                 gsub('_target_genes', '',colnames(AP2.targets)))))

AP2.targets <- AP2.targets %>% 
  dplyr::select(-c('AP2X5_ikD_upregulated_genes')) %>%
  dplyr::select(-c('AP2X5_ikD_downregulated_genes')) %>%
  dplyr::select(-c('AP2IX.9_target_genes')) %>%
  dplyr::select(-c("AP2IX9_target_genes_up")) ## they play no role
colnames(AP2.targets) <- gsub("AP2IX9_down", "AP2IX9", 
                              gsub('\\.', '', gsub('_target_genes', '',colnames(AP2.targets))))


AP2.targets <- AP2.targets %>% gather(key = GeneSet, value = targets, -GeneName)
AP2.targets <- AP2.targets %>% group_by(GeneSet) %>% 
  summarise(genes = list(as.character(GeneName[!is.na(targets)])), total = sum(!is.na(targets)))

AP2.targets$GeneSet <- factor(AP2.targets$GeneSet, levels = unique(AP2.targets$GeneSet))

AP2.enrich.up <- getEnrichment(AP2.targets, GeneSet.list.up)
AP2.enrich.down <- getEnrichment(AP2.targets, GeneSet.list.down)

ap2.up.file <- paste(abs.path, "AP2_enrichment_up_genes.xlsx")
write.xlsx(x = AP2.enrich.up, file = ap2.up.file, sheetName = "overlaps", row.names = F, col.names = T)
ap2.down.file <- paste(abs.path, "AP2_enrichment_down_genes.xlsx")
write.xlsx(x = AP2.enrich.down, file = ap2.down.file, sheetName = "overlaps", row.names = F, col.names = T)


p1 <- ggplot(subset(AP2.enrich.up, pvalue < 0.05), aes(x = GeneSet, y = Contrast)) + 
  geom_point(aes(color = pvalue, size = -log(pvalue))) +
  theme_bw(base_size = 12) +
  scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + 
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2, size = 8)) + facet_grid(Category ~ ., scales='free') + 
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1, linetype="solid"))

plot(p1)

out.pic <- paste(abs.path, "AP2_enrichment_up_genes.pdf", sep = '')
ggsave(filename=out.pic, plot=p1,
       width = 12, height = 12, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)



p2 <- ggplot(subset(AP2.enrich.down, pvalue < 0.01), aes(x = GeneSet, y = Contrast)) + 
  geom_point(aes(color = pvalue, size = -log(pvalue))) +
  theme_bw(base_size = 12) +
  scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + 
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2, size = 8)) + facet_grid(Category ~ ., scales='free') + 
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1, linetype="solid"))

plot(p2)

out.pic <- paste(abs.path, "AP2_enrichment_down_genes.pdf", sep = '')
ggsave(filename=out.pic, plot=p2,
       width = 12, height = 12, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

### Perform a Three Way overlap between DEGs, AP2 targets, and GeneSets
XX <- GeneSet.190 %>% unnest(genes) %>% dplyr::select(-c('total'))
colnames(XX) <- c('GeneSet', 'GeneName')
XX <- full_join(XX, toxo.fc.qval.down.reg, by = 'GeneName')

YY <- AP2.targets %>% unnest(genes) %>% dplyr::select(-c('total'))
colnames(YY) <- c('AP2', 'GeneName')
ThreeWay.down <- full_join(XX, YY, by = 'GeneName')
ThreeWay.down$Product.Description <- toxo.tab$Product.Description[match(ThreeWay.down$GeneName, toxo.tab$GeneName)]
ThreeWay.down <- ThreeWay.down[!is.na(ThreeWay.down$fc), ]

## repeate fro up genes
XX <- GeneSet.190 %>% unnest(genes) %>% dplyr::select(-c('total'))
colnames(XX) <- c('GeneSet', 'GeneName')
XX <- full_join(XX, toxo.fc.qval.up.reg, by = 'GeneName')
YY <- AP2.targets %>% unnest(genes) %>% dplyr::select(-c('total'))
colnames(YY) <- c('AP2', 'GeneName')
ThreeWay.up <- full_join(XX, YY, by = 'GeneName')
ThreeWay.up$Product.Description <- toxo.tab$Product.Description[match(ThreeWay.up$GeneName, toxo.tab$GeneName)]
ThreeWay.up <- ThreeWay.up[!is.na(ThreeWay.up$fc), ]
ThreeWay <- rbind(ThreeWay.down, ThreeWay.up)

ThreeWay.spread.by.AP2 <- ThreeWay %>% 
  mutate(is.targ = 'yes') %>% distinct() %>% 
  spread(key = AP2, value = is.targ, fill = 'no') %>% dplyr::select(-c('<NA>')) %>%
  dplyr::filter(!is.na(GeneSet))

ThreeWay.spread.by.GeneSet <- ThreeWay %>% 
  mutate(is.targ = 'yes') %>% distinct() %>% 
  spread(key = GeneSet, value = is.targ, fill = 'no') %>% dplyr::select(-c('<NA>')) %>%
  dplyr::filter(!is.na(AP2))

ThreeWay.spread.by.both <- ThreeWay %>% 
  mutate(is.targ = 'yes') %>% distinct() %>% 
  spread(key = GeneSet, value = is.targ, fill = 'no') %>% dplyr::select(-c('<NA>')) %>%
  mutate(is.targ = 'yes') %>% distinct() %>% 
  spread(key = AP2, value = is.targ, fill = 'no') %>% dplyr::select(-c('<NA>'))

threeway.file <- paste(abs.path, "ThreeWayDEGoverlaps.xlsx", sep = '')
write.xlsx(x = ThreeWay, file = threeway.file , sheetName = "overlaps", row.names = F, col.names = T)
threeway.file.ap2 <- paste(abs.path, "ThreeWayDEGoverlaps_Spread_by_AP2.xlsx", sep = '')
write.xlsx(x = ThreeWay.spread.by.AP2, file = threeway.file.ap2, sheetName = "overlaps", row.names = F, col.names = T)
threeway.file.gse <- paste(abs.path, "ThreeWayDEGoverlaps_Spread_by_GeneSet.xlsx", sep = '')
write.xlsx(x = ThreeWay.spread.by.GeneSet , file = threeway.file.gse, sheetName = "overlaps", row.names = F, col.names = T)
threeway.file.both <- paste(abs.path, "ThreeWayDEGoverlaps_Spread_by_both.xlsx", sep = '')
write.xlsx(x = ThreeWay.spread.by.both, file = threeway.file.both, sheetName = "overlaps", row.names = F, col.names = T)
