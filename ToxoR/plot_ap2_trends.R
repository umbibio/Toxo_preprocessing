###
library(dplyr)
library(ggplot2)
library(tidyverse)

toxo.tab <- read.xlsx('../toxo_table_batch_corrected_logCPM_expression_edgeR_DEGs_ap2_targs.xlsx')

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
                                as.character(toxo.fc$Contrast)))))

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


## Get expressions for extra cellular
toxo.tab.extra <- toxo.tab %>% 
  dplyr::select(contains('extra_'))
colnames(toxo.tab.extra) <- gsub('extra_', '', colnames(toxo.tab.extra))

toxo.tab.extra.mean <- toxo.tab.extra %>% dplyr::select(contains('mean'))
colnames(toxo.tab.extra.mean) <- gsub('\\.mean', '', colnames(toxo.tab.extra.mean))
toxo.tab.extra.mean$GeneName <- toxo.tab$GeneName
toxo.tab.extra.mean <- toxo.tab.extra.mean %>% gather(key=Passage, value= mean, -GeneName)

toxo.tab.extra.sd <- toxo.tab.extra %>% dplyr::select(contains('sd'))
colnames(toxo.tab.extra.sd) <- gsub('\\.sd', '', colnames(toxo.tab.extra.sd))
toxo.tab.extra.sd$GeneName <- toxo.tab$GeneName
toxo.tab.extra.sd <- toxo.tab.extra.sd %>% gather(key=Passage, value= sd, -GeneName)

toxo.tab.extra.mean.sd <- left_join(toxo.tab.extra.mean, toxo.tab.extra.sd, by = c('GeneName', 'Passage'))

toxo.tab.extra.mean.sd$cond <- 'extra'

## Same for Intra
toxo.tab.intra <- toxo.tab %>% 
  dplyr::select(contains('intra_'))
colnames(toxo.tab.intra) <- gsub('intra_', '', colnames(toxo.tab.intra))

toxo.tab.intra.mean <- toxo.tab.intra %>% dplyr::select(contains('mean'))
colnames(toxo.tab.intra.mean) <- gsub('\\.mean', '', colnames(toxo.tab.intra.mean))
toxo.tab.intra.mean$GeneName <- toxo.tab$GeneName
toxo.tab.intra.mean <- toxo.tab.intra.mean %>% gather(key=Passage, value= mean, -GeneName)

toxo.tab.intra.sd <- toxo.tab.intra %>% dplyr::select(contains('sd'))
colnames(toxo.tab.intra.sd) <- gsub('\\.sd', '', colnames(toxo.tab.intra.sd))
toxo.tab.intra.sd$GeneName <- toxo.tab$GeneName
toxo.tab.intra.sd <- toxo.tab.intra.sd %>% gather(key=Passage, value= sd, -GeneName)

toxo.tab.intra.mean.sd <- left_join(toxo.tab.intra.mean, toxo.tab.intra.sd, by = c('GeneName', 'Passage'))

toxo.tab.intra.mean.sd$cond <- 'intra'


toxo.tab.mean.sd <- rbind(toxo.tab.extra.mean.sd, toxo.tab.intra.mean.sd)


toxo.tab.means.sd.fc.qval <- left_join(toxo.tab.mean.sd, toxo.fc.qval, 
                                        by = c('GeneName', 'Passage' = "Passage2")) %>%
  dplyr::select(-c('Passage1'))  %>% distinct() %>% na.omit()



AP2s <- toxo.tab %>% dplyr::filter(str_detect(Product.Description, 'AP2')) %>%
  dplyr::select('GeneName', 'Product.Description')
AP2s$Product.Description <- gsub('AP2 domain transcription factor ', '', AP2s$Product.Description)
my.AP2s <- c('AP2XI-5', 'AP2X-5','AP2IX-9', 'AP2IX-4', 'AP2IV-4', 'AP2IV-3', 'AP2XI-4', 'AP2Ib-1', 'AP2XII-6')
#my.AP2s <- c('AP2Ib-1', 'AP2IV-3', 'AP2VI-3', 'AP2VIIa-1', 'AP2VIII-4', 'AP2IX-9', 
#             'AP2XI-5', 'AP2X-5', 'AP2IX-4','AP2IV-4', 'AP2XI-4')
#my.AP2s <- c('AP2Ib-1', 'AP2IV-3', 'AP2IX-9', 
#             'AP2XI-5', 'AP2X-5', 'AP2IX-4','AP2IV-4', 'AP2XI-4')

## AP2s associated with Bradyzoites
#activators <- c('AP2XII-6', 'AP2XI-4','AP2IV-3', 'AP2Ib-1')
#repressors <- c('AP2IX-9','AP2IX-4','AP2IV-4')
#others <- c('AP2XI-5', 'AP2X-5')
#stress_induced <- c('AP2XI-5', 'AP2X-5','AP2IX-9', 'AP2IX-4', 'AP2IV-4', 'AP2IV-3', 'AP2XI-4')
#my.AP2s <- c(others, activators, repressors)

my.AP2s <- AP2s %>% dplyr::filter(Product.Description %in% my.AP2s)
colnames(my.AP2s) <- c('GeneName', 'AP2')

Ap2.tab.means.sd.qval <- right_join(toxo.tab.means.sd.fc.qval, my.AP2s, by = 'GeneName')

Ap2.tab.means.sd.qval$GeneName <- factor(Ap2.tab.means.sd.qval$GeneName,
                                         levels = unique(Ap2.tab.means.sd.qval$GeneName))



#### Filtering out RH and P7 and compare extra vs intra
Ap2.tab.means.sd.qval <- Ap2.tab.means.sd.qval %>% 
  dplyr::filter(!(Passage %in% c('RH', 'P7')))

#Passage <- as.numeric(gsub('P', '', gsub('RH', 'P300', as.character(Ap2.tab.means.sd.qval$Passage))))
Passage <- as.numeric(gsub('P', '', as.character(Ap2.tab.means.sd.qval$Passage)))

ind <- sort(unique(Passage), index.return = T)$ix
Ap2.tab.means.sd.qval$Passage <- factor(Ap2.tab.means.sd.qval$Passage,
                                        levels = unique(Ap2.tab.means.sd.qval$Passage)[ind])

Ap2.tab.means.sd.qval$cond <- factor(Ap2.tab.means.sd.qval$cond,
                                     levels = unique(Ap2.tab.means.sd.qval$cond))

Ap2.tab.means.sd.qval$AP2 <- factor(Ap2.tab.means.sd.qval$AP2,
                                    levels = unique(Ap2.tab.means.sd.qval$AP2))
Ap2.tab.means.sd.qval <- Ap2.tab.means.sd.qval %>% 
  mutate(sig =ifelse(Category == 'extra.vs.intra' & qval < 0.0005, '***', 
                     ifelse(Category == 'extra.vs.intra' & qval < 0.005, '**', 
                            ifelse(Category == 'extra.vs.intra' & qval < 0.05, '*', ''))))
Ap2.tab.means.sd.qval$sig[Ap2.tab.means.sd.qval$cond == 'intra'] = ''

str.pos <- Ap2.tab.means.sd.qval %>% group_by(GeneName) %>% summarise(max.expr = max(mean + sd))
Ap2.tab.means.sd.qval <- left_join(Ap2.tab.means.sd.qval, str.pos, by = 'GeneName')
XX <- Ap2.tab.means.sd.qval %>% dplyr::select('GeneName', 'Passage', 'mean', 'sd', 'AP2', 'cond', 'sig', 'max.expr')

p1 <- ggplot(XX, aes(x = Passage, y = mean, group = cond, color = cond)) + 
  geom_line(size=1) + 
  geom_point(size = 1) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +  
  geom_text(aes(x = Passage, y = max.expr + 1, label=sig), color = 'black') +
  theme_bw() + 
  facet_wrap(AP2~., scales='free') + 
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1, linetype="solid"))

plot(p1)



ggsave(filename="../AP2_trends.pdf", plot=p1,
       width = 10, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)
