library(fda)
library(openxlsx)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(sme)


## Phenotypes
abs.path <- '~/work/ToxoplasmaGondii/'
pheno.file <- paste(abs.path, 'new_phenotypes.xlsx', sep = '')
phenotypes <- read.xlsx(pheno.file)
phenotypes.B2 <- phenotypes %>% dplyr::filter(str_detect(passage, 'B2'))
phenotypes.B2$passage <- gsub('B2 ', '', phenotypes.B2$passage)
phenotypes.B2$phenotype <- factor(phenotypes.B2$phenotype, levels = unique(phenotypes.B2$phenotype))
phenotypes.B2.plaq <- phenotypes.B2 %>% dplyr::filter(phenotype == 'plaquesize')
phenotypes.B2.reinv <- phenotypes.B2 %>% dplyr::filter(phenotype == 'reinvasion')
phenotypes.B2.rep <- phenotypes.B2 %>% dplyr::filter(phenotype == 'replication')
phenotypes.B2.plaq$passage <- factor(phenotypes.B2.plaq$passage, levels = unique(phenotypes.B2.plaq$passage))
phenotypes.B2.reinv$passage <- factor(phenotypes.B2.reinv$passage, levels = unique(phenotypes.B2.reinv$passage))
phenotypes.B2.rep$passage <- factor(phenotypes.B2.rep$passage, levels = unique(phenotypes.B2.rep$passage))


## Get the expression of each replicate
toxo.file <- paste(abs.path, 'toxo_table_batch_corrected_logCPM_expression_edgeR_DEGs_ap2_targs.xlsx', sep = '')
sep_rep <- read.xlsx(toxo.file)

## For correlation exclude P7 and RH and look at extra only
sep_rep.extra <- sep_rep %>% 
  dplyr::select(matches('extra.*rep')) %>% dplyr::select(-contains('RH')) %>%
  dplyr::select(-contains('P7'))
colnames(sep_rep.extra) <- gsub('.extra', '', colnames(sep_rep.extra))
sep_rep.extra$GeneName <- sep_rep$GeneName

sep_rep.extra <- sep_rep.extra %>% gather(key=PassageRep, value = Expr, -GeneName)
sep_rep.extra <- sep_rep.extra %>% 
  mutate(Passage = gsub("\\.rep.*", '', PassageRep), Rep = gsub(".*rep\\.", '', PassageRep)) %>%
  dplyr::select(-c('PassageRep'))

#sep_rep.extra$Passage <- gsub('RH', 'P300', sep_rep.extra$Passage)
Passage <- as.numeric(gsub('P', '', sep_rep.extra$Passage))
sep_rep.extra$Passage <- Passage

sep_rep.extra$Rep <- factor(sep_rep.extra$Rep, levels = unique(sort(sep_rep.extra$Rep)))
sep_rep.extra$Passage <- factor(sep_rep.extra$Passage, levels = unique(sort(sep_rep.extra$Passage)))

## Use the fda smoothing splines to fit the curves
getFit <- function(gene, sep_rep.extra, tm.grid = T){
  W.mat <- sep_rep.extra %>% dplyr::select(c('GeneName', 'Passage', 'Rep', 'Expr')) %>%
    dplyr::filter(GeneName == gene) %>%
    spread(key = Passage, value = Expr) %>%
    dplyr::select(-c('GeneName', 'Rep'))  %>% t()
  passages <- sort(as.numeric(as.character(unique(sep_rep.extra$Passage))))
  p.start <- passages[1]
  p.end <- passages[length(passages)]
  B.basis <- create.bspline.basis(rangeval = c(p.start, p.end), nbasis = length(passages))
  exprPar <- fdPar(fdobj = B.basis,lambda = 2)
  exprSmooth <- smooth.basis(passages, W.mat, exprPar)
  if(tm.grid){
    tm <- p.start:p.end
  }else{
    tm <- passages
  }
  fitted.mean <- eval.fd(evalarg = tm, fdobj = mean(exprSmooth$fd))
  return(fitted.mean)
}


sms.fits <- lapply(unique(sep_rep.extra$GeneName), function(x) getFit(x, sep_rep.extra))
sms.fits <- do.call('cbind', sms.fits)
sms.fits <- as.data.frame(sms.fits)
colnames(sms.fits) <- unique(sep_rep.extra$GeneName)

passages <- sort(as.numeric(as.character(unique(sep_rep.extra$Passage))))
p.start <- passages[1]
p.end <- passages[length(passages)]

sms.fits$Passage <- p.start:p.end

sms.fits <- sms.fits %>% gather(key = GeneName, value = Expr, -Passage)


## Fit the phenotypes
fitPheno <- function(Pheno){
  P.Mat <- Pheno %>% dplyr::select(c('passage', 'mean')) %>% 
    mutate(passage = as.numeric(gsub('P', '', passage)))
  Passages <- P.Mat$passage
  P.Mat <- P.Mat %>% dplyr::select('mean') %>% as.matrix()
  rownames(P.Mat) <- Passages
  B.basis <- create.bspline.basis(rangeval = c(Passages[1], Passages[length(Passages)]), 
                                  nbasis = length(Passages) + 1)
  phenoPar <- fdPar(fdobj = B.basis,lambda = 2)
  phenoSmooth <- smooth.basis(Passages, P.Mat, phenoPar)
  fitted.mean <- eval.fd(evalarg = Passages[1]:Passages[length(Passages)], fdobj = mean(phenoSmooth$fd))
  XX <- data.frame(Passage = c(Passages[1]:Passages[length(Passages)]), 
                   fit = fitted.mean)
  colnames(XX) <- c('Passage', 'Pheno')
  return(XX)
}

sms.Pheno.reinv <- fitPheno(phenotypes.B2.reinv)
colnames(sms.Pheno.reinv) <- c('Passage', 'reinv')
sms.Pheno.plaq <- fitPheno(phenotypes.B2.plaq)
colnames(sms.Pheno.plaq) <- c('Passage', 'plaq')
sms.Pheno.rep <- fitPheno(phenotypes.B2.rep)
colnames(sms.Pheno.rep) <- c('Passage', 'rep')

sms.fit.pheno <- sms.fits %>% inner_join(sms.Pheno.reinv, by = 'Passage') %>% 
  inner_join(sms.Pheno.plaq, by = 'Passage') %>% 
  inner_join(sms.Pheno.rep, by = 'Passage')

cor.mat <- sms.fit.pheno %>% group_by(GeneName) %>% 
  summarise(cor.p.reinv = cor(reinv, Expr, method = "pearson"),
            cor.s.reinv = cor(reinv, Expr, method = "spearman"),
            cor.p.plaq = cor(plaq, Expr, method = "pearson"),
            cor.s.plaq = cor(plaq, Expr, method = "spearman"),
            cor.p.rep = cor(rep, Expr, method = "pearson"),
            cor.s.rep = cor(rep, Expr, method = "spearman"))

expr_tab <- read.xlsx(toxo.file)

expr_tab <- left_join(expr_tab, cor.mat, by = 'GeneName')
cor.file <- paste(abs.path, 
                  'toxo_table_batch_corrected_logCPM_expression_edgeR_DEGs_ap2_targs_phenotype_correlations.xlsx', 
                  sep = '')
write.xlsx(expr_tab, cor.file)


### Make a 4 Way matirx with correlations, AP2 targets, gene sets and DEGs.
threewayfile <- paste(abs.path, 'ThreeWayDEGoverlaps.xlsx', sep = '')
ThreeWay <- read.xlsx(threewayfile)
FourWay <- left_join(ThreeWay, cor.mat, by = 'GeneName')
fourway.out <- paste(abs.path, 'FourWayDEGoverlaps.xlsx', sep = '')
write.xlsx(FourWay, fourway.out)
