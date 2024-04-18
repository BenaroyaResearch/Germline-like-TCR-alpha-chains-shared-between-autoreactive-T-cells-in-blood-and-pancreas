rm(list = ls())
require(Biostrings)
require(limma)
require(plyr)
require(gdata)
require(ggplot2)
library(stringdist)
library(ggsignif)
library(Peptides)
#library(brigg)
library(tcrGraph)
library(tidyverse)
library(igraph)
library(graphlayouts)
library(ggraph)
library(viridis)
library(gridExtra)
library(plyr)
library(scales)
library(Peptides)

#################################
## load vdjdb database - test different epitope per junction

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data/")

db <- read.delim("vdjdb.slim.txt", stringsAsFactors=FALSE) # 

colnames(db)[2] = c("junction")

hu = subset(db, species == "HomoSapiens") #4085; 0900417 version 8827; 62753 07312020 version; 67871 11242022 version
huNOT0 = subset(hu, !complex.id == 0)
colnames(huNOT0) = gsub("complex.id", "libid", colnames(huNOT0))
colnames(huNOT0) = gsub("v.segm", "v_gene", colnames(huNOT0))
colnames(huNOT0) = gsub("j.segm", "j_gene", colnames(huNOT0))
huNOT0$reference.id = gsub("https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#", "https://pages.10xgenomics.com/rs/446-PBO-704/images/10x_AN047_IP_A_New_Way_of_Exploring_Immunity_Digital.pdf", huNOT0$reference.id)
huNOT0$juncLen = nchar(huNOT0$junction)

noChains =  ddply(huNOT0,.(libid), plyr::summarize, sum = length(junction)) # 
pairs = subset(noChains, sum == 2)
#huNOT0 = subset(huNOT0, libid %in% pairs$libid)

hu10X = subset(huNOT0, huNOT0$reference.id %in% grep("10x", huNOT0$reference.id, value = T))
#hu10X = hu10X[c("libid", "Gene", "v_gene", "junction", "j_gene", "Epitope.gene")]
#colnames(hu10X) = c("libid", "chainType", "v_gene", "junction", "j_gene", "study_group")

###############################
## hydrophobicity of VDJdb single and multiple epitope TCR chains

huTra = subset(huNOT0, gene == "TRA")
randTrav= subset(huTra, huTra$juncLen == sample(huTra$juncLen, size = length(huTra$juncLen)))
huTrb = subset(huNOT0, gene == "TRB")

## select chain for statistics

vdjdbVar = c("huTra")
vdjdb = get(vdjdbVar)

vdjdb$hydro = hydrophobicity( vdjdb$junction, scale = "Eisenberg")

noEpiPerJunc =  ddply(vdjdb,.(gene, junction, juncLen, hydro), plyr::summarize, sum = length(antigen.epitope)) # 20786

one = subset(noEpiPerJunc, sum ==1)
gtOne = subset(noEpiPerJunc, sum >1)
one$NoEpitopes = c("one")
gtOne$NoEpitopes = c("multiple")

summary(one$juncLen)
summary(gtOne$juncLen)

## tests
## for TRA
ks.test(one$juncLen, gtOne$juncLen) #p-value = <2.2e-16 for TRA,****
ks.test(one$hydro, gtOne$hydro) #p-value = 0.004015 for TRA, ***

## for TRB
ks.test(one$juncLen, gtOne$juncLen) #p-value = 0.2212 for TRB
ks.test(one$hydro, gtOne$hydro) #p-value = 0.8615 for TRB

summary(one$hydro) # median = 0.19571 for TRA, 0.114 for TRB
summary(gtOne$hydro) # median = 0.21168 fro TRA, 0.115 for TRB

###################
## plots for TRA

## No. specificities density

## select chain

vdjdbVar = c("huTra")
vdjdb = get(vdjdbVar)

toPlot2 = vdjdb

## aggregate by numbers of antigen specificities

noEpiPerJunc =  ddply(vdjdb,.(gene, junction, vdjdb.score, juncLen), plyr::summarize, sum = length(antigen.gene)) # 19435

one = subset(noEpiPerJunc, sum ==1)
gtOne = subset(noEpiPerJunc, sum >1)
one$NoEpitopes = c("one")
gtOne$NoEpitopes = c("multiple")
one$hydro = hydrophobicity( one$junction, scale = "Eisenberg")
gtOne$hydro = hydrophobicity( gtOne$junction, scale = "Eisenberg")

toPlot = rbind(one, gtOne)
toPlot$NoEpitopes = factor(toPlot$NoEpitopes, levels = c("one", "multiple"))
toPlot$group = ifelse(toPlot$sum == 1, "oneSpecificity","multiSpecificity")
toPlot$group = factor(toPlot$group, levels = c("oneSpecificity", "multiSpecificity"))
#toPlot$hydro = hydrophobicity( toPlot$junction, scale = "Eisenberg")

if(dev.cur()>1)dev.off()

quartz(width=14,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 2))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(toPlot, aes(x = sum)) + geom_density(aes(fill = NoEpitopes), alpha = 0.5, adjust = 1)
last_plot() + scale_x_continuous(limits = c(0,6), breaks = c(0, 2,4, 6)) 
last_plot() + scale_fill_manual(values=cbPalette)

xlab = "\nNumber VDJdb specificities"
ylab = "TCR chain density\n"
last_plot() + labs(x = xlab, y = ylab)

p = last_plot()
setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("Fig5A_", vdjdbVar, "_chain_density_number_of_specificities_in_VDJdb.pdf")
ggsave(filename, p)

################ 
## Shubham's Stitchr-derived IMGT parameters

## CDR3 length density

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data") # folder containing data required to generate figures

vdjMeta = read.csv("petervdjTRATRBmeta.csv")
hu2 = subset(vdjMeta, Species == "HomoSapiens") #60621
huNOT02 = subset(hu2, !complex.id == 0)
colnames(huNOT02) = gsub("complex.id", "libid", colnames(huNOT02))
colnames(huNOT02) = gsub("v.segm", "v_gene", colnames(huNOT02))
colnames(huNOT02) = gsub("j.segm", "j_gene", colnames(huNOT02))
huNOT02$Reference = gsub("https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#", "https://pages.10xgenomics.com/rs/446-PBO-704/images/10x_AN047_IP_A_New_Way_of_Exploring_Immunity_Digital.pdf", huNOT02$Reference)
huNOT02$juncLen = nchar(huNOT02$CDR3)

imgt2 = read.delim("VDJdb_3_Nt-sequences.txt")
imgt2 = imgt2 %>% mutate_all(str_replace_all, "Homsap", "") 

## add junction protein sequence to IMGT nucleotide data

table(imgt2$Sequence.ID %in% hu2$TCR_name)
#FALSE  TRUE 
#   91 60530

imgt2$junction = hu2$CDR3[match(imgt2$Sequence.ID, hu2$TCR_name)]

## calculate lengths of sequence regions

imgt2$CDR3 = nchar(imgt2$CDR3.IMGT) 
imgt2$Nregion = nchar(imgt2$N.REGION) + nchar(imgt2$N1.REGION)

## select chain and prepare toPlot

vdjdbVar = c("huTra")
vdjdb = get(vdjdbVar)

toPlot2 = vdjdb

## aggregate by numbers of antigen specificities

noEpiPerJunc =  ddply(vdjdb,.(gene, junction, vdjdb.score, juncLen), plyr::summarize, sum = length(antigen.gene)) # 19435

one = subset(noEpiPerJunc, sum ==1)
gtOne = subset(noEpiPerJunc, sum >1)
one$NoEpitopes = c("one")
gtOne$NoEpitopes = c("multiple")
one$hydro = hydrophobicity( one$junction, scale = "Eisenberg")
gtOne$hydro = hydrophobicity( gtOne$junction, scale = "Eisenberg")

toPlot = rbind(one, gtOne)
toPlot$NoEpitopes = factor(toPlot$NoEpitopes, levels = c("one", "multiple"))
toPlot$group = ifelse(toPlot$sum == 1, "oneSpecificity","multiSpecificity")
toPlot$group = factor(toPlot$group, levels = c("oneSpecificity", "multiSpecificity"))
#toPlot$hydro = hydrophobicity( toPlot$junction, scale = "Eisenberg")

## add sequence region lengths to toPlot

toPlot$CDR3.IMGT = imgt2$CDR3.IMGT[match(toPlot$junction, imgt2$junction)]
toPlot$N.REGION = imgt2$N.REGION[match(toPlot$junction, imgt2$junction)]
toPlot$N1.REGION = imgt2$N1.REGION[match(toPlot$junction, imgt2$junction)]

toPlot$CDR3 = imgt2$CDR3[match(toPlot$junction, imgt2$junction)]
toPlot$Nregion = imgt2$Nregion[match(toPlot$junction, imgt2$junction)]

toPlotSub1 = subset(toPlot, NoEpitopes == "one")
toPlotSub2 = subset(toPlot, NoEpitopes == "multiple")

summary(toPlotSub1$CDR3)
summary(toPlotSub2$CDR3)

table(toPlotSub1$NoEpitopes, toPlotSub1$chainType)

if(dev.cur()>1)dev.off()

quartz(width=12,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 2))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(toPlot, aes(x = CDR3 )) + geom_density(aes(fill = NoEpitopes), alpha = 0.5, adjust = 0.5) # y= after_stat(scaled)
last_plot() + scale_x_continuous(limits = c(20,60), breaks = c(20, 30,40, 50, 60)) 
last_plot() + geom_vline(xintercept = median(toPlotSub1$CDR3, na.rm = T), lty = "dashed")
last_plot() + geom_vline(xintercept = median(toPlotSub2$CDR3, na.rm = T), lty = "solid")
last_plot() + scale_fill_manual(values=cbPalette)

xlab = "\nCDR3 length (nt)"
ylab = "Density\n"
last_plot() + labs(x = xlab, y = ylab)

p = last_plot()
setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("Fig5B_", vdjdbVar, "_CDR3_length_by_number_of_specificities_in_VDJdb.pdf")
ggsave(filename, p)

################ 
## Shubham's Stitchr-derived IMGT parameters

## N region length density

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data") # folder containing data required to generate figures

vdjMeta = read.csv("petervdjTRATRBmeta.csv")
hu2 = subset(vdjMeta, Species == "HomoSapiens") #60621
huNOT02 = subset(hu2, !complex.id == 0)
colnames(huNOT02) = gsub("complex.id", "libid", colnames(huNOT02))
colnames(huNOT02) = gsub("v.segm", "v_gene", colnames(huNOT02))
colnames(huNOT02) = gsub("j.segm", "j_gene", colnames(huNOT02))
huNOT02$Reference = gsub("https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#", "https://pages.10xgenomics.com/rs/446-PBO-704/images/10x_AN047_IP_A_New_Way_of_Exploring_Immunity_Digital.pdf", huNOT02$Reference)
huNOT02$juncLen = nchar(huNOT02$CDR3)

imgt2 = read.delim("VDJdb_3_Nt-sequences.txt")
imgt2 = imgt2 %>% mutate_all(str_replace_all, "Homsap", "") # 60621

## add junction protein sequence to IMGT nucleotide data

table(imgt2$Sequence.ID %in% hu2$TCR_name)
#FALSE  TRUE 
#   91 60530

imgt2$junction = hu2$CDR3[match(imgt2$Sequence.ID, hu2$TCR_name)] #44870 unique, 60621 total
imgt2 = imgt2[!duplicated(imgt2$junction), ] # 44870

## calculate lengths of sequence regions

imgt2$CDR3 = nchar(imgt2$CDR3.IMGT) 
imgt2$Nregion = nchar(imgt2$N.REGION) + nchar(imgt2$N1.REGION)

## select chain and prepare toPlot

vdjdbVar = c("huTra")
vdjdb = get(vdjdbVar) # 22100

toPlot2 = vdjdb

## aggregate by numbers of antigen specificities

noEpiPerJunc =  ddply(vdjdb,.(gene, junction, vdjdb.score, juncLen), plyr::summarize, sum = length(antigen.gene)) # 19435

one = subset(noEpiPerJunc, sum ==1)
gtOne = subset(noEpiPerJunc, sum >1)
one$NoEpitopes = c("one")
gtOne$NoEpitopes = c("multiple")
one$hydro = hydrophobicity( one$junction, scale = "Eisenberg")
gtOne$hydro = hydrophobicity( gtOne$junction, scale = "Eisenberg")

toPlot = rbind(one, gtOne)
toPlot$NoEpitopes = factor(toPlot$NoEpitopes, levels = c("one", "multiple"))
toPlot$group = ifelse(toPlot$sum == 1, "oneSpecificity","multiSpecificity")
toPlot$group = factor(toPlot$group, levels = c("oneSpecificity", "multiSpecificity"))
#toPlot$hydro = hydrophobicity( toPlot$junction, scale = "Eisenberg")

## add sequence region lengths to toPlot

toPlot$CDR3.IMGT = imgt2$CDR3.IMGT[match(toPlot$junction, imgt2$junction)]
toPlot$N.REGION = imgt2$N.REGION[match(toPlot$junction, imgt2$junction)]
toPlot$N1.REGION = imgt2$N1.REGION[match(toPlot$junction, imgt2$junction)]

toPlot$CDR3 = imgt2$CDR3[match(toPlot$junction, imgt2$junction)]
toPlot$Nregion = imgt2$Nregion[match(toPlot$junction, imgt2$junction)]

toPlotSub1 = subset(toPlot, NoEpitopes == "one")
toPlotSub2 = subset(toPlot, NoEpitopes == "multiple")

summary(toPlotSub1$CDR3) # median = 36
summary(toPlotSub2$CDR3) # median = 33

summary(toPlotSub1$Nregion) # median = 5
summary(toPlotSub2$Nregion) # median = 4

ks.test(toPlotSub1$CDR3, toPlotSub2$CDR3) # p-value <=2.2e-16, or ****
ks.test(toPlotSub1$Nregion, toPlotSub2$Nregion) # p-value <=2.2e-16, or ****
ks.test(toPlotSub1$hydro, toPlotSub2$hydro) # p-value 0.01045, or *
p.adjust(c(2.2e-16,2.2e-16,0.1045), method = "BH")
# 3.300e-16 3.300e-16 1.045e-01

if(dev.cur()>1)dev.off()

quartz(width=12,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 2))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(toPlot, aes(x = Nregion, y= after_stat(scaled))) + geom_density(aes(fill = NoEpitopes), alpha = 0.5, adjust = 0.5) # 
last_plot() + scale_x_continuous(limits = c(0,20)) #, breaks = c(0, 5, 10, )) 
last_plot() + geom_vline(xintercept = median(toPlotSub1$Nregion, na.rm = T), lty = "dashed")
last_plot() + geom_vline(xintercept = median(toPlotSub2$Nregion, na.rm = T), lty = "solid")
last_plot() + scale_fill_manual(values=cbPalette)

xlab = "\nNregion length (nt)"
ylab = "Density\n"
last_plot() + labs(x = xlab, y = ylab)

p = last_plot()
setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("Fig5C_", vdjdbVar, "_Nregion_length_by_number_of_specificities_in_VDJdb.pdf")
ggsave(filename, p)

## hydrophobicity density plot

if(dev.cur()>1)dev.off()

quartz(width=11,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 2))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(toPlot, aes(x = hydro)) + geom_density(aes(fill = NoEpitopes), alpha = 0.5, adjust = 1)

xlab = "\nHydrophobicity"
ylab = "Density\n"
last_plot() + labs(x = xlab, y = ylab)
last_plot() + geom_vline(xintercept = median(one$hydro), lty = "dashed")
last_plot() + geom_vline(xintercept = median(gtOne$hydro), lty = "solid")
last_plot() + scale_fill_manual(values=cbPalette)
last_plot() + scale_x_continuous(limits = c(-0.6, 0.6))

p = last_plot()

forLab = ggplot_build(p)

test = ks.test(one$hydro, gtOne$hydro) #p-value = <2.2e-16 for TRA,****

Y = 0.9*(max(forLab$data[[1]]$density))
X = -0.9*(max((forLab$data[[1]]$x)))

label = paste("KS test, \np-value =", scientific(test$p.value, 2), sep = " ")
#last_plot() + annotate("text", X, Y, label = label, size = 8, hjust = 0)

p = last_plot()
setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("Fig5D_", vdjdbVar, "_hydrophobicity_density_by Nos_specificities_in_VDJdb.pdf")
ggsave(filename, p)

########################
## numbers of junctions

length((toPlot2$junction))
#[1] 22100
length((one$junction))
#[1] 17921
length((gtOne$junction))
#[1] 1664

length(unique(toPlot2$junction))
#[1] 19435
length(unique(one$junction))
#[1] 17826
length(unique(gtOne$junction))
#[1] 1664



