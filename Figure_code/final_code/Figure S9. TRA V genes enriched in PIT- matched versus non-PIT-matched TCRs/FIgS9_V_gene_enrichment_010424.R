rm(list = ls())
require(Biostrings)
require(limma)
require(plyr)
require(gdata)
library(reshape2)
library(psych)
library(ggplot2)
library(r2symbols)
library(stringdist)
library(scales)
library(Peptides)

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data") # folder containing data required to generate figures

## load TCRs

## load PIT matched TCRs

levIAR1 = read.csv("Levenshtein_index_IAR_CD4_with_islet_TCRS_lv.lt9.csv", stringsAsFactors = F)
levIAR1$index = NULL
levIAR1$set = c("IAR1")

filename = ("Levenshtein_index_P324_P474_IAR_CD4_with_islet_TCRS_lv.lt6.csv")

levIAR2 = read.csv(filename, stringsAsFactors = F)
levIAR2$index1 = NULL
levIAR2$index2 = NULL
levIAR2$set = NULL
levIAR2$set = c("IAR2")

levSubIAR1 = subset(levIAR1, levIAR1$lv <2); nrow(levSubIAR1) # 1389 (573 unique) with IAR1; 1343 (555 unique(with IAR2)
levSubIAR2 = subset(levIAR2, levIAR2$lv <2); nrow(levSubIAR2) # 1389 (573 unique) with IAR1; 1343 (555 unique(with IAR2)

levSub = rbind(levSubIAR1, levSubIAR2) # 2732

#levSub = subset(levSub, set == "IAR1") # for troubleshooting. 1389. 
###########################
## compile IMGT parameters

## compare segment sequence features for PIT matching and non-matching TRA chains

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data") # folder containing data required to generate figures

imgt1 = read.delim("P91_3_Nt-sequences.txt", stringsAsFactors = F)
imgt1$set = c('IAR1')
seqId1 = data.frame(strsplit2(imgt1$Sequence.ID, split = "_"))
colnames(seqId1) = c("libid", "junction")
imgt1 = cbind(seqId1, imgt1)

imgt2 = read.delim("P325_P474_3_Nt-sequences.txt", stringsAsFactors = F)
imgt2$set = c('IAR2')

## add junction and libid to imgt2

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data") # folder containing data required to generate figures
imgt2Tcrs = read.delim("P325_P474_CD4_unique_TCRs.txt", stringsAsFactors = F)
seqId2 = data.frame(imgt2Tcrs$libid, imgt2Tcrs$junction)
colnames(seqId2) = c("libid", "junction")
imgt2 = cbind(seqId2, imgt2)

imgt1AA =  read.delim("P91_5_AA-sequences.txt", stringsAsFactors = F)
imgt1AA$set = c('IAR1')
imgt2AA =  read.delim("P325_P474_5_AA-sequences.txt", stringsAsFactors = F)
imgt2AA$set = c('IAR2')
imgt2AA$X = NULL

imgt = rbind(imgt1, imgt2)
imgtAA = rbind(imgt1AA, imgt2AA)

## check on id order

all(imgt$Sequence.number == imgtAA$Sequence.number)
#[1] TRUE
all(imgt$Sequence.ID == imgtAA$ID)
#[1] TRUE

table(imgt$set)
table(imgtAA$set)

## both give
#IAR1 IAR2 
#3168 3217 

## add protein sequence to imgt

imgt$junction = imgtAA$JUNCTION[match(imgt$Sequence.number, imgtAA$Sequence.number)]
imgt$PITmatch = imgt$junction %in% levSub$aJunc1

imgt$chainType = ifelse(imgt$V.GENE.and.allele %in% grep("TRA", imgt$V.GENE.and.allele, value = T), "TRA", 
					ifelse(imgt$V.GENE.and.allele %in% grep("TRB", imgt$V.GENE.and.allele, value = T), "TRB", "other"))

imgt = subset(imgt, !chainType == "other") #3226 T, 3141 F

#imgt = subset(imgt, set == "IAR2") # for troubleshooting

imgt = subset(imgt, !(set == "IAR2" & chainType == "TRB" & PITmatch == "TRUE")) # these shouldn't be here. remove from analysis
      
## calculate lengths

imgt$NRegion = nchar(imgt$N.REGION) + nchar(imgt$N1.REGION)
imgt$juncLen = nchar(imgt$JUNCTION) 
imgt$threeVRegion = nchar(imgt$X3.V.REGION) 
imgt$fiveJRegion = nchar(imgt$X5.J.REGION) 
imgt$JRegion = nchar(imgt$J.REGION) 
imgt$cdr1 = nchar(imgt$CDR1.IMGT)
imgt$cdr2 = nchar(imgt$CDR2.IMGT)
imgt$cdr3 = nchar(imgt$CDR3.IMGT)
imgt$fr1 = nchar(imgt$FR1.IMGT)
imgt$fr2 = nchar(imgt$FR2.IMGT)
imgt$fr3 = nchar(imgt$FR3.IMGT)
imgt$fr4 = nchar(imgt$FR4.IMGT)

imgt$chainType = ifelse(imgt$V.GENE.and.allele %in% grep("TRA", imgt$V.GENE.and.allele, value = T), "TRA", 
					ifelse(imgt$V.GENE.and.allele %in% grep("TRB", imgt$V.GENE.and.allele, value = T), "TRB", "other"))


imgt$hydro = hydrophobicity( imgt$junction, scale = "Eisenberg")

imgtTra = subset(imgt, chainType == "TRA") # 3226
imgtTrb = subset(imgt, chainType == "TRB") # 2863

table(imgt$PITmatch, imgt$chainType)

# IAR1 only
#table(imgt$PITmatch, imgt$chainType)
       
#         TRA  TRB
#  FALSE 1033 1562
#  TRUE   573    0
  
# IAR2 only

#table(imgt$PITmatch, imgt$chainType)
       
#         TRA  TRB
#  FALSE 1330 1301
#  TRUE   290  278
  
  
# both
 
#table(imgt$PITmatch, imgt$chainType)
       
#         TRA  TRB
#  FALSE 2356 2863
#  TRUE   870  278

# after removal of erroneous TRB chains

table(imgt$PITmatch, imgt$chainType)
       
#         TRA  TRB
#  FALSE 2356 2863
#  TRUE   870    0

## identify TRB chains that pair with PITmatched TRA chains & add PIT matching status for TRA chains.

traMatchT = subset(imgtTra,PITmatch == "TRUE") # 870
traMatchF = subset(imgtTra,PITmatch == "FALSE") # 2356

trbMatchT = subset(imgtTrb,PITmatch == "TRUE") # 0
imgtTrb$PITmatch = ifelse(imgtTrb$libid %in% traMatchT$libid, "TRUE", "FALSE") # 739 TRUE, 2124 FALSE

trbMatchT = subset(imgtTrb,PITmatch == "TRUE") # 739
trbMatchF = subset(imgtTrb,PITmatch == "FALSE") # 2124

#toPlot = imgtTra

ks.test(traMatchT$NRegion, traMatchF$NRegion, simulate.p.value = T, B = 100000, alternative = "greater") #p-value = 1e-5

summary(traMatchT$NRegion) # median = 3
summary(traMatchF$NRegion) # median = 4

x = imgtTrb

colnames(x) = gsub("V.GENE.and.allele", "v_gene", colnames(x))
xnames = unique(x$v_gene)
ncx = length(xnames)

DF = data.frame(matrix(nrow= 0, ncol = 7))

for(i in 1:ncx){
q = xnames[i]
qsub = subset(x, v_gene == q)

tSub = subset(qsub, PITmatch == "TRUE")
fSub = subset(qsub, PITmatch == "FALSE")

t = nrow(traMatchT)
f = nrow(traMatchF)

forTestA = c(nrow(tSub), t - nrow(tSub))
forTestB = c(nrow(fSub), f - nrow(fSub))
c = rbind(forTestA, forTestB)

test = fisher.test(c)

result = c(q, nrow(tSub), t - nrow(tSub), nrow(fSub), f-nrow(fSub), test$p.value, test$estimate)

DF[i,] = result
}

colnames(DF) = c("v_gene", "noT", "totT", "noF", "totF", "pVal", "OR")
DF$pVal = as.numeric(DF$pVal)
DF$fdr = p.adjust(DF$pVal, method = "BH")
DF$PITmatch = x$PITmatch[match(DF$v_gene, x$v_gene)]
DF$PITmatch = gsub("TRUE", "PIT-matched", DF$PITmatch)
DF$PITmatch = gsub("FALSE", "non-PIT-matched", DF$PITmatch)

DFup = subset(DF, fdr <0.05 & OR>1)
DFdn = subset(DF, fdr <0.05 & OR<1)

DF$fdrAdj = DF$fdr
#DF$fdrAdj = gsub(1, 0.6, DF$fdrAdj)

DF$fdrAdj = ifelse(DF$fdrAdj == 1.0, 0.7, DF$fdrAdj)
DF$fdrAdj = as.numeric(DF$fdrAdj)

## intersect with long and short CDR1 V genes

v = data.frame(strsplit2(imgt$V.GENE.and.allele, split = " "))
imgt$v_gene = v$X2

imgtTra = subset(imgt, chainType == "TRA") # 1606

imgt15 = subset(imgtTra, nchar(imgtTra$CDR1.IMGT) == 15)
imgt18 = subset(imgtTra, nchar(imgtTra$CDR1.IMGT) == 18)
imgt21 = subset(imgtTra, nchar(imgtTra$CDR1.IMGT) == 21)

intersect(imgt15$v_gene, DFup$v_gene) # TRAV41*01
intersect(imgt18$v_gene, DFup$v_gene) # TRAV12-2*01
intersect(imgt21$v_gene, DFup$v_gene) # 0

intersect(imgt15$v_gene, DFdn$v_gene) # 0
intersect(imgt18$v_gene, DFdn$v_gene) # 0
intersect(imgt21$v_gene, DFdn$v_gene) # TRAV26-2*01, TRAV4*01

## plot v gene enrichment

if(dev.cur() >1) dev.off()

quartz(width=13,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 2))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(DF, aes(x = v_gene, y = -log10(fdrAdj), fill = PITmatch)) + geom_bar(position = position_dodge(preserve = "single"), stat="identity", width=1)
last_plot() + geom_hline(yintercept = -log10(0.05), lty = "dotted")
last_plot() + theme(axis.text.x = element_blank()) 
#last_plot() + theme(axis.text.x=element_text(angle=90, hjust=1))
xlab = "\nTRA V gene segment "
ylab = "Enrichment, -log10(pAdj)\n"
last_plot() + scale_y_continuous(limits = c(0, 4), breaks = c(0, 1, 2, 3, 4))
last_plot() + labs(x = xlab, y = ylab)
last_plot() + scale_fill_manual(values=cbPalette)

p = last_plot()

forLab = ggplot_build(p)

Y = 0.95*(max(max(-log10(p$data$fdr))))
X =30

label = c("TRAV41*01")
last_plot() + annotate("text", X, Y, label = label, size = 8, hjust = 0)

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename1 = paste0("FigS9A_TRA_V_gene_enrichment_PITmatched-TRCs.pdf")
filename2 = paste0("FigS9B_TRB_V_gene_enrichment_PITmatched-TRCs.pdf")
ggsave(filename2, p)


















