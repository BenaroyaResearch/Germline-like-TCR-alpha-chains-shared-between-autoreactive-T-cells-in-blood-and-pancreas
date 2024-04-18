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
## load TCRs
## load IAR1 TCRs

IAR1Tcrs = read.csv("201512_TCR_MasterList_w_CloneIDs.csv", stringsAsFactors = F) # 5417
colnames(IAR1Tcrs) = gsub("tcrGraph_sharing_level", "sharing_level", colnames(IAR1Tcrs))

IAR1Tcrs = data.frame(IAR1Tcrs, set = "IAR1")

colNames = c("libid", "v_gene", "j_gene", "junction", "project", "donor_id", "set", "study_group")

IAR1Tcrs = IAR1Tcrs[c(colNames)]
IAR1Tcrs$study_group = gsub("early onset T1D", "newT1D", IAR1Tcrs$study_group)

## create IAR1 hla

hla = read_excel("subject_char_w_HLA.xlsx")

## correct donor id

hla$'Subject ID' = gsub('TID', "T1D", hla$'Subject ID')

## load ID key

idKey = read_excel("P91_P168 Sample_ID_key.xlsx")

idKey$tcrId = gsub("HC10, Ctrl10", "CTRL10", idKey$tcrId)

## remove whitespace

hla$'Subject ID' = stringr::str_trim(hla$'Subject ID')
idKey$suppTabId = stringr::str_trim(idKey$suppTabId)
IAR1Tcrs$donor_id = stringr::str_trim(IAR1Tcrs$donor_id)

## add suppTabId and then match HLA

IAR1Tcrs$suppTabId = idKey$suppTabId[match(IAR1Tcrs$donor_id, idKey$tcrId)]
IAR1Tcrs$hla = hla$DRB1[match(IAR1Tcrs$suppTabId, hla$"Subject ID")]

## check to see that HLA is present for most donors

temp = subset(IAR1Tcrs, is.na(IAR1Tcrs$hla)) # 0/5417. 
unique(temp$donor_id) #0

tempa = subset(hla, hla$"Subject ID" == "CTRL10")
tempb = subset(IAR1Tcrs, IAR1Tcrs$donor_id == "CTRL10")

tempa$'Subject ID'
unique(tempb$donor_id)

colNames1 = c("libid", "v_gene", "j_gene", "junction", "project", "donor_id", "set", "study_group", "hla")

IAR1Tcrs = IAR1Tcrs[c(colNames1)] # 5417

## load IAR2 TCRs

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data")

anno.filename = "P325_P474_comb_CD4+_anno_w_hla.txt"
tcrs.filename = "P325_P474_comb_CD4+_TCR_w_hla.txt"

test1 = read.delim(anno.filename, stringsAsFactors = F)
test2 = read.delim(tcrs.filename, stringsAsFactors = F)

test2$donor_id = test1$donorId[match(test2$libid, test1$libid)]
test2$study_group = test1$studyGroup2[match(test2$libid, test1$libid)]
test2$study_group = gsub("roT1D", "newT1D", test2$study_group)

test2 = subset(test2, !study_group == "estT1D")

IAR2Tcrs = data.frame(test2, set = "IAR2")
IAR2Tcrs$hla = test1$hla[match(IAR2Tcrs$donor_id, test1$donorId)]

colNames = c("libid", "v_gene", "j_gene", "junction", "project", "donor_id", "set", "study_group", "hla")

IAR2Tcrs = IAR2Tcrs[c(colNames)]

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_code")

## combine set 1 amd set2

tcrsComb = rbind(IAR1Tcrs, IAR2Tcrs) # 7727

tcrsComb$studyGroup = tcrsComb$study_group

tcrsComb$studyGroup = factor(tcrsComb$studyGroup, levels = c("HC", "AAbNeg", "1AAb", "2AAb",  "newT1D", "T1D" ))

## save tcrsComb as Supplemental Table for manuscript

toSave = tcrsComb
toSave$set = gsub("IAR1", "Cohort1", toSave$set)
toSave$set = gsub("IAR2", "Cohort2", toSave$set)
table(toSave$set )

#setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data")
#filename = "Table_S4_Compiled_and_filtered_TCR_sequences_used_in_this_study.csv"
#write.csv(toSave, filename)

## subset tcrsComb by HLA, if desired

#tcrsComb = subset(tcrsComb, tcrsComb$hla %in% grep("04", tcrsComb$hla, value = T)) # 3626 for "04"; 4101 for NOT 04

## determine expanded TCRs

no =  ddply(tcrsComb,.(junction), plyr::summarize, sum = length(libid)) # 6529 junctions, sum(no$sum) = 7727

cut <- 2
no.sub = subset(no, sum>=cut) # 474

E = subset(tcrsComb, junction %in% no.sub$junction) # 855 E junctions, 135 unique
libs = E$libid
E.cell = subset(tcrsComb, tcrsComb$libid %in% libs) # 956

E.cell$study_group = tcrsComb$study_group[match(E.cell$libid, tcrsComb$libid)] # 1954
	
frxn.e = length(unique(E.cell$libid))/length(unique(tcrsComb$libid)) # 24.4%
	
table(E.cell$studyGroup)

#HC AAbNeg   1AAb   2AAb newT1D    T1D 
# 284     23     56    101    812    678 # all tcrs, not subsetted
#   267      8     21     31    119    510 # hla  "04" subset
   
## add expanded cells to tcrsComb

tcrsComb$E = tcrsComb$libid %in% E.cell$libid # 1954 TRUE 5773 FALSE for junctions; 
tcrsComb$expanded = tcrsComb$E 
tcrsComb$expanded = gsub("TRUE", "E", tcrsComb$expanded)
tcrsComb$expanded = gsub("FALSE", "NE", tcrsComb$expanded)

tcrsComb$chainType = ifelse(tcrsComb$v_gene %in% grep("TRA", tcrsComb$v_gene, value = T), "TRA",
						ifelse(tcrsComb$v_gene %in% grep("TRB", tcrsComb$v_gene, value = T), "TRB","other"))

## remove iNKT and MAIT cell sequencces

iNkt1 = subset(tcrsComb, junction == "CVVSDRGSTLGRLYF")
iNkt2 = subset(tcrsComb, libid %in% iNkt1$libid)
mait1 = subset(tcrsComb, v_gene == "TRAV1-2" & (j_gene == "TRAJ33" | j_gene == "TRAJ20" | j_gene == "TRAJ12")) # 18
mait2 = subset(tcrsComb, libid %in% mait1$libid) # 45

tcrsCombSub3 = subset(tcrsComb, !junction %in% iNkt2$junction) # 
tcrsCombSub3 = subset(tcrsCombSub3, !junction %in% mait2$junction) # 7633

## modify combined TCRs

tcrsCombSub3$IARmatch = tcrsCombSub3$junction %in% levSubIAR1$aJunc1 | tcrsCombSub3$junction %in% levSubIAR2$aJunc1 # 1235 T, 6398 F
tcrsCombSub3$expanded = factor(tcrsCombSub3$expanded, levels = c("E", "NE")) # 

tcrsCombSub3Tra = subset(tcrsCombSub3, tcrsCombSub3$chainType == "TRA") # 1814
tcrsCombSub3Tra.u = tcrsCombSub3Tra[!duplicated(tcrsCombSub3Tra$junction),] # 1512

tcrsCombSub3.u = tcrsCombSub3[!duplicated(tcrsCombSub3$junction),] # 2967
tcrsTra = subset(tcrsCombSub3.u, chainType == "TRA") # 1512
tcrsTrb = subset(tcrsCombSub3.u, chainType == "TRB") # 1433

###########################
## compile IMGT parameters

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data") # folder containing data required to generate figures

imgtAA = read.delim("Combined_P91_P325_P474_AB_TCRs_5_AA-sequences.txt")
imgtNt = read.delim("Combined_P91_P325_P474_AB_TCRs_3_Nt-sequences.txt")

all(imgtAA$Sequence.ID == imgtNt$Sequence.ID)
#[1] TRUE

## add protein sequence to imgt

imgt = imgtNt

libids = data.frame(strsplit2(imgt$Sequence.ID, split = "_"))
colnames(libids) = c("libid", "junction")
all(imgt$junction == libids$junction)
#[1] TRUE # OK to combine

imgt$libid = libids$libid
imgt$junction = imgtAA$JUNCTION[match(imgt$Sequence.ID, imgtAA$Sequence.ID)]
imgt$PITmatch = imgt$junction %in% levSubIAR1$aJunc1 | imgt$junction %in% levSubIAR2$aJunc1 # 

imgt$chainType = ifelse(imgt$V.GENE.and.allele %in% grep("TRA", imgt$V.GENE.and.allele, value = T), "TRA", 
					ifelse(imgt$V.GENE.and.allele %in% grep("TRB", imgt$V.GENE.and.allele, value = T), "TRB", "other"))

table(imgt$PITmatch, imgt$chainType)

imgtTra = subset(imgt, chainType == "TRA")
aMatchT = subset(imgtTra, PITmatch == "TRUE") # 942
imgtTrb = subset(imgt, chainType == "TRB")
imgtTrb$PITmatch = imgtTrb$libid %in% aMatchT$libid

imgt = rbind(imgtTra, imgtTrb)
      
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

imgt$hydro = hydrophobicity( imgt$junction, scale = "Eisenberg")

imgtTra = subset(imgt, chainType == "TRA") # 3226
imgtTrb = subset(imgt, chainType == "TRB") # 2863

table(imgt$PITmatch, imgt$chainType)

#       TRA  TRB
# FALSE 2322 2355
# TRUE   942  832 
  
## identify TRB chains that pair with PITmatched TRA chains & add PIT matching status for TRA chains.

traMatchT = subset(imgtTra,PITmatch == "TRUE") # 942
traMatchF = subset(imgtTra,PITmatch == "FALSE") # 2322

trbMatchT = subset(imgtTrb,PITmatch == "TRUE") # 832
trbMatchF = subset(imgtTrb,PITmatch == "FALSE") # 2355

## summarize junction length and hydrophobicity

summary(traMatchT$juncLen) # 39 median 
summary(traMatchF$juncLen) # 42 median

summary(traMatchT$hydro) # 0.2321
summary(traMatchF$hydro) # 0.17786

summary(trbMatchT$juncLen) # 42
summary(trbMatchF$juncLen) # 42

summary(trbMatchT$hydro) # 0.08728 
summary(trbMatchF$hydro) # 0.09077

## KS tests
test1 = ks.test(traMatchT$juncLen, traMatchF$juncLen, alternative = "two.sided"); test1 #p-value <2.2e-16 for junction TRA; 
test2 = ks.test(traMatchT$hydro, traMatchF$hydro, alternative = "two.sided"); test2 #p-value = 4.745e-11 for hydro TRA;
test3 = ks.test(trbMatchT$juncLen, trbMatchF$juncLen, alternative = "two.sided"); test3 #p-value = 0.9668 for junction TRB; 
test4 = ks.test(trbMatchT$hydro, trbMatchF$hydro, alternative = "two.sided"); test4 #p-value = 0.3768 for hydro TRB; 

p.adjust(c(2.2e-16,4.745e-11,0.9668, 0.3768), method = "BH")
#[1] 8.800e-16 9.490e-11 9.668e-01 5.024e-01
# 
## use **** adjusted p-values for both TRA juncLen and hydro; NS for both TRB juncLen and hydro

###################### 
## for TRA junctions

toPlotName = c("imgtTra")
toPlot = get(toPlotName) # select TRA or TRB chains
toPlot$PITmatch = gsub("TRUE", "PIT-matched", toPlot$PITmatch)
toPlot$PITmatch = gsub("FALSE", "non-PIT-matched", toPlot$PITmatch)


## Junction length

if(dev.cur()>1) dev.off()
quartz(width=13,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 2))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(toPlot, aes(x = juncLen)) + geom_density(aes(fill = PITmatch), alpha = 0.5, adjust = 1)
#last_plot() + theme(legend.spacing.y = unit(0.5, 'cm')) + guides(fill = guide_legend(byrow = TRUE))

xlab = "\nJunc length (nt)"
ylab = "Density\n"
last_plot() + labs(x = xlab, y = ylab)
last_plot() + geom_vline(xintercept = median(traMatchT$juncLen), lty = "solid")
last_plot() + geom_vline(xintercept = median(traMatchF$juncLen), lty = "dashed")
last_plot() + scale_fill_manual(values=cbPalette)

p = last_plot()

forLab = ggplot_build(p)

Y = 0.9*(max(forLab$data[[1]]$density))
X = 1*(min((forLab$data[[1]]$x)))

label = paste("KS test, \np-value =", scientific(test1$p.value, 2), sep = " ")
#last_plot() + annotate("text", X, Y, label = label, size = 8, hjust = 0)

p = last_plot()
setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("Fig2A_", toPlotName, "_junction_length_by_PITmatch.pdf") #TRA

ggsave(filename, p)

## Hydrophobicity

if(dev.cur()>1) dev.off()
quartz(width=13,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 2))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(toPlot, aes(x = hydro)) + geom_density(aes(fill = PITmatch), alpha = 0.5, adjust = 1)

xlab = "\nHydrophobicity"
ylab = "Density\n"
last_plot() + labs(x = xlab, y = ylab)
last_plot() + geom_vline(xintercept = median(traMatchT$hydro), lty = "solid")
last_plot() + geom_vline(xintercept = median(traMatchF$hydro), lty = "dashed")
last_plot() + scale_fill_manual(values=cbPalette)

p = last_plot()

forLab = ggplot_build(p)

Y = 0.9*(max(forLab$data[[1]]$density))
X = 1*(min((forLab$data[[1]]$x)))

label = paste("KS test, \np-value =", scientific(test2$p.value, 2), sep = " ")
#last_plot() + annotate("text", X, Y, label = label, size = 8, hjust = 0)

p = last_plot()

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("Fig2C_", toPlotName, "_hydrophobicity_by_PITmatch.pdf") #TRA

ggsave(filename, p)

 ######################
 ## plots for TRB junctions

 toPlotName = c("imgtTrb")
 toPlot = get(toPlotName) # select TRA or TRB chains
 toPlot$PITmatch = gsub("TRUE", "PIT-matched", toPlot$PITmatch)
 toPlot$PITmatch = gsub("FALSE", "non-PIT-matched", toPlot$PITmatch)

 ## Junction length

 if(dev.cur()>1) dev.off()
 quartz(width=13,height=8, dpi=72)  ### open plotting window
 update_geom_defaults("line", aes(size = 2))
 update_geom_defaults("point", aes(size = 2))

 theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

 cbPalette = c('#66c2a5','#fc8d62', 'gray')

 ggplot(toPlot, aes(x = juncLen)) + geom_density(aes(fill = PITmatch), alpha = 0.5, adjust = 1)
 #last_plot() + theme(legend.spacing.y = unit(0.5, 'cm')) + guides(fill = guide_legend(byrow = TRUE))

 xlab = "\nJunc length (nt)"
 ylab = "Density\n"
 last_plot() + labs(x = xlab, y = ylab)
 last_plot() + geom_vline(xintercept = median(trbMatchT$juncLen), lty = "solid")
 last_plot() + geom_vline(xintercept = median(trbMatchF$juncLen), lty = "dashed")
 last_plot() + scale_fill_manual(values=cbPalette)

 p = last_plot()

 forLab = ggplot_build(p)

 Y = 0.9*(max(forLab$data[[1]]$density))
 X = 1*(min((forLab$data[[1]]$x)))

 label = paste("KS test, \np-value =", scientific(test1$p.value, 2), sep = " ")
 #last_plot() + annotate("text", X, Y, label = label, size = 8, hjust = 0)

 p = last_plot()
 setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
 filename = paste0("Fig2B_", toPlotName, "_junction_length_by_PITmatch.pdf") #TRA

 ggsave(filename, p)

 ## Hydrophobicity

 if(dev.cur()>1) dev.off()
 quartz(width=13,height=8, dpi=72)  ### open plotting window
 update_geom_defaults("line", aes(size = 2))
 update_geom_defaults("point", aes(size = 2))

 theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

 cbPalette = c('#66c2a5','#fc8d62', 'gray')

 ggplot(toPlot, aes(x = hydro)) + geom_density(aes(fill = PITmatch), alpha = 0.5, adjust = 1)

 xlab = "\nHydrophobicity"
 ylab = "Density\n"
 last_plot() + labs(x = xlab, y = ylab)
 last_plot() + geom_vline(xintercept = median(trbMatchT$hydro), lty = "solid")
 last_plot() + geom_vline(xintercept = median(trbMatchF$hydro), lty = "dashed")
 last_plot() + scale_fill_manual(values=cbPalette)

 p = last_plot()

 forLab = ggplot_build(p)

 Y = 0.9*(max(forLab$data[[1]]$density))
 X = 1*(min((forLab$data[[1]]$x)))

 label = paste("KS test, \np-value =", scientific(test2$p.value, 2), sep = " ")
 #last_plot() + annotate("text", X, Y, label = label, size = 8, hjust = 0)

 p = last_plot()

 setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
 filename = paste0("Fig2D_", toPlotName, "_hydrophobicity_by_PITmatch.pdf") #TRA

 ggsave(filename, p)






