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
library(readxl)

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data") # folder containing data required to generate figures

## load matches

filename1 = ("Levenshtein_index_IAR_CD4_with_islet_TCRS_lv.lt9.csv")
filename2 = ("Levenshtein_index_P324_P474_IAR_CD4_with_islet_TCRS_lv.lt6.csv")

levIAR1 = read.csv(filename1, stringsAsFactors = F)
levIAR2 = read.csv(filename2, stringsAsFactors = F)

levIAR1$index = NULL

levIAR2$index1 = NULL
levIAR2$index2 = NULL
levIAR2$set = NULL

levIAR1$set = c("IAR1")
levIAR2$set = c("IAR2")

levSubIAR1 = subset(levIAR1, levIAR1$lv <2); nrow(levSubIAR1) # 1389 (573 unique) with IAR1; 1343 (555 unique(with IAR2)
levSubIAR2 = subset(levIAR2, levIAR2$lv <2); nrow(levSubIAR2) # 1389 (573 unique) with IAR1; 1343 (555 unique(with IAR2)

## load TCRs

## load IAR1 TCRs

IAR1Tcrs = read.csv("201512_TCR_MasterList_w_CloneIDs.csv", stringsAsFactors = F) # 5729
colnames(IAR1Tcrs) = gsub("tcrGraph_sharing_level", "sharing_level", colnames(IAR1Tcrs))

IAR1Tcrs = data.frame(IAR1Tcrs, set = "IAR1")

colNames = c("libid", "v_gene", "j_gene", "junction", "project", "donor_id", "set", "study_group", "full_nt_sequence")

IAR1Tcrs = IAR1Tcrs[c(colNames)]
IAR1Tcrs$study_group = gsub("early onset T1D", "newT1D", IAR1Tcrs$study_group)

## add HLA
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

## check to see that HLA is present for most donos

temp = subset(IAR1Tcrs, is.na(IAR1Tcrs$hla)) # 0/5417. 
unique(temp$donor_id) #0

tempa = subset(hla, hla$"Subject ID" == "CTRL10")
tempb = subset(IAR1Tcrs, IAR1Tcrs$donor_id == "CTRL10")

tempa$'Subject ID'
unique(tempb$donor_id)

colNames1 = c("libid", "v_gene", "j_gene", "junction", "project", "donor_id", "set", "study_group", "hla", "full_nt_sequence")

IAR1Tcrs = IAR1Tcrs[c(colNames1)]

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

colNames = c("libid", "v_gene", "j_gene", "junction", "project", "donor_id", "set", "study_group", "hla", "full_nt_sequence")

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

tcrsCombSub3$PITmatch = tcrsCombSub3$junction %in% levSubIAR1$aJunc1 | tcrsCombSub3$junction %in% levSubIAR2$aJunc1 # 1235 T, 6398 F
tcrsCombSub3$expanded = factor(tcrsCombSub3$expanded, levels = c("E", "NE")) # 
tcrsCombSub3 = subset(tcrsCombSub3, !chainType == "other") # 7601

tcrsCombSub3Tra = subset(tcrsCombSub3, tcrsCombSub3$chainType == "TRA") # 3907
tcrsCombSub3Tra.u = tcrsCombSub3Tra[!duplicated(tcrsCombSub3Tra$junction),] # 3264

tcrsCombSub3.u = tcrsCombSub3[!duplicated(tcrsCombSub3$junction),] # 6451
tcrsTra = subset(tcrsCombSub3.u, chainType == "TRA") # 3264
tcrsTrb = subset(tcrsCombSub3.u, chainType == "TRB") # 3187

tcrsCombSub3Tra$studyGroup = gsub("AAbNeg", "HC", tcrsCombSub3Tra$studyGroup)

lev = rbind(levIAR1, levIAR2)

table(lev$set)

#   IAR1    IAR2 
#3002403  518011  

levSub = subset(lev, lv <2); nrow(levSub) # 2732
levSub = subset(levSub, aJunc1 %in% tcrsCombSub3$junction) # 2345, 942 unique

toCalc = subset(levSub, lv ==1) # 2258, 927 unique
toCalc = toCalc[!duplicated(toCalc$aJunc1, fromLast = T),]
toCalc$len1 = nchar(toCalc$aJunc1)
toCalc$len2 = nchar(toCalc$aJunc2)

toPlotName = c("Levenshtein_1")

## identify region of mismatch between single mismatched IAR TRA and PIT chains.. Use randomized subset for comparison
res1 <- mapply(function(x,y) first(which(bitwXor(utf8ToInt(x),utf8ToInt(y))>0)), toCalc$aJunc1, toCalc$aJunc2, USE.NAMES = FALSE) # 
res2 <- mapply(function(x,y) first(which(bitwXor(utf8ToInt(x),utf8ToInt(y))>0)), sample(lev$aJunc1, size = length(toCalc$aJunc1)), toCalc$aJunc2, USE.NAMES = FALSE) # random 

toCalc$res1 = as.numeric(res1)
toCalc$res2 = as.numeric(res2)

toCalc$res1Norm = toCalc$res1/toCalc$len1
toCalc$res2Norm = toCalc$res2/toCalc$len1

###########################
## Compare IMGT parameters

toPlot = data.frame(aJunc1 = toCalc$aJunc1, aJunc2 = toCalc$aJunc2, observed = toCalc$res1, random = toCalc$res2)
mToPlot = melt(toPlot)
#mToPlot = subset(mToPlot, variable == "random")
colnames(mToPlot) = gsub("variable", "Mismatch", colnames(mToPlot))

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

imgtNt = rbind(imgt1, imgt2)
imgtAA = rbind(imgt1AA, imgt2AA)

imgt = imgtNt

## check on id order

all(imgt$Sequence.number == imgtAA$Sequence.number)
#[1] TRUE
all(imgt$Sequence.ID == imgtAA$ID)
#[1] TRUE

imgt$lenN = nchar(imgt$N.REGION) + nchar(imgt$N1.REGION)
imgt$lenJunc = nchar(imgt$JUNCTION) 
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

#seqId = data.frame(strsplit2(imgt$Sequence.ID, split = "_"))
#colnames(seqId) = c("libid", "junction")
#imgt$junction = seqId$junction

toPlot$vEnd = imgt$V.REGION.end[match(toPlot$aJunc1, imgt$junction)]
toPlot$jStart = imgt$J.REGION.start[match(toPlot$aJunc1, imgt$junction)]
toPlot$juncStart = imgt$JUNCTION.start[match(toPlot$aJunc1, imgt$junction)]
toPlot$juncEnd = imgt$JUNCTION.end[match(toPlot$aJunc1, imgt$junction)]
toPlot$juncStartAA = ((toPlot$juncStart-toPlot$juncStart)/3 + 1)
toPlot$vEndAA = floor((toPlot$vEnd-toPlot$juncStart)/3 + 1)
toPlot$jStartAA = ceiling((toPlot$jStart-toPlot$juncStart)/3 + 1)

toPlotSub = data.frame(aJunc1 = toPlot$aJunc1, aJunc2 = toPlot$aJunc2, mismatch = toPlot$observed, vEnd = toPlot$vEndAA, jStart = toPlot$jStartAA)
toPlotSub$gte = toPlot$observed >= toPlot$vEndAA & toPlot$observed < toPlot$jStartAA # 74%
toPlotSub$gt = toPlot$observed > toPlot$vEndAA & toPlot$observed < toPlot$jStartAA # 41%
toPlotSub$eq = toPlot$observed == toPlot$vEndAA & toPlot$observed < toPlot$jStartAA # 41%

## calculate mismatches betwen V region and J region
## gte >=
table(toPlotSub$gte)

#FALSE  TRUE 
#  247   680 
# 680/(680+247) # 927 total
#[1] 0.7335491 

## gt >

table(toPlotSub$gt)

#FALSE  TRUE 
#  546   381 
# 381/(546+381) # 927 total
#[1] 0.4110032

## eq =

table(toPlotSub$gte)

#FALSE  TRUE 
#  628   299 
# 291/(628+299) # 927 total
#[1] 0.3139159

mToPlotSub = melt(toPlotSub)
mToPlotSub$variable = gsub("observed", "mismatch", mToPlotSub$variable)
mToPlotSub$variable = factor(mToPlotSub$variable, levels = c("mismatch", "vEnd", "jStart"))

#######################
## distribution of mismatches

if(dev.cur()>1) dev.off()

quartz(width=16,height=9, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 0.5))
update_geom_defaults("point", aes(size = 2))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#1b9e77','#d95f02', '#7570b3')

ggplot(mToPlotSub, aes(x = value)) + geom_density(aes(fill = variable), alpha = 0.5)
last_plot() + scale_x_continuous(limits = c(0, 8), breaks = seq(from = 0, to = 8, by = 1))
last_plot() + scale_fill_manual(values=cbPalette)
xlab = "\nPosition (AA)"
ylab = "Density\n"
last_plot() + labs(x = xlab, y = ylab)
last_plot() + facet_wrap(~variable)
last_plot() + geom_vline(xintercept = median(toPlot$observed))
#last_plot() + stat_function(fun = dnorm, n = length(toPlot$observed), args = list(mean = mean(toPlot$random), sd = 1/3*mean(toPlot$random))) + ylab("")

p = last_plot()
setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("FigS7A_", toPlotName, "_mismatch_junction_positions.pdf") #

ggsave(filename, p)

#####################
## statistics
library(psych)

var1 = subset(mToPlotSub, variable == "mismatch")
var2 = subset(mToPlotSub, variable == "vEnd")
var3 = subset(mToPlotSub, variable == "jStart")

median(var1$value) # 4
median(var2$value) # 4
median(var3$value) # 6

mean(var1$value) # 4.272923
mean(var2$value) # 3.737864
mean(var3$value) # 5.648328

ks.test(var1$value, var2$value, alternative = c("greater")) # 1
ks.test(var1$value, var3$value, alternative = c("greater")) # < 2.2e-16
ks.test(var1$value, var2$value, alternative = c("less")) # 2.109e-15
ks.test(var1$value, var3$value, alternative = c("less")) # 1

## bargraph of individual mismatches

a = c(680, 381, 299)
b = c(rep(927, 3))
c = data.frame(a,b)
c$set = c("greater or equal", "greater", "equal")
colnames(c) = c("match", "total", "set")
c$frxn = a/b

if(dev.cur()>1) dev.off()

quartz(width=8,height=12, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 0.5))
update_geom_defaults("point", aes(size = 2))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#1b9e77','#d95f02', '#7570b3')

ggplot(c, aes(x = set, y = frxn)) +  geom_bar(stat="identity", fill = "#1b9e77")
xlab = "\nMismatch position \nversus V gene C-term"
ylab = "Fraction mismatches\n"
last_plot() + labs(x = xlab, y = ylab)
last_plot() + theme(axis.text.x=element_text(angle=45, hjust=1))

p = last_plot()
setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("FigS7B_", toPlotName, "_mismatch_junction_positions_relative_to_V_gene_end.pdf") #

ggsave(filename, p)
