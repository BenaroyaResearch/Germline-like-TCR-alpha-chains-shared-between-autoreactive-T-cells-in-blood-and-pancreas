rm(list = ls())
require(Biostrings)
require(limma)
require(plyr)
require(readxl)
library(reshape2)
library(psych)
library(ggplot2)
library(r2symbols)
library(stringdist)

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data") # folder containing data required to generate figures

## load matches

levIAR1 = read.csv("Levenshtein_index_IAR_CD4_with_islet_TCRS_lv.lt9.csv", stringsAsFactors = F)

filename = ("Levenshtein_index_P324_P474_IAR_CD4_with_islet_TCRS_lv.lt6.csv")

levIAR2 = read.csv(filename, stringsAsFactors = F)

levSubIAR1 = subset(levIAR1, levIAR1$lv <2); nrow(levSubIAR1) # 1389 (573 unique) with IAR1; 1343 (555 unique(with IAR2)
levSubIAR2 = subset(levIAR2, levIAR2$lv <2); nrow(levSubIAR2) # 1389 (573 unique) with IAR1; 1343 (555 unique(with IAR2)

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
tcrsCombSub3 = subset(tcrsCombSub3, !chainType == "other") # 7601

## modify combined TCRs

tcrsCombSub3$PITmatch = tcrsCombSub3$junction %in% levSubIAR1$aJunc1 | tcrsCombSub3$junction %in% levSubIAR2$aJunc1 # 1235 T, 6366 F
tcrsCombSub3$expanded = factor(tcrsCombSub3$expanded, levels = c("E", "NE")) # 
tcrsCombSub3 = subset(tcrsCombSub3, !chainType == "other")

tcrsCombSub3.u = tcrsCombSub3[!duplicated(tcrsCombSub3$junction),] # 6451

## determine junction sharing

publicity =  ddply(tcrsCombSub3,.(junction), plyr::summarize, sum = length(unique(donor_id))) # 6451 junctions, sum(publicity$sum) = 7601

tcrsCombSub3.u$publicity = publicity$sum[match(tcrsCombSub3.u$junction, publicity$junction)]

tcrsCombSub3.u$pubT = tcrsCombSub3.u$publicity>1
table(tcrsCombSub3.u$pubT )

#FALSE  TRUE 
# 6347   104 

tcrsCombSub3.u$privT = tcrsCombSub3.u$publicity ==1 & tcrsCombSub3.u$expanded == "E"
table(tcrsCombSub3.u$privT)

#FALSE  TRUE 
# 5831   620

a = subset(tcrsCombSub3.u, chainType == "TRA") # 3264
aMatchT = subset(a, PITmatch == "TRUE") # 942
b = subset(tcrsCombSub3.u, chainType == "TRB") #3187
b$PITmatch = b$libid %in% aMatchT$libid

df = table(a$pubT, a$PITmatch)
df
         
#  FALSE  2302  886
#  TRUE     20   56
fisher.test(df) # 3.475e-16

df = table(b$pubT, b$PITmatch)
df

#          FALSE TRUE
#  FALSE  2339  820
#  TRUE     16   12

fisher.test(df) # 0.05153

###################
## hypergeometric p-values

## TRA vs public

aUniv = length(unique(a$junction))# 3264
a.u = a[!duplicated(a$junction),] # 3264

aPub = subset(a.u, pubT == "TRUE") # 76
aPriv = subset(a.u, privT == "TRUE") # 281
aPrivS = subset(a.u, a.u$junction %in% sample(a.u$junction, size = length(aPub$junction))) # 76

aPitT = subset(a.u, PITmatch == "TRUE") # 942
aPitF = subset(a.u, PITmatch == "FALSE") # 2322

lenApub = length(aPub$junction) # 76
lenApriv = length(aPriv$junction) # 281
lenAprivS = length(aPrivS$junction) # 76

lenApitT = length(aPitT$junction) # 942
lenApitF = length(aPitF$junction) # 2322

pubI = length(intersect(aPub$junction, aPitT$junction)) # 56
privI = length(intersect(aPriv$junction, aPitT$junction)) # 94
privIs = length(intersect(aPrivS$junction, aPitT$junction)); privIs # 23.3 (mean of 10 determinations each)

phyper(pubI, lenApub, aUniv-lenApub, lenApitT, lower.tail = F) # 4.642685e-17
phyper(privI, lenApriv, aUniv-lenApriv, lenApitT, lower.tail = F) # 0.03374707
phyper(privIs, lenAprivS, aUniv-lenAprivS, lenApitT, lower.tail = F) # 0.7298442

## Venn diagram

library(VennDiagram)
library(ggplot2)

## public PITmatched

if(dev.cur() >1) dev.off()
quartz(height = 10, width = 10, dpi =72);
theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("errorbar", aes(size = 2))
update_geom_defaults("point", aes(size = 8))
theme_update(plot.title = element_text(hjust = 0.5))

A = "Public TRA"
B = "PIT-matched junctions"

#A = paste("islet", "\nintersect", "AgSp", sep = " ")
#B = paste("islet", "\nintersect", "Shingrix", sep = " ")

venn.plot <- draw.pairwise.venn(
  area1 = lenApub,
  area2 = lenApitT,
  cross.area = pubI,
  category = c(A, B),
  fill = c("mediumseagreen", "red3"),
  col = c("mediumseagreen", "red3"),
  cex = c(2,2,2),
  cat.cex = c(2,2),
  cat.pos = c(-20, 10),
  cat.dist = 0.09,
  fontfamily = "sans",
  cat.fontfamily = "sans"
)

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs")
filename = "FigS4A. PIT-matched TRA junctions are enriched with public TRA junctions.pdf"

quartz.save(filename, type = "pdf")
dev.off()

###################
## hypergeometric p-values

## TRA vs private

aUniv = length(unique(a$junction))# 3264
a.u = a[!duplicated(a$junction),] # 3264

aPub = subset(a.u, pubT == "TRUE") # 76
aPriv = subset(a.u, privT == "TRUE") # 281
aPrivS = subset(a.u, a.u$junction %in% sample(a.u$junction, size = length(aPub$junction))) # 76

aPitT = subset(a.u, PITmatch == "TRUE") # 942
aPitF = subset(a.u, PITmatch == "FALSE") # 2322

lenApub = length(aPub$junction) # 76
lenApriv = length(aPriv$junction) # 281
lenAprivS = length(aPrivS$junction) # 76

lenApitT = length(aPitT$junction) # 942
lenApitF = length(aPitF$junction) # 2322

pubI = length(intersect(aPub$junction, aPitT$junction)) # 56
privI = length(intersect(aPriv$junction, aPitT$junction)) # 94
privIs = length(intersect(aPrivS$junction, aPitT$junction)) # # 23.3 (mean of 10 determinations each)

phyper(pubI, lenApub, aUniv-lenApub, lenApitT, lower.tail = F) # 4.64e-17
phyper(privI, lenApriv, aUniv-lenApriv, lenApitT, lower.tail = F) # 0.034
phyper(privIs, lenAprivS, aUniv-lenAprivS, lenApitT, lower.tail = F) # 0.6375176

p.adjust(c(4.64e-17, 0.034, 0.6375176))
# 4.640000e-17 3.400000e-02 6.375176e-01
## Venn diagram

library(VennDiagram)
library(ggplot2)

## public PITmatched
if(dev.cur() >1) dev.off()
quartz(height = 10, width = 10, dpi =72);
theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("errorbar", aes(size = 2))
update_geom_defaults("point", aes(size = 8))
theme_update(plot.title = element_text(hjust = 0.5))

A = "Private TRA"
B = "PIT-matched junctions"

#A = paste("islet", "\nintersect", "AgSp", sep = " ")
#B = paste("islet", "\nintersect", "Shingrix", sep = " ")

venn.plot <- draw.pairwise.venn(
  area1 = lenAprivS,
  area2 = lenApitT,
  cross.area = privIs,
  category = c(A, B),
  fill = c("mediumseagreen", "red3"),
  col = c("mediumseagreen", "red3"),
  cex = c(2,2,2),
  cat.cex = c(2,2),
  cat.pos = c(-20, 10),
  cat.dist = 0.09,
  fontfamily = "sans",
  cat.fontfamily = "sans"
)

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs")
filename = "FigS4B. PIT-matched TRA junctions are enriched with public TRA junctions.pdf"

quartz.save(filename, type = "pdf")
dev.off()

###################
## hypergeometric p-values

## TRB vs public

bUniv = length(unique(b$junction))# 3187
b.u = b[!duplicated(b$junction),] # 3187

bPub = subset(b.u, pubT == "TRUE") # 28
bPriv = subset(b.u, privT == "TRUE") # 339
bPrivS = subset(b.u, b.u$junction %in% sample(b.u$junction, size = length(bPub$junction))) # 28

bPitT = subset(b.u, PITmatch == "TRUE") # 832
bPitF = subset(b.u, PITmatch == "FALSE") # 2355

lenBpub = length(bPub$junction) # 28
lenBpriv = length(bPriv$junction) # 339
lenBprivS = length(bPrivS$junction) # 28

lenBpitT = length(bPitT$junction) # 832
lenBpitF = length(bPitF$junction) # 2355

pubI = length(intersect(bPub$junction, bPitT$junction)) # 12
privI = length(intersect(bPriv$junction, bPitT$junction)) # 102
privIs = length(intersect(bPrivS$junction, bPitT$junction)) # 6 

phyper(pubI, lenBpub, bUniv-lenBpub, lenBpitT, lower.tail = F) # 0.01583616
phyper(privI, lenBpriv, bUniv-lenBpriv, lenBpitT, lower.tail = F) # 0.03484931
phyper(privIs, lenBprivS, bUniv-lenBprivS, lenBpitT, lower.tail = F) # 0.6248287

p.adjust(c(0.01583616, 0.03484931, 0.6248287))
# [1] 0.04750848 0.06969862 0.62482870

## Venn diagram

library(VennDiagram)
library(ggplot2)

## public PITmatched
if(dev.cur() >1) dev.off()
quartz(height = 10, width = 10, dpi =72);
theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("errorbar", aes(size = 2))
update_geom_defaults("point", aes(size = 8))
theme_update(plot.title = element_text(hjust = 0.5))

A = "Public TRB"
B = "PIT-matched junctions"

venn.plot <- draw.pairwise.venn(
  area1 = lenBpub,
  area2 = lenBpitT,
  cross.area = pubI,
  category = c(A, B),
  fill = c("mediumseagreen", "red3"),
  col = c("mediumseagreen", "red3"),
  cex = c(2,2,2),
  cat.cex = c(2,2),
  cat.pos = c(-20, 10),
  cat.dist = 0.09,
  fontfamily = "sans",
  cat.fontfamily = "sans"
)

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs")
filename = "FigS4C. PIT-matched TRB junctions were enriched with public TRB junctions.pdf"

quartz.save(filename, type = "pdf")
dev.off()

###################
## hypergeometric p-values

## TRB vs private

bUniv = length(unique(b$junction))# 3187
b.u = b[!duplicated(b$junction),] # 3187

bPub = subset(b.u, pubT == "TRUE") # 28
bPriv = subset(b.u, privT == "TRUE") # 339
bPrivS = subset(b.u, b.u$junction %in% sample(b.u$junction, size = length(bPub$junction))) # 28

bPitT = subset(b.u, PITmatch == "TRUE") # 832
bPitF = subset(b.u, PITmatch == "FALSE") # 2355

lenBpub = length(bPub$junction) # 28
lenBpriv = length(bPriv$junction) # 339
lenBprivS = length(bPrivS$junction) # 28

lenBpitT = length(bPitT$junction) # 832
lenBpitF = length(bPitF$junction) # 2355

pubI = length(intersect(bPub$junction, bPitT$junction)) # 12
privI = length(intersect(bPriv$junction, bPitT$junction)) # 102
privIs = length(intersect(bPrivS$junction, bPitT$junction)) # 6 

phyper(pubI, lenBpub, bUniv-lenBpub, lenBpitT, lower.tail = F) # 0.01583616
phyper(privI, lenBpriv, bUniv-lenBpriv, lenBpitT, lower.tail = F) # 0.03484931
phyper(privIs, lenBprivS, bUniv-lenBprivS, lenBpitT, lower.tail = F) # 0.6248287

## Venn diagram

library(VennDiagram)
library(ggplot2)

## public PITmatched
if(dev.cur() >1) dev.off()
quartz(height = 10, width = 10, dpi =72);
theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("errorbar", aes(size = 2))
update_geom_defaults("point", aes(size = 8))
theme_update(plot.title = element_text(hjust = 0.5))

A = "Private TRB"
B = "PIT-matched junctions"

#A = paste("islet", "\nintersect", "AgSp", sep = " ")
#B = paste("islet", "\nintersect", "Shingrix", sep = " ")

venn.plot <- draw.pairwise.venn(
  area1 = lenBprivS,
  area2 = lenBpitT,
  cross.area = privIs,
  category = c(A, B),
  fill = c("mediumseagreen", "red3"),
  col = c("mediumseagreen", "red3"),
  cex = c(2,2,2),
  cat.cex = c(2,2),
  cat.pos = c(-20, 10),
  cat.dist = 0.09,
  fontfamily = "sans",
  cat.fontfamily = "sans"
)

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs")
filename = "FigS4D. PIT-matched TRB junctions were enriched with public TRB junctions.pdf"

quartz.save(filename, type = "pdf")
dev.off()

