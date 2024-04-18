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
library(readxl)

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data") # folder containing data required to generate figures

## load TCRs

levIAR1 = read.csv("Levenshtein_index_IAR_CD4_with_islet_TCRS_lv.lt9.csv", stringsAsFactors = F)

filename = ("Levenshtein_index_P324_P474_IAR_CD4_with_islet_TCRS_lv.lt6.csv")

levIAR2 = read.csv(filename, stringsAsFactors = F)

levSubIAR1 = subset(levIAR1, levIAR1$lv <2); nrow(levSubIAR1) # 1389 (573 unique) with IAR1; 1343 (555 unique(with IAR2)
levSubIAR2 = subset(levIAR2, levIAR2$lv <2); nrow(levSubIAR2) # 1389 (573 unique) with IAR1; 1343 (555 unique(with IAR2)

## load IAR1 TCRs

IAR1Tcrs = read.csv("201512_TCR_MasterList_w_CloneIDs.csv", stringsAsFactors = F) # 5729
colnames(IAR1Tcrs) = gsub("tcrGraph_sharing_level", "sharing_level", colnames(IAR1Tcrs))

IAR1Tcrs = data.frame(IAR1Tcrs, set = "IAR1")

colNames = c("libid", "v_gene", "j_gene", "junction", "project", "donor_id", "set", "study_group")

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

colNames1 = c("libid", "v_gene", "j_gene", "junction", "project", "donor_id", "set", "study_group", "hla")

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

tcrsCombSub3$PITmatch = tcrsCombSub3$junction %in% levSubIAR1$aJunc1 | tcrsCombSub3$junction %in% levSubIAR2$aJunc1 # 1235 T, 6398 F
tcrsCombSub3$expanded = factor(tcrsCombSub3$expanded, levels = c("E", "NE")) # 

tcrsCombSub3Tra = subset(tcrsCombSub3, tcrsCombSub3$chainType == "TRA") # 3907
tcrsCombSub3Tra.u = tcrsCombSub3Tra[!duplicated(tcrsCombSub3Tra$junction),] # 3264

tcrsCombSub3Trb = subset(tcrsCombSub3, tcrsCombSub3$chainType == "TRB") # 3694
tcrsCombSub3Trb.u = tcrsCombSub3Tra[!duplicated(tcrsCombSub3Trb$junction),] # 3367

tcrsCombSub3.u = tcrsCombSub3[!duplicated(tcrsCombSub3$junction),] # 6481
tcrsTra = subset(tcrsCombSub3.u, chainType == "TRA") # 1512
tcrsTrb = subset(tcrsCombSub3.u, chainType == "TRB") # 1433

tcrsCombSub3Tra$studyGroup = gsub("AAbNeg", "HC", tcrsCombSub3Tra$studyGroup)

## select TCR subset to use

all = subset(tcrsCombSub3Trb, set == "IAR2") # 2725 IAR1, 1182 IAR2, 1104 IAR2 TRB
all.u = all[!duplicated(all$junction),] # 2174 IAR1, 1117 IAR2, 1053 IAR2 TRB

allSub = all.u

all$study_group = gsub("early onset T1D", "newT1D", all$study_group)

all.HC = subset(all, study_group == "HC") # 259
all.newT1D = subset(all, study_group == "newT1D") # 497
all.T1D = subset(all, study_group == "T1D") # 597

## load healthy (and COVID) repertoires and filter sequencces

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data")

su = read.csv("Su_et_al_filtered_TCRs.csv", stringsAsFactors = F) # 326782
su$X = NULL

suNkt1 = subset(su, junction == "CVVSDRGSTLGRLYF") #45
suNkt2 = subset(su, barcode %in% suNkt1$barcode) # 90
sumait1 = subset(su, v_gene == "TRAV1-2" & (j_gene == "TRAJ33" | j_gene == "TRAJ20" | j_gene == "TRAJ12")) # 2533
sumait2 = subset(su, barcode %in% sumait1$barcode) # 5106

su = subset(su, !junction %in% suNkt2$junction & !junction %in% sumait2$junction) # 320372

suTRA = subset(su, chainType == "TRA")
suTRB = subset(su, chainType == "TRB")

suCd4 = subset(su, su$cellType == "cd4")
suCd8 = subset(su, su$cellType == "cd8")

su.u = su[!duplicated(su$junction),]
suHc = subset(su.u, group == "HC")
suCovid = subset(su.u, group == "COVID")

suHcCd4 = subset(suHc, suHc$cellType == "cd4")
suHcCd8 = subset(suHc, suHc$cellType == "cd8") # 5825

suHcTra = subset(suHcCd8, chainType == "TRA") # 2528
suHcTrb = subset(suHcCd8, chainType == "TRB") # 3297

suHcTra$chainLen = nchar(suHcTra$junction)
suHcTrb$chainLen = nchar(suHcTrb$junction)

summary(suHcTra$chainLen)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   7.00   13.00   14.00   13.91   15.00   22.00 

summary(suHcTrb$chainLen)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   8.00   13.00   15.00   14.62   16.00   21.00 

## get pancreatic TCR seqs from Maki Nakayama, MAKI.NAKAYAMA@CUANSCHUTZ.EDU

mn = read_excel("220519_scPCR_Illumina.xlsx")

colnames(mn) = gsub("Cell ID", "Cell.ID", colnames(mn))

colnames(mn) = gsub("Vgene...6", "Vgene", colnames(mn))
colnames(mn) = gsub("Jgene...7", "Jgene", colnames(mn))
colnames(mn) = gsub("Junction...8", "Junction", colnames(mn))
colnames(mn) = gsub("frame...9", "frame", colnames(mn))

colnames(mn) = gsub("Vgene...16", "Vgene.2", colnames(mn))
colnames(mn) = gsub("Jgene...17", "Jgene.2", colnames(mn))
colnames(mn) = gsub("Junction...18", "Junction.2", colnames(mn))
colnames(mn) = gsub("frame...19", "frame.2", colnames(mn))

mn$barcode = seq(1, nrow(mn), by = 1)
mn$alphaChain = paste(mn$Vgene, mn$Junction, mn$Jgene, sep = "_")
mn$betaChain = paste(mn$Vgene.2, mn$Junction.2, mn$Jgene.2, sep = "_")

comId = c("barcode", "Illumina", "Group", "Case_Tissue", "Subset", "Cell.ID" )
traId = c("Vgene", "Jgene", "Junction", "frame")
trbId = c("Vgene.2", "Jgene.2", "Junction.2", "frame.2")

mnA = mn[,c(comId, traId)]
mnB = mn[,c(comId, trbId)]

colnames(mnA) = c(comId, "v_gene", "j_gene", "junction", "frame")
colnames(mnB) = c(comId, "v_gene", "j_gene", "junction", "frame")

mnTcrs = rbind(mnA, mnB) # 18684
mnTcrs = subset(mnTcrs, frame == "in-frame") # 14049

length(unique(mnTcrs$junction)) # 9798  

mnTcrs.u = mnTcrs[!duplicated(mnTcrs$junction),] # 9798 

## load HLA for pancreatic TCR seqs

mnHla = read_excel("HLA_Typing.xlsx")

## bring patients IDs into harmony
colnames(mnHla) = gsub("Sample ID", "Sample.ID", colnames(mnHla))
colnames(mnHla) = gsub("DRB1...4", "DRB1", colnames(mnHla))
colnames(mnHla) = gsub("DRB1...5", "DRB1.1", colnames(mnHla))

mnHla$Sample.ID = gsub(" ", "", mnHla$Sample.ID)
mnHla$Sample.ID = gsub("IIDP", "iidp", mnHla$Sample.ID)
mn$Case_Tissue = gsub("_islets", "", mn$Case_Tissue)
mn$Case_Tissue = gsub("_Islets", "", mn$Case_Tissue)

mnTcrs$drb1 = mnHla$DRB1[match(mnTcrs$Case_Tissue, mnHla$Sample.ID)]
mnTcrs$drb1.1 = mnHla$DRB1.1[match(mnTcrs$Case_Tissue, mnHla$Sample.ID)]

hla0401 = subset(mnTcrs, mnTcrs$drb1 %in% grep("04:01", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("04:01", mnTcrs$drb1.1, value = T))
hla0301 = subset(mnTcrs, mnTcrs$drb1 %in% grep("03:01", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("03:01", mnTcrs$drb1.1, value = T))
hla0404 = subset(mnTcrs, mnTcrs$drb1 %in% grep("04:04", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("04:04", mnTcrs$drb1.1, value = T))
hla0101 = subset(mnTcrs, mnTcrs$drb1 %in% grep("01:01", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("01:01", mnTcrs$drb1.1, value = T))

mnTcrs$hla = ifelse(mnTcrs$drb1 %in% grep("03:01", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("03:01", mnTcrs$drb1.1, value = T), "03:01", 
				ifelse(mnTcrs$drb1 %in% grep("07:01", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("07:01", mnTcrs$drb1.1, value = T),
"07:01", 				
					ifelse(mnTcrs$drb1 %in% grep("04:01", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("04:01", mnTcrs$drb1.1, value = T), "04:01",
						ifelse(mnTcrs$drb1 %in% grep("01:01", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("01:01", mnTcrs$drb1.1, value = T), "01:01", "other"))))
						
mn$hla = mnTcrs$hla[match(mn$Case_Tissue, mnTcrs$Case_Tissue)]

## reformant Maki's TCRs

comId = c("barcode", "Illumina", "Group", "Case_Tissue", "Subset", "Cell.ID", "hla" )
traId = c("Vgene", "Jgene", "Junction", "frame")
trbId = c("Vgene.2", "Jgene.2", "Junction.2", "frame.2")

mnA = mn[,c(comId, traId)]
mnB = mn[,c(comId, trbId)]

colnames(mnA) = c(comId, "v_gene", "j_gene", "junction", "frame")
colnames(mnB) = c(comId, "v_gene", "j_gene", "junction", "frame")

mnTcrs = rbind(mnA, mnB) # 18684
mnTcrs = subset(mnTcrs, frame == "in-frame") # 14049

length(unique(mnTcrs$junction)) # 9798  

mnTcrs.u = mnTcrs[!duplicated(mnTcrs$junction),] # 9798 

mn.u = mn[!duplicated(mn$alphaChain),] # 5039

## overlap of unique junctions

var1 = mnTcrs.u$junction %in% intersect(mnTcrs$junction, all.u$junction)

var1 = all.u$junction %in% intersect(mnTcrs$junction, all.u$junction)
table(var1)

#FALSE  TRUE 
# 1083    34 #IAR2 TRA; 1049, 4 for TRB

#####################
## plot Lv macthes IAR T cells versus Hc Rand

## load levenshtein index comparisons

levTra = read.csv("Levenshtein_index_P324_P474_IAR_CD4_with_islet_TCRS_lv.lt6.csv", stringsAsFactors = F) # 518011
levTrb = read.csv("Levenshtein_index_P325_P474_IAR_CD4_TRB_junctions_with_islet_TCRS_lv.lt6.csv", stringsAsFactors = F)# 887901

#lev = subset(levTra, aJunc1 %in% all.u$junction) # 354074
lev = subset(levTrb, aJunc1 %in% all.u$junction) # 354074


lev = subset(lev, lv<9)
levRand = read.csv("Levenshtein_index_IAR_CD4_with_islet_TCRS_lv.lt9_random_set.csv", stringsAsFactors = F) 
levRand = read.csv("Levenshtein_index_TRB_chains_downsampled_randHC.csv", stringsAsFactors = F) 

levRand = subset(levRand, lv<10)

iarAg =  ddply(lev,.(lv), plyr::summarize, sum = length(lv)) # 
hcAg =  ddply(levRand,.(lv), plyr::summarize, sum = length(lv)) # 

agMerge = merge(iarAg, hcAg, by = "lv")
colnames(agMerge)[1] = "LvIndex"

if(dev.cur()>1) dev.off()

quartz(width=13,height=10, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 6))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(agMerge, aes(y = log10(sum.x), x = log10(sum.y))) + geom_point(aes(size = LvIndex))
last_plot() + scale_x_continuous(limits = c(0,7))
last_plot() + scale_y_continuous(limits = c(0,7))

last_plot() + geom_abline(intercept = 0, slope = 1)
last_plot() + geom_smooth(method = lm)
ylab = "\nIAR vs. PIT TRB, \nlog10(No. matches)"
xlab = "HC vs. PIT, \nlog10(No. mismatches)\n"
last_plot() + labs(x = xlab, y = ylab)
#last_plot() + scale_x_continuous(limits = c(0,5))

p = last_plot()
setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs")
filename = "Fig_2SG_IAR T cell and HC TRA junctions comparisons over a range_of_LI_values_IAR2_TRB.pdf"
ggsave(filename, p)

########################
## test for significant difference between line slope and 1
## https://www.researchgate.net/post/Does-anyone-know-how-to-test-the-significant-difference-between-a-line-slopeeg-081-and-1-by-using-R-or-SPSS

## method 1
r.x <- lm(sum.x ~ sum.y, data = agMerge); summary(r.x)
r1 <- lm(sum.x ~1 + offset(sum.y), data = agMerge); summary(r1)
anova(r.x, r1)
#Analysis of Variance Table

#Model 1: sum.x ~ sum.y
#Model 2: sum.x ~ 1 + offset(sum.y)
#   Res.Df        RSS Df   Sum of Sq     F    Pr(>F)    
#1      5 7.6233e+06                                   
#2      6 4.0234e+10 -1 -4.0226e+10 26384 1.678e-10 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## method 2
summary(r.x)
plot(agMerge$sum.x,agMerge$sum.y)
summary(r.x)$coefficients
a=summary(r.x)$coefficients[1,1]
b=summary(r.x)$coefficients[2,1]
abline(a,b)
confint(r.x)
#                   2.5 %       97.5 %
#(Intercept) -868.59867 2027.67457
#sum.y          2.54227    2.59187

# thus, confidence interval for slope does NOT contain 1. Significant.

## calculate fraction of PIT matches

levComb = unique(c(levSubIAR1$aJunc1, levSubIAR2$aJunc1)) # 1106

tcrsCombSub3.u = tcrsCombSub3[!duplicated(tcrsCombSub3$junction),] # 6481
tcrsTra = subset(tcrsCombSub3.u, chainType == "TRA") # 3264
tcrsTrb = subset(tcrsCombSub3.u, chainType == "TRB") # 3187

#fracLevTra = 1106/3264 = 33.9%

## no. donors by cohort

noDonors =  ddply(tcrsComb,.(set, study_group), plyr::summarize, sum = length(unique(donor_id))) #
# noDonors
#   set study_group sum
#1 IAR1          HC  11
#2 IAR1      newT1D  26
#3 IAR1         T1D  16
#4 IAR2        1AAb   8
#5 IAR2        2AAb   6
#6 IAR2      AAbNeg   6
#7 IAR2      newT1D  11

