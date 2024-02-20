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

colNames = c("libid", "full_nt_sequence", "v_gene", "j_gene", "junction", "project", "donor_id", "set", "study_group")

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

colNames1 = c("libid", "full_nt_sequence", "v_gene", "j_gene", "junction", "project", "donor_id", "set", "study_group", "hla")

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

colNames = c("libid", "full_nt_sequence", "v_gene", "j_gene", "junction", "project", "donor_id", "set", "study_group", "hla")

IAR2Tcrs = IAR2Tcrs[c(colNames)]

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_code")

## combine set 1 amd set2

tcrsComb = rbind(IAR1Tcrs, IAR2Tcrs) # 7727

tcrsComb$studyGroup = tcrsComb$study_group

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
tcrsComb = subset(tcrsComb, !chainType == "other")

tcrsComb$study_group = gsub("AAbNeg", "HC", tcrsComb$study_group)
tcrsComb$study_group = gsub("1AAb", "AAb+", tcrsComb$study_group)
tcrsComb$study_group = gsub("2AAb", "AAb+", tcrsComb$study_group)

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

table(tcrsCombSub3$PITmatch, tcrsCombSub3$chainType)
       
#         TRA  TRB
#  FALSE 2672 3694
#  TRUE  1235    0

toSaveTra = subset(tcrsCombSub3, chainType == "TRA")
aMatchT = subset(toSaveTra, PITmatch == "TRUE") # 1235
toSaveTrb = subset(tcrsCombSub3, chainType == "TRB")
toSaveTrb$PITmatch = toSaveTrb$libid %in% aMatchT$libid
tcrsCombSub3 = rbind(toSaveTra, toSaveTrb)
table(tcrsCombSub3$PITmatch, tcrsCombSub3$chainType)
       
#         TRA  TRB
#  FALSE 2672 2590
#  TRUE  1235 1104

tcrsCombSub3Tra = subset(tcrsCombSub3, tcrsCombSub3$chainType == "TRA") # 3907
tcrsCombSub3Tra.u = tcrsCombSub3Tra[!duplicated(tcrsCombSub3Tra$junction),] # 3264

tcrsCombSub3.u = tcrsCombSub3[!duplicated(tcrsCombSub3$junction),] # 6451

table(tcrsCombSub3.u$PITmatch, tcrsCombSub3.u$chainType)
       
#         TRA  TRB
#  FALSE 2322 2287
#  TRUE   942  900

tcrsTra = subset(tcrsCombSub3.u, chainType == "TRA") # 3264
tcrsTrb = subset(tcrsCombSub3.u, chainType == "TRB") # 3187

## load healthy (and COVID) repertoires and filter sequencces

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data")

su = read.csv("Su_et_al_filtered_TCRs.csv") # 326782
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

## get pancreatic TCR seqs from Maki Nakayama, MAKI.NAKAYAMA@CUANSCHUTZ.EDU

mn = readxl::read_excel("220519_scPCR_Illumina.xlsx")
mn$barcode = seq(1, nrow(mn), by = 1)

comId = c("barcode", "Illumina", "Group", "Case_Tissue", "Subset", "Cell ID" )
traId = c("Vgene...6", "Jgene...7", "Junction...8", "frame...9")
trbId = c("Vgene...16", "Jgene...17", "Junction...18", "frame...19")

mnA = mn[,c(comId, traId)]
mnB = mn[,c(comId, trbId)]

colnames(mnA) = c(comId, "v_gene", "j_gene", "junction", "frame")
colnames(mnB) = c(comId, "v_gene", "j_gene", "junction", "frame")

mnA$chain = paste(mnA$v_gene, mnA$junction, mnA$j_gene, sep = "_")
mnB$chain = paste(mnB$v_gene, mnB$junction, mnB$j_gene, sep = "_")

mnTcrs = rbind(mnA, mnB) # 18684
mnTcrs = subset(mnTcrs, frame == "in-frame") # 14049
mnTcrs$chainType = ifelse(mnTcrs$v_gene %in% grep("TRA", mnTcrs$v_gene, value = T), "TRA",
						ifelse(mnTcrs$v_gene %in% grep("TRB", mnTcrs$v_gene, value = T), "TRB", "other"))
mnTcrs = subset(mnTcrs, !chainType == "other")

length(unique(mnTcrs$junction)) # 9798  

mn.u = mnTcrs[!duplicated(mnTcrs$chain),] # 

## load HLA for pancreatic TCR seqs

mnHla = readxl::read_excel("HLA_Typing.xlsx")
colnames(mnHla)[1] = c("Sample.ID")

## bring patients IDs into harmony
mnHla$Sample.ID = gsub(" ", "", mnHla$Sample.ID)
mnHla$Sample.ID = gsub("IIDP", "iidp", mnHla$Sample.ID)
mn$Case_Tissue = gsub("_islets", "", mn$Case_Tissue)
mn$Case_Tissue = gsub("_Islets", "", mn$Case_Tissue)

mnTcrs$Case_Tissue = gsub("_islets", "", mnTcrs$Case_Tissue)
mnTcrs$Case_Tissue = gsub("_Islets", "", mnTcrs$Case_Tissue)

mnTcrs$drb1 = mnHla$DRB1...4[match(mnTcrs$Case_Tissue, mnHla$Sample.ID)]
mnTcrs$drb1.1 = mnHla$DRB1...5[match(mnTcrs$Case_Tissue, mnHla$Sample.ID)]

hla0401 = subset(mnTcrs, mnTcrs$drb1 %in% grep("04:01", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("04:01", mnTcrs$drb1.1, value = T))
hla0301 = subset(mnTcrs, mnTcrs$drb1 %in% grep("03:01", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("03:01", mnTcrs$drb1.1, value = T))
hla0404 = subset(mnTcrs, mnTcrs$drb1 %in% grep("04:04", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("04:04", mnTcrs$drb1.1, value = T))
hla0101 = subset(mnTcrs, mnTcrs$drb1 %in% grep("01:01", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("01:01", mnTcrs$drb1.1, value = T))

mnTcrs$hla = ifelse(mnTcrs$drb1 %in% grep("03:01", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("03:01", mnTcrs$drb1.1, value = T), "03:01", 
				ifelse(mnTcrs$drb1 %in% grep("07:01", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("07:01", mnTcrs$drb1.1, value = T),
"07:01", 				
					ifelse(mnTcrs$drb1 %in% grep("04:01", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("04:01", mnTcrs$drb1.1, value = T), "04:01",
						ifelse(mnTcrs$drb1 %in% grep("01:01", mnTcrs$drb1, value = T) | mnTcrs$drb1.1 %in% grep("01:01", mnTcrs$drb1.1, value = T), "01:01", "other"))))

mnTcrs$chainType = ifelse(mnTcrs$v_gene %in% grep("TRA", mnTcrs$v_gene, value = T), "TRA",
						ifelse(mnTcrs$v_gene %in% grep("TRB", mnTcrs$v_gene, value = T), "TRB", "other"))
mnTcrs = subset(mnTcrs, !chainType == "other")

mnTcrs$cellType = mnTcrs$Subset
mnTcrs$cellType = ifelse(mnTcrs$cellType %in% grep("CD4", mnTcrs$cellType, value = T), "CD4",
						ifelse(mnTcrs$cellType %in% grep("CD8", mnTcrs$cellType, value = T), "CD8", "other"))

mnTcrs.u = mnTcrs[!duplicated(mnTcrs$junction),] # 9758 
						
mn$hla = mnTcrs$hla[match(mn$Case_Tissue, mnTcrs$Case_Tissue)]

toSave = mnTcrs.u
toSave$chainType = ifelse(toSave$v_gene %in% grep("TRA", toSave$v_gene, value = T), "TRA",
						ifelse(toSave$v_gene %in% grep("TRB", toSave$v_gene, value = T), "TRB", "other"))
toSave = subset(toSave, !chainType == "other")

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data")
filename = "Table_S2_Unique_compiled_and_filtered_PIT_TCR_sequences_used_in_this_study.csv"
#write.csv(toSave, filename)

## modify combined TCRs

tcrsCombSub3$PITmatch = tcrsCombSub3$junction %in% levSubIAR1$aJunc1 | tcrsCombSub3$junction %in% levSubIAR2$aJunc1 # 1235 T, 6398 F
tcrsCombSub3$expanded = factor(tcrsCombSub3$expanded, levels = c("E", "NE")) # 

tcrsCombSub3Tra = subset(tcrsCombSub3, tcrsCombSub3$chainType == "TRA") # 1814
tcrsCombSub3Tra.u = tcrsCombSub3Tra[!duplicated(tcrsCombSub3Tra$junction),] # 1512

tcrsCombSub3.u = tcrsCombSub3[!duplicated(tcrsCombSub3$junction),] # 2967
tcrsTra = subset(tcrsCombSub3.u, chainType == "TRA") # 1512
tcrsTrb = subset(tcrsCombSub3.u, chainType == "TRB") # 1433

tcrsCombSub3Tra$studyGroup = gsub("AAbNeg", "HC", tcrsCombSub3Tra$studyGroup)

## save tcrsComb as Supplemental Table for manuscript

toSave = tcrsCombSub3.u
toSave$set = gsub("IAR1", "Cohort1", toSave$set)
toSave$set = gsub("IAR2", "Cohort2", toSave$set)
toSave = subset(toSave, !chainType == "other") # 7601
table(toSave$PITmatch, toSave$chainType)
       
#         TRA  TRB
#  FALSE 2672 3694
#  TRUE  1235    0

toSaveTra = subset(toSave, chainType == "TRA")
aMatchT = subset(toSaveTra, PITmatch == "TRUE") # 942
toSaveTrb = subset(toSave, chainType == "TRB")
toSaveTrb$PITmatch = toSaveTrb$libid %in% aMatchT$libid
toSave = rbind(toSaveTra, toSaveTrb)
toSave$set = tcrsCombSub3$set[match(toSave$libid, tcrsCombSub3$libid)]

table(toSave$PITmatch, toSave$chainType)
       
#         TRA  TRB
#  FALSE 2322 2355
#  TRUE   942  832

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data")
filename = "Table_S2_Unique_compiled_and_filtered_TCR_sequences_used_in_this_study.csv"
#write.csv(toSave, filename)

## select TCR subset to use

all = tcrsCombSub3
all.u = tcrsCombSub3.u

allSub = all

## overlap of unique junctions

var1 = mnTcrs.u$junction %in% intersect(mnTcrs$junction, all.u$junction)
table(var1)
#var1
#FALSE  TRUE 
# 6355    96 

var1 = all.u$junction %in% intersect(mnTcrs$junction, all.u$junction)
table(var1)

#FALSE  TRUE 
# 6355    96 

##################################
## violin plots junctions per donor

tcrsCombSub3$PITmatch = tcrsCombSub3$junction %in% levSubIAR1$aJunc1 | tcrsCombSub3$junction %in% levSubIAR2$aJunc1 # 1235 T, 6266 F
tcrsCombSub3$PITmatch = gsub("TRUE", "PIT-matched", tcrsCombSub3$PITmatch)
tcrsCombSub3$PITmatch = gsub("FALSE", "non-PIT-matched", tcrsCombSub3$PITmatch)
tcrsCombSub3$PITmatch = factor(tcrsCombSub3$PITmatch, levels = c("PIT-matched", "non-PIT-matched"))

tcrsCombSub3$expanded = factor(tcrsCombSub3$expanded, levels = c("E", "NE")) # 

## focus on TRA chains

tcrsCombSub3Tra = subset(tcrsCombSub3, tcrsCombSub3$chainType == "TRA") # 3907
tcrsCombSub3Tra.u = tcrsCombSub3Tra[!duplicated(tcrsCombSub3Tra$junction),] # 3264
#tcrsCombSub3$study_group = factor(tcrsCombSub3$studyGroup, levels = c("HC", "AAb+", "newT1D", "T1D"))
tcrsCombSub3$expanded = factor(tcrsCombSub3$expanded, levels = c("E", "NE")) # 1874 E, 5727 NE

table(tcrsCombSub3$chainType, tcrsCombSub3$expanded)
     
#         E   NE
#  TRA 1000 2907
#  TRB  874 2820

table(tcrsCombSub3Tra$chainType, tcrsCombSub3Tra$study_group)
     
#      AAb+   HC newT1D  T1D
#  TRA  632  575   1859  841

dr4 = subset(tcrsComb, tcrsComb$hla %in% grep("04", tcrsComb$hla, value = T))

#####################
## select PIT match and expanded groups if desired

toUse = subset(tcrsCombSub3Tra, expanded == "E"); nrow(toUse) #1000 E, 2907 NE

## tally numbers per donor

totNoPerDonor =  ddply(toUse,.(donor_id, PITmatch, study_group), plyr::summarize, sum = length(donor_id))

## add total cells per donor and calculate fractions

total = ddply(totNoPerDonor, .(donor_id), plyr::summarize, sum = sum(sum))
totNoPerDonor$total = total$sum[match(totNoPerDonor$donor_id, total$donor_id)]

totNoPerDonor$frxnPITmatch = totNoPerDonor$sum/totNoPerDonor$total

## for plots
donorSub = totNoPerDonor
donorSub$DRB1 = ifelse(donorSub$donor_id %in% dr4$donor_id, "DRB1-04", "other")

donorSub$PITmatch = gsub("TRUE", "PIT-matched", donorSub$PITmatch)
donorSub$PITmatch = gsub("FALSE", "non-PIT-matched", donorSub$PITmatch)

donorSub$study_group = factor(donorSub$study_group, levels = c("HC", "AAb+", "newT1D", "T1D"))
donorSub$PITmatch = factor(donorSub$PITmatch, levels = c("PIT-matched", "non-PIT-matched"))

## check whether small cell yields affect results

#gt10 = subset(donorSub, total >=10)
#donorSub = subset(donorSub, donor_id %in% gt10$donor_id) # to 49 from 81

library(ggsignif)

if(dev.cur()>1) dev.off()
quartz(width=16,height=9, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 0.5))
update_geom_defaults("point", aes(size = 4))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(donorSub, aes(x = study_group, y = frxnPITmatch)) + geom_violin(draw_quantiles = c(0.5), stat = "ydensity", position = "dodge", fill = "#eeeeee", lwd = 1.5) + facet_grid(~PITmatch)
#last_plot() + geom_jitter(data = donorSub, aes(x = study_group, y = frxnPITmatch, colour = DRB1))#
last_plot() + geom_dotplot(data = donorSub, aes(x = study_group, y = frxnPITmatch, fill = DRB1), position = "identity", method = "histodot", binaxis = "y", stackdir = "center", dotsize = 0.75)
last_plot() + scale_y_continuous(limits = c(0.0, 1.125), breaks = c(0.0, 0.25, 0.5, 0.75, 1.0))

last_plot() + scale_fill_manual(values=cbPalette)
last_plot() + theme(axis.text.x=element_text(angle=45, hjust=1))

last_plot() + theme(axis.text.x=element_text(angle=45, hjust=1))

xlab = "\nGroup"
ylab = "PIT matches, \nfraction IAR cells per donor\n"
last_plot() + labs(x = xlab, y = ylab)

ymax = max(donorSub$frxnPITmatch + 0.125) # 0.125
ymin = min(donorSub$frxnPITmatch)

p = last_plot()

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("FigS3_Fraction", "_PIT_matches_per_donor_E_TCRs.pdf")
ggsave(filename, p)

## statistics- adjust for groups tested

var1 = subset(donorSub, donorSub$study_group == "HC" & donorSub$PITmatch == "PIT-matched")
var2 = subset(donorSub, donorSub$study_group == "AAb+" & donorSub$PITmatch == "PIT-matched")
var3 = subset(donorSub, donorSub$study_group == "newT1D" & donorSub$PITmatch == "PIT-matched")
var4 = subset(donorSub, donorSub$study_group == "T1D" & donorSub$PITmatch == "PIT-matched")

test1 = wilcox.test(var1$frxnPITmatch, var2$frxnPITmatch); test1
test2 = wilcox.test(var1$frxnPITmatch, var3$frxnPITmatch); test2
test3 = wilcox.test(var1$frxnPITmatch, var4$frxnPITmatch); test3
test4 = wilcox.test(var2$frxnPITmatch, var3$frxnPITmatch); test4
test5 = wilcox.test(var2$frxnPITmatch, var4$frxnPITmatch); test5
test6 = wilcox.test(var3$frxnPITmatch, var4$frxnPITmatch); test6

padj = p.adjust(c(test1$p.value, test2$p.value, test3$p.value, test4$p.value), method = "BH"); padj
#[1] 0.09783087 0.46465182 0.91303871 0.21375656 # E PIT-matched
#[1] 0.2033694 0.4640637 0.8289894 0.8163858 # E non-PIT-matched
#[1] 0.094756480 0.002122289 0.533739322 0.094756480 # NE PIT-matched
#[1] 0.30104874 0.04889583 0.55778895 0.14213472 # NE non-PIT-matched


