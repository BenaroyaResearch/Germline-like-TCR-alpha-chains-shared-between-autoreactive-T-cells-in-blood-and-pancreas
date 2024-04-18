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
library(readxl)

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data") # folder containing data required to generate figures

## load TCRs with pgen scores

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data") # folder containing data required to generate figures
pgen = read.csv("Table_S4_Compiled_and_filtered_TCR_sequencces_used_in_this_study.csv")

## load PIT matched TCRs

levIAR1 = read.csv("Levenshtein_index_IAR_CD4_with_islet_TCRS_lv.lt9.csv", stringsAsFactors = F)

filename = ("Levenshtein_index_P324_P474_IAR_CD4_with_islet_TCRS_lv.lt6.csv")

levIAR2 = read.csv(filename, stringsAsFactors = F)

levSubIAR1 = subset(levIAR1, levIAR1$lv <2); nrow(levSubIAR1) # 1389 (573 unique) with IAR1; 1343 (555 unique(with IAR2)
levSubIAR2 = subset(levIAR2, levIAR2$lv <2); nrow(levSubIAR2) # 1389 (573 unique) with IAR1; 1343 (555 unique(with IAR2)

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

tcrsCombSub3Tra = subset(tcrsCombSub3, tcrsCombSub3$chainType == "TRA") # 3907
tcrsCombSub3Tra.u = tcrsCombSub3Tra[!duplicated(tcrsCombSub3Tra$junction),] # 3264

tcrsCombSub3.u = tcrsCombSub3[!duplicated(tcrsCombSub3$junction),] # 6481
tcrsTra = subset(tcrsCombSub3.u, chainType == "TRA") # 1512
tcrsTrb = subset(tcrsCombSub3.u, chainType == "TRB") # 1433

tcrsCombSub3Tra$studyGroup = gsub("AAbNeg", "HC", tcrsCombSub3Tra$studyGroup)

## select TCR subset to use

all = subset(tcrsCombSub3, set == "IAR1") # 5338
all.u = all[!duplicated(all$junction),] # 4331

## add pgen scores to all

toPlot = all.u
toPlot$pgen = pgen$pgen[match(toPlot$junction, pgen$junction)]

toPlotTra = subset(toPlot, chainType == "TRA") # 2174 unique
#toPlotTrb = subset(toPlot, chainType == "TRB") # 2136

## select TCRss to use

traMatchT = subset(toPlotTra,PITmatch == "TRUE") # 582
traMatchF = subset(toPlotTra,PITmatch == "FALSE") # 1592

## identify subset of TRB chains that pair with PITmatched TRA chains & add PIT matching status for TRA chains

toPlotTrb = subset(toPlot, libid %in% toPlotTra$libid & chainType == "TRB") # 1785
toPlotTrb$PITmatch = ifelse(toPlotTrb$libid %in% traMatchT$libid, "TRUE", "FALSE") # 500 TRUE, 1285 FALSE

## subset to PITmatch T and F

trbMatchT = subset(toPlot,PITmatch == "TRUE") # 492
trbMatchF = subset(toPlot,PITmatch == "FALSE") # 1070

trbMatchT = subset(toPlotTrb,PITmatch == "TRUE") # 500
trbMatchF = subset(toPlotTrb,PITmatch == "FALSE") # 1285

###################### 
## 
toPlotName = c("toPlotTra")
toPlot2 = get(toPlotName) # select TRA or TRB chains

toPlot2$PITmatch = gsub("TRUE", "PIT-matched", toPlot2$PITmatch)
toPlot2$PITmatch = gsub("FALSE", "non-PIT-matched", toPlot2$PITmatch)

t = subset(toPlot2, PITmatch == "PIT-matched")
f = subset(toPlot2, PITmatch == "non-PIT-matched")
ks.test(t$pgen, f$pgen)

## pgen

if(dev.cur()>1) dev.off()
quartz(width=13,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 2))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(toPlot2, aes(x = -log10(pgen))) + geom_density(aes(fill = PITmatch), alpha = 0.5, adjust = 1) 
last_plot() + theme(legend.spacing.y = unit(0.5, 'cm')) + guides(fill = guide_legend(byrow = TRUE))

xlab = "\n-log10(Pgen)"
ylab = "Density\n"
last_plot() + labs(x = xlab, y = ylab)
last_plot() + geom_vline(xintercept = -log10(median(traMatchT$pgen, na.rm = T)), lty = "solid")
last_plot() + geom_vline(xintercept = -log10(median(traMatchF$pgen, na.rm = T)), lty = "dashed")
last_plot() + scale_fill_manual(values=cbPalette, name = "PIT match")

last_plot() + scale_x_continuous(limits = c(5, 25))

p = last_plot()

forLab = ggplot_build(p)

Y = 0.9*(max(forLab$data[[1]]$density))
X = 1*(min((forLab$data[[1]]$x)))

label = paste("KS test, \np-value =", scientific(test1$p.value, 2), sep = " ")
#last_plot() + annotate("text", X, Y, label = label, size = 8, hjust = 0)

p = last_plot()
setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("FigS5A_", toPlotName, "_Pgen_by_PITmatch.pdf") #TRA

ggsave(filename, p)

###################### 
## for TRB junctions

toPlotName = c("toPlotTrb")
toPlot2 = get(toPlotName) # select TRA or TRB chains
toPlot2$PITmatch = gsub("TRUE", "PIT-matched", toPlot2$PITmatch)
toPlot2$PITmatch = gsub("FALSE", "non-PIT-matched", toPlot2$PITmatch)

t = subset(toPlot2, PITmatch == "PIT-matched")
f = subset(toPlot2, PITmatch == "non-PIT-matched")
ks.test(t$pgen, f$pgen)

## pgen

if(dev.cur()>1) dev.off()
quartz(width=13,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 2))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(toPlot2, aes(x = -log10(pgen))) + geom_density(aes(fill = PITmatch), alpha = 0.5, adjust = 1) 
last_plot() + theme(legend.spacing.y = unit(0.5, 'cm')) + guides(fill = guide_legend(byrow = TRUE))

xlab = "\n-log10(Pgen)"
ylab = "Density\n"
last_plot() + labs(x = xlab, y = ylab)
last_plot() + geom_vline(xintercept = -log10(median(trbMatchT$pgen, na.rm = T)), lty = "solid")
last_plot() + geom_vline(xintercept = -log10(median(trbMatchF$pgen, na.rm = T)), lty = "dashed")
last_plot() + scale_fill_manual(values=cbPalette, name = "PIT match")
last_plot() + scale_x_continuous(limits = c(5, 25))

p = last_plot()

forLab = ggplot_build(p)

Y = 0.9*(max(forLab$data[[1]]$density))
X = 1*(min((forLab$data[[1]]$x)))

label = paste("KS test, \np-value =", scientific(test1$p.value, 2), sep = " ")
#last_plot() + annotate("text", X, Y, label = label, size = 8, hjust = 0)

p = last_plot()
setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("FigS5B_", toPlotName, "_Pgen_by_PITmatch.pdf") #TRA

ggsave(filename, p)


p.adjust(c(2.2e-16, 0.01376))
# [1] 4.400e-16 1.376e-02