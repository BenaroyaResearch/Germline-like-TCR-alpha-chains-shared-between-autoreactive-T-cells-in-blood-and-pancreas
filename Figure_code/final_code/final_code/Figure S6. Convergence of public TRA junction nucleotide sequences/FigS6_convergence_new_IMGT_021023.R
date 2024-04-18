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
tcrsCombSub3 = subset(tcrsCombSub3, !chainType == "other") # 7601

tcrsCombSub3Tra = subset(tcrsCombSub3, tcrsCombSub3$chainType == "TRA") # 3907
tcrsCombSub3Tra.u = tcrsCombSub3Tra[!duplicated(tcrsCombSub3Tra$junction),] # 3264

tcrsCombSub3.u = tcrsCombSub3[!duplicated(tcrsCombSub3$junction),] # 6451
tcrsTra = subset(tcrsCombSub3.u, chainType == "TRA") # 3264
tcrsTrb = subset(tcrsCombSub3.u, chainType == "TRB") # 3187

tcrsCombSub3Tra$studyGroup = gsub("AAbNeg", "HC", tcrsCombSub3Tra$studyGroup)

## select TCR subset to use

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data") # folder containing data required to generate figures

#load all TCRs

tcrs = read.csv("combined_TCRs_IAR1_IAR2.csv", stringsAsFactors = F)

## load IMGT data all TCRs

imgt = read.delim("IMGT_all_filtered_3_Nt-sequences.txt", stringsAsFactors = F)
colnames(imgt) = gsub("JUNCTION", "ntJunction", colnames(imgt))

imgt$chainType = ifelse(imgt$V.GENE.and.allele %in% grep("TRA", imgt$V.GENE.and.allele, value = "TRUE"), "TRA",
						ifelse(imgt$V.GENE.and.allele %in% grep("TRB", imgt$V.GENE.and.allele, value = "TRUE"), "TRB", "other")) # 3907 TRA, 3694 TRB, 32 other
imgt = subset(imgt, !chainType == "other") # 7601

ids= data.frame(strsplit2(imgt$Sequence.ID, split = "_"))
ids$X1 = NULL
colnames(ids) = c("libid", "protJunction", "set", "PITmatch", "expanded", "chainType")

imgt = cbind(ids, imgt) # 7601

imgt$chainType = imgt$chainType.1 # corrects corruption string split

imgt$donor_id = tcrs$donor_id[match(imgt$libid, tcrs$libid)]

## correct TRB PIT matches

table(imgt$chainType, imgt$PITmatch)
     
#      FALSE TRUE
#  TRA  2672 1235
#  TRB  3694    0

imgtTra = subset(imgt, chainType == "TRA") # 3907
imgtTraPitT = subset(imgtTra, PITmatch == "TRUE") # 1235

imgtTrb = subset(imgt, chainType == "TRB") # 3694

imgtTrb$PITmatch = ifelse(imgtTrb$libid %in% imgtTraPitT$libid, "TRUE", "FALSE") # 1104 TRUE, 2590 FALSE

imgt = rbind(imgtTra, imgtTrb) # 7601

table(imgt$chainType, imgt$PITmatch)
     
#      FALSE TRUE
#  TRA  2672 1235
#  TRB  2590 1104

## calculate publicity

publicity =  ddply(imgt,.(protJunction,  PITmatch, expanded, chainType), plyr::summarize, sum = length(unique(donor_id))) # 6474 total, 6451 unique
publicity =  ddply(imgt,.(protJunction, chainType), plyr::summarize, sum = length(unique(donor_id))) # 6451 total
colnames(publicity)[3] = c("publicity")

convergence =  ddply(imgt,.(protJunction,  PITmatch, expanded, chainType), plyr::summarize, sum = length(unique(ntJunction))) # 6474 total, 6451 unique
convergence =  ddply(imgt,.(protJunction,  chainType), plyr::summarize, sum = length(unique(ntJunction))) # 6451

all(publicity$protJunction == convergence$protJunction)
all(publicity$chainType == convergence$chainType)
#[1] both TRUE # thus, OK to combine by column

## combine publicity and convergence data frames
toPlot = cbind(publicity, convergence$sum) # 6451
colnames(toPlot) = c("protJunction", "chainType", "publicity", "convergence")

toPlot$PITmatch = imgt$PITmatch[match(toPlot$protJunction, imgt$protJunction)] # 942 TRA, 900 TRB
toPlot$expanded = imgt$expanded[match(toPlot$protJunction, imgt$protJunction)] # 

#toPlot$PITmatch = tcrsCombSub3.u$PITmatch[match(toPlot$protJunction, tcrsCombSub3.u$junction)] # 942 TRA, 900 TRB
#toPlot$expanded = tcrsCombSub3.u$expanded[match(toPlot$protJunction, tcrsCombSub3.u$junction)] # 

toPlot = subset(toPlot, protJunction %in% tcrsCombSub3.u$junction)
toPlot = toPlot[!duplicated(toPlot$protJunction),] # 6451

#colnames(toPlot) = c("protJunction", "chainType", "publicity", "convergence")

#toPlot = tcrsCombSub3.u
#tcrsCombSub3.u$publicity = publicity$publicity[match(tcrsCombSub3.u$junction, publicity$protJunction)]

#tcrsCombSub3.u$pubT = tcrsCombSub3.u$publicity>1
#table(tcrsCombSub3.u$pubT )

## convert publicity adn convergence to binary variables
toPlot$pubT = toPlot$publicity >1
table(toPlot$pubT)

#non-Public     Public 
#      6347        104 

toPlot$privT = toPlot$publicity ==1 & toPlot$expanded == "E"
table(toPlot$privT)

#FALSE  TRUE 
# 5831   620 

toPlot$convT = toPlot$convergence >1
table(toPlot$convT)

#FALSE  TRUE 
# 6358    93 

toPlot$PITmatch = gsub("TRUE", "PIT-matched", toPlot$PITmatch)
toPlot$PITmatch = gsub("FALSE", "non-PIT-matched", toPlot$PITmatch)

toPlot$convT = gsub("TRUE", "Convergent", toPlot$convT)
toPlot$convT = gsub("FALSE", "non-Convergent", toPlot$convT)

toPlot$pubT = gsub("TRUE", "Public", toPlot$pubT)
toPlot$pubT = gsub("FALSE", "non-Public", toPlot$pubT)
toPlot$privT = gsub("TRUE", "Private", toPlot$privT)
toPlot$privT = gsub("FALSE", "non-Private", toPlot$privT)

##########################
## make plot- Convergence by chain Type

df = as.data.frame.matrix(table(toPlot$chainType, toPlot$convT))
fisher.test(df) # p-value = 4.804e-06

df$chainType = row.names(df)
xNames = data.frame(table(toPlot$convT))

mdf = melt(df)
mdf$tot = xNames$Freq[match(mdf$variable, xNames$Var1)]
mdf$Freq = mdf$value/mdf$tot

if(dev.cur()>1) dev.off()
quartz(width=8,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 6))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(mdf, aes(x = variable, y = Freq)) + geom_col(position = "stack", aes(fill = chainType), alpha = 01) + scale_fill_manual(values= cbPalette)
last_plot() + theme(axis.text.x=element_text(angle=45, hjust=1))
xlab = ""
ylab = "Fraction\n"
last_plot() + labs(x = xlab, y = ylab)
last_plot() + scale_y_continuous(limits = c(0, 1.10), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00))

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("FigS6A_convergence_by_chainType.pdf") #TRA

p = last_plot()

ggsave(filename, p)

###################### 
## make plot- Convergence by Public

toPlotSub = subset(toPlot, chainType == "TRA")

df = as.data.frame.matrix(table(toPlotSub$pubT, toPlotSub$convT))
fisher.test(df) # p-value = <2.2e-16

df$share = row.names(df)
xNames = data.frame(table(toPlotSub$convT))

mdf = melt(df)
mdf$tot = xNames$Freq[match(mdf$variable, xNames$Var1)]
mdf$Freq = mdf$value/mdf$tot

if(dev.cur()>1) dev.off()
quartz(width=8,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 6))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(mdf, aes(x = variable, y = Freq)) + geom_col(position = "stack", aes(fill = share), alpha = 01) + scale_fill_manual(values= cbPalette)
last_plot() + theme(axis.text.x=element_text(angle=45, hjust=1))
xlab = ""
ylab = "Fraction\n"
last_plot() + labs(x = xlab, y = ylab)
last_plot() + scale_y_continuous(limits = c(0, 1.10), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00))

p = last_plot()

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("FigS6B_TRA_convergence_by_Public.pdf") #TRA

ggsave(filename, p)

###################### 
## make plot- Convergence by Private

toPlotSub = subset(toPlot, chainType == "TRA")

df = as.data.frame.matrix(table(toPlotSub$privT, toPlotSub$convT))
fisher.test(df) # p-value = 0.083

df$share = row.names(df)
xNames = data.frame(table(toPlotSub$convT))

mdf = melt(df)
mdf$tot = xNames$Freq[match(mdf$variable, xNames$Var1)]
mdf$Freq = mdf$value/mdf$tot

if(dev.cur()>1) dev.off()
quartz(width=8,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 6))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(mdf, aes(x = variable, y = Freq)) + geom_col(position = "stack", aes(fill = share), alpha = 01) + scale_fill_manual(values= cbPalette)
last_plot() + theme(axis.text.x=element_text(angle=45, hjust=1))
xlab = ""
ylab = "Fraction\n"
last_plot() + labs(x = xlab, y = ylab)
last_plot() + scale_y_continuous(limits = c(0, 1.10), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00))

p = last_plot()

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("FigS6C_TRA_convergence_by_Private.pdf") #TRA

ggsave(filename, p)

###################### 
## make plot- boxplots of expanded

## determine expanded TCRs

no =  ddply(imgt,.(protJunction), plyr::summarize, sum = length(libid)) # 6446 junctions, sum(no$sum) = 7596 or nrow in imgt

toPlot$Exp = no$sum[match(toPlot$protJunction, no$protJunction)]

toPlotSub = subset(toPlot, chainType == "TRA")

df1 = as.data.frame.matrix(table(toPlotSub$pubT, toPlotSub$convT))
fisher.test(df1) # p-value = <2.2e-16

df2 = as.data.frame.matrix(table(toPlotSub$privT, toPlotSub$convT))
fisher.test(df2) # p-value = 0.083

var1 = subset(toPlotSub, pubT == "Public")
var2 = subset(toPlotSub, privT == "Private")
var1$set = c("public")
var2$set = c("private")

wilcox.test(log2(var1$Exp), log2(var2$Exp)) #  6.581e-05

toPlotSub2 = rbind(var1, var2)

library(ggsignif)

if(dev.cur()>1) dev.off()
quartz(width=6,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 6))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(toPlotSub2, aes(x = set, y = log2(Exp))) + geom_boxplot(fill = "#eeeeee", outlier.alpha = 0) #+ scale_fill_manual(values= cbPalette)
last_plot() + geom_dotplot(data = toPlotSub2, aes(x = set, y = log2(Exp)), position = "identity", method = "histodot", binaxis = "y", stackdir = "center", dotsize = 0.3)

last_plot() + theme(axis.text.x=element_text(angle=45, hjust=1))
xlab = ""
ylab = "Degree expansion,\n no. cells\n"
last_plot() + labs(x = xlab, y = ylab)
last_plot() + scale_y_continuous(limits = c(0, 6), breaks = c(0, 1, 2, 3,4, 5, 6))

last_plot() + geom_signif(comparisons = list(c("private", "public")), textsize = 20, y_position = 5, map_signif_level = T, test = "wilcox.test", test.args=list(alternative = "two.sided", var.equal = T, paired=F)) 

p = last_plot()

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("FigS6D_boxplot_public_private_expansion.pdf") #TRA

ggsave(filename, p)

###################### 
## make plot- Convergence TRA by PITmatch

toPlotSub = subset(toPlot, chainType == "TRA") # 3264

df = as.data.frame.matrix(table(toPlotSub$PITmatch, toPlotSub$convT))
fisher.test(df) # p-value = <2.2e-16

df$share = row.names(df)
xNames = data.frame(table(toPlotSub$convT))

mdf = melt(df)
mdf$tot = xNames$Freq[match(mdf$variable, xNames$Var1)]
mdf$Freq = mdf$value/mdf$tot

if(dev.cur()>1) dev.off()
quartz(width=9,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 6))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(mdf, aes(x = variable, y = Freq)) + geom_col(position = "stack", aes(fill = share), alpha = 01) + scale_fill_manual(values= cbPalette)
last_plot() + theme(axis.text.x=element_text(angle=45, hjust=1))
xlab = ""
ylab = "Fraction\n"
last_plot() + labs(x = xlab, y = ylab)
last_plot() + scale_y_continuous(limits = c(0, 1.10), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00))

p = last_plot()

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("FigS6E_TRA_convergence_by_PITmatch.pdf") #TRA

ggsave(filename, p)

###################### 
## make plot- Convergence TRB by PITmatch

toPlotSub = subset(toPlot, chainType == "TRB") # 3187

df = as.data.frame.matrix(table(toPlotSub$PITmatch, toPlotSub$convT))
fisher.test(df) # p-value = 0.01001

df$share = row.names(df)
xNames = data.frame(table(toPlotSub$convT))

mdf = melt(df)
mdf$tot = xNames$Freq[match(mdf$variable, xNames$Var1)]
mdf$Freq = mdf$value/mdf$tot

if(dev.cur()>1) dev.off()
quartz(width=9,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 6))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(mdf, aes(x = variable, y = Freq)) + geom_col(position = "stack", aes(fill = share), alpha = 01) + scale_fill_manual(values= cbPalette)
last_plot() + theme(axis.text.x=element_text(angle=45, hjust=1))
xlab = ""
ylab = "Fraction\n"
last_plot() + labs(x = xlab, y = ylab)
last_plot() + scale_y_continuous(limits = c(0, 1.10), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00))

p = last_plot()

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = paste0("FigS6F_TRB_convergence_by_PITmatch.pdf") #TRA

ggsave(filename, p)

p.adjust(c(4.80e-6,2.2e-16, 0.083, 6.581e-5, 2.2e-16, .01001))
# [1] 1.9200e-05 1.3200e-15 8.3000e-02 1.9743e-04 1.3200e-15
#[6] 2.0020e-02



