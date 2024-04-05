rm(list = ls())
library(plyr)
library(edgeR)
library(compareGroups)
library(dbscan)
library(MAST)
library(ggsignif)
library(monocle3)
library(dplyr)
library(geosphere)
library(RColorBrewer)
library(gtools)
library(psych)
library(gdata)
library(ggplot2); library(reshape2); theme_set(theme_bw(26) + theme(panel.grid.major = element_blank(), 
                                                 panel.grid.minor = element_blank()) +
                              theme(legend.key = element_blank()))
update_geom_defaults("point", aes(size = 4))

##################
## TCR CDR peptide contacts- all 

## Prepare and load data file
## Use all PIT-matched IAR1 TCRs with known specificity. Construct models using TCRmodel (https://tcrmodel.ibbr.umd.edu/). Assumed ##same peptide bound as by the exmplar with known specificity. View models in Mol* viewer (https://molstar.org/viewer/). Select peptide chain, the select residues, followed by Manipulate Selection, Surrounding Residues 5A of selection. Record number of residues highlighted for each chain in sequence viewer- repeat counts at least twice or until getting a consensus result. Identify CDR1, CDR2, CDR3 and other regions using IMGT. Manually tabulate results in CDR regions. Verify results by checking at least twice.

## load tabulated and annotated Structural results. 

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data/")

contacts = read.delim("PIT_matched_TCRs_w_chain_seqs_manual_All_080123.txt")

contacts$traR13 = contacts$chain1.cdr1/contacts$chain.1.cdr3
contacts$trbR13 = contacts$chain2.cdr1/contacts$chain.2.cdr3

mCont = melt(contacts, measure.vars = c("chain1.cdr1",	"chain.1.cdr2", "chain.1.cdr3", "chain2.cdr1",	"chain.2.cdr2",	"chain.2.cdr3", "otherTRA",	"otherTRB"))

colnames(mCont) = gsub("value","contacts", colnames(mCont))
mCont$variable = gsub("chain1.cdr1","TraCDR1", mCont$variable, fixed = T)
mCont$variable = gsub("chain.1.cdr1","TraCDR2", mCont$variable, fixed = T)
mCont$variable = gsub("chain.1.cdr2","TraCDR2", mCont$variable, fixed = T)

mCont$variable = gsub("chain2.cdr1","TrbCDR1", mCont$variable, fixed = T)
mCont$variable = gsub("chain.2.cdr1","TrbCDR1", mCont$variable, fixed = T)
mCont$variable = gsub("chain.2.cdr2","TrbCDR2", mCont$variable, fixed = T)
mCont$variable = gsub("chain.2.cdr3","TrbCDR3", mCont$variable, fixed = T)
mCont$variable = gsub("otherTRA","TraOther", mCont$variable, fixed = T)
mCont$variable = gsub("otherTRB","TrbOther", mCont$variable, fixed = T)

mCont$variable = gsub("chain.1.cdr3","TraCDR3", mCont$variable, fixed = T)

#mCont$variable = factor(mCont$variable, levels = )
mCont$chain = ifelse(mCont$variable %in% grep("Tra", mCont$variable, value = T), "TRA", 
				ifelse(mCont$variable %in% grep("Trb", mCont$variable, value = T), "TRB","other"))
				
mCont$PITmatch = gsub("TRUE", "PIT-matched", mCont$PITmatch)
mCont$PITmatch = gsub("FALSE", "non-PIT-matched", mCont$PITmatch)

toPlot2 = subset(mCont, !variable %in% grep("other", mCont$variable, value = T))
toPlot2 = subset(toPlot2, !variable %in% grep("cdr2", mCont$variable, value = T))
toPlot2 = na.omit(toPlot2)

###################

## facet plot

if(dev.cur()>1)dev.off()
quartz(width=16,height=8, dpi=72)  ### open plotting window

theme_set(theme_bw(32) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))
update_geom_defaults("point", aes(size = 4))
update_geom_defaults("density2d", aes(size = 4))

#pal = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffed6f','#b15928', "gray")	
pal = c("#1b9e77", "#d95f02", "#7570b3")
#pal = c("#d3d3d3", "#1b9e77", "other") # all dots gray palette

ggplot(mCont, aes( y = contacts, x = traLen)) + geom_jitter(aes(colour = PITmatch, group = variable), width = 0.5, height = 0.5) + facet_wrap(~variable, ncol = 4, scales = "free_y")
last_plot() + geom_smooth(method='lm',formula=y ~ x, aes(), colour = "black");
last_plot() + scale_color_manual(values = pal, name = "PITmatch")# custom colors

xlab = "\nTRA junction length (AA)"
ylab = "Peptide contacts\n"
last_plot() + labs(x = xlab, y = ylab)
#last_plot() + scale_y_continuous(limits = c(0, 13), breaks = c(0, 5, 10))

p = last_plot()
setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/Figure_PDFs/")
filename = "FigS10_peptide_contacts_by_CDRs_5A.pdf"
ggsave(filename, p)

## linear model

toModel = mCont

var1 = subset(toModel, variable == "contactsCDR1A")
var2 = subset(toModel, variable == "contactsCDR2A")
var3 = subset(toModel, variable == "contactsCDR3A")
var4 = subset(toModel, variable == "contacts5JA")
var5 = subset(toModel, variable == "contactsOthersA")

var6 = subset(toModel, variable == "contactsCDR1B")
var7 = subset(toModel, variable == "contactsCDR2B")
var8 = subset(toModel, variable == "contactsCDR3B")
var9 = subset(toModel, variable == "contacts5JB")
var10 = subset(toModel, variable == "contactsOthersB")

DF = data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors = F)
xnames = unique(toModel$variable)
ncx = length(xnames)

## model

for(i in 1:ncx){
	q = as.character(xnames[i])
	qsub = subset(mCont, variable == q)
	m = lm(contacts~traLen, data = qsub)
	estimate = summary(m)$coefficients[2]
	pval = summary(m)$coefficients[8]
	result = c(q, estimate, pval)
	DF[i, ] = result
	
	}
colnames(DF) = c("variable", "estimate", "pval")
DF$pAdj = p.adjust(DF$pval)

DF
#  variable            estimate                 pval        pAdj
#1  TraCDR1               -0.25     0.12512485406065 0.750749124
#2  TraCDR2 -0.0326086956521738    0.223666500707135 0.894666003
#3  TraCDR3   0.532608695652174 0.000105305495491008 0.000842444
#4  TrbCDR1 -0.0652173913043473    0.503990726222478 1.000000000
#5  TrbCDR2  -0.032608695652174    0.782835382730538 1.000000000
#6  TrbCDR3  -0.271739130434783    0.149431260262739 0.750749124
#7 TraOther  0.0217391304347826    0.662519556328938 1.000000000
#8 TrbOther  -0.130434782608695   0.0134871534032227 0.094410074


