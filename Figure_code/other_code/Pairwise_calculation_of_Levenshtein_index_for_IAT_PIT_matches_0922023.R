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
library(foreach)

setwd("/Users/peterlinsley/Desktop/PIT_TCR_paper_code/data") # 

## load HLA info

hla = read.xls("subject_char_w_HLA.xlsx", stringsAsFactors = F)

hla1 = subset(hla, DRB1 %in% grep("0401", hla$DRB1, value = T))
hla2 = subset(hla, DRB1 %in% grep("03", hla$DRB1, value = T))


## load public and private TCRs

KfTcrs = read.csv("201512_TCR_MasterList_w_CloneIDs.csv", stringsAsFactors = F) # 5729
colnames(KfTcrs) = gsub("tcrGraph_sharing_level", "sharing_level", colnames(KfTcrs))

KfTcrs$chainType = ifelse(KfTcrs$v_gene %in% grep("TRA", KfTcrs$v_gene, value = T), "TRA",
                          ifelse(KfTcrs$v_gene %in% grep("TRB", KfTcrs$v_gene, value = T), "TRB","other"))

kfSub1 = subset(KfTcrs, !chainType == "other") # 4124
kfSub1 = subset(kfSub1, CD45RA == "FALSE")
kfSub2 = subset(kfSub1, sharing_level %in% c("public", "private")) # 1412
iNkt1 = subset(KfTcrs, junction == "CVVSDRGSTLGRLYF")
iNkt2 = subset(KfTcrs, libid %in% iNkt1$libid)
#mait1 = subset(KfTcrs, junction == "CAVKDSNYQLIW" )
mait1 = subset(KfTcrs, v_gene == "TRAV1-2" & (j_gene == "TRAJ33" | j_gene == "TRAJ20" | j_gene == "TRAJ12")) # 12
mait2 = subset(KfTcrs, libid %in% mait1$libid) # 30

kfSub2 = subset(kfSub2, !junction %in% iNkt2$junction) # 1363 by juncton 1365 by libid
kfSub3 = subset(kfSub1, !junction %in% iNkt2$junction) # 4165
kfSub3 = subset(kfSub3, !junction %in% mait2$junction) # 4135

kfSub2.u = kfSub2[!duplicated(kfSub2$junction),] # 455
kfSub3.u = kfSub3[!duplicated(kfSub3$junction),] # 3168

chainsPerLib = ddply(kfSub3.u,.(libid), plyr::summarize, sum = length(junction)) # 1731
chainsPerLib = subset(chainsPerLib, sum == 2) # 1140

clonesTra = subset(kfSub3.u, chainType == "TRA" & libid %in% chainsPerLib$libid) # 1144
clonesTrb = subset(kfSub3.u, chainType == "TRB" & libid %in% chainsPerLib$libid) # 1136

clonesTra = subset(clonesTra, libid %in% clonesTrb$libid) # 1114
clonesTrb = subset(clonesTrb, libid %in% clonesTra$libid) # 1114

length(setdiff(clonesTra$libid, clonesTrb$libid))
#[1] 0
length(setdiff(clonesTrb$libid, clonesTra$libid))
#[1] 0

clones = merge(clonesTra, clonesTrb, by = "libid") # 1114
clones$alphaChain = paste(clones$v_gene.x, clones$junction.x, clones$j_gene.x, sep = "_")
clones$betaChain = paste(clones$v_gene.y, clones$junction.y, clones$j_gene.y, sep = "_")

noTra =  ddply(clones,.(libid), plyr::summarize, sum = length(junction.x)) # 1114
noTrb =  ddply(clones,.(libid), plyr::summarize, sum = length(junction.y)) # 1114

clones.u = clones[!duplicated(clones$alphaChain),] # 1114

## Su et al TCRs

su = read.csv("Su_et_al_filtered_TCRs.csv", stringsAsFactors = F) # 326782
su$X = NULL

suNkt1 = subset(su, junction == "CVVSDRGSTLGRLYF") #45
suNkt2 = subset(su, barcode %in% suNkt1$barcode) # 90
sumait1 = subset(su, v_gene == "TRAV1-2" & (j_gene == "TRAJ33" | j_gene == "TRAJ20" | j_gene == "TRAJ12")) # 2533
sumait2 = subset(su, barcode %in% sumait1$barcode) # 5106

su = subset(su, !junction %in% suNkt2$junction & !junction %in% sumait2$junction) # 320372

suCd4 = subset(su, su$cellType == "cd4")
suCd8 = subset(su, su$cellType == "cd8")

su.u = su[!duplicated(su$junction),]
suHc = subset(su.u, group == "HC")
suCovid = subset(su.u, group == "COVID")

myTcrs = kfSub2

all = kfSub2

all$study_group = gsub("early onset T1D", "newT1D", all$study_group)

allSub = all

all.HC = subset(all, study_group == "HC") # 259
all.newT1D = subset(all, study_group == "newT1D") # 497
all.T1D = subset(all, study_group == "T1D") # 597

all.TRA = subset(all, chainType == "TRA")

## get islet TCR seqs from Maki Nakayama, MAKI.NAKAYAMA@CUANSCHUTZ.EDU

mn = read.xls("220519_scPCR_Illumina.xlsx", stringsAsFactors = F)
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

mn.u = mn[!duplicated(mn$alphaChain),]

length(intersect(mnTcrs.u$junction, kfSub3.u$junction)) # 50/3181 = 0.01571833
length(intersect(mnTcrs.u$junction, suHc$junction)) # 185/25153 = 0.007354987
length(intersect(mnTcrs.u$junction, suCovid$junction)) # 1342/198753 = 0.006752099

## load Mitchell et al TRB junctions that overlap with known INS TCR TRB chains (JCI Insight. 2022;7(18):e161885). https://doi.org/10.1172/jci.insight.161885.

ins = read.delim("Mitchell_Michels_insulin_TRB_TCRs.txt", stringsAsFactors = F, header = F)

## load 0,1 TRA mismatches

#setwd("/mnt/bioinformatics/workspace/plinsley/IAR_T_cell_TCRs") # 
#mMatch = read.csv("Mismatcvh_0_1_IARvsIslet_TCRs.csv", stringsAsFactors = F) # 1389

#randMm = read.csv("Mismatcvh_0_1_randomSuvsIslet_TCRs.csv", stringsAsFactors = F) # 1389

## load all random matches

#randAll = read.csv("Islet_matches_randomSuvsIslet_TCRs.csv")

## find TRB chains of mismatches

#IarMmSub = subset(kfSub3.u, junction %in% mMatch$aJunc1)
#IarMmSub2 = subset(kfSub3.u, libid %in% IarMmSub$libid)
#IarTrbSub3 = subset(IarMmSub2, chainType == "TRB")

## find TRB chains of random matches

#randHcSub = subset(suHc, suHc$junction %in% sample(suHc$junction, size = nrow(IarMmSub2)))
#randHcTrb = subset(randHcSub, chainType == "TRB")
#randHcTrb = subset(randHcTrb, randHcTrb$junction %in% sample(randHcTrb$junction, size = nrow(IarTrbSub3)))

## find TRB chains of random matches from non-matched 

#notMmSub = subset(kfSub3.u, !junction %in% mMatch$aJunc1)
#notMmSub2 = subset(kfSub3.u, libid %in% notMmSub$libid)
#notMmSub3 = subset(IarMmSub2, chainType == "TRB")
#notMmSub4 = subset(notMmSub3, junction %in% sample(notMmSub3$junction, size = length(IarTrbSub3$junction)))

## add Fari's antigen specificity

ag = read.csv("TCRs_specificity_accounting.csv", stringsAsFactors = F)
ag = subset(ag, !primary_specificity == "glucosylseramide") # remove iNKT

## overlap of unique junctions

var1 = mnTcrs.u$junction %in% intersect(mnTcrs$junction, kfSub3$junction)
table(var1)

## Levenshtein distances in AgSp versus islet overlapping TRA chains

start = Sys.time()

xnames = unique(all.TRA$junction)
ynames = mn.u$Junction.2
egs <- expand.grid(xnames, ynames, stringsAsFactors=FALSE)
ncx = nrow(egs)
ncx = 100 # for testing

DF = (matrix(nrow = 2479188, ncol = 4))
colnames(DF) = c("index",  "aJunc1", "aJunc2", "lv")
timeDf = matrix(nrow = 2479188, ncol = 2)

foreach(i = 1:ncx) %do% {
  r  = egs[i, 1]
  s  = egs[i, 2]
  d = stringdistmatrix(a = r, b = s, method = "lv", useNames = F, nthread =20)
  timeResult = c(i, Sys.time()-start)
  timeDf[i, ] = timeResult
  if (as.numeric(d)>20) {
  } else {
    result = as.character(c(i, r, s, d))
    DF[i, ] = result
  }
}

DF = data.frame(DF)
DF = na.omit(DF)
DF$lv = as.numeric(DF$lv)
filename1 = "Pairwise_Levenshtein_index_calcs.csv"
write.csv(DF, filename1)

timeDf = na.omit(data.frame(timeDf))

#quartz(width=10,height=8, dpi=72)  ### open plotting window
update_geom_defaults("line", aes(size = 2))
update_geom_defaults("point", aes(size = 2))

theme_set(theme_bw(36) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.key = element_blank()))

cbPalette = c('#66c2a5','#fc8d62', 'gray')

ggplot(timeDf, aes(x = X1/1000, y = X2/60)) + geom_point()

end = Sys.time()
elapsed = end-start; elapsed # 1.35953 hours

