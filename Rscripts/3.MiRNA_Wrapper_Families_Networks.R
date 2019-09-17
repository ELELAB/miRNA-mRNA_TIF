# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Copyright (C) 2018, Thilde Bagger Terkelsen <thilde@cancer.dk> <thildebate@gmail.com>
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Working directory
my.wd <- "~/Desktop/Thilde/MS_MS_TIF_analysis_2014_2015/TIF_miRNAmRNA/TIF_miRNA"



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load Datasets MiRNAs
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd(paste0(my.wd,"/Data"))

# All TIF miRNAs, no normals
AllmiR <- read.delim("TIF.txt", header = TRUE)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd(paste0(my.wd,"/Results/DE_Tables"))

# TIF_NIF
TIFNIFmiR_up <- read.delim("TIF_NIF_up.txt", header = TRUE)
TIFNIFmiR_down <- read.delim("TIF_NIF_down.txt", header = TRUE)
TIFNIFmiR <- rbind(TIFNIFmiR_up, TIFNIFmiR_down)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Other vs TNBC
OtherTNBCmiR_up <-  read.delim("Other_TNBC_up.txt", header = TRUE)
OtherTNBCmiR_down <-  read.delim("Other_TNBC_down.txt", header = TRUE)
OtherTNBCmiR <- rbind(OtherTNBCmiR_up, OtherTNBCmiR_down)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Cluster 1 vs Cluster 2
C1C2miR_up <- read.delim("C1_C2_up.txt", header = TRUE)
C1C2miR_down <- read.delim("C1_C2_down.txt", header = TRUE)
C1C2miR <- rbind(C1C2miR_up, C1C2miR_down)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ER+ vs ER-
ERpERnmiR_up <- read.delim("ERp_ERn_up.txt")
ERpERnmiR_down <- read.delim("ERp_ERn_down.txt")
ERpERnmiR <- rbind(ERpERnmiR_up, ERpERnmiR_down)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# PgR+ vs PgR-
PGRpPGRnmiR_up <- read.delim("PGRp_PGRn_up.txt")
PGRpPGRnmiR_down <- read.delim("PGRp_PGRn_down.txt")
PGRpPGRnmiR <- rbind(PGRpPGRnmiR_up, PGRpPGRnmiR_down)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#  High TILs vs Low TILs
HTLTmiR_up <- read.delim("HT_LT_up.txt", header = TRUE)
HTLTmiR_down <- read.delim("HT_LT_down.txt", header = TRUE)
HTLTmiR <- rbind(HTLTmiR_up, HTLTmiR_down)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#  High Grade vs Low/Medium Grade
HGLGmiR_up <- read.delim("HG_LG_up.txt", header = TRUE) 
HGLGmiR_down <- read.delim("HG_LG_down.txt", header = TRUE)
HGLGmiR <- rbind(HGLGmiR_up, HGLGmiR_down)





# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




                                                                                                            ### MiRNA FAMILIES ###




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# TargetScan miRNA families DB
setwd(paste0(my.wd,"/Backgrounds_and_Databases/"))

# Load miRNA families
mirFam <- read.delim("miRFamTargetScan.txt", header = TRUE)

# Simplify annotation
mirFamSimple <- mirFam[,c(1,4)]
mirFamSimple$Family <- gsub("-3p|-5p|-2-3p|-1-3p|-2-3p|-1-5p|-2-5p|[.]1|[.]2", "",mirFamSimple$Family)

# List of all sets
miRNAsets <- list(rownames(AllmiR), rownames(TIFNIFmiR), rownames(OtherTNBCmiR), rownames(C1C2miR), rownames(ERpERnmiR), rownames(PGRpPGRnmiR), rownames(HTLTmiR), rownames(HGLGmiR))
miRNANames <- list("All", "TIFNIF", "Subtype", "Clusters", "ER", "PGR", "TILs", "Grade")

setwd(paste0(my.wd,"/Results/Plots"))

# For all sets map miRNAs to miRNA family
families <- map_miRFam(miRNAsets, mirFamSimple, miRNANames)





# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




                                                                                                                 ### MiRNA-GENE NETWORKS ###




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




# Load miRNA-Gene interactions
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd(paste0(my.wd,"/Backgrounds_and_Databases/"))
TargetScan <- read.delim("TargetScanHuman.txt", header = TRUE)
InBioMapDE <- read.delim("InBioMapDB.txt", header = TRUE)




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load mRNA sets
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load Datasets mRNAs


setwd("~/Desktop/Thilde/MS_MS_TIF_analysis_2014_2015/TIF_miRNAmRNA/TIF_mRNA/Results/DE_Tables")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Other vs TNBC
HER2TNBCmR_up <-  read.delim("HER2_TNBC_up.txt", header = TRUE)
HER2TNBCmR_down <-  read.delim("HER2_TNBC_down.txt", header = TRUE)
HER2TNBCmR <- rbind(HER2TNBCmR_up, HER2TNBCmR_down)

LumATNBCmR_up <-  read.delim("LumA_TNBC_up.txt", header = TRUE)
LumATNBCmR_down <-  read.delim("LumA_TNBC_down.txt", header = TRUE)
LumATNBCmR <- rbind(LumATNBCmR_up, LumATNBCmR_down)

LumBTNBCmR_up <-  read.delim("LumB_TNBC_up.txt", header = TRUE)
LumBTNBCmR_down <-  read.delim("LumB_TNBC_down.txt", header = TRUE)
LumBTNBCmR <- rbind(LumBTNBCmR_up, LumBTNBCmR_down)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Cluster 1 vs Cluster 2
C1C2mR_up <- read.delim("C1_C2_up.txt", header = TRUE)
C1C2mR_down <- read.delim("C1_C2_down.txt", header = TRUE)
C1C2mR <- rbind(C1C2mR_up, C1C2mR_down)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ER+ vs ER-
ERpERnmR_up <- read.delim("ERp_ERn_up.txt")
ERpERnmR_down <- read.delim("ERp_ERn_down.txt")
ERpERnmR <- rbind(ERpERnmR_up, ERpERnmR_down)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# High HER2 vs Low HER2
HHHLmR_up <- read.delim("HH_HL_up.txt")
HHHLmR_down <- read.delim("HH_HL_down.txt")
HHHLmR <- rbind(HHHLmR_up, HHHLmR_down)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#  High TILs vs Low TILs
HTLTmR_up <- read.delim("HT_LT_up.txt", header = TRUE)
HTLTmR_down <- read.delim("HT_LT_down.txt", header = TRUE)
HTLTmR <- rbind(HTLTmR_up, HTLTmR_down)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#  High Grade vs Low/Medium Grade
HGLGmR_up <- read.delim("HG_LG_up.txt", header = TRUE) 
HGLGmR_down <- read.delim("HG_LG_down.txt", header = TRUE)
HGLGmR <- rbind(HGLGmR_up, HGLGmR_down)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Construct Networks
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


setwd("~/Desktop/Thilde/MS_MS_TIF_analysis_2014_2015/TIF_miRNAmRNA/Joint/Tables/Network_Tables")



HER2TNBCNetwork <- DENetWork(rownames(OtherTNBCmiR_up), rownames(OtherTNBCmiR_down), HER2TNBCmR_up$GeneName, HER2TNBCmR_down$GeneName, TargetScan, InBioMapDE, FALSE)
HER2TNBCNetworkLFC <- DENetWorkFC(HER2TNBCNetwork, OtherTNBCmiR, HER2TNBCmR)


# Write out tables
#write.table(HER2TNBCNetworkLFC[[1]], "HER2TNBCNetworkLFC.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(HER2TNBCNetworkLFC[[2]], "HER2TNBCNodeInfo.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

LumATNBCNetwork <- DENetWork(rownames(OtherTNBCmiR_up), rownames(OtherTNBCmiR_down), LumATNBCmR_up$GeneName, LumATNBCmR_down$GeneName, TargetScan, InBioMapDE, FALSE)
LumATNBCNetworkLFC <- DENetWorkFC(LumATNBCNetwork, OtherTNBCmiR, LumATNBCmR)


# Write out tables
#write.table(LumATNBCNetworkLFC[[1]], "LumATNBCNetworkLFC.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(LumATNBCNetworkLFC[[2]], "LumATNBCNodeInfo.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


LumBTNBCNetwork <- DENetWork(rownames(OtherTNBCmiR_up), rownames(OtherTNBCmiR_down), LumBTNBCmR_up$GeneName, LumBTNBCmR_down$GeneName, TargetScan, InBioMapDE, FALSE)
LumBTNBCNetworkLFC <- DENetWorkFC(LumBTNBCNetwork, OtherTNBCmiR, LumBTNBCmR)


# Write out tables
#write.table(LumBTNBCNetworkLFC[[1]], "LumBTNBCNetworkLFC.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(LumBTNBCNetworkLFC[[2]], "LumBTNBCNodeInfo.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

C1C2Network <- DENetWork(rownames(C1C2miR_up), rownames(C1C2miR_down), C1C2mR_up$GeneName, C1C2mR_down$GeneName, TargetScan, InBioMapDE, FALSE)
C1C2NetworkLFC <- DENetWorkFC(C1C2Network, C1C2miR, C1C2mR)

# Write out tables
#write.table(C1C2NetworkLFC[[1]], "C1C2NetworkLFC.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(C1C2NetworkLFC[[2]], "C1C2NodeInfo.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

ERpERnNetwork <- DENetWork(rownames(ERpERnmiR_up), rownames(ERpERnmiR_down), ERpERnmR_up$GeneName, ERpERnmR_down$GeneName, TargetScan, InBioMapDE, FALSE)
ERpERnNetworkLFC <- DENetWorkFC(ERpERnNetwork, ERpERnmiR, ERpERnmR)


# Write out tables
#write.table(ERpERnNetworkLFC[[1]], "ERpERnNetworkLFC.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(ERpERnNetworkLFC[[2]], "ERpERnNodeInfo.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


HTLTNetwork <- DENetWork(rownames(HTLTmiR_up), rownames(HTLTmiR_down), HTLTmR_up$GeneName, HTLTmR_down$GeneName, TargetScan, InBioMapDE, FALSE)
HTLTNetworkLFC <- DENetWorkFC(HTLTNetwork, HTLTmiR, HTLTmR)


# Write out tables
#write.table(HTLTNetworkLFC[[1]], "HTLTNetworkLFC.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(HTLTNetworkLFC[[2]], "HTLTNodeInfo.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


HGLGNetwork <- DENetWork(rownames(HGLGmiR_up), rownames(HGLGmiR_down), HGLGmR_up$GeneName, HGLGmR_down$GeneName, TargetScan, InBioMapDE, FALSE)
HGLGNetworkLFC <- DENetWorkFC(HGLGNetwork, HGLGmiR, HGLGmR)


# Write out tables
#write.table(HGLGNetworkLFC[[1]], "HGLGNetworkLFC.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(HGLGNetworkLFC[[2]], "HGLGNodeInfo.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)





# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Network Plots
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd("~/Desktop/Thilde/MS_MS_TIF_analysis_2014_2015/TIF_miRNAmRNA/Joint/Plots/Network_Plots")
#setwd("~/Desktop/Thilde/MS_MS_TIF_analysis_2014_2015/TIF_miRNAmRNA/Joint/Plots/Network_Plots2")

Arcplot(HER2TNBCNetworkLFC[[1]], 12, 7, "HER2TNBC")
Arcplot(LumATNBCNetworkLFC[[1]], 12, 7, "LumATNBC")
Arcplot(LumBTNBCNetworkLFC[[1]], 12, 7, "LumBTNBC")
Arcplot(C1C2NetworkLFC[[1]], 22, 12, "C1C2")
Arcplot(ERpERnNetworkLFC[[1]], 10, 5, "ERpERn")
Arcplot(HTLTNetworkLFC[[1]], 14, 8, "HTLT")
Arcplot(HGLGNetworkLFC[[1]], 10, 5, "HGLG")



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Compare Networks and Perform Enrichment Analysis
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# List of all networks
list.full <- list(HER2TNBCNetworkLFC[[2]], LumATNBCNetworkLFC[[2]], LumBTNBCNetworkLFC[[2]], ERpERnNetworkLFC[[2]], HTLTNetworkLFC[[2]], HGLGNetworkLFC[[2]], C1C2NetworkLFC[[2]])
names(list.full) <- c("HER2TNBC", "LumATNBC", "LumBTNBC", "ERpERn", "HTLT", "HGLG", "C1C2")

# Extract common miRNAs and mRNAs from networks
full <- ExtractmiRmRNA(list.full)

# Extract unique miRNAs and mRNAs from networks
Ufull <- UniquemiRmRNA(list.full)


# List of networks with subtypes
list.sub <- list(HER2TNBCNetworkLFC[[2]], LumATNBCNetworkLFC[[2]], LumBTNBCNetworkLFC[[2]])
names(list.sub) <- c("HER2TNBC", "LumATNBC", "LumBTNBC")

# Extract common miRNAs and mRNAs from networks
sub <- ExtractmiRmRNA(list.sub)

# Extract unique miRNAs and mRNAs from networks
Usub <- UniquemiRmRNA(list.sub)


# List of networks with clusters and TILs
list.clus <- list(HTLTNetworkLFC[[2]], C1C2NetworkLFC[[2]])
names(list.clus) <- c("HTLT", "C1C2")

# Extract common miRNAs and mRNAs from networks
clus <- ExtractmiRmRNA(list.clus)

# Extract unique miRNAs and mRNAs from networks
Uclus <- UniquemiRmRNA(list.clus)




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



 
                                                                                                       ### PATHWAY ENRICHMENT ANALYSIS ###




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd(paste0(my.wd,"/Backgrounds_and_Databases/"))

# Background genes, all genes (mRNAs) from set
BGGenes <- as.character(unique(sort(read.delim("GeneBGforEnrich.txt", header = TRUE)$GeneName)))

# List of KEGG pathways with genes assigned to them
load("pwGenesHGNC.Rdata")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# HER2 vs TNBC
HER2TNBCpw <-  mirnaPwEnrich(pwGenesHGNC, BGGenes, Usub$mRNAs$HER2TNBC, simple = TRUE)
HER2TNBCpw <- HER2TNBCpw[HER2TNBCpw$Pval < 0.01,]

# LumA vs TNBC
LumATNBCpw <-  mirnaPwEnrich(pwGenesHGNC, BGGenes, Usub$mRNAs$LumATNBC, simple = TRUE)
LumATNBCpw <- LumATNBCpw[LumATNBCpw$FDR < 0.05,]

# LumB vs TNBC
LumBTNBCpw <-  mirnaPwEnrich(pwGenesHGNC, BGGenes, Usub$mRNAs$LumBTNBC, simple = TRUE)
LumBTNBCpw <- LumBTNBCpw[LumBTNBCpw$Pval < 0.01,]


# Common for subtypes
Subcommon <-  mirnaPwEnrich(pwGenesHGNC, BGGenes, Usub$mRNAs$Common, simple = TRUE)
Subcommon <- Subcommon[Subcommon$Pval < 0.05,]

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Cluster 1 vs Cluster 2
C1C2pw <- mirnaPwEnrich(pwGenesHGNC, BGGenes, Uclus$mRNA$C1C2, simple = TRUE)
C1C2pw <- C1C2pw[C1C2pw$FDR < 0.12,]

C1C2fullpw <- mirnaPwEnrich(pwGenesHGNC, BGGenes, clus$mRNA$C1C2, simple = TRUE)
C1C2fullpw <- C1C2fullpw[C1C2fullpw$FDR < 0.05,]

# High TILs vs Low TILs
TILspw <- mirnaPwEnrich(pwGenesHGNC, BGGenes, Uclus$mRNA$HTLT, simple = TRUE)
TILspw<- TILspw[TILspw$FDR < 0.05,]

TILsfullpw <- mirnaPwEnrich(pwGenesHGNC, BGGenes, clus$mRNA$HTLT, simple = TRUE)
TILsfullpw <- TILsfullpw[TILsfullpw$FDR < 0.05,]

# Common clusters and TILs
C1C2TILsCommon <- mirnaPwEnrich(pwGenesHGNC, BGGenes, Uclus$mRNA$Common, simple = TRUE)
C1C2TILsCommon <- C1C2TILsCommon[C1C2TILsCommon$FDR < 0.05,]
