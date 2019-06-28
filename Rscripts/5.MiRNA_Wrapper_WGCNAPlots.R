# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Copyright (C) 2018, Thilde Bagger Terkelsen <thilde@cancer.dk> <thildebate@gmail.com>
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Working directory

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
my.wd <- ""
setwd(paste0(my.wd, "/DE_Tables/mRNA_Tables"))


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load mRNA Sets
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




                                                                                                      ### WGCNA CO-EXPRESSED GENES IN MODULES ###




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load Modules
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load("genes_modules_colors.RData")


# subset by colors of interest
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Blue <- genes_colors[genes_colors$mergedColors == "blue",]
Yellow <- genes_colors[genes_colors$mergedColors == "yellow",]
Lightgreen <- genes_colors[genes_colors$mergedColors == "lightgreen",]
Red <- genes_colors[genes_colors$mergedColors == "red",]
Green <- genes_colors[genes_colors$mergedColors == "green",]
Grey60 <-  genes_colors[genes_colors$mergedColors == "grey60",]



# Extract genes from modules
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TILsModules <- table(genes_colors[genes_colors$V1 %in% HTLTmR$GeneName, ]$mergedColors)
ClusModules <- table(genes_colors[genes_colors$V1 %in% C1C2mR$GeneName, ]$mergedColors)
GradeModules <- table(genes_colors[genes_colors$V1 %in% HGLGmR$GeneName, ]$mergedColors)
ERModules <- table(genes_colors[genes_colors$V1 %in% ERpERnmR$GeneName, ]$mergedColors)
HER2TNBCModules <- table(genes_colors[genes_colors$V1 %in% HER2TNBCmR$GeneName, ]$mergedColors)
LumATNBCModules <- table(genes_colors[genes_colors$V1 %in% LumATNBCmR$GeneName, ]$mergedColors)
LumBTNBCModules <- table(genes_colors[genes_colors$V1 %in% LumBTNBCmR$GeneName, ]$mergedColors)



# Venn Diagrams quantify the genes in modules 
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pdf("vennTILs.pdf")
vennTILs <- venn.diagram(list(A=HTLTmR$GeneName, B=Yellow$V1, C=Lightgreen$V1), category.names = c("DE TILs", "Yellow Module","LightGreen Module"), filename=NULL, lwd = 0, cat.pos=0, sub.cex = 2, cat.cex= 1.5, cex=1.5, fill=c("#FFF4EC", "#DDEA7C", "#91F5AD"))
grid.draw(vennTILs)
dev.off()

pdf("vennGR.pdf")
vennGR <- venn.diagram(list(A=HGLGmR$GeneName, B=Blue$V1, C=Green$V1), category.names = c("DE Grade", "Blue Module", "Green Module"), filename=NULL, lwd = 0, cat.pos=0, sub.cex = 2, cat.cex= 1.5, cex=1.5, fill=c("#FFF4EC", "#445E93", "#4FB286"))
grid.draw(vennGR)
dev.off()

pdf("vennClus.pdf")
vennClus <- venn.diagram(list(A=C1C2mR$GeneName, B=Yellow$V1, C=Green$V1), category.names = c("DE Clusters", "Yellow Module", "Green Module"), filename=NULL, lwd = 0, cat.pos=0, sub.cex = 2, cat.cex= 1.5, cex=1.5, fill=c("#FFF4EC", "#DDEA7C", "#4FB286"))
grid.draw(vennClus)
dev.off()

pdf("vennER.pdf")
vennER <- venn.diagram(list(A=ERpERnmR$GeneName, B=Green$V1), category.names = c("DE Estrogen Receptor", "Green Module"), filename=NULL, lwd = 0, cat.pos=0, sub.cex = 2, cat.cex= 1.5, cex=1.5, fill=c("#FFF4EC", "#4FB286"))
grid.draw(vennER)
dev.off()

pdf("vennHER2.pdf")
vennHER2 <- venn.diagram(list(A=HER2TNBCmR$GeneName, B=Grey60$V1), category.names = c("DE Her2 vs TNBC", "Grey60 Module"), filename=NULL, lwd = 0, cat.pos=0, sub.cex = 2, cat.cex= 1.5, cex=1.5, fill=c("#FFF4EC", "grey60"))
grid.draw(vennHER2)
dev.off()

pdf("vennLumA.pdf")
vennLumA <- venn.diagram(list(A=LumATNBCmR$GeneName, B=Grey60$V1, C=Green$V1, D=Red$V1), category.names = c("DE LumA vs TNBC", "Grey60 Module", "Green Module", "Red Module"), filename=NULL, lwd = 0, cat.pos=0, sub.cex = 2, cat.cex= 1.5, cex=1.5, fill=c("#FFF4EC", "grey60", "#4FB286", "#D35269"))
grid.draw(vennLumA)
dev.off()

pdf("vennLumB.pdf")
vennLumB <- venn.diagram(list(A=LumBTNBCmR$GeneName, B=Grey60$V1, C=Green$V1, D=Red$V1), category.names = c("DE LumB vs TNBC", "Grey60 Module", "Green Module", "Red Module"), filename=NULL, lwd = 0, cat.pos=0, sub.cex = 2, cat.cex= 1.5, cex=1.5, fill=c("#FFF4EC", "grey60", "#4FB286", "#D35269"))
grid.draw(vennLumB)
dev.off()



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extract Genesets of Interest Across Modules Correlated with Clinical Variables.
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

TILsWGCNA <- intersect(HTLTmR$GeneName, c(as.character(Yellow$V1), as.character(Lightgreen$V1)))
GRWGCNA <- intersect(HGLGmR$GeneName, c(as.character(Blue$V1), as.character(Green$V1)))
ClusWGCNA <- intersect(C1C2mR$GeneName, c(as.character(Yellow$V1), as.character(Green$V1)))
ERWGCNA <- intersect(ERpERnmR$GeneName, as.character(Green$V1))
HER2TNBCWGCNA <- intersect(HER2TNBCmR$GeneName, as.character(Grey60$V1))
LumATNBCWGCNA <- intersect(LumATNBCmR$GeneName, c(as.character(Grey60$V1),as.character(Green$V1), as.character(Red$V1)))
LumBTNBCWGCNA <- intersect(LumBTNBCmR$GeneName, c(as.character(Grey60$V1),as.character(Green$V1), as.character(Red$V1)))
Subtypes <- unique(sort(c(HER2TNBCWGCNA, LumATNBCWGCNA,LumBTNBCWGCNA )))


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Make Heatmaps 
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd(paste0(my.wd, "/Data"))

# Load full genes expression dataframe.
load("yfilt.avg.Rdata")


# Extract genes of interest
TILsWGCNA <- yfilt.avg$E[rownames(yfilt.avg$E) %in% TILsWGCNA,]
GRWGCNA <- yfilt.avg$E[rownames(yfilt.avg$E) %in% GRWGCNA,]
#ClusWGCNA <- yfilt.avg$E[rownames(yfilt.avg$E) %in% ClusWGCNA,]
ERWGCNA <- yfilt.avg$E[rownames(yfilt.avg$E) %in% ERWGCNA,]
Subtypes <- yfilt.avg$E[rownames(yfilt.avg$E) %in% Subtypes,]


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Set color scheme for colorbar
my.ER <- get_colors(yfilt.avg$targets$ER, c("#96C3CE","#CC978E"))
my.GR <- get_colors(yfilt.avg$targets$GR, c("#0071AA","#0071AA","#032A63"))
my.TILS <- get_colors(yfilt.avg$targets$TILs, c("#D9EAFC","#AECFE8","#0071AA","#032A63"))
my.TS <- get_colors(yfilt.avg$targets$Tumor_subtype_corrected_2015_11_20, c("#93C0C1", "#C6FFE8", "#C6FFE8", "#C6FFE8", "#C6FFE8", "#3B8476"))

spacer <- as.matrix(replicate(ncol(yfilt.avg$E), "white"))
ER.cols <- cbind(spacer, my.ER)
GR.cols <- cbind(spacer, my.GR)
TILS.cols <- cbind(spacer, my.TILS)
TS.cols <- cbind(spacer, my.TS)

# Heat colors
heat.cols <- viridis(option = "magma", n=20, direction = -1)




#pdf("BestCandGrade.pdf", height = 6, width = 8)
heatmap.plus(as.matrix(GRWGCNA), col=heat.cols, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", Rowv = NA, labCol="", ColSideColors=GR.cols, margins = c(14,8), cexCol=1.2, cexRow = 1)
#dev.off()

#pdf("BestCandTILs.pdf", height = 19, width = 8)
heatmap.plus(as.matrix(TILsWGCNA), col=heat.cols, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", Rowv = NA, labCol="", ColSideColors=TILS.cols, margins = c(14,8), cexCol=1.2, cexRow = 0.9)
#dev.off()

#pdf("BestCandER.pdf", height = 8, width = 8)
heatmap.plus(as.matrix(ERWGCNA), col=heat.cols, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", Rowv = NA, labCol="", ColSideColors=ER.cols, margins = c(14,8), cexCol=1.2, cexRow = 1)
#dev.off()

#pdf("BestCandTNBC.pdf",height = 10, width = 8)
heatmap.plus(as.matrix(scale(Subtypes, scale = FALSE)), col=heat.cols, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", Rowv = NA, labCol="", ColSideColors=TS.cols, margins = c(14,8), cexCol=1.2, cexRow = 0.9)
#dev.off()




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load Networks
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd(paste0(my.wd, "/Network_Tables"))


HER2TNBCNetworkLFC <- read.delim("HER2TNBCNetworkLFC.txt", header = TRUE)
LumATNBCNetworkLFC <- read.delim("LumATNBCNetworkLFC.txt", header = TRUE)
LumBTNBCNetworkLFC <- read.delim("LumBTNBCNetworkLFC.txt", header = TRUE)
C1C2NetworkLFC <- read.delim("C1C2NetworkLFC.txt", header = TRUE)
ERpERnNetworkLFC <- read.delim("ERpERnNetworkLFC.txt", header = TRUE)
HTLTNetworkLFC <- read.delim("HTLTNetworkLFC.txt", header = TRUE)
HGLGNetworkLFC <- read.delim("HGLGNetworkLFC.txt", header = TRUE)


# Extract best gene-miRNA pairs
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

TILsWGCNAmiR <- BestPairs(HTLTNetworkLFC, TILsWGCNA)
#write.table(TILsWGCNAmiR, "TILsWGCNAmiR.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

GRWGCNAmiR <- BestPairs(HGLGNetworkLFC, GRWGCNA)
#write.table(GRWGCNAmiR, "GRWGCNAmiR.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

ClusWGCNAmiR <- BestPairs(C1C2NetworkLFC, ClusWGCNA)
#write.table(ClusWGCNAmiR, "ClusWGCNAmiR.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

LumATNBCWGCNAmiR <- BestPairs(LumATNBCNetworkLFC, LumATNBCWGCNA)
#write.table(LumATNBCWGCNAmiR, "LumATNBCWGCNAmiR.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#HER2TNBCWGCNAmiR <- BestPairs(HER2TNBCNetworkLFC, HER2TNBCWGCNA)
#write.table(HER2TNBCWGCNAmiR, "HER2TNBCWGCNAmiR.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

LumBTNBCWGCNAmiR <- BestPairs(LumBTNBCNetworkLFC, LumBTNBCWGCNA)
#write.table(LumBTNBCWGCNAmiR, "LumBTNBCWGCNAmiR.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Make Summary Plots
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load dataframe with best candidates
BestCand <- read.delim("BestMiRInfoConcatnatedSmall.txt", header = TRUE)

# Concatnate miRNA name with miRNA family
BestCand$miRNA <- paste0(BestCand$miRNA, " (", BestCand$Family, ")")

InformationPlot(4, BestCand, c("green", "red"), "Luminal and HER2 vs TNBC", 4, 13)
InformationPlot(5, BestCand, c("green"), "ER+ vs ER-", 4, 11)
InformationPlot(8, BestCand, c("green",  "yellow"), "High TILs vs Low TILs and High Gr vs Low Gr", 10, 14)
InformationPlot(9, BestCand, c("green"), "High Gr vs Low Gr",4, 11)
InformationPlot(10, BestCand, c("green", "yellow"), "Metastasis", 6, 12)
InformationPlot(11, BestCand, c("green", "red"), "Poor outcome", 3, 11)











