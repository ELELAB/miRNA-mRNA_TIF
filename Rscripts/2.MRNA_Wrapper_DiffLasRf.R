# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Copyright (C) 2018, Thilde Bagger Terkelsen <thilde@cancer.dk> <thildebate@gmail.com>
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# Working Directory and Load/Read
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
my.wd <- "~/Desktop/Thilde/MS_MS_TIF_analysis_2014_2015/TIF_miRNAmRNA/TIF_mRNA"
setwd(paste0(my.wd,"/Data/mRNARawFiles"))

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read datafiles
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
targets <- readTargets("MergedTarFile.txt")
x <- read.maimages(targets, source="agilent", green.only=TRUE, other.columns="gIsWellAboveBG")



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Normalization and Removal of Controls and Outliers
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Background correction
y <- backgroundCorrect(x, method="normexp")
y <- normalizeBetweenArrays(y, method="quantile")

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Remove controls and genes not expressed above background
Control <- y$genes$ControlType!=0
IsExpr <- rowSums(y$other$gIsWellAboveBG > 0) >= 9

yfilt <- y[!Control & IsExpr, ]

# Low tumour percentage, Replicates, Missing Info on TILs and Grade
LowTP <- y$targets$tp <= 30
LowTP <- which(ifelse(is.na(LowTP), FALSE, LowTP))
NoGR <- which(is.na(y$targets$GR))
NoTILS <- which(is.na(y$targets$TILs))
IsApo <- which(y$targets$Tumor_subtype_corrected_2015_11_20 == "Apocrine")
IsOther <- which(y$targets$ID %in% c("TIF54", "TIF106", "TIF290", "TIF109.1", "TIF123.1", "TIF123.2", "TIF46.1", "TIF49.1", "TIFkontroll", "TIFkontroll.1")) 
rmsamp <- unique(sort(c(LowTP,NoGR,NoTILS,IsApo,IsOther)))

yfilt <- yfilt[,-rmsamp]

yfilt.avg <- avereps(yfilt, ID=yfilt$genes$GeneName)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Map missing ensembl transcript IDs

yfilt.avg <- MapENST(yfilt.avg)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Factor Vectors for DEA, arrays and grade/immuno information
Array <- as.factor(as.character(yfilt.avg$targets$Array))
TILS <- factor(as.character(yfilt.avg$targets$TILs), levels=c("T3", "T2", "T1", "T0"))
TILSS <- as.factor(ifelse(as.character(yfilt.avg$targets$TILs) %in% c("T0", "T1"), "LT", "HT"))
GR <- factor(paste0("G", as.character(yfilt.avg$targets$GR)), levels=c("G3", "G2", "G1"))
GRS <- as.factor(ifelse(paste0("G", as.character(yfilt.avg$targets$GR)) %in% c("G1", "G2"), "LG", "HG"))
ER <- factor(as.character(yfilt.avg$targets$ER), levels = c("ERp", "ERm"))
PGR <- factor(as.character(yfilt.avg$targets$PGR),levels = c("PGRp", "PGRm"))
AR <- factor(as.character(yfilt.avg$targets$AR),levels = c("ARp", "ARRm"))
HER2 <- factor(as.character(yfilt.avg$targets$HER2),levels = c("H3","H2","H1","H0"))

TS <- as.factor(as.character(yfilt.avg$targets$Tumor_subtype_corrected_2015_11_20))
TS <- factor(ifelse(TS == "Luminal", "LumA", as.character(TS)), levels = c("LumA", "LumB", "LumB_HER2_enriched", "HER2", "TNBC"))
TS.cols <- c("#04724D","#BDCC9F","#D2AB99","#DD614A","#310A31")





# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Combat Correction for MDS Plotting
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# Combat correction for array with model design (tumour subtype information)
mod_design <-  model.matrix(~TS)
mrna_combat_TS <- ComBat(dat=yfilt.avg$E, mod= mod_design, batch = Array)
colnames(mrna_combat_TS) <- yfilt.avg$targets$ID

# Combat correction for array NO model design 
mrna_combat <- ComBat(dat=yfilt.avg$E, batch = Array)
colnames(mrna_combat) <- yfilt.avg$targets$ID


# MDS plotting
myMDSplot(yfilt.avg$E, TS, "", TS.cols)
myMDSplot(mrna_combat_ST, TS, "", TS.cols)
myMDSplot(mrna_combat, TS, "", TS.cols)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Clustering Analysis
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Optimal number of clusters, this is a slow process
#optimal_nc(mrna_combat)

# Kmeans with 2 clusters
set.seed(10)
K2.clus <- kmeans(t(mrna_combat), 2, iter.max = 10000)
clusters2 <- ifelse(K2.clus$cluster == 1, 3, K2.clus$cluster)
clusters2 <- ifelse(clusters2 == 2, 1, 2)
CLUS <- as.factor(paste0("C",clusters2))

# Kmeans with 3 clusters
set.seed(10)
K3.clus <- kmeans(t(mrna_combat), 3,iter.max = 10000)
clusters3 <- as.factor(as.character(K3.clus$cluster))


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Dendogram of MRNA Samples
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Dendogram
dend <- mrna_combat
#colnames(dend) <- yfilt.avg$targets$ID
dend <- as.dendrogram(hclust(dist(t(dend)), method = "ward.D2"))

# Colors
my.CLUS2 <- get_colors(clusters2, c("#B2B09B","#43AA8B"))
my.CLUS3 <- get_colors(clusters3,  c("#B2B09B", "#B2B09B", "#43AA8B"))
my.ER <- get_colors(ER, c("#FFFCF7","grey60"))
my.PGR <- get_colors(PGR, c("#FFFCF7","grey60"))
my.GR <- get_colors(GR, c("#0071AA","#0071AA","#032A63"))
my.TILS <- get_colors(TILS, c("#F7D76F" ,"#EFCA4F", "#EDAB49" ,"#EE964B"))
my.TS <- get_colors(TS, c("#FFAAA3", "#6F73E2", "#BFDDFF", "#BFDDFF", "#FF5465"))



# Plot dendogram

#setwd(paste0(my.wd,"/Results/Plots"))
#pdf("Dendogram_mRNA.pdf", height = 16, width = 16)

par(cex=0.6, mar = c(8,3,3,8))
nodePar <- list(lab.cex = 1.5, pch = c(NA, 19), col = "black")
plot(dend, horiz = TRUE, nodePar = nodePar)
colored_bars(cbind(my.TS, my.ER, my.PGR, my.TILS, my.GR, my.CLUS2), dend, rowLabels = c("TS","ER", "PGR", "TILs", "GR", "CLUST"), horiz = TRUE, cex.rowLabels = 1.1)

#dev.off()


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Scatter Plot and FVIZ Cluster Plot
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# FVIZ cluster plot
fviz_cluster(K2.clus, t(mrna_combat), palette = c("#43AA8B","#B2B09B"), ggtheme = theme_minimal())

# Distances, corrected for patient effects
d<-dist(t(mrna_combat))
fit <- cmdscale(d,eig=TRUE, k=2)
res<-data.frame(M1=fit$points[,1],M2=fit$points[,2])


# Scatterplot
scatterplot3d(res, pch = 16, cex.symbols = 2, color = my.CLUS2, grid=TRUE, angle = 120)






# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



                                                                                                ### Differential Expression Analysis ###



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differnetial Expression Analysis with Subtypes
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


mod_design <-  model.matrix(~0+TS+Array+TILSS)

# Create all combinations of BC subtypes for LIMMA DA analysis
combinations <- data.frame(t(combn(paste0("TS", levels(TS)), 2)))
combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")

# Making group contrasts 
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))

# Apply DE_limma function to all comparisons
my.mRNAsST <- DE_mRNA_apply(contrast.matrix, yfilt.avg, mod_design, 1, 0.05, FALSE, NULL)


#setwd(paste0(my.wd,"/Results/DE_Tables"))
#write.table(my.mRNAsST$`STHER2-STLumA`[[1]], "HER2_LumA_up.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(my.mRNAsST$`STHER2-STLumA`[[2]], "HER2_LumA_down.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#write.table(my.mRNAsST$`STHER2-STLumB`[[1]], "HER2_LumB_up.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(my.mRNAsST$`STHER2-STLumB`[[2]], "HER2_LumB_down.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#write.table(my.mRNAsST$`STHER2-STTNBC`[[1]], "HER2_TNBC_up.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(my.mRNAsST$`STHER2-STTNBC`[[2]], "HER2_TNBC_down.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#write.table(my.mRNAsST$`STLumA-STTNBC`[[1]], "LumA_TNBC_up.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(my.mRNAsST$`STLumA-STTNBC`[[2]], "LumA_TNBC_down.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#write.table(my.mRNAsST$`STLumB-STTNBC`[[1]], "LumB_TNBC_up.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(my.mRNAsST$`STLumB-STTNBC`[[2]], "LumB_TNBC_down.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differnetial Expression Analysis with Kmeans Clusters
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


mod_design <-  model.matrix(~0+CLUS+Array+ER+TILSS)

# Create all combinations of BC subtypes for LIMMA DA analysis
combinations <- data.frame(t(combn(paste0("CLUS", levels(CLUS)), 2)))
combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")

# Making group contrasts 
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))

# Apply DE_limma function to all comparisons
my.mRNAsCLUS <- DE_mRNA_apply(contrast.matrix, yfilt.avg, mod_design, 1, 0.05, FALSE, NULL)

#write.table(my.mRNAsCLUS$`CLUSC1-CLUSC2`[[1]], "C1_C2_up.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(my.mRNAsCLUS$`CLUSC1-CLUSC2`[[2]], "C1_C2_down.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differnetial Expression Analysis with Immune-Infiltration 
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DEA with TILS - All levels

mod_design <-  model.matrix(~0+TILS+Array+ER+PGR)

# Create all combinations of BC subtypes for LIMMA DA analysis
combinations <- data.frame(t(combn(paste0("TILS", levels(TILS)), 2)))
combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")

# Making group contrasts 
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))

# Apply DE_limma function to all comparisons
my.mRNAsTILS <- DE_mRNA_apply(contrast.matrix, yfilt.avg, mod_design, 1, 0.05, FALSE, NULL)


#write.table(my.mRNAsTILS$`TILST3-TILST1`[[1]], "T3_T1_up.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(my.mRNAsTILS$`TILST3-TILST1`[[2]], "T3_T1_down.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#write.table(my.mRNAsTILS$`TILST2-TILST1`[[1]], "T2_T1_up.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(my.mRNAsTILS$`TILST2-TILST1`[[2]], "T2_T1_down.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DEA with TILS - High Low level


mod_design <-  model.matrix(~0+TILSS+Array+ER+PGR)

# Create all combinations of BC subtypes for LIMMA DA analysis
combinations <- data.frame(t(combn(paste0("TILSS", levels(TILSS)), 2)))
combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")

# Making group contrasts 
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))

# Apply DE_limma function to all comparisons
my.mRNAsTILSS <- DE_mRNA_apply(contrast.matrix, yfilt.avg, mod_design, 1, 0.05, FALSE, NULL)

#write.table(my.mRNAsTILSS$`TILSSHT-TILSSLT`[[1]], "HT_LT_up.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(my.mRNAsTILSS$`TILSSHT-TILSSLT`[[2]], "HT_LT_down.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differnetial Expression Analysis with Grade
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## DEA with Grade - All levels

mod_design <-  model.matrix(~0+GR+Array+ER+TILS)

# Create all combinations of BC subtypes for LIMMA DA analysis
combinations <- data.frame(t(combn(paste0("GR", levels(GR)), 2)))
combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")

# Making group contrasts 
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))

# Apply DE_limma function to all comparisons
my.mRNAsGR <- DE_mRNA_apply(contrast.matrix, yfilt.avg, mod_design, 1, 0.05, FALSE, NULL)

#write.table(my.mRNAsGR$`GRG3-GRG1`[[1]], "HG_LG_up.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(my.mRNAsGR$`GRG3-GRG1`[[2]], "HG_LG_down.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## DEA with Grade - High Low level

#mod_design <-  model.matrix(~0+GRS+Array+ER+TILS)

## Create all combinations of BC subtypes for LIMMA DA analysis
#combinations <- data.frame(t(combn(paste0("GRS", levels(GRS)), 2)))
#combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")

## Making group contrasts 
#contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))

## Apply DE_limma function to all comparisons
#my.mRNAsGRS <- DE_mRNA_apply(contrast.matrix, yfilt.avg, mod_design, 1, 0.05, FALSE, NULL)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




                                                                                                              ### LASSO Regression ###




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LASSO Regression with Kmeans Clusters, Tumour Grade and Immune-Infiltration
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Clusters
LASSO_C <- Reduce(intersect, list(LASSO_mRNA(1, mrna_combat, clusters2, FALSE),LASSO_mRNA(5, mrna_combat, clusters2, FALSE),LASSO_mRNA(55, mrna_combat, clusters2, FALSE),LASSO_mRNA(555, mrna_combat, clusters2, FALSE),LASSO_mRNA(5050, mrna_combat, clusters2, FALSE)))


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GR
group_GR <- as.integer(GR)
LASSO_GR <- Reduce(intersect, list(LASSO_mRNA(5, mrna_combat, group_GR, TRUE),LASSO_mRNA(5, mrna_combat, group_GR, TRUE),LASSO_mRNA(55, mrna_combat, group_GR, TRUE),LASSO_mRNA(5555, mrna_combat, group_GR, TRUE),LASSO_mRNA(5050, mrna_combat, group_GR, TRUE)))

group_GR2 <- as.integer(GR2)
LASSO_GR2 <- Reduce(intersect, list(LASSO_mRNA(5, mrna_combat, group_GR2, FALSE),LASSO_mRNA(5, mrna_combat, group_GR2, FALSE),LASSO_mRNA(55, mrna_combat, group_GR2, FALSE),LASSO_mRNA(5555, mrna_combat, group_GR2, FALSE),LASSO_mRNA(5050, mrna_combat, group_GR2, FALSE)))


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TILs
group_TILS <- as.integer(TILS)
LASSO_TILS <- Reduce(intersect, list(LASSO_mRNA(1, mrna_combat, group_TILS, TRUE),LASSO_mRNA(5, mrna_combat, group_TILS, TRUE),LASSO_mRNA(55, mrna_combat, group_TILS, TRUE),LASSO_mRNA(555, mrna_combat, group_TILS, TRUE),LASSO_mRNA(5555, mrna_combat, group_TILS, TRUE)))

group_TILS2 <- as.integer(TILS2)
LASSO_TILS2 <- Reduce(intersect, list(LASSO_mRNA(1, mrna_combat, group_TILS2, FALSE),LASSO_mRNA(5, mrna_combat, group_TILS2, FALSE),LASSO_mRNA(55, mrna_combat, group_TILS2, FALSE),LASSO_mRNA(555, mrna_combat, group_TILS2, FALSE),LASSO_mRNA(5555, mrna_combat, group_TILS2, FALSE)))




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Random Forest with Kmeans Clusters, Tumour Grade and Immune-Infiltration
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Clusters
RFCLUS <- Reduce(intersect, list(my_forest(1, t(mrna_combat),clusters2, 660), my_forest(5, t(mrna_combat), clusters2, 660), my_forest(55, t(mrna_combat), clusters2, 660), my_forest(555, t(mrna_combat), clusters2, 660), my_forest(5050, t(mrna_combat), clusters2, 660)))
my_forest_conver(5, t(mrna_combat), clusters2)

# Convergence -TRUE

#set.seed(30)
#VSCLUS  <- varSelRF(t(mrna_combat), clusters2, ntree = 1000, ntreeIterat = 5000, vars.drop.frac = 0.2)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Grade

RFGR <- Reduce(intersect, list(my_forest(1, t(mrna_combat),GR2, 70), my_forest(5, t(mrna_combat), GR2, 70), my_forest(55, t(mrna_combat), GR2, 70), my_forest(555, t(mrna_combat), GR2, 70), my_forest(5050, t(mrna_combat), GR2, 70)))
my_forest_conver(5, t(mrna_combat), GR2)

# Convergence - FALSE

#set.seed(30)
#VSGR  <- varSelRF(t(mrna_combat), GR2, ntree = 1000, ntreeIterat = 5000, vars.drop.frac = 0.2)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TILS

# RF
RFTILS <- Reduce(intersect, list(my_forest(1, t(mrna_combat), TILS2, 360), my_forest(5, t(mrna_combat), TILS2, 360), my_forest(55, t(mrna_combat), TILS2, 360), my_forest(555, t(mrna_combat), TILS2, 360), my_forest(5050, t(mrna_combat), TILS2, 360)))
my_forest_conver(5, t(mrna_combat), TILS2)

# Convergence - SOMEWHAT

# Var. selction
#set.seed(30)
#VSTILS  <- varSelRF(t(mrna_combat), TILS, ntree = 1000, ntreeIterat = 5000, vars.drop.frac = 0.2)




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




                                                                                            ### ANALYSIS OF HORMONE RECEPTOR STATUS RELATED MRNA PROFILES ###




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differnetial Expression Analysis with Hormone Receptor Statuses
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# Estrogen Receptor

mod_design <-  model.matrix(~0+ER+Array+PGR+TILS)

# Create all combinations of BC subtypes for LIMMA DA analysis
combinations <- data.frame(t(combn(paste0("ER", levels(ER)), 2)))
combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")

# Making group contrasts 
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))

# Apply DE_limma function to all comparisons
my.mRNAsER <- DE_mRNA_apply(contrast.matrix, yfilt.avg, mod_design, 1, 0.05, FALSE, NULL)


# Write out tables
#write.table(my.mRNAsER$`ERERp-ERERm`[[1]], "ERp_ERn_up.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(my.mRNAsER$`ERERp-ERERm`[[2]], "ERp_ERn_down.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Progesterone Receptor

mod_design <-  model.matrix(~0+PGR+Array+ER+TILS)

# Create all combinations of BC subtypes for LIMMA DA analysis
combinations <- data.frame(t(combn(paste0("PGR", levels(PGR)), 2)))
combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")

# Making group contrasts 
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))

# Apply DE_limma function to all comparisons
my.mRNAsPGR <- DE_mRNA_apply(contrast.matrix, yfilt.avg, mod_design, 1, 0.05, FALSE, NULL)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# HER2 receptor

mod_design <-  model.matrix(~0+HER2+Array)

# Create all combinations of BC subtypes for LIMMA DA analysis
combinations <- data.frame(t(combn(paste0("HER2", levels(HER2)), 2)))
combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")

# Making group contrasts 
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(mod_design))))

# Apply DE_limma function to all comparisons
my.mRNAsHER2 <- DE_mRNA_apply(contrast.matrix, yfilt.avg, mod_design, 1, 0.05, FALSE, NULL)


HH_HL_up <- rbind(my.mRNAsHER2$`HER2H3-HER2H0`[[1]], my.mRNAsHER2$`HER2H2-HER2H0`[[1]])
HH_HL_up  <- HH_HL_up[order(HH_HL_up$adj.P.Val),]

HH_HL_down <- rbind(my.mRNAsHER2$`HER2H3-HER2H0`[[2]], my.mRNAsHER2$`HER2H2-HER2H0`[[2]])
HH_HL_down  <- HH_HL_down[order(HH_HL_down$adj.P.Val),]

# Write out tables
#write.table(HH_HL_up, "HH_HL_up.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(HH_HL_down, "HH_HL_down.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LASSO with Hormone Receptor Statuses
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# PGR
group_PGR <- as.integer(as.factor(PGR))
LASSO_PGR <- Reduce(intersect, list(LASSO_mRNA(1, mrna_combat, group_PGR, FALSE),LASSO_mRNA(5, mrna_combat, group_PGR, FALSE),LASSO_mRNA(55, mrna_combat, group_PGR, FALSE),LASSO_mRNA(555, mrna_combat, group_PGR, FALSE),LASSO_mRNA(5050, mrna_combat, group_PGR, FALSE)))

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# ER
group_ER <- as.integer(as.factor(ER))
LASSO_ER <- Reduce(intersect, list(LASSO_mRNA(1, mrna_combat, group_ER, FALSE),LASSO_mRNA(5, mrna_combat, group_ER, FALSE),LASSO_mRNA(55, mrna_combat, group_ER, FALSE),LASSO_mRNA(555, mrna_combat, group_ER, FALSE),LASSO_mRNA(5050, mrna_combat, group_ER, FALSE)))



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Random Forest with Hormone Receptor Statuses
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# ER
RFER <-  Reduce(intersect, list(my_forest(1, t(mrna_combat),ER, 320), my_forest(5, t(mrna_combat), ER, 320), my_forest(55, t(mrna_combat), ER, 320), my_forest(555, t(mrna_combat), ER, 320), my_forest(5050, t(mrna_combat), ER, 320)))
my_forest_conver(5, t(mrna_combat), ER)

# Convergence - FALSE

#set.seed(30)
#VSPGR  <- varSelRF(t(mrna_combat), ER, ntree = 1000, ntreeIterat = 5000, vars.drop.frac = 0.2)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PGR
RFPGR <-  Reduce(intersect, list(my_forest(1, t(mrna_combat),PGR, 10), my_forest(5, t(mrna_combat), PGR, 10), my_forest(55, t(mrna_combat), PGR, 10), my_forest(555, t(mrna_combat), PGR, 10), my_forest(5050, t(mrna_combat), PGR, 10)))
my_forest_conver(5, t(mrna_combat), PGR)

# Convergence - SOMEWHAT

#set.seed(30)
#VSPGR  <- varSelRF(t(mrna_combat), PGR, ntree = 1000, ntreeIterat = 5000, vars.drop.frac = 0.2)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




                                                                                                              ### UP-SET PLOTS OF DE MRNAS ###




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load Saved DE Sets From Above
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd(paste0(my.wd,"/Results/DE_Tables"))

C1C2up <- read.delim("C1_C2_up.txt", header = TRUE)$GeneName
C1C2down <- read.delim("C1_C2_down.txt", header = TRUE)$GeneName

HTLTup <- read.delim("HT_LT_up.txt", header = TRUE)$GeneName
HTLTdown <- read.delim("HT_LT_down.txt", header = TRUE)$GeneName

HGLGup <- read.delim("HG_LG_up.txt", header = TRUE)$GeneName
HGLGdown <- read.delim("HG_LG_down.txt", header = TRUE)$GeneName

ERpERnup <- read.delim("ERp_ERn_up.txt", header = TRUE)$GeneName
ERpERndown <- read.delim("ERp_ERn_down.txt", header = TRUE)$GeneName

HHHLup <- read.delim("HH_HL_up.txt", header = TRUE)$GeneName
HHHLdown <- read.delim("HH_HL_down.txt", header = TRUE)$GeneName

#PGRmPGRpup <- read.delim("PGRm_PGRp_up.txt", header = TRUE)$GeneName
#PGRmPGRpdown<- read.delim("PGRm_PGRp_down.txt", header = TRUE)$GeneName

HER2TNBCup <-  read.delim("HER2_TNBC_up.txt", header = TRUE)$GeneName
HER2TNBCdown <-  read.delim("HER2_TNBC_down.txt", header = TRUE)$GeneName

HER2LumAup <-  read.delim("HER2_LumA_up.txt", header = TRUE)$GeneName
HER2LumAdown <-  read.delim("HER2_LumA_down.txt", header = TRUE)$GeneName

HER2LumBup <-  read.delim("HER2_LumB_up.txt", header = TRUE)$GeneName
HER2LumBdown <-  read.delim("HER2_LumB_down.txt", header = TRUE)$GeneName

LumATNBCup <-  read.delim("LumA_TNBC_up.txt", header = TRUE)$GeneName
LumATNBCdown <-  read.delim("LumA_TNBC_down.txt", header = TRUE)$GeneName

LumBTNBCup <-  read.delim("LumB_TNBC_up.txt", header = TRUE)$GeneName
LumBTNBCdown <-  read.delim("LumB_TNBC_down.txt", header = TRUE)$GeneName


# All Subtypes Together
SubtypesAll <- unique(sort(c(as.character(HER2TNBCup) , as.character(HER2TNBCdown), as.character(HER2LumAup), as.character(HER2LumAdown), as.character(HER2LumBup), as.character(HER2LumBdown), as.character(LumATNBCup), as.character(LumATNBCdown), as.character(LumBTNBCup), as.character(LumBTNBCdown))))



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
setwd(paste0(my.wd,"/Results/Plots"))


# Colors
colov <- c("#08605F", "#177E89","#598381", "#99AD93", "#8E936D","#CAD1D1")

# List of Datasets
intsec.list <- list(as.character(SubtypesAll), c(as.character(C1C2up),as.character(C1C2down)), c(as.character(HTLTup),as.character(HTLTdown)), c(as.character(ERpERnup),as.character(ERpERndown)), c(as.character(HHHLup),as.character(HHHLdown)), c(as.character(HGLGup),as.character(HGLGdown))) 
names(intsec.list) <- c("Subtypes","Clus1_Clus2", "TILSH_TILSL", "ERp_ERn", "Her2H_Her2H", "GRH_GRL")

# Plot 
myplot <- plot_upsetR(intsec.list, names(intsec.list), "Overlap_ConsensusSets", colov, TRUE, FALSE)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




                                                                                                ### INTEGRATION OF MUTATIONS, CNVs AND GENE ROLE INFORMATION ###


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load All COSMIC Sets and CancerMine Set
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set WD
my.wd <- "~/Desktop/Thilde/MS_MS_TIF_analysis_2014_2015/TIF_miRNAmRNA/Joint"
setwd(my.wd)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BC specific pathogenic mutations from COSMIC.
CosmicMuts <-  read.delim("/Databases/Cosmic_Muts_clean.txt", header = TRUE)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Cosmic Census Genes
CosmicCancerGenes <- read.delim("/Databases/Cosmic_cancer_gene_sensus_clean.txt", header = TRUE)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Cosmic Copy Number Variations
CosmicCNVs <- read.delim("/Databases/Cosmic_CNVs_clean.txt", header = TRUE)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Cancermine oncogenes and tumour suppressors:
CancerMine <-  read.delim("/Databases/CancerMine_clean.txt", header = TRUE)
CancerMine <- CancerMine[CancerMine$Citation.Sum >= 5, ]






# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Integrate Database Information with DE genes from Modules; Green, Red and Yellow.
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd(paste0(my.wd, "/MiRSytem/Module_Genes"))

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
LuminalTNBCGenes <- as.character(read.delim("ModuleGenes_LuminalTNBC.txt", header = TRUE)$Gene)

LuminalTNBC <- GenesImpact("ModuleGenes_LuminalTNBC.txt", CosmicCancerGenes, CosmicMuts, CosmicCNVs, CancerMine, FALSE)
LuminalTNBC <- LuminalTNBC[match(LuminalTNBCGenes[-6], LuminalTNBC$Gene.Symbol),]

PlotMuts(LuminalTNBC, "LuminalTNBC")


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
LuminalBTNBCGenes <- as.character(read.delim("ModuleGenes_LumBTNBC.txt", header = TRUE)$Gene)

LuminalBTNBC <- GenesImpact("ModuleGenes_LumBTNBC", CosmicCancerGenes, CosmicMuts, CosmicCNVs, CancerMine, FALSE)
LuminalBTNBC <- LuminalBTNBC[match(LuminalBTNBCGenes, LuminalBTNBC$Gene.Symbol),]

PlotMuts(LuminalBTNBC, "LuminalBTNBC")


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TILsGenes <- as.character(read.delim("ModuleGenes_TILs.txt", header = TRUE)$Gene)

TILs <- GenesImpact("ModuleGenes_TILs", CosmicCancerGenes, CosmicMuts, CosmicCNVs, CancerMine, FALSE)
TILs <- TILs[match(TILsGenes, TILs$Gene.Symbol),]

PlotMuts(TILs, "TILs")


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
GradeTILsGenes <- as.character(read.delim("ModuleGenes_TILsGrade.txt", header = TRUE)$Gene)

GradeTILs <- GenesImpact("ModuleGenes_TILsGrade", CosmicCancerGenes, CosmicMuts, CosmicCNVs, CancerMine, FALSE)
GradeTILs <- GradeTILs[match(GradeTILsGenes[-c(4:5)], GradeTILs$Gene.Symbol),]

PlotMuts(GradeTILs, "GradeTILs")




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Overlap all genes from networks with information from COSMIC and CancerMine.
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set WD to load network files:
setwd(paste0(my.wd, "/Tables/Network_Tables"))


# Call function "GenesImpact" on each set, using the COSMIC datasets.
C1C2 <- GenesImpact("C1C2NodeInfo", CosmicCancerGenes, CosmicMuts, CosmicCNVs, CancerMine)
ERpERn <- GenesImpact("ERpERnNodeInfo", CosmicCancerGenes, CosmicMuts, CosmicCNVs, CancerMine)
Her2TNBC <- GenesImpact("Her2TNBCNodeInfo", CosmicCancerGenes, CosmicMuts, CosmicCNVs, CancerMine)
LumATNBC <- GenesImpact("LumATNBCNodeInfo", CosmicCancerGenes, CosmicMuts, CosmicCNVs, CancerMine)
LumBTNBC <- GenesImpact("LumBTNBCNodeInfo", CosmicCancerGenes, CosmicMuts, CosmicCNVs, CancerMine)
HGLG <- GenesImpact("HGLGNodeInfo", CosmicCancerGenes, CosmicMuts, CosmicCNVs, CancerMine)
HTLT <- GenesImpact("HTLTNodeInfo", CosmicCancerGenes, CosmicMuts, CosmicCNVs, CancerMine)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Make Figure 9 - Overlap of Oncogenes, TSGs, Dual Role Genes etc. with Subsets of Genes from Modules.
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in all differentially expressed genes

fl <- list.files(path=paste0(wd,"DE_Tables"))
fl <- fl[grep(".txt", fl)]

setwd(paste0(my.wd,"/DE_Tables"))

DEGs <- lapply(fl, function(x) read.delim(x, header = TRUE))
DEGs <- data.frame(unique(sort(do.call("rbind", DEGs)$GeneName)))
colnames(DEGs) <- "Gene.Symbol"

DEGs$DE.gene <- rep(2, nrow(DEGs))

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in all genes from miRNA-mRNA interaction networks

fl <- list.files(path=paste0(wd,"/Network_Tables"))
fl <- fl[grep("NetworkLFC", fl)]

setwd(paste0(my.wd,"/Network_Tables"))

nws <- lapply(fl, function(x) read.delim(x, header = TRUE))
nws <- do.call("rbind", nws)

# Extract genes from each node and remove miRNAs
nwgenes <- unique(sort(c(as.character(nws$node1), as.character(nws$node2))))
nwgenes <- data.frame(nwgenes[-grep("hsa", nwgenes)])
colnames(nwgenes) <- "Gene.Symbol"

hasMiR <- as.character(unique(nws[grep("hsa", nws$node1),]$node2))
nwgenes$From.NW <- rep(2, nrow(nwgenes))
nwgenes$MiR.Pair <- ifelse(nwgenes$Gene.Symbol %in% hasMiR, 2, 1)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Merge DE genes with genes from networks
MyGenes <- merge(DEGs, nwgenes, by ="Gene.Symbol", all.x = TRUE, all.y = TRUE)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load genes from co-expression modules.
setwd(my.wd)
load("genes_modules_colors.RData")
colnames(genes_colors) <- c("Gene.Symbol", "Module")


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extract information from CancerMine

# Is the gene BC-related
is.BC <- unique(CancerMine[grep("breast", CancerMine$Cancer.Normalized),]$Gene.Symbol)
CancerMine$Breast.Cancer <- ifelse(CancerMine$Gene.Symbol %in% is.BC, 2, 1)

# Conver Driver to Fusion for easier convergence with COSMIC
CancerMine$Role.in.Cancer <- ifelse(CancerMine$Role.in.Cancer == "Driver", "Fusion", as.character(CancerMine$Role.in.Cancer))

CancerMine <- CancerMine[, c(1,5:7)]
CancerMine$In.CancerMine <- rep(2, nrow(CancerMine))


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extract information from COSMIC database

# Is the gene BC-related
is.BC <- CosmicCancerGenes[unique(sort(c(grep("breast", CosmicCancerGenes$Tumour.Types.Somatic.), grep("breast", CosmicCancerGenes$Tumour.Types.Germline.)))),]$Gene.Symbol
CosmicCancerGenes$Breast.Cancer <- ifelse(CosmicCancerGenes$Gene.Symbol %in% is.BC, 2, 1)

CosmicCancerGenes$In.Cosmic <- rep(2, nrow(CosmicCancerGenes))
CosmicCancerGenesSmall <- CosmicCancerGenes[, c(1,6,10,11,15,21,22)]



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Merge COSMIC and CancerMine

Cosmic.CancerMine <- merge(CosmicCancerGenesSmall, CancerMine, by = "Gene.Symbol", all.x = TRUE, all.y = TRUE)

# Consensus for gene Role.in.Cancer and is Breast.Cancer related
Cosmic.CancerMine$Role.in.Cancer <- ifelse(is.na(Cosmic.CancerMine$Role.in.Cancer.x), as.character(Cosmic.CancerMine$Role.in.Cancer.y), as.character(Cosmic.CancerMine$Role.in.Cancer.x))
Cosmic.CancerMine$Breast.Cancer <- ifelse(is.na(Cosmic.CancerMine$Breast.Cancer.x), as.character(Cosmic.CancerMine$Breast.Cancer.y), as.character(Cosmic.CancerMine$Breast.Cancer.x))


# Make key to change role to numeric
KeyNumeric <- data.frame("Role.in.Cancer" = levels(as.factor(Cosmic.CancerMine$Role.in.Cancer)), "Role" = as.numeric(c(4,5,6,3)))
Cosmic.CancerMine <- merge(Cosmic.CancerMine, KeyNumeric, by ="Role.in.Cancer", all.x = TRUE, all.y =FALSE)

# Small and ordered set
Cosmic.CancerMine <- Cosmic.CancerMine[, c(2,1,8,9,12:14)]
Cosmic.CancerMine <- Cosmic.CancerMine[order(Cosmic.CancerMine$Gene.Symbol),]



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Merge Cosmic.CancerMine with set of DE network genes
Cosmic.CancerMine <- merge(Cosmic.CancerMine, MyGenes, by = "Gene.Symbol", all.x = TRUE, all.y =TRUE)

# Replace NA with nuemric 1 for plotting
Cosmic.CancerMine [is.na(Cosmic.CancerMine)] <- 1

# Filter dataset to strict version, either in COSMIC or BC-related or From network or Citation.Sum > 10.
Cosmic.CancerMine  <- Cosmic.CancerMine[which(Cosmic.CancerMine$In.Cosmic == 2 | Cosmic.CancerMine$Breast.Cancer == 2 | Cosmic.CancerMine$From.NW == 2 | Cosmic.CancerMine$Citation.Sum > 10),]

Cosmic.CancerMine <- Cosmic.CancerMine[, c(1,7,3,5,6,8,9)]


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Merge Cosmic.CancerMine with modules

InModules <- merge(Cosmic.CancerMine, genes_colors, by = "Gene.Symbol")
InModules$Breast.Cancer <- as.numeric(InModules$Breast.Cancer)
InModules$Role <- as.numeric(InModules$Role)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot

# Melt for ggplot
InModules <- melt(InModules)

# Split by module
GeneMod <- split(InModules, f = InModules$Module)

# Module names
ModName <- names(GeneMod)

# Color for plotting
my.colors <- data.frame("lev" = c(1:6), "my.color" = c("white", "#C0CAAD", "#3FB5A5", "#7B8BD8", "#CED67A", "#FCB97D"))

# Plot all modules
for(idx in 1:length(GeneMod)) {
  my.df <- GeneMod[[idx]]
  my.name <- ModName[[idx]]
  OncoTSP(my.df, my.name, 3, 14)
}




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Fishers Exact Test - Enrichment of Oncogenes and TSGs in Modules.
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Get all oncogenes and TSGs
OncoTSGs <- data.frame("Gene.Symbol" = as.character(Cosmic.CancerMine$Gene.Symbol))

# Merge oncogenes and TSGs with modules
IsOncoTSG <- merge(OncoTSGs, genes_colors, by = "Gene.Symbol")

# Split by modules
IsOncoTSG <- split(IsOncoTSG, f = IsOncoTSG$Module)

# Count number of oncogenes and TSGs in each modules
OT <- data.frame(unlist(lapply(IsOncoTSG, function(x) nrow(x))))
OT$mod <- rownames(OT)
colnames(OT) <- c("OncoTGs", "mod")
rownames(OT) <- NULL

# Count number of genes in each modules
ModulesSplit <- split(genes_colors, f = genes_colors$Module)
M <- data.frame(unlist(lapply(ModulesSplit, function(x) nrow(x))))
M$mod <- rownames(M)
colnames(M) <- c("Total", "mod")
rownames(M) <- NULL


# Merge into one dataframe
Fishers.res <- merge(OT, M, by = "mod")

# Set up contingency tables
Fishers.res$IsNOT <- Fishers.res$Total-Fishers.res$OncoTGs
Fishers.res$totalwitout <- sum(Fishers.res$Total)-Fishers.res$Total
Fishers.res$totalwitoutOnco <- sum(Fishers.res$OncoTGs)-Fishers.res$OncoTGs

# Open lists for loop
oRs <- list()
ps <- list()

for (idx in 1:nrow(Fishers.res)) {
  f <- fisher.test(data.frame(c(Fishers.res$OncoTGs[[idx]], Fishers.res$totalwitoutOnco[[idx]]), c(Fishers.res$IsNOT[[idx]], Fishers.res$totalwitout[[idx]])))
  oRs[[idx]] <- f$estimate
  ps[[idx]] <- f$p.value
}

oRs <- as.numeric(unlist(oRs))
ps <- unlist(ps)
ps.adj <- p.adjust(ps, method = "fdr")

# Final dataset
Fishers.res$ps.adj <- ps.adj
Fishers.res$OR <- oRs
