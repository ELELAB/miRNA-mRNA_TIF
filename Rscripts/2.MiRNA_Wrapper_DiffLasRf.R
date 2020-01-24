# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Copyright (C) 2018, Thilde Bagger Terkelsen <thilde@cancer.dk> <thildebate@gmail.com>
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Working directory
my.wd <- ""

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





                                                                                                            ### LOAD, PREPROCESS AND FILTER ###





# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load MiRNA Data and Metadata
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd(paste0(my.wd,"/Data/"))

NIFTIF <- read.delim("miRNA_all_samples_NIFTIF.txt", header = TRUE)
NIFTIFinfo <- read.table("miRNA_all_samples_NIFTIFinfo.txt", header = TRUE)

setwd(paste0(my.wd,"/Backgrounds_and_Databases/"))
Converter <- read.delim("TIFmiRNAConverter.txt", header = TRUE)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Filtering and Preprocessing of MiRNA Data
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Add correct IDs
NIFTIF$ID <- Converter$finalname
rownames(NIFTIF) <- NULL

# Average duplicates
NIFTIF <- data.frame(data.table(NIFTIF)[, lapply(.SD, mean), by=ID])
rownames(NIFTIF) <- NIFTIF$ID
NIFTIF$ID <- NULL

# Remove samples < 60% tumor, remove technical replicates, remove Apocrine type and outlier samples TIF78 
remove <- sort(c(which(NIFTIFinfo$tp %in% c("10", "30", "40")), which(is.na(NIFTIFinfo$tp)), which(NIFTIFinfo$ID %in% c("TIF302_2","TIF81a","TIF234_2", "TIF78")), which(NIFTIFinfo$Type == "Apocrine")))
NIFTIF <- NIFTIF[,-remove]
NIFTIFinfo <- NIFTIFinfo[-remove,]


# Number of tumor samples
Tn <- (ncol(NIFTIF)-(length(grep("tumor", NIFTIFinfo$TN))))+1


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Filtering out rows where less TIF samples than number of samples in (second)-smallest group (HER2, 8 samples) have values => 0 (0.006)
greater_than_background <- data.frame(apply(NIFTIF[,Tn:ncol(NIFTIF)], 1, function(x) sum(x > 0.006)))
greater_than_TNBC <- which(greater_than_background[,1] < 8)
NIFTIF <- NIFTIF[-greater_than_TNBC,]


# Sub missing values with min per miRNA
min_per_miRNA <- as.vector(apply(NIFTIF, 1, function(x) min(x[x > 0.006])))

for(i in 1:nrow(NIFTIF)){
  NIFTIF[i, NIFTIF[i,] == 0.006] <- min_per_miRNA[i]
}


 


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TIF samples only

TIF <- NIFTIF[, Tn:ncol(NIFTIF)]
TIFinfo <- NIFTIFinfo[Tn:ncol(NIFTIF),]

#setwd(paste0(my.wd,"/Data/"))
#write.table(TIF, "TIF.txt", sep ="\t", quote =FALSE, col.names =TRUE, row.names =TRUE)
#write.table(TIFinfo, "TIFinfo.txt", sep ="\t", quote =FALSE, col.names =TRUE, row.names =TRUE)




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Factor vectors

TN <- as.factor(as.character(NIFTIFinfo$TN))
TSN <- as.factor(as.character(NIFTIFinfo$Tumor_subtype_corrected_2015_11_20))
TSN <- factor(as.character(TSN), levels = c("normal", "LumA", "LumB", "LumB_HER2_enriched", "HER2", "TNBC"))
TS <- as.factor(as.character(TIFinfo$Tumor_subtype_corrected_2015_11_20))
TS <- factor(as.character(TS), levels = c("LumA", "LumB", "LumB_HER2_enriched", "HER2", "TNBC"))

TNBC <-  as.factor(as.character(ifelse(!(TSN %in% c("TNBC", "normal")), "Other", as.character(TSN))))
patient <- as.factor(paste0("p", NIFTIFinfo$patient))


ER <- factor(as.character(ifelse(TIFinfo$ER=="ER+", "ERp", "ERm")),levels = c("ERp", "ERm"))
PGR <- factor(as.character(ifelse(TIFinfo$PGR=="PGR+", "PGRp", "PGRm")),levels = c("PGRp", "PGRm"))
HER2 <- as.factor(as.character(paste0("H", gsub("[+]", "", TIFinfo$HER2))))
AR <- factor(as.character(ifelse(TIFinfo$AR=="AR+", "ARp", "ARm")),levels = c("ARp", "ARm"))

GR <- factor(paste0("G", ifelse(as.character(TIFinfo$Gr) == "1", "2", as.character(TIFinfo$Gr))), levels = c("G3", "G2"))
GRN <- factor(as.character(ifelse(NIFTIFinfo$Gr == "1", "2", as.character(NIFTIFinfo$Gr))),  levels = c("3", "2","0"))
TILS <- factor(as.character(ifelse(TIFinfo$TILS == "T3_outside_tumor", "T3", as.character(TIFinfo$TILS))), levels = c("T3", "T2", "T1", "T0"))
TILSS <- as.factor(as.character(ifelse(TILS %in% c("T0", "T1"), "LT", "HT")))
TILSN <- factor(as.character(ifelse(NIFTIFinfo$TILS == "T3_outside_tumor", "T3", as.character(NIFTIFinfo$TILS))), levels = c("T3", "T2", "T1", "T0", "TN"))



# Color vectors
TN.cols <- c("grey60","violetred4")
#TSN.cols <- c("purple", "aquamarine3", "darkblue", "pink", "orange", "grey60")
TSN.cols <- c("grey60","#04724D","#BDCC9F","#D2AB99","#DD614A","#310A31")
TS.cols  <- c("#04724D","#BDCC9F","#D2AB99","#DD614A","#310A31")

colov <- c("#08605F", "#177E89", "#CAD1D1", "#598381","#8E936D")
Cl.cols <- c("#D86654", "#254441", "#43AA8B", "#B2B09B")




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







                                                                                                            ###  VISUALIZATION ###







# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Multidimensional Scaling Plots
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
setwd(paste0(my.wd,"/Results/"))


# Combat batch correction for patient effects - only for visualization
datacorr <- ComBat(as.matrix(log2(NIFTIF)), batch = patient)

#pdf("MDSTN.pdf", height = 9, width = 10)

# Corrected, color by tumour or normal
TNcomb <- myMDSplot(datacorr, TN, "", TN.cols)
# Not corrected, color by tumour or normal
TNnone <- myMDSplot(log2(NIFTIF), TN, "", TN.cols)


# Not corrected, color by tumour subtypes
mod_design <-  model.matrix(~TS)
TSnone <- myMDSplot(log2(TIF), TS, "", TS.cols)

#grid.arrange(grobs = list(TNnone, TNcomb, TSnone), ncol = 2, main = "Main title")
#dev.off()

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Combat batch correction for patient effects - only for visualization
mod_design <-  model.matrix(~TSN)
datacorr <- ComBat(as.matrix(log2(NIFTIF)), mod = mod_design, batch = patient)

# Corrected, color by tumour or normal
TSNcomb <- myMDSplot(datacorr, TSN, "", TSN.cols)

# Not corrected, color by tumour subtypes or normal
TSNnone <- myMDSplot(log2(NIFTIF), TSN, "", TSN.cols)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Clustering Analysis
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Determine optimal number of clusters

optimal_nc(log2(TIF))

# Kmeans clustering - 2 clusters
set.seed(10)
K2.clus <- kmeans(t(log2(TIF)), 2)
clusters2 <- as.factor(as.character(K2.clus$cluster))
CLUS <- as.factor(paste0("C", clusters2))

# Kmeans clustering - 3 clusters
set.seed(10)
K3.clus <- kmeans(t(CP1), 3)
clusters3 <- as.factor(as.character(K3.clus$cluster))




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Dendogram of MiRNA Samples
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Plot dendogram
dend <- as.dendrogram(hclust(dist(t(log2(TIF))), method = "ward.D2"))

# Color scheme
my.ER <- get_colors(ER, c("#FFFCF7","grey60"))
my.PGR <- get_colors(PGR, c("#FFFCF7","grey60"))
my.GR <- get_colors(TIFinfo$Gr, c("#0071AA","#0071AA","#032A63"))
my.TILS <- get_colors(TIFinfo$TILS, c("#F7D76F" ,"#EFCA4F", "#EDAB49" ,"#EE964B", "#EE964B"))
my.CLUS2 <- get_colors(clusters2, c("#43AA8B", "#B2B09B"))
my.CLUS3 <- get_colors(clusters3, c("#43AA8B","#696773", "#B2B09B"))
my.MFS <- get_colors(TIFinfo$Meta_free_surv, c("white", "black"))
my.TS <- get_colors(TS, c("#FFAAA3", "#6F73E2", "#BFDDFF", "#BFDDFF", "#FF5465"))

#pdf("Dendogram_miRNA.pdf", height = 10, width = 14)

par(cex=0.6, mar = c(8,3,3,8))
nodePar <- list(lab.cex = 1.5, pch = c(NA, 19), col = "black")
plot(dend, horiz = TRUE, nodePar = nodePar)
colored_bars(cbind(my.TS, my.ER, my.PGR, my.TILS, my.GR, my.CLUS2), dend, rowLabels = c("TS", "ER", "PGR", "TILs", "GR", "Clusters"), horiz = TRUE, cex.rowLabels = 1.1)

#dev.off()



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Scatter Plot and FVIZ Cluster Plot
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Kmeans clustering - 2 clusters
set.seed(10)
K3.clus <- kmeans(t(CP1), 3)

# FVIZ cluster plot
fviz_cluster(K3.clus, t(CP1), palette = c("#43AA8B","#696773", "#B2B09B"), ggtheme = theme_minimal())

# Distances, corrected for patient effects
d<-dist(t(CP1))
fit <- cmdscale(d,eig=TRUE, k=3)
res<-data.frame(M1=fit$points[,1],M2=fit$points[,2],M3=fit$points[,3])

# Scatterplot
scatterplot3d(res, pch = 16, cex.symbols = 2, color = my.CLUS3, grid=TRUE, angle = 110)




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







                                                                                                            ###   DIFFERENTIAL EXPRESSION ANALYSIS ###








# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential Expression Analysis with TIF vs NIF and TIF Subtypes
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Tumor normal
my_design_TN <- model.matrix(~0+TN)
my_contrast <- makeContrasts("TNtumor-TNnormal", levels = my_design_TN)
DE_TN <- DE_miRNA(my_contrast, log2(NIFTIF), my_design_TN, 1, 0.05, patient)


#setwd(paste0(my.wd,"/Results/DE_Tables"))
#write.table(DE_TN[[1]], "TIF_NIF_up.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
#write.table(DE_TN[[2]], "TIF_NIF_down.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Subtypes with normal
my_design_TSN <- model.matrix(~0+TSN)
DE_TSN <- DE_all_contrasts(log2(NIFTIF), my_design_TSN, TSN, "TSN", 1, 0.05, patient)

my_design_TS <- model.matrix(~0+TS+TILSS)
DE_TS <- DE_all_contrasts(log2(TIF), my_design_TS, TS, "TS", 1, 0.05)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Subtypes no normal
#my_design_TS <- model.matrix(~0+TS)
#DE_TS <- DE_all_contrasts(log2(TIF), my_design_TS, TS, "TS", 1, 0.05)
#DE_subtypes_TS <- unique(sort(c(rownames(DE_TS$`TSHER2-TSLumA`[[1]]), rownames(DE_TS$`TSHER2-TSLumA`[[2]]), rownames(DE_TS$`TSHER2-TSTNBC`[[1]]),  rownames(DE_TS$`TSLumA-TSTNBC`[[1]]),  rownames(DE_TS$`TSLumA-TSTNBC`[[2]]), rownames(DE_TS$`TSLumB-TSTNBC`[[1]]),  rownames(DE_TS$`TSLumB-TSTNBC`[[2]]))))


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TNBC normal
my_design_TNBCN <- model.matrix(~0+TNBC)
DE_TNBCN <- DE_all_contrasts(log2(NIFTIF), my_design_TNBCN, TNBC, "TNBC", 1, 0.05, patient)

#write.table(DE_TNBCN$`TNBCOther-TNBCTNBC`[[1]], "Other_TNBC_up.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
#write.table(DE_TNBCN$`TNBCOther-TNBCTNBC`[[2]], "Other_TNBC_down.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


TNBC <- as.factor(as.character(TNBC[Tn:length(TNBC)])) 
my_design_TNBC <- model.matrix(~0+TNBC)
DE_TNBC<- DE_all_contrasts(log2(TIF), my_design_TNBC, TNBC, "TNBC", 1, 0.05)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TNBC no normal
#TNBC <- as.factor(as.character(TNBC[Tn:ncol(NIFTIF)]))
#my_design_TNBC <- model.matrix(~0+TNBC)
#my_contrast <- makeContrasts("TNBCOTHER-TNBCTNBC", levels = my_design_TNBC)
#DE_TNBC_NN <- DE_miRNA(my_contrast, log2(TIF), my_design_TNBC, 1, 0.05)
#DE_TNBC_NN <-  unique(sort(c(rownames(DE_TNBC_NN[[1]]), rownames(DE_TNBC_NN[[2]]))))





# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential Expression Analysis with Kmeans Clusters
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#CLUSN <- as.factor(c(rep("CN", Tn-1), CLUS))
#my_design_CLUSN <- model.matrix(~0+CLUSN+patient)
#DE_CLUSN <- DE_all_contrasts(log2(NIFTIF), my_design_CLUSN, CLUSN, "CLUSN", 1, 0.05)


my_design_CLUS <- model.matrix(~0+CLUS+TILS+ER+PGR)
DE_CLUS <- DE_all_contrasts(log2(TIF), my_design_CLUS, CLUS, "CLUS", 1, 0.05)

# Write out Tables

#write.table(DE_CLUS$`CLUSC1-CLUSC2`[[1]], "C1_C2_up.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
#write.table(DE_CLUS$`CLUSC1-CLUSC2`[[2]], "C1_C2_down.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential Expression Analysis with Hormone Receptor Statuses
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ER
#ERN  <- as.factor(c(rep("N",Tn-1), as.character(ER)))
#my_design_ERN <- model.matrix(~0+ERN)
#DE_ERN <- DE_all_contrasts(log2(NIFTIF), my_design_ERN, ERN, "ERN", 1, 0.05, patient)


my_design_ER <- model.matrix(~0+ER+PGR+TILSS)
DE_ER <- DE_all_contrasts(log2(TIF), my_design_ER, ER, "ER", 1, 0.05)

# Write out Tables

#write.table(DE_ER$`ERERp-ERERm`[[1]], "ERp_ERn_up.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
#write.table(DE_ER$`ERERp-ERERm`[[2]], "ERp_ERn_down.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PGR
#PGRN  <- as.factor(c(rep("N",Tn-1), as.character(PGR)))
#my_design_PGRN <- model.matrix(~0+PGRN)
#DE_PGRN <- DE_all_contrasts(log2(NIFTIF), my_design_PGRN, PGRN, "PGRN", 1, 0.05, patient)

my_design_PGR <- model.matrix(~0+PGR+ER+TILSS)
DE_PGR <- DE_all_contrasts(log2(TIF), my_design_PGR, PGR, "PGR", 1, 0.05)

# Write out Tables

#write.table(DE_PGR$`PGRPGRp-PGRPGRm`[[1]], "PGRp_PGRn_up.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
#write.table(DE_PGR$`PGRPGRp-PGRPGRm`[[2]], "PGRp_PGRn_down.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# HER2
#HER2N  <- as.factor(c(rep("N",Tn-1), as.character(HER2)))
#my_design_HER2N <- model.matrix(~0+HER2N)
#DE_HER2N <- DE_all_contrasts(log2(NIFTIF), my_design_HER2N, HER2N, "HER2N", 1, 0.05, patient)

#my_design_HER2 <- model.matrix(~0+HER2)
#DE_HER2 <- DE_all_contrasts(log2(TIF), my_design_HER2, HER2, "HER2", 1, 0.05)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# AR
#ARN  <- as.factor(c(rep("N",Tn-1), as.character(AR)))
#my_design_ARN <- model.matrix(~0+ARN)
#DE_ARN <- DE_all_contrasts(log2(NIFTIF), my_design_ARN, ARN, "ARN", 1, 0.05, patient)

#my_design_AR <- model.matrix(~0+AR+ER+PGR)
#DE_AR <- DE_all_contrasts(log2(TIF), my_design_AR, AR, "AR", 1, 0.05)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Receptor
#RC <- as.factor(paste0(ER, PGR))
#my_design_RC <- model.matrix(~0+RC)
#DE_RC <- DE_all_contrasts(log2(NIFTIF), my_design_RC, RC, "RC", 1, 0.05, patient)







# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Differential Expression Analysis with Tumour Grade and Immune-Infiltration
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TILS

TILSNS <- as.character(ifelse(TILSN %in% c("T0", "T1"), "LT", as.character(TILSN)))
TILSNS <- as.factor(as.character(ifelse(TILSNS %in% c("T2", "T3"), "HT", as.character(TILSNS))))

my_design_TILSN <- model.matrix(~0+TILSN)
DE_TILSN <- DE_all_contrasts(log2(NIFTIF), my_design_TILSN, TILSN, "TILSN", 1, 0.05, patient)

#my_design_TILSNS <- model.matrix(~0+TILSNS)
#DE_TILSNS <- DE_all_contrasts(log2(NIFTIF), my_design_TILSNS, TILSNS, "TILSNS", 1, 0.05, patient)

my_design_TILS <- model.matrix(~0+TILS+ER+PGR)
DE_TILS <- DE_all_contrasts(log2(TIF), my_design_TILS, TILS, "TILS", 1, 1)

#my_design_TILSS <- model.matrix(~0+TILSS+ER+PGR)
#DE_TILSS <- DE_all_contrasts(log2(TIF), my_design_TILSS, TILSS, "TILSS", 1, 1)


intup <- intersect(rownames(DE_TILSN$`TILSNT3-TILSNT1`[[1]]), rownames(DE_TILS$`TILST3-TILST1`[[1]]))
intdown <- intersect(rownames(DE_TILSN$`TILSNT3-TILSNT1`[[2]]), rownames(DE_TILS$`TILST3-TILST1`[[2]]))

TILSup <- DE_TILSN$`TILSNT3-TILSNT1`[[1]][rownames(DE_TILSN$`TILSNT3-TILSNT1`[[1]]) %in% intup, ]
TILSdown <- DE_TILSN$`TILSNT3-TILSNT1`[[2]][rownames(DE_TILSN$`TILSNT3-TILSNT1`[[2]]) %in% intdown, ]

# Write out Tables

#write.table(TILSup, "HT_LT_up.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
#write.table(TILSdown, "HT_LT_down.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GRADE

Gradeset1 <- NIFTIF[,-which(is.na(NIFTIFinfo$Gr))]
GRpatient <- as.factor(as.character(patient[-which(is.na(NIFTIFinfo$Gr))]))
Gradeset2 <- TIF[,-which(is.na(TIFinfo$Gr))]
GR <- factor(as.character(GR[-which(is.na(TIFinfo$Gr))]), levels = c("G3","G2"))
GRER <- as.factor(as.character(ER[-which(is.na(TIFinfo$Gr))])) 
GRPGR <- as.factor(as.character(PGR[-which(is.na(TIFinfo$Gr))])) 
GRTILSS <- as.factor(as.character(TILSS[-which(is.na(TIFinfo$Gr))]))

GRN  <- factor(c(rep("N",Tn-1), as.character(GR)), levels = c("G3", "G2", "N"))
my_design_GRN <- model.matrix(~0+GRN)
DE_GRN <- DE_all_contrasts(log2(Gradeset1), my_design_GRN, GRN, "GRN", 1, 0.05, GRpatient)

my_design_GR <- model.matrix(~0+GR+GRER+GRPGR+GRTILSS)
DE_GR <- DE_all_contrasts(log2(Gradeset2), my_design_GR, GR, "GR", 1, 1)

intup <- intersect(rownames(DE_GRN$`GRNG3-GRNG2`[[1]]), rownames(DE_GR$`GRG3-GRG2`[[1]]))
intdown <- intersect(rownames(DE_GRN$`GRNG3-GRNG2`[[2]]), rownames(DE_GR$`GRG3-GRG2`[[2]]))

GRup <- DE_GRN$`GRNG3-GRNG2`[[1]][rownames(DE_GRN$`GRNG3-GRNG2`[[1]]) %in% intup, ]
GRdown <- DE_GRN$`GRNG3-GRNG2`[[2]][rownames(DE_GRN$`GRNG3-GRNG2`[[2]]) %in% intdown, ]

# Write out Tables

#write.table(GRup, "HG_LG_up.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
#write.table(GRdown, "HG_LG_down.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)





# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                                                              




                                                                                                                      ### UP-SET PLOTS OF DE MIRNAS ###





# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Load DE Datasets written out from analysis above.

setwd(paste0(my.wd,"/Results/DE_Tables"))

# TIF_NIF
TIFNIFup <- rownames(read.delim("TIF_NIF_up.txt", header = TRUE))
TIFNIFdown <- rownames(read.delim("TIF_NIF_down.txt", header = TRUE))
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Other vs TNBC
OtherTNBCup <- rownames(read.delim("Other_TNBC_up.txt", header = TRUE))
OtherTNBCdown <- rownames(read.delim("Other_TNBC_down.txt", header = TRUE))
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Cluster 1 vs Cluster 2
C1C2up <- rownames(read.delim("C1_C2_up.txt", header = TRUE))
C1C2down <- rownames(read.delim("C1_C2_down.txt", header = TRUE))
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ER+ vs ER-
ERpERnup <- rownames(read.delim("ERp_ERn_up.txt"))
ERpERndown <- rownames(read.delim("ERp_ERn_down.txt"))
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# PgR+ vs PgR-
PGRpPGRnup <- rownames(read.delim("PGRp_PGRn_up.txt"))
PGRpPGRndown <- rownames(read.delim("PGRp_PGRn_down.txt"))
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#  High TILs vs Low TILs
HTLTup <- rownames(read.delim("HT_LT_up.txt", header = TRUE))
HTLTdown <- rownames(read.delim("HT_LT_down.txt", header = TRUE))
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#  High Grade vs Low/Medium Grade
HGLGup <- rownames(read.delim("HG_LG_up.txt", header = TRUE))
HGLGdown <- rownames(read.delim("HG_LG_down.txt", header = TRUE))




setwd(paste0(my.wd,"/Results/"))

# Color scheme
colov <- c("#496377", "#177E89", "#08605F","#598381", "#B5CA8D", "#CAD1D1", "#99AD93")

# List of all sets
intsec.list <- list(c(as.character(TIFNIFup), as.character(TIFNIFdown)), c(as.character(C1C2up),as.character(C1C2down)), c(as.character(OtherTNBCup), as.character(OtherTNBCdown)), c(as.character(HTLTup),as.character(HTLTdown)), c(as.character(PGRpPGRnup), as.character(PGRpPGRndown)), c(as.character(HGLGup),as.character(HGLGdown)),c(as.character(ERpERnup),as.character(ERpERndown))) 

# Names of sets
names(intsec.list) <- c("TIFNIF", "Clus1_Clus2", "Subtypes","TILSH_TILSL", "PgRp_PgRn", "GRH_GRL", "ERp_ERn")

# Plotting
myplot <- plot_upsetR(intsec.list, names(intsec.list), "Overlap_ConsensusSets", colov, TRUE, FALSE)





# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







                                                                                                            ### RANDOM FOREST ###









# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Random Forest TIF vs NIF
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###  Seeds 1, 5, 55, 555, 5050

# Tumor and Normal
#Vset <- c(as.vector(window(which(NIFTIFinfo$TN == "normal"), deltat=3)), as.vector(window(which(NIFTIFinfo$TN == "tumor"), deltat=3)))
#Tset <- seq(1:ncol(NIFTIF))[!seq(1:ncol(NIFTIF)) %in% Vset] 

# Validation set
#Vsetdat <- data.frame(t(log2(NIFTIF[,Vset])))
#VsetY <- TN[Vset]

# Test set
#Tsetdat <- data.frame(t(log2(NIFTIF[,Tset])))
#TsetY <- TN[Tset]

# RF
#RFTIFNIF <- my_forest(30, Tsetdat, TsetY, 286)
# Var. selection
#set.seed(30)
#VSTIFNIF  <- varSelRF(Tsetdat, TsetY, ntree = 10000, ntreeIterat = 5000, vars.drop.frac = 0.2)

RFTIFNIF <- Reduce(intersect, list(my_forest(1, t(log2(NIFTIF)), TN, 286), my_forest(5, t(log2(NIFTIF)), TN, 286), my_forest(55, t(log2(NIFTIF)), TN, 286), my_forest(555, t(log2(NIFTIF)), TN, 286), my_forest(5050, t(log2(NIFTIF)), TN, 286)))
my_forest_conver(1, t(log2(NIFTIF)), TN)

# Convergence - TRUE



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Random Forest Tumour Grade Aand Immune-Infiltration
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TILS
TILSY <- as.character(TILS[Tn:length(TILS)])
TILSY <- as.factor(ifelse(TILSY %in% c("T0", "T1"), "low", "high"))

# RF
RFTILS <- Reduce(intersect, list(my_forest(1, t(log2(TIF)), TILSY, 108), my_forest(5, t(log2(TIF)), TILSY, 108), my_forest(55, t(log2(TIF)), TILSY, 108), my_forest(555, t(log2(TIF)), TILSY, 108), my_forest(5050, t(log2(TIF)), TILSY, 108)))
my_forest_conver(1, t(log2(TIF)), TILSY)

# Convergence - SOMEWHAT

# Var. selction
#set.seed(30)
#VSTILS  <- varSelRF(t(log2(TIF)), TILSY, ntree = 10000, ntreeIterat = 5000, vars.drop.frac = 0.2)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GRADE
TIFGR <- TIF[,-which(GR == "GNA")]
GRY <- as.factor(as.character(GR[-which(GR == "GNA")]))

RFGR <- Reduce(intersect, list(my_forest(1, t(log2(TIFGR)), GRY, 56), my_forest(5, t(log2(TIFGR)), GRY, 56), my_forest(55, t(log2(TIFGR)), GRY, 56), my_forest(555, t(log2(TIFGR)), GRY, 56), my_forest(5050, t(log2(TIFGR)), GRY, 56)))
my_forest_conver(1, t(log2(TIFGR)), GRY)


# Convergence - SOMEWHAT

#set.seed(30)
#VSGR  <- varSelRF(t(log2(TIF)), GR, ntree = 10000, ntreeIterat = 5000, vars.drop.frac = 0.2)




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Random Forest Kmeans Clusters
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#  CLUSTERS
RFCLUS <- Reduce(intersect, list(my_forest(1, t(log2(TIF)), clusters2, 196),my_forest(5, t(log2(TIF)), clusters2, 196),my_forest(55, t(log2(TIF)), clusters2, 196),my_forest(555, t(log2(TIF)), clusters2, 196),my_forest(5050, t(log2(TIF)), clusters2, 196)))
my_forest_conver(1, t(log2(TIF)), clusters2)

# Convergence - TRUE

#set.seed(30)
#VSCLUS  <- varSelRF(t(log2(TIF)), clusters, ntree = 10000, ntreeIterat = 5000, vars.drop.frac = 0.2)







# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Random Forest Hormone Receptor Statuses
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# ER
RFER <- Reduce(intersect, list(my_forest(1, t(log2(TIF)), ER, 110),my_forest(5, t(log2(TIF)), ER, 110),my_forest(55, t(log2(TIF)), ER, 110),my_forest(555, t(log2(TIF)), ER, 110),my_forest(5050, t(log2(TIF)), ER, 110)))
my_forest_conver(1, t(log2(TIF)), ER)

# Convergence - FALSE

#set.seed(30)
#VSPGR  <- varSelRF(t(log2(TIF)), ER, ntree = 10000, ntreeIterat = 5000, vars.drop.frac = 0.2)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PGR
RFPGR <- Reduce(intersect, list(my_forest(1, t(log2(TIF)), PGR, 110),my_forest(5, t(log2(TIF)), PGR, 110),my_forest(55, t(log2(TIF)), PGR, 110),my_forest(555, t(log2(TIF)), PGR, 110),my_forest(5050, t(log2(TIF)), PGR, 110)))
my_forest_conver(1, t(log2(TIF)), PGR)

# Convergence - SOMEWHAT

#set.seed(30)
#VSPGR  <- varSelRF(t(log2(TIF)), PGR, ntree = 10000, ntreeIterat = 5000, vars.drop.frac = 0.2)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------








                                                                                                            ### LASSO REGRESSION ###












# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LASSO Regression
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# Tumor Normal
group_TN <- as.integer(as.factor(TN))
LASSO_TN <- Reduce(intersect, list(LASSO_miRNA(1, log2(NIFTIF), group_TN, FALSE), LASSO_miRNA(5, log2(NIFTIF), group_TN, FALSE), LASSO_miRNA(55, log2(NIFTIF), group_TN, FALSE), LASSO_miRNA(555, log2(NIFTIF), group_TN, FALSE), LASSO_miRNA(5050, log2(NIFTIF), group_TN, FALSE)))

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Clusters
LASSO_C <- Reduce(intersect, list(LASSO_miRNA(1, log2(TIF), clusters, FALSE),LASSO_miRNA(5, log2(TIF), clusters, FALSE),LASSO_miRNA(55, log2(TIF), clusters, FALSE),LASSO_miRNA(555, log2(TIF), clusters, FALSE),LASSO_miRNA(5050, log2(TIF), clusters, FALSE)))

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PGR
group_PGR <- as.integer(as.factor(PGR))
LASSO_PGR <- Reduce(intersect, list(LASSO_miRNA(1, log2(TIF), group_PGR, FALSE), LASSO_miRNA(5, log2(TIF), group_PGR, FALSE), LASSO_miRNA(55, log2(TIF), group_PGR, FALSE), LASSO_miRNA(555, log2(TIF), group_PGR, FALSE), LASSO_miRNA(5050, log2(TIF), group_PGR, FALSE)))

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GR
group_GR <- as.integer(GRY)
LASSO_GR <- Reduce(intersect, list(LASSO_miRNA(1, log2(TIFGR), group_GR , FALSE), LASSO_miRNA(5, log2(TIFGR), group_GR , FALSE), LASSO_miRNA(55, log2(TIFGR), group_GR, FALSE), LASSO_miRNA(555, log2(TIFGR), group_GR, FALSE), LASSO_miRNA(5050, log2(TIFGR), group_GR, FALSE)))

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TILS
group_TILS <- as.integer(TILSY)
LASSO_TILS <- Reduce(intersect, list(LASSO_miRNA(1, log2(TIF), group_TILS, FALSE), LASSO_miRNA(5, log2(TIF), group_TILS, FALSE), LASSO_miRNA(55, log2(TIF), group_TILS, FALSE), LASSO_miRNA(555, log2(TIF), group_TILS, FALSE), LASSO_miRNA(5050, log2(TIF), group_TILS, FALSE)))



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




                                                                                                          ### CONSENSUS SETS ###







# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#setwd(paste0(my.wd,"/Results/consensus_DE_RF_LASSO/"))


# TIF vs NIF
TIF_NIF_up_consensus <- unique(c(intersect(RFTIFNIF, rownames(DE_TN[[1]])), intersect(LASSO_TN, rownames(DE_TN[[1]]))))
TIF_NIF_down_consensus <- unique(c(intersect(RFTIFNIF, rownames(DE_TN[[2]])), intersect(LASSO_TN, rownames(DE_TN[[2]]))))

# Write out Tables

#write.table(data.frame(TIF_NIF_up_consensus), "TIF_NIF_up_consensus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(data.frame(TIF_NIF_down_consensus), "TIF_NIF_down_consensus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)



# High TILS vs Low TILS
TILS_high_low_up_consensus <- unique(c(intersect(RFTILS, rownames(DE_TILSN$`TILSNhigh-TILSNlow`[[1]])),  intersect(LASSO_TILS, rownames(DE_TILSN$`TILSNhigh-TILSNlow`[[1]]))))
TILS_high_low_down_consensus <- unique(c(intersect(RFTILS, rownames(DE_TILSN$`TILSNhigh-TILSNlow`[[2]])), intersect(LASSO_TILS, rownames(DE_TILSN$`TILSNhigh-TILSNlow`[[2]]))))

# Write out Tables

#write.table(data.frame(TILS_high_low_up_consensus), "TILS_high_low_up_consensus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(data.frame(TILS_high_low_down_consensus), "TILS_high_low_down_consensus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


# PGR- vs PGR +
GR3_GR2_up_consensus <- unique(c(intersect(RFGR, rownames(DE_GR$`GRNG3-GRNG2`[[1]])), intersect(LASSO_GR, rownames(DE_GR$`GRNG3-GRNG2`[[1]]))))
GR3_GR2_down_consensus <- unique(c(intersect(RFGR, rownames(DE_GR$`GRNG3-GRNG2`[[2]])), intersect(LASSO_GR, rownames(DE_GR$`GRNG3-GRNG2`[[2]]))))

# Write out Tables

#write.table(data.frame(GR3_GR2_up_consensus), "GR3_GR2_up_consensus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(data.frame(GR3_GR2_down_consensus), "GR3_GR2_down_consensus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


# PGR- vs PGR +
PGRm_PGRp_up_consensus <- unique(c(intersect(RFPGR, rownames(DE_PGR$`PGRNPGRm-PGRNPGRp`[[1]])), intersect(LASSO_PGR, rownames(DE_PGR$`PGRNPGRm-PGRNPGRp`[[1]]))))
PGRm_PGRp_down_consensus <- unique(c(intersect(RFPGR, rownames(DE_PGR$`PGRNPGRm-PGRNPGRp`[[2]])), intersect(LASSO_PGR, rownames(DE_PGR$`PGRNPGRm-PGRNPGRp`[[2]]))))

# Write out Tables

#write.table(data.frame(PGRm_PGRp_up_consensus), "PGRm_PGRp_up_consensus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(data.frame(PGRm_PGRp_down_consensus), "PGRm_PGRp_down_consensus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


# Cluster 1 vs Cluster 2
C1_C2_up_consensus <- unique(c(intersect(RFCLUS, rownames(DE_clusters$`clustersNC1-clustersNC2`[[1]])), intersect(LASSO_C, rownames(DE_clusters$`clustersNC1-clustersNC2`[[1]]))))
C1_C2_down_consensus <- unique(c(intersect(RFCLUS, rownames(DE_clusters$`clustersNC1-clustersNC2`[[2]])), intersect(LASSO_C, rownames(DE_clusters$`clustersNC1-clustersNC2`[[2]]))))

# Write out Tables

#write.table(data.frame(C1_C2_up_consensus), "C1_C2_up_consensus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(data.frame(C1_C2_down_consensus), "C1_C2_down_consensus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# UPsetR plot



# Overlap of Datasets
intsec.list <- list(c(as.character(TIF_NIF_up_consensus),as.character(TIF_NIF_down_consensus)), c(as.character(C1_C2_up_consensus),as.character(C1_C2_down_consensus)),c(as.character(PGRm_PGRp_up_consensus),as.character(PGRm_PGRp_down_consensus)), c(as.character(TILS_high_low_up_consensus),as.character(TILS_high_low_down_consensus)), c(as.character(GR3_GR2_up_consensus), as.character(GR3_GR2_down_consensus))) 
names(intsec.list) <- c("TIF_NIF", "C1_C2", "PGRm_PGRp", "HTILS_LTILS", "HGR_LGR")
myplot <- plot_upsetR(intsec.list, names(intsec.list), "Overlap_ConsensusSets", colov, TRUE, FALSE)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

