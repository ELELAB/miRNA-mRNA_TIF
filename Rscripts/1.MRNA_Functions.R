# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Copyright (C) 2018, Thilde Bagger Terkelsen <thilde@cancer.dk> <thildebate@gmail.com>
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






### FUNCTIONS ###








# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD PACKAGES
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(cluster)
library(data.table)
library(dendextend)
library(ggplot2)
library(glmnet)
library(heatmap.plus)
library(pamr)
library(plyr)
library(RColorBrewer)
library(reshape)
library(rms)
library(scales)
library(squash)
library(VennDiagram)

library(biomaRt)
library(limma)
library(RnAgilentDesign028282.db)
library(sva)
library(UpSetR)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# MULTIDIMENSIONAL SCALING PLOT
# Takes as arguments:
# my.data = matrix of countes/expression/abundance
# my.groups = vector of group IDs for coloring
# my.labels = vector of group IDs for labeling

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


myMDSplot <- function(my.data, my.group, my.labels) {
  d<-dist(t(my.data))
  fit <- cmdscale(d,eig=TRUE, k=2)
  res<-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
  p <- ggplot(data=res) + geom_point(aes(x=M1,y=M2,color=my.group)) + geom_text(aes(x=M1,y=M2, label= my.labels, color=my.group)) +
    coord_cartesian(xlim=c(min(res$M1)*1.4,max(res$M1)*1.4)) + theme_minimal() + theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold"), axis.text.x=element_blank()) + guides(colour = guide_legend(override.aes = list(size=6))) + theme(legend.position = "top") + theme(axis.text=element_text(size=16, face="bold")) + theme(axis.text = element_text(colour = "black"))
  return(p)
}




# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# MAPPING MISSING ENSEMBL IDS TO GENE NAME
# Takes as arguments:
# dataset.list = Object from limma read.maimages, list object with dataset and metadata
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

MapENST <- function(dataset.list) {
  # Get ENST IDs
  ENST <- grep("ENST", dataset.list$genes$GeneName)
  ENST <- unique(dataset.list$genes$GeneName[ENST])
  ensembl <- useDataset("hsapiens_gene_ensembl",mart = useMart("ensembl"))
  my.entrez <- unique(getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol"), filters = "ensembl_transcript_id", values = ENST, mart = ensembl))
  missing <- grep("^[[:space:]]*$", my.entrez$hgnc_symbol)
  my.entrez <- my.entrez[-missing,]
  colnames(my.entrez) <- c("GeneName", "symbol")
  
  # Merge IDs unto dataset
  dataset.list$genes$order <- 1:nrow(dataset.list$genes)
  dataset.list$genes <- merge(dataset.list$genes, my.entrez, by = "GeneName", all.x = TRUE)
  dataset.list$genes$GeneName <- ifelse(!is.na(dataset.list$genes$symbol), as.character(dataset.list$genes$symbol), as.character(dataset.list$genes$GeneName))
  dataset.list$genes$symbol <- NULL
  dataset.list$genes <- dataset.list$genes[order(dataset.list$genes$order),]
  rownames(dataset.list$E) <- dataset.list$genes$GeneName
  return(dataset.list)
}




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR CLUSTER ANALYSIS
# optimal_nc function takes as arguments:
# my.dataframe = a dataframe of expression/abundance values
# plot_clusters function takes as arguments:
# my.dataframe = a dataframe of expression/abundance values
# my.clusters = number of clusters to plot
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Calculating number of optimal clusters - clustGap 


optimal_nc <- function (my.dataframe) {
  pam1 <- function(my.dataframe,k) list(cluster = pam(my.dataframe,k, cluster.only=TRUE))
  optimal_number <- clusGap(t(my.dataframe), FUN = pam1, K.max = 15, B = 500)
  plot(optimal_number)
  return(optimal_number)
}




# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR GETTING HEATMAP COLORS
# Takes as arguments:
# my.truestatus = a vetcor of groups/labels (a character vector, length of ncol in the matrix to be plotted)
# my.cols = a vector with colors to use (a character vector with the length of the number of groups/levels).
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


get_colors <- function(my.truestatus, my.cols) {
  hm_col <- data.frame(levels(as.factor(as.character(my.truestatus))), my.cols)
  colnames(hm_col) <- c("status", "mycolor")
  true_status <- data.frame(my.truestatus)
  myorder <- 1:nrow(true_status)
  true_status$order <- myorder
  colnames(true_status) <- c("status", "order")
  col <- merge(true_status, hm_col, by="status", all.x =TRUE)
  col <- col[order(col$order),]
  col$mycolor <- ifelse(is.na(col$mycolor), "white", as.character(col$mycolor))
  return(as.matrix(col$mycolor))
}



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION TO OBTAIN DIFFERENTIALLY ABUNDANT HITS:
# Takes as arguments;
# my.contrast = a contrast between groups of interest
# my.data = a dataframe with expression/abundance counts
# my.design = a design matrix with all comparisons
# my.coLFC, my.coFDR = cutoffs for logFC and FDR
# If blocking than a vector of patient IDs
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


DE_mRNA <- function(my.contrast, my.data, my.design, my.coLFC, my.coFDR, my.block=NULL) {
  arrayw <- arrayWeights(my.data)
  if(is.null(my.block)) {
    fit3 <- eBayes(contrasts.fit(lmFit(object = my.data, design = my.design, weights = arrayw), my.contrast))
  }
  else {
    corfit <- duplicateCorrelation(my.data, my.design, block=my.block) 
    fit3 <- eBayes(contrasts.fit(lmFit(object = my.data, design = my.design, block = my.block, weights = arrayw, correlation=corfit$consensus), my.contrast))
  }
  tt <- topTable(fit3, coef=1, adjust='fdr', number=nrow(my.data))
  
  up <- tt[tt$logFC >= my.coLFC & tt$adj.P.Val < my.coFDR, ]
  down <- tt[tt$logFC <= -my.coLFC & tt$adj.P.Val < my.coFDR, ]
  
  up$dir <- rep("up", nrow(up))
  down$dir <- rep("down", nrow(down))
  
  final <- list(up, down)
  return(final)
}




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTIONS TO APPLY DIFFERENTIALLY ABUNDANCE ANALYSIS TO ALL COMPARISONS AT ONCE:
# Takes as arguments;
# my.contrasts= all contrasts between groups of interest
# my.data = a dataframe with expression/abundance counts
# my.design = a design matrix with all comparisons
# my.coLFC, my.coFDR = cutoffs for logFC and FDR
# If blocking than a vector of patient IDs
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


DE_mRNA_apply <- function(my.contrasts, my.data, my.design, my.coLFC, my.coFDR, my.vector, my.block=NULL) {
  my.mRNAs.l <- apply(my.contrasts, 2, function(x) DE_mRNA(x, my.data, my.design, my.coLFC, my.coFDR, my.block)) 
  if(my.vector == TRUE) {
    my.mRNAs <- do.call(rbind, lapply(my.mRNAs.l, function(x) do.call(rbind, x)))
    my.mRNAs <- unique(do.call(rbind, strsplit(rownames(my.mRNAs), "[.]"))[,2])
    return(my.mRNAs)
  }
  else {
    return(my.mRNAs.l)
  }
}




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTION FOR DA ANALYSIS WITH CLINICAL PARAMETERS. THE FUNCTION CALLS "DA_miRNA_apply" FROM ABOVE.
# Takes as arguments;
# my.data = a dataframe with expression/abundance counts
# my.design = a design matrix with all comparisons
# a cutoff for logFC and FDR
# my.coLFC, my.coFDR = cutoffs for logFC and FDR
# If blocking than a vector of patient IDs
# If remove is different from NULL, a vector of indices to remove must be supplied

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


DE_all_contrasts <- function(my.data, my.design, my.group, my.group.name, my.logFC, my.FDR, my.block=NULL, my.remove=NULL) {
  if (!is.null(my.remove)) {
    my.data <- my.data[, -my.remove]
    my.group <- my.group[-my.remove]
  }
  combinations<- data.frame(t(combn(paste0(my.group.name, levels(my.group)), 2)))
  combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")
  contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(my.design))))
  my.DE <- DE_mRNA_apply(contrast.matrix, my.data, my.design, my.logFC, my.FDR, FALSE, my.block)
  return(my.DE)
}





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR LASSO REGRESSION
# Takes as arguments:
# my.seed = A random seed
# my.data = matrix of countes/expression/abundance
# my.groups = vector of group IDs, must be as.integer()
# If my.multinorm=TRUE the my.groups vector has > 2 groups, else my.multinorm=FALSE which will result in binomial regression.

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



LASSO_mRNA <- function(my.seed, my.data, my.groups, my.multinorm=TRUE) {
  if(my.multinorm == TRUE) {
    set.seed(my.seed)
    my.fit <- cv.glmnet(x = t(my.data), y = my.groups, family="multinomial", type.multinomial = "grouped", nfolds = 10)
    my.coef <- coef(my.fit, s=my.fit$lambda.min)
    my.ma <- as(my.coef$`1`, "matrix")
    rm(my.fit)
    rm(my.coef)
  }
  else {
    set.seed(my.seed)
    my.fit <- cv.glmnet(x = t(my.data), y = my.groups, family = "binomial", type.measure = "class", nfolds = 10)
    my.coef <- coef(my.fit, s=my.fit$lambda.min)
    my.ma <- as(my.coef, "matrix")
    rm(my.fit)
    rm(my.coef)
  }
  my.ma <- names(my.ma[my.ma[,1] != 0, ])
  #my.ma <- my.data[rownames(my.data) %in% my.ma, ]
  return(my.ma)
}





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR RANDOM FOREST
# Takes as arguments:
# my.seed = A random seed
# my.data = matrix of countes/expression/abundance
# my.groups = vector of group IDs, must be as.factor()
# my.nhits = Number of hits to return from RF ranking.
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



my_forest_conver <- function(my.seed, my.data, my.groups) {
  set.seed(my.seed)
  rf <- randomForest(x=my.data, y=my.groups, ntree=300)
  plot(rf)
  return(rf)
}



my_forest <- function(my.seed, my.data, my.groups, my.nhits) {
  set.seed(my.seed)
  rf <- randomForest(x=my.data, y=my.groups, ntree=300, replace = TRUE)
  # plot(rf)
  var.imp <- data.frame(importance(rf,type=2))
  # make row names as columns
  var.imp$Variables <- row.names(var.imp)
  var.imp <- var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),]
  rf <- gsub(pattern = "[.]", replacement = "-", var.imp$Variables[1:my.nhits])
  return(rf)
}





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot UpSetR
# Takes as arguments;
# list.of.sets = list with sets to plot
# my.intersection = the names of the sets to intersect
# my.name = name of output plot
# my.cols = colors vector with as, one color per set
# if my.plot= TRUE a pdf is writen out, og write.ids = TRUE, write out intersection of all sets
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


plot_upsetR <- function(list.of.sets, my.intersection, my.name, my.cols, my.plot, write.ids) {
  full.set <- data.frame(unique(sort(c(unlist(list.of.sets)))))
  colnames(full.set) <- "Accession"
  for (name in  names(list.of.sets)) {
    full.set <- data.frame(full.set, ifelse(full.set$Accession %in% as.character(list.of.sets[[name]]), 1, 0))
  }
  colnames(full.set) <- c("Accession", names(list.of.sets))
  metadata <- data.frame("sets" = colnames(full.set)[-1], "sets2" = colnames(full.set)[-1])
  if (my.plot==TRUE) {
    pdf(paste0(my.name, ".pdf"), height = 7, width = 11)
    my.combination <- my.cols
    names(my.combination) <- my.intersection
    upset(full.set, sets=colnames(full.set)[2:ncol(full.set)], sets.bar.color = my.cols, set.metadata = list(data = metadata, plots = list(list(type="matrix_rows", column = "sets", colors = my.combination, alpha = 0.5))), order.by = "freq", text.scale = 1.7, keep.order = TRUE) 
    dev.off()
  }
  if (write.ids == TRUE) {
    idx <- which(names(list.of.sets) %in% my.intersection)
    write_out(Reduce(intersect, list.of.sets[idx]), my.name)
  } 
  return(full.set)
}





