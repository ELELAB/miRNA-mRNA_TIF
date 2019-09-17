# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Copyright (C) 2018, Thilde Bagger Terkelsen <thilde@cancer.dk> <thildebate@gmail.com>
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






                                                                                                        ### FUNCTIONS ###








# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD PACKAGES
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




library(caTools)
library(data.table)
library(dendextend)
library(factoextra)
library(ggplot2)
library(glmnet)
library(gridExtra)
library(heatmap.plus)
library(pamr)
library(RColorBrewer)
library(randomForest)
library(reshape)
library(rms)
library(scales)
library(varSelRF)
library(VennDiagram)
library(viridis)
library(WGCNA)


library(arcdiagram)
library(biomaRt)
library(limma)
library(pathview)
library(sva)
library(UpSetR)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Multidimensional Scaling Plot:
# Takes as arguments;
    # my.data = a dataframe of expression/abundance counts
    # my.group, my.labels = a vector of IDs for coloring and labeling (may be the same or different, length should be equal to ncol(dataframe))
    # my.cols = a vector of colors (one color for each group)
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


myMDSplot <- function(my.data, my.group, my.labels, my.cols) {
    d<-dist(t(my.data))
    fit <- cmdscale(d,eig=TRUE, k=2)
    res<-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
    p <- ggplot(data=res)
    p + geom_point(aes(x=M1,y=M2,color=my.group), size=2) + geom_text(aes(x=M1,y=M2, label= my.labels, color=my.group)) + scale_color_manual(values  = my.cols) +
    coord_cartesian(xlim=c(min(res$M1)*1.4,max(res$M1)*1.4)) + theme_bw() + theme(legend.title=element_blank()) + theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold")) + guides(colour = guide_legend(override.aes = list(size=6))) + theme(legend.position = "top") + theme(axis.text=element_text(size=16, face="bold")) + theme(axis.text = element_text(colour = "black"))
}






# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION TO OBTAIN DIFFERENTIALLY ABUNDANT HITS:
# Takes as arguments;
    # my.contrast = a contrast between groups of interest
    # my.data = a dataframe with expression/abundance counts
    # my.design = a design matrix with all comparisons
    # my.coLFC, my.coFDR = cutoffs for logFC and FDR
    # If blocking than a vector of patient IDs
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



DE_miRNA <- function(my.contrast, my.data, my.design, my.coLFC, my.coFDR, my.block=NULL) {
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




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTIONS TO APPLY DIFFERENTIALLY ABUNDANCE ANALYSIS TO ALL COMPARISONS AT ONCE:
# Takes as arguments;
    # my.contrasts= all contrasts between groups of interest
    # my.data = a dataframe with expression/abundance counts
    # my.design = a design matrix with all comparisons
    # my.coLFC, my.coFDR = cutoffs for logFC and FDR
    # If blocking than a vector of patient IDs
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



DE_miRNA_apply <- function(my.contrasts, my.data, my.design, my.coLFC, my.coFDR, my.vector, my.block=NULL) {
  my.miRNAs.l <- apply(my.contrasts, 2, function(x) DE_miRNA(x, my.data, my.design, my.coLFC, my.coFDR, my.block)) 
  if(my.vector == TRUE) {
    my.miRNAs <- do.call(rbind, lapply(my.miRNAs.l, function(x) do.call(rbind, x)))
    my.miRNAs <- unique(do.call(rbind, strsplit(rownames(my.miRNAs), "[.]"))[,2])
    return(my.miRNAs)
  }
  else {
    return(my.miRNAs.l)
  }
}



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR DA ANALYSIS WITH CLINICAL PARAMETERS. THE FUNCTION CALLS "DA_miRNA_apply" FROM ABOVE.
# Takes as arguments;
    # my.data = a dataframe with expression/abundance counts
    # my.design = a design matrix with all comparisons
    # a cutoff for logFC and FDR
    # my.coLFC, my.coFDR = cutoffs for logFC and FDR
    # If blocking than a vector of patient IDs
    # If remove is different from NULL, a vector of indices to remove must be supplied
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

DE_all_contrasts <- function(my.data, my.design, my.group, my.group.name, my.logFC, my.FDR, my.block=NULL, my.remove=NULL) {
  if (!is.null(my.remove)) {
    my.data <- my.data[, -my.remove]
    my.group <- my.group[-my.remove]
  }
  combinations<- data.frame(t(combn(paste0(my.group.name, levels(my.group)), 2)))
  combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")
  contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(my.design))))
  my.DE <- DE_miRNA_apply(contrast.matrix, my.data, my.design, my.logFC, my.FDR, FALSE, my.block)
  return(my.DE)
}






# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR GETTING HEATMAP COLORS
# Takes as arguments:
    # my.truestatus = a vetcor of groups/labels (a character vector, length of ncol in the matrix to be plotted)
    # my.cols = a vector with colors to use (a character vector with the length of the number of groups/levels).
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


get_colors <- function(my.truestatus, my.cols) {
    hm_col <- data.frame(levels(as.factor(as.character(my.truestatus))), my.cols)
    colnames(hm_col) <- c("status", "mycolor")
    true_status <- data.frame(my.truestatus)
    myorder <- 1:nrow(true_status)
    true_status$order <- myorder
    colnames(true_status) <- c("status", "order")
    col <- merge(true_status, hm_col, by="status", all.x =TRUE)
    col <- col[order(col$order),]
    col$mycolor <- ifelse(is.na(col$mycolor), "black", as.character(col$mycolor))
    return(as.matrix(col$mycolor))
}




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR CLUSTER ANALYSIS
# optimal_nc function takes as arguments:
    # my.dataframe = a dataframe of expression/abundance values
# plot_clusters function takes as arguments:
    # my.dataframe = a dataframe of expression/abundance values
    # my.clusters = number of clusters to plot
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Calculating number of optimal clusters - clustGap 

optimal_nc <- function (my.dataframe) {
  pam1 <- function(my.dataframe,k) list(cluster = pam(my.dataframe,k, cluster.only=TRUE))
  optimal_number <- clusGap(t(my.dataframe), FUN = pam1, K.max = 15, B = 500)
  plot(optimal_number)
}


# Plotting clusters with patient IDs

plot_clusters <- function(my.dataframe, my.clusters) {
  d <- dist(t(my.dataframe))
  clust <- kmeans(d, my.clusters)
  library(fpc)
  library(cluster)
  clusplot(as.matrix(d), clust$cluster, color=TRUE, shade=TRUE, col.txt=col, col.clus = c("dodgerblue3", "magenta","springgreen2"), col.p=col, labels=2, cex=0.6, lines=0, main="CLUSPLOT")
}




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR OBTAINING miRNA CLUSTERS - get expression/abundance matrix for clusters, remove subtypes with too few samples
# takes arguments:
    # my.clusterIDs = a vector clusters (example with three clusters: c(1,2,1,2,2,3,3,2,3,3,1,2,1,3))
    # my.datainfo = metadata (must contain columns: TN (tumor or normal sample), ID (sample ID) and Tumor_subtype_corrected_2015_11_20 (breast cancer subtypes))
    # my.data = dataframe of expression/abundance values
    # my.remove = vector of subtype groups to remove (example: c("TNBC", "HER2"))
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


get_clusters <- function(my.clusterIDs, my.datainfo, my.data, my.remove) {
  NN <- length(which(my.datainfo$TN == "normal"))
  CLUS <- c(as.character(my.datainfo$ID[1:NN]), my.clusterIDs)
  CLUS <- my.datainfo[my.datainfo$ID %in% CLUS,]
  CLUS <- CLUS[-which(CLUS$Tumor_subtype_corrected_2015_11_20 %in% my.remove),]
  CLUSdata <- my.data[colnames(my.data) %in% as.character(CLUS$ID)]
  CLUSTS <- as.factor(as.character(CLUS$Tumor_subtype_corrected_2015_11_20))
  CLUSdataTS <- list(CLUSTS, CLUSdata)
  return(CLUSdataTS)
}




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR RANDOM FOREST
# Takes as arguments:
    # my.seed = A random seed
    # my.data = matrix of countes/expression/abundance
    # my.groups = vector of group IDs, must be as.factor()
    # my.nhits = Number of hits to return from RF ranking.
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



my_forest_conver <- function(my.seed, my.data, my.groups) {
  set.seed(my.seed)
  rf <- randomForest(x=my.data, y=my.groups, ntree=1000)
  plot(rf)
  return(rf)
}



my_forest <- function(my.seed, my.data, my.groups, my.nhits) {
  set.seed(my.seed)
  rf <- randomForest(x=my.data, y=my.groups, ntree=1000, replace = TRUE)
  # plot(rf)
  var.imp <- data.frame(importance(rf,type=2))
  # make row names as columns
  var.imp$Variables <- row.names(var.imp)
  var.imp <- var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),]
  rf <- gsub(pattern = "[.]", replacement = "-", var.imp$Variables[1:my.nhits])
  return(rf)
}



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR LASSO REGRESSION
# Takes as arguments:
    # my.seed = A random seed
    # my.data = matrix of countes/expression/abundance
    # my.groups = vector of group IDs, must be as.integer()
    # If my.multinorm=TRUE the my.groups vector has > 2 groups, else my.multinorm=FALSE which will result in binomial regression.

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


LASSO_miRNA <- function(my.seed, my.data, my.groups, my.multinorm=TRUE) {
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



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# KEGG FUNCTION
# Takes as arguments:
    # my.pathways = A pathway list object. Each pathway as name and genes or miRNAs assigned to name as vector.
    # my.background = vector of background genes or miRNAs
    # my.query = queryset of genes or miRNAs
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


mirnaPwEnrich <- function(my.pathways, my.background, my.query, simple=NULL) {
  BGSetcount <- data.frame(as.numeric(unlist(lapply(my.pathways, function(x) length(intersect(x, my.background))))))
  colnames(BGSetcount) <- c("inPWBG")
  BGSetcount$notPWBG <- length(my.background)-BGSetcount$inPWBG
  BGSetcount$isPWSet  <- as.numeric(unlist(lapply(my.pathways, function(x) length(intersect(x, my.query)))))
  BGSetcount$notPWSet <- length(my.query)-BGSetcount$isPWSet
  rownames(BGSetcount) <- names(my.pathways)
  my.fishers <- apply(BGSetcount, 1, function(x) fisher.test(matrix(x,nr=2))[1:3])
  if (!is.null(simple)) {
    my.fishers <- do.call(rbind, lapply(my.fishers, function(x) data.frame(x[1][[1]], x[3][[1]], x[2][[1]][1], x[2][[1]][2])))
    colnames(my.fishers) <- c("Pval", "OR", "Lower", "Upper")
    my.fishers$FDR <- p.adjust(my.fishers$Pval, method = "fdr")
  }
  return(my.fishers)
}



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot UpSetR
# Takes as arguments;
    # list.of.sets = list with sets to plot
    # my.intersection = the names of the sets to intersect
    # my.name = name of output plot
    # my.cols = colors vector with as, one color per set
    # if my.plot= TRUE a pdf is writen out, og write.ids = TRUE, write out intersection of all sets
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


plot_upsetR <- function(list.of.sets, my.intersection, my.name, my.cols, my.plot, write.ids) {
  full.set <- data.frame(unique(sort(c(unlist(list.of.sets)))))
  colnames(full.set) <- "Accession"
  for (name in  names(list.of.sets)) {
    full.set <- data.frame(full.set, ifelse(full.set$Accession %in% as.character(list.of.sets[[name]]), 1, 0))
  }
  colnames(full.set) <- c("Accession", names(list.of.sets))
  metadata <- data.frame("sets" = colnames(full.set)[-1], "sets2" = colnames(full.set)[-1])
  if (my.plot==TRUE) {
    pdf(paste0(my.name, ".pdf"), height = 6, width = 10)
    my.combination <- my.cols
    names(my.combination) <- my.intersection
    p <- upset(full.set, sets=colnames(full.set)[2:ncol(full.set)], sets.bar.color = my.cols, set.metadata = list(data = metadata, plots = list(list(type="matrix_rows", column = "sets", colors = my.combination, alpha = 0.5))), order.by = "freq", text.scale = 1.7, keep.order = TRUE) 
    print(p)
    dev.off()
  }
  if (write.ids == TRUE) {
    idx <- which(names(list.of.sets) %in% my.intersection)
    write_out(Reduce(intersect, list.of.sets[idx]), my.name)
  } 
  return(full.set)
}






# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# miRNA families list + plot, multiple datasets
# Takes as arguments;
        # A list containing N vectors of miRNAs to be compared.
        # A dataframe with miRNAs assigned to miRNA families (must contain at least two columns, Name (name of miRNA) and Family (miRNA family))
        # A vector of names (same length as length of list above yo name each set)
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


map_miRFam <- function(my.list.of.lists, my.miRFams, vec.of.names) {
  alist <- list()
  for (idx in 1:length(my.list.of.lists)) {
    my.list <- my.list.of.lists[[idx]]
    df <- unique(my.miRFams[my.miRFams$Name %in% my.list,])
    dfplot <- data.frame(table(as.character(df$Family)))
    colnames(dfplot) <- c("Family", as.character(vec.of.names[[idx]]))
    alist[[idx]] <-  dfplot
  }
  
  # Merge all dataframes
  mirfams <- Reduce(function(...) merge(..., by="Family", all=TRUE), alist)

  # Clean up data for plotting
  mirfams[is.na(mirfams)] <- 0
  mirfams <- mirfams[rowSums(mirfams[,2:ncol(mirfams)]) > 3, ]
  mirfams <- melt(mirfams, id.vars = "Family")
  colnames(mirfams) <- c("Family", "DEset", "miRNAs")
  families <- unique(as.character(mirfams$Family))
  
  # Plot
  mirfams$Family <- factor(mirfams$Family, levels=c(families[order(nchar(families), families, decreasing = TRUE)]))
  p1 <- ggplot(mirfams , aes(x = DEset, y = Family)) + geom_tile(aes(fill = miRNAs), color = "white") +
    scale_fill_gradient2(low="white", mid="#F3E8EE", high="#475B63", midpoint=0, limits=range(mirfams$miRNAs)) + scale_y_discrete(limits=levels(as.factor(mirfams$Family))) + theme_classic() +
    coord_fixed(ratio = 0.2) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10), axis.text.y = element_text(size=7)) + xlab("") + ylab("")
  pdf("miRFamPlot.pdf", height = 14, width = 6)
  print(p1)
  dev.off()
  return(mirfams)
}




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Make Network of differentially expreassed miRNA-Gene targets and Gene-Gene interactions.
# Takes as arguments;
      # miRNAup = list of up-regulated miRNAs
      # miRNAdown = list of down-regulated miRNAs
      # mRNAup = list of up-regulated mRNAs
      # mRNAdown = list of down-regulated mRNAs
      # Database with miRNA-mRNA interactions
      # Database with gene-gene interactions
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



DENetWork <- function(miRNAup, miRNAdown, mRNAup, mRNAdown, miGeIntDB, GeGeIntDB, strict = TRUE) {
  
  # miRNA-Gene Interactions
  miGeup <- miGeIntDB[miGeIntDB$miRNAs %in% miRNAup,]
  miGedown <- miGeIntDB[miGeIntDB$miRNAs %in% miRNAdown,]
  miRNAs <- rbind(miGeup, miGedown)
  miRNAs$miRNAdir <- c(as.character(rep("up", nrow(miGeup))), as.character(rep("down", nrow(miGedown))))
  
  # DE genes only
  mRNAs <- data.frame(c(as.character(mRNAup), as.character(mRNAdown)), c(as.character(rep("up", length(mRNAup))), as.character(rep("down", length(mRNAdown)))))
  colnames(mRNAs) <- c("genes", "mRNAdir")
  
  # Merge datasets
  pairs <- merge(miRNAs, mRNAs, by="genes", all.x=FALSE, all.y=FALSE)
  
  # Oposite directionality only
  pairsOposite <- pairs[pairs$miRNAdir != pairs$mRNAdir,]
  pairsOposite <- pairsOposite[,match(c("miRNAs", "genes", "miRNAdir", "mRNAdir"), colnames(pairsOposite))]
  colnames(pairsOposite) <- c("node1", "node2", "dirnode1", "dirnode2" )
  
  # Gene-Gene interactions
  PPIs <- GeGeIntDB[GeGeIntDB$node1 %in% c(as.character(mRNAup), as.character(mRNAdown)),]
  DEPPIs <- PPIs[PPIs$node2 %in% c(as.character(mRNAup), as.character(mRNAdown)),]
  
  # Remove self-Interactions
  DEPPIs <- DEPPIs[as.character(DEPPIs$node1) != as.character(DEPPIs$node2),]
  
  DEPPIs$dirnode1 <- ifelse(DEPPIs$node1 %in% mRNAup, "up", "down")
  DEPPIs$dirnode2 <- ifelse(DEPPIs$node2 %in% mRNAup, "up", "down")
  
  if(strict == TRUE) {
    IntNode1 <- intersect(pairsOposite$node2, DEPPIs$node1)
    IntNode1 <- DEPPIs[DEPPIs$node1 %in% IntNode1, ]
    IntNode2 <- intersect(pairsOposite$node2, DEPPIs$node2)
    IntNode2 <- DEPPIs[DEPPIs$node2 %in% IntNode2, ]
    DEPPIs <- unique(rbind(IntNode1, IntNode2))
  }
  
  # Final Network
  finalpairs <- rbind(pairsOposite, DEPPIs)
  finalpairs$dir <- rep("o", nrow(finalpairs))
  finalpairs$dir <- ifelse(finalpairs$dirnode1 == "up" & finalpairs$dirnode2 == "up", "u", as.character(finalpairs$dir))
  finalpairs$dir <- ifelse(finalpairs$dirnode1 == "down" & finalpairs$dirnode2 == "down", "d", as.character(finalpairs$dir))
  return(unique(finalpairs))
}








# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Adds and averages logFCs and adjusted p.values of differentially expreassed miRNA-Gene targets and Gene-Gene interactions.
        # miGePairs = miRNA-mRNA and mRNA-mRNA pairs (Object from DENetwork)
        # Dataframe with miRNA IDs and logFCs, p-values and directionality
        # Dataframe with mRNA IDs and loFCs, p-values and directionality
        # Conversion file with miRNA names (can be omitted)
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



DENetWorkFC <- function(miGePairs, milogFCs, GelogFCs) {
  
  milogFCs$node1 <- rownames(milogFCs)
  milogFCs <- milogFCs[, c(1, 5, 7, 8)]
  milogFCs$full <- rownames(milogFCs)
  rownames(milogFCs) <- NULL
  
  # Merge miRNA logFCs with pairs
  Network2 <- merge(miGePairs, milogFCs, by="node1", all.x =TRUE, all.y=FALSE)
  Network2$dir.y <- NULL
  Network2$pairs <- paste0(Network2$node1,"|",Network2$node2)
  
  # Average duplicated IDs
  Network3 <- data.table(Network2[,c(6:7,9)])
  Network3 <- Network3[, lapply(.SD, mean), by=pairs]
  Network3 <- data.frame(Network3)
  
  # Final merge and correted names
  miGePairs$pairs <- paste0(miGePairs$node1,"|",miGePairs$node2)
  Network4 <- merge(miGePairs, Network3, by="pairs", all.x = TRUE, all.y = FALSE)
  colnames(Network4) <- c("pairs", "node1", "node2", "dirnode1",  "dirnode2", "dir", "logFCnode1", "fdrnode1")
  
  # mRNAsets
  
  # Extract logFCs for mRNAs
  node1 <- GelogFCs$GeneName
  GelogFCs <- GelogFCs[, c(13, 17, 19)]
  GelogFCs$node1 <- node1
  
  Network5 <- merge(Network4, GelogFCs, by="node1", all.x =TRUE, all.y=FALSE)
  Network5$logFCnode1 <- ifelse(is.na(Network5$logFCnode1), Network5$logFC, Network5$logFCnode1)
  Network5$fdrnode1 <- ifelse(is.na(Network5$fdrnode1), Network5$adj.P.Val, Network5$fdrnode1)
  Network5 <- Network5[, c(1:8)]
  
  colnames(GelogFCs)[4] <- "node2"
  Network6 <- merge(Network5, GelogFCs, by="node2", all.x =TRUE, all.y=FALSE)
  Network6 <- Network6[,c(2,1,4,5,6,7,9,8,10,3)]
  colnames(Network6) <- c("node1", "node2", "dirnode1", "dirnode2", "dir", "logFCnode1", "logFCnode2", "fdrnode1", "fdrnode2",  "pairs")
  
  dups <- which(duplicated(Network6$pairs))
  if(length(dups) > 0) {
    NWlogFCs <- Network6[-dups,] 
  } else {
    NWlogFCs <- Network6
  }
  
  # Make nodeattribute file
  nodeattributes <- data.frame(c(as.character(NWlogFCs$node1), as.character(NWlogFCs$node2)), c(as.character(NWlogFCs$dirnode1), as.character(NWlogFCs$dirnode2)), c(NWlogFCs$logFCnode1, NWlogFCs$logFCnode2), c(NWlogFCs$fdrnode1, NWlogFCs$fdrnode2))
  colnames(nodeattributes) <- c("node","dir","logFC","fdr")
  nodeattributes <- unique(nodeattributes)
  nodeattributes$AbslogFC <- abs(nodeattributes$logFC)
  nodeattributes$Log2IvFdr <- log2(1/nodeattributes$fdr)
  
  return(list(NWlogFCs, nodeattributes))
}


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Archplot
# Takes arguments:
        # p.nodes = Output of function DENetWorkFC above
        # p.width = width of .tif plot
        # p.height = height of ,tif plot
        # p.name = Name of .tif plot
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


Arcplot <- function(p.nodes, p.width, p.height, p.name) {
  p.info <- data.table(c(as.character(p.nodes$node1), as.character(p.nodes$node2)), c(p.nodes$logFCnode1, p.nodes$logFCnode2))
  colnames(p.info) <- c("name", "logFC")
  p.info<- data.frame(p.info[, .(Freq = .N), by = .(name, logFC)])
  p.info <- p.info[order(p.info$Freq, decreasing = TRUE),]
  p.info$group <- 1:nrow(p.info)
  vircls <- c("#C6614D","#484A47")
  #vircls <- viridis(2, end = 0.6, direction = -1, option = "magma")
  p.info$calfb <- ifelse(p.info$logFC > 0, vircls[1], "grey50")
  p.info$calfb <- ifelse(p.info$logFC < 0, vircls[2], as.character(p.info$calfb))
  
  myorder <- unique(as.vector(t(data.frame(p.nodes[,1:2]))))
  p.info <- p.info[match(myorder, p.info$name),]
  
  final.nodes <- as.matrix(p.nodes[, 1:2])
  degrees <- as.numeric(p.info$Freq)
  names(degrees) <- p.info$name
  meta <- data.frame(p.info$group, degrees, p.info$name, ind=1:nrow(p.info))
  my.order <- as.integer(meta[order(meta$degrees, decreasing = TRUE),]$ind)
  
  tiff(paste0(p.name,"_Interacions.tiff"), width = p.width, height = p.height, units = 'in', res = 300)
  arcplot(final.nodes, ordering=my.order, labels=p.info$name, cex.labels=0.8,
          show.nodes=TRUE, col.nodes=p.info$calfb, bg.nodes=p.info$calfb,
          cex.nodes = log2(abs(p.info$logFC)), pch.nodes=21,
          lwd.nodes = 2, line=-0.5,
          col.arcs = hsv(0, 0, 0.2, 0.25), lwd.arcs = log2(p.nodes$score/100)*1.5)
  dev.off()
}



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extract unique and common miRNAs and mRNAs from two or more networks
# Takes as arguments;
        # List of node df output as the second element from DENetWorkFC
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ExtractmiRmRNA <- function(list.of.df) {
  lname <- names(list.of.df)
  miRlist <- list()
  mRlist <- list()
  for (idx in 1:length(list.of.df)) {
    df <- list.of.df[[idx]]
    mir <- as.character(df[grep("miR|let", df$node),]$node)
    mr <- as.character(df[-grep("miR|let", df$node),]$node)
    miRlist[[idx]] <- mir
    mRlist[[idx]] <- mr
  }
  names(miRlist) <- lname
  names(mRlist) <- lname
  l <- list(miRlist , mRlist)
  names(l) <- c("miRNA", "mRNA")
  return(l)
}



UniquemiRmRNA <- function(list.of.df) {
  lname <- c(names(list.of.df), "Common")
  miRlist <- list()
  mRlist <- list()
  UmiRs <- list()
  UmRs <- list()
  for (idx in 1:length(list.of.df)) {
    df <- list.of.df[[idx]]
    mir <- as.character(df[grep("miR|let", df$node),]$node)
    mr <- as.character(df[-grep("miR|let", df$node),]$node)
    miRlist[[idx]] <- mir
    mRlist[[idx]] <- mr
  }
  for (idx in 1:length(miRlist)) {
    other <- unique(as.character(unlist(miRlist[-idx])))
    Umir <- setdiff(miRlist[[idx]], other)
    UmiRs[[idx]] <- Umir
  }
  UmiRs[[length(UmiRs)+1]] <- Reduce(intersect, miRlist)
  names(UmiRs) <- lname
  
  for (idx in 1:length(mRlist)) {
    other <- unique(as.character(unlist(mRlist[-idx])))
    Umr <- setdiff(mRlist[[idx]], other)
    UmRs[[idx]] <- Umr
  }
  UmRs[[length(UmRs)+1]] <- Reduce(intersect, mRlist)
  names(UmRs) <- lname
  
  l <- list(UmiRs, UmRs)
  names(l) <- c("miRNAs", "mRNAs")
  return(l)
}




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extract miRNA-gene pairs based on genes in module.
# Takes as arguments;
    # my.network = Network of miRNA-gene and gene-gene pairs, must have atleast two columns (node1 and node2)
    # my.geneset = vector of gene names 
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#BestPairs <- function(my.network, my.geneset) {
#  WGCNA <- my.network[my.network$node1 %in%  my.geneset | my.network$node2 %in%  my.geneset ,]
#  WGCNAmiR <- WGCNA[grep("hsa|let", WGCNA$node1),]
#  WGCNAmiR <- data.frame(setDT(WGCNAmiR[,1:2])[, list(id=paste(node2, collapse=",")), by = node1])
#  colnames(WGCNAmiR) <- c("miRNA", "Gene")
#  return(WGCNAmiR)
#}

BestPairs <- function(my.network, my.geneset) {
  WGCNA <- my.network[my.network$node1 %in%  my.geneset | my.network$node2 %in%  my.geneset ,]
  WGCNAmiR <- WGCNA[grep("hsa|let", WGCNA$node1),1:4]
  colnames(WGCNAmiR) <- c("miRNA", "Gene", "Dir miRNA", "Dir Gene")
  return(WGCNAmiR)
}

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plor co-expressed genes from modules with miRNA interaction partner
# Takes as arguments;
    # my.f = dataframe where column 1 = genes, 2 = miRNA, 3 = direction of either miRNA OR Gene, 4 = miRNA family and 5 = module colors
    # my.set = name of set.
    # my.h, my.w =  integers specifying hight and width of plot.
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


InformationPlot <- function(my.df, my.modname, my.set,  my.h, my.w) {
  colnames(my.df) <- c("Gene", "miRNA", "Dir", "Family", "Module")
  my.df$Dir <- factor(my.df$Dir, levels = c("up", "down"))
  my.df <- my.df[with(my.df, order(Dir,Family)),]
  my.df$miRNA <- paste0(as.character(my.df$miRNA), " (", as.character(my.df$Family), ")")
  my.df$miRNA <- factor(as.character(my.df$miRNA), levels =  c(unique(as.character(my.df$miRNA))))
  if ("up" %in% my.df$Dir & "down" %in% my.df$Dir) {
    cols <- c("#FF7777", "#242F40")
  } else if ("up" %in% my.df$Dir) {
    cols <- c("#FF7777")
  } else {
    cols <-  c("#242F40")
  }
  m <- ggplot(my.df, aes(x=Gene, y=miRNA,  fill=Dir, colour="white")) + geom_tile(colour="white") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = cols) + labs(title=paste0(my.modname, " module")) + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size=.1, color="black")) 
  pdf(paste0(my.set, ".pdf"), height = my.h, width = my.w)
  m
  dev.off()
}

