# miRNA-mRNA-TIF

'''Created by Thilde Bagger Terkelsen, 27-06-2019, thilde@cancer.dk, thildebate@gmail.com'''
Computational Biology Laboratory, Danish Cancer Society Research Center, Strandboulevarden 49, 2100, Copenhagen, Denmark.

Repository associated with the publication:

###

Corresponding author: Elena Papaleo, elenap@cancer.dk

R-scripts were created by: Thilde Bagger Terkelsen, thilde@cancer.dk

This repository contains tumour interstitial fluid miRNA expression data and paired solid tumour tissue mRNA expression data from patients with breast cancer. The repository was made with intent of openly sharing both data and R-scripts used for analysis in relation to the publication.

The repository contains X folders:

    (1) Data: miRNA and mRNA expression datasets and non-sensitive patient metadata. 
    (2) Backgrounds_and_Databases: Databases and files used for the analysis. 
    (3) R-scripts: A collection of R-scripts that recapitulate the work.
                                

Requirements:

    R version 3.5.1 or higher
    Rstudio version 1.1.463 or higher        

Although R-packages should automatically be installed and errors raised if they cannot be, we here provide the user with the list of required packages for manual installation:

CRAN:

    caTools
    data.table
    dendextend
    factoextra
    ggplot2
    glmnet
    heatmap.plus
    pamr
    RColorBrewer
    randomForest
    reshape
    rms
    scales
    varSelRF
   
   
   
Bioconductor:   
    
    arcdiagram
    biomaRt
    limma
    sva
    UpSetR
    WGCNA



