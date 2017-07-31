################################################################
#                                                              #
#      BENCHMARK SCRIPT FOR CLUSTERING ALGORITHMS              #
#                                                              #
################################################################

## Packages to load for PCAREDUCE ##
library(SC3)
library(pcaReduce)
library(pcaMethods)
library(scater)
library(mclust)
## PCAREDUCE ##
library(Rtsne)
library(SINCERA)
library(Seurat)
library(rPython)
library(pheatmap)
library(reshape2)
library(e1071)
library(dbscan)
library(NMF)
library(ggplot2)
library(RColorBrewer)

##############
# PARAMETERS #
##############

options(expressions = 5e5)
nCores = 3

#########
# UTILS #
#########

## Functions to load for PCAREDUCE ##
#' Get gene expression matrix from dataset -sceset
#'
#' Gets the gene expression matrix from the proper slot of the dataset -sceset.
#'
#' @param dataset input SCESet 
#' @return gene expression matrix
#' 
#' @importFrom scater counts
#' @importFrom scater tpm
#' @importFrom scater fpkm
#' @importFrom scater cpm
#' 
#' @export
getData <- function(dataset) {
  if (!is.null(counts(dataset))) return(counts(dataset))
  if (!is.null(tpm(dataset))) return(tpm(dataset))
  if (!is.null(fpkm(dataset))) return(fpkm(dataset))
  if (!is.null(cpm(dataset))) return(cpm(dataset))
}

#' Get "reference" clustering from dataset
#'
#' Gets the "reference" clustering from the input dataset.
#'
#' @param dataset input SCESet 
#' @param nbhier such as "cell_type@nbhier" is the slot where the reference clustering is contained
#' @return the "reference clustering" from the dataset
#' 
#' @importFrom scater pData
#' 
#' @export
getClustering <- function(dataset, nbhier) {
  ctypes <- pData(dataset)[,sprintf("cell_type%d", nbhier)]
  hir <- unique(ctypes)
  clusters <- data.frame(1:length(hir), row.names = hir)
  cl <- sapply(ctypes, function(celltype) clusters[celltype, 1])
  return(cl)
}
## PCAREDUCE ##

#' Import dataset from CSV file
#'
#' Imports the dataset -features as rows, conditions as columns- from the 
#' specified CSV file name, contained in @pathFile.
#'
#' @param dataset name of the CSV file
#' @param removeAllNA boolean set to TRUE iff all missing values should be removed
#' @param pathFile absolute path to the CSV file
#' @return the dataset contained in the input CSV file
#' 
#' @importFrom utils read.csv
#' 
#' @export
importDataset <- function(dataset, removeAllNA = FALSE, pathFile="~/Documents/dataaa/") {
  data <- read.csv(paste0(pathFile, dataset, ".csv", sep = ""), header=TRUE, na.strings=-1)
  checkDuplicated <- duplicated(data[,1])
  indices <- c(unlist(sapply(1:length(checkDuplicated), function(i) if (!checkDuplicated[i]) i)))
  data <- data[indices,]
  if (removeAllNA) data <- data[apply(data, 1, function(y) !all(is.na(y))),]
  row.names(data) <- data[,1]
  data <- data[,2:dim(data)[2]]
  if (dataset == "yan" | dataset == "treutlein" | dataset == "deng") data <- as.data.frame(t(data.matrix(data)))
  filterSpikeIn <- row.names(data)[!grepl("ERCC", row.names(data))]
  return(data[filterSpikeIn,])
}

#' Construct a SCESet object from a gene expression matrix
#'
#' Constructs a SCESet object from a gene expression matrix.
#'
#' @param dataset gene expression matrix
#' @param unitcount unit for the counts in @dataset: RPKM, TPM, ...
#' @return the SCESet associated with the input matrix
#' 
#' @importFrom scater newSCESet
#' @importFrom scater calculateQCMetrics
#' @importFrom methods new
#' 
#' @export
turn_into_sceset <- function(dataset, unitcount) {
  varMetadataT <- data.frame(labelDescription = NA, row.names = "cell")
  phenodata <- new("AnnotatedDataFrame", data = data.frame(cell = colnames(dataset)), 
                   varMetadata = varMetadataT)
  row.names(phenodata) <- colnames(dataset)
  if (unitcount == "fpkm" | unitcount == "rpkm") sceset <- newSCESet(fpkmData = data.matrix(dataset), 
                                                                     phenoData = phenodata, 
                                                                     logExprsOffset = 1)
  if (unitcount == "tpm") sceset <- newSCESet(tpmData = data.matrix(dataset), phenoData = phenodata)
  if (unitcount == "cpm") sceset <- newSCESet(cpmData = data.matrix(dataset), phenoData = phenodata)
  else sceset <- newSCESet(countData = data.matrix(dataset), phenoData = phenodata)
  sceset <- calculateQCMetrics(sceset)
  return(sceset)
}

##################
# TEST FUNCTIONS #
##################

## These functions write in slots (for each algorithm and number of clusters) of the input SCESet ##
## the runtime, the parameter values and the clustering labels                                    ##

#_______________#
#      SC3      #
#_______________#

#ref: https://github.com/hemberg-lab/SC3/blob/master/vignettes/my-vignette.Rmd
#ref: https://hemberg-lab.github.io/scRNA.seq.course/clust-methods.html#sc3-1
#' Test function for SC3
#'
#' Applies SC3 algorithm to the dataset.
#'
#' @param dataset input SCESet
#' @param nbclusters expected number of clusters
#' @return the input SCESet with the newly-written slots
#' 
#' @importFrom SC3 sc3
#' @importFrom SC3 sc3_prepare
#' @importFrom SC3 sc3_estimate_k
#' @importFrom scater pData
#' 
#' @export
applySC3 <- function(dataset, nbclusters) {
  start <- Sys.time()
  dataset <- sc3_prepare(dataset, ks = 2:length(sampleNames(dataset)), n_cores = nCores)
  dataset <- sc3_estimate_k(dataset)
  end <- Sys.time()
  timeEst <- end-start
  k <- if (is.null(nbclusters)) dataset@sc3$k_estimation else nbclusters
  start <- Sys.time()
  dataset <- sc3(dataset, ks = k, biology = TRUE, 
                 gene_filter = TRUE, pct_dropout_min = 10,
                 pct_dropout_max = 90, d_region_min = 0.04, d_region_max = 0.07,
                 svm_num_cells = NULL, svm_train_inds = NULL, svm_max = 5000,
                 n_cores = nCores, kmeans_nstart = NULL, kmeans_iter_max = 1e+09,
                 k_estimator = FALSE, rand_seed = 1)
  end <- Sys.time()
  if (is.null(nbclusters)) {
    pData(dataset)$runtime_sc3_est <- as.character(end-start)
    pData(dataset)$timeEst_sc3 <- as.character(timeEst)
    pData(dataset)$parameters_sc3_est <- list(k)
  }
  else {
    pData(dataset)$runtime_sc3 <- as.character(end-start)
    pData(dataset)$parameters_sc3 <- list(k)
  }
  pData(dataset)$clusters <- pData(dataset)[, sprintf("sc3_%d_clusters", k)]
  return(dataset)
}

#_____________________#
#      PCAREDUCE      #
#_____________________#

#ref: https://hemberg-lab.github.io/scRNA.seq.course/clust-methods.html#sc3-1
#' Test function for pcaReduce
#'
#' Applies pcaReduce algorithm to the dataset.
#'
#' @param dataset input SCESet
#' @param nbclusters expected number of clusters
#' @return the input SCESet with the newly-written slots
#' 
#' @importFrom pcaReduce PCAreduce
#' @importFrom scater exprs
#' @importFrom scater pData
#' 
#' @export
applyPCAREDUCE <- function(dataset, nbclusters, nbt = 100) {
  input <- exprs(dataset)
  start <- Sys.time()
  ## "run pcaReduce nbt = 'no of samples' times creating hierarchies from 1 to q clusters" ##
  ## PCAreduce centers the data matrix                                                     ##
  pcared <- PCAreduce(t(input), nbt = nbt, q = nbclusters, method = 'S')
  end <- Sys.time()
  pData(dataset)[, paste0("pcareduce_", nbclusters, "_clusters", sep = "")] <- as.character(pcared[[1]][,"Cl_mat"])
  pData(dataset)$runtime_pcareduce <- as.character(end-start)
  pData(dataset)$parameters_pcareduce <- list(list(nbclusters, nbt))
  pData(dataset)$clusters <- as.character(pcared[[1]][,"Cl_mat"])
  return(dataset)
}

#________________________#
#      tSNE+K-means      #
#________________________#

#' Test function for tSNE+K-means
#'
#' Applies tSNE+K-means algorithm to the dataset.
#'
#' @param dataset input SCESet
#' @param nbclusters expected number of clusters
#' @param perplexity perplexity parameter for tSNE, between 2 and min 50, no. samples
#' @param theta theta parameter for tSNE, between 0 and 1
#' @param itmax maximum number of iterations for tSNE
#' @return the input SCESet with the newly-written slots
#' 
#' @importFrom Rtsne Rtsne
#' @importFrom stats kmeans
#' @importFrom scater pData
#' 
#' @export
applyTSNE <- function(dataset, nbclusters, perplexity = 50, theta = 1, itmax = 1000) {
  ## samples in rows, dimensions in columns ##
  data <- getData(dataset)
  counts <- t(data.matrix(data))
  start <- Sys.time()
  ## 2d dimension reduction ##
  reddata <- Rtsne(counts, perplexity = perplexity, theta = theta, max_iter = itmax,
                   initial_dims = dim(data)[2], check_duplicates = FALSE, 
                   verbose = T, 
                   dims = 2, pca = TRUE, is_distance = FALSE,
                   pca_center = TRUE, pca_scale = FALSE, 
                   momentum = 0.5, final_momentum = 0.8, eta = 200, exaggeration_factor = 12)
  end <- Sys.time()
  timeEst <- end-start
  start <- Sys.time()
  clustObject <- kmeans(reddata$Y, centers = nbclusters, 
                        iter.max = 10, nstart = 1, algorithm = "Hartigan-Wong", trace=FALSE)
  end <- Sys.time()
  pData(dataset)$runtime_tSNE_kmeans <- as.character(end-start)
  pData(dataset)$timeEst_tSNE_kmeans <- as.character(timeEst)
  pData(dataset)[, paste0("tSNE_kmeans_", nbclusters, "_clusters", sep = "")] <- as.character(clustObject$cluster)
  pData(dataset)$parameters_tSNE_kmeans <- list(list(nbclusters, perplexity, theta, itmax))
  pData(dataset)$clusters <- as.character(clustObject$cluster)
  return(dataset)
}

#_______________________#
#      tSNE+DBSCAN      #
#_______________________#

#' Test function for tSNE+DBSCAN
#'
#' Applies tSNE+DBSCAN algorithm to the dataset.
#'
#' @param dataset input SCESet
#' @param nbclusters expected number of clusters
#' @param perplexity perplexity parameter for tSNE, between 2 and min 50, no. samples
#' @param theta theta parameter for tSNE, between 0 and 1
#' @param itmax maximum number of iterations for tSNE
#' @param eps epsilon parameter for DBSCAN
#' @return the input SCESet with the newly-written slots
#' 
#' @importFrom Rtsne Rtsne
#' @importFrom DBSCAN dbscan
#' @importFrom scater pData
#' 
#' @export
applyDBSCAN <- function(dataset, nbclusters, perplexity = 34, theta = 0.5, itmax = 1000, eps = 1) {
  ## samples in rows, dimensions in columns ##
  data <- getData(dataset)
  counts <- t(data.matrix(data))
  start <- Sys.time()
  ## 2d dimension reduction ##
  reddata <- Rtsne(counts, perplexity = perplexity, theta = theta, max_iter = itmax,
                   initial_dims = dim(data)[2], check_duplicates = FALSE, 
                   verbose = T, 
                   dims = 2, pca = TRUE, is_distance = FALSE,
                   pca_center = TRUE, pca_scale = FALSE,
                   momentum = 0.5,
                   final_momentum = 0.8, eta = 200, exaggeration_factor = 12)
  end <- Sys.time()
  pData(dataset)$timeEst_dbscan <- as.character(end-start)
  start <- Sys.time()
  clustObject <- dbscan(reddata$Y, eps = eps, minPts = 5, weights = NULL, borderPoints = TRUE)
  end <- Sys.time()
  pData(dataset)$runtime_tSNE_dbscan <- as.character(end-start)
  pData(dataset)[, paste0("tSNE_dbscan_", nbclusters, "_clusters", sep = "")] <- as.character(clustObject$cluster)
  pData(dataset)$parameters_tSNE_dbscan <- list(list(perplexity, theta, itmax, eps))
  pData(dataset)$clusters <- as.character(clustObject$cluster)
  return(dataset)
}

#___________________#
#      K-means      #
#___________________#

#' Test function for K-means
#'
#' Applies K-means algorithm to the dataset.
#'
#' @param dataset input SCESet
#' @param nbclusters expected number of clusters
#' @param itmax maximum number of iterations for K-means
#' @param nstart number of random sets that should be chosen at first
#' @return the input SCESet with the newly-written slots
#' 
#' @importFrom stats kmeans
#' @importFrom scater pData
#' @importFrom scater exprs
#' 
#' @export
applyKMEANS <- function(dataset, nbclusters, itmax = 10000, nstart = 1) {
  data <- t(scale(exprs(dataset))[,])
  start <- Sys.time()
  kmeansObject <- kmeans(data, centers = nbclusters, 
                         iter.max = itmax, nstart = nstart, algorithm = "Hartigan-Wong")
  end <- Sys.time()
  pData(dataset)$runtime_kmeans <- as.character(end-start)
  pData(dataset)[, paste0("kmeans_", nbclusters, "_clusters", sep = "")] <- as.character(kmeansObject$cluster)
  pData(dataset)$parameters_kmeans <- list(list(nbclusters, itmax, nstart))
  pData(dataset)$clusters <- as.character(kmeansObject$cluster)
  return(dataset)
}

#____________________#
#      SNN-Cliq      #
#____________________#

#ref: https://hemberg-lab.github.io/scRNA.seq.course/clust-methods.html#sc3-1
#ref: http://bioinfo.uncc.edu/SNNCliq/#manual
#' Test function for SNN-Cliq
#'
#' Applies SNN-Cliq algorithm to the dataset.
#'
#' @param dataset input SCESet
#' @param nbclusters expected number of clusters
#' @param r 0 <= r < 1, "density threshold of quasi-cliques"
#' @param m 0 <= m < 1, "threshold on the overlapping rate for merging"
#' @return the input SCESet with the newly-written slots
#' 
#' @importFrom scater pData
#' @importFrom utils read.table
#' 
#' @export
applySNNCLIQ <- function(dataset, nbclusters, r = 0.8, m = 0.8) {
  folderSNNCLIQ = "./"
  source(paste(folderSNNCLIQ, "snn.R", sep = ""))
  start <- Sys.time()
  data <- getData(dataset)
  SNN(data = t((data)), outfile="snn-cliq.txt", k=nbclusters, distance="euclidean")
  snn.res <- system(paste0("python ", folderSNNCLIQ, "Cliq.py ", 
                           "-i snn-cliq.txt ","-o res-snn-cliq.txt ",
                           "-r ", r," -m ", m), intern = TRUE)
  end <- Sys.time()
  snn.res <- read.table("res-snn-cliq.txt")
  system("rm snn-cliq.txt res-snn-cliq.txt")
  pData(dataset)$runtime_SNNCliq <- as.character(end-start)
  pData(dataset)[, paste0("SNNCliq_", nbclusters, "_clusters", sep = "")] <- as.character(snn.res[,1])
  pData(dataset)$parameters_SNNCliq <- list(list(nbclusters, r, m))
  pData(dataset)$clusters <- as.character(snn.res[,1])
  return(dataset)
}

#__________________#
#      SINCERA     #
#__________________#

#ref: https://github.com/minzheguo/SINCERA/blob/master/demo/humanIPF.R
#ref: https://hemberg-lab.github.io/scRNA.seq.course/clust-methods.html
#' Test function for SINCERA
#'
#' Applies SINCERA algorithm to the dataset.
#'
#' @param dataset input SCESet
#' @param nbclusters expected number of clusters
#' @return the input SCESet with the newly-written slots
#' 
#' @importFrom SINCERA construct
#' @importFrom SINCERA normalization.zscore
#' @importFrom SINCERA getGenesForClustering
#' @importFrom SINCERA getCellMeta
#' @importFrom SINCERA cluster.assignment
#' @importFrom stats cutree
#' @importFrom stats cor 
#' @importFrom stats hclust
#' @importFrom cluster clusGap
#' @importFrom scater pData
#' 
#' @export
applySINCERA <- function(dataset, nbclusters) {
  data <- as.data.frame(getData(dataset))
  data <- data[apply(data, 1, function(x) !all(x==0 | is.na(x)) & !(var(x) == 0)),]
  data <- data[apply(data, 2, function(x) !(var(x) == 0)),]
  sc <- construct(exprmatrix=data, samplevector=colnames(data))
  sc <- normalization.zscore(sc, pergroup=FALSE)
  cordist <- function(y) as.dist((1-cor(t(y), method = "pearson"))/2)
  if (is.null(nbclusters)) {
    hclustForGap <- function(y, k) list(cluster=cutree(hclust(cordist(y), method="average"), k=k))
    x <- getExpression(sc, scaled=T, genes=getGenesForClustering(sc))
    start <- Sys.time()
    gapstats <- clusGap(t(x), FUN=hclustForGap, K.max=dim(data)[2])
    end <- Sys.time()
    pData(dataset)$timeEst_sincera <- as.character(end-start)
    pData(dataset)$gapstats_sincera <- list(gapstats)
    plot(1:dim(data)[2], gapstats$Tab[,3])
    #Should restart the algorithm with the elbow point
    return(dataset)
  }
  start <- Sys.time()
  sc <- cluster.assignment(sc, feature.type = "gene", k=nbclusters, verbose = F)
  end <- Sys.time()
  pData(dataset)$runtime_sincera <- as.character(end-start)
  pData(dataset)$parameters_sincera <- list(nbclusters)
  pData(dataset)[, paste0("sincera_", nbclusters, "_clusters", sep = "")] <- getCellMeta(sc, "CLUSTER")
  pData(dataset)$clusters <- getCellMeta(sc, "CLUSTER")
  return(dataset)
}

#_________________#
#      SEURAT     #
#_________________#

#ref: https://hemberg-lab.github.io/scRNA.seq.course/clust-methods.html#sc3-1
#ref: http://satijalab.org/seurat/pbmc-tutorial.html
#' Test function for SEURAT
#'
#' Applies SEURAT algorithm to the dataset.
#'
#' @param dataset input SCESet
#' @param nbclusters expected number of clusters
#' @param resolution SEURAT parameter: between 0.6 and 1.2 according to the authors
#' @return the input SCESet with the newly-written slots
#' 
#' @importFrom seurat Setup
#' @importFrom seurat MeanVarPlot
#' @importFrom seurat RegressOut
#' @importFrom seurat PCAFast
#' @importFrom seurat RunTSNE
#' @importFrom seurat FindClusters
#' @importFrom scater pData
#' 
#' @export
applySEURAT <- function(dataset, nbclusters, resolution = 1.34) {
  data <- getData(dataset)
  data <- data[apply(data, 1, function(x) !all(x==0 | is.na(x))),]
  seuratdata <- Setup(new("seurat", raw.data = data), project = "Data")
  start <- Sys.time()
  seuratdata <- MeanVarPlot(seuratdata, do.plot = FALSE)
  seuratdata <- RegressOut(seuratdata, latent.vars = c("nUMI"), genes.regress = seuratdata@var.genes)
  seuratdata <- PCAFast(seuratdata)
  seuratdata <- RunTSNE(seuratdata)
  seuratdata <- FindClusters(seuratdata, resolution = resolution)
  end <- Sys.time()
  samples <- colnames(seuratdata@data)
  v <- rep(0, dim(dataset)[2])
  names(v) <- colnames(dataset)
  v[samples] <- seuratdata@ident
  pData(dataset)$runtime_seurat <- as.character(end-start)
  pData(dataset)$parameters_seurat <- list(resolution)
  pData(dataset)[, paste0("seurat_", nbclusters, "_clusters", sep = "")] <- v
  pData(dataset)$clusters <- v
  return(dataset)
}

#______________#
#      NMF     #
#______________#

#ref: https://github.com/naikai/sake/blob/master/R/nmf_utils.R
#ref: https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf
#ref: https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btw607
#ref: http://nimfa.biolab.si
#' Test function for NMF
#'
#' Applies NMF algorithm -using NIMFA initialization- to the dataset.
#'
#' @param dataset input SCESet
#' @param nbclusters expected number of clusters
#' @param maxiter maximum number of iterations for NMF
#' @param nrun number of runs for NMF
#' @return the input SCESet with the newly-written slots
#' 
#' @importFrom NMF nmf
#' @importFrom scater pData
#' 
#' @export
applyNMF <- function(dataset, nbclusters, maxIter = 5000, nrun = 100) {
  data <- getData(dataset)
  data <- data[apply(data, 1, function(x) !all(x==0 | is.na(x))),]
  start <- Sys.time()
  res <- nmf(data, rank = nbclusters, method = "snmf/r", 
             nrun = nrun, seed = "nndsvd")
  end <- Sys.time()
  cluster <- sapply(1:dim(data)[2], function(i) which.max(coef(res)[,i]))
  pData(dataset)$runtime_nmf <- as.character(end-start)
  pData(dataset)$parameters_nmf <- list(list(nbclusters, maxIter, nrun, "nndsvd", "snmf/r"))
  pData(dataset)[, paste0("nmf_", nbclusters, "_clusters", sep = "")] <- as.character(cluster)
  pData(dataset)$clusters <- as.character(cluster)
  return(dataset)
}

#___________________#
#      F-CMEANS     #
#___________________#

#ref: https://cran.r-project.org/web/packages/e1071/e1071.pdf
#' Test function for Fuzzy C-means
#'
#' Applies NMF algorithm -using NIMFA initialization- to the dataset.
#'
#' @param dataset input SCESet
#' @param nbclusters expected number of clusters
#' @param itmax maximum number of iterations for Fuzzy C-means
#' @param m fuzzication degree
#' @return the input SCESet with the newly-written slots
#' 
#' @importFrom e1071 cmeans
#' @importFrom scater exprs
#' @importFrom scater pData
#' 
#' @export
applyFKMEANS <- function(dataset, nbclusters, itmax = 28, m = 18) {
  data <- t(scale(exprs(dataset))[,])
  start <- Sys.time()
  cl <- cmeans(data, centers = nbclusters, 
               iter.max = itmax, verbose = TRUE, method = "cmeans", m = m)
  end <- Sys.time()
  pData(dataset)$runtime_fkmeans <- as.character(end-start)
  pData(dataset)$parameters_fkmeans <- list(list(nbclusters, itmax, m))
  pData(dataset)[, paste0("fkmeans_", nbclusters, "_clusters", sep = "")] <- as.character(cl$cluster)
  pData(dataset)$clusters <- as.character(cl$cluster)
  return(dataset)
}

#_______________#
#      MINE     #
#_______________#

#' Test function for Mine
#'
#' Applies Mine algorithm to the dataset.
#'
#' @param dataset input SCESet
#' @param nbclusters expected number of clusters
#' @param d probability threshold for cell similarity
#' @param f frequency gene trimming criterion
#' @param K number of neighbors that should agree to the merge of two clusters to make it happen
#' @param normalized boolean set to TRUE iff the gene expression matrix is already normalized
#' @return the input SCESet with the newly-written slots
#' 
#' @importFrom scater pData
#' 
#' @export
applyMINE <- function(dataset, nbclusters, d = 0.5, f = 0.7, K = NULL, normalized = T) {
  data <- getData(dataset)
  start <- Sys.time()
  clusters <- mine(data, d = d, f = f, K = K, normalized = normalized)
  end <- Sys.time()
  pData(dataset)$runtime_mine <- as.character(end-start)
  pData(dataset)[, paste0("mine_", nbclusters, "_clusters", sep = "")] <- as.character(clusters) ##TODO
  pData(dataset)$parameters_mine <- list(list(d, f, K))
  pData(dataset)$clusters <- as.character(clusters)
  return(dataset)
}

###################
# FUNCTION ARRAYS #
###################

#Cleaning: R.utils::gcDLLs() 

#_________________________#
#      FUNCTION ARRAY     #
#_________________________#

## To load for PCAREDUCE ##
a <- list(applySC3, 
  #applyPCAREDUCE#, 
  applyNMF,
  applyTSNE, 
  applyKMEANS, 
  applySNNCLIQ,
  applySINCERA, 
  applySEURAT, 
  #applyMINE, 
  applyDBSCAN, applyFKMEANS
)
names(a) <- c("sc3", 
  #"pcareduce"#, 
  "nmf", 
  "tSNE_kmeans", 
  "kmeans", 
  "SNNCliq",
  "sincera", 
  "seurat", 
  #"mine",
  "tSNE_dbscan", "fkmeans"
)
## PCAREDUCE ##

#________________________#
#      DATASET ARRAY     #
#________________________#

unitcounts <- list("fpkm", "rpkm",
                   "rawcounts", "rawcounts", "fpkm",
                   #"rawcounts", "tpm", 
                   "rpkm", #"rawcounts",
                   #"rawcounts", 
                   #"?", 
                   "rawcounts", "rawcounts"#, 
                   #"rawcounts", "rawcounts", "rawcounts"
                   )
names(unitcounts) <- c("biase", "yan",
                       "goolam", "deng", "treutlein",
                       #"ting", "patel", 
                       "usoskin", #"klein",
                       #"zeisel", 
                       #"macosko", 
                       "tintori", "ciona"#, 
                       #"green1", "green2", "pollen"
                       )

#__________________________#
#      SCESet FUNCTION     #
#__________________________#

#' Get SCESet from CSV file
#'
#' Gets the SCESet object associated with the specified CSV file.
#'
#' @param dataname name of the CSV file
#' @return the associated SCESet
#' 
#' @export
getSceset <- function(dataname) 
  return(turn_into_sceset(importDataset(dataname, removeAllNA = F), unitcounts[[dataname]]))

############################
# ITERATE CLUSTERING TESTS #
############################

#' Iterate all algorithms on dataset
#'
#' Iterates all the algorithms on the input dataset.
#'
#' @param dataset the SCESet object associated to the dataset
#' @param nbclusters the expected number of clusters
#' @param nbhier the number of the slot "cell_type" where the reference clustering can be found
#' @return the SCESet with newly-written slots, for each algorithm
#' 
#' @export
allIterate <- function(dataset, nbclusters, nbhier = 1) {
  dataname <- deparse(substitute(dataset))
  for (name in names(a)) {
    ## Seurat, SNNCliq and SINCERA are deterministic ##
    if (!(name == "sincera" | name == "SNNCliq" | name == "seurat")) {
      dataset <- iterateAlgo(name, dataset, nbclusters, nbhier = nbhier)
    }
    else dataset <- iterateAlgo(name, dataset, nbclusters, nbhier = nbhier, maxit = 1)
    save(dataset, file = paste0("clust_", dataname, ".Rdata", sep = ""))
  }
  return(dataset)
}

## Function to load for PCAREDUCE ##
#' Iterate one algorithm on dataset
#'
#' Iterates the specified algorithm on the input dataset.
#'
#' @param name string for the name of the algorithm -matching those in names@a
#' @param dataset the SCESet object associated to the dataset
#' @param nbclusters the expected number of clusters
#' @param nbhier the number of the slot "cell_type" where the reference clustering can be found
#' @param maxit the number of iterations for the algorithm
#' @param ref_cl reference clustering label vector, if there is none in the input SCESet
#' @return the SCESet with newly-written slots
#' 
#' @importFrom scater pData
#' @importFrom mclust adjustedRandIndex
#' 
#' @export
iterateAlgo <- function(name, dataset, nbclusters, nbhier = 1, maxit = 100, ref_cl = NULL) {
  field <- if (is.null(nbclusters)) sprintf("ari_%s", name) else sprintf("ari_%s_%dcl", name, nbclusters)
  ari <- tryCatch({pData(dataset)[, field]}, error=function(e){NULL})
  if (is.null(ari)) pData(dataset)[, field] <- -Inf
  algo <- a[[name]]
  if (is.null(pData(dataset)$ref_cl) & is.null(ref_cl)) pData(dataset)$ref_cl <- getClustering(dataset, nbhier)
  if (is.null(pData(dataset)$ref_cl)) pData(dataset)$ref_cl <- ref_cl
  i <- 0
  while(i < maxit) {
    tmp <- algo(dataset, nbclusters)
    nari <- adjustedRandIndex(pData(tmp)$clusters, pData(tmp)$ref_cl)
    
    print(paste0(name, " -- #iteration : ", i, " -- ari : ", nari))
    
    if (nari > pData(dataset)[1, field]) {
      dataset <- tmp
      pData(dataset)[, field] <- nari
    }
    i <- i + 1
  }
  return(dataset)
}
## PCAREDUCE ##

#########
# PLOTS #
#########

#_______________________#
#      COMPUTATIONS     #
#_______________________#

#' Load datasets
#'
#' Loads datasets.
#'
#' @param path the absolute path where the datasets can be found
#' @return a list of the loaded datasets
#' 
#' @export
loadDatasets <- function(path="~/Documents/benchmark/datasets/") {
  dNames <- c("biase", "goolam", "deng", #"treutlein", "klein", 
              #"ciona",
              "tintori", "kolodziejczyk", #"usoskin", 
              "yan")
  ld <- function(name) return(paste(path, "clust_", name, ".Rdata", sep=""))
  for (d in dNames) load(ld(d))
  datasets <- c(biase, goolam, deng, #treutlein,
                #klein,
                #ciona,
                tintori, 
                kolodziejczyk, #usoskin, 
                yan)
  names(datasets) <- dNames
  return(datasets)
}

#' Create summary of benchmark and store it in a file
#'
#' Creates a summary of benchmark and stores it in a file called plots.Rdata.
#'
#' @param datasets list of datasets
#' @param path the absolute path where the datasets can be found
#' @return a data frame summarizing the benchmark results
#' 
#' @importFrom scater pData
#' 
#' @export
scriptPlots <- function(datasets=NULL, path="~/Documents/benchmark/") {
  if (is.null(datasets)) datasets <- loadDatasets()
  dNames <- names(datasets)
  kNumbers <- c(6, 8, 9, #12, 16, 
                #5, 
                5, #11, 
                11, 7)
  algos <- c("sc3", "pcareduce", "tSNE_kmeans", "tSNE_dbscan", "kmeans", 
             "SNNCliq", 
             "sincera", "seurat",
             "nmf",
             "fkmeans")
  getRuntime <- function(algoName) {
    if (algoName == "tSNE_kmeans") return(unlist(lapply(datasets, function(d) as.numeric(pData(d)$runtime_tSNE_kmeans[1])
                                                        +as.numeric(pData(d)$timeEst_tSNE_kmeans[1]))))
    if (algoName == "tSNE_dbscan") return(unlist(lapply(datasets, function(d) as.numeric(pData(d)$runtime_tSNE_dbscan[1])
                                                        +as.numeric(pData(d)$timeEst_dbscan[1]))))
    return(unlist(lapply(datasets, function(d) as.numeric(pData(d)[, sprintf("runtime_%s", algoName)][1]))))
  }
  for (n in algos) assign(paste(n, "RUNTIME", sep=""), getRuntime(n))
  getARI <- function(algoName) {
    print(algoName)
    return(unlist(lapply(1:length(datasets), 
                         function(i) as.numeric(pData(datasets[[i]])[, sprintf("ari_%s_%dcl", algoName, kNumbers[i])][[1]][1]))))
  }
  for (n in algos) assign(paste(n, "ARI", sep=""), getARI(n))
  df <- data.frame(pcareduceRUNTIME=pcareduceRUNTIME, tSNE_kmeansRUNTIME=tSNE_kmeansRUNTIME, 
                   tSNE_dbscanRUNTIME=tSNE_dbscanRUNTIME, kmeansRUNTIME=kmeansRUNTIME, SNNCliqRUNTIME=SNNCliqRUNTIME,
                   sc3RUNTIME=sc3RUNTIME, 
                   sinceraRUNTIME=sinceraRUNTIME, 
                   seuratRUNTIME=seuratRUNTIME, 
                   nmfRUNTIME=nmfRUNTIME, 
                   fkmeansRUNTIME=fkmeansRUNTIME,
                   pcareduceARI=pcareduceARI, tSNE_kmeansARI=tSNE_kmeansARI, 
                   tSNE_dbscanARI=tSNE_dbscanARI, kmeansARI=kmeansARI, SNNCliqARI=SNNCliqARI,
                   sc3ARI=sc3ARI, 
                   sinceraARI=sinceraARI, 
                   seuratARI=seuratARI, 
                   nmfARI=nmfARI, 
                   fkmeansARI=fkmeansARI)
  getStability <- function(algoName)
    return(c(as.numeric(pData(datasets[["biase"]])[, sprintf("ariDataset_%s", algoName)][[1]]), 
             as.numeric(pData(datasets[["biase"]])[, sprintf("stability_%s", algoName)][1])))
  for (n in algos) if (!(n == "SNNCliq" | n == "sincera" | n == "seurat")) assign(paste("ariDataset", n, sep=""), getStability(n))
  ddf <- data.frame(ariDatasetsc3=ariDatasetsc3, 
                    ariDatasetpcareduce=ariDatasetpcareduce,
                    ariDatasettSNE_kmeans=ariDatasettSNE_kmeans,
                    ariDatasettSNE_dbscan=ariDatasettSNE_dbscan,
                    ariDatasetkmeans=ariDatasetkmeans,
                    #ariDatasetnmf=ariDatasetnmf,
                    ariDatasetfkmeans=ariDatasetfkmeans)
  nbCells <- lapply(datasets, function(d) dim(d)[2])
  nbGenes <- lapply(datasets, function(d) dim(d)[1])
  df$nbCells <- nbCells
  df$nbGenes <- nbGenes
  df$maxiter <- 100
  df$kNumbers <- kNumbers
  save(df, ddf, file=paste(path, "plots.Rdata", sep = ""))
  return(df)
}

#________________#
#      PLOTS     #
#________________#

#' Post-clustering plots
#'
#' Plots different useful plots -PCA, tSNE, histograms, ...- to analyze the results of the benchmark.
#'
#' @param plotName plot to draw
#' @param algoName specify a single name
#' @return a plot object
#' 
#' @importFrom scater plotPCA
#' @importFrom scater plotTSNE
#' @importFrom RColorBrewer palette
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 qplot
#' @importFrom graphics axis
#' @importFrom graphics plot
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics abline
#' @importFrom graphics hist
#' @importFrom graphics points
#' @importFrom stats median
#' @mportFrom stats mean
#' 
#' @export
library(RColorBrewer)
library(ggplot2)
plots <- function(plotName, algoName=NULL, datasets=NULL, path="~/Documents/benchmark/") {
  load(paste(path, "plots.Rdata", sep=""))
  if (is.null(df)) df <- scriptPlots()
  algos <- c("sc3", "pcareduce", "tSNE_kmeans", "tSNE_dbscan", "kmeans", 
             "SNNCliq", 
             "sincera", "seurat", 
             "nmf", 
             "fkmeans")
  colours = palette(rainbow(max(length(algos), length(datasets))))
  prettyAlgos <- c("SC3", "PCAReduce", "+Kmeans", "+DBSCAN", "K-means", 
                   "SNN-Cliq", 
                   "SINCERA", 
                   "Seurat", 
                   "NMF",
                   "Fuzzy C-means")
  names(prettyAlgos) <- algos
  prettyUndetAlgos <- prettyAlgos[!(prettyAlgos == "SINCERA" | prettyAlgos == "Seurat" | prettyAlgos == "SNN-Cliq")]
  undetAlgos <- algos[!(algos == "sincera" | algos == "seurat" | algos == "SNNCliq")]
  if (is.null(datasets)) datasets <- loadDatasets()
  ncells <- unlist(df$nbCells)
  ngenes <- unlist(df$nbGenes)
  
  switch(plotName,
         
         "time one cell gene" = {
           if (is.null(algoName)) return(NULL)
           main = paste("Runtime=f(#cells, #genes): ", prettyAlgos[algoName], sep = "")
           ylim = c(min(ncells), max(ncells)+1e2)
           ylim2 = c(min(ngenes), max(ngenes)+1e4)
           y = df[[paste(algoName,"RUNTIME",sep="")]]
           par(mar = c(5,5,2,5))
           plot(y, ncells, main=main, type="p", pch = 0, col="red3", xlab = "Runtime in seconds", ylab="#cells", ylim=ylim)
           abline(v=mean(y), col="red", lwd = 3)
           par(new = T)
           plot(y, ngenes, pch=16, type = "p", axes=F, xlab = NA, ylab=NA, ylim=ylim2, cex=1.2)
           axis(side = 4)
           mtext(side = 4, line = 3, '#genes')
           legend("topleft", legend=c("time=f(#cells)", "time=f(#genes)", sprintf("mean runtime: %2.2f sec.", mean(y))), 
                  lty=c(0,0,1), pch=c(0, 16, NA), col=c("red3", "black", "red")
           )
         },
         
         "time all cell" = {
           main = "Runtime=f(#cells): all algorithms"
           y = df[, paste(algos[1], "RUNTIME", sep="")]
           my <- c(mean(y))
           gi = 1:length(datasets)
           xlim = c(min(ncells), max(ncells)+1e2)
           ylim = c(min(y), max(y)+10)
           plot(ncells, y, main=main, type="p", col=colours[1], xlab = "#cells", xaxt='n', ylab = "Runtime in seconds", xlim=xlim, ylim=ylim)
           abline(h=mean(y), col=colours[1], lwd = 3)
           for (i in 2:length(algos)) {
             y = df[, paste(algos[i], "RUNTIME", sep="")]
             points(ncells, y, col=colours[i])
             my <- c(my, mean(y))
             abline(h=my[i], col=colours[i], lwd = 3)
           }
           legend("topleft", legend = sapply(1:length(algos), function(i) paste0(prettyAlgos[i], ": ", round(my[i], 2), "s")), fill = colours)
           axis(side=1, at=unique(ncells), labels=sapply(gi, function(i) paste(names(datasets)[i], "\n", ncells[i], sep="")))
         },
         
         "time all gene" = {
           main = "Runtime=f(#genes): all algorithms"
           y = df[, paste(algos[1], "RUNTIME", sep="")]
           my <- c(mean(y))
           gi = 1:length(datasets)
           xlim = c(min(ngenes), max(ngenes)+1e4)
           ylim = c(min(y), max(y)+10)
           plot(ngenes, y, main=main, type="p", col=colours[1], xlab = "#genes", xaxt='n', ylab = "Runtime in seconds", xlim=xlim, ylim=ylim)
           abline(v=mean(y), col=colours[1], lwd = 3)
           for (i in 2:length(algos)) {
             y = df[, paste(algos[i], "RUNTIME", sep="")]
             points(ngenes, y, col=colours[i])
             my <- c(my, mean(y))
             abline(h=mean(y), col=colours[i], lwd = 3)
           }
           legend("topleft", legend = sapply(1:length(algos), function(i) paste0(prettyAlgos[i], ": ", round(my[i], 2), "s")), fill = colours)
           axis(side=1, at=unique(ngenes), labels=sapply(gi, function(i) paste(names(datasets)[i], "\n", ngenes[i], sep="")))
         },
         
         "ari average all" = {
           aries <- sapply(algos, function(name) mean(df[, paste(name, "ARI", sep="")]))
           plot(1:length(algos), aries, main = "Average ARI value: all algorithms", type="h", col=colours,
                xlab = "Algorithm", ylab = "Average ARI", ylim=c(0,1), xaxt='n', lwd=30)
           lines(1:length(algos), rep(0.5, length(algos)), col="red")
           legend("topleft", legend = "ARI=0.5", lty = 1, pch = NA, fill = "red")
           axis(side=1, at=unique(1:length(algos)), labels=prettyAlgos)
         },
         
         "ari average all datasets" = {
           aries <- sapply(1:length(datasets), function(i) mean(as.numeric(df[names(datasets)[i], grepl(".*ARI$", colnames(df))])))
           plot(1:length(datasets), aries, main = "Average ARI value: all datasets", type="h", col=colours, 
                xlab = "Dataset", ylab = "Average ARI", ylim=c(min(min(aries), 0),1), xaxt='n', lwd=30)
           lines(1:length(datasets), rep(0.5, length(datasets)), col="red")
           legend("topleft", legend = "ARI=0.5", lty = 1, pch = NA, fill = "red")
           axis(side=1, at=unique(1:length(datasets)), labels=names(datasets))
         },
         
         "ari all" = {
           Algorithms <- NULL
           for (name in prettyAlgos) Algorithms <- c(Algorithms, rep(name, length(datasets)))
           Datasets=rep(names(datasets), length(algos))
           ARI <- NULL
           for (name in algos) ARI <- c(ARI, df[, paste(name, "ARI", sep="")])
           dat=data.frame(Algorithms, Datasets, ARI)
           ggplot(dat, aes(Algorithms, ARI, fill=Datasets)) + geom_bar(position="dodge", stat="identity") + ggtitle("Adjusted Rand Index")
         },
         
         "ari one stability hist" = {
           if (is.null(algoName)) return(NULL)
           y <- ddf[1:100, paste("ariDataset", algoName, sep="")]
           hist(y, col="grey", main=sprintf("Histogram of the ARI for 100 iterations of %s on Biase dataset", 
                                            prettyAlgos[algoName]), xlab="ARI values", ylab="Frequency")
           abline(v=median(y), col="blue", lwd=3)
           legend("topleft", legend = "median ARI value", fill = "blue", lty = 1, pch = NA)
         },

         "stability one" = {
           stabilities <- as.numeric(ddf[101, ])
           plot(1:length(stabilities), stabilities, main = "Stability measure values: all algorithms on Biase dataset", type="h", col=colours, 
                xlab = "Algorithm", ylab = "Stability measure", ylim=c(0,1), xaxt='n', lwd=30)
           axis(side=1, at=unique(1:length(undetAlgos)), labels=prettyUndetAlgos)
         },
         
         "ari all boxplot" = {
           aries <- NULL
           for (name in undetAlgos) aries <- c(aries, ddf[1:100, paste("ariDataset", name, sep="")])
           algorithms <- NULL
           for (name in prettyUndetAlgos) algorithms <- c(algorithms, rep(name, 100))
           ariesData <- data.frame(algorithms, aries)
           qplot(x=algorithms, y=aries, data = ariesData, geom=c("boxplot"), fill=algorithms, 
                 main="Adjusted Rand Indices boxplot on Biase dataset", xlab = "Algorithm", ylab = "ARI")
         }
  )
}

###################
# STABILITY TESTS #
###################

#' Update solution stability vector
#'
#' Updates the solution count vector for stability tests according to the newly-computed clustering solution.
#'
#' @param solution newly-computed clustering solution
#' @param solutions vector containing the previously seen solutions
#' @param countSolutions solution count vector -of the number of time this solution appeared
#' @return a list containing the updated values of @solutions and @countSolutions
#' 
#' @importFrom mclust adjustedRandIndex
#' 
#' @export
updateCountSolution <- function(solution, solutions, countSolutions) {
  ind <- NULL
  if (!is.null(solutions)) ind <- unlist(sapply(1:length(solutions), function(i) if (adjustedRandIndex(solution, solutions[[i]])==1) return(i)))
  if (is.null(ind)) {
    solutions <- append(solutions, list(solution))
    countSolutions <- c(countSolutions, 1)
  }
  else countSolutions[ind] <- countSolutions[ind] + 1
  return(list(solutions, countSolutions))
}

#' Perform stability test
#'
#' Performs a stability test on a given algorithm, on a given dataset.
#'
#' @param name the name of the algorithm -should match those in names@a
#' @param dataset a SCESet object associated with the dataset
#' @param nbclusters the expected number of clusters
#' @param maxit the maximum number of iterations
#' 
#' @importFrom mclust adjustedRandIndex
#' @importFrom scater pData
#' 
#' @export
testStability <- function(name, dataset, nbclusters, maxit = 100, randDelete=F, delRate=0.2) {
  solutions <- NULL
  countSolutions <- NULL
  ariDataset <- NULL
  algo <- a[[name]]
  i <- 0
  while(i < maxit) {
    print(i)
    if (randDelete) {
      d <- deparse(substitute(dataset))
      d <- turn_into_sceset(getData(dataset), unitcounts[[d]])
      keep <- sample(1:dim(d)[2], floor((1-delRate)*dim(d)[2]))
      keep <- colnames(d)[keep]
      datasetT <- d[, keep] 
      nbclusters <- length(table(pData(dataset)$ref_cl))
    }
    else datasetT <- dataset
    tmp <- algo(datasetT, nbclusters)
    if (randDelete) {
      tmpp <- adjustedRandIndex(pData(tmp)[, paste0(name, "_", nbclusters, "_clusters", sep = "")], pData(tmp)$clusters)
      ariDataset <- c(ariDataset, tmpp)
    }
    else {
      ariDataset <- c(ariDataset, adjustedRandIndex(pData(tmp)$clusters, pData(tmp)$ref_cl))
      tmp <- updateCountSolution(pData(tmp)$clusters, solutions, countSolutions)
      solutions <- tmp[[1]]
      countSolutions <- tmp[[2]]
    }
    i <- i + 1
  }
  if (randDelete) pData(dataset)[, paste0("ariDataset_", name, "DEL", sep = "")] <- list(list(ariDataset))
  else {
    pData(dataset)[, paste0("stability_", name, sep="")] <- max(countSolutions)/maxit
    pData(dataset)[, paste0("ariDataset_", name, sep = "")] <- list(list(ariDataset))
  }
  return(dataset)
}

###############
# GET RESULTS #
###############

#' Return pairs of ARI, nb of clusters for each algorithm on a given dataset
#'
#' Returns pairs of ARI, nb of clusters for each algorithm on a given dataset.
#'
#' @param dataset a SCESet object associated with the dataset
#' @return pairs of ARI, nb of clusters for all algorithms on the given dataset
#' 
#' @importFrom mclust adjustedRandIndex
#' @importFrom scater pData
#' 
#' @export
readResults <- function(dataset) {
  nn1 <- sort(colnames(pData(dataset))[grepl("ari_", colnames(pData(dataset)))])
  print(nn1)
  l1 <- as.numeric(pData(dataset)[1, nn1])
  nn <- sort(colnames(pData(dataset))[grepl("_clusters", colnames(pData(dataset)))])
  l <- as.numeric(apply(pData(dataset)[, nn], 2, min))
  l2 <- as.numeric(apply(pData(dataset)[, nn], 2, max))
  l2 <- l2-l+1
  l3 <- sapply(1:length(l2), function(i) sprintf("(%2.2f, %1.0f)", l1[i], l2[i]))
  refNumber <- length(table(pData(dataset)$ref_cl))
  names(l3) <- nn
  print(refNumber)
  return(l3)
}

#' Return pairs of ARI, nb of clusters for each algorithm on all datasets
#'
#' Returns pairs of ARI, nb of clusters for each algorithm on all datasets.
#'
#' @param datasetsNames list of dataset names
#' @param datasets a SCESet object associated with the datasets in @datasetsNames
#' @return pairs of ARI, nb of clusters for all datasets
#' 
#' @export
readAllResults <- function(datasetsNames, datasets) 
  return(lapply(datasetsNames, function(datasetNm) {
    print(datasetNm)
    readResults(datasets[datasetNm])
    }))

###############
# RANKING     #
###############

#' Compute weights for the benchmark ranking
#'
#' Computes weights for the benchmark ranking.
#'
#' @param df result of function scriptPlots
#' @param ddf result of function scriptPlots
#' @return weights for the ranking
#' 
#' @importFrom stats median
#' 
#' @export
computeWeights <- function(df, ddf=NULL) {
  dnames <- row.names(df)
  weights <- sapply(dnames, function(dname) median(as.numeric(df[dname, colnames(df)[grepl(".*ARI$", colnames(df))]])))
  weights[weights < 0] <- 0
  weights <- weights/sum(weights)
  names(weights) <- dnames
  return(weights)
}

#' Compute weights for the benchmark ranking
#'
#' Computes weights for the benchmark ranking.
#'
#' @param df result of function scriptPlots
#' @param q quantity in c("ARI", "RUNTIME")
#' @param ddf result of function scriptPlots
#' @param datasets list of datasets
#' @return scores for all algorithms for quantity @q
#' 
#' @export
computeScores <- function(df, q, ddf=NULL, datasets=NULL, undetAlgos = NULL) {
  weights <- computeWeights(df)
  algos <- c("sc3", "pcareduce", "tSNE_kmeans", "tSNE_dbscan", "kmeans", 
             #"SNNCliq", 
             "sincera", "seurat", 
             #"nmf", 
             "fkmeans")
  if (q %in% c("ARI", "RUNTIME")) values <- lapply(algos, function(algo) df[, paste0(algo, q)])
  if (q == "STABILITY") values <- as.numeric(ddf[101, colnames(ddf)[grepl("ariDataset.*", colnames(ddf))]])
  else {
    values <- rep(list(df[, "kNumbers"]), length(algos))
    v <- lapply(algos, function(a) as.numeric(sapply(datasets, 
                                          function(d) 
                                            length(
                                                unique(
                                                  pData(d)[, 
                                                           colnames(pData(d))[grepl(paste0(a, "_.*_clusters$"), colnames(pData(d)))]
                                                           ]
                                                  )
                                              )
                                          )
    )
                )
    values <- lapply(1:length(algos), function(i) abs(values[[i]]-v[[i]]))
  }
  if (q != "STABILITY") {
    values <- sapply(1:length(values), function(i) weights %*% values[[i]])
    names(values) <- algos
  }
  else names(values) <- undetAlgos
  return(values)
}

#' Compute ranking for the benchmark ranking
#'
#' Computes ranking for the benchmark ranking.
#'
#' @param df result of function scriptPlots
#' @param q quantity in c("ARI", "RUNTIME")
#' @param ddf result of scriptPlots
#' @param datasets list of datasets
#' @param increasing is q == "RUNTIME" or "NUMBER"?
#' @return ranking of algorithms, from the best one, respect to q, to the worst one, with values
#' 
#' @export
computeRanking <- function(df, q, ddf=NULL, datasets=NULL, increasing=F, undetAlgos=NULL) {
  values <- sort(computeScores(df, q, ddf, datasets, undetAlgos=undetAlgos), decreasing=!increasing)
  return(values)
}