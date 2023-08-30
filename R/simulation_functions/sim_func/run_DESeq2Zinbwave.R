###########################
# DESeq2 Default Pipeline #
###########################

#####################################
# Install or Load Required Packages #
#####################################

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('dplyr', 'reshape2')
if(! require("DESeq2")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
  }
library(DESeq2)
library(zinbwave)
library(SummarizedExperiment)
library(dplyr)
library(readr)
library(tibble)
library(BiocParallel)
library(scran)
library(GUniFrac)
###########################
# Fit DESeq2ZinbwaveRobust To A Dataset #
###########################

fit.DESeq2Zinbwave<-function(features, metadata, libSize, ID, transformation, MultTestCorrection, weights){
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default DESeq2 model. Use NONE.')
  
  
  # Random Effect Adjustment
  if(!length(ID)==length(unique(ID))){
    stop('edgeR random effect model is currently not implemented.')
  }
  
  formula <- as.formula(paste('~', paste(colnames(metadata), collapse = "+"), sep=''))
  
  
  # Fit Model
  x <- DESeqDataSetFromMatrix(countData = t(features), colData = metadata, design = formula)
  # gm_mean = function(x, na.rm=TRUE){
  #   exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  # }
  # geoMeans = apply(counts(x), 1, gm_mean)
  # x = estimateSizeFactors(x, geoMeans = geoMeans)
  # fit <- DESeq(x)
  
  #scrnew=GMPR(OTUmatrix=t(features), min_ct = 2, intersect_no = 4)
  #sizeFactors(x) <- scrnew
  
  x <- estimateSizeFactors(x, type="poscounts")
  
  #scrnew=GMPR(OTUmatrix=t(features)) #, min_ct = 2, intersect_no = 4)
  #sizeFactors(x) <- scrnew
  
  #scr <- computeSumFactors(x)         
        ## Warning in FUN(...): encountered negative size factor estimates
  #sizeFactors(x) <- sizeFactors(scr)
    #Note that we are returned a warning that some of the values are negative. See the description of this issue in the scran manual for computeSumFactors. Typically, I have found that more extensive filtering resolves this issue. Here we will proceed with the estimated positive values, but aware of the limitation/warning.
  # dds_zinb <- DESeq(dds_zinb, test="LRT", reduced=~1, sfType="poscounts",
  #                   minmu=1e-6, minReplicatesForReplace=3, fitType = "local") #minReplicatesForReplace=Inf
  
  # countsp=as.data.frame(features)
  # #dim(counts)
  # pdata = metadata
  # countsp$Run=row.names(countsp)
  # pdata$Run=row.names(countsp)
  # row.names(countsp) <- NULL
  # row.names(pdata) <- NULL
  # pdata <- column_to_rownames(pdata, "Run")
  # countsp <- column_to_rownames(countsp, "Run")
  # physeqSE <- SummarizedExperiment(assays=SimpleList(counts=t(countsp)), colData = pdata)
  # #physeqf
  # zinbmodel <- zinbFit(physeqSE, 
  #                      K = 0,
  #                      epsilon = 1e10, commondispersion = TRUE, verbose = FALSE, BPPARAM = SerialParam())
  # weights <- computeExactWeights(model = zinbmodel,x = t(features))
  # 
  weights[which(weights<1e-6)] <- 1e-06
  assays(x)[["weights"]] = weights
  
  fit<- DESeq(x, test="LRT", reduced=~1, #sfType="poscounts", #useT=TRUE,
              minmu=1e-6, minReplicatesForReplace=Inf) #, fitType = "local") #minReplicatesForReplace=Inf,
  
  # fit<- DESeq(x, test="Wald", sfType="poscounts",
  #             minmu=1e-6, minReplicatesForReplace=Inf, fitType = "local") #
  # 
  #res_zinb <- results(dds_zinb, alpha = 0.05, cooksCutoff = TRUE,independentFiltering=FALSE)
  # coef<-coef(fit)[,-1]
  # pval<-results(fit,  name=resultsNames(fit)[2])$pvalue  #alpha = 0.05, cooksCutoff = FALSE,
  # DD<-cbind.data.frame(coef,pval)
  # DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
  # DD<-DD[order(DD$qval, decreasing=FALSE),]
  # DD = subset(DD, qval < 0.05)
  # DeSeq2_Zinbwave = rownames(DD)
  # 
  
  
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(coef(fit)[,-1], 'coef')
    pvalMatrix<-get_pval_DESeq2(fit)
    pvalue.vector<-rename.features(pvalMatrix[,-1], 'pval')
    DD<-cbind.data.frame(coef.vector, pvalue.vector)
    DD<-DD[, !duplicated(colnames(DD))]
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD$relevance<- (DD$qval < 0.05)*1
    #DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
    rownames(DD)<-NULL;
  }
  else{
    coef<-coef(fit)[,-1]
    pval<-results(fit,  name=resultsNames(fit)[2], independentFiltering=FALSE)$pvalue #, alpha = 0.05, cooksCutoff = FALSE)$pvalue  #
    DD<-cbind.data.frame(coef,pval)
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD$relevance<- (DD$qval < 0.05)*1
    #DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD$feature<-rownames(DD)
    DD$metadata<- names(metadata)
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
    rownames(DD)<-NULL;
  }
  return(DD)
}

# Get P-values from DESeq2 Fit
get_pval_DESeq2<-function(fit){
  List <- list()
  for(i in 1:length(resultsNames(fit))){
    List[[i]] <- results(fit,name=resultsNames(fit)[i])$pvalue   #,alpha = 0.05, cooksCutoff = FALSE
  }
  Matrix = do.call(cbind, List)
  rownames(Matrix)<-names(fit)
  colnames(Matrix)<-resultsNames(fit)
  return(Matrix)
}

computeExactWeights <- function (model, x) 
{
  mu <- getMu(model)
  pi <- getPi(model)
  theta <- getTheta(model)
  theta <- matrix(rep(theta, each = ncol(x)), ncol = nrow(x))
  nb_part <- dnbinom(t(x), size = theta, mu = mu)
  zinb_part <- pi * ( t(x) == 0 ) + (1 - pi) *  nb_part
  zinbwg <- ( (1 - pi) * nb_part ) / zinb_part
  zinbwg <- t(zinbwg)
  zinbwg[x > 0] <- 1
  zinbwg[zinbwg < 1e-15] <- 1e-15
  zinbwg
}
####################################
# Fit DESeq2 To A List of Datasets #
####################################

list.DESeq2Zinbwave<-function(physeq, transformation='NONE', MultTestCorrection = 'BH'){
  foreach(physeq=physeq, 
          .export=c("fit.DESeq2Zinbwave", "rename.features", "get_pval_DESeq2",
                    "computeExactWeights"),
          .packages = c("DESeq2", "dplyr", "reshape2", "zinbwave", "SummarizedExperiment",
                        "readr", "tibble", "BiocParallel", "scran","GUniFrac"), 
          .errorhandling = 'remove') %dopar% 
          {
            start.time <- Sys.time()
            features<-physeq$features; 
            metadata<-physeq$metadata;
            libSize <- physeq$libSize;
            weights <- physeq$weights;
            ID<-physeq$ID;
            DD<-fit.DESeq2Zinbwave(features, metadata, libSize, ID, transformation, MultTestCorrection, weights)
            DD$pairwiseAssociation<-paste('pairwiseAssociation', 1:nrow(DD), sep='')
            wh.TP = intersect(grep("[[:print:]]+\\_TP$", DD$metadata), grep("[[:print:]]+\\_TP$", DD$feature))
            newname = paste0(DD$pairwiseAssociation[wh.TP], "_TP")
            DD$pairwiseAssociation[wh.TP] <- newname;
            DD<-dplyr::select(DD, c('pairwiseAssociation', 'feature', 'metadata'), everything())
            stop.time <- Sys.time()
            time<-as.numeric(round(difftime(stop.time, start.time, units="min"),3), units = "mins")
            DD$time<-time
            return(DD)
          }
}

#xy=list.DESeq2ZinbwaveRobust(fitlist, transformation='NONE', MultTestCorrection = 'BH')
#physeqf=fitlist[[2]]
# features<-physeqf$features; 
# metadata<-physeqf$metadata;
# libSize <- physeqf$libSize;
# ID<-physeqf$ID;
#physeq=fitlist

