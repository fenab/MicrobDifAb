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
library(scran)

library(GUniFrac)
###########################
# Fit DESeq2 To A Dataset #
###########################

fit.DESeq2<-function(features, metadata, libSize, ID, transformation, MultTestCorrection){
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default DESeq2 model. Use NONE.')
  
  
  # Random Effect Adjustment
  if(!length(ID)==length(unique(ID))){
    stop('edgeR random effect model is currently not implemented.')
  }
  
  formula <- as.formula(paste('~', paste(colnames(metadata), collapse = "+"), sep=''))
  #sel=c("Feature30_TP","Feature92","Feature103","Feature109","Feature142","Feature165")
  #all=colnames(features)
  #subs=which(all %in% sel)
  #featuresR=features[,-c(21, 53, 63, 69, 79, 90)]
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
  
  #scr <- computeSumFactors(x)         
  # ## Warning in FUN(...): encountered negative size factor estimates
  #sizeFactors(x) <- sizeFactors(scr)
  #Note that we are returned a warning that some of the values are negative. See the description of this issue in the scran manual for computeSumFactors. Typically, I have found that more extensive filtering resolves this issue. Here we will proceed with the estimated positive values, but aware of the limitation/warning.
  # dds_zinb <- DESeq(dds_zinb, test="LRT", reduced=~1, sfType="poscounts",
  #                   minmu=1e-6, minReplicatesForReplace=3, fitType = "local") #minReplicatesForReplace=Inf
  
  fit<- DESeq(x, test="LRT", reduced=~1, sfType="poscounts", useT = TRUE,
                    minmu=1e-6, minReplicatesForReplace=Inf, fitType = "local") #

  # fit<- DESeq(x, test="Wald", sfType="poscounts", useT = TRUE,
  #             minmu=1e-6, minReplicatesForReplace=Inf, fitType = "local") #
  #res_zinb <- results(dds_zinb, alpha = 0.05, cooksCutoff = TRUE,independentFiltering=FALSE)
 
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
    pval<-results(fit, name=resultsNames(fit)[2])$pvalue  #alpha = 0.05, cooksCutoff = FALSE,
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
get_pval_DESeq2<-function(fit1){
  List <- list()
  for(i in 1:length(resultsNames(fit1))){
    List[[i]] <- results(fit1,name=resultsNames(fit1)[i])$pvalue #alpha = 0.05, cooksCutoff = FALSE,
  }
  Matrix = do.call(cbind, List)
  rownames(Matrix)<-names(fit1)
  colnames(Matrix)<-resultsNames(fit1)
  return(Matrix)
}

####################################
# Fit DESeq2 To A List of Datasets #
####################################

list.DESeq2<-function(physeq, transformation='NONE', MultTestCorrection = 'BH'){
  foreach(physeq=physeq, 
          .export=c("fit.DESeq2", "rename.features", "get_pval_DESeq2"),
          .packages = c("DESeq2", "dplyr", "reshape2", "scran","GUniFrac"), 
          .errorhandling = 'remove') %dopar% 
          {
            start.time <- Sys.time()
            features<-physeq$features; 
            metadata<-physeq$metadata;
            libSize <- physeq$libSize;
            ID<-physeq$ID;
            DD<-fit.DESeq2(features, metadata, libSize, ID, transformation, MultTestCorrection)
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

#xy=list.DESeq2(simlist, transformation='NONE', MultTestCorrection = 'BH')
#physeqf=fitlist[[2]]
# features<-physeqf$features; 
# metadata<-physeqf$metadata;
# libSize <- physeqf$libSize;
# ID<-physeqf$ID;
#physeq=simlist

