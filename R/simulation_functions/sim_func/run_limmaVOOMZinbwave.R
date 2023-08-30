##############################
# limmaVOOM Default Pipeline #
##############################

#####################################
# Install or Load Required Packages #
#####################################

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('dplyr', 'reshape2')
if(! require("limma")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("limma")
}
library(limma)
library(edgeR)
##############################
# Fit limmaVOOM To A Dataset #
##############################

fit.limmaVOOMZinbwave = function(features, metadata, libSize, ID, transformation, MultTestCorrection, weights) {
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default limmaVOOM model. Use NONE.')
  
  # Fit limmaVOOM
   x<-t(as.matrix(features)+1) # Convert to matrix, round up to nearest integer, and transpose
  design <- model.matrix(~., data=metadata)
  
  #counts2p=t(features)
  #counts2p=t(countsp)  #simulation based
  #dge <- DGEList(counts2p, group=metadata$Metadata1_TP) 
  #dge <- edgeR::calcNormFactors(dge)
  #dge$samples$norm.factors <- NFs
  #design = as.formula("~ SolA_log")
  #design <- model.matrix(design, data.frame(sample_data(physeqt)))
  #v <- voom(dge, design, plot=FALSE, lib.size = dge$samples[,2]*dge$samples[,3])
  #v <- voom(dge, design, plot = FALSE, weights = weights,lib.size = dge$samples[,2]*dge$samples[,3])
  #v$weights <- v$weights * weights
  #fit <- lmFit(v, design, weights = v$weights , method="robust")  #observation outliers
  
  #fit <- lmFit(v,design)
  #fit <- eBayes(fit, robust=TRUE)  ##dispersion outliers
  #tt <- topTable(fit, coef = 2, n = nrow(dge), sort.by="p")
  
  y <- voom(x,design,plot=FALSE)
  
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
  # weights[which(weights<1e-6)] <- 1e-06
  
  
  y$weights <- y$weights * weights
  
  # Fit limma 
  # Random Effect Adjustment
  if (!length(ID)==length(unique(ID))){
    dupcor <-  limma::duplicateCorrelation(v, design, block = ID)
    fit <- limma::lmFit(y,design, weights = y$weights, block = ID, correlation = dupcor$cor)
  } else{ 
    fit <- limma::lmFit(y,design, weights = y$weights)}
  
  # Empirical Bayes Adjustment
  fit <- limma::eBayes(fit)
  
  # coef<-fit$coefficients[,-1]
  # pval<-fit$p.value[,-1]
  # DD<-cbind.data.frame(coef,pval)
  # DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
  # DD<-DD[order(DD$qval, decreasing=FALSE),]
  # DD = subset(DD, qval < 0.05)
  # limmavoom_Zinbwave = rownames(DD)
  
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(fit$coefficients[,-1], 'coef')
    pvalue.vector<-rename.features(fit$p.value[,-1], 'pval')
    DD<-cbind.data.frame(coef.vector, pvalue.vector)
    DD<-DD[, !duplicated(colnames(DD))]
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD$relevance<- (DD$qval < 0.05)*1
    #DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
    rownames(DD)<-NULL;
  } 
  else{
    coef<-fit$coefficients[,-1]
    pval<-fit$p.value[,-1]
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
#######################################
# Fit limmaVOOM To A List of Datasets #
#######################################

list.limmaVOOMZinbwave<-function(physeq, transformation='NONE', MultTestCorrection= 'BH'){
  foreach(physeq=physeq, 
          .export=c("fit.limmaVOOMZinbwave", "rename.features","computeExactWeights"),
          .packages = c("limma", "dplyr", "reshape2","edgeR", "zinbwave", "SummarizedExperiment",
                        "readr", "tibble", "BiocParallel"), 
          .errorhandling = 'remove') %dopar% 
          {
            start.time <- Sys.time()
            features<-physeq$features; 
            metadata<-physeq$metadata;
            libSize <- physeq$libSize;
            weights <- physeq$weights;
            ID<-physeq$ID;
            DD<-fit.limmaVOOMZinbwave(features, metadata, libSize, ID, transformation, MultTestCorrection, weights)
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

#warnings()
