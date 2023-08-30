##########################
# edgeR Default Pipeline #
##########################

#####################################
# Install or Load Required Packages #
#####################################

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('dplyr', 'reshape2')
if(! require("edgeR")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")
}
library(edgeR)
library(zinbwave)
library(SummarizedExperiment)
library(dplyr)
library(readr)
library(tibble)
library(BiocParallel)

##########################
# Fit edgeR To A Dataset #
##########################

fit.edgeRZinbwave = function(features, metadata, libSize, ID, transformation, MultTestCorrection, weights) {
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default edgeR model. Use NONE.')
  
  d <- DGEList(counts=t(features))
  d <- edgeR::calcNormFactors(d, method='TMM')
  
  # Random Effect Adjustment
  if(!length(ID)==length(unique(ID))){
    stop('edgeR weighted for random effect model is currently not implemented.')
  }
  design <- model.matrix(~., data=metadata)
  
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
  d$weights <- weights
  d <- estimateDisp(d, design)
  # d <- estimateGLMTrendedDisp(d,design)
  # d <- estimateGLMTagwiseDisp(d,design)
  fit <- glmFit(d, design)
  
  # fit<- glmWeightedF(fit, 2) #, ZI = TRUE, independentFiltering = FALSE)
  # coef<-fit$coefficients[,-1]
  # pval<-fit$table$PValue
  # DD<-cbind.data.frame(coef,pval)
  # DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
  # DD<-DD[order(DD$qval, decreasing=FALSE),]
  # DD = subset(DD, qval < 0.05)
  # edgeR_Zinbwave = rownames(DD)
  
  
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(fit$coefficients[,-1], 'coef')
    pvalMatrix<-get_pval_edgeR(fit)
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
    fit<- glmWeightedF(fit, 2) #, ZI = TRUE, independentFiltering = FALSE)
    coef<-fit$coefficients[,-1]
    pval<-fit$table$PValue
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


#f0=rownames(DD)
#f1=cbind(DD$feature, DD$relevance)

# Get P-values from An edgeR Fit

get_pval_edgeR<-function(fit){
  mat<-as.matrix(fit$coefficients)
  rownames<-rownames(mat)
  colnames<-colnames(mat)
  n<-dim(fit$coefficients)[2]
  List <- list()
  for(i in 1:n){
    List[[i]] <-  glmWeightedF(fit, i)$table$PValue
  }
  Matrix = do.call(cbind, List)
  rownames(Matrix)<-rownames
  colnames(Matrix)<-colnames
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

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


glmWeightedF2 <- 
  function (glmfit=glmfit, coef = ncol(glmfit$design), contrast = NULL, 
            ZI = TRUE, independentFiltering = FALSE, filter = NULL) 
  {
    if (!is(glmfit, "DGEGLM")) {
      if (is(glmfit, "DGEList") && is(coef, "DGEGLM")) {
        stop("First argument is no longer required. Rerun with just the glmfit and coef/contrast arguments.")
      }
      stop("glmfit must be an DGEGLM object (usually produced by glmFit).")
    }
    if (is.null(glmfit$AveLogCPM)) 
      glmfit$AveLogCPM <- aveLogCPM(glmfit)
    nlibs <- ncol(glmfit)
    design <- as.matrix(glmfit$design)
    nbeta <- ncol(design)
    if (nbeta < 2) 
      stop("Need at least two columns for design, usually the first is the intercept column")
    coef.names <- colnames(design)
    if (is.null(contrast)) {
      if (length(coef) > 1) 
        coef <- unique(coef)
      if (is.character(coef)) {
        check.coef <- coef %in% colnames(design)
        if (any(!check.coef)) 
          stop("One or more named coef arguments do not match a column of the design matrix.")
        coef.name <- coef
        coef <- match(coef, colnames(design))
      }
      else coef.name <- coef.names[coef]
      logFC <- glmfit$coefficients[, coef, drop = FALSE]/log(2)
    }
    else {
      contrast <- as.matrix(contrast)
      qrc <- qr(contrast)
      ncontrasts <- qrc$rank
      if (ncontrasts == 0) 
        stop("contrasts are all zero")
      coef <- 1:ncontrasts
      if (ncontrasts < ncol(contrast)) 
        contrast <- contrast[, qrc$pivot[coef]]
      logFC <- drop((glmfit$coefficients %*% contrast)/log(2))
      if (ncontrasts > 1) {
        coef.name <- paste("LR test of", ncontrasts, 
                           "contrasts")
      }
      else {
        contrast <- drop(contrast)
        i <- contrast != 0
        coef.name <- paste(paste(contrast[i], coef.names[i], 
                                 sep = "*"), collapse = " ")
      }
      Dvec <- rep.int(1, nlibs)
      Dvec[coef] <- diag(qrc$qr)[coef]
      Q <- qr.Q(qrc, complete = TRUE, Dvec = Dvec)
      design <- design %*% Q
    }
    if (length(coef) == 1) 
      logFC <- as.vector(logFC)
    design0 <- design[, -coef, drop = FALSE]
    fit.null <- glmFit(glmfit$counts, design = design0, offset = glmfit$offset, 
                       weights = glmfit$weights, dispersion = glmfit$dispersion, 
                       prior.count = 0)
    LR <- fit.null$deviance - glmfit$deviance
    if (ZI) 
      fit.null$df.residual <- rowSums(fit.null$weights) - ncol(design0) #
    if (ZI) 
      glmfit$df.residual <- rowSums(glmfit$weights) - ncol(design)  #
    df.test <- fit.null$df.residual - glmfit$df.residual
    LRT.pvalue <- {
      phi <- quantile(glmfit$dispersion, p = 0.5)  #glmfit$dispersion #
      mu <- quantile(glmfit$fitted.values, p = 0.5) #rowMeans(glmfit$fitted.values) #
      gamma.prop <- (phi * mu/(1 + phi * mu))^2
      prior.df <- glmfit$prior.df
      if (is.null(prior.df)) 
        prior.df <- 20
      
      bartlett.corc = 1-(2*dim(glmfit$weights)[1]+5)/(6*(dim(glmfit$weights)[2]-1)) 
      glmfit$df.total <- glmfit$df.residual + prior.df/gamma.prop
      pf((LR/df.test), df1 = df.test, df2 = glmfit$df.total,   #0.25*glmfit$df.residual,   #glmfit$df.total, 
         lower.tail = FALSE, log.p = FALSE)
    }
    LRT.pvalue[glmfit$df.residual <= 0] = NA
    rn <- rownames(glmfit)
    if (is.null(rn)) 
      rn <- 1:nrow(glmfit)
    else rn <- make.unique(rn)
    tab <- data.frame(logFC = logFC, logCPM = glmfit$AveLogCPM, 
                      LR = LR, PValue = LRT.pvalue, row.names = rn)
    glmfit$counts <- NULL
    glmfit$table <- tab
    glmfit$comparison <- coef.name
    glmfit$df.test <- df.test
    res <- new("DGELRT", unclass(glmfit))
    if (independentFiltering) {
      if (is.null(filter)) 
        filter = rowMeans(glmfit$fitted.values)
      res <- independentFiltering(res, filter = filter, objectType = "edgeR")
    }
    else return(res)
  }
###################################
# Fit edgeR To A List of Datasets #
###################################

list.edgeRZinbwave<-function(physeq, transformation='NONE', MultTestCorrection = 'BH'){
  foreach(physeq=physeq, 
          .export=c("fit.edgeRZinbwave", "rename.features", "get_pval_edgeR","computeExactWeights","glmWeightedF2"),
          .packages = c("edgeR", "dplyr", "reshape2","zinbwave","SummarizedExperiment",
                        "readr","tibble","BiocParallel"), 
          .errorhandling = 'remove') %dopar% 
    {
      start.time <- Sys.time()
      features<-physeq$features; 
      metadata<-physeq$metadata;
      libSize <- physeq$libSize;
      weights <- physeq$weights;
      weights[which(weights<1e-6)] <- 1e-6
      ID<-physeq$ID;
      DD<-fit.edgeRZinbwave(features, metadata, libSize, ID, transformation, MultTestCorrection, weights)
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

#xy=list.edgeRZinbwave(fitlist, transformation='NONE', MultTestCorrection = 'BH')
#physeq=fitlist[[2]]
# features<-physeq$features; 
# metadata<-physeq$metadata;
# libSize <- physeq$libSize;
# ID<-physeq$ID;
# weights<-physeq$weights;
#physeq=fitlist


# a1=c('a', 'b', 'c','d', 'e', 'f')
# a2=c(1,0,1,1,0,0)
# 
# b1=c('a','c','f')
# b2=c(1,0,1)
# df1=data.frame(a1,a2)
# colnames(df1)=c("feat","rel")
# df2=data.frame(b1,b2)
# colnames(df2)=c("feat","rel")
# 
# df3 = merge(df1,df2, by="feat",all=TRUE)
# df3[is.na(df3)] <- 0
# 
# df4=df3[-1]
# df3$Rj=rowMeans(df4)
# Mallj=mean(Rj)
# df3$accept= (df3$Rj>Mallj)
# df3<-df3[order(df3$accept, decreasing=TRUE),]
# wh_predicted=subset(df3,accept=="TRUE")



