######################
## RFA Normalization #
######################
list.RFA<- function(physeq){
  foreach(physeq = physeq, .packages = "phyloseq") %dopar% {
    physeq$features = rarefy_even_depth(otu_table(physeq$features, taxa_are_rows = FALSE))
  return(physeq)
  }
}

######################
## CSS Normalization #
######################

# Apply CSS Normalization To A Dataset

CSSnorm = function(physeq) {

  # Extract ID and Metadata and Keep Them Unchanged
  ID<-physeq$ID
  metadata<-physeq$metadata

  # Extract OTUs
  features = as.matrix(physeq$features)
  dd<-colnames(features)

  # CSS Normalizing the Data
  # Create the metagenomeSeq object
  MGS = newMRexperiment(t(features), featureData=NULL, libSize=NULL, normFactors=NULL)
  # Trigger metagenomeSeq to calculate its Cumulative Sum scaling factor.
  MGS = cumNorm(MGS, p=cumNormStat(MGS))
  # Calculate scaling factors
  libSize <- normFactors(MGS)
  # Save the normalized data as data.frame
  features_norm = as.data.frame(t(MRcounts(MGS, norm=TRUE, log=FALSE)))

  # Rename the True Positive Features - Same Format as Before
  colnames(features_norm) <- dd;

  # Return as list
  return(list(metadata=metadata, features=features_norm, ID=ID, libSize=libSize))
}

# Apply CSS Normalization To A List of Datasets

list.CSS<-function(physeq){
  foreach(physeq = physeq, .packages = "metagenomeSeq", .export="CSSnorm", .errorhandling = 'remove') %dopar% {
    return(CSSnorm(physeq))
  }
}

######################
## TSS Normalization #
######################

# Apply TSS Normalization To A Dataset

TSSnorm = function(physeq) {

  # Extract ID and Metadata and Keep Them Unchanged
  ID<-physeq$ID
  metadata<-physeq$metadata

  # Extract OTUs
  features = as.matrix(physeq$features)
  dd<-colnames(features)

  # TSS Normalizing the Data
  features_norm <- decostand(features, method="total", MARGIN=1)

  # Convert back to data frame
  features_norm<-as.data.frame(features_norm)

  # Rename the True Positive Features - Same Format as Before
  colnames(features_norm) <- dd;

  # Force the Sequencing Depth to be the Original Library Size (Same as CLR)
  libSize<-rowSums(features)

  # Return as list
  return(list(metadata=metadata, features=features_norm, ID=ID, libSize=libSize))
}

# Apply TSS Normalization To A List of Datasets

list.TSS<-function(physeq){
  foreach(physeq = physeq, .packages = "vegan", .export="TSSnorm", .errorhandling = 'remove') %dopar% {
    return(TSSnorm(physeq))
  }
}

######################
## CLR Normalization #
######################

# Apply CLR Normalization To A Dataset

CLRnorm = function(physeq) {

  # Extract ID and Metadata and Keep Them Unchanged
  ID<-physeq$ID
  metadata<-physeq$metadata

  # Extract OTUs
  features = as.matrix(physeq$features)
  dd<-colnames(features)

  # Force the Sequencing Depth to be the Original Library Size (Same as TSS)
  libSize<-rowSums(features)

  # CLR
  features_norm_clr<-clr(features+1)

  # Convert back to data frame
  features_norm_clr<-as.data.frame(features_norm_clr)

  # Rename the True Positive Features - Same Format as Before
  colnames(features_norm_clr) <- dd;

  # Return as list
  return(list(metadata=metadata, features= features_norm_clr, ID=ID, libSize=libSize))
}

# Apply CLR Normalization To A List of Datasets

list.CLR<-function(physeq){
  foreach(physeq = physeq,
          .packages = c("vegan", "chemometrics"),
          .export="CLRnorm", .errorhandling = 'remove') %dopar% {
            return(CLRnorm(physeq))
          }
}

######################
## TMM Normalization #
######################

# Apply TMM Normalization To A Dataset

TMMnorm = function(physeq) {

  # Extract ID and Metadata and Keep Them Unchanged
  ID<-physeq$ID
  metadata<-physeq$metadata

  # Extract OTUs
  features = as.matrix(physeq$features)
  dd<-colnames(features)

  # TMM Normalizing the Data
  X<-t(features);
  libSize = edgeR::calcNormFactors(X,method="TMM") #Calculate normaization factors
  eff.lib.size = colSums(X)*libSize;
  ref.lib.size = mean(eff.lib.size); #Use the mean of the effective library sizes as a reference library size
  X.output = sweep(X,MARGIN=2,eff.lib.size,"/")*ref.lib.size; #Normalized read counts

  # Convert back to data frame
  features_norm<-as.data.frame(t(X.output))

  # Rename the True Positive Features - Same Format as Before
  colnames(features_norm) <- dd;

  # Return as list
  return(list(metadata=metadata, features=features_norm, ID=ID, libSize=libSize))
}

# Apply TMMnorm Normalization To A List of Datasets

list.TMM<-function(physeq){
  foreach(physeq = physeq, .packages = "edgeR", .export="TMMnorm", .errorhandling = 'remove') %dopar% {
    return(TMMnorm(physeq))
  }
}


######################
## RLE Normalization #
######################

# Apply RLE Normalization To A Dataset

RLEnorm = function(physeq) {

  # Extract ID and Metadata and Keep Them Unchanged
  ID<-physeq$ID
  metadata<-physeq$metadata

  # Extract OTUs
  features = as.matrix(physeq$features)
  dd<-colnames(features)

  # RLE Normalizing the Data
  formula <- as.formula(paste('~', paste(colnames(metadata), collapse = "+"), sep=''))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = t(features), colData = metadata, design = formula)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(dds), 1, gm_mean)
  dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  X.output<-counts(dds, normalized=TRUE)

  # Convert back to data frame
  features_norm<-as.data.frame(t(X.output))

  # Rename the True Positive Features - Same Format as Before
  colnames(features_norm) <- dd;

  # Extract Scaling Factors
  libSize<-sizeFactors(dds)

  # Return as list
  return(list(metadata=metadata, features=features_norm, ID=ID, libSize=libSize))
}

# Apply RLE Normalization To A List of Datasets

list.RLE<-function(physeq){
  foreach(physeq = physeq, .packages = "DESeq2", .export="RLEnorm", .errorhandling = 'remove') %dopar% {
    return(RLEnorm(physeq))
  }
}


# Convert A Limma/edgeR Output Matrix into Long Format
rename.features<-function(x, name){
  x<-reshape2::melt(x);
  colnames(x)<-c('feature', 'metadata', name);
  return(x);
}




# Calculated Evaluation Metris for A Single Data Frame

eval_res_list = function(resi, alpha=0.05) { #, simdataModel=simdataModel) {

  # #resi=methOutput[[1]][[2]]
  # if (simdataModel=="NULLmodel"){
  #   # Replace Missing Q-values to the Highest Possible Value 1.0
  #   resi[is.na(resi[, "pval"]), "pval"] <- 1
  #   # Evaluate Detection Performance
  #   time= mean(resi[,"time"], na.rm=TRUE)
  #   wh.pred = (resi[, "pval"] < alpha)
  # }
  #
  # else{
      # Replace Missing Q-values to the Highest Possible Value 1.0
    resi[is.na(resi[, "qval"]), "qval"] <- 1
      # Evaluate Detection Performance
    time= mean(resi[,"time"], na.rm=TRUE)
    wh.pred = (resi[, "qval"] < alpha)
  #}
  wh.pos = which(wh.pred)
  wh.neg = which(!wh.pred)
  wh.TP = grep("[[:print:]]+\\_TP$", resi[, "pairwiseAssociation"])
  FPs = sum(!wh.pos %in% wh.TP)
  TPs = sum(wh.pos %in% wh.TP)
  TNs = sum(!wh.neg %in% wh.TP)
  FNs = sum(wh.neg %in% wh.TP)

  # Sensitivity: True Positives Divided by All Positives (Sum of True
  # Positives and False Negatives)
  Sensitivity = TPs/(TPs + FNs)

  # Specificity: True Negatives Divided by All Negatives (Sum of True
  # Negatives and False Positives)
  Specificity = TNs/(TNs + FPs)

  # False Discovery Rate: False Positives Divided by All Detected Positives
  FDR = if ((TPs + FPs) == 0)
    0 else FPs/(TPs + FPs)
  # If no true positives, return NA's for irrelevant measures

  # FScore
  FScore<-2*TPs/(2*TPs + FNs + FPs)

  # Matthew's Correlation Coefficient
  numerator <- (TPs * TNs - FPs * FNs)
  denominator <- sqrt((TPs + FPs)*(TPs + FNs)) * sqrt((TNs + FPs)*(TNs + FNs))
  if(denominator == 0 ) denominator <- 1
  MCC <- numerator/denominator

  # AUC, pAUC (FPR < 0.20)
  wh.truth = (1:nrow(resi) %in% wh.TP)
  if (all(wh.truth=='TRUE')) {
    AUC=1
    pAUC=1
    fAUC=1
  } else if (all(wh.truth=='FALSE')) {
    AUC=0
    pAUC=0
  } else {
    pred <- prediction(as.numeric(wh.pred), factor(wh.truth, levels=c("TRUE", "FALSE")))
    AUC = performance(pred, "auc")@y.values[[1]]
    pAUC = performance(pred, "auc", fpr.stop=0.20)@y.values[[1]][[1]]
  }

  # df.res = data.frame(FPs,TPs,TNs,FNs,length(wh.TP),dim(resi)[1],Sensitivity,Specificity,FDR,MCC,AUC,"DESeq2")
  # colnames(df.res)=c("FP","TP","TN","FN","Spiked","nFeatures","Sensitivity","Specificity","FDR","MCC","AUC","Method")
  #

  # Departure from Uniformity Under the Null
  AucAocVals <- AucAocFun(resi[!wh.truth, "pval"], plotIt = FALSE, pch = 20, type = "l")
  totalArea = AucAocVals["conservArea"] + AucAocVals["liberalArea"]
  names(totalArea) = "totalArea"

  # Return
  return(c(Sensitivity = Sensitivity,
           Specificity=Specificity,
           FDR=FDR,
           FScore=FScore,
           MCC=MCC,
           AUC=AUC,
           pAUC=pAUC,
           AucAocVals["conservArea"],
           AucAocVals["liberalArea"],
           totalArea,
           time=time))
}


# Extract Parameters from A List of Results

make_power_df = function(reslist, simparamslabels) {
  powerdf = ldply(reslist, .parallel=TRUE)
  colnames(powerdf)[1] <- "Combinations"
  paramdf = ldply(strsplit(powerdf[, "Combinations"], "_"), .parallel=TRUE)
  colnames(paramdf) <- simparamslabels
  powerdf = cbind(powerdf, paramdf)
  return(powerdf)
}


# Organize All Lists into A Coherent Data Frame

list.perfdf<-function(resultslist, simparams, simparamslabels){
  foreach(resultslist = resultslist, .export=c("make_power_df", "eval_res_list", "AucAocFun"),
          .packages = c("ROCR", "plyr"), .errorhandling = 'remove') %dopar% {
  perflist <- lapply(resultslist, eval_res_list)
  if (is.null(names(resultslist))) {
    names(perflist) <- simparams
  } else {
    names(perflist) <- names(resultslist)
  }
  perfdf = make_power_df(perflist, simparamslabels)
  return(perfdf)
}}



#######################################
# Arc-Sine Square Root Transformation #
#######################################

# Arc Sine Square Root Transformation
AST<-function(x){
  return(sign(x)*asin(sqrt(abs(x))))
}

########################
# Logit Transformation #
########################

# LOGIT Transformation
LOGIT<-function(x){
  y<-car::logit(x, adjust=0)
  y[!is.finite(y)]<-0
  return(y)
}

# Shifted Logit Transformation (Lukens et al, 2014, Nature)
# LOGIT_S<-function(x){
#   y<-0.5*log(x/(1-x)) + 10
#   y[!is.finite(y)]<-0
#   return(y)
# }

######################
# Log Transformation #
######################

# LOG Transformation
LOG<-function(x){
  return(log(x+1))
}

# Apply SimpleFilter (Prevalence and Abundance Filtering) To A Dataset

# SimpleFilter = function(physeq, Threshold_Abundance, Threshold_Prevalence) {
#   features<-physeq$features
#   filtered_features<-features[,colSums(features > Threshold_Abundance) > nrow(features)*Threshold_Prevalence]
#   physeq$features<-filtered_features
#   return(physeq)
# }

SimpleFilter = function(physeq, Threshold_Reads, Threshold_Samples) {
  features<-physeq$features
  #filtered_features<-features[,colSums(features > Threshold_Abundance) > nrow(features)*Threshold_Prevalence]
  #physeq$features<-filtered_features
  metadata=physeq$metadata
  group=factor(metadata[,1])
  ugroup=unique(group)
  filt_rep=list()
  df_counts = cbind(features,group)
  for(i in 1:length(ugroup)){
    sub_filt=subset(df_counts, group == ugroup[i])
    dd=dim(sub_filt)[2]
    sub_filt = sub_filt[,-dd]
    kp = colSums(sub_filt >= Threshold_Reads) >= Threshold_Samples
    sub_filt=sub_filt[,kp]
    colnames_filt=colnames(sub_filt)
    filt_rep[[i]]=colnames_filt
  }
  filtered_cols <- which(colnames(features) %in% unique(unlist(filt_rep)))
  physeq$features <- features[,filtered_cols]
  return(physeq)
}


# Apply SimpleFilter (Prevalence and Abundance Filtering) To A List of Datasets

# list.SimpleFilter<- function(simlist, Threshold_Abundance, Threshold_Prevalence){
#   return(lapply(simlist, SimpleFilter, Threshold_Abundance, Threshold_Prevalence))
# }

list.SimpleFilter<- function(simlist, Threshold_Reads, Threshold_Samples){
  return(lapply(simlist, SimpleFilter, Threshold_Reads, Threshold_Samples))
}

# Calculate No. of TP Features in A Dataset

TruePositives<-function(physeq, breakitdown) {
  features<-physeq$features
  d<-colnames(features)[grep("[[:print:]]+\\_TP$", colnames(features))]
  if(breakitdown==TRUE) {
    return(d)
  } else {
    return(length(d))
  }
}

# Calculate No. of TP Features in A List of Datasets

list.TruePositives<- function(simlist, breakitdown){
  return(lapply(simlist, TruePositives, breakitdown))
}

# Calculate Prevalence of All Features in A Dataset

calculatePrevalence<-function(physeq) {
  features<-physeq$features
  d<-colSums(features > Threshold_Abundance)/nrow(features)
  return(d)
}

# Calculate Prevalence of All Features in A List of Datasets

list.calculatePrevalence<- function(simlist){
  return(lapply(simlist, calculatePrevalence))
}

# Check if a metadata is binary

is.binary <- function(v) {
  x <- unique(v)
  dim(x)[1] - sum(is.na(x)) == 2L
}


# Function from Hawinkel et al. (2017)
## A function to look at the distribution of the P-values under the null
## distribution

### Compute AUC and AOC of P-value distribution to be applied on the Monte
### Carlo distribution of each OTU separately
AucAocFun <- function(pVals, maxVal = 0.25, plotIt = FALSE, ...) {
  ## sum of height differences between the curve and the y=x line, half maximum
  ## = 0.5 * max(pVals) * length(pVals) this is the case of 100% rejection with
  ## any alpha value _onlyTotArea_ controls if only the total area is returned,
  ## and not also the Area Over and Under the y=x line.  halfMaxArea <- 0.5 *
  ## max(maxVal) * length(pVals)
  halfMaxArea <- 0.5 * length(pVals)

  pVals <- pVals[!is.na(pVals)]
  estimated <- sort(pVals)
  theoretic <- seq_along(pVals)/length(pVals)

  if (plotIt) {
    plot(theoretic, estimated, ...)
    abline(0, 1)
  } else {
  }

  diffPVals <- theoretic - estimated
  indConserv <- theoretic <= estimated
  conservArea <- sum(-diffPVals[indConserv])/halfMaxArea
  liberalArea <- sum(diffPVals[!indConserv])/halfMaxArea

  c(conservArea = conservArea, liberalArea = liberalArea, totalArea = liberalArea +
      conservArea)

}  # END - function: aucAocFun, AUC/AOC calculation for p-values distribution

##Ensemble
## relevance matrix formation

resRelevanceMatrix=function(itern,resMerge){
  resMerge.df=data.frame(resMerge[[itern]][[1]],resMerge[[itern]][[2]], resMerge[[itern]][[3]],
                         resMerge[[itern]][[6]],resMerge[[itern]][[7]], resMerge[[itern]][[14]],resMerge[[itern]][[15]],
                         resMerge[[itern]][[22]],resMerge[[itern]][[23]])
  colnames(resMerge.df)=c("pairwiseAssociation","feature","metadata","qvalLimmaVoomZinb","relevanceLimmaVoomZinb",
                          "qvalDESeq2Zinb","relevanceDESeq2Zinb","qvalEdgeRZinb","relevanceEdgeRZinb")
  resMerge.df[is.na(resMerge.df)] <- 0
  resMerge.df
}


#methEnsemble=methEnseOutput[[2]]
#str(methEnsemble)
#ensm_res_list(methEnsemble)

ensm_res_list = function(methEnsemble) {

  methEnsembleMatrix=methEnsemble[,-c(1,2,3)]
  methEnsemble$RelevanceScore= rowMeans(methEnsembleMatrix[,c(2,4,6)])
  #methEnsemble$pvalEnsemble=(methEnsemble$qvalLimmaVoomZinb^0.2 * methEnsemble$qvalDESeq2Zinb^0.3 * methEnsemble$qvalEdgeRZinb^0.5)
  #methEnsemble$pvalEnsemble=(methEnsemble$qvalLimmaVoomZinb^0.002 * methEnsemble$qvalDESeq2Zinb^0.002 * methEnsemble$qvalEdgeRZinb^0.996)
  #methEnsemble$pvalEnsemble=-2*( log(methEnsemble$qvalLimmaVoomZinb) + log(methEnsemble$qvalDESeq2Zinb) + log(methEnsemble$qvalEdgeRZinb))

  ensemble_cutoff=1  #0.05 #quantile(methEnsemble$RelevanceScore, 0.65, 0.75) #mean(methEnsemble$RelevanceScore) #
  #ensemble_cutoff=0.05
  #ensemble_cutoff=qchisq(0.95,6)
  # methEnsemble$important= (methEnsemble$RelevanceScore > ensemble_cutoff)
  # methEnsemble <- methEnsemble[order(methEnsemble$important, decreasing=TRUE),]
  # wh_predicted=subset(methEnsemble,important=="TRUE")
  #
  #time= mean(methEnsemble[,"time"], na.rm=TRUE)
  wh.pred = (methEnsemble[, "RelevanceScore"] >= ensemble_cutoff)
  #wh.pred = (methEnsemble[, "pvalEnsemble"] < ensemble_cutoff)
  #wh.pred = (methEnsemble[, "pvalEnsemble"] > ensemble_cutoff)
  wh.pos = which(wh.pred)
  wh.neg = which(!wh.pred)
  wh.TP = grep("[[:print:]]+\\_TP$", methEnsemble[, "pairwiseAssociation"])
  #methEnsemble$truth= (methEnsemble$pairwiseAssociation %in% methEnsemble[wh.TP,"pairwiseAssociation"])*1
  # sum(methEnsemble$relevanceLimmaVoomZinb)
  #glmfit=glm(truth~relevanceLimmaVoomZinb + relevanceDESeq2Zinb + relevanceEdgeRZinb, data=methEnsemble)
  #summary(glmfit)
  #glmfit$coefficients
  # glmfit=glm(truth~ relevanceEdgeRZinb, data=methEnsemble)
  # glmfit=glm(truth~ qvalLimmaVoomZinb+qvalDESeq2Zinb + qvalEdgeRZinb, data=methEnsemble)
  # summary(glmfit)
  FPs = sum(!wh.pos %in% wh.TP)
  TPs = sum(wh.pos %in% wh.TP)
  TNs = sum(!wh.neg %in% wh.TP)
  FNs = sum(wh.neg %in% wh.TP)
  # plot(sort(methEnsemble$qvalEdgeRZinb))
  # points(sort(methEnsemble$qvalDESeq2Zinb))
  # points(sort(methEnsemble$qvalLimmaVoomZinb))
  #0.03*methEnsemble$qvalDESeq2Zinb+
  #methEnsemble$RelevanceScore=(0.70*methEnsemble$qvalLimmaVoomZinb + 0.05*methEnsemble$qvalEdgeRZinb)/(0.7+0.05)
  # Sensitivity: True Positives Divided by All Positives (Sum of True
  # Positives and False Negatives)
  Sensitivity = TPs/(TPs + FNs)

  # Specificity: True Negatives Divided by All Negatives (Sum of True
  # Negatives and False Positives)
  Specificity = TNs/(TNs + FPs)

  # False Discovery Rate: False Positives Divided by All Detected Positives
  FDR = if ((TPs + FPs) == 0)
    0 else FPs/(TPs + FPs)
  # If no true positives, return NA's for irrelevant measures

  # FScore
  FScore<-2*TPs/(2*TPs + FNs + FPs)

  # Matthew's Correlation Coefficient
  numerator <- (TPs * TNs - FPs * FNs)
  denominator <- sqrt((TPs + FPs)*(TPs + FNs)) * sqrt((TNs + FPs)*(TNs + FNs))
  if(denominator == 0 ) denominator <- 1
  MCC <- numerator/denominator

  # AUC, pAUC (FPR < 0.20)
  wh.truth = (1:nrow(methEnsemble) %in% wh.TP)
  if (all(wh.truth=='TRUE')) {
    AUC=1
    pAUC=1
    fAUC=1
  } else if (all(wh.truth=='FALSE')) {
    AUC=0
    pAUC=0
  } else {
    pred <- prediction(as.numeric(wh.pred), factor(wh.truth, levels=c("TRUE", "FALSE")))
    AUC = performance(pred, "auc")@y.values[[1]]
    pAUC = performance(pred, "auc", fpr.stop=0.20)@y.values[[1]][[1]]
  }

  # Departure from Uniformity Under the Null
  #AucAocVals <- AucAocFun(methEnsemble[!wh.truth, "pval"], plotIt = FALSE, pch = 20, type = "l")
  #totalArea = AucAocVals["conservArea"] + AucAocVals["liberalArea"]
  #names(totalArea) = "totalArea"

  # Return
  return(c(Sensitivity = Sensitivity,
           Specificity=Specificity,
           FDR=FDR,
           FScore=FScore,
           MCC=MCC,
           AUC=AUC,
           pAUC=pAUC))
}


# features_importance_list = function(ensm_res_list) {
#
#   # identify significant features
#   importance_matrix = cbind(ensm_res_list$egdeRZinbwave, ensm_res_list$DEDeq2Zinbwave, ensm_res_list$limmavoomZinbwave)
#   feature_importance = rowMeans(importance_matrix)
#   mean_cut_off = mean(feature_importance)
#
#
#   if(feature_importance > mean_cut_off)
#   # Evaluate Detection Performance
#
#   time= mean(resi[,"time"], na.rm=TRUE)
#   wh.pred = (resi[, "qval"] < alpha)
#   wh.pos = which(wh.pred)
#   wh.neg = which(!wh.pred)
#
#
#   wh.TP = grep("[[:print:]]+\\_TP$", resi[, "pairwiseAssociation"])
#
#
#   # Return
#   return(c( wh.pred =  wh.pred,
#             wh.pos=wh.pos,
#             time=time))
# }



