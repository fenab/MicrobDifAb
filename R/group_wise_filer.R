group_wise_filer  <- function(group, counts, min_counts, min_replicates){
  ugroup <- sort(unique(group))
  if(length(ugroup)==2){
    df_counts <- data.frame(counts, group)

    filtAll <- lapply(1:length(ugroup), function(i) {
      bhasf <- subset(df_counts, group == ugroup[i])
      bhasfs <- bhasf[,-dim(df_counts)[2]]
      kp <- colSums(bhasfs >= min_counts) >= min_replicates
      d2k <- bhasfs[, kp]
      frep <- colnames(d2k)
      return(frep)
    })

    fcolst <- unique(unlist(filtAll))
    length(fcolst)
    features <- counts[, fcolst]

    df_counts1 = data.frame(features,group)
    filtZero <- lapply(1:length(ugroup), function(i) {
      bhasf <- subset(df_counts1, group == ugroup[i])
      bhasfs <- bhasf[,-dim(df_counts1)[2]]
      kp <- colSums(bhasfs > 0) == 0
      d2k <- bhasfs[, kp]
      frep <- colnames(d2k)
      return(frep)
    })

    filtZero_fcols <- unique(unlist(filtZero))
    length(filtZero_fcols)
    if(length(filtZero_fcols)==0) ex_counts <- NULL
    else  ex_counts <- features[,filtZero_fcols]

    ffeatures1 = features[,!(colnames(features) %in% filtZero_fcols)]

    return=list(taxa_All_filt=features, taxa_NonZeroGroup=ffeatures1, taxa_ZeroGroup=ex_counts)
  }

  else{
    #ugroup <- unique(group, incomparables = TRUE)
    df_counts <- data.frame(counts, group)

    filtAll <- lapply(1:length(ugroup), function(i) {
      bhasf <- subset(df_counts, group == ugroup[i])
      bhasfs <- bhasf[,-dim(df_counts)[2]]
      kp <- colSums(bhasfs >= min_counts) >= min_replicates
      d2k <- bhasfs[, kp]
      frep <- colnames(d2k)
      return(frep)
    })

    fcolst <- unique(unlist(filtAll))
    length(fcolst)

    features <- counts[, fcolst]

    df_counts1 = data.frame(features,group)
    filtZero <- lapply(1:length(ugroup), function(i) {
      bhasf <- subset(df_counts1, group == ugroup[i])
      bhasfs <- bhasf[,-dim(df_counts1)[2]]
      kp <- colSums(bhasfs > 0) == 0
      d2k <- bhasfs[, kp]
      frep <- colnames(d2k)
      return(frep)
    })

    filtZero_fcols12 <- unique(c(filtZero[[1]],filtZero[[2]]))
    length(filtZero_fcols12)
    ex_counts12 <- features[,filtZero_fcols12]
    ffeatures12 = features[,!(colnames(features) %in% filtZero_fcols12)]

    filtZero_fcols13 <- unique(c(filtZero[[1]],filtZero[[3]]))
    length(filtZero_fcols13)
    ex_counts13<- features[,filtZero_fcols13]
    ffeatures13 = features[,!(colnames(features) %in% filtZero_fcols13)]


    return=list(taxa_All_filt=features, taxa_NonZeroGroup12=ffeatures12, taxa_ZeroGroup12=ex_counts12,
                taxa_NonZeroGroup13=ffeatures13, taxa_ZeroGroup13=ex_counts13)
  }
}
