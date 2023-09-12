clean_data_sim <- function(physeqLOAD)
{

  physeqtr <- filter_taxa(physeqLOAD, function(x) sum(x > 0) >= (0.01*length(x)), TRUE)
  filtered_real_counts = (otu_table(physeqtr))
  dim(filtered_real_counts)
  readDepth=round(median(colSums(filtered_real_counts)),0)
  nMicrobes=nrow(filtered_real_counts)


  featuress = t(otu_table(physeqLOAD, 'matrix'))
  metadatas = sample_data(physeqLOAD)
  group=factor(metadatas$group)
  ugroup=unique(group)
  df_counts = data.frame(featuress,group)

  frepi=list()
  for(i in 1:length(ugroup)){
    bhasf=subset(df_counts, group == ugroup[i])
    dim(bhasf)
    bhasfs=bhasf[,-dim(df_counts)[2]] #-c(11794,11795)]
    kp = colSums(bhasfs >= 2) >= round(0.10*length(group),0)  #20%
    d2k=bhasfs[,kp]
    dim(d2k)
    frep=colnames(d2k)
    frepi[[i]]=frep
  }
  fcolst <- unique(unlist(frepi))
  length(fcolst)
  features_otu1 <- featuress[,fcolst]
  dim(features_otu1)

  #remove taxa with group-wise structured zeros from the spike-ins

  df_counts1 = data.frame(features_otu1,group)

  filti=list()
  for(i in 1:length(ugroup)){
    cbhasf=subset(df_counts1, group == ugroup[i])
    dflast=dim(cbhasf)[2]
    cbhasff=cbhasf[,-dflast]
    filt=colSums(cbhasff) == 0
    filt_d2k=cbhasff[,filt]
    dim(filt_d2k)
    filtr=colnames(filt_d2k)
    filti[[i]]=filtr
  }
  filt_fcols <- unique(unlist(filti))
  length(filt_fcols)

  #ex_counts1 <- features[,filti[[1]]]
  drop1 <- filt_fcols #[[1]]
  if(length(drop1)==0) {
    ex_features <- NULL
    }
  else {
    ex_features <- features_otu1[, drop1]
    }
  #features retained for further analysis
  features_otu = features_otu1[,!(colnames(features_otu1) %in% drop1)]
  dim(features_otu)
  #number of spike-in features

  #features_otu <- features_otu1

  dim(features_otu)

  return(list(filtered_real_counts=filtered_real_counts,features_otu=features_otu,zeroGroup_features=ex_features,
              readDepth=readDepth,nMicrobes=nMicrobes))

}
