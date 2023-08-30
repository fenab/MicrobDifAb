generateMetadata<-function(metadataType,
                           nSubjects,
                           nPerSubject,
                           nMetadata,
                           spikeMetadata){
  
  # Calculate Number of Samples
  nSamples = round(nSubjects*nPerSubject)
  
  # Create Blocking Variable
  if (nPerSubject==1){  # NO RANDOM EFFECTS
    subjectRandomEffects <- as.matrix(rnorm(nSubjects,mean=0,sd=0))
  }
  if (nPerSubject>1){  # SUBJECT-SPECIFIC RANDOM EFFECTS
    subjectRandomEffects <- as.matrix(rnorm(nSubjects,mean=0,sd=1))
  }
  BLOCK <- as.vector(matrix(subjectRandomEffects,nrow=nPerSubject,ncol=length(subjectRandomEffects),byrow=TRUE))
  
  # Specify Mean and Covariance Structure
  mu<-rep(0,nMetadata)
  cov<-diag(1,nMetadata, nMetadata)
  
  if (metadataType == 'MVB'){
    for (i in 1:nMetadata){
      for (j in 1:nMetadata){
        if(i!=j) cov[i,j]=0.5**(abs(i-j)) # AR(1)
      }
    }
  }
  
  # Generate from MVN
  fakeMetadata<-as.matrix(MASS::mvrnorm(n=nSamples, mu,cov))
  
  # Transpose and Add Blocking Structure
  finalMetadata<-apply(fakeMetadata, 2, function(x) x+BLOCK)
  
  #############################
  # Modularize Specific Cases #
  #############################
  
  # Multivariable Scenario - Dichotomize Half of the Features
  if (metadataType %in% c('MVA', 'MVB')){
    t_UserMetadata<-apply(finalMetadata, 2, function(x) ifelse(x>median(x), 1, 0))
    columns_not_to_binarize<-sample(1:nMetadata, nMetadata/2)
    t_UserMetadata[,columns_not_to_binarize]<-finalMetadata[, columns_not_to_binarize]
    UserMetadata<-t(t_UserMetadata)
  }
  
  # Univariate Binary
  else if (metadataType == 'UVB'){
    UserMetadata<-t(apply(finalMetadata, 2, function(x) ifelse(x>median(x), 1, 0)))
  }
  # Univariate Continuous
  else {
    UserMetadata<-t(finalMetadata)
  }
  
  # Collect Relevant Spike-in Information
  spikeCount<- round(nMetadata*spikeMetadata)
  Metadatafrozenidx<-sample(1:nMetadata, spikeCount, replace=FALSE)
  significant_metadata<-paste('Metadata', Metadatafrozenidx, sep='')
  spikeCount<-as.character(spikeCount)
  
  # Return
  return(list(UserMetadata=UserMetadata,
              Metadatafrozenidx=Metadatafrozenidx,
              significant_metadata=significant_metadata,
              spikeCount=spikeCount))
}
