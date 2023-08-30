#Code is modified from https://github.com/biobakery/maaslin2_benchmark
run_simulator_SparseDOSSA2 <- function(es, ns, nIterations, workingDirectory,
                                       nMicrobes,readDepth,filtered_real_counts)

{
  SparseDOSSA2_fit <- SparseDOSSA2::fit_SparseDOSSA2(data = filtered_real_counts ,
                                                     control = list(verbose = TRUE))

  for(b in 1:length(es)){
    for(k in 1:length(ns)){

      # Set Parameters for Testing Purposes
      noZeroInflate<- FALSE # High-level parameter
      RandomEffect<-FALSE # High-level parameter
      metadataType <- 'UVB' #'MVA' # High-level parameter
      nSubjects <- ns[k] # # Low-level parameter
      nPerSubject<-1 # Low-level parameter
      nMicrobes <- nMicrobes #3413 #945 #758 # 239, 168 #nMicrobes # Low-level parameter
      spikeMicrobes <- 0.1 # Low-level parameter
      nMetadata<- 1 #5 # Low-level parameter
      spikeMetadata<- 1 #0.2 # Low-level parameter
      effectSize<- es[b] # Low-level parameter
      readDepth<- readDepth #248434 #247294 #246709 #1708 #1438 #readDepth # Default parameter
      noParallel<- FALSE # Default parameter
      nIterations<- nIterations # Default parameter
      rSeed<- 1234 # Default parameter
      nCores<- 4 # Default parameter
      # workingDirectory <- '/Users/hmallick/Dropbox (Huttenhower Lab)/Maaslin2' # Default parameter


      if (noZeroInflate==TRUE && RandomEffect==TRUE){
        inputSubString<-'/noZeroInflate_RandomEffect'
      }

      if (noZeroInflate==TRUE && RandomEffect==FALSE){
        inputSubString<-'/noZeroInflate_noRandomEffect'
      }

      if (noZeroInflate==FALSE && RandomEffect==TRUE){
        inputSubString<-'/ZeroInflate_RandomEffect'
      }

      if (noZeroInflate==FALSE && RandomEffect==FALSE){
        inputSubString<-'/ZeroInflate_noRandomEffect'
      }

      inputString<-paste(inputSubString, metadataType, nSubjects, nPerSubject, nMicrobes, spikeMicrobes, nMetadata, spikeMetadata, effectSize, nIterations, readDepth, sep='_')


      # Create Input Directory
      inputDirectory <- file.path(workingDirectory, 'Input')

      if (!dir.exists(inputDirectory)){
        print("Input directory created!")
        dir.create(inputDirectory)
      } else {
        print("Input directory already exists!")
      }

      inputStringR<-paste(inputDirectory, paste(inputString, '.RData', sep=''), sep='')


      # If input does not exist, then only generate data. Otherwise, skip.

      if (!file.exists(inputStringR)){
        # Do some operation based on user input
        simlist<-newtrigger_sparseDOSSA_Simulator(noZeroInflate=noZeroInflate,
                                                  RandomEffect=RandomEffect,
                                                  metadataType=metadataType,
                                                  nSubjects=nSubjects,
                                                  nPerSubject=nPerSubject,
                                                  nMicrobes=nMicrobes,
                                                  spikeMicrobes = spikeMicrobes,
                                                  nMetadata = nMetadata,
                                                  spikeMetadata=spikeMetadata,
                                                  effectSize = effectSize,
                                                  #df_spike_metadata=df_spike_metadata,
                                                  features_otu=features_otu,
                                                  readDepth = readDepth,
                                                  SparseDOSSA2_fit=SparseDOSSA2_fit,
                                                  #physeqtr=physeqtr,
                                                  noParallel = noParallel,
                                                  nIterations=nIterations,
                                                  nspike=nspike,
                                                  rSeed=rSeed,
                                                  nCores=nCores)


        #simlist

        save(simlist, file=inputStringR)
      } else{
        print("Input file already exists. No new data generated!")
        load(inputStringR)
      }


      #######################################
      # Delete Temporary sparseDOSSA Files  #
      #######################################

      if (file.exists("SyntheticMicrobiome-Counts.pcl")) file.remove("SyntheticMicrobiome-Counts.pcl")
      if (file.exists("SyntheticMicrobiome.pcl")) file.remove("SyntheticMicrobiome.pcl")
      if (file.exists("SyntheticMicrobiomeParameterFile.txt")) file.remove("SyntheticMicrobiomeParameterFile.txt")

    }
  }

}
