#Code is modified from https://github.com/biobakery/maaslin2_benchmark
summarizing_performance <- function(nMicrobes,readDepth, nIterations, methodName, ns, es)
  {

    ss=length(methodName)
    noZeroInflate<- FALSE # High-level parameter
    RandomEffect<-FALSE # High-level parameter
    metadataType <- 'UVB' #'MVA' # High-level parameter
    nSubjects <- ns #c(10,20,50,100,200) # # Low-level parameter
    nPerSubject<-1 # Low-level parameter
    nMicrobes <- nMicrobes #3413 #945 #758 # 239, 168 #nMicrobes # Low-level parameter
    spikeMicrobes <- 0.1 # Low-level parameter
    nMetadata<- 1 #5 # Low-level parameter
    spikeMetadata<- 1 #0.2 # Low-level parameter
    effectSize<- es #c(0.5,1,2) # Low-level parameter
    readDepth<- readDepth #248434 #247294 #246709 #1708 #1438 #readDepth # Default parameter
    noParallel<- FALSE # Default parameter
    nIterations<- nIterations # Default parameter
    rSeed<- 1234 # Default parameter
    nCores<- 4 # Default parameter

    rResults=NULL

    for(m in 1:length(nSubjects)){
      nResults=NULL

      for(t in 1:length(effectSize)){
        kResults=NULL
        outputDirectory <- file.path(workingDirectory, 'Output')
        for(s in 1:ss){


          outputSubString1 <-'ZeroInflate_noRandomEffect'

          #outputString1 <-paste(methodName1[s], outputSubString1, metadataType[1], nSubjects[m], nPerSubject[1], nMicrobes[1], spikeMicrobes[1], nMetadata[1], spikeMetadata[1], effectSize[t], readDepth[1], sep='_')
          outputString1 <-paste(methodName[s], outputSubString1, metadataType[1], nSubjects[m], nPerSubject[1], nMicrobes[1], spikeMicrobes[1], nMetadata[1], spikeMetadata[1], effectSize[t], readDepth[1], sep='_')
          outputStringR1<-paste(outputDirectory, paste(outputString1, '.RData', sep=''), sep='/')
          # }

          fao=get(load(outputStringR1))

          kResults=rbind(kResults, fao)
        }

        nResults=rbind(nResults, kResults)

      }
      rResults=rbind(rResults, nResults)
    }
    dim(rResults)
    rResultsAll=rResults #rbind(rResults,rResultsEns,rResultsEnsRob)
    return(rResultsAll)
}

