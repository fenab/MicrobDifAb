
##Code is modified from https://github.com/biobakery/maaslin2_benchmark

run_performance_simulations <- function(es, ns, nIterations,
                                        nMicrobes,readDepth) #, simdataModel=c("spikeinModel", "NULLmodel"))
{
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


      max_percent_outliers = 0
      percent_outlier_spikins = 0



      # Set Up Clustering Environment
      no_cores <- nCores  # detectCores() - 1
      cl <- makeCluster(no_cores)
      registerDoParallel(cl)


      # Track Start Time
      cat(c("Job started at:",date()), "\n")
      start.time <- Sys.time()


      #a few selected methods
      model_iter<-c('DESeq2','edgeR',
                    'limmaVOOM',
                    'limmaVOOMZinbwave',
                    'DESeq2Zinbwave',
                    'edgeRZinbwave',
                    'MaAsLin2')

     #Normalization and transformation methods
      normParams<-c('CSS', 'CLR', 'RLE', 'TMM', 'TSS')
      transfParams<-c('AST', 'LOG', 'LOGIT')

      modelParams <- model_iter

      methOutput=list()

      for (j in 1:length(model_iter)){

        #readDepth<- 50010

        # Set Parameters for Testing Purposes
        methodName =  model_iter[j]  #'edgeRZinbwave' 'limmaVOOMZinbwave' #LM.CLR' # High-level parameter

        ###################################################
        # Extract Method + Normalization + Transformation #
        ###################################################

        inputMethodParams<-unlist(strsplit(methodName, "[.]"))
        modelName<-inputMethodParams[inputMethodParams %in% modelParams]
        normMethod<- inputMethodParams[inputMethodParams %in% normParams]
        transfMethod<- inputMethodParams[inputMethodParams %in% transfParams]

        if(is.character(modelName) & length(modelName) == 0) stop('The input model is invalid!')
        if(is.character(transfMethod) & length(transfMethod) == 0) transfMethod='NONE'
        if(is.character(normMethod) & length(normMethod) == 0) normMethod='NONE'


          # Track Start Time
        cat(c("Job started at:",date()), "\n")
        start.time1 <- Sys.time()

        # Extract Pre-generated Dataset
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

        # Load Dataset
        if (!file.exists(inputStringR)) stop('The input file does not exist. Generate the dataset first.')
        load(inputStringR)

        # Set Up Clustering Environment
        no_cores <- nCores  # detectCores() - 1
        cl <- makeCluster(no_cores)
        registerDoParallel(cl)

        ####################
        # Simple Filtering #
        ####################

        Threshold_Abundance = 0
        Threshold_Prevalence = 0.01
        simlist.filtered <- list.SimpleFilter(simlist, Threshold_Abundance, Threshold_Prevalence)


        # Threshold_Reads = 2
        # Threshold_Samples =2
        # simlist.filtered <- list.SimpleFilter(simlist, Threshold_Reads, Threshold_Samples)

        # Choose Appropriate Normalized Dataset
        if (normMethod=='NONE') {
          fitlist<-simlist.filtered
        }
        #str(fitlist[[1]])
        #####################
        # RFA Normalization #
        #####################

        if (normMethod=='RFA') {
          RFAlist<-list.RFA(simlist.filtered); names(RFAlist) <- names(simlist.filtered)
          fitlist<-RFAlist
        }

        #####################
        # TSS Normalization #
        #####################

        if (normMethod=='TSS') {
          TSSlist<-list.TSS(simlist.filtered); names(TSSlist) <- names(simlist.filtered)
          fitlist<-TSSlist
        }

        #####################
        # CLR Normalization #
        #####################
        if (normMethod=='CLR') {
          CLRlist<-list.CLR(simlist.filtered); names(CLRlist) <- names(simlist.filtered)
          fitlist<-CLRlist
        }

        #####################
        # CSS Normalization #
        #####################

        if (normMethod=='CSS') {
          CSSlist<-list.CSS(simlist.filtered); names(CSSlist) <-  names(simlist.filtered)
          fitlist<-CSSlist
        }

        #####################
        # TMM Normalization #
        #####################

        if (normMethod=='TMM') {
          TMMlist<-list.TMM(simlist.filtered); names(TMMlist) <-  names(simlist.filtered)
          fitlist<-TMMlist
        }

        #####################
        # RLE Normalization #
        #####################

        if (normMethod=='RLE') {
          RLElist<-list.RLE(simlist.filtered); names(simlist.filtered) <-  names(simlist.filtered)
          fitlist<-RLElist
        }


        # How Many Datasets to Run In A List
        reps<-1:nIterations
        fitlist<-fitlist[reps]

        # Assign Consistent Names to Results
        if (noZeroInflate==TRUE && RandomEffect==TRUE){
          outputSubString<-'noZeroInflate_RandomEffect'
        }

        if (noZeroInflate==TRUE && RandomEffect==FALSE){
          outputSubString<-'noZeroInflate_noRandomEffect'
        }

        if (noZeroInflate==FALSE && RandomEffect==TRUE){
          outputSubString<-'ZeroInflate_RandomEffect'
        }

        if (noZeroInflate==FALSE && RandomEffect==FALSE){
          outputSubString<-'ZeroInflate_noRandomEffect'
        }
        #readDepth=50005
        outputString<-paste(methodName, outputSubString, metadataType, nSubjects, nPerSubject, nMicrobes, spikeMicrobes, nMetadata, spikeMetadata, effectSize, readDepth, sep='_')

        # Create Output Directory
        outputDirectory <- file.path(workingDirectory, 'Output')

        if (!dir.exists(outputDirectory)){
          print("Output directory created!")
          dir.create(outputDirectory)
        } else {
          print("Output directory already exists!")
        }

        # Save the Results in A .RData File
        outputStringR<-paste(outputDirectory, paste(outputString, '.RData', sep=''), sep='/')

        # Run the Model Only if the Output Does Not Exist
        if (!file.exists(outputStringR)){
          print(methodName)
          # Model
          list.function<- eval(parse(text=paste('list', modelName, sep ='.')))
          output<-eval(parse(text="list.function(fitlist, transformation=transfMethod)"))

          # Discard Failed Methods
          if (length(output)<1) stop('Consistent error in the model fitting. No output returned.')
          names(output)<-paste(outputString, 1:length(output), sep='_')


          # Organize into A Coherent Data Frame
          dflist <- lapply(output, eval_res_list)
          names(dflist)<- names(output)
          simparamslabels = c("methodName", "ZeroInflate", "RandomEffect","metadataType",
                              "nSubjects", "nPerSubject", "nMicrobes", "spikeMicrobes", "nMetadata",
                              "spikeMetadata", "effectSize", "readDepth", "rep")

          df = make_power_df(dflist, simparamslabels)
          save(df, file=outputStringR)




        } else{
          print("Output file already exists. No new results generated!")
          load(outputStringR)
        }


        stopCluster(cl)

      }
      # # Track End Time
      stop.time <- Sys.time()
      time<-round(difftime(stop.time, start.time, units="min"),3)
      cat(c("Job finished at:",date()), "\n");
      cat("Computational time:",time,"minutes \n")
      cat("The output is in:", outputDirectory, fill=TRUE)

    }
  }

}

