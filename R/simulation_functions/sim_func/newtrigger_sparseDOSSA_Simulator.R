#Code is modified from https://github.com/biobakery/maaslin2_benchmark
newtrigger_sparseDOSSA_Simulator<-function(noZeroInflate=FALSE,
                                           RandomEffect=FALSE,
                                           metadataType= "UVB",
                                           nSubjects=nSubjects,
                                           nPerSubject=1,
                                           nMicrobes=nMicrobes,
                                           spikeMicrobes=0.1,
                                           nMetadata=1,
                                           spikeMetadata=1,
                                           effectSize=effectSize,
                                           #df_spike_metadata=df_spike_metadata,
                                           features_otu=features_otu,
                                           readDepth = readDepth,
                                           SparseDOSSA2_fit=SparseDOSSA2_fit,
                                           #physeqtr=physeqtr,
                                           nIterations = nIterations,
                                           nspike=nspike,
                                           noParallel = FALSE,
                                           rSeed = 1234,
                                           nCores = 4){

  # Create Replicates
  reps = 1:nIterations

  effect_Size=effectSize

  ########################
  # Catch Obvious Errors #
  ########################

  # Check Character Values
  if (!metadataType %in% c('UVA', 'UVB', 'MVA', 'MVB'))
    stop('Must be one of the following: UVA, UVB, MVA, or MVB.')

  # Check Positive Integer Values
  if (round(nSubjects) != nSubjects ||
      nSubjects<0 ||
      round(nPerSubject) != nPerSubject ||
      nPerSubject<0 ||
      round(nMicrobes) != nMicrobes ||
      nMicrobes<0 ||
      round(nMetadata) != nMetadata ||
      nMetadata<0 ||
      round(readDepth) != readDepth ||
      readDepth<0)
    stop('nSubjects/nPerSubject/nMicrobes/nMetadata/readDepth must be positive integers.')

  # Check Proportion Values
  if (spikeMicrobes>1 || spikeMicrobes<=0 || spikeMetadata<=0 || spikeMetadata>1)
    stop('spikeMicrobes/spikeMetadata must be in (0, 1].')

  # Check Illegal Combinations
  if(RandomEffect==TRUE && nPerSubject==1)
    stop('nPerSubject must be greater 1 when RandomEffect is TRUE.')

  if(RandomEffect==FALSE && nPerSubject>1)
    stop('nPerSubject must be equal to  1 when RandomEffect is FALSE.')

  if(metadataType %in% c('UVA', 'UVB') && (nMetadata!=1 || spikeMetadata!=1))
    stop('Both nMetadata and spikeMetadata must be equal to 1 when metadataType is UVA or UVB.')

  if(!metadataType %in% c('UVA', 'UVB') && nMetadata==1)
    stop('nMetadata must be greater than 1 when metadataType is MVA or MVB')

  # Define the Simulation Parameters Combinations
  simparams = apply(expand.grid(metadataType,
                                nSubjects,
                                nPerSubject,
                                nMicrobes,
                                spikeMicrobes,
                                nMetadata,
                                spikeMetadata,
                                effectSize,
                                readDepth,
                                reps), 1, paste, collapse = '_')

  # Define the Labels to Go with Each Element of the Simulation Parameter
  simparamslabels = c("metadataType","nSubjects", "nPerSubject", "nMicrobes", "spikeMicrobes", "nMetadata", "spikeMetadata", "effectSize", "readDepth", "rep")

  # Track Start Time
  cat(c("Job started at:",date()), "\n")
  start.time <- Sys.time()

  # Set Reproducibility Seed
  set.seed(rSeed)

  # Call Grid Computing Only When Specified

  # if (noParallel){
  #
  #   # Call SparseDOSSA Wrapper (noParallel)
  #   simlist <- sparseDOSSA_Wrapper_noParallel(simparams, simparamslabels, SparseDOSSA2_fit, noZeroInflate=noZeroInflate)
  # }
  # else{

  # Set Up Clustering Environment
  no_cores <- nCores
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)


  ####################
  # Data Generation #
  ###################
  keep_filtered = colnames(features_otu)
  keep_fil= sample(keep_filtered,  nspike)


  keep_fil1=keep_fil[1:(nspike/2)]
  keep_fil2=keep_fil[(nspike/2 + 1):length(keep_fil)]

  df_spike_metadata <-
    rbind(
      tibble::tibble(metadata_datum = 1,
                     feature_spiked = keep_fil1,
                     associated_property = "abundance",
                     effect_size = effect_Size),
      tibble::tibble(metadata_datum = 1,
                     feature_spiked = keep_fil2,
                     associated_property = "abundance",
                     effect_size = - effect_Size)
    )

  # Call SparseDOSSA Wrapper
  simlist <- newsparseDOSSA_Wrapper(simparams, simparamslabels, noZeroInflate=noZeroInflate,SparseDOSSA2_fit=SparseDOSSA2_fit,
                                    df_spike_metadata=df_spike_metadata) #, physeqtr=physeqtr)

  # Stop the Cluster
  stopCluster(cl)
  #}

  # Set Names
  if (noZeroInflate==TRUE && RandomEffect==TRUE) {
    simnames<- paste('noZeroInflate_RandomEffect', simparams, sep='_')}
  if (noZeroInflate==TRUE && RandomEffect==FALSE) {
    simnames<- paste('noZeroInflate_noRandomEffect', simparams, sep='_')}
  if (noZeroInflate==FALSE && RandomEffect==TRUE) {
    simnames<- paste('ZeroInflate_RandomEffect', simparams, sep='_')}
  if (noZeroInflate==FALSE && RandomEffect==FALSE) {
    simnames<- paste('ZeroInflate_noRandomEffect', simparams, sep='_')}
  names(simlist) <- simnames

  # Track End Time
  cat(c("Job ended at:",date()), "\n")
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units="min"),3)
  cat("Computational time:", minutes, "minutes \n")

  # Return
  return(simlist)
}
