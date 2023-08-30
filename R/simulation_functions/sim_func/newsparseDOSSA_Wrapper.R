#Code is modified from https://github.com/biobakery/maaslin2_benchmark
# Trigger sparseDOSSA
newsparseDOSSA_Wrapper<-function(simparams, simparamslabels, noZeroInflate, SparseDOSSA2_fit, df_spike_metadata){   #, physeqtr
  f<-foreach(i = simparams, .packages = c("SparseDOSSA2", "MASS", "stringi","zinbwave","SummarizedExperiment",
                                          "dplyr","readr","tibble","BiocParallel","edgeR","DESeq2","phyloseq","ape"),
             .export = c("generateMetadata","SummarizedExperiment"), .errorhandling = 'remove') %dopar% {

               # Extract Parameter Strings
               params = strsplit(i, '_')[[1]]
               names(params) <- simparamslabels

               # Extract Relevant Parameters
               metadataType = as.character(params["metadataType"]) # Type of Metadata
               nSubjects <- as.numeric(params["nSubjects"])  # Number of Subjects
               nPerSubject <- as.numeric(params["nPerSubject"])  # Number of Samples Per Subject
               nSamples<-round(nSubjects*nPerSubject) # Number of Samples
               nMicrobes <- as.numeric(params["nMicrobes"])  # Number of Microbes
               spikeMicrobes <- as.numeric(params["spikeMicrobes"]) # Proportion of Spiked-in Microbes
               nMetadata<-as.numeric(params["nMetadata"])  # Number of Metadata
               spikeMetadata<-as.numeric(params["spikeMetadata"])  # Proportion of Spiked-in Metadata
               effectSize<-as.numeric(params["effectSize"]) # Effect Size
               readDepth<-as.numeric(params["readDepth"]) # Library Size


               # # Initialize
               DSim = NULL
               #
               # # sparseDOSSA Error Control
               tryAgain = TRUE
               infiniteloopcounter = 1
               while (tryAgain & infiniteloopcounter < 5) {
               #nSubjects=nSamples=100
               # Generate Metadata
                 FF<-generateMetadata(metadataType=metadataType,
                                      nSubjects=nSubjects,
                                      nPerSubject=nPerSubject,
                                      nMetadata=nMetadata,
                                      spikeMetadata=spikeMetadata)

                 # Extract Relevant Information
                 UserMetadata<-FF$UserMetadata;
                 Metadatafrozenidx<-FF$Metadatafrozenidx;
                 significant_metadata<-FF$significant_metadata;
                 spikeCount<-FF$spikeCount


                 mat_metadata <- t(UserMetadata)
                 # mat_metadata <- matrix(c(rep(1, nSubjects/2),
                 #                            rep(0,nSubjects/2)),ncol = 1)


                 DSim <- SparseDOSSA2::SparseDOSSA2(template = SparseDOSSA2_fit,
                                                    new_features = FALSE,
                                                    n_sample = nSamples,
                                                    #spike_metadata = "abundance",
                                                    metadata_matrix = mat_metadata,
                                                    spike_metadata = df_spike_metadata,
                                                    median_read_depth= readDepth, #1438,   #250000,
                                                    verbose = TRUE)


                   if (is.null(DSim) | inherits(DSim, "try-error")) {
                     tryAgain = TRUE
                     infiniteloopcounter = infiniteloopcounter + 1
                   } else {
                     tryAgain = FALSE
                   }
               }
               if (infiniteloopcounter >= 5) {
                 stop("Consistent error found during simulation. Need to investigate cause.")
               }

               # Gather sparseDOSSA Outputs
               sparsedossa_results <- as.data.frame(DSim$simulated_data)
               #dim(sparsedossa_results)
               colnames(sparsedossa_results)<-paste('Sample', 1:ncol(sparsedossa_results), sep='')

               if(abs(effectSize) > 0){
                 truth<-c(unlist(DSim$spike_metadata$spike_metadata$feature_spiked))
                 truth<-truth[!stringi::stri_detect_fixed(truth,":")]
                 #truth<-truth[(5+nMetadata):length(truth)]
                 #truth<-as.data.frame(truth)
                 significant_features<- truth #as.vector(truth) #[seq(1, (as.numeric(spikeCount)+1)*(nMicrobes*spikeMicrobes), (as.numeric(spikeCount)+1)),])
               } else{
                 truth<- "NULL" #c(unlist(DD$spike_metadata$spike_metadata$feature_spiked))
                 truth<-truth[!stringi::stri_detect_fixed(truth,":")]
                 #truth<-truth[(5+nMetadata):length(truth)]
                 #truth<-as.data.frame(truth)
                 significant_features<- truth #as.vector(truth) #[seq(1, (as.numeric(spikeCount)+1)*(nMicrobes*spikeMicrobes), (as.numeric(spikeCount)+1)),])
              }


               #rownames(sparsedossa_results)<-sparsedossa_results$X1
               #sparsedossa_results<-sparsedossa_results[-1,-1]
               data<- mat_metadata #as.matrix(sparsedossa_results[-c((nMetadata+1):(2*nMicrobes+nMetadata)),])
               colnames(data)="Metadata1"
               data<-data.matrix(data)
               class(data) <- "numeric"


               # Separate Metadata and Taxa

               # Extract Metadata
               if (metadataType %in% c('UVA', 'UVB')){
                 metadata<- as.data.frame(data) #as.data.frame(data[1,])
                 rownames(metadata)<-colnames(sparsedossa_results)
                 #colnames(metadata)<-rownames(data)[1]
               } else{
                 metadata<-data #as.data.frame(t(data[(1:nMetadata),]))
               }

               # Mark True Positive Metadata - Same Format at Mcmurdie and Holmes (2014)
               which.TP = colnames(metadata) %in% significant_metadata
               meta_newname = paste0(colnames(metadata)[which.TP], "_TP")
               colnames(metadata)[which.TP] <- meta_newname

               # Extract Features
               features<- as.data.frame(t(sparsedossa_results)) #as.data.frame(t(data[-c(1:nMetadata),]))

               # Mark True Positive Features - Same Format at Mcmurdie and Holmes (2014)
               wh.TP = colnames(features) %in% significant_features
               nMicrobesSim=length(colnames(features))
               colnames(features)<-paste("Feature", 1:nMicrobesSim, sep = "")
               newname = paste0(colnames(features)[wh.TP], "_TP")
               colnames(features)[wh.TP] <- newname;

               # Add Sample ID
               ID<-rep(paste('Subject', 1:nSubjects, sep=''), each = nPerSubject)

               # Add Library Size (or Sequencing Depth)
               #libSize<-rowSums(features)

               #dim(features)

               group=factor(metadata$Metadata1_TP)
               ugroup=unique(group)
               df_counts = data.frame(features,group)

               frepi=list()
               for(i in 1:length(ugroup)){
                 bhasf=subset(df_counts, group == ugroup[i])
                 dim(bhasf)
                 bhasfs=bhasf[,-dim(df_counts)[2]]
                 kp = colSums(bhasfs >= 2) >= 2
                 d2k=bhasfs[,kp]
                 dim(d2k)
                 frep=colnames(d2k)
                 frepi[[i]]=frep
               }
               fcolst <- unique(unlist(frepi))
               length(fcolst)
               features1 <- features[,fcolst]


               group=factor(metadata$Metadata1_TP)
               ugroup=unique(group)
               df_counts1 = data.frame(features1,group)

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
               drop1 <- filt_fcols

               if(length(drop1)==0){
                 ex_counts1 <- NULL
                 } else {
                 ex_counts1 <- features1[,drop1]
                 }

               #features retained for further analysis
               retain_counts = features1[,!(colnames(features1) %in% drop1)]



               design <- model.matrix(~Metadata1_TP, data = metadata)

               #dds_zinbwave <- DESeqDataSetFromMatrix(countData = t(filtered_counts), colData = metadata, design = design)
               dds_zinbwave <- DESeqDataSetFromMatrix(countData = t(retain_counts), colData = metadata, design = design)

               dds_zinbwave <- zinbwave(dds_zinbwave,
                                        X="~ Metadata1_TP", #Metadata1_TP orobanchol_pmol_g",
                                        epsilon = 1e10,
                                        verbose = FALSE,
                                        K = 0,
                                        observationalWeights = TRUE,
                                        BPPARAM = BiocParallel::SerialParam())

               weights <- assay(dds_zinbwave, "weights")



               featuresEx=as.data.frame(ex_counts1)
               features=as.data.frame(retain_counts)

               # Add Library Size (or Sequencing Depth)
               libSize<-rowSums(features)

               return(list(metadata=metadata, features=features, featuresEx=featuresEx, ID=ID, libSize=libSize,weights=weights))
             }
  return(f)
}
