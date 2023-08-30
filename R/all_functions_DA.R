group_wise_filer  <- function(group, counts, min_counts, min_replicates){
  
  ugroup <- unique(group)
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
  ex_counts <- features[,filtZero_fcols] 
  
  ffeatures1 = features[,!(colnames(features) %in% filtZero_fcols)]
  
  return=list(taxa_All_filt=features, taxa_NonZeroGroup=ffeatures1, taxa_ZeroGroup=ex_counts)
  
}

DA_DESeq2_ZINBWaVE_DESeq2 <- function( taxa_All_filt,taxa_NonZeroGroup, taxa_ZeroGroup, group,
                                       taxadata, reduced=~1, normalization = c("GPRM","poscounts")) 
{ 
  
  #Normalization
  group = data.frame(group)
  rownames(group) <- rownames(taxa_All_filt)
  SAM1f <- sample_data(group) 
  OTU = otu_table(as.matrix(t(taxa_All_filt)), taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(taxadata))
  
  TAX1 <- TAX[rownames(OTU),]
  all(rownames(OTU)==rownames(TAX1))
  physeq1f = phyloseq(OTU, TAX1)
  
  # merge the data
  ps_all = merge_phyloseq(physeq1f, SAM1f) 
  ps_all
  
  
  sample_data(ps_all)$group <- as.factor(sample_data(ps_all)$group) 
  
  #normalizing data
  if(normalization == "GPRM") {
    size.factor=GMPR(OTUmatrix=data.frame(taxa_All_filt))
  }
  else {
    dds_all <- phyloseq_to_deseq2(ps_zs, ~ group ) 
    size.factor <- estimateSizeFactors(dds, type="poscounts")
    
  }
  
  ### Taxa in the NonZero group
  
  ps_zs <- prune_taxa(names(taxa_ZeroGroup), ps_all)
  
  dds <- phyloseq_to_deseq2(ps_zs, ~ group ) 
  sizeFactors(dds) <- size.factor
  
  dds <- DESeq(dds, test="LRT", reduced=reduced,    
               minmu=1e-6, minReplicatesForReplace=Inf) 
  
  coef_tested= resultsNames(dds)
  
  resP <- results(dds, name =  coef_tested[2]) 
  #summary(resP)
  if(!is.null(resP)){
      resP = resP[order(resP$padj, na.last=NA), ]
      resPP=resP
      resP = subset(resP, padj < 0.05)
      taxa_sigP = rownames(resP) 
      
      
      df_res_P <- cbind(as(resP, "data.frame"), as(tax_table(ps_zs)[rownames(resP), ], "matrix"))
      df_res_P <- df_res_P %>%
        #rownames_to_column(var = "OTU") %>%
        arrange(padj)
      
      #rownames(df_res_P)= df_res_P$OTU
      
      df_res_PP <- df_res_P %>%
        filter(abs(log2FoldChange) > 0 & padj < 0.05) %>%
        droplevels() %>%
        arrange(padj) #desc(log2FoldChange)
  }
  ## Taxa in the NonZero group
  
  else{
    df_res_P = NULL 
    df_res_PP = NULL 
  }
  ps_nzs <- prune_taxa(names(taxa_NonZeroGroup), ps_all)
  
  
  #generating observation weights using ZINBWaVE.
  
  dds_zinb <- phyloseq_to_deseq2(ps_nzs, ~ group )
  dds_zinb <- zinbwave(dds_zinb,
                       X="~ group", 
                       epsilon = 1e10,
                       verbose = FALSE,
                       K = 0,
                       observationalWeights = TRUE,
                       BPPARAM = BiocParallel::SerialParam())
  
  weights = assay(dds_zinb, "weights")
  
  weights[weights <= 10^-6] = 10^-6
  dds_zinb <- DESeqDataSet(dds_zinb, design = ~group) 
  #dds_zinb <- estimateSizeFactors(dds_zinb, type="poscounts")
  sizeFactors(dds_zinb) <- size.factor
  
  dds_zinb <- DESeq(dds_zinb, test="LRT", reduced=reduced,     
                    minmu=1e-6, minReplicatesForReplace=Inf) 
  
  
  coef_tested_zinb= resultsNames(dds_zinb)
  resW <- results(dds_zinb, name =  coef_tested_zinb[2]) 
  #summary(resW)
  if(!is.null(resW)){
      resW = resW[order(resW$padj, na.last=NA), ]
      resWW = resW 
      resW = data.frame(subset(resW, padj < 0.05))
      taxa_sigW = rownames(resW) 
      
      df_res_W <- cbind(as(resW, "data.frame"), as(tax_table(ps_nzs)[rownames(resW), ], "matrix"))
      df_res_W <- df_res_W %>%
        #rownames_to_column(var = "OTU") %>%
        arrange(padj)
      
      #rownames(df_res_W)= df_res_W$OTU
      
      df_res_WW <- df_res_W %>%
        filter(abs(log2FoldChange) > 0 & padj < 0.05) %>%
        droplevels() %>%
        arrange(padj)
  }
  
  else{
    df_res_W = NULL 
    df_res_WW = NULL 
  }
  
  return=list(sign_taxa_ZeroGroup=df_res_PP, sign_taxa_NonZeroGroup=df_res_WW, results_taxa_ZeroGroup=df_res_P,
              results_taxa_NonZeroGroup=df_res_W, ps_all=ps_all)
  
  
}

#Plotingresults


plot_results <- function(taxa_sigPWLR, results_taxa_ZeroGroup, results_taxa_NonZeroGroup, ps_all )
{
  
  ps.taxa.pse.sub <- results_D2ZD2$ps_all
  
  ddsall <- phyloseq_to_deseq2(ps.taxa.pse.sub, ~ group  ) 
  #ddsall <- estimateSizeFactors(ddsall, type="poscounts")
  
  #scr0all <- computeSumFactors(ddsall)         
  ## Warning in FUN(...): encountered negative size factor estimates
  
  #Size_Factor = sizeFactors(scr0all)
  #sizeFactors(ddsall) <- sizeFactors(scr0all) #Size_Factor #
  
  ps.taxa.rel <- transform_sample_counts(ps.taxa.pse.sub, function(x) x/sum(x)*100)
  ps.taxa.rel.sig <- prune_taxa(taxa_sigPWLR, ps.taxa.rel) #sig_taxa selected_sigtaxa
  
  # Only keep N and C samples
  ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps.taxa.pse.sub)), ps.taxa.rel.sig)
  #ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps.taxa)), ps.taxa.rel.sig)
  
  
  matrix1 <- otu_table(ps.taxa.rel.sig, matrix) #as.matrix(data.frame())
  #dim(matrix1)
  #rownames(matrix1) <- paste(as.character(tax_table(ps.taxa.rel.sig)[, "Species"]),rownames(matrix1))
  metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig))
  
  # Define the annotation color for columns and rows
  annotation_col = data.frame(
    group = as.factor(metadata_sub$group), 
    check.names = FALSE
  )
  rownames(annotation_col) = rownames(metadata_sub)
  
  annotation_row = data.frame(
    Phylum = as.factor(tax_table(ps.taxa.rel.sig)[, "Phylum"])
  )
  rownames(annotation_row) = rownames(matrix1)
  
  
  
  
  set.seed(123456)
  #palette <- distinctColorPalette(n)
  phylum_col =  distinctColorPalette(length(levels(annotation_row$Phylum)))
  
  names(phylum_col) = levels(annotation_row$Phylum)
  group_col =  distinctColorPalette(length(levels(annotation_col$group)))
  names(group_col) = levels(annotation_col$group)
  ann_colors = list(
    group = group_col, # c(levels(group)[1] = "green", levels(group)[2]  = "purple"), 
    Phylum = phylum_col
  ) 
  #ComplexHeatmap
  ph <- pheatmap::pheatmap(matrix1, scale= "row", cluster_cols = F,
                           annotation_col = annotation_col, 
                           annotation_row = annotation_row, 
                           annotation_colors = ann_colors,
                           fontsize = 10, fontsize_row = 10, fontsize_col = 10
  )
  
  
  #Count plots
  top10_sigOE_genes = taxa_sigPWLR 
  
  ## normalized counts for top 20 significant genes
  
  normalized_counts <- DESeq2::counts(ddsall, normalized=F)
  #top10_sigOE_genes = colnames(ffeatures1) rownames(normalized_counts)[1:83]
  top10_sigOE_norm <- normalized_counts[top10_sigOE_genes, ]
  #which(top10_sigOE_genes %in% rownames(normalized_counts))
  
  
  ## use melt to modify the format of the data frame
  melted_top10_sigOE <- data.frame(melt(top10_sigOE_norm))
  
  ## add column names that make sense
  colnames(melted_top10_sigOE) <- c("Feature", "samplename", "normalized_counts")
  
  ## add metadata to "melted" dataframe
  metadata_sub <- data.frame(sample_data(ps.taxa.pse.sub))
  metadata_sub$samplename <- rownames(metadata_sub)
  
  melted_top10_sigOE <- merge(melted_top10_sigOE, metadata_sub)
  
  
  pt_counts <- ggplot(melted_top10_sigOE, aes(x = as.factor(Feature), y = normalized_counts)) +
    geom_boxplot(aes(fill = group), position = position_dodge(0.9)) +
    geom_point(aes(fill = group),position = position_jitterdodge(0.2)) +
    # theme(axis.text.x=element_blank(),
    #       axis.ticks.x=element_blank() ) +
    facet_wrap(~ Feature, scales = "free") +
    stat_summary(fun=mean, geom="point", aes(group=group), position=position_dodge(.9), color="red", size=2) +
    #scale_fill_manual(values = c(BoraForest = "purple", KeyGene = "blue", Enza='green')) + #c("#09E359", "#E31009", "#E31009")) + 
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), text = element_text(size = 14), panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + #c(0.8, 0.2)) +
    xlab("Soil type") +
    theme(legend.position = "top") 
  
  #Logfold change plots
  
  df_resPW = rbind(results_taxa_ZeroGroup, results_taxa_NonZeroGroup)
  dim(df_resPW)
  if(max(df_resPW$baseMean)<50){
    df_resPW <- within(df_resPW, {   
      averageCounts <- NA # need to initialize variable
      averageCounts[baseMean < 10] <- "<10"
      averageCounts[baseMean >= 10 & baseMean <= 30] <- "10-30"
      averageCounts[baseMean > 30] <- ">30"
    } )
    df_resPW$averageCounts=factor(df_resPW$averageCounts,levels=c("<10","10-30",">30"))
  }
  else{
    df_resPW <- within(df_resPW, {
      averageCounts <- NA # need to initialize variable
      averageCounts[baseMean < round(quantile(baseMean,0.25),0)] <- paste0(c("<",round(quantile(baseMean,0.33),0)), collapse="")
      averageCounts[baseMean >= round(quantile(baseMean,0.33),0) & baseMean <= round(quantile(baseMean,0.66),0)] <- paste0(c(round(quantile(baseMean,0.33),0),round(quantile(baseMean,0.66),0)), collapse="-")
      averageCounts[baseMean > round(quantile(df_resPW$baseMean,0.66),0)] <- paste0(c(">",round(quantile(df_resPW$baseMean,0.66),0)), collapse="")
    } )
    df_resPW$averageCounts=factor(df_resPW$averageCounts,levels=c(paste0(c("<",round(quantile(df_resPW$baseMean,0.33),0)), collapse=""),paste0(c(round(quantile(df_resPW$baseMean,0.33),0),round(quantile(df_resPW$baseMean,0.66),0)), collapse="-"),paste0(c(">",round(quantile(df_resPW$baseMean,0.66),0)), collapse="")))
    
  }
  
  
  # df_resPW = df_resPW %>%
  #    mutate(Genus = recode(Genus, 'Burkholderia-Caballeronia-Paraburkholderia' = 'Burkholderia-C-P'))
  #  df_resPW = df_resPW %>%
  #    mutate(Genus = recode(Genus, 'Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = 'Allorhizobium-N-Par-Rh'))
  #  
  
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  # Phylum order
  x = tapply(df_resPW$log2FoldChange, df_resPW$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  df_resPW$Phylum = factor(as.character(df_resPW$Phylum), levels=names(x))
  # Genus order
  x = tapply(df_resPW$log2FoldChange, df_resPW$Genus, function(x) max(x))
  x = sort(x, TRUE)
  df_resPW$Genus = factor(as.character(df_resPW$Genus), levels=names(x))
  pt_genus <- ggplot(df_resPW, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
    theme(text = element_text(size = 16, face="bold", family="serif"), axis.line = element_line(colour = "black"))  +
    theme(axis.text.x = element_text(face="bold", size=14),
          axis.text.y = element_text(face="bold", size=14)) +
    geom_point(aes(size = averageCounts, color=Phylum))
  
  
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  # Phylum order
  x = tapply(df_resPW$log2FoldChange, df_resPW$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  df_resPW$Phylum = factor(as.character(df_resPW$Phylum), levels=names(x))
  # Genus order
  x = tapply(df_resPW$log2FoldChange, df_resPW$Family, function(x) max(x))
  x = sort(x, TRUE)
  df_resPW$Family = factor(as.character(df_resPW$Family), levels=names(x))
  pt_family <- ggplot(df_resPW, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
    theme(text = element_text(size = 16, face="bold", family="serif"), axis.line = element_line(colour = "black"))  +
    theme(axis.text.x = element_text(face="bold", size=14),
          axis.text.y = element_text(face="bold", size=14)) +
    geom_point(aes(size = averageCounts, color=Phylum))
  
  return(list(plot_heatmap=ph, plot_family=pt_family, plot_genus=pt_genus, plot_counts=pt_counts))
}
