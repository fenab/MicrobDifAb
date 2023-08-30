plot_results <- function(taxa_sigPWLR=taxa_sigPWLR, sign_taxa_ZeroGroup=sign_taxa_ZeroGroup,
                         sign_taxa_NonZeroGroup=sign_taxa_NonZeroGroup, ps_all=ps_all, group=group )
{

  ps.taxa.pse.sub <- ps_all

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
          axis.ticks.x=element_blank(), text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + #c(0.8, 0.2)) +
    xlab("Group") +
    theme(legend.position = "top")

  #Logfold change plots

  df_resPW = rbind(sign_taxa_ZeroGroup, sign_taxa_NonZeroGroup)
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
      averageCounts[baseMean < 100] <- "<100"
      averageCounts[baseMean >= 100 & baseMean < 500] <- "100-500"
      averageCounts[baseMean >= 500 & baseMean < 1000] <- "500-1000"
      averageCounts[baseMean > 1000] <- ">1000"
    } )
    df_resPW$averageCounts=factor(df_resPW$averageCounts,levels=c("<100","100-500","500-1000",">1000"))
    # df_resPW <- within(df_resPW, {
    #   averageCounts <- NA # need to initialize variable
    #   averageCounts[baseMean < round(quantile(baseMean,0.25),0)] <- paste0(c("<",round(quantile(baseMean,0.33),0)), collapse="")
    #   averageCounts[baseMean >= round(quantile(baseMean,0.33),0) & baseMean <= round(quantile(baseMean,0.66),0)] <- paste0(c(round(quantile(baseMean,0.33),0),round(quantile(baseMean,0.66),0)), collapse="-")
    #   averageCounts[baseMean > round(quantile(df_resPW$baseMean,0.66),0)] <- paste0(c(">",round(quantile(df_resPW$baseMean,0.66),0)), collapse="")
    # } )
    # df_resPW$averageCounts=factor(df_resPW$averageCounts,levels=c(paste0(c("<",round(quantile(df_resPW$baseMean,0.33),0)), collapse=""),paste0(c(round(quantile(df_resPW$baseMean,0.33),0),round(quantile(df_resPW$baseMean,0.66),0)), collapse="-"),paste0(c(">",round(quantile(df_resPW$baseMean,0.66),0)), collapse="")))

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
  pt_genus <- ggplot(df_resPW, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
    geom_point(size=3) +
    #scale_x_discrete(labels = label_wrap(10)) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
    theme(text = element_text(size = 12, face="bold", family="serif"), axis.line = element_line(colour = "black"))  +
    theme(axis.text.x = element_text(face="bold", size=10),
          axis.text.y = element_text(face="bold", size=10)) +
    geom_point(aes(size = averageCounts, color=Phylum))

  #pt_genus <- pt_genus +  aes(stringr::str_wrap(Genus, 10), log2FoldChange) + xlab(NULL)



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
  pt_family <- ggplot(df_resPW, aes(x=Family, y=log2FoldChange, color=Phylum)) +
    geom_point(size=3) +
    #scale_x_discrete(labels = label_wrap(10)) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
    theme(text = element_text(size = 12, face="bold", family="serif"), axis.line = element_line(colour = "black"))  +
    theme(axis.text.x = element_text(face="bold", size=10),
          axis.text.y = element_text(face="bold", size=10)) +
    geom_point(aes(size = averageCounts, color=Phylum))

  #pt_family <- pt_family +  aes(stringr::str_wrap(Family, 10), log2FoldChange) + xlab(NULL)

  return(list(plot_heatmap=ph, plot_family=pt_family, plot_genus=pt_genus, plot_counts=pt_counts))
}
