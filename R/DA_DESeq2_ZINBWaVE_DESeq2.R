DA_DESeq2_ZINBWaVE_DESeq2 <- function( taxa_All_filt,taxa_NonZeroGroup, taxa_ZeroGroup, group,
                                       taxadata, reduced=~1, normalization = c("GPRM","poscounts"))
{
  #normalization = "GPRM"
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
  } else {
    dds_all <- phyloseq_to_deseq2(ps_all, ~ group )
    dds_all <- estimateSizeFactors(dds_all, type="poscounts")
    size.factor <- dds_all@colData@listData$sizeFactor
  }

  ### Taxa in the NonZero group
  if(!is.null(taxa_ZeroGroup)){
    ps_zs <- prune_taxa(colnames(taxa_ZeroGroup), ps_all)

    dds <- phyloseq_to_deseq2(ps_zs, ~ group )
    sizeFactors(dds) <- size.factor

    dds <- DESeq(dds, test="LRT", reduced=reduced,
                 minmu=1e-6, minReplicatesForReplace=Inf)

    coef_tested= resultsNames(dds)

    resP <- results(dds, name =  coef_tested[2])
    resP_all <- resP
    #summary(resP)
    if(!is.null(resP)){
      resP = resP[order(resP$padj, na.last=NA), ]
      resPP=resP
      resP = subset(resP, padj < 0.05)
      taxa_sigP = rownames(resP)


      df_res_P <- cbind(as(resPP, "data.frame"), as(tax_table(ps_zs)[rownames(resPP), ], "matrix"))
      df_res_P <- df_res_P %>%
        rownames_to_column(var = "OTU") %>%
        arrange(padj)

      rownames(df_res_P)= df_res_P$OTU

      df_res_PP <- df_res_P %>%
        dplyr::filter(abs(log2FoldChange) > 0 & padj < 0.05) %>%
        droplevels() %>%
        arrange(padj) #%>%
        #column_to_rownames("OTU")
    } else{
      df_res_P <- NULL
      df_res_PP <- NULL
    }

  } else{
    df_res_P <- NULL
    df_res_PP <- NULL
    resP_all <- NULL
  }

  ## Taxa in the NonZero group
  ps_nzs <- prune_taxa(colnames(taxa_NonZeroGroup), ps_all)


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

  #weights[weights <= 10^-6] = 10^-6
  dds_zinb <- DESeqDataSet(dds_zinb, design = ~group)
  #dds_zinb <- estimateSizeFactors(dds_zinb, type="poscounts")
  sizeFactors(dds_zinb) <- size.factor

  dds_zinb <- DESeq(dds_zinb, test="LRT", reduced=reduced,
                    minmu=1e-6, minReplicatesForReplace=Inf)


  coef_tested_zinb= resultsNames(dds_zinb)
  resW <- results(dds_zinb, name =  coef_tested_zinb[2])
  resW_all <- resW
  #summary(resW)
  if(!is.null(resW)){
    resW = resW[order(resW$padj, na.last=NA), ]
    resWW = resW
    resW = data.frame(subset(resW, padj < 0.05))
    taxa_sigW = rownames(resW)

    df_res_W <- cbind(as(resWW, "data.frame"), as(tax_table(ps_nzs)[rownames(resWW), ], "matrix"))
    df_res_W <- df_res_W %>%
      rownames_to_column(var = "OTU") %>%
      arrange(padj)

    rownames(df_res_W)= df_res_W$OTU

    df_res_WW <- df_res_W %>%
      dplyr::filter(abs(log2FoldChange) > 0 & padj < 0.05) %>%
      droplevels() %>%
      arrange(padj) #%>%
      #column_to_rownames("OTU")

  } else{
    df_res_W <- NULL
    df_res_WW <- NULL
  }

  rownames(df_res_PP) <- df_res_PP$OTU
  rownames(df_res_WW) <- df_res_WW$OTU

  return=list(sign_taxa_ZeroGroup=df_res_PP, sign_taxa_NonZeroGroup=df_res_WW, results_taxa_ZeroGroup=resP_all,
              results_taxa_NonZeroGroup=resW_all, ps_all=ps_all)


}
