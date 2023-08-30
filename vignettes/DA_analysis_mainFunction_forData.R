---
title: "DA analysis of microbiome data"
author: "Yazew"
date: "29/06/2023"
output: 
#BiocManager::install("BiocStyle")
  BiocStyle::html_document:
    toc: yes
    number_section: yes
    fig_caption: yes
vignette: >
  %\VignetteEcoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Intro}
  \usepackage[utf8]{inputenc}
#bibliography: bib_intro.json
#csl: bioinformatics.csl
---

```{r options, include=FALSE, echo=FALSE}
knitr::opts_chunk$set(
    warning = FALSE, 
    error = FALSE, 
    message = FALSE)
```

#Load libraries
```{r, echo=FALSE}
library(phyloseq)          
library(DESeq2)             
library(zinbwave)
library(apeglm)    
library(parallel)
BiocParallel::register(BiocParallel::SerialParam())
library(readr)
library(tibble)
library(dplyr)
library(tidyverse)  
library(ape)
library(GUniFrac)
library(vegan)
#library(devtools)
#install_github("lichen-lab/GMPR")
library(GMPR)
library(ggplot2)
library(RColorBrewer)
library(randomcoloR)
library(pheatmap)
library(reshape2)
#library(ggarrange)
```

#First run the followings three functions
```{r, echo=FALSE}
source("")
```   

#Load phyloseq data 
```{r}
Data_phyloseq <- readRDS("C:\\Users\\fenta\\Documents\\OneDrive\\Documents\\2021\\Amsterdam\\MicrobiomeProject\\CodesR\\Simulation_ageJohFred\\maaslin2_benchmark-master\\mocksamples\\fredMock.RDS")
Data_phyloseq
counts <- t(otu_table(Data_phyloseq, 'matrix'))
group <- sample_data(Data_phyloseq)$Soil_Type 
taxadata = tax_table(Data_phyloseq)

```

#Groupwise filtering
```{r}
g_filtered <- group_wise_filer(group=group, counts=counts, min_counts=2, min_replicates=2)
taxa_All_filt <- g_filtered$taxa_All_filt
cat("Total number of filtered taxa = ", dim(taxa_All_filt)[2],"\n")
taxa_NonZeroGroup <- g_filtered$taxa_NonZeroGroup
cat("Number of taxa in the Non-Zero Group = ", dim(taxa_NonZeroGroup)[2],"\n")
taxa_ZeroGroup <- g_filtered$taxa_ZeroGroup
cat("Number of taxa in the Zero Group = ",dim(taxa_ZeroGroup)[2],"\n")
```

#Differential abundance analysis
```{r}


results_D2ZD2 <- DA_DESeq2_ZINBWaVE_DESeq2(taxa_All_filt,taxa_NonZeroGroup, taxa_ZeroGroup, group, 
                                           taxadata, reduced=~1, normalization="GPRM") 


#Table of significant taxa
taxa_sig_All_names = union(rownames(results_D2ZD2$sign_taxa_ZeroGroup),
                     rownames(results_D2ZD2$sign_taxa_NonZeroGroup))
results_taxa_ZeroGroup <-results_D2ZD2$results_taxa_ZeroGroup
results_taxa_NonZeroGroup <-results_D2ZD2$results_taxa_NonZeroGroup

results_taxa_ALL <- rbind(results_taxa_ZeroGroup, results_taxa_NonZeroGroup)
results_taxa_ALL
```
#Plots (all)

```{r,fig.width=15, fig.height=10}
# table of significant taxa

taxa_sigPWLR = union(rownames(results_D2ZD2$sign_taxa_ZeroGroup),
                     rownames(results_D2ZD2$sign_taxa_NonZeroGroup))
ps_all <- results_D2ZD2$ps_all
results_taxa_ZeroGroup <-results_D2ZD2$results_taxa_ZeroGroup
results_taxa_NonZeroGroup <-results_D2ZD2$results_taxa_NonZeroGroup

Da_plots <- plot_results(taxa_sigPWLR, results_taxa_ZeroGroup, results_taxa_NonZeroGroup, ps_all )
Da_plots
# Da_plots$plot_family
# Da_plots$plot_genus
# Da_plots$plot_counts

```
#Extracting and adjusting width and height of individaul plots: heatmap
```{r,fig.width=15, fig.height=10}
hmap <- ggpubr::ggarrange(Da_plots$plot_heatmap$gtable,  nrow = 1)
hmap
```
#Extracting and adjusting width and height of individaul plots:logfold change plot
```{r,fig.width=10, fig.height=6}
Logfold_family <- Da_plots$plot_family 
Logfold_family
```
