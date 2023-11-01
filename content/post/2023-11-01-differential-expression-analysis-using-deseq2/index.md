---
title: Differential expression analysis using DEseq2
author: Package Build
date: '2023-11-01'
slug: []
categories: []
tags: []
---




I borrow the data from this article ( https://www.pnas.org/doi/10.1073/pnas.2008762117#data-availability) to analyze Alzheimer’s disease (AD) data, aiming to show the standard analysis during my work. The previous blog mainly focused on the quality control, dimensionality reduction, clustering, and cell type annotation of single-nuclei RNA seq data, and further analyzes the data for differential abundance analysis. Here, I apply DEseq2 for differential expression analysis to identify significant genes and explore their functions.


### Libraries

```r
library(Seurat)
library(scuttle)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(annotables)
library(clusterProfiler)
library(org.Hs.eg.db)
```

### Data

```r
seurat <- readRDS("../../../../../transcriptome_protein_AD_plat/sc_RNA/Seu_PCA_UMAP_annotated.RDS")
DimPlot(seurat, reduction = "umap",label = TRUE)+NoLegend()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-1.png" width="672" />

```r
DimPlot(seurat, reduction = "umap", group.by = "group_id")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-2.png" width="672" />


```r
DimPlot(seurat, reduction = "umap", split.by = "orig.ident")+NoLegend()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-3-1.png" width="960" />

From above plots, we can see that cells in NC11 is much lesser than other samples. Here, I assume NC11 as a outlier.


### Creating pseudo-bulk samples


```r
sce <- as.SingleCellExperiment(seurat)
summed <- aggregateAcrossCells(sce, 
    ids=colData(sce)[,c("orig.ident", "cell_type")])
summed <- summed[,summed$orig.ident != "NC11"]
```



### DE analysis and visualization

```r
dir.create("DE_res")
ht_list <- list()
pca_list <- list()
for (i in unique(summed$cell_type)) {
  current <- summed[,summed$cell_type == i]
  ## Remove lowly expressed genes which have less than 6 cells with any counts
  current <- current[rowSums(counts(current) > 1) >= 5,]
  cluster_counts <- counts(current)
  cluster_metadata <- data.frame(colData(current))
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ group_id)
  # Transform counts for data visualization
  # rld <- rlog(dds, blind=TRUE)
  

  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  compare_list <- c("AD","NC")
 
    contrast <- c("group_id", compare_list)
    res.all <- results(dds, contrast = contrast) %>% data.frame()
    write.csv(res.all, paste0("DE_res/",i,".csv"))
    
    res.other <- res.all[res.all$pvalue >= 0.001,]
    
  # Here select 20% insignificant genes randomly for plot. 
    set.seed(123)
    res.sel <- res.other[sample(1:nrow(res.other), nrow(res.other) * 0.2),]
    res <- rbind(res.all[res.all$pvalue < 0.001,], res.sel)
    ## volcanoplot
    p <- EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim = c(min(res$log2FoldChange, na.rm = TRUE) - 1,max(res$log2FoldChange, na.rm = TRUE) + 1),
    ylim = c(0, max(-log10(res$pvalue), na.rm = TRUE) + 1),
    title = i,
    pCutoff = 0.001,
    FCcutoff = 1,
    pointSize = 1.4,
    labSize = 2.3,
    axisLabSize = 8,
    titleLabSize = 10,
    subtitleLabSize = 8,
    captionLabSize = 8,
    legendLabSize = 7,
    legendIconSize = 2,
    # max.overlaps = Inf,
    drawConnectors = TRUE,
    widthConnectors = 0.3)
    
    print(p)
  
}
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-1.png" width="672" /><img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-2.png" width="672" /><img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-3.png" width="672" /><img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-4.png" width="672" /><img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-5.png" width="672" /><img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-6.png" width="672" /><img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-7.png" width="672" /><img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-8.png" width="672" />

Volcano plots in each cell type. Here cutoff for the significant genes is: P value < 0.001, log2FC > 1.


### Functional analysis

Here, I explore the significant genes in "Interneuron; GABAergic neuron"  cell type. Because the  previous post in https://xiwen514.github.io/post/2023-10-30-diff-abun-analysis/ shows that the p value of this cell type is the lowerst, 0.1.


```r
head(grch37)
```

```
## # A tibble: 6 × 9
##   ensgene         entrez symbol   chr       start     end strand biotype descr…¹
##   <chr>            <int> <chr>    <chr>     <int>   <int>  <int> <chr>   <chr>  
## 1 ENSG00000000003   7105 TSPAN6   X      99883667  9.99e7     -1 protei… tetras…
## 2 ENSG00000000005  64102 TNMD     X      99839799  9.99e7      1 protei… tenomo…
## 3 ENSG00000000419   8813 DPM1     20     49551404  4.96e7     -1 protei… dolich…
## 4 ENSG00000000457  57147 SCYL3    1     169818772  1.70e8     -1 protei… SCY1-l…
## 5 ENSG00000000460  55732 C1orf112 1     169631245  1.70e8      1 protei… chromo…
## 6 ENSG00000000938   2268 FGR      1      27938575  2.80e7     -1 protei… feline…
## # … with abbreviated variable name ¹​description
```

```r
res.all <- read.csv("DE_res/Interneuron; GABAergic neuron.csv")
colnames(res.all)[1] <- "gene"
res_tableOE_tb <- res.all 


idx <- grch37$symbol%in% res_tableOE_tb$gene
ids <- grch37[idx, ]

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates <- which(duplicated(ids$symbol) == FALSE)
ids <- ids[non_duplicates, ] 

## Merge the IDs with the results 
res_ids <- inner_join(res_tableOE_tb, ids, by=c("gene"="symbol"))   


## Create background dataset for hypergeometric testing using all genes tested for significance in the results
allOE_genes <- as.character(res_ids$gene)
## Extract significant results
sigOE <- res_ids %>% dplyr::filter(padj< 0.05 & abs(log2FoldChange) > 1)
sigOE.order<- sigOE[order(sigOE$pvalue,decreasing = F),]
head(sigOE.order)
```

```
##         gene baseMean log2FoldChange     lfcSE      stat       pvalue
## 8  MTRNR2L12 330.4722      -1.862599 0.3153567 -5.906324 3.498266e-09
## 29       UBC 331.1639      -1.187903 0.2354932 -5.044320 4.551377e-07
## 21     CEND1 274.6389      -1.208004 0.2421301 -4.989070 6.067061e-07
## 7    FAM19A1 266.9323       1.383870 0.2950563  4.690190 2.729509e-06
## 13    EEF1A1 288.4591      -1.175869 0.2524100 -4.658568 3.184170e-06
## 31    DYNLL2 244.6369      -1.125122 0.2439659 -4.611800 3.991969e-06
##            padj         ensgene    entrez chr     start       end strand
## 8  2.287166e-05 ENSG00000269028 100463486   3  96335981  96337000     -1
## 29 9.916611e-04 ENSG00000150991      7316  12 125396150 125401914     -1
## 21 9.916611e-04 ENSG00000184524     51286  11    787104    790123     -1
## 7  2.974015e-03 ENSG00000183662    407738   3  68053359  68594776      1
## 13 2.974015e-03 ENSG00000156508      1915   6  74225473  74233520     -1
## 31 3.262437e-03 ENSG00000264364    140735  17  56160776  56172897      1
##           biotype
## 8  protein_coding
## 29 protein_coding
## 21 protein_coding
## 7  protein_coding
## 13 protein_coding
## 31 protein_coding
##                                                                   description
## 8                                                MT-RNR2-like 12 (pseudogene)
## 29                                                                ubiquitin C
## 21                             cell cycle exit and neuronal differentiation 1
## 7  family with sequence similarity 19 (chemokine (C-C motif)-like), member A1
## 13                         eukaryotic translation elongation factor 1 alpha 1
## 31                                            dynein, light chain, LC8-type 2
```



```r
## Run KEGG enrichment analysis 
ekegg <- enrichKEGG(gene = sigOE$entrez, 
                organism = "hsa",
                keyType = "kegg",
                pvalueCutoff = 1,
                pAdjustMethod = "fdr",
                universe = res_ids$entrez,
                minGSSize = 10,
                maxGSSize = 500,
                qvalueCutoff = 1,
                use_internal_data = FALSE)
                

## Dotplot 
dotplot(ekegg, color = "pvalue",showCategory=20)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-1.png" width="672" />



```r
## Run GO enrichment analysis 
ego <- enrichGO(gene = sigOE$gene, 
                universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "fdr", 
                pvalueCutoff = 1, 
                qvalueCutoff = 1,
                readable = TRUE)
                

## Dotplot 
dotplot(ego, color = "pvalue",showCategory=20)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-8-1.png" width="672" />
