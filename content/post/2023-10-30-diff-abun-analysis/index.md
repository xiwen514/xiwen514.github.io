---
title: Differential abundance analysis 
author: Package Build
date: '2023-10-30'
slug: []
categories: []
tags: []
---




I borrow the data from this article (https://www.pnas.org/doi/10.1073/pnas.2008762117#data-availability) to analyze Alzheimer’s disease (AD) data, aiming to show the standard analysis during my work. The previous blog mainly focused on the quality control, dimensionality reduction, clustering, and cell type annotation of single-nuclei RNA seq data. This blog further analyzes the data for differential abundance analysis.


### Libraries

```r
library(Seurat)
library(dplyr)
library(ggplot2)
library(edgeR)
```



### Abundance visualization

```r
seu <- readRDS("../../../../../transcriptome_protein_AD_plat/sc_RNA/Seu_PCA_UMAP_annotated.RDS")
Seurat::DimPlot(seu, reduction="umap",label = TRUE)+NoLegend()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-1.png" width="672" />

```r
Seurat::DimPlot(seu, reduction="umap", group.by="group_id")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-2.png" width="672" />


```r
# Create single cell experiment object
sce <- as.SingleCellExperiment(seu, assays="RNA")
```

```
## Warning: The following arguments are not used: assays
```

```r
n_cells <- table(sce$cell_type, sce$orig.ident) 
freqs <- prop.table(n_cells, margin = 2)

coldata <- data.frame(sce@colData)
coldata <- coldata[!duplicated(coldata$orig.ident),]
colnames(coldata)[1] <- "Sample"

df <- as.data.frame(t(freqs))
colnames(df)[1:3] <- c("Sample", "cluster_id", "frequency")
df <- left_join(df, coldata[,c("Sample", "group_id")])
```

```
## Joining with `by = join_by(Sample)`
```

```r
colnames(df)[4] <- "diagnosis_id"

df$diagnosis_id <- factor(df$diagnosis_id, levels = c("AD","NC"))

ggplot(df, aes(x = Sample, y = frequency, fill = cluster_id)) +
    geom_bar(stat = "identity", position = "fill") + 
    facet_wrap(~ diagnosis_id, scales = "free_x", ncol = 4) +
    theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-3-1.png" width="672" />

```r
ggplot(df, aes(x = diagnosis_id, y = frequency, fill = cluster_id)) +
    geom_bar(stat = "identity", position = "fill") + 
    theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-3-2.png" width="672" />

```r
par(mfrow = c(4, 2))

ggplot(df, aes(x = diagnosis_id, y = frequency, color = diagnosis_id)) +
    geom_boxplot(outlier.colour = NA) +  geom_jitter() +
    facet_wrap(~ cluster_id, scales = "free_y", ncol = 3) +
    theme_classic() +
    theme(legend.position = "none")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-3-3.png" width="672" />


### Use EdgeR to perform differential abundance analysis


```r
abundances <- unclass(n_cells) 
abundances <-  abundances[,-4] # delete NC11, assume NC11 as a outlier
#Attaching some column metadata.
extra.info <- sce@colData[match(colnames(abundances), sce$orig.ident),]

# Creates a DGEList object from a table of counts
y.ab <- DGEList(abundances, samples=extra.info)
summary(y.ab$samples$lib.size)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    3244    4237    4640    6438    5503   14568
```

```r
# Filter out low-abundance labels
keep <- filterByExpr(y.ab, group=y.ab$samples$seurat_clusters)
y.ab <- y.ab[keep,]
```

For a DA analysis of cluster abundances, filtering is generally not required as most clusters will not be of low-abundance (otherwise there would not have been enough evidence to define the cluster in the first place).



```r
sample <- as.factor(y.ab$samples$orig.ident)
plotMDS(y.ab, pch=16, col=c(2:8)[sample], main="MDS")
 legend("topright", legend= levels(sample),pch=16, col=2:8, cex=0.8)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-1.png" width="672" />

```r
# Formulate the design matrix
design <- model.matrix(~factor(group), y.ab$samples)


# Estimate the NB dispersion for each cluster
y.ab <- estimateDisp(y.ab, design, trend="none")
y.ab$common.dispersion
```

```
## [1] 0.3259443
```

```r
plotBCV(y.ab)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-2.png" width="672" />

Biological coefficient of variation (BCV) for each label with respect to its average abundance. BCVs are defined as the square root of the NB dispersion. Common dispersion estimates are shown in red.


```r
# QL dispersion
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
plotQLDisp(fit.ab)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-6-1.png" width="672" />

QL dispersion estimates for each label with respect to its average abundance. Quarter-root values of the raw estimates are shown in black while the shrunken estimates are shown in red. Shrinkage is performed towards the common dispersion in blue.



```r
# Test for differences in abundance
res <- glmQLFTest(fit.ab, coef=ncol(design))
res_df <- topTags(res, n = Inf)$table %>% round(2)
res_df
```

```
##                                logFC logCPM    F PValue  FDR
## Interneuron; GABAergic neuron   0.75  16.39 4.32   0.10 0.74
## Oligodendrocyte                -0.62  18.07 1.39   0.29 0.74
## Endothelial cell               -1.65  13.89 1.26   0.32 0.74
## Oligodendrocyte precursor cell  0.63  16.40 0.99   0.37 0.74
## Glutamatergic neuron            0.16  18.14 0.35   0.58 0.93
## Macrophage; Microglia cell      0.18  16.23 0.10   0.76 0.98
## Bergmann glial cell; Astrocyte -0.04  16.78 0.00   0.95 0.98
## GABAergic neuron                0.02  15.87 0.00   0.98 0.98
```

```r
summary(decideTests(res))
```

```
##        factor(group)NC
## Down                 0
## NotSig               8
## Up                   0
```
None of clusters are significant between AD and NC.





