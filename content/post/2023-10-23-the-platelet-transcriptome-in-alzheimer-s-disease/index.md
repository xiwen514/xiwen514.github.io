---
title: The platelet transcriptome in Alzheimer’s disease and aging data analysis
author: Package Build
date: '2023-10-23'
slug: []
categories: []
tags: []
---





I borrow the data from this article ( https://www.frontiersin.org/articles/10.3389/fmolb.2023.1196083/full ) to analyze Alzheimer’s disease (AD) data, aiming to show the standard analysis during my work. This report mainly analyze Bulk RNA data by comparing AD VS old to identify differential genes and explore their functions. 


This article provides DESeq2-normalized counts using Median of ratios method by 4 steps (https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html):

1. creates a pseudo-reference sample (row-wise geometric mean)
2. calculates ratio of each sample to the reference
3. calculate the normalization factor for each sample (size factor)
4. calculate the normalized count values using the normalization factor


This data has 27 samples, 3 groups(AD, old and young), each group has 9 samples. Here, I compare AD vs old only.

- AD : AD
- CO: old
- CY: young

### libraries

```r
library(readxl)
library(DESeq2)
library(apeglm)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)
library(annotables)
library(RColorBrewer)
library(ggplot2)
library(scales) 
library(viridis)
library(ComplexHeatmap)
```



### dataset pre-filterring

Before data exploration, quality control the data. While it is not strictly necessary, it is good to do some preliminary filtering of the data before running the differential expression analysis. This will reduce the size of the DEseq2DataSet object and speed up the runtime of the algorithm. The filtering criteria:

1. genes should have more than a sum total of 5 reads of support in all samples
2. number of 0 of each gene have to less than 9


```r
dat <- read_excel("Table3_The platelet transcriptome and proteome in Alzheimer_s disease and aging_ an exploratory cross-sectional study.XLSX", sheet = "DESeq2_normalized_counts") %>% data.frame()
dat1<- apply(dat[,-1],2, as.numeric)
rownames(dat1) <- dat$Ensembl_ID
dat1 <- dat1[,-c(grep("CY", colnames(dat1)))]

sumzero <- rowSums(dat1==0)
dat1 <- dat1[sumzero < 9,]
dat1 <- dat1[rowSums(dat1) > 5,] 
colnames(dat1)
```

```
##  [1] "AD2"  "AD3"  "AD4"  "AD5"  "AD6"  "AD7"  "AD8"  "AD9"  "AD10" "CO1" 
## [11] "CO2"  "CO3"  "CO4"  "CO5"  "CO9"  "CO10" "CO13" "CO14"
```

```r
dim(dat1)
```

```
## [1] 17656    18
```

```r
dat1 <- data.frame(dat1)
#write.csv(dat1,"bulk_RNA_deseq2_normalized_afterQC.csv",col.names = F)

meta <- read_excel("Table3_The platelet transcriptome and proteome in Alzheimer_s disease and aging_ an exploratory cross-sectional study.XLSX")
colnames(meta)[1] <- "sample"
meta1 <- meta[meta$sample %in% colnames(dat1),]
meta1$Condition <- as.factor(meta1$Condition)
meta1$Condition <- factor(meta1$Condition,labels=c("Old", "AD"), levels= c("Old", "AD"))
rownames(meta1) <- meta1$sample 
### Check that sample names match in both data
all(colnames(dat1) %in% rownames(meta1))
```

```
## [1] TRUE
```

```r
all(colnames(dat1) == rownames(meta1))
```

```
## [1] TRUE
```




### data exploration
#### PCA


```r
### data has been normalized, so it has non-integer values. Use round() to create a DESeq2 subject
dds <- DESeqDataSetFromMatrix(countData = round(dat1), colData = meta1, design = ~ Condition)
### normalization, here omit this step
#dds <- estimateSizeFactors(dds)
#normalized_counts <- counts(dds, normalized=TRUE)

### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
### Plot PCA 
plotPCA(rld, intgroup="Condition")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-3-1.png" width="672" />

```r
plotPCA(rld, intgroup="Sex")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-3-2.png" width="672" />

The variation on PC1 and PC2 cannot be explained by condition or sex factor.


```r
# Input is a matrix of log transformed values
rld_mat <- assay(rld)
write.csv(rld_mat,"bulk_RNA_deseq2_normalized_afterQC_rlogged.csv",col.names = F)
pca <- prcomp(t(rld_mat))
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta1, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = Condition))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-4-1.png" width="672" />

```r
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = Sex))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-4-2.png" width="672" />

The variation on PC3 and PC4 still cannot be explained by  condition or sex factor.


#### hierarchical tree


```r
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    
### Plot heatmap
pheatmap(rld_cor)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-1.png" width="672" />

```r
### PCA for samples
ggplot(df,aes(x=PC1, y=PC2, color = sample)) + geom_point()+geom_label(label=df$sample)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-2.png" width="672" />

Since the majority of genes are not differentially expressed, samples generally have high correlations with each other (values higher than 0.80). In the heatmap, CO1, AD7, CO14, AD4, CO9 have lower correlations by each other. In the PCA plot, these samples are not in the same cluster. I think they are not outliers.



### differential analysis

**The data has been normalized. If you want to use normalized count data, you need to normalize manually by this: normalized_counts <- counts(dds, normalized=TRUE).**

I would perform the differential analysis by these steps, which are integrated by DESeq function.

1. estimation of size factors: estimateSizeFactors 
2. estimation of dispersion: estimateDispersions
3. Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest


```r
dds <- DESeq(dds)
plotDispEsts(dds)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-6-1.png" width="672" />

Genes with low dispersion estimates are shrunken towards the curve, and the more accurate, higher shrunken values are output for fitting of the model and differential expression testing.

The plot shows data generally scatter around the curve, with the dispersion decreasing with increasing mean expression levels. I think this dataset is a good fit for the DESeq2 model.



#### building the results table

The shrinkage method in DESeq2 is particularly important to reduce false positives in the differential expression analysis. In the most recent versions of DESeq2, the shrinkage of LFC estimates is not performed by default. Shrinkage of the LFC estimates toward zero when the information for a gene is low, which could include:

- Low counts
- High dispersion values

**Here omit lfcShrink() step for selecting more differential genes.**


```r
contrast_oe <- c("Condition", "AD","Old")

res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.1)
res_tableOE <- lfcShrink(dds, coef = "Condition_AD_vs_Old", res=res_tableOE_unshrunken)

plotMA(res_tableOE_unshrunken)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-1.png" width="672" />

```r
plotMA(res_tableOE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-2.png" width="672" />

```r
res_tableOE_tb <- res_tableOE_unshrunken %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
```

NOTE: on p-values set to NA

- If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates, p-value and adjusted p-value will all be set to NA.
- If a row contains a sample with an extreme count outlier then the p-value and adjusted p-value will be set to NA. These outlier counts are detected by Cook’s distance.
- If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p-value will be set to NA.



```r
# Set a boolean column for significance
res_tableOE_tb$significant <- ifelse(res_tableOE_tb$pvalue < .1, "Significant", NA)

# Plot the results similar to DEseq2
ggplot(res_tableOE_tb, aes(baseMean, log2FoldChange, colour=significant)) +
  geom_point(size=1) + 
  scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() +
  geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + 
  scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-8-1.png" width="672" />

```r
# Let's add some more detail
ggplot(res_tableOE_tb, aes(baseMean, log2FoldChange, colour=pvalue)) +
  geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + 
  labs(x="mean of normalized counts", y="log fold change") +
  scale_colour_viridis(direction=-1, trans='sqrt') + 
  theme_bw() +
  geom_density_2d(colour="black", size=2)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-8-2.png" width="672" />

In general we would expect the expression of genes to remain consistent between conditions and so the MA plot should be similar to the shape of a trumpet with most points residing on a y intercept of 0. We can see from the above plots that they are in the characteristic trumpet shape of MA plots. Further we have overlayed density contours in the second plot and, as expected, these density contours are centered around a y-intercept of 0. We can further see that as the average counts increase there is more power to call a gene as differentially expressed based on the fold change.

### visualization

Here cutoff for significant genes
- pvalue< 0.05
- abs(log2FoldChange) > 1

#### heatmap

```r
sigOE <-  res_tableOE_tb %>% dplyr::filter(pvalue< 0.05 & abs(log2FoldChange) > 1)
dat2 <- rld_mat[rownames(rld_mat) %in% sigOE$gene,]
dat3 <- t(apply(dat2,1,scale))
colnames(dat3) <- colnames(dat2)
rownames(dat3) <- NULL
### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap
Heatmap(dat3,
        column_split = meta1$Condition,
        cluster_column_slices = F,
        column_names_gp = gpar(fontsize = 4),
        row_names_gp = gpar(fontsize = 6))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-9-1.png" width="672" />

286 significant genes are showed in the heatmap.

### functional analysis

```r
## Return the IDs for the gene symbols in the DE results
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
idx <- grch37$ensgene %in% res_tableOE_tb$gene
ids <- grch37[idx, ]

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates <- which(duplicated(ids$symbol) == FALSE)
ids <- ids[non_duplicates, ] 

## Merge the IDs with the results 
res_ids <- inner_join(res_tableOE_tb, ids, by=c("gene"="ensgene"))   


## Create background dataset for hypergeometric testing using all genes tested for significance in the results
allOE_genes <- as.character(res_ids$gene)
## Extract significant results
sigOE <- res_ids %>% dplyr::filter(pvalue< 0.05 & abs(log2FoldChange) > 1)
sigOE.order<- sigOE[order(sigOE$pvalue,decreasing = F),]
head(sigOE.order)
```

```
## # A tibble: 6 × 16
##   gene    baseM…¹ log2F…² lfcSE  stat  pvalue   padj signi…³ entrez symbol chr  
##   <chr>     <dbl>   <dbl> <dbl> <dbl>   <dbl>  <dbl> <chr>    <int> <chr>  <chr>
## 1 ENSG00…   19.9    -1.51 0.327 -4.62 3.81e-6 0.0440 Signif…   9283 GPR37… 1    
## 2 ENSG00…   23.1     1.39 0.305  4.57 4.98e-6 0.0440 Signif…  11067 C10or… 10   
## 3 ENSG00…   39.8    -3.12 0.715 -4.36 1.29e-5 0.0570 Signif…  85455 DISP2  15   
## 4 ENSG00…   17.2    -1.83 0.436 -4.20 2.73e-5 0.0962 Signif… 619279 ZNF704 8    
## 5 ENSG00…   23.7     2.07 0.525  3.95 7.81e-5 0.222  Signif…   5069 PAPPA  9    
## 6 ENSG00…    9.38   -3.00 0.765 -3.92 8.79e-5 0.222  Signif…  55089 SLC38… 12   
## # … with 5 more variables: start <int>, end <int>, strand <int>, biotype <chr>,
## #   description <chr>, and abbreviated variable names ¹​baseMean,
## #   ²​log2FoldChange, ³​significant
```

```r
sigOE_genes <- as.character(sigOE$gene)
```



#### KEGG terms

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

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-11-1.png" width="672" />


#### GO terms


```r
## Run GO enrichment analysis 
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "fdr", 
                pvalueCutoff = 1, 
                qvalueCutoff = 1,
                readable = TRUE)
                

## Dotplot 
dotplot(ego, color = "pvalue",showCategory=20)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-12-1.png" width="672" />




#### GSEA 

```r
## Remove any NA values
res_entrez <- dplyr::filter(res_ids, entrez != "NA")

## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrez) == F), ]

## Extract the foldchanges
foldchanges <- res_entrez$log2FoldChange

## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$entrez
## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)

head(foldchanges)
```

```
##    54996     1770     5475    79820     6549     2842 
## 3.043742 2.867605 2.703314 2.683883 2.641870 2.584722
```

```r
## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes 
              organism = "hsa", # supported organisms listed below
              nPerm = 1000, # default number permutations
              minGSSize = 20, # minimum gene set size (# genes in set) 
              pvalueCutoff = 1, # padj cutoff value
              pAdjustMethod = "fdr",
              verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result

gseaplot(gseaKEGG, geneSetID = 'hsa04810') # Regulation of actin cytoskeleton
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-13-1.png" width="672" />

```r
### gseGO
gseaGO <- gseGO(geneList = foldchanges, 
              OrgDb = org.Hs.eg.db, 
              ont = 'BP', 
              nPerm = 1000, 
              minGSSize = 20, 
              pvalueCutoff = 1,
              verbose = FALSE) 

gseaGO_results <- gseaGO@result
```

geneSetID = 'hsa04810' is Regulation of actin cytoskeleton, which is enriched in the KEGG proven by the article.

Platelet function strongly relies on a dynamic cytoskeleton, as platelet shape change is involved in granule release, platelet adhesion, and aggregation. In AD, previous platelet proteomic studies showed dysregulation of several actin cytoskeleton binding proteins.

Here, the gseaplot shows that this pathway is actually meaningful, although its p value is 0.24.




### WGCNA
Next shows WGCNA for bulk RNA data. To reduce the run time, I choose 5000 genes with highest variance among all genes.

#### Load Libraries

```r
library(WGCNA)
library(tibble)
library(dplyr)
library(rstatix)
library(corrplot)
library(ComplexHeatmap)
library(ggpubr)
library(stringr)
library(readxl)
library(ggpubr)
library(rstatix)
library(skimr)
```

#### Data preprocessing

```r
var <- apply(dat1,1,var)
var <- sort(var,decreasing = T)
datExpr <- dat1[rownames(dat1) %in% names(var)[1:5000],] %>% t()
skim(datExpr[,1:10])
```


Table: <span id="tab:unnamed-chunk-15"></span>Table 1: Data summary

|                         |                |
|:------------------------|:---------------|
|Name                     |datExpr[, 1:10] |
|Number of rows           |18              |
|Number of columns        |10              |
|_______________________  |                |
|Column type frequency:   |                |
|numeric                  |10              |
|________________________ |                |
|Group variables          |None            |


**Variable type: numeric**

|skim_variable   | n_missing| complete_rate|    mean|      sd|      p0|     p25|     p50|     p75|    p100|hist  |
|:---------------|---------:|-------------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|:-----|
|ENSG00000000460 |         0|             1|  272.27|  121.36|   94.89|  216.22|  241.52|  314.75|  491.48|▅▇▅▁▅ |
|ENSG00000000938 |         0|             1|  876.55|  200.42|  619.91|  753.09|  822.05|  963.21| 1404.22|▆▇▅▁▁ |
|ENSG00000001461 |         0|             1|  299.29|   71.72|  202.52|  246.94|  278.86|  334.22|  451.29|▇▇▃▁▃ |
|ENSG00000001629 |         0|             1|  492.39|  218.38|  188.62|  336.75|  488.04|  594.56| 1134.29|▇▇▅▁▁ |
|ENSG00000001631 |         0|             1|  254.12|   82.70|   88.58|  217.85|  258.38|  302.66|  400.88|▃▂▇▅▂ |
|ENSG00000002586 |         0|             1| 3028.06| 1840.10| 1157.23| 1763.71| 2695.92| 3950.94| 8800.71|▇▅▂▁▁ |
|ENSG00000002834 |         0|             1| 1597.10|  327.60| 1089.40| 1384.16| 1637.35| 1798.25| 2420.96|▂▅▇▁▁ |
|ENSG00000003056 |         0|             1|  541.26|   70.53|  426.13|  507.32|  522.66|  564.46|  727.09|▃▇▅▁▁ |
|ENSG00000003402 |         0|             1| 2101.26|  306.63| 1790.41| 1851.98| 1999.45| 2288.40| 2956.85|▇▂▃▁▁ |
|ENSG00000003436 |         0|             1| 1809.30|  780.81|  396.38| 1209.75| 1844.49| 2406.72| 2960.53|▅▃▇▆▆ |

```r
dim(datExpr)
```

```
## [1]   18 5000
```

```r
rownames(meta1) <- meta1$sample
all(rownames(meta1)==rownames(datExpr))
```

```
## [1] TRUE
```

#### Construction of co-expression network


```r
enableWGCNAThreads()
```

```
## Allowing parallel execution with up to 7 working processes.
```

```r
powers <- c(1:10,seq(12,20,2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType = "signed")
```

```
## pickSoftThreshold: will use block size 5000.
##  pickSoftThreshold: calculating connectivity for given powers...
##    ..working on genes 1 through 5000 of 5000
##    Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
## 1      1  0.87900  4.0800          0.920  2670.0    2730.0   3120
## 2      2  0.00432  0.0563         -0.034  1590.0    1600.0   2230
## 3      3  0.54800 -0.8330          0.481  1040.0     990.0   1740
## 4      4  0.67500 -1.0900          0.645   715.0     644.0   1420
## 5      5  0.73700 -1.1700          0.731   517.0     437.0   1190
## 6      6  0.74600 -1.2200          0.762   388.0     309.0   1020
## 7      7  0.76700 -1.2300          0.799   299.0     225.0    888
## 8      8  0.76700 -1.2600          0.814   236.0     169.0    779
## 9      9  0.78000 -1.2700          0.842   190.0     128.0    690
## 10    10  0.78400 -1.2900          0.857   155.0      99.0    616
## 11    12  0.82300 -1.3100          0.900   108.0      61.6    498
## 12    14  0.82700 -1.3400          0.913    77.6      40.3    410
## 13    16  0.84500 -1.3600          0.934    57.7      26.8    343
## 14    18  0.85000 -1.3800          0.943    44.0      18.6    290
## 15    20  0.85600 -1.4100          0.948    34.3      13.1    247
```

```r
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-16-1.png" width="672" />

```r
#Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-16-2.png" width="672" />

```r
sft$powerEstimate
```

```
## [1] 1
```




To achieve a highest R^2 and mean connectivity，assume Soft power is 7. Meanwhile, the "ward.D2" method for distance calculation shows a better performance.

#### Network construction and module detection

**softPower = 7; minModuleSize = 20; hclust : ward.D2**

- Ward: the original Ward's method essentially tried to minimize the increase in total within-cluster variance at each step using the Euclidean distance itself, not its square. 
- Ward.D2: the D2 refers to the square of the Euclidean distance. This method considers the square of the distance rather than just the distance itself.


```r
## adjacency matrix
softPower <- 7
adjacency <- adjacency(datExpr, power = softPower, type = "signed")

## TOM matrix
TOM = TOMsimilarity(adjacency, TOMType = "signed")
```

```
## ..connectivity..
## ..matrix multiplication (system BLAS)..
## ..normalization..
## ..done.
```

```r
dissTOM <- 1-TOM

## Create a dendogram using a hierarchical clustering tree
TaxaTree <- hclust(as.dist(dissTOM), method = "ward.D2")

## Module identification using dynamic tree cut
minModuleSize <- 20

dynamicMods <- cutreeDynamic(dendro = TaxaTree, distM = dissTOM, 
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
```

```
##  ..cutHeight not given, setting it to 19.2  ===>  99% of the (truncated) height range in dendro.
##  ..done.
```

```r
## Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods)
plotDendroAndColors(TaxaTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    main = "Dendrogram and module colors")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-17-1.png" width="1536" />





The hierarchical tree shows 12 modules.

Find the relationships modules and clinical traits


```r
nTaxa = ncol(datExpr)
nSamples = nrow(datExpr)

traits <- meta1
traits$group<- as.numeric(factor(traits$Condition, levels = c("AD", "Old")))
traits$Gender  <-  as.numeric(factor(traits$Sex,levels = c("F","M")))

moduleColors <- dynamicColors
MEs0 <- moduleEigengenes(datExpr, moduleColors,softPower = softPower)$eigengenes
MEs <- orderMEs(MEs0)

moduleTraitCor <- cor(MEs, traits[,5:6], use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

## visualize it
textMatrix <- paste(signif(moduleTraitCor, 2), 
                    "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = colnames(traits[5:6]), yLabels = names(MEs), 
               ySymbols = names(MEs), colorLabels = FALSE, 
               colors = blueWhiteRed(50), textMatrix = textMatrix, 
               setStdMargins = T, cex.text = 0.7, 
               main = paste("Module-trait relationships"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-18-1.png" width="1344" />






```r
df.2 <- cbind(traits["group"], MEs0)
rownames(df.2) <- traits$New.ID

palettes = c("#00aba9", "#fa6800", "#647687")

p_list <- list()
for (i in colnames(df.2)[-1]) {
  stat.test <- df.2 %>%
  pairwise_wilcox_test(as.formula(paste(i, "group", sep = " ~ ")), p.adjust.method = "fdr") %>% 
  add_xy_position(x = "group")

  p_list[[i]] <- ggboxplot(df.2, 
          x = "group", y = i,
          color = "black", 
         palette = palettes,
          add = "jitter",
          xlab = "group",
          ylab = "engine genes",
          title = i,
          width = 0.3,
          add.params = list(color = "group", shape=20, size = 1.5,jitter = 0.07)) +
  theme(legend.position='none') +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = TRUE, size = 7) +
  theme(plot.title = element_text(hjust = 0.5))
  
}

ggarrange(plotlist = p_list)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-19-1.png" width="768" />




From above heatmap and boxplot, we can see that the correlation between MEturquoise and group is significant, which p value is 0.03.

Next apply enrichment analysis for proteins in this module.


**KEGG terms**


```r
mod_conn <- intramodularConnectivity(adjacency, moduleColors, scaleByMax = FALSE)
mod_con <- add_column(mod_conn, gene=rownames(mod_conn),moduleColors=moduleColors,
                      .before = colnames(mod_conn)[1])

mod.sig <- mod_con[mod_con$moduleColors =="turquoise",]$gene
eg <- res_ids[res_ids$gene %in% mod.sig,]

uni <- res_ids[res_ids$gene %in% mod_con$gene,]

ekegg <- enrichKEGG(gene = eg$entrez, 
                organism = "hsa",
                keyType = "kegg",
                pvalueCutoff = 0.2,
                pAdjustMethod = "fdr",
                universe = uni$entrez,
                minGSSize = 10,
                maxGSSize = 500,
                qvalueCutoff = 1,
                use_internal_data = FALSE)
                

## Dotplot 
dotplot(ekegg, color = "p.adjust",showCategory=20)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-20-1.png" width="672" />

From the plot, we can see AD is enriched.

**GO terms**


```r
## Run GO enrichment analysis 
ego <- enrichGO(gene = eg$entrez, 
                universe = uni$entrez,
                keyType = "ENTREZID",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "fdr", 
                pvalueCutoff = 0.05, 
                qvalueCutoff = 0.2,
                readable = TRUE)
                

## Dotplot 
dotplot(ego, color = "p.adjust",showCategory=20)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-21-1.png" width="672" />


