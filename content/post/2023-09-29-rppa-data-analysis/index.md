---
title: RPPA data analysis
author: Package Build
date: '2023-09-29'
slug: []
categories: []
tags: []
---
 
### What is RPPA?

Reverse phase protein array (RPPA) is a high-throughput antibody-based technique with the procedures similar to that of Western blots. Proteins are extracted from tumor tissue or cultured cells, denatured by SDS, printed on nitrocellulose-coated slides followed by antibody probe. 

### What are the advantages of RPPA?

- Inexpensive, high-throughput method utilizing automation for increased quality and reliability
- Sample preparation requirements are similar to that of Western blots
- Complete assay requires only 40 microliters of each sample for 150 antibodies
- Robust quantification due to serial dilution of samples

Here, I compare the differentially expressed genes based on the patients' status(aLive or dead) and their associated functions.


## libraries


```r
library(utils)
library(TCGAbiolinks)
library(tidyverse)
library(skimr)
library(factoextra)
library(FactoMineR)
library(factoextra)
library(limma)
library(tidyverse)
library(modelbased)
library(patchwork)
library(EnhancedVolcano)
library(WGCNA)
library(ggpubr)
library(rstatix)
library(circlize)
library(ComplexHeatmap)
library(clusterProfiler)
library(fgsea)
library(enrichplot)
```


## Data extraction from TCGA

```r
query.rppa <- GDCquery(
  project = "TCGA-ESCA", 
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification"
)
GDCdownload(query.rppa) 
rppa <- GDCprepare(query.rppa) 
```

## Normalization using limma

```r
data <- rppa[,-c(1:5)]
colnames(data) <- gsub("-01A","", colnames(data))
rownames(data) <- make.names(rppa$peptide_target)
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
data1 <- data.frame(t(data))
# delete the proteins which expresssion values are NA
na.count <- apply(is.na(data1),2,sum)
data2 <- data1[,na.count !=(nrow(data1))]
dim(data2)
```

```
## [1] 126 457
```

```r
skim(data2[,1:10])
```


Table: <span id="tab:unnamed-chunk-3"></span>Table 1: Data summary

|                         |              |
|:------------------------|:-------------|
|Name                     |data2[, 1:10] |
|Number of rows           |126           |
|Number of columns        |10            |
|_______________________  |              |
|Column type frequency:   |              |
|numeric                  |10            |
|________________________ |              |
|Group variables          |None          |


**Variable type: numeric**

|skim_variable  | n_missing| complete_rate|  mean|   sd|    p0|   p25|   p50|   p75| p100|hist  |
|:--------------|---------:|-------------:|-----:|----:|-----:|-----:|-----:|-----:|----:|:-----|
|X1433BETA      |         0|             1|  0.25| 0.18| -0.02|  0.14|  0.21|  0.32| 0.85|▅▇▂▁▁ |
|X1433EPSILON   |         0|             1|  0.09| 0.12| -0.24|  0.04|  0.10|  0.16| 0.64|▁▇▆▁▁ |
|X1433ZETA      |         0|             1|  0.00| 0.28| -0.58| -0.22| -0.05|  0.16| 0.72|▂▇▆▃▂ |
|X4EBP1         |         0|             1|  0.13| 0.51| -1.08| -0.18|  0.10|  0.37| 1.63|▂▇▇▂▁ |
|X4EBP1_pS65    |         0|             1| -0.13| 0.24| -0.73| -0.30| -0.14|  0.04| 0.71|▁▇▇▂▁ |
|X4EBP1_pT37T46 |         0|             1| -0.27| 0.63| -2.09| -0.68| -0.21|  0.09| 1.39|▁▅▇▆▁ |
|X4EBP1_pT70    |         0|             1|  0.12| 0.14| -0.13|  0.03|  0.11|  0.19| 0.67|▅▇▃▁▁ |
|X53BP1         |         0|             1| -0.94| 0.65| -3.00| -1.31| -0.79| -0.50| 0.36|▁▂▃▇▂ |
|ACC_pS79       |         0|             1| -0.19| 0.55| -1.51| -0.57| -0.22|  0.20| 1.10|▂▆▇▆▂ |
|ACC1           |         0|             1| -0.04| 0.60| -2.10| -0.42| -0.06|  0.31| 1.48|▁▂▇▆▂ |

```r
data.norm <- normalizeQuantiles(data2)
skim(data.norm[,1:10])
```


Table: <span id="tab:unnamed-chunk-3"></span>Table 1: Data summary

|                         |                  |
|:------------------------|:-----------------|
|Name                     |data.norm[, 1:10] |
|Number of rows           |126               |
|Number of columns        |10                |
|_______________________  |                  |
|Column type frequency:   |                  |
|numeric                  |10                |
|________________________ |                  |
|Group variables          |None              |


**Variable type: numeric**

|skim_variable  | n_missing| complete_rate|  mean|   sd|    p0|  p25|   p50|  p75| p100|hist  |
|:--------------|---------:|-------------:|-----:|----:|-----:|----:|-----:|----:|----:|:-----|
|X1433BETA      |         0|             1| -0.01| 0.45| -1.12| -0.3| -0.04| 0.24| 1.54|▂▇▇▂▁ |
|X1433EPSILON   |         0|             1| -0.01| 0.45| -1.12| -0.3| -0.04| 0.24| 1.54|▂▇▇▂▁ |
|X1433ZETA      |         0|             1| -0.01| 0.45| -1.12| -0.3| -0.04| 0.24| 1.54|▂▇▇▂▁ |
|X4EBP1         |         0|             1| -0.01| 0.45| -1.12| -0.3| -0.04| 0.24| 1.54|▂▇▇▂▁ |
|X4EBP1_pS65    |         0|             1| -0.01| 0.45| -1.12| -0.3| -0.04| 0.24| 1.54|▂▇▇▂▁ |
|X4EBP1_pT37T46 |         0|             1| -0.01| 0.45| -1.12| -0.3| -0.04| 0.24| 1.54|▂▇▇▂▁ |
|X4EBP1_pT70    |         0|             1| -0.01| 0.45| -1.12| -0.3| -0.04| 0.24| 1.54|▂▇▇▂▁ |
|X53BP1         |         0|             1| -0.01| 0.45| -1.12| -0.3| -0.04| 0.24| 1.54|▂▇▇▂▁ |
|ACC_pS79       |         0|             1| -0.01| 0.45| -1.12| -0.3| -0.04| 0.24| 1.54|▂▇▇▂▁ |
|ACC1           |         0|             1| -0.01| 0.45| -1.12| -0.3| -0.04| 0.24| 1.54|▂▇▇▂▁ |

## meta info

```r
met <- read.delim("GDCdata/clinical.tsv")
met_fil <-met[match(colnames(data),met$case_submitter_id),] 
met_fil[met_fil=="'--"] <- NA
na_count <- apply(is.na(met_fil), 2, sum)
met_filNA <- met_fil[,na_count!=126]
met1 <- met_filNA[,c("case_submitter_id", "age_at_index" ,"gender", 
               "vital_status" )]
rownames(met1) <- met1$case_submitter_id
met1$vital_status <- factor(met1$vital_status)
head(met1)
```

```
##              case_submitter_id age_at_index gender vital_status
## TCGA-LN-A49Y      TCGA-LN-A49Y           77   male        Alive
## TCGA-Q9-A6FW      TCGA-Q9-A6FW           61   male        Alive
## TCGA-LN-A49L      TCGA-LN-A49L           44   male         Dead
## TCGA-LN-A49N      TCGA-LN-A49N           50   male        Alive
## TCGA-LN-A49P      TCGA-LN-A49P           71   male        Alive
## TCGA-L5-A4OJ      TCGA-L5-A4OJ           70 female        Alive
```

## Dimentionally reduction (PCA) 

```r
data.pca <- PCA(data.norm, graph = FALSE) 
fviz_pca_ind(data.pca, label="none", habillage=met1$vital_status,
             addEllipses=TRUE, ellipse.level=0.95)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-1.png" width="672" />

## Differential analysis
#### linear regression

```r
data3 <- cbind(met1,data.norm) %>% mutate(across(gender, as.factor))
lm_list <- NULL
for (i in names(data3)[5:461]) {
  model <- lm(as.formula(paste0(i, "~ vital_status")), data = data3)
  # Recompute contrasts with a higher precision (for a smoother plot)
  contrasts <- estimate_contrasts(model, contrast ="vital_status", length = 20, standardize = TRUE)
  contrasts$pro <- i
  lm_list <- rbind(lm_list,contrasts)
}
```


```r
# filter the significant proteins
sig.res <- lm_list %>%
dplyr::mutate(asterisk_label = case_when(
  p > 0.1 ~ "ns",
  p < 0.1 & p > 0.05 ~ "*",
  p < 0.05 & p > 0.01 ~ "**",
  p < 0.01 & p > 0.001 ~ "***",
  p < 0.001 ~ "****"
)) %>%
dplyr::filter(asterisk_label != "ns")
```

#### limma

```r
design <- model.matrix(~ met1$vital_status)
limma.model <- lmFit(t(data.norm), design = design)
data.eb <- eBayes(limma.model)
limma.res <- topTable(data.eb, number = Inf)
```

```
## Removing intercept from test coefficients
```

```r
rppa$symbol <- make.names(rppa$peptide_target)
limma.res$symbol <-rownames(limma.res) 
limma.res <- left_join(limma.res,rppa[,c("peptide_target", "symbol")])
```

```
## Joining with `by = join_by(symbol)`
```

```r
sig.res.limma <- limma.res %>% filter(P.Value < 0.1)
intersect(sig.res$pro, rownames(sig.res.limma))
```

```
## character(0)
```

## Visualization 
#### boxplot for significant proteins  (p.value < 0.05)


```r
# 
# pst <- NULL
# for (i in sig.res$pro) {
#   p <- ggplot(data3,
#       aes_string(x = "vital_status", y = i)) +
#       geom_violin(aes(colour = vital_status, fill = vital_status), trim = FALSE) +
#       geom_boxplot(aes(fill = vital_status), width = 0.2, colour = "black")
#   p <- p+ ggplot2::theme_minimal() + theme(legend.position="none")
#   pst[[i]] <- p
# }
#   
```

#### Volcano plot

```r
EnhancedVolcano(limma.res,
    lab = limma.res$peptide_target,
    x = 'logFC',
    y = 'P.Value',
    pCutoff = 0.05,
    FCcutoff = 0.01,
    xlim = c(-0.5, 0.5),
    ylim = c(0, 3),
    #xlab = "difference",
    #ylab = "p.value",
    labSize = 4.0,
    #legendLabSize = 8,
  legendIconSize = 5,
    legendLabels = c("NS", "Coefficient", "FDR", "FDR and Coefficient"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-10-1.png" width="672" />
  
  
#### heatmap

```r
data_sig <- data.norm[,sig.res$pro]
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

Heatmap(data_sig,
        col= col_fun,
        row_split =met1$vital_status,
        cluster_column_slices = F,
        column_names_gp = gpar(fontsize = 6),
        row_names_gp = gpar(fontsize = 4))
```

```
## Warning: The input is a data frame-like object, convert it to a matrix.
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-11-1.png" width="960" />


## kegg

```r
## to acquire more pathways, I set p.value cutoff is 0.2.
eg = bitr(sig.res.limma$peptide_target, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
```

```
## 
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## Warning in bitr(sig.res.limma$peptide_target, fromType = "SYMBOL", toType =
## "ENTREZID", : 68.97% of input gene IDs are fail to map...
```

```r
kegg <- enrichKEGG(eg$ENTREZID, organism = 'hsa',
                   pvalueCutoff = 0.2)
```

```
## Reading KEGG annotation online: "https://rest.kegg.jp/link/hsa/pathway"...
```

```
## Reading KEGG annotation online: "https://rest.kegg.jp/list/pathway/hsa"...
```

```r
barplot(kegg, drop = F, showCategory = 12)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-12-1.png" width="672" />

```r
heatplot(kegg, foldChange=sig.res.limma[sig.res.limma$P.Value<0.2&abs(sig.res.limma$logFC)>1,]$log2fc,showCategory = 10) + ggtitle("Heatplot")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-12-2.png" width="672" />


## GO

```r
## to acquire more pathways, I set p.value cutoff is 0.2.
go <- enrichGO(
  eg$ENTREZID,
  "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.2,
  pAdjustMethod = "BH",
  #universe,
  qvalueCutoff = 0.2,
  #minGSSize = 10,
  #maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)
barplot(go, drop = F, showCategory = 12)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-13-1.png" width="672" />

```r
enrichplot::cnetplot(go,categorySize = "pvalue", 
                    foldChange=sig.res.limma[sig.res.limma$P.Value<0.2&abs(sig.res.limma$logFC)>1,]$log2fc,
                    colorEdge= TRUE)
```

```
## Warning in cnetplot.enrichResult(x, ...): Use 'color.params = list(foldChange = your_value)' instead of 'foldChange'.
##  The foldChange parameter will be removed in the next version.
```

```
## Warning in cnetplot.enrichResult(x, ...): Use 'color.params = list(edge = your_value)' instead of 'colorEdge'.
##  The colorEdge parameter will be removed in the next version.
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-13-2.png" width="672" />

```r
heatplot(go, foldChange=sig.res.limma[sig.res.limma$P.Value<0.2&abs(sig.res.limma$logFC)>1,]$log2fc,showCategory = 10) + ggtitle("Heatplot")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-13-3.png" width="672" />


## fgsea

```r
ranks <-sig.res.limma[sig.res.limma$peptide_target %in% eg$SYMBOL,]$logFC
names(ranks) <-eg$ENTREZID
ranks  <- sort(ranks,decreasing = T)

fgsea <- fgsea(pathways= kegg@geneSets,
               stats = ranks)
fgsea[fgsea$pval <0.1,]
```

```
##     pathway       pval      padj   log2err      ES       NES size leadingEdge
## 1: hsa04066 0.05172414 0.8277135 0.2765006 -0.9375 -1.716359    2   2597,7037
## 2: hsa04216 0.05172414 0.8277135 0.2765006 -0.9375 -1.716359    2   2180,7037
## 3: hsa04664 0.04619249 0.8277135 0.3217759  0.9375  1.733384    2   9846,7409
## 4: hsa04666 0.04619249 0.8277135 0.3217759  0.9375  1.733384    2   9846,7409
```

```r
#plotEnrichment(  kegg@geneSets[[head(fgsea[order(pval), ], 1)$pathway]],
#                ranks)
```

## GSEA KEGG

```r
eg = bitr(limma.res$peptide_target, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## Warning in bitr(limma.res$peptide_target, fromType = "SYMBOL", toType =
## "ENTREZID", : 66.74% of input gene IDs are fail to map...
```

```r
limma.res$SYMBOL <- limma.res$peptide_target
df <- merge(limma.res,eg,"SYMBOL")
df1 <- df[order(df$logFC,decreasing = T),]
ranks <- df1$logFC
names(ranks) <- df1$ENTREZID
head(ranks)
```

```
##      9846      2886      7409      7019       472       331 
## 0.2248719 0.1903437 0.1900618 0.1883713 0.1770109 0.1728815
```

```r
gse.KEGG <- gseKEGG(ranks, 
                    organism = "hsa", 
                    pvalueCutoff = 0.5,
                    pAdjustMethod = "BH") 
```

```
## preparing geneSet collections...
```

```
## GSEA analysis...
```

```
## leading edge analysis...
```

```
## done...
```

```r
ges_kegg_res <- gse.KEGG@result
head(gse.KEGG[order(gse.KEGG@result$p.adjust,decreasing = F),],10)
```

```
##                ID              Description setSize enrichmentScore      NES
## hsa01524 hsa01524 Platinum drug resistance      10       0.7182333 1.775986
##               pvalue  p.adjust    qvalue rank                   leading_edge
## hsa01524 0.004948142 0.1533924 0.1533924   17 tags=40%, list=11%, signal=38%
##            core_enrichment
## hsa01524 472/331/4292/7507
```

```r
dotplot(gse.KEGG, showCategory = 10, title = "Enriched Pathways", split=".sign") + facet_grid(.~.sign)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-15-1.png" width="672" />

```r
x2 = pairwise_termsim(gse.KEGG)
emapplot(x2, showCategory = 20)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-15-2.png" width="672" />

```r
ridgeplot(gse.KEGG,showCategory = 15) + labs(x = "enrichment distribution")
```

```
## Picking joint bandwidth of 0.0177
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-15-3.png" width="672" />

## GSEA GO

```r
barplot(sort(ranks, decreasing = T))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-16-1.png" width="672" />

```r
gse.GO <- gseGO(
  ranks, #geneList
  ont = "ALL",  
  OrgDb = "org.Hs.eg.db", 
  keyType = "ENTREZID",
  pvalueCutoff = 0.5,
  pAdjustMethod = "BH",
)
```

```
## preparing geneSet collections...
```

```
## GSEA analysis...
```

```
## leading edge analysis...
```

```
## done...
```

```r
gse_go_res <- gse.GO@result
head(gse.GO[order(gse.GO@result$p.adjust,decreasing = F),],10)
```

```
##            ONTOLOGY         ID
## GO:0031347       BP GO:0031347
## GO:0006325       BP GO:0006325
## GO:0050727       BP GO:0050727
## GO:0005615       CC GO:0005615
## GO:0032101       BP GO:0032101
## GO:0051649       BP GO:0051649
## GO:0009986       CC GO:0009986
## GO:0008585       BP GO:0008585
## GO:0046545       BP GO:0046545
## GO:0046660       BP GO:0046660
##                                                     Description setSize
## GO:0031347                       regulation of defense response      18
## GO:0006325                               chromatin organization      13
## GO:0050727                  regulation of inflammatory response      12
## GO:0005615                                  extracellular space      38
## GO:0032101          regulation of response to external stimulus      28
## GO:0051649                establishment of localization in cell      24
## GO:0009986                                         cell surface      12
## GO:0008585                             female gonad development      10
## GO:0046545 development of primary female sexual characteristics      10
## GO:0046660                           female sex differentiation      10
##            enrichmentScore       NES       pvalue  p.adjust    qvalue rank
## GO:0031347       0.6766236  1.941582 0.0008383577 0.2448004 0.2420942   32
## GO:0006325       0.7373285  1.937396 0.0007421146 0.2448004 0.2420942   32
## GO:0050727       0.7403734  1.899500 0.0005116069 0.2448004 0.2420942   25
## GO:0005615      -0.4741697 -1.759906 0.0043141588 0.3153585 0.3118722   51
## GO:0032101       0.5461552  1.735473 0.0057316234 0.3153585 0.3118722   47
## GO:0051649       0.5604075  1.725255 0.0062627097 0.3153585 0.3118722   14
## GO:0009986      -0.6206356 -1.721293 0.0072921276 0.3153585 0.3118722   27
## GO:0008585       0.7062552  1.707106 0.0062689631 0.3153585 0.3118722    9
## GO:0046545       0.7062552  1.707106 0.0062689631 0.3153585 0.3118722    9
## GO:0046660       0.7062552  1.707106 0.0062689631 0.3153585 0.3118722    9
##                              leading_edge
## GO:0031347 tags=56%, list=21%, signal=50%
## GO:0006325 tags=62%, list=21%, signal=53%
## GO:0050727 tags=58%, list=16%, signal=53%
## GO:0005615 tags=61%, list=34%, signal=54%
## GO:0032101 tags=54%, list=31%, signal=45%
## GO:0051649  tags=33%, list=9%, signal=36%
## GO:0009986 tags=50%, list=18%, signal=45%
## GO:0008585  tags=30%, list=6%, signal=30%
## GO:0046545  tags=30%, list=6%, signal=30%
## GO:0046660  tags=30%, list=6%, signal=30%
##                                                                                                          core_enrichment
## GO:0031347                                                            7409/472/331/6714/6647/27250/3717/2625/253943/5058
## GO:0006325                                                                      9557/29072/5080/64783/3717/2625/546/5058
## GO:0050727                                                                             472/331/6714/6647/27250/3717/2625
## GO:0005615 3486/780/4313/581/942/2194/5052/4323/26227/2539/1938/6275/952/1956/1977/983/1445/960/6648/4893/3485/7037/2597
## GO:0032101                                   7409/472/331/6714/6647/27250/3717/2625/253943/5058/5728/2113/8805/83481/673
## GO:0051649                                                                      9846/472/6714/6647/29072/4292/83548/5080
## GO:0009986                                                                                     94/952/1956/4851/960/7037
## GO:0008585                                                                                                 472/6714/6647
## GO:0046545                                                                                                 472/6714/6647
## GO:0046660                                                                                                 472/6714/6647
```

```r
dotplot(gse.GO, showCategory = 10, title = "Enriched Pathways", split=".sign") + facet_grid(.~.sign)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-16-2.png" width="672" />

```r
x2 = pairwise_termsim(gse.GO)
emapplot(x2, showCategory = 20)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-16-3.png" width="672" />

```r
ridgeplot(gse.GO,showCategory = 15) + labs(x = "enrichment distribution")
```

```
## Picking joint bandwidth of 0.017
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-16-4.png" width="672" />


## WGCNA
### Construction of co-expression network


```r
enableWGCNAThreads()
```

```
## Allowing parallel execution with up to 7 working processes.
```

```r
powers <- c(1:10,seq(12,20,2))
sft <- pickSoftThreshold(data.norm, powerVector = powers, verbose = 5,networkType = "signed")
```

```
## pickSoftThreshold: will use block size 457.
##  pickSoftThreshold: calculating connectivity for given powers...
##    ..working on genes 1 through 457 of 457
##    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
## 1      1    0.244 -67.40          0.939 228.000   228.000 232.00
## 2      2    0.509 -19.30          0.885 120.000   119.000 130.00
## 3      3    0.728  -7.94          0.858  66.000    64.700  80.40
## 4      4    0.779  -5.00          0.898  37.900    36.500  53.80
## 5      5    0.822  -3.72          0.908  22.500    21.200  38.00
## 6      6    0.873  -3.02          0.938  13.900    12.700  28.00
## 7      7    0.878  -2.61          0.908   8.850     7.790  21.30
## 8      8    0.912  -2.35          0.928   5.810     4.910  16.70
## 9      9    0.880  -2.20          0.871   3.920     3.190  13.30
## 10    10    0.896  -2.07          0.886   2.710     2.100  10.80
## 11    12    0.906  -1.92          0.881   1.400     0.966   7.69
## 12    14    0.303  -3.51          0.194   0.779     0.473   5.94
## 13    16    0.940  -1.80          0.924   0.467     0.246   4.80
## 14    18    0.929  -1.77          0.909   0.297     0.129   4.00
## 15    20    0.319  -2.79          0.202   0.199     0.071   3.42
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

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-17-1.png" width="672" />

```r
#Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-17-2.png" width="672" />

```r
sft$powerEstimate
```

```
## [1] 6
```


### Network construction and module detection

**softPower = 6; minModuleSize = 20; hclust : average**


```r
## adjacency matrix
softPower <- 6
adjacency <- adjacency(data.norm, power = softPower, type = "signed")

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
TaxaTree <- hclust(as.dist(dissTOM), method = "average")

## Module identification using dynamic tree cut
minModuleSize <- 20

dynamicMods <- cutreeDynamic(dendro = TaxaTree, distM = dissTOM, 
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
```

```
##  ..cutHeight not given, setting it to 0.977  ===>  99% of the (truncated) height range in dendro.
##  ..done.
```

```r
## Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods)

## Plot the dendrogram with module colors underneath
plotDendroAndColors(TaxaTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    main = "Dendrogram and module colors")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-18-1.png" width="1536" />

### Find the relationships modules and clinical met1

```r
nTaxa = ncol(data.norm)
nSamples = nrow(data.norm)

met1$age_at_index <- as.numeric(met1$age_at_index)
met1$gender<- as.numeric(factor(met1$gender))
met1$vital_status  <-  as.numeric(met1$vital_status)
#met1 <- met1[,-1]
#met1 <- apply(met1, 2, as.numeric)
  
moduleColors <- dynamicColors
MEs0 <- moduleEigengenes(data.norm, moduleColors,softPower = softPower)$eigengenes
MEs <- orderMEs(MEs0)

moduleTraitCor <- cor(MEs, met1[,2:4], use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

## visualize it
textMatrix <- paste(signif(moduleTraitCor, 2), 
                    "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = colnames(met1[2:4]), yLabels = names(MEs), 
               ySymbols = names(MEs), colorLabels = FALSE, 
               colors = blueWhiteRed(50), textMatrix = textMatrix, 
               setStdMargins = T, cex.text = 1, 
               main = paste("Module-trait relationships"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-19-1.png" width="1344" />


### Output hubgenes

```r
## output connectivity
mod_conn <- intramodularConnectivity(adjacency, moduleColors, scaleByMax = FALSE)
mod_con <- add_column(mod_conn, protein=rownames(mod_conn),moduleColors=moduleColors,.before = colnames(mod_conn)[1])
yellow <- mod_con[mod_con$moduleColors=="yellow",]$protein
yellow.df <- limma.res[limma.res$symbol %in% yellow,]
```


## kegg

```r
## to acquire more pathways, I set p.value cutoff is 0.2.
eg = bitr(yellow.df$peptide_target, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## Warning in bitr(yellow.df$peptide_target, fromType = "SYMBOL", toType =
## "ENTREZID", : 73.68% of input gene IDs are fail to map...
```

```r
kegg <- enrichKEGG(eg$ENTREZID, organism = 'hsa',
                   pvalueCutoff = 0.2)
barplot(kegg, drop = F, showCategory = 12)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-21-1.png" width="672" />


## GO

```r
## to acquire more pathways, I set p.value cutoff is 0.2.
go <- enrichGO(
  eg$ENTREZID,
  "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.2,
  pAdjustMethod = "BH",
  #universe,
  qvalueCutoff = 0.2,
  #minGSSize = 10,
  #maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)
barplot(go, drop = F, showCategory = 12)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-22-1.png" width="672" />




```r
## membership
datKME = signedKME(data.norm, MEs,outputColumnName="MM")
kme <- datKME


## connectivity to hub genes
top10 <- mod_con %>% group_by(moduleColors) %>% top_n(10,kWithin)
hub_dat <- add_column(mod_con, hub=0, .after = colnames(mod_con)[1])
hub_dat$hub[which(hub_dat$protein %in% top10$protein)] <- 1


cor_mat <- t(moduleTraitCor)
corp_mat <- t(moduleTraitPvalue)
module_rename <- data.frame(cbind(colnames(cor_mat), paste("M",1:ncol(cor_mat),sep = "")),stringsAsFactors = F)
module_colors <- gsub("ME","",module_rename$X1)

node_names <- hub_dat$moduleColors
for (i in 1:length(module_colors)) {
  node_names[which(hub_dat$moduleColors==module_colors[i])] <- module_rename$X2[i]
}
hub_out <- add_column(hub_dat, moduleNames = node_names, .after = colnames(hub_dat)[1])

kme1 <- kme
colnames(kme1) <- paste("M",1:ncol(cor_mat),sep = "")
hub_out_file <- cbind(hub_out,kme1)


cor_data <- cbind(data.norm,met1[2:4])
name_1 <- colnames(data.norm)
name_2 <- colnames(met1[2:4])
dat_sig <- cor_test(cor_data,vars = name_1, vars2 = name_2)
```

```
## Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
## ℹ Please use `all_of()` or `any_of()` instead.
##   # Was:
##   data %>% select(name_1)
## 
##   # Now:
##   data %>% select(all_of(name_1))
## 
## See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

```
## Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
## ℹ Please use `all_of()` or `any_of()` instead.
##   # Was:
##   data %>% select(name_2)
## 
##   # Now:
##   data %>% select(all_of(name_2))
## 
## See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

```r
sig_mat <- as_cor_mat(dat_sig)
datSIG <- sig_mat[,-1]
rownames(datSIG) <- sig_mat$rowname
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
sig <- datSIG
```
