---
title: Intergration of the platelet transcriptome and proteome in Alzheimer’s disease
  and aging data analysis
author: Package Build
date: '2023-10-23'
slug: []
categories: []
tags: []
---






I borrow the data from this article ( https://www.frontiersin.org/articles/10.3389/fmolb.2023.1196083/full ) to analyze Alzheimer’s disease (AD) data, aiming to show the standard analysis during my work. This report mainly integrate transccriptomes and proteomics data to identify a highly correlated multi-omics signature that discriminates AD and old group.



### Load libraries

```r
library(annotables)
library(mixOmics)
library(dplyr)
library(DescTools)
library(readxl)
library(tibble)
library(stringr)
library(BiocParallel)
library(tidyselect)
```


### Data preparation

Regarding the preprocessing, all of methods in mixOmics assume the right preprocessing method has been applied.

Here, for transcriptome data, I use rlog data by DESeq2.For proteomic data, I use log2 data. Additionally, 200 features with highest variance in each dataset are selected for saving runtime.


```r
met1 <- read.csv("APPA-participants.csv")
met2 <- read_excel("Table3_The platelet transcriptome and proteome in Alzheimer_s disease and aging_ an exploratory cross-sectional study.XLSX")
colnames(met2)[1] <- "New.ID"
met <- inner_join(met1,met2)
met <- met[-grep('CY', met$New.ID),]

rna <- read.csv("bulk_RNA_deseq2_normalized_afterQC_rlogged.csv")
rna1 <- t(rna[,-1]) 
colnames(rna1) <- rna$X
rna1 <- rna1[rownames(rna1) %in% met$New.ID,]
rna1 <- rna1[match(met$New.ID, rownames(rna1)),]
```





```r
## select 100 genes with the highest variance among all genes
rna.var <- apply(rna1,2,var)
rna.var <- sort(rna.var,decreasing = T)
rna2 <- rna1[,colnames(rna1) %in% names(rna.var)[1:500]]


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
idx <- grch37$ensgene %in% colnames(rna2)
ids <- grch37[idx, ]

## The gene names can map to more than one Ensembl ID (some genes change ID over time),
## so we need to remove duplicate IDs 
non_duplicates <- which(duplicated(ids$symbol) == FALSE)
ids <- ids[non_duplicates, ]

rna2 <- rna2[,colnames(rna2) %in% ids$ensgene]
all(colnames(rna2)== ids$ensgene)
```

```
## [1] TRUE
```

```r
colnames(rna2) <- ids$symbol
dim(rna2)
```

```
## [1]  15 476
```

```r
pro <- read.csv("protein_after_QC_log2_KNNimputation.csv")
rownames(pro) <- gsub("LFQ.intensity.","", pro$X)
pro1 <- pro[rownames(pro) %in% met$New.ID, -1]
pro1 <- pro1[match(met$New.ID,rownames(pro1)),]

## select 100 proteins with the highest variance among all proteins
pro.var <- apply(pro1,2,var)
pro.var <- sort(pro.var,decreasing = T)
pro2 <- pro1[,colnames(pro1) %in% names(pro.var)[1:500]]
dim(pro2)
```

```
## [1]  15 421
```

```r
all(rownames(rna2)==rownames(pro2))
```

```
## [1] TRUE
```



```r
## data design
data = list(transcriptome = rna2,
            proteome = pro2)
Y = factor(met$Condition)

## design matrix
design = matrix(0.1, ncol = length(data), nrow = length(data), 
              dimnames = list(names(data), names(data)))
diag(design) = 0

## Tuning the number of components
sgccda.res = block.splsda(X = data, Y = Y, ncomp = 5, 
                         design = design)
set.seed(123) 
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 5)
plot(perf.diablo)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-4-1.png" width="672" />

From above plot, 2 PC has the lowerst error rate of mahalanobis.dist calculation. But I choose centroids distance calculation for more understandable comprehension and faster calculation. In addition, the data does not have large samples and features so that I set folds is 5 for validation.



```r
# set the optimal ncomp value
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
# show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote 
```

```
##             max.dist centroids.dist mahalanobis.dist
## Overall.ER         1              1                1
## Overall.BER        1              1                1
```

```r
## Tuning keepX
set.seed(123)
test.keepX = list(transcriptome = c(seq(10,50,5)),
                proteome = c(seq(10,50,5)))

Y
```

```
##  [1] Old Old Old Old Old Old AD  AD  AD  AD  AD  AD  AD  AD  AD 
## Levels: AD Old
```

```r
tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = 2,
                             test.keepX = test.keepX, design = design,
                             validation = 'Mfold', folds = 5, nrepeat = 1,
                            dist = "centroids.dist")
## list.keepX
list.keepX = tune.TCGA$choice.keepX # set the optimal values of features to retain
list.keepX$transcriptome
```

```
## [1] 20 10
```

```r
list.keepX$proteome
```

```
## [1] 10 30
```

```r
## Final model
sgccda.res = block.splsda(X = data, Y = Y, ncomp = 2, 
                        keepX = list.keepX, design = design)
```

We can see from above results that the optimal features number is 10 in PC1 and PC2 in each omics data.



```r
# the features selected to form the first component
selectVar(sgccda.res, block = 'transcriptome', comp = 1)$transcriptome$value
```

```
##                 value.var
## HMOX2          0.46792977
## CSNK2B         0.40781994
## PBX2          -0.31405520
## SLC39A7       -0.30168835
## RP11-958N24.1 -0.28663687
## ABHD16A        0.27505066
## CXCL5          0.27310683
## THRB           0.21602778
## CATSPERB       0.20611244
## GMPR2          0.20277692
## DISP2         -0.17867785
## REEP1         -0.10297366
## PDGFRA         0.09515494
## SIK1           0.06325613
## RP11-203J24.9 -0.05524005
## FAM60A         0.04428573
## LY6E          -0.02420327
## PKIB           0.02372416
## NEXN           0.01968624
## LONRF1        -0.01656365
```

```r
sel.rna <- selectVar(sgccda.res, block = 'transcriptome', comp = 1)$transcriptome$name
selectVar(sgccda.res, block = 'proteome', comp = 1)$proteome$value
```

```
##                  value.var
## S100A6         -0.60610829
## A30             0.41131318
## KNG1            0.38445985
## RBP4            0.32364220
## CFI.1           0.30044365
## HBA2           -0.26009559
## DKFZp686C11235  0.16921710
## ALAD           -0.09620698
## HBB            -0.09251638
## HEL.S.2a       -0.08083907
```

```r
sel.pro <- selectVar(sgccda.res, block = 'proteome', comp = 1)$proteome$name
```
Here shows the selected features on PC1 from two data sets.



### Plot results

```r
plotDiablo(sgccda.res, ncomp = 1)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-1.png" width="672" />

```r
plotDiablo(sgccda.res, ncomp = 2)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-2.png" width="672" />

As we can see, the first components from two data set are highly correlated to each other.The colours and ellipses related to the sample subtypes indicate the discriminative power of each component to separate AD and old group.  


```r
plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-8-1.png" width="672" />

Clustering of the samples can be better assessed with this plot. It seems that the quality of clustering of the RNA data and proteomics  is similar. This suggests that two data sets ikely to hold similar discriminative power within the model.


```r
plotVar(sgccda.res, var.names =T, style = 'graphics', cutoff=0.5, legend = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-9-1.png" width="864" />

These first two components correlate highly with the selected variables from the proteomics and transcriptomes dataset.


```r
circosPlot(sgccda.res,  cutoff = 0.5, line = F,size.variables= 0.6, size.labels = 1,size.legend = 1)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-10-1.png" width="672" />

We can see that each RNA and proteins are correlated whether they are negative or postiive correlated.



```r
plotLoadings(sgccda.res, 1,comp = 1, contrib = 'max', method = 'median')
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-11-1.png" width="672" />

```r
plotLoadings(sgccda.res, 2,comp = 1, contrib = 'max', method = 'median')
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-11-2.png" width="672" />



```r
par(mar = c(5, 4, 2, 2))
cimDiablo(sgccda.res,legend.position = "topright",size.legend = 0.7,cluster = "column",comp=1,margins = c(5, 5))
```

```
## 
## trimming values to [-3, 3] range for cim visualisation. See 'trim' arg in ?cimDiablo
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-12-1.png" width="960" />

```r
cimDiablo(sgccda.res,legend.position = "topright",size.legend = 0.7,cluster = "column", comp=c(1,2),margins = c(5, 5))
```

```
## 
## trimming values to [-3, 3] range for cim visualisation. See 'trim' arg in ?cimDiablo
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-12-2.png" width="960" />


### Performance of the model

```r
# run repeated CV performance evaluation
perf.diablo = perf(sgccda.res, validation = 'Mfold', 
                   M = 10, nrepeat = 10, 
                   dist = 'centroids.dist') 

perf.diablo$MajorityVote.error.rate
```

```
## $centroids.dist
##                 comp1     comp2
## AD          0.6555556 0.6555556
## Old         0.7500000 0.8000000
## Overall.ER  0.6933333 0.7133333
## Overall.BER 0.7027778 0.7277778
```

```r
perf.diablo$WeightedVote.error.rate
```

```
## $centroids.dist
##                 comp1     comp2
## AD          0.4888889 0.5444444
## Old         0.5666667 0.6166667
## Overall.ER  0.5200000 0.5733333
## Overall.BER 0.5277778 0.5805556
```



```r
auc.splsda = auroc(sgccda.res, roc.block = "transcriptome", 
                   roc.comp = 2, print = FALSE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-14-1.png" width="672" />

```r
auc.splsda = auroc(sgccda.res, roc.block = "proteome", 
                   roc.comp = 2, print = FALSE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-14-2.png" width="672" />







