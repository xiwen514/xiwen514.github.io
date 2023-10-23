---
title: The platelet proteome in Alzheimer’s disease and aging data analysis
author: Package Build
date: '2023-10-10'
slug: []
categories: []
tags: []
---





I borrow the data from this article ( https://www.frontiersin.org/articles/10.3389/fmolb.2023.1196083/full ) to analyze Alzheimer’s disease (AD) data, aiming to show the standard analysis during my work. This report mainly analyze proteomics data by comparing AD VS old to identify differential proteins and explore their functions. 


### Libraries

```r
library(utils)
library(PCAtools)
library(factoextra)
library(FactoMineR)
library(parameters)
library(tibble)
library(ggplot2)
library(ggvenn)
library(skimr)
library(dplyr)
library(ggpubr)
library(ggprism)
library(VIM)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(gridExtra)
# library(reactome.db)
# library(ReactomePA)
```

### Proteomics data introduction

This article provides two proteomics data tables(AD and old, old and young) after MAXQUANT preprocessing, but without AD and young group protein abundance table. 

- AD: AD, 10samples
- CO : old,9 samples
- CY : young, 10 samples

### input preprocessing

I cannot compare AD and young group abundance from two tables directly due to the nonrepeatability of MAXQUANT. To validate this, PCA is applied to see the patterns of same old samples from two tables.

First, quality control is a essential step to exclude data noise. Some criteria show here.

1. gene names is not NULL
2. Reverse is NULL( not "+")
3. Only.identified.by.site is NULL( not "+")
4. Potential.contaminant is NULL( not "+")
5. Q.value < 0.01
6. delete the proteins those of number of missing values more than 50% samples

For missing values, I use KNN algorithm to impute. Before imputation, transform the data using log2.

```r
data.ad <- read.delim("proteinGroups_AD_patients_2.txt") # AD and old samples,2031 proteins
## quality control
data.ad[data.ad==""] <- NA
data.ad1 <- data.ad %>%
  filter(!is.na(Gene.names)) %>%
  filter(is.na(Reverse)) %>% 
  filter(is.na(Only.identified.by.site)) %>% 
  filter(is.na(Potential.contaminant)) %>%
  filter(Q.value < 0.01) # 1810 proteins

data.ad1 <- data.ad1[,c(grep("Gene.names", colnames(data.ad1)),
                  grep("LFQ.intensity",colnames(data.ad1)))]
colnames(data.ad1)
```

```
##  [1] "Gene.names"         "LFQ.intensity.AD6"  "LFQ.intensity.AD7" 
##  [4] "LFQ.intensity.AD5"  "LFQ.intensity.AD4"  "LFQ.intensity.AD9" 
##  [7] "LFQ.intensity.AD1"  "LFQ.intensity.AD2"  "LFQ.intensity.AD3" 
## [10] "LFQ.intensity.AD8"  "LFQ.intensity.AD10" "LFQ.intensity.CO6" 
## [13] "LFQ.intensity.CO5"  "LFQ.intensity.CO7"  "LFQ.intensity.CO9" 
## [16] "LFQ.intensity.CO3"  "LFQ.intensity.CO8"  "LFQ.intensity.CO1" 
## [19] "LFQ.intensity.CO4"  "LFQ.intensity.CO2"
```

```r
## drop the proteins those of number of missing values more than 50% samples
data.ad1[data.ad1 == 0] <- NA
na_per <- rowSums(is.na(data.ad1))/(ncol(data.ad1)-1)
per_cut <- 0.5
data.ad1 <- data.ad1[na_per <= per_cut,] # 421 proteins

## log transformation and KNN impute missing values
data.ad2 <- log2(data.ad1[,-1]) %>% t()
ad.imp <- VIM::kNN(data.ad2, k=10, imp_var=F) %>% t() %>% data.frame()
colnames(ad.imp) <- colnames(data.ad1)[-1]
ad.imp$Gene.names <- data.ad1$Gene.names


#skim(ad.imp)


################################################################################
data.cy <- read.delim("proteinGroups_Drerup_2.txt") # old and young samples, 1115 proteins

## quality control
data.cy[data.cy==""] <- NA
data.cy1 <- data.cy %>%
  filter(!is.na(Gene.names)) %>%
  filter(is.na(Reverse)) %>% 
  filter(is.na(Only.identified.by.site)) %>% 
  filter(is.na(Potential.contaminant))  %>% 
  filter(Q.value < 0.01)# 911 proteins


data.cy1 <- data.cy1[,c(grep("Gene.names", colnames(data.cy1)),
                  grep("LFQ.intensity",colnames(data.cy1)))]
colnames(data.cy1)
```

```
##  [1] "Gene.names"         "LFQ.intensity.CO6"  "LFQ.intensity.CO5" 
##  [4] "LFQ.intensity.CO7"  "LFQ.intensity.CO9"  "LFQ.intensity.CO3" 
##  [7] "LFQ.intensity.CO8"  "LFQ.intensity.CO1"  "LFQ.intensity.CO4" 
## [10] "LFQ.intensity.CO2"  "LFQ.intensity.CY1"  "LFQ.intensity.CY10"
## [13] "LFQ.intensity.CY2"  "LFQ.intensity.CY3"  "LFQ.intensity.CY4" 
## [16] "LFQ.intensity.CY5"  "LFQ.intensity.CY6"  "LFQ.intensity.CY7" 
## [19] "LFQ.intensity.CY8"  "LFQ.intensity.CY9"
```

```r
## drop the proteins those of number of missing values more than 50% samples,
data.cy1[data.cy1 == 0] <- NA
na_per <- rowSums(is.na(data.cy1))/(ncol(data.cy1)-1)
per_cut <- 0.5
data.cy1 <- data.cy1[na_per <= per_cut,] # 218 proteins
## log transformation and KNN impute missing values
data.cy2 <- log2(data.cy1[,-1]+0.00001) %>% t()
cy.imp <- VIM::kNN(data.cy2, k=10, imp_var=F) %>% t() %>% data.frame()
colnames(cy.imp) <- colnames(data.cy1)[-1]
cy.imp$Gene.names <- data.cy1$Gene.names
#skim(t(cy.imp))
```

Find the shared proteins between old samples from two table, and apply PCA analysis.


```r
common.names <- intersect(ad.imp$Gene.names,cy.imp$Gene.names) # 197 common genes
ad.lfq <- ad.imp[match(common.names,ad.imp$Gene.names),grep("LFQ.intensity",colnames(ad.imp))]
rownames(ad.lfq) <- common.names

cy.lfq <- cy.imp[match(common.names,cy.imp$Gene.names),grep("LFQ.intensity",colnames(cy.imp))]
rownames(cy.lfq) <- common.names

lfq <- cbind(ad.lfq,cy.lfq)
colnames(lfq) <- gsub("LFQ.intensity.","",colnames(lfq))
## pca between old samples from two tables
old <- lfq[,grep("CO", colnames(lfq))]
old.met <- data.frame(sample= colnames(old),
                      group = c(rep("g1",9),rep("g2",9)))
rownames(old.met) <- old.met$sample

old.pca <- prcomp(old, scale = TRUE)
pca.res <- PCAtools::pca(old, metadata = old.met, removeVar = 0.1, scale = T)
biplot(pca.res,
       lab = pca.res$metadata$sample,
       colby = "group",
       hline = 0, vline = 0,
       legendPosition = 'right')
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-3-1.png" width="672" />


After QA, AD and old table remains 421 proteins; old and young table remains 218 proteins.

There are 197 common proteins in two tables and the PCA plot demonstrate a shift between the two groups, yet the intrinsic relationship among the samples within each group remains unchanged. I assume that the nonrepeatability of MAXQUANT is because of the algorithm of MAXQUANT. If I have to compare AD and young group protein abundance, I would like to use "combat" function to correct two batch table.



### differential analysis using linear regression

1. Compare AD and old(Control group) by linear regression to identify differential proteins( cutoff: p < 0.05)
2. Compare old and young(Control group) by linear regression to identify differential proteins( cutoff: p < 0.05)
3. And detect the common differential proteins.

#### AD vs old

```r
genes <- sub(";.*", "", ad.imp$Gene.names)
df1 <- ad.imp[,grep("LFQ.intensity",colnames(ad.imp))]
tdf1 <- data.frame(t(df1))
colnames(tdf1) <- make.names(genes)
#write.csv(tdf1,"protein_after_QC_log2_KNNimputation.csv")
tdf1$group <- c(rep("AD",10), rep("CO",9))
ad.genes.df <- data.frame(gene = genes, colnames=colnames(tdf1)[-ncol(tdf1)])
ad.genes.df$raw <- data.ad1$Gene.names


## PCA
ad.met <- data.frame(group=tdf1$group)
rownames(ad.met) <- rownames(tdf1)
ad.met$sample <- gsub("LFQ.intensity.", "", rownames(ad.met))
df.pca <- prcomp(df1, scale = TRUE)
pca.res <- PCAtools::pca(df1, metadata = ad.met, removeVar = 0.1, scale = T)
biplot(pca.res,
       lab = pca.res$metadata$sample,
       colby = "group",
       hline = 0, vline = 0,
       legendPosition = 'right')
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-4-1.png" width="672" />

```r
## linear regression
lm_list <- list()
for (i in colnames(tdf1)[-ncol(tdf1)]) {
  lm_list[[i]] <- lm(paste0(i, "~ group"), data = tdf1)
}
lm_fit <- lapply(lm_list, parameters)
lm_dat <- do.call(rbind,lm_fit)
lm_df <- lm_dat[lm_dat$Parameter == "groupCO",]
lm_df$fdr <- p.adjust(lm_df$p, method = "fdr")
ad.lm_df <- add_column(lm_df, features = names(lm_list), .before = 1)
ad.sig <- ad.lm_df %>%
  filter(p < 0.05)
```




#### old vs young

```r
genes1 <- sub(";.*", "", cy.imp$Gene.names)
df1 <- cy.imp[,grep("LFQ.intensity",colnames(cy.imp))]

tdf1 <- data.frame(t(df1))
colnames(tdf1) <- make.names(genes1)
tdf1$group <- c(rep("CO",9), rep("CY",10))
cy.genes.df <- data.frame(gene = genes1, colnames=colnames(tdf1)[-ncol(tdf1)])
cy.genes.df$raw <- data.cy1$Gene.names
## PCA
cy.met <- data.frame(group=tdf1$group)
rownames(cy.met) <- rownames(tdf1)
cy.met$sample <- gsub("LFQ.intensity.", "", rownames(cy.met))

df.pca <- prcomp(df1, scale = TRUE)
pca.res <- PCAtools::pca(df1, metadata = cy.met, removeVar = 0.1, scale = T)
biplot(pca.res,
       lab = pca.res$metadata$sample,
       colby = "group",
       hline = 0, vline = 0,
       legendPosition = 'right')
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-1.png" width="672" />

```r
## linear regression
lm_list <- list()
for (i in colnames(tdf1)[-ncol(tdf1)]) {
  lm_list[[i]] <- lm(paste0(i, "~ group"), data = tdf1)
}
lm_fit <- lapply(lm_list, parameters)
lm_dat <- do.call(rbind,lm_fit)
lm_df <- lm_dat[lm_dat$Parameter == "groupCY",]
lm_df$fdr <- p.adjust(lm_df$p, method = "fdr")
cy.lm_df <- add_column(lm_df, features = names(lm_list), .before = 1)
cy.sig <- cy.lm_df %>%
  filter(p < 0.05)
```



```r
sets_list <- list(AD_old = ad.sig$features, old_young = cy.sig$features)
ggvenn(sets_list)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-6-1.png" width="672" />

```r
intersect(ad.sig$features,cy.sig$features) ## common differential proteins
```

```
## [1] "HEL111"      "CFHR5"       "APOC4.APOC2" "ACTC1"
```



There are 35 differential proteins compared AD and old group, and there are 50 differential proteins compared old and young group. 4 Shared proteins are HEL111, CFHR5, APOC4-APOC2, ACTC1. The article show 4 shared platelet DEPs are STOM, MMRN1, TMSB4X, and GMPR. I think of some reasons why I have different results.

1. I do not normalize the logged data.I think logged data doesn't have to be in a specific range such as (0-1 or -1-1). Logged data has the same unit and satisfies the normal distribution for following linear regression.
2. GMPR is dropped during QA step.

### visualization

Here, I do not analyze old VS young because
1. The data quality is not good. After QA, only remains 218 proteins per sample.
2. Only 4 shared proteins between AD vs old and old VS young.

boxplot of shared proteins

I use original value of proteins for boxplot, which is more understandable easily.


```r
diff.shared <- c("HEL111", "CFHR5", "APOC4-APOC2", "ACTC1", "STOM", "MMRN1", "TMSB4X", "GMPR")
genes.df <- rbind(ad.genes.df, cy.genes.df)
shared <- unique(genes.df[genes.df$gene %in% diff.shared,]$raw)

diff.shared.ad <- data.ad1[data.ad1$Gene.names %in% shared,]
diff.shared.ad1 <- data.frame(t(diff.shared.ad[,-1]))
colnames(diff.shared.ad1) <- make.names(diff.shared.ad$Gene.names)
diff.shared.ad1$group <- c(rep("AD",10), rep("CO",9))

palette = c("#00aba9", "#fa6800", "#647687")
par(mfrow = c(4, 2))

plot_list <- list()

for (i in colnames(diff.shared.ad1)[-ncol(diff.shared.ad1)]) {
  base <- ggplot(diff.shared.ad1,
                 aes_string(x = "group", y = i)) +
    geom_violin(aes(colour = group, fill = group), trim = FALSE) +
    geom_boxplot(aes(fill = group), width = 0.2, colour = "black") +
    stat_compare_means(label = "p.signif", 
                       method = "t.test")+
    theme(legend.position = "none") 
  plot_list[[length(plot_list) + 1]] <- base
  #print(base)

}
grid.arrange(grobs = plot_list, ncol = 3)  
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-1.png" width="672" />

```r
#########################################################################
diff.shared.cy <- data.cy1[data.cy1$Gene.names %in% shared,]
diff.shared.cy1 <- data.frame(t(diff.shared.cy[,-1]))
colnames(diff.shared.cy1) <- make.names(diff.shared.cy$Gene.names)
diff.shared.cy1$group <- c(rep("CO",9), rep("CY",10))

plot_list <- list()

for (i in colnames(diff.shared.cy1)[-ncol(diff.shared.cy1)]) {
  base <- ggplot(diff.shared.cy1,
                 aes_string(x = "group", y = i)) +
    geom_violin(aes(colour = group, fill = group), trim = FALSE) +
    geom_boxplot(aes(fill = group), width = 0.2, colour = "black") +
    stat_compare_means(label = "p.signif", 
                       method = "t.test")+
    theme(legend.position = "none") 
 # print(base)
    plot_list[[length(plot_list) + 1]] <- base

}
grid.arrange(grobs = plot_list, ncol = 3)  
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-2.png" width="672" />

heatmap for 35 significant proteins (AD VS old)


```r
genes.df$unique <- make.names(genes.df$gene,unique = T)
rownames(ad.imp) <- make.names(genes,unique = T)
ad.sig.imp <- ad.imp[rownames(ad.imp) %in% ad.sig$features,-ncol(ad.imp)]
ad.sig.scale <- apply(ad.sig.imp,1,scale)
rownames(ad.sig.scale) <- colnames(ad.sig.imp)
ad.sig.scale1 <- t(ad.sig.scale)
Heatmap(ad.sig.scale1,
        column_split = ad.met$group,
        cluster_column_slices = F,
        column_names_gp = gpar(fontsize = 4),
        row_names_gp = gpar(fontsize = 6))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-8-1.png" width="672" />


The heatmap shows that the significant differences between AD and old.


volcano plot (AD VS old)


```r
ad.imp$logFC <- NULL
for (i in c(1:nrow(ad.imp))){
  ad.imp$logFC[i] <- mean(as.numeric(as.character(ad.imp[i,grep("AD", colnames(ad.imp))])))-
    mean(as.numeric(as.character(ad.imp[i,grep("CO", colnames(ad.imp))])))
}

ad.imp$features <- rownames(ad.imp)  
ad.lm_df <- left_join(ad.lm_df, ad.imp[,c("features", "logFC")])

EnhancedVolcano(ad.lm_df,
    lab = ad.lm_df$features,
    x = 'logFC',
    y = 'p',
    pCutoff = 0.05,
    FCcutoff = 0.01,
    xlim = c(-3.5, 3.5),
    ylim = c(0, 6),
    labSize = 4.0,
    legendLabSize = 8,
    legendIconSize = 5,
    legendLabels = c("NS", "Coefficient", "FDR", "FDR and Coefficient"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-9-1.png" width="672" />


The Volcano plot shows 414 proteins because of replicated proteins exist in 421 qualified proteins.

### functional analysis

To get more pathways , use the proteins those of p < 0.1.


```r
## KEGG
ad.sig <- ad.lm_df %>%
  filter(p < 0.1)


ad.sig.pro <- unique(genes.df[genes.df$colnames %in% ad.sig$features,]$gene)
## to acquire more pathways, I set p.value cutoff is 0.2.
eg = bitr(ad.sig.pro, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kegg <- enrichKEGG(eg$ENTREZID, organism = 'hsa',
                   pvalueCutoff = 0.2)
barplot(kegg, drop = F, showCategory = 12)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-10-1.png" width="672" />

```r
## GO
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

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-10-2.png" width="672" />




### WGCNA
Next shows WGCNA for proteomics data.

#### Load Library

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
#library(ggstatsplot)
library(rstatix)
#library(propr)
#library(robCompositions)
```


#### data preprocessing

```r
datExpr <- t(ad.imp[,grep("LFQ",colnames(ad.imp))])
rownames(datExpr) <- gsub("LFQ.intensity.","",rownames(datExpr))
skim(datExpr[,1:10])
```


Table: <span id="tab:unnamed-chunk-12"></span>Table 1: Data summary

|                         |                |
|:------------------------|:---------------|
|Name                     |datExpr[, 1:10] |
|Number of rows           |19              |
|Number of columns        |10              |
|_______________________  |                |
|Column type frequency:   |                |
|numeric                  |10              |
|________________________ |                |
|Group variables          |None            |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate|  mean|   sd|    p0|   p25|   p50|   p75|  p100|hist  |
|:-------------|---------:|-------------:|-----:|----:|-----:|-----:|-----:|-----:|-----:|:-----|
|ALDOC         |         0|             1| 25.94| 0.91| 23.71| 25.71| 26.01| 26.39| 28.00|▂▁▇▃▁ |
|VCL           |         0|             1| 31.16| 1.08| 29.16| 30.74| 31.42| 31.80| 33.41|▃▃▇▇▁ |
|HIST1H2AH     |         0|             1| 26.22| 1.73| 23.42| 25.10| 26.07| 27.14| 30.18|▅▇▇▂▂ |
|C9            |         0|             1| 28.74| 0.98| 26.14| 28.49| 28.95| 29.25| 30.03|▂▁▂▇▅ |
|GNB1          |         0|             1| 26.11| 0.84| 24.62| 25.65| 25.67| 26.54| 28.03|▂▇▇▁▂ |
|CAPZA1        |         0|             1| 26.78| 1.09| 24.53| 26.41| 26.75| 27.38| 29.00|▂▂▇▆▁ |
|TGFB1         |         0|             1| 24.13| 1.26| 22.37| 23.67| 23.78| 24.36| 27.80|▃▇▂▁▂ |
|NAPA          |         0|             1| 27.44| 1.08| 25.82| 26.87| 27.41| 27.76| 30.38|▃▇▂▂▁ |
|APOC1         |         0|             1| 31.08| 1.02| 28.89| 30.39| 31.11| 31.68| 32.95|▁▇▇▃▅ |
|VASP          |         0|             1| 27.33| 0.81| 25.37| 26.85| 27.30| 27.69| 29.34|▁▅▇▃▁ |

```r
traits <- read.csv("APPA-participants.csv")
traits <- traits[traits$New.ID %in% rownames(datExpr),]
traits <-traits[match(rownames(datExpr), traits$New.ID),]
traits$group <- ifelse(traits$AD=="Yes", "AD", "old")
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
## pickSoftThreshold: will use block size 421.
##  pickSoftThreshold: calculating connectivity for given powers...
##    ..working on genes 1 through 421 of 421
##    Power SFT.R.sq     slope truncated.R.sq mean.k. median.k. max.k.
## 1      1 6.20e-02  1.550000         0.8710  226.00    248.00  264.0
## 2      2 8.38e-02  1.190000         0.6710  150.00    169.00  202.0
## 3      3 3.16e-02  0.335000         0.0923  111.00    123.00  167.0
## 4      4 5.40e-02  0.222000        -0.2030   85.60     92.10  142.0
## 5      5 4.21e-02  0.131000        -0.1030   68.40     70.70  123.0
## 6      6 1.15e-06 -0.000518         0.4290   55.80     55.50  108.0
## 7      7 1.00e-01 -0.146000         0.8180   46.30     44.20   96.9
## 8      8 3.06e-01 -0.269000         0.9290   39.00     36.00   87.5
## 9      9 4.11e-01 -0.380000         0.8280   33.10     29.80   79.8
## 10    10 5.18e-01 -0.452000         0.9290   28.50     24.80   73.4
## 11    12 6.67e-01 -0.674000         0.8420   21.60     18.10   63.5
## 12    14 7.24e-01 -0.818000         0.8190   16.90     13.40   56.4
## 13    16 7.96e-01 -0.969000         0.7970   13.50     10.30   50.9
## 14    18 8.55e-01 -1.040000         0.8360   11.10      7.88   46.5
## 15    20 8.81e-01 -1.110000         0.8500    9.22      6.13   42.9
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

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-13-1.png" width="672" />

```r
#Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-13-2.png" width="672" />

```r
sft$powerEstimate
```

```
## [1] 18
```


To achieve a highest R^2 and mean connectivity，assume Soft power is 18. Meanwhile, the "complete" method for distance calculation shows a better performance.

#### Network construction and module detection

**softPower = 18; minModuleSize = 20; hclust : complete**


```r
## adjacency matrix
softPower <- 18
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
TaxaTree <- hclust(as.dist(dissTOM), method = "complete")

## Module identification using dynamic tree cut
minModuleSize <- 20

dynamicMods <- cutreeDynamic(dendro = TaxaTree, distM = dissTOM, 
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
```

```
##  ..cutHeight not given, setting it to 0.995  ===>  99% of the (truncated) height range in dendro.
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

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-14-1.png" width="1536" />


The hierarchical tree shows 10 modules.


Find the relationships modules and clinical traits


```r
nTaxa = ncol(datExpr)
nSamples = nrow(datExpr)

traits <- read.csv("APPA-participants.csv")
traits <- traits[traits$New.ID %in% rownames(datExpr),]
traits <-traits[match(rownames(datExpr), traits$New.ID),]
traits$group <- ifelse(traits$AD=="Yes", "AD", "old")
traits$group<- as.numeric(factor(traits$group, levels = c("AD", "old")))
traits$Gender  <-  as.numeric(factor(traits$Gender,levels = c("Female","Male")))

moduleColors <- dynamicColors
MEs0 <- moduleEigengenes(datExpr, moduleColors,softPower = softPower)$eigengenes
MEs <- orderMEs(MEs0)

moduleTraitCor <- cor(MEs, traits[,3:5], use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

## visualize it
textMatrix <- paste(signif(moduleTraitCor, 2), 
                    "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = colnames(traits[3:5]), yLabels = names(MEs), 
               ySymbols = names(MEs), colorLabels = FALSE, 
               colors = blueWhiteRed(50), textMatrix = textMatrix, 
               setStdMargins = T, cex.text = 0.5, 
               main = paste("Module-trait relationships"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-15-1.png" width="1344" />


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
          ylab = "engine proteins",
          title = i,
          width = 0.3,
          add.params = list(color = "group", shape=20, size = 1.5,jitter = 0.07)) +
  theme(legend.position='none') +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = TRUE, size = 7) +
  theme(plot.title = element_text(hjust = 0.5))
  
}

ggarrange(plotlist = p_list)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-16-1.png" width="672" />

In term of the correlation with group, MEred module is significant, which p value is 0.05.
Next apply enrichment analysis for proteins in MEred module.




```r
mod_conn <- intramodularConnectivity(adjacency, moduleColors, scaleByMax = FALSE)
mod_con <- add_column(mod_conn, protein=rownames(mod_conn),moduleColors=moduleColors,
                      .before = colnames(mod_conn)[1])

red.pro <- mod_con[mod_con$moduleColors =="red",]$protein
red.pro1 <- unique(genes.df[genes.df$colnames %in% red.pro,]$gene)

eg = bitr(red.pro1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kegg <- enrichKEGG(eg$ENTREZID, organism = 'hsa',
                   pvalueCutoff = 0.2)
barplot(kegg, drop = F, showCategory = 12)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-17-1.png" width="672" />

```r
## GO
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

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-17-2.png" width="672" />

The kegg pathways are African Trypanosomiasis and Malaria, which are also exist in functional analysis after linear regression. However, only 2 proteins are on the pathways.







