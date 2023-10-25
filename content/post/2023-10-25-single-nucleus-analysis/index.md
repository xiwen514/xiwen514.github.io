---
title: Single nucleus analysis
author: Package Build
date: '2023-10-25'
slug: []
categories: []
tags: []
---







I borrow the data from this article (https://www.pnas.org/doi/10.1073/pnas.2008762117#data-availability) to analyze Alzheimer’s disease (AD) data, aiming to show the standard analysis during my work. This report mainly analyze single nucleus RNA seq data of prefrontal cortex from 12 AD and 9 NC to identify differential signatures and explore these signatures' functions.


The data is prepossessed by Cell Ranger (version 3.0.1) with the default settings.

### Libraries

```r
library(Seurat)
library(ggplot2)
library(gridExtra)
library(ggprism)
library(SeuratObject)
library(future)
library(future.apply)
library(patchwork)
#library(scMayoMap)
```



### GEO Data preprocessing

Create directories for each sample and each directory should have three files: show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz.


```r
# For output from CellRanger >= 3.0 with multiple data types
files <- list.files() # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz


sample_names <- unique(sub("_(barcodes|features|matrix).+$", "", files))
sample_names <- sample_names[-22]
# create directories for each sample
lapply(sample_names, dir.create)

# move corresponding files to each dir and rename each file.
for (file in files) {
  sample_name <- sub("_(barcodes|features|matrix).+$", "", file)
  new.file <- sub(paste0(sample_name,"_"), "",file)
  file.rename(file, file.path(sample_name, new.file))
}
```


There are 21 samples totally, here I use the first three samples per group.

#### Initialize the Seurat object with the raw (non-normalized data).


```r
# Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
files <- list.files() 
files <- c("GSM4775561_AD1", "GSM4775562_AD2", "GSM4775563_AD4", 
           "GSM4775573_NC3", "GSM4775574_NC7", "GSM4775575_NC11" )

# filenames extraction
sample_names <- unique(sub("_(barcodes|features|matrix).+$", "", files))

sce <- list()
for (sample in sample_names){
  s <- Read10X(sample)
  sce[sample] = CreateSeuratObject(counts =s, min.cells = 3, min.features = 200)
}


pbmc  = merge(x = sce[[1]], y = sce[2:6], add.cell.ids = names(sce), merge.data = TRUE)

pbmc@meta.data$group_id <- NA
pbmc@meta.data$group_id[grep("AD", rownames(pbmc@meta.data))] <- "AD"
pbmc@meta.data$group_id[grep("NC", rownames(pbmc@meta.data))] <- "NC"
pbmc@meta.data$orig.ident[grep("AD", rownames(pbmc@meta.data))] <- gsub(".*_(AD\\d+)_.*", "\\1", rownames(pbmc@meta.data)[grep("AD", rownames(pbmc@meta.data))])
pbmc@meta.data$orig.ident[grep("NC", rownames(pbmc@meta.data))] <- gsub(".*_(NC\\d+)_.*", "\\1", rownames(pbmc@meta.data)[grep("NC", rownames(pbmc@meta.data))])
#saveRDS(pbmc,file= "../GSE157827_seurat.RDS")
```

### Quality control



```r
pbmc <- readRDS("GSE157827_seurat.RDS")
pbmc$percent.mt <- PercentageFeatureSet(pbmc, pattern = "^MT-")
## nFeature_RNA
qc.p1 <- VlnPlot(pbmc, 
        features = "nFeature_RNA", 
        ncol = 3,
        pt.size = 0,
        group.by = "orig.ident",
        raster=FALSE) +
  scale_fill_prism(palette = "floral") +
  theme(legend.position="none",
        axis.title.x = element_blank())

qc.p2 <- VlnPlot(pbmc, 
        features = "nFeature_RNA", 
        ncol = 3,
        pt.size = 0,
        group.by = "group_id") +
  scale_fill_prism(palette = "prism_dark") +
  theme(legend.position="none",
        axis.title.x = element_blank())

grid.arrange(qc.p1, qc.p2, ncol = 2)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-4-1.png" width="672" />





```r
## nCount_RNA
qc.p3 <- VlnPlot(pbmc, 
        features = "nCount_RNA", 
        ncol = 3,
        pt.size = 0,
        group.by = "orig.ident",
        raster=FALSE) +
  scale_fill_prism(palette = "floral") +
  theme(legend.position="none",
        axis.title.x = element_blank())

qc.p4 <- VlnPlot(pbmc, 
        features = "nCount_RNA", 
        ncol = 3,
        pt.size = 0,
        group.by = "group_id") +
  scale_fill_prism(palette = "prism_dark") +
  theme(legend.position="none",
        axis.title.x = element_blank())

grid.arrange(qc.p3, qc.p4, ncol = 2)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-1.png" width="672" />



```r
## percent_mt
qc.p5 <- VlnPlot(pbmc, 
        features = "percent.mt", 
        ncol = 3,
        pt.size = 0,
        group.by = "orig.ident",
        raster=FALSE) +
  scale_fill_prism(palette = "floral") +
  theme(legend.position="none",
        axis.title.x = element_blank())

qc.p6 <- VlnPlot(pbmc, 
        features = "percent.mt", 
        ncol = 3,
        pt.size = 0,
        group.by = "group_id") +
  scale_fill_prism(palette = "prism_dark") +
  theme(legend.position="none",
        axis.title.x = element_blank())

grid.arrange(qc.p5, qc.p6, ncol = 2)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-6-1.png" width="672" />


From above plots, we can see that there are around 10000 features, 10000 counts in each cell. Cells of NC11 sample have  higher mitochondrial percents, which means the sample quality is lower than others.


```r
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-1.png" width="768" />

```r
dim(pbmc)
```

```
## [1] 27231 44137
```





Select cells for further analysis with three criteria:

1. filter cells that more than three median absolute deviations (MADs) of the log10 read counts below the median.
2. filter cells that more than three MADs of the log10 genes detected below the median.
3. filtered out nuclei that have more than 2% mitochondrial unique molecular identifiers (UMI). Here the median is 0.8 and 3*MAD mito percent is 2.18.



```r
min.nfeatures.thr <- median(log10(pbmc$nFeature_RNA)) - 3*mad(log10(pbmc$nFeature_RNA))
max.nfeatures.thr <- median(log10(pbmc$nFeature_RNA)) + 3*mad(log10(pbmc$nFeature_RNA))

min.ncount.thr <- median(log10(pbmc$nCount_RNA)) - 3*mad(log10(pbmc$nCount_RNA))
max.ncount.thr <- median(log10(pbmc$nCount_RNA)) + 3*mad(log10(pbmc$nCount_RNA))

nfeatures.passed <- log10(pbmc$nFeature_RNA) > min.nfeatures.thr & 
  log10(pbmc$nFeature_RNA) < max.nfeatures.thr
ncount.passed <- log10(pbmc$nCount_RNA) > min.ncount.thr & 
  log10(pbmc$nCount_RNA) < max.ncount.thr
mito.passed <- pbmc$percent.mt < 2

pbmc.sel <- pbmc[,which(nfeatures.passed & ncount.passed & mito.passed)]
dim(pbmc.sel)
```

```
## [1] 27231 32376
```

```r
## nFeature_RNA
qc.p1 <- VlnPlot(pbmc.sel, 
        features = "nFeature_RNA", 
        ncol = 3,
        pt.size = 0,
        group.by = "orig.ident",
        raster=FALSE) +
  scale_fill_prism(palette = "floral") +
  theme(legend.position="none",
        axis.title.x = element_blank())

qc.p2 <- VlnPlot(pbmc.sel, 
        features = "nFeature_RNA", 
        ncol = 3,
        pt.size = 0,
        group.by = "group_id") +
  scale_fill_prism(palette = "prism_dark") +
  theme(legend.position="none",
        axis.title.x = element_blank())

grid.arrange(qc.p1, qc.p2, ncol = 2)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-8-1.png" width="672" />



```r
## nCount_RNA
qc.p3 <- VlnPlot(pbmc.sel, 
        features = "nCount_RNA", 
        ncol = 3,
        pt.size = 0,
        group.by = "orig.ident",
        raster=FALSE) +
  scale_fill_prism(palette = "floral") +
  theme(legend.position="none",
        axis.title.x = element_blank())

qc.p4 <- VlnPlot(pbmc.sel, 
        features = "nCount_RNA", 
        ncol = 3,
        pt.size = 0,
        group.by = "group_id") +
  scale_fill_prism(palette = "prism_dark") +
  theme(legend.position="none",
        axis.title.x = element_blank())

grid.arrange(qc.p3, qc.p4, ncol = 2)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-9-1.png" width="672" />



```r
## percent.mt
qc.p5 <- VlnPlot(pbmc.sel, 
        features = "percent.mt", 
        ncol = 3,
        pt.size = 0,
        group.by = "orig.ident",
        raster=FALSE) +
  scale_fill_prism(palette = "floral") +
  theme(legend.position="none",
        axis.title.x = element_blank())

qc.p6 <- VlnPlot(pbmc.sel, 
        features = "percent.mt", 
        ncol = 3,
        pt.size = 0,
        group.by = "group_id") +
  scale_fill_prism(palette = "prism_dark") +
  theme(legend.position="none",
        axis.title.x = element_blank())

grid.arrange(qc.p5, qc.p6, ncol = 2)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-10-1.png" width="672" />


From above plots, NC11 sample is not perform well. Only reamins less than 5000 fearures.

### Dimensionally reduction

Here, Normalization the data use Log. 


```r
pbmc.sel <- NormalizeData(pbmc.sel, normalization.method = "LogNormalize", scale.factor = 10000)

### Identification of highly variable features (modeling the mean-variance relationship)
pbmc.sel <- FindVariableFeatures(pbmc.sel, selection.method = "vst", nfeatures = 2000)
#### Identify the 10 most highly variable genes
```


```r
pbmc.sel <- readRDS("Seu_PCA_UMAP_annotated.RDS")
top10 <- head(VariableFeatures(pbmc.sel), 10)
#### plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc.sel)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-12-1.png" width="672" />

```r
plot2
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-12-2.png" width="672" />




```r
### Scaling the data
all.genes <- rownames(pbmc.sel)
#plan("multiprocess", workers = 10)
pbmc.sel <- ScaleData(pbmc.sel, features = all.genes)
pbmc.sel <- RunPCA(pbmc.sel, features = VariableFeatures(object = pbmc.sel), npcs = 50)
```


```r
all.genes <- rownames(pbmc.sel)
VizDimLoadings(pbmc.sel, dims = 1:5, reduction = "pca")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-14-1.png" width="768" />


```r
pca.p1 <- DimPlot(object = pbmc.sel, 
        reduction = "pca", 
        group.by = "orig.ident",
        pt.size = 1,
        cols = ggprism_data$colour_palettes$floral2[1:12]) 

pca.p2 <- DimPlot(object = pbmc.sel, 
        reduction = "pca", 
        group.by = "group_id",
        pt.size = 1,
        cols = ggprism_data$colour_palettes$prism_dark2[1:4]) 


pca.p1 | pca.p2 
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-15-1.png" width="960" />


```r
DimHeatmap(pbmc.sel, dims = 1:15, cells = 1000, balanced = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-16-1.png" width="768" />



```r
pbmc.sel <- FindNeighbors(pbmc.sel, dims = 1:30)

pbmc.sel <- FindClusters(pbmc.sel, resolution = 0.04)
pbmc.sel <- RunUMAP(pbmc.sel, dims = 1:30)
#saveRDS(pbmc.sel, file = "Seu_PCA_UMAP.RDS")
```


```r
umap.p1 <- DimPlot(pbmc.sel, reduction = "umap", pt.size = 1)

umap.p2 <- DimPlot(pbmc.sel, reduction = "umap",
        group.by = "group_id", pt.size = 1)

umap.p3 <- DimPlot(pbmc.sel, reduction = "umap",group.by = "group_id",
        split.by = "group_id",pt.size = 1)

(umap.p1 | umap.p2) /
      umap.p3
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-18-1.png" width="672" />

Show in New Window
Error in DimPlot(pbmc.mer, reduction = "umap", group.by = "orig.ident",  : 
  could not find function "DimPlot"

R version 4.3.1 (2023-06-16 ucrt) -- "Beagle Scouts"
Copyright (C) 2023 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64 (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with

The UMAP shows 10 clusters, which all present in AD and NC group.


### Cell type annotation

```r
# pbmc.sel <- readRDS("Seu_PCA_UMAP.RDS")
# seurat.markers <- FindAllMarkers(pbmc.sel,  method = 'MAST')

seurat.markers <- readRDS("seurat_markers.RDS")
scMayoMap.obj <- scMayoMap(data = seurat.markers, 
                           database=scMayoMapDatabase, 
                           tissue = 'brain')
plt <- scMayoMap.plot(scMayoMap.object = scMayoMap.obj)


res <- scMayoMap.obj$res
pbmc.mer <- pbmc.sel
pbmc.mer$cell_type <- ""

clus <- unique(pbmc.sel$seurat_clusters)
for (i in 1:length(clus)) {
  pbmc.mer$cell_type[pbmc.mer$seurat_clusters == clus[i]] <- res$celltype[res$cluster == clus[i]]
}



new.cluster.ids <- res$celltype
names(new.cluster.ids) <- levels(pbmc.mer)
pbmc.mer <- RenameIdents(pbmc.mer, new.cluster.ids)
```


```r
pbmc.mer <- pbmc.sel
DimPlot(pbmc.mer, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-20-1.png" width="672" />



```r
DimPlot(pbmc.mer, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident", label = F, pt.size = 0.5)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-21-1.png" width="960" />


```r
DimPlot(pbmc.mer, reduction = "umap", group.by = "group_id", split.by ="group_id",   label = F, pt.size = 0.5)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-22-1.png" width="768" />

```r
#saveRDS(pbmc.mer, file = "Seu_PCA_UMAP_annotated.RDS")
```

From above UMAP, all clusters are annotated by scMayoMapDatabse.

