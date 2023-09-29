---
title: Introduction to PCA
author: Package Build
date: '2023-09-14'
slug: []
categories: []
tags: []
---


PCA is an unsupervised analysis that reflects the true trend of the data itself. 

PCA is one of the dimension reduction methods.The purpose of PCA (Principal Component Analysis) is to transform the original data into a new coordinate system, where the variance of the data  along the axes is maximized. The first principal component (PC1) represents the direction of the greatest variance in the data, while the second principal component (PC2) is orthogonal to PC1 and represents the direction of the second greatest variance, and so on.

When you plot the results of PCA on a two-dimensional plane，you are essentially looking at the projection of the data onto these two main directions. The distance between points reflects their relative positions in these two directions, and indeed, it is the Euclidean distance.

However, remember that PCA involves a linear transformation of the data, implying that distances in the original data space may have been distorted. Thus, the Euclidean distances between points on a PCA plot might not fully represent the actual relationships in the original data, especially when the total variance explained by the principal components is relatively low.

### Generally, the Principal Component Analysis(PCA) steps are:
 
1) Scaling our data (if the columns don’t scale the same). This is important because PCA is an algorithm that is strongly influenced by the size of each column.
2) Calculate Covariance Matrix
3) Calculate Eigenvalues and Eigenvector
4) Sort Eigenvalues and Eigenvector
5) Pick top-2 or top-3 (or any amount of Principal Components that you want) eigenvalues
6) Transform the original data

### How to calculate covariance? 
Here provides a formula for you.

<img src="covariance.png" width="50%" />



### Example codes in R

```r
library(factoextra)
library(FactoMineR)
library(ggpubr) 
data("iris")
head(iris)
iris.pca <- PCA(iris[,-5], graph = FALSE) 
# The quality of representation of the variables of the principal components are called the cos2.
print(iris.pca)
fviz_pca_ind(iris.pca, label="none")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-1.png" width="672" />

```r
fviz_pca_ind(iris.pca,  label="none", habillage=iris$Species)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-2.png" width="672" />

```r
fviz_pca_ind(iris.pca, label="none", habillage=iris$Species,
             addEllipses=TRUE, ellipse.level=0.95)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-3.png" width="672" />

```r
fviz_pca_biplot(iris.pca, 
                habillage = iris$Species, addEllipses = TRUE,
                col.var = "red", alpha.var ="cos2",
                label = "var") +
  scale_color_brewer(palette="Dark2")+
  theme_minimal()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-4.png" width="672" />

```r
fviz_screeplot(iris.pca, ncp=10)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-5.png" width="672" />

```r
plot(iris.pca, choix = "var")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-6.png" width="672" />

```r
fviz_pca_var(iris.pca, col.var="contrib")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-7.png" width="672" />

```r
fviz_pca_var(iris.pca, alpha.var="contrib")+
  theme_minimal()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-8.png" width="672" />

```r
# Coordinates of individuals on the principal components
head(iris.pca$ind$coord)
head(iris.pca$ind$cos2)
head(iris.pca$ind$contrib)
plot(iris.pca, choix = "ind")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-9.png" width="672" />

```r
fviz_pca_ind(iris.pca)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-10.png" width="672" />

```r
fviz_pca_biplot(iris.pca,  geom = "text")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-11.png" width="672" />



### References 

https://bioconductor.org/packages/release/bioc/vignettes/roFFF/inst/doc/ropls-vignette.html
https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html
https://towardsdatascience.com/the-mathematics-behind-principal-component-analysis-fff2d7f4b643
https://www.mathsisfun.com/algebra/eigenvalue.html
http://www.sthda.com/english/wiki/wiki.php?id_contents=7866
http://blog.codinglabs.org/articles/pca-tutorial.html
https://www.youtube.com/watch?v=GEn-_dAyYME
https://towardsdatascience.com/the-mathematics-behind-principal-component-analysis-fff2d7f4b643
http://blog.codinglabs.org/articles/pca-tutorial.html
http://www.sthda.com/english/wiki/wiki.php?id_contents=7851
https://peterbloem.nl/blog/pca





