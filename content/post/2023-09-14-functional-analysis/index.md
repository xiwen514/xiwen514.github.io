---
title: Functional analysis
author: Sivan
date: '2023-09-14'
slug: []
categories: []
tags: []
---

Here mainly introduce over representation analysis(ORA), gene set enrichment analysis(GSEA) and gene set variation analysis(GSVA).

### 1. Over representation analysis(ORA)

Hypergeometric test is a statistical method used to determine if the frequency of a specific subset in a sample significantly deviates from what would be expected under random sampling. To be specific, we classify the genes we've selected to see which ones have functions relevant to our research. Then, we conduct a test and use the p-value to determine if the association is random or significant. To be more accurate, we also perform a false discovery rate (FDR) test and use the q-value to further assess which gene functions are most related to our study.

**a four-step approach to map gene lists onto pathways**

- Calculate a local (gene-level) statistic
- Calculate a global (gene set or pathway-level) statistic
- Determine significance of the global statistic
- Adjust for multiple testing

In my work, I usually apply ORA to find functional pathways after differential analysis by comparing the ratio between **genes number in a specific pathway / Foreground genes number** and **genes number in a specific pathway / Background genes number**.

- Foreground genes: These are typically the differentially expressed genes between the experimental and control groups, showing significant changes in expression.
- Background genes: For organisms with a reference genome, it's generally recommended to use all genes from the reference genome. For organisms without a reference, the assembled unigenes are used as the background genes.

#### ORA approaches have three major limitations.
- First, the inclusion criteria for input gene lists are rather arbitrary and typically involves selecting genes that exceed some user-defined statistical cutoff.  This risks excluding potentially important genes that for whatever reason fail to reach statistical significance.  
- Second, ORA approaches use gene names but not any of the rich quantitative information associated with gene expression experiments.  In this way, equal importance is assigned to each an every gene.  
- Third, many of the ORA procedures uses statistical procedures that assume independence among genes: Changes in any gene do not affect or are not influenced by any others.  Clearly, this is unrealistic for biological systems and has the effect of making ORA procedures more prone to erroneous discoveries or false-positives. Single-gene analysis may miss important effects on pathways. 

#### More on gene set collections

Gene Ontology (GO)

- Cellular components (CC)
- Biological processes (BP)
- Molecular functions (MF)

Well curated pathway database

- KEGG pathway
- Biocarta
- GenMAPP
- IPA pathway database

Gene set collections

- MSigDB
- GAzer

#### Example codes in R
clusterProfiler takes as input a significant gene list and a background gene list and performs statistical enrichment analysis using hypergeometric testing.


### 2. Gene set enrichment analysis(GSEA)

GSEA is a method designed to assess whether a predefined set of genes shows significant distribution patterns within a gene list ranked by its association with a phenotype. Rather than focusing solely on individual genes, it evaluates the collective prominence of a gene set at the top or bottom of the ranked list, inferring the set's overall contribution to the phenotype.

After conducting differential gene expression analysis and obtaining logFC values for each gene, these logFC values offer a ranking based on expression differences. By employing GSEA, you can determine if a predefined gene set, such as those from a specific biological process or pathway, is significantly enriched within your ranked list of differentially expressed genes. This assists in discerning functional and pathway differences between the two conditions.

**Why all genes?** The hypothesis of FCS(Functional class scoring) methods is that although large changes in individual genes can have significant effects on pathways (and will be detected via ORA methods), weaker but coordinated changes in sets of functionally related genes (i.e., pathways) can also have significant effects. Thus, rather than setting an arbitrary threshold to identify ‘significant genes’, all genes are considered in the analysis, regardless of whether they are deemed differentially expressed based on arbitrary thresholds. All these genes are compared and enriched against predefined GSEA gene sets (similar to pathways, GO, which represent gene-function associations) to determine their contributions to phenotypes (functions). Hence, from the perspective of gene set enrichment, GSEA is not limited to differentially expressed genes and is theoretically more adept at identifying subtle impacts on biological pathways/functions, even when gene fold-changes are minor.

#### GSEA rationale

a significance analysis of function and expression.
Methods that fall under the SAFE framework use a four-step approach to map gene lists onto pathways
- Calculate a local (gene-level) statistic
- Calculate a global (gene set or pathway-level) statistic
- Determine significance of the global statistic
- Adjust for multiple testing

<img src="gsea.png" width="500" height="300">

The Leading-Edge Subset: Gene sets can be defined by using a variety of methods, but not all of the members of a gene set will typically participate in a biological process. Often it is useful to extract the core members of high scoring gene sets that contribute to the ES. We define the leading-edge subset to be those genes in the gene set S that appear in the ranked list L at, or before, the point where the running sum reaches its maximum deviation from zero (Fig. 1B). The leading-edge subset can be interpreted as the core of a gene set that accounts for the enrichment signal.

- Step 1: Calculation of enrichment score:

An enrichment score for a particular gene set is calculated by walking down the list of log2 fold changes and increasing the running-sum statistic every time a gene in the gene set is encountered and decreasing it when genes are not part of the gene set. The size of the increase/decrease is determined by magnitude of the log2 fold change. Larger (positive or negative) log2 fold changes will result in larger increases or decreases. The final enrichment score is where the running-sum statistic is the largest deviation from zero.

<img src="gseastep1.png" width="400" height="300"style="margin-right:30px;"/><img src="gseastep11.png" width="400" height="300"/>


- Step 2: Estimation of significance:

The significance of the enrichment score is determined using permutation testing, which performs rearrangements of the data points to determine the likelihood of generating an enrichment score as large as the enrichment score calculated from the observed data. Essentially, for this step, the first permutation would reorder the log2 fold changes and randomly assign them to different genes, reorder the gene ranks based on these new log2 fold changes, and recalculate the enrichment score. The second permutation would reorder the log2 fold changes again and recalculate the enrichment score again, and this would continue for the total number of permutations run. Therefore, the number of permutations run will increase the confidence in the signficance estimates.

<img src="gseastep2.png" width="500" height="400"/>

- Step 3: Adjust for multiple test correction

After all gene sets are tested, the enrichment scores are normalized for the size of the gene set, then the p-values are corrected for multiple testing.
The GSEA output will yield the core genes in the gene sets that most highly contribute to the enrichment score. The genes output are generally the genes at or before the running sum reaches its maximum value (eg. the most influential genes driving the differences between conditions for that gene set).



The predefined gene sets originate from the MSigDB database (https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). These sets predominantly consist of functional genes. If one only has non-functional genes (e.g., miRNA, lncRNA, circRNA), GSEA analysis might not be feasible due to the absence of appropriate gene sets.

<img src="gseadatabase.png" width="300" height="400"/>

### 3. Gene set variation analysis(GSVA)

GSVA is a non-parametric, unsupervised analytical method primarily for evaluating gene set enrichment in microarrays and transcriptomes. It works by transforming gene expression matrices across samples into gene set expression matrices to assess the enrichment of metabolic pathways. GSVA quantifies gene set enrichment scores per sample, facilitating subsequent statistical analyses. Using the limma package, one can identify differentially expressed genes between samples. Similarly, applying limma to GSVA results (still a matrix) identifies significantly different gene sets. These differentially expressed gene sets, compared to individual genes, offer richer biological insights and interpretability, which can further aid in areas like tumor subtyping and other biologically-relevant inquiries



### References and Codes

online analysis 

- http://impala.molgen.mpg.de/
- Goprofiler: https://biit.cs.ut.ee/gprofiler/gost
- shinygo: http://bioinformatics.sdstate.edu/go/

visualization 

- https://github.com/jokergoo/KeywordsEnrichment
- https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
- Goprofiler:https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html#gene-list-functional-enrichment-analysis-with-gost

ORA: https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/10_FA_over-representation_analysis.html

GSEA 

- guide: https://www.pathwaycommons.org/guide/primers/data_analysis/gsea/
- paper:https://www.pnas.org/doi/10.1073/pnas.0506580102
- code: https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/11_FA_functional_class_scoring.html

GSVA: https://towardsdatascience.com/decoding-gene-set-variation-analysis-8193a0cfda3

fGSEA: https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

ConsensusPathDB documentation: http://cpdb.molgen.mpg.de/CPDB/tutorial#pathwaya.visconcepts

In Chinese : https://zhuanlan.zhihu.com/p/339046340

