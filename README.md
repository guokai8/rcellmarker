# rcellmarker
# richR <a href="https://travis-ci.org/guokai8/rcellmarker"><img src="https://travis-ci.org/guokai8/rcellmarker.svg" alt="Build status"></a>  [![Project Status:](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)  [![](https://img.shields.io/badge/devel%20version-0.0.9-green.svg)](https://github.com/guokai8/rcellmarker)  ![Code Size:](https://img.shields.io/github/languages/code-size/guokai8/rcellmarker)
## Description
_rcellmarker_ provides method to identify cell type based on single cell sequencing data. Since most methods try to annotate cell types manually after clustering the single-cell RNA-seq data. Such methods are labor-intensive and heavily rely on user expertise, which may lead to inconsistent results. Here, we present _rcellmarker_ package -- an automatic tool to annotate cell types from single-cell RNA-seq data. Now only support human and mouse.
## Installation
```
library(devtools)
install_github("guokai8/rcellmarker")
``` 
## Quick tour
```{r}
set.seed(123)   
library(rcellmarker)   
gene=sample(unique(human$SYMBOL),20)
res<-cells(gene,species = "human",keytype ="SYMBOL")
head(res)
#if you have the result include the cluster information, padj value and avg_logFC named as "single"
head(getdetail(res,single))
## if yuo have the results from Seurat of cellranger named as 'df'
res <- cellMarker(df,type='seurat',species='human',keytype='SYMBOL',weight=1)
# or
res <- cellMarker(df,type='cellranger',species='human',keytype='SYMBOL',weight=100)
# or you just have two columns include cluster and gene name
res <- cellMarker(df,type='custom',species='human',keytype='SYMBOL')
```   
## Note
The _rcellmarker_ package use the __CellMarker__ and __PanglaoDB__ database as the reference with ID cleaning and merging. We also add information from the newest public papers related with single cell field from NCBI. The package is still under development. Other species will be supported soon.

## Contact information

For any questions please contact guokai8@gmail.com
