# rcellmarker
# richR <a href="https://travis-ci.org/guokai8/rcellmarker"><img src="https://travis-ci.org/guokai8/rcellmarker.svg" alt="Build status"></a>  [![Project Status:](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)  [![](https://img.shields.io/badge/devel%20version-0.0.8-green.svg)](https://github.com/guokai8/rcellmarker)  ![Code Size:](https://img.shields.io/github/languages/code-size/guokai8/rcellmarker)
## Description
_rcellmarker_ provides method to identify cell type based on single cell sequencing data.  Now only support human and mouse.
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
```   
## Note
The _rcellmarker_ package use the __CellMarker__ and __PanglaoDB__ database as the reference with ID cleaning and merging. We also add information from the newest public papers related with single cell field from NCBI. The package is still under development. Other species will be supported soon.

## Contact information

For any questions please contact guokai8@gmail.com
