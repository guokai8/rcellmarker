# rcellmarker
# richR <a href="https://travis-ci.org/guokai8/rcellmarker"><img src="https://travis-ci.org/guokai8/rcellmarker.svg" alt="Build status"></a>  [![Project Status:](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)  [![](https://img.shields.io/badge/devel%20version-0.0.7-green.svg)](https://github.com/guokai8/rcellmarker) 
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
```       
