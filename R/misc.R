##' @method as.data.frame cellResult
##' @export
as.data.frame.cellResult <- function(x, ...) {
    as.data.frame(x@result, ...)
}
##' @method row.names cellResult
##' @export
row.names.cellResult <- function(x, ...) {
    row.names(x@result)
}
##' @method names cellResult
##' @export
names.cellResult <- function(x, ...) {
    names(x@result)
}
##' @importFrom utils head
##' @method head cellResult
##' @export
head.cellResult <- function(x, n=6L, ...) {
    head(x@result, n, ...)
}
##' @importFrom utils tail
##' @method tail cellResult
##' @export
tail.cellResult <- function(x, n=6L, ...) {
    tail(x@result, n, ...)
}
##' @method dim cellResult
##' @export
dim.cellResult <- function(x) {
    dim(x@result)
}
##' @method [ cellResult
##' @export
`[.cellResult` <- function(x, i, j) {
    x@result[i,j]
}
##' @method $ cellResult
##' @export
`$.cellResult` <-  function(x, name) {
    x@result[, name]
}
##' @method detail cellResult
##' @export
detail.cellResult<-function(x){
    as.data.frame(x@detail)
}
##' @method result cellResult
##' @export
result.cellResult<-function(x){
    as.data.frame(x@result)
}
##' @importFrom dplyr left_join
getdetail<-function(rese,resd,sep){
    if(!is.data.frame(resd)){
        resd=data.frame(gene=resd)
    }
    if(!("gene"%in%colnames(resd))){
        resd$gene=rownames(resd)
    }
    gene<-strsplit(as.vector(rese$GeneID),split=sep)
    names(gene)<-rese$Annot
    gened<-data.frame("TERM"=rep(names(gene),times=unlist(lapply(gene,length))),
                      "Annot"=rep(resultFis$Term,times=unlist(lapply(gene,length))),
                      "GeneID"=unlist(gene),row.names=NULL,
                      "Pvalue"=rep(resultFis$Pvalue,times=unlist(lapply(gene,length))),
                      "Padj"=rep(resultFis$Padj,times=unlist(lapply(gene,length)))
    )
    gened$GeneID<-as.character(gened$GeneID)
    res<-left_join(gened,resd,by=c("GeneID"="gene"))
    return(res)
}
##'
##'
setAs(from = "cellResult", to = "data.frame", def = function(from){
    resultFis <- from@result
    resultFis
})

.getdata <-function(species){
    species = tolower(species)
    if(species=='human'){
        data(human)
        dat <- human
    }else{
        data(mouse)
        dat <- mouse
    }
}

#' reverse List
#' @param lhs list with names
#' @export
#' @author Kai Guo
reverseList<-function(lhs){
    lhs_n<-rep(names(lhs),times=lapply(lhs,function(x)length(x)))
    res<-sf(as.data.frame(cbind(lhs_n,unlist(lhs))))
    return(res)
}
