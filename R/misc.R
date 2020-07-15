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
##' get detail from with cellResult and combine with other information
##' @importFrom dplyr left_join
##' @param x cellResult object from cells function
##' @param y  a data frame with gene name and other information
##' @param sep cellResult object sep
##' @author Kai Guo
##' @export
getdetail<-function(x,y,sep=","){
    if(!is.data.frame(y)){
        y=data.frame(gene=y)
    }
    if(!("gene"%in%colnames(y))){
        y$gene=rownames(y)
    }
    if(is.data.frame(x)){
        sep=sep
        result <- x
    }else{
        sep = x@sep
        result <- x@result
    }
    gene<-strsplit(as.vector(x$GeneID),split=sep)
    gened<-data.frame("cellType"=rep(result$cellType,times=unlist(lapply(gene,length))),
                      "GeneID"=unlist(gene),row.names=NULL,
                      "Pvalue"=rep(result$Pvalue,times=unlist(lapply(gene,length))),
                      "Padj"=rep(result$Padj,times=unlist(lapply(gene,length)))
    )
    gened$GeneID<-as.character(gened$GeneID)
    res<-left_join(gened,y,by=c("GeneID"="gene"))
    return(res)
}
##' Functions to coerce cellResult to data.frame
##' @method as.data.frame cellResult
##' @export
as.data.frame.cellResult <- function(x, ...) {
    as.data.frame(x@result, ...)
}

#' load the data based on the species name
#' @param species species name
#' @author Kai Guo
.getdata <-function(species){
    species = tolower(species)
    if(species=='human'){
        data(human)
        dat <- humancells
    }else{
        data(mouse)
        dat <- mousecells
    }
    dat
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
.empty_class <- function(){
    new("cellResult",
        result = data.frame(),
        detail = data.frame(),
        species = "", 
        pvalueCutoff   = 0.05,
        pAdjustMethod  = "BH",
        padjCutoff   = 0.05,
        gene           = "",
        keytype        = "",
        sep = ""
    )
}
##' get the unique cluster markers for each cluster 
##' @importFrom dplyr distinct
##' @param x cellResult object or result from cellMarker
##' @param sep character string used to separate the genes in GeneID column
##' @param .unique filter duplicate gene id or not
##' @author Kai Guo
##' @export
marker <- function(x,sep=",",.unique=TRUE){
        x <- as.data.frame(x)
        gene<-strsplit(as.vector(x$GeneID),split=sep)
        res <- data.frame("Cluster"= rep(x[,1],times=unlist(lapply(gene,length))),
                          'cellType'=rep(x$cellType,times=unlist(lapply(gene,length))),
                          "GeneID" = unlist(gene))
        if(isTRUE(.unique)){
            res <- res%>%distinct(GeneID,.keep_all = T)
        }
        res
}
