##' detail generic
##' @param x cellResult object
##' @return detail return detial for these significant genes
##' @export
detail<-function(x){
    UseMethod("detail",x)
}
##' result generic
##' @param x cellResult object
##' @return result return dataframe and print summary
##' @export
result<-function(x){
    UseMethod("result",x)
}