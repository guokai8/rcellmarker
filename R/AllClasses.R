##' Class "cellResult"
##' This class represents the result of enrichment analysis.
##'
##'
##' @name cellResult-class
##' @aliases cellResult-class
##'   show,cellResult-method plot,cellResult-method
##'   summary,cellResult-method
##'
##' @docType class
##' @slot result enrichment analysis results
##' @slot detail genes included in significant terms and original information
##' @slot species species name
##' @slot pvalueCutoff cutoff pvalue
##' @slot pAdjustMethod pvalue adjust method
##' @slot padjCutoff pvalue adjust cutoff value
##' @slot gene Gene IDs
##' @slot keytype Gene ID type
##' @slot sep character string used to separate the genes when concatenating 
##' @exportClass cellResult
##' @author Kai Guo
##' @keywords classes
setClass("cellResult",
         representation=representation(
             result         = "data.frame",
             detail         = "data.frame",
             species       = "character",
             pvalueCutoff   = "numeric",
             pAdjustMethod  = "character",
             padjCutoff   = "numeric",
             gene           = "character",
             keytype        = "character",
             sep = "character"
         )
)