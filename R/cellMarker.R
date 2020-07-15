#' assign cell type based on cell cluster results
#' @importFrom magrittr %>%
#' @importFrom purrr safely
#' @importFrom purrr map
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom tidyr gather
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr top_n
#' @param x input file for marker annotation
#' @param species species for annotation
#' @param type source of marker genes (seurat(default), cellranger, custom)
#' @param keytype keytype for input genes
#' @param weight weight threshold for marker filtering 
#' @param format file format for user supplied data.frame(long: Cluster, gene;
#'  wide: gene, following cluster name as column name)
#' @param cluster clutser number (default: NULL for annotate all clusters)
#' @param tissue tissue for annotation (default: NULL to use all tissues)
#' @param topn the number of cell type to list fro each cluster
#' @param padj adjust p value threshold
#' @param minSize minimal number of genes included in significant cell type
#' @param maxSize maximum number of genes included in significant cell type
#' @param p.adjust.methods pvalue adjust method(default: "BH")
#' @export
#' @author Kai Guo
cellMarker <- function(x, type = 'seurat', species="human", keytype = 'SYMBOL', 
                       weight = NULL, format="long",
                       cluster = NULL,tissue = NULL, topn = 3,
                       padj = 0.05, minSize=2,maxSize=500,
                       p.adjust.methods = "BH"){
    options(warn = -1)
        cells_ <- safely(cells, otherwise = .empty_class())
        if(type == 'cellranger'){
            if(is.null(weight)) weight <- 100
            colnames(x)[1:2]<-c('GeneID','GeneName')
            x <- x%>%select(GeneName,contains("Weight"))
            colnames(x) <- sub('_Weight','',gsub('[\\.| ]','_',colnames(x)))
            x <- x%>%gather(Cluster,val,-GeneName)%>%filter(val>=weight)%>%
                select(Cluster,GeneName)%>%
                group_by(Cluster)%>%nest()
            if(!is.null(cluster)){
                cluster <- paste0("Cluster_",cluster)
                x <- x%>%filter(Cluster%in%cluster)
            }
            x <- x %>%mutate(cellType=map(data,
                function(y)result(cells_(y$GeneName,species=species,
                            keytype=keytype,minSize=minSize,padj=padj,
                            maxSize=maxSize,
                            p.adjust.methods=p.adjust.methods)$result))) 
            x <- x%>%select(Cluster,cellType)%>%unnest(cellType)%>%
                group_by(Cluster)%>%top_n(topn,wt=Padj)
        }else if(type == "seurat"){
            if(is.null(weight)) weight <- 1
            x <-x%>%filter(avg_logFC >= weight,p_val_adj<padj)%>%
                select(cluster,gene)%>%group_by(cluster)%>%nest()
            if(!is.null(cluster)){
                cl <- cluster
                x <- x%>%filter(cluster%in%cl)
            }
            x <- x %>%mutate(cellType=map(data,
                    function(y)result(cells_(y$gene,species=species,
                    keytype=keytype,minSize=minSize,padj=padj,
                    maxSize=maxSize,
                    p.adjust.methods=p.adjust.methods)$result))) 
            x <- x%>%select(cluster,cellType)%>%unnest(cellType)%>%
                group_by(cluster)%>%top_n(topn,wt=Padj)
        }else{
            colnames(x)[1]<-'gene'
            if(format=="wide"){
                x <- x%>%gather(Cluster,val,-gene)%>%select(Cluster,gene)
            }
            x <- x%>%group_by(Cluster)%>%nest()
            if(!is.null(cluster)){
                cl <- cluster
                x <- x%>%filter(Cluster%in%cl)
            }
            x <- x %>%mutate(cellType=map(data,
                    function(y)result(cells_(y$gene,species=species,
                    keytype=keytype,minSize=minSize,padj=padj,
                    maxSize=maxSize,
                    p.adjust.methods=p.adjust.methods)$result))) 
            x <- x%>%select(Cluster,cellType)%>%unnest(cellType)%>%
                group_by(Cluster)%>%top_n(topn,wt=Padj)
        }
    as.data.frame(x)
}