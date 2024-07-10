#' cells function for cell type identification 
#' @importFrom magrittr %>%
#' @importFrom stats p.adjust
#' @importFrom methods new
#' @param x vector contains gene names 
#' @param species species name
#' @param tissue tissue type (default NULL)
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant cell type
#' @param maxSize maximum number of genes included in significant cell type
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @export
#' @author Kai Guo
cells<-function(x,species='human',db='default',keytype="SYMBOL", tissue = NULL, padj=0.05, pvalue=NULL,
                minSize=3,maxSize=500,
                padj.method="BH",sep = ","){
    annot <- .getdata(species=species,db=db)
    annot <- na.omit(annot)
    keytype <- toupper(keytype)
    if(!is.null(tissue)){
        annot <- annot[grepl(tissue,annot$tissueType,ignore.case = T),]
    }
    annot <- annot[,c(keytype,'cellType')]
    if(sum(x%in%annot[,1])==0){
        return(.empty_class())
    }
    ao2gene<-sf(annot)
    ao2gene_num<-name_table(ao2gene)
    gene2ao<-sf(annot[,c(2,1)])
    input=as.vector(x)
    fgene2ao=gene2ao[input]
    fao2gene=reverseList(fgene2ao)
    k=name_table(fao2gene)
    n=length(unique(unlist(fao2gene)))
    M=ao2gene_num[names(k)]
    N=length(unique(annot[,1]))
    rhs<-hyper_bench_vector(k,M,N,n)
    lhs<-p.adjust(rhs,method=padj.method)
    rhs_gene<-unlist(lapply(fao2gene, function(x)paste(unique(x),sep="",collapse = sep)))
    resultFis<-data.frame("cellType"=names(rhs),"Annotated"=M[names(rhs)],
                          "Significant"=k[names(rhs)],"Pvalue"=as.vector(rhs),"Padj"=lhs,
                          "GeneID"=rhs_gene[names(rhs)])
    resultFis<-resultFis[order(resultFis$Pvalue),]
    if(!is.null(pvalue)){
        resultFis<-resultFis[resultFis$Pvalue<pvalue,]
        padj <- numeric()
    }else{
        resultFis<-resultFis[resultFis$Padj<padj,]
        pvalue <- numeric()
    }
    resultFis<-resultFis[resultFis$Significant<=maxSize,]
    resultFis<-resultFis[resultFis$Significant>=minSize,]
    rownames(resultFis)<-resultFis$cellType
    gene<-strsplit(as.vector(resultFis$GeneID),split=sep)
       # names(gene)<-resultFis$cellType
    gened<-data.frame("cellType"=rep(resultFis$cellType,times=unlist(lapply(gene,length))),
                          "GeneID"=unlist(gene),row.names=NULL,
                          "Pvalue"=rep(resultFis$Pvalue,times=unlist(lapply(gene,length))),
                          "Padj"=rep(resultFis$Padj,times=unlist(lapply(gene,length)))
        )
    gened$GeneID<-as.character(gened$GeneID)
    result<-new("cellResult",
                result = resultFis,
                detail = gened,
                species = species, 
                pvalueCutoff   = pvalue,
                pAdjustMethod  = padj.method,
                padjCutoff   = padj,
                gene           = input,
                keytype        = keytype,
                sep = sep
    )
    return(result)
}
#' show support tissues
#' @importFrom utils data
#' @param species species name
#' @export
#' @author Kai Guo
list.tissue <- function(species='human'){
    dat <- .getdata(species=species)
    data.frame('Tissue'=sort(unique(humancells$tissueType)))
}





