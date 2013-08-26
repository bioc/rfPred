#' Motor of \code{rfPred_scores}
#'
#' Motor of \code{rfPred_scores}
#'
#' @param variant_list Variants list in a \code{data.frame} containing 4 or 5 columns:
#' chromosome number, hg19 genomic position on the chromosome, reference nucleotid, variant nucleotid and uniprot protein identifier (optional)
#' @param data Path to the compressed TabixFile, either on the server (default) or on the user's computer
#' @param index Path to the index of the TabixFile, either on the server (default) or on the user's computer
#' @param all.col \code{TRUE} to return all available information, \code{FALSE} to return a more compact result (the most informative columns, see Value)
#' @param file.export Optional, name of the CSV file in which export the results (default is \code{NULL})
#' @param n.cores number of cores to use when scaning the TabixFile, can be efficient for large request (default is 1)
#' @note This function is called by the \code{rfPred_scores} S4 method
#' @return see the \code{\link{rfPred_scores}} function

# last updated: july 17, 2013

rfPred_scores_motor=function(variant_list,data,index,all.col,file.export,n.cores){
  # warning message if request on the server
  if (nrow(variant_list)>1000 & data=="http://www.sbim.fr/rfPred/all_chr_rfPred.txtz"){
    stop("For variants lists with more than 1000 rows, please download the TabixFile/index on http://www.sbim.fr/rfPred/ and run rfPred locally")
  }
  
  # formating columns in the right class
  variant_list[,c(1,3:ncol(variant_list))]=apply(variant_list[,c(1,3:ncol(variant_list))],2,as.character)
  variant_list[,2]=as.numeric(as.character(variant_list[,2]))
  names(variant_list)=c("chromosome","position_hg19","reference","alteration",if(ncol(variant_list)==5){"proteine"}else{NULL})
  variant_list=as.data.table(variant_list)

  # output columns
  names=c("chromosome","position_hg19","reference","alteration","proteine","aaref","aaalt",
          "rfPred_score","SIFT_score","MutationTaster_score","Polyphen2_score","PhyloP_score","LRT_score","aapos",
          "Uniprot_id","genename","position_hg18","Polyphen2_HDIV_score","Polyphen2_HDIV_pred","Polyphen2_HVAR_score",
          "Polyphen2_HVAR_pred","MutationTaster_pred","phyloP","LRT_Omega","LRT_pred")
  if (all.col){names.out=names}else{names.out=names[1:14]}

  # request on TabixFile (chr+pos)
  param=reduce(GRanges(variant_list$chromosome,IRanges(start=variant_list$position_hg19,width=1)),drop.empty.ranges=TRUE)
  if (n.cores==1){
    res.rst=scanTabix(TabixFile(data,index),param=param)  
  } else{
    Clust=makeCluster(n.cores)
    res.rst=unlist(parLapply(Clust,clusterSplit(Clust,param), 
                     function(variants,file){require(Rsamtools);return(scanTabix(file,param=variants))},TabixFile(data,index)))  
    stopCluster(Clust)
  }
  
  res.rst=.Call("outstring",unlist(res.rst),as.integer(length(names)),"\t",PACKAGE="rfPred")
  colnames(res.rst)=names
  res.rst=as.data.frame(res.rst,stringsAsFactors=FALSE)[,names.out]
  res.rst[,2]=as.numeric(res.rst[,2]) 
  res.rst=as.data.table(res.rst)
  
  # out for no matching on chr+pos
  out=merge(unique(variant_list[,1:2,with=FALSE]),res.rst,all.x=TRUE,by=c("chromosome","position_hg19"),allow.cartesian=TRUE)
  out[apply(out[,3:ncol(out),with=FALSE],1,function(x){all(is.na(x))}),3:ncol(out)]="no matching"
  out=out[apply(out[,3:ncol(out),with=FALSE],1,function(x){all(x=="no matching")}),]
  variant_list=merge(variant_list,unique(res.rst[,1:2,with=FALSE]),by=c("chromosome","position_hg19"),allow.cartesian=TRUE)
   
  # matching on columns 3 and 4 (and possibly 5)
  out=rbind(out,merge(variant_list,res.rst,all.x=TRUE,by=names(variant_list),allow.cartesian=TRUE))
  if (ncol(variant_list)==5){
    variant_list[,5]="."
    out=rbind(out,merge(variant_list,res.rst,by=names(variant_list),allow.cartesian=TRUE))
  }
  out[apply(out[,(ncol(variant_list)+1):ncol(out),with=FALSE],1,function(x){all(is.na(x))}),(ncol(variant_list)+1):ncol(out)]="no matching"

  # sorting the output
  out=as.data.frame(out[order(factor(out$chromosome,levels=c(1:22,"X","Y")),as.numeric(out$position_hg19)),])
  rownames(out)=1:nrow(out)

  # writing results in a CSV file
  if (!is.null(file.export)){write.csv2(out,file=file.export,row.names=FALSE)}

  return(out)
}
