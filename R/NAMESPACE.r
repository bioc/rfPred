#' @exportPattern ^[a-zA-Z]
#' @useDynLib rfPred
#' @import methods
#' @importFrom utils read.table write.csv2
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @importFrom data.table as.data.table
#' @importFrom IRanges reduce IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom Rsamtools scanTabix TabixFile
#' @importFrom parallel makeCluster stopCluster parLapply clusterSplit
#' @exportMethod rfPred_scores
NULL
