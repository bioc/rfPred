#' Assign functional prediction rfPred scores to human missense variants
#'
#' rfPred is a statistical method which combines 5 algorithms predictions in a
#' random forest model: SIFT, Polyphen2, LRT, PhyloP and MutationTaster. These scores
#' are available in the dbNFSP database for all the possible missense variants in hg19
#' version, and the package rfPred gives a composite score more reliable than each of
#' the isolated algorithms.
#'
#' @param variant_list A variants list in a \code{data.frame} containing 4 or 5 columns:
#' chromosome number, hg19 genomic position on the chromosome, reference nucleotid, variant nucleotid and uniprot protein identifier (optional);
#' or a character string of the path to a VCF (Variant Call Format) file; or a \code{GRanges} object with metadata containing textually \code{reference}, 
#' \code{alteration} and \code{proteine} (optional) columns names
#' for reference and alteration
#' @param data Path to the compressed TabixFile, either on the server (default) or on the user's computer
#' @param index Path to the index of the TabixFile, either on the server (default) or on the user's computer
#' @param all.col \code{TRUE} to return all available information, \code{FALSE} to return a more compact result (the most informative columns, see Value)
#' @param file.export Optional, name of the CSV file in which export the results (default is \code{NULL})
#' @param n.cores number of cores to use when scaning the TabixFile, can be efficient for large request (default is 1)
#' @return The variants list with the assigned rfPred scores, as well as the scores used to build rfPred meta-score: SIFT, phyloP, MutationTaster, LRT (transformed) and Polyphen2 (corresponding to Polyphen2_HVAR_score). The data frame returned contains these columns:
#'   \item{chromosome}{chromosome number}
#'   \item{position_hg19}{physical position on the chromosome as to hg19 (1-based coordinate)}
#'   \item{reference}{reference nucleotide allele (as on the + strand)}
#'   \item{alteration}{alternative nucleotide allele (as on the + strand)}
#'   \item{proteine}{Uniprot accession number}
#'   \item{aaref}{reference amino acid}
#'   \item{aaalt}{alternative amino acid}
#'   \item{aapos}{amino acid position as to the protein}
#'   \item{rfPred_score}{rfPred score betwen 0 and 1 (higher it is, higher is the probability of pathogenicity)}
#'   \item{SIFT_score}{SIFT score between 0 and 1 (higher it is, higher is the probability of pathogenicity contrary to the original SIFT score) = 1-original SIFT score}
#'   \item{Polyphen2_score}{Polyphen2 (HVAR one) score between 0 and 1, used to calculate rfPred (higher it is, higher is the probability of pathogenicity)}
#'   \item{MutationTaster_score}{MutationTaster score between 0 and 1 (higher it is, higher is the probability of pathogenicity)}
#'   \item{PhyloP_score}{PhyloP score between 0 and 1 (higher it is, higher is the probability of pathogenicity): PhyloP_score=1-0.5x10^phyloP if phyloP>0 or PhyloP_score=0.5x10^-phyloP if phyloP<0}
#'  \item{LRT_score}{LRT score between 0 and 1 (higher it is, higher is the probability of pathogenicity): LRT_score=1-LRToriginalx0.5 if LRT_Omega<1 or LRT_score=LRToriginalx0.5 if LRT_Omega>=1}
#' The following columns are also returned if \code{all.col} is \code{TRUE}:
#'   \item{Uniprot_id}{Uniprot ID number}
#'   \item{genename}{gene name}
#'   \item{position_hg18}{physical position on the chromosome as to hg18 (1-based coordinate)}
#'   \item{Polyphen2_HDIV_score}{Polyphen2 score based on HumDiv, i.e. hdiv_prob. The score ranges from 0 to 1: the corresponding prediction is "probably damaging" if it is in [0.957,1]; "possibly damaging" if it is in [0.453,0.956]; "benign" if it is in [0,0.452]. Score cut-off for binary classification is 0.5, i.e. the prediction is "neutral" if the score is lower than 0.5 and "deleterious" if the score is higher than 0.5. Multiple entries separated by ";"}
#'   \item{Polyphen2_HDIV_pred}{Polyphen2 prediction based on HumDiv: \code{D} (probably damaging), \code{P} (possibly damaging) and \code{B} (benign). Multiple entries separated by ";"}
#'   \item{Polyphen2_HVAR_score}{Polyphen2 score based on HumVar, i.e. hvar_prob. The score ranges from 0 to 1, and the corresponding prediction is "probably damaging" if it is in [0.909,1]; "possibly damaging" if it is in [0.447,0.908]; "benign" if it is in [0,0.446]. Score cut-off for binary classification is 0.5, i.e. the prediction is "neutral" if the score is lower than 0.5 and "deleterious" if the score is higher than 0.5. Multiple entries separated by ";"}
#'   \item{Polyphen2_HVAR_pred}{Polyphen2 prediction based on HumVar: \code{D} (probably damaging), \code{P} (possibly damaging) and \code{B} (benign). Multiple entries separated by ";"}
#'   \item{MutationTaster_pred}{MutationTaster prediction: \code{A} (disease_causing_automatic), \code{D} (disease_causing), \code{N} (polymorphism) or \code{P} (polymorphism_automatic)}
#'   \item{phyloP}{original phyloP score}
#'   \item{LRT_Omega}{estimated nonsynonymous-to-synonymous-rate ratio}
#'   \item{LRT_pred}{LRT prediction, \code{D}(eleterious), \code{N}(eutral) or \code{U}(nknown)}
#' @author Fabienne Jabot-Hanin, Hugo Varet and Jean-Philippe Jais
#' @references Jabot-Hanin F, Varet H, Tores F and Jais J-P. 2013. rfPred: a new meta-score for functional prediction of missense variants in human exome (submitted).
#' 
#' @export
#' @docType methods
#' @rdname rfPred_scores-methods
#'
#' @examples
#' # from a data.frame without uniprot protein identifier
#' data(variant_list_Y)
#' res=rfPred_scores(variant_list = variant_list_Y[,1:4],
#'         data = system.file("extdata", "chrY_rfPred.txtz", package="rfPred",mustWork=TRUE),
#'         index = system.file("extdata", "chrY_rfPred.txtz.tbi", package="rfPred",mustWork=TRUE))
#' # from a data.frame with uniprot protein identifier
#' res2=rfPred_scores(variant_list = variant_list_Y,
#'         data = system.file("extdata", "chrY_rfPred.txtz", package="rfPred",mustWork=TRUE),
#'         index = system.file("extdata", "chrY_rfPred.txtz.tbi", package="rfPred",mustWork=TRUE))
#' # from a VCF file
#' res3=rfPred_scores(variant_list = system.file("extdata", "example.vcf", package="rfPred",mustWork=TRUE),
#'         data = system.file("extdata", "chrY_rfPred.txtz", package="rfPred",mustWork=TRUE),
#'         index = system.file("extdata", "chrY_rfPred.txtz.tbi", package="rfPred",mustWork=TRUE))
#' # from a GRanges object
#' data(example_GRanges)
#' res4=rfPred_scores(variant_list = example_GRanges,
#'         data = system.file("extdata", "chrY_rfPred.txtz", package="rfPred",mustWork=TRUE),
#'         index = system.file("extdata", "chrY_rfPred.txtz.tbi", package="rfPred",mustWork=TRUE))
#' 

setGeneric(
 name = "rfPred_scores",
 def = function(variant_list,data="http://www.sbim.fr/rfPred/all_chr_rfPred.txtz",
                index="http://www.sbim.fr/rfPred/all_chr_rfPred.txtz.tbi",
                all.col=FALSE,file.export=NULL,n.cores=1){standardGeneric("rfPred_scores")}
)

#' @rdname rfPred_scores-methods
#' @aliases rfPred_scores,data.frame-method
#' @usage rfPred_scores(variant_list,data="http://www.sbim.fr/rfPred/all_chr_rfPred.txtz",
#'                 index="http://www.sbim.fr/rfPred/all_chr_rfPred.txtz.tbi",
#'                 all.col=FALSE,file.export=NULL)
setMethod(
  f = "rfPred_scores",
  signature = c("data.frame"),
  definition = function(variant_list,data,index,all.col,file.export,n.cores){
    rfPred_scores_motor(variant_list,data,index,all.col,file.export,n.cores)
  }
)

#' @rdname rfPred_scores-methods
#' @aliases rfPred_scores,character-method
setMethod(
  f = "rfPred_scores",
  signature = c("character"),
  definition = function(variant_list,data,index,all.col,file.export,n.cores){
    d=read.table(variant_list,stringsAsFactors=FALSE)
    d=d[,c(1,2,4,5)]
    d[,1]=sub("chr","",d[,1])
    d[,2]=sub(" ","",d[,2])
    rfPred_scores_motor(d,data,index,all.col,file.export,n.cores)
  }
)

#' @rdname rfPred_scores-methods
#' @aliases rfPred_scores,GRanges-method
setMethod(
  f = "rfPred_scores",
  signature = c("GRanges"),
  definition = function(variant_list,data,index,all.col,file.export,n.cores){
    seqlevels(variant_list)=sub("chr","",seqlevels(variant_list))
    d=as.data.frame(variant_list,stringsAsFactors=FALSE)
    d=d[,c("seqnames","start","reference","alteration",if(any(names(d)=="proteine")){"proteine"}else{NULL})]
    rfPred_scores_motor(d,data,index,all.col,file.export,n.cores)
  }
)
