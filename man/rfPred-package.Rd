\docType{package}
\name{rfPred-package}
\alias{rfPred-package}
\title{Assign functional prediction rfPred scores to human missense variants
(random forest method based on SIFT, Polyphen2, PhyloP, LRT and Mutation Taster)}
\description{
  The package provides a function which returns the rfPred
  score for a list of non-synonymous missense variants. All
  the rfPred scores are pre-calculated and stored in a
  \code{TabixFile} available on a server and which can be
  downloaded for using the package while not connected on
  the Internet. The package does not work without an access
  to the \code{TabixFile}. However, a toy example on the
  chromosome Y is available within the package to test the
  \code{rfPred_scores} function. curves with numbers of
  subjects at risk, compare data sets, display
  spaghetti-plot, build multi-contingency tables...
}
\author{
  Fabienne Jabot-Hanin, Hugo Varet and Jean-Philippe Jais
}
\references{
  dbNSFP database: Liu X, Jian X and Boerwinkle E. 2011.
  dbNSFP: a lightweight database of human non-synonymous
  SNPs and their functional predictions. Human Mutation.
  32:894-899.

  rfPred method: Jabot-Hanin F, Varet H, Tores F and Jais
  J-P. 2013. rfPred: a new meta-score for functional
  prediction of missense variants in human exome
  (submitted).
}

