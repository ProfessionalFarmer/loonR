% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methylation.R
\name{ChAMP_QC_Pipeline_Frome_Beta_Value}
\alias{ChAMP_QC_Pipeline_Frome_Beta_Value}
\title{ChAMP QC Pipeline for beta matrix}
\usage{
ChAMP_QC_Pipeline_Frome_Beta_Value(
  Sample.beta.df = NULL,
  Sample.Group = "",
  Slide = "",
  arraytype = c("450K", "EPIC"),
  CpG.GUI = FALSE,
  QC.GUI = FALSE,
  combat = FALSE,
  batchname = c("Slide"),
  cores = 50
)
}
\arguments{
\item{Sample.beta.df}{Column is sample}

\item{arraytype}{450K or EPIC}

\item{CpG.GUI}{If use CpG.GUI}

\item{combat}{If perform combat to remove batch effect}

\item{batchname}{Batch variable name. Default is c("Slide")}

\item{cores}{}
}
\description{
ChAMP QC Pipeline for beta matrix
}
\examples{

library(ChAMP)
testDir=system.file("extdata",package="ChAMPdata")
myImport <- champ.import(testDir)

beta.df = myImport$beta
group = myImport$pd$Sample_Group

########### Step 1
res = loonR::ChAMP_QC_Pipeline_Frome_Beta_Value(Sample.beta.df=beta.df, Sample.Group=group)
myLoad = res$myLoad

########### Step 2
res = res$followAnalysisImputeQCNormCombat(arraytype="450K")
head(res$myNorm)
}
