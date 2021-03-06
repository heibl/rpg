% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compareMSAs.R
\name{compareMSAs}
\alias{compareMSAs}
\title{Compare reference MSAs with alternative MSAs}
\usage{
compareMSAs(ref, com, dir_path)
}
\arguments{
\item{ref}{of class data.frame, is the reference MSA ('BASE MSA') with sequences as columns}

\item{com}{like ref, but 1 alternative MSA}

\item{dir_path}{directory with multiple alternative MSA files}
}
\value{
list containing following scores:

mean_scores residue pair score and mean column score

column_score

residue_column_score

residue_pair_residue_score

residual_pair_sequence_pair_score

residual_pair_sequence_score

residue_pair_score
}
\description{
MSA reliability scores (Penn et al. 2010)

Wrapper function for program msa_set_score v2.01 of the GUIDANCE program (see reference). Copyright: To modify the code, or use parts of it for other purposes, permission should be requested. Please contact Tal Pupko: talp@post.tau.ac.il. Please note that the use of the GUIDANCE program is for academic use only.

C code computing basic MSA comparision. The most basic is the residue pairs residue score, which checks if residue pairs combinations are correctly aligned in both MSAs. From this the residue score, residue column score (GUIDANCE score), residue sequence score are computed. It also calulates the column score (CS), which simply checks if a column is identically aligned in the alternative MSA.
}
\references{
Penn et al. (2010). An alignment confidence score capturing
  robustness to guide tree uncertainty. Molecular Biology and Evolution
  27:1759--1767
}
\seealso{
\code{\link{guidance}}, \code{\link{guidance2}}, \code{\link{HoT}}
}
\author{
Franz-Sebastian Krah
}
