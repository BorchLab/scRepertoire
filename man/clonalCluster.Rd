% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clonalCluster.R
\name{clonalCluster}
\alias{clonalCluster}
\title{Clustering adaptive receptor sequences by edit distance}
\usage{
clonalCluster(
  input.data,
  chain = "TRB",
  sequence = "aa",
  samples = NULL,
  threshold = 0.85,
  group.by = NULL,
  exportGraph = FALSE
)
}
\arguments{
\item{input.data}{The product of \code{\link[=combineTCR]{combineTCR()}},
\code{\link[=combineBCR]{combineBCR()}} or \code{\link[=combineExpression]{combineExpression()}}.}

\item{chain}{Indicate if both or a specific chain should be used -
e.g. "both", "TRA", "TRG", "IGH", "IGL".}

\item{sequence}{Clustering based on either \strong{"aa"} or
\strong{"nt"}.}

\item{samples}{The specific samples to isolate for visualization.}

\item{threshold}{The normalized edit distance to consider.
The higher the number the more similarity of sequence will be
used for clustering.}

\item{group.by}{The column header used for to group contigs.
If (\strong{NULL}), clusters will be calculated across samples.}

\item{exportGraph}{Return an igraph object of connected
sequences (\strong{TRUE}) or the amended input with a
new cluster-based variable (\strong{FALSE}).}
}
\value{
Either amended input with edit-distanced clusters added
or igraph object of connect sequences
}
\description{
This function uses edit distances of either the nucleotide or amino acid
sequences of the CDR3 and V genes to cluster similar TCR/BCRs together.
As a default, the function takes the input from \code{\link[=combineTCR]{combineTCR()}},
\code{\link[=combineBCR]{combineBCR()}} or \code{\link[=combineExpression]{combineExpression()}} and amends a
cluster to the data frame or meta data. If \strong{exportGraph} is set
to TRUE, the function returns an igraph object of the connected sequences.
If multiple sequences per chain are present, this function only compares
the first sequence.
}
\examples{
# Getting the combined contigs
combined <- combineTCR(contig_list, 
                        samples = c("P17B", "P17L", "P18B", "P18L", 
                                    "P19B","P19L", "P20B", "P20L"))

sub_combined <- clonalCluster(combined[c(1,2)], 
                              chain = "TRA", 
                              sequence = "aa")

}
\concept{Visualizing_Clones}
