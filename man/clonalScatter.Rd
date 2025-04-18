% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clonalScatter.R
\name{clonalScatter}
\alias{clonalScatter}
\title{Scatter plot comparing the clonal expansion of two samples}
\usage{
clonalScatter(
  input.data,
  cloneCall = "strict",
  x.axis = NULL,
  y.axis = NULL,
  chain = "both",
  dot.size = "total",
  group.by = NULL,
  graph = "proportion",
  exportTable = FALSE,
  palette = "inferno"
)
}
\arguments{
\item{input.data}{The product of \code{\link[=combineTCR]{combineTCR()}},
\code{\link[=combineBCR]{combineBCR()}}, or \code{\link[=combineExpression]{combineExpression()}}.}

\item{cloneCall}{How to call the clone - VDJC gene (\strong{gene}),
CDR3 nucleotide (\strong{nt}), CDR3 amino acid (\strong{aa}),
VDJC gene + CDR3 nucleotide (\strong{strict}) or a custom variable
in the data.}

\item{x.axis}{name of the list element to appear on the x.axis.}

\item{y.axis}{name of the list element to appear on the y.axis.}

\item{chain}{indicate if both or a specific chain should be used -
e.g. "both", "TRA", "TRG", "IGH", "IGL".}

\item{dot.size}{either total or the name of the list element to
use for size of dots.}

\item{group.by}{The variable to use for grouping.}

\item{graph}{graph either the clonal "proportion" or "count".}

\item{exportTable}{Returns the data frame used for forming the graph.}

\item{palette}{Colors to use in visualization - input any
\link[grDevices:palettes]{hcl.pals}.}
}
\value{
ggplot of the relative clone numbers between two sequencing runs or groups
}
\description{
This function produces a scatter plot directly comparing
the specific clones between two samples. The clones will
be categorized by counts into singlets or expanded, either
exclusive or shared between the selected samples.
}
\examples{
#Making combined contig data
combined <- combineTCR(contig_list, 
                        samples = c("P17B", "P17L", "P18B", "P18L", 
                                    "P19B","P19L", "P20B", "P20L"))
clonalScatter(combined, 
              x.axis = "P17B", 
              y.axis = "P17L",
              graph = "proportion")

}
\concept{Visualizing_Clones}
