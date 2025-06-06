% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clonalProportion.R
\name{clonalProportion}
\alias{clonalProportion}
\title{Examining the clonal space occupied by specific clones}
\usage{
clonalProportion(
  input.data,
  clonalSplit = c(10, 100, 1000, 10000, 30000, 1e+05),
  cloneCall = "strict",
  chain = "both",
  group.by = NULL,
  order.by = NULL,
  exportTable = FALSE,
  palette = "inferno"
)
}
\arguments{
\item{input.data}{The product of \code{\link[=combineTCR]{combineTCR()}},
\code{\link[=combineBCR]{combineBCR()}}, or \code{\link[=combineExpression]{combineExpression()}}.}

\item{clonalSplit}{The cut points for the specific clones}

\item{cloneCall}{How to call the clone - VDJC gene (\strong{gene}),
CDR3 nucleotide (\strong{nt}), CDR3 amino acid (\strong{aa}),
VDJC gene + CDR3 nucleotide (\strong{strict}) or a custom variable
in the data}

\item{chain}{indicate if both or a specific chain should be used -
e.g. "both", "TRA", "TRG", "IGH", "IGL"}

\item{group.by}{The variable to use for grouping}

\item{order.by}{A vector of specific plotting order or "alphanumeric"
to plot groups in order}

\item{exportTable}{Exports a table of the data into the global.
environment in addition to the visualization}

\item{palette}{Colors to use in visualization - input any
\link[grDevices:palettes]{hcl.pals}}
}
\value{
ggplot of the space occupied by the specific rank of clones
}
\description{
This function calculates the relative clonal space occupied by the
clones. The grouping of these clones is based on the parameter
\strong{clonalSplit}, at default, \strong{clonalSplit} will group the clones
into bins of 1:10, 11:100, 101:1001, etc. To adjust the clones
selected, change the numbers in the variable split. If a matrix output
for the data is preferred, set \strong{exportTable} = TRUE.
}
\examples{
#Making combined contig data
combined <- combineTCR(contig_list, 
                        samples = c("P17B", "P17L", "P18B", "P18L", 
                                    "P19B","P19L", "P20B", "P20L"))
clonalProportion(combined, cloneCall = "gene")

}
\concept{Visualizing_Clones}
