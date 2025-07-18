% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clonalSizeDistribution.R
\name{clonalSizeDistribution}
\alias{clonalSizeDistribution}
\title{Plot powerTCR Clustering Based on Clonal Size}
\usage{
clonalSizeDistribution(
  input.data,
  cloneCall = "strict",
  chain = "both",
  method = "ward.D2",
  threshold = 1,
  group.by = NULL,
  exportTable = FALSE,
  palette = "inferno",
  ...
)
}
\arguments{
\item{input.data}{The product of \code{\link[=combineTCR]{combineTCR()}},
\code{\link[=combineBCR]{combineBCR()}}, or \code{\link[=combineExpression]{combineExpression()}}.}

\item{cloneCall}{Defines the clonal sequence grouping. Accepted values
are: \code{gene} (VDJC genes), \code{nt} (CDR3 nucleotide sequence), \code{aa} (CDR3 amino
acid sequence), or \code{strict} (VDJC + nt). A custom column header can also be used.}

\item{chain}{The TCR/BCR chain to use. Use \code{both} to include both chains
(e.g., TRA/TRB). Accepted values: \code{TRA}, \code{TRB}, \code{TRG}, \code{TRD}, \code{IGH}, \code{IGL}
(for both light chains), \code{both}.}

\item{method}{The clustering parameter for the dendrogram.}

\item{threshold}{Numerical vector containing the thresholds
the grid search was performed over.}

\item{group.by}{A column header in the metadata or lists to group the analysis
by (e.g., "sample", "treatment"). If \code{NULL}, data will be analyzed as
by list element or active identity in the case of single-cell objects.}

\item{exportTable}{If \code{TRUE}, returns a data frame or matrix of the results
instead of a plot.}

\item{palette}{Colors to use in visualization - input any
\link[grDevices:palettes]{hcl.pals}.}

\item{...}{Additional arguments passed to the ggplot theme}
}
\value{
A ggplot object visualizing dendrogram of clonal size distribution
or a data.frame if \code{exportTable = TRUE}.
}
\description{
This function produces a hierarchical clustering of clones by sample
using discrete gamma-GPD spliced threshold model. If using this
model please read and cite powerTCR (more info available at
\href{https://pubmed.ncbi.nlm.nih.gov/30485278/}{PMID: 30485278}).
}
\details{
The probability density function (pdf) for the \strong{Generalized Pareto Distribution (GPD)} is given by:
\deqn{f(x|\mu, \sigma, \xi) = \frac{1}{\sigma} \left( 1 + \xi \left( \frac{x - \mu}{\sigma} \right) \right)^{-\left( \frac{1}{\xi} + 1 \right)}}

Where:
\itemize{
\item{\eqn{\mu} is a location parameter}
\item{\eqn{\sigma > 0} is a scale parameter}
\item{\eqn{\xi} is a shape parameter}
\item{\eqn{x \ge \mu} if \eqn{\xi \ge 0} and \eqn{\mu \le x \le \mu - \sigma/\xi} if \eqn{\xi < 0}}
}

The probability density function (pdf) for the \strong{Gamma Distribution} is given by:
\deqn{f(x|\alpha, \beta) = \frac{x^{\alpha-1} e^{-x/\beta}}{\beta^\alpha \Gamma(\alpha)}}

Where:
\itemize{
\item{\eqn{\alpha > 0} is the shape parameter}
\item{\eqn{\beta > 0} is the scale parameter}
\item{\eqn{x \ge 0}}
\item{\eqn{\Gamma(\alpha)} is the gamma function of \eqn{\alpha}}
}
}
\examples{
# Making combined contig data
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B","P19L", "P20B", "P20L"))

# Using clonalSizeDistribution()
clonalSizeDistribution(combined, 
                       cloneCall = "strict", 
                       method="ward.D2")

}
\author{
Hillary Koch
}
\concept{Visualizing_Clones}
