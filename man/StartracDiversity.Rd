% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/startracDiversity.R
\name{StartracDiversity}
\alias{StartracDiversity}
\title{Startrac-based diversity indices for single-cell RNA-seq}
\usage{
StartracDiversity(
  sc.data,
  cloneCall = "strict",
  chain = "both",
  type = NULL,
  group.by = NULL,
  exportTable = FALSE,
  palette = "inferno"
)
}
\arguments{
\item{sc.data}{The single-cell object after \code{\link[=combineExpression]{combineExpression()}}.
For SCE objects, the cluster variable must be in the meta data under
"cluster".}

\item{cloneCall}{How to call the clone - VDJC gene (\strong{gene}),
CDR3 nucleotide (\strong{nt}), CDR3 amino acid (\strong{aa}),
VDJC gene + CDR3 nucleotide (\strong{strict}) or a custom variable
in the data.}

\item{chain}{indicate if both or a specific chain should be used -
e.g. "both", "TRA", "TRG", "IGH", "IGL".}

\item{type}{The variable in the meta data that provides tissue type.}

\item{group.by}{The variable in the meta data to group by, often samples.}

\item{exportTable}{Returns the data frame used for forming the graph.}

\item{palette}{Colors to use in visualization - input any \link[grDevices:palettes]{hcl.pals}.}
}
\value{
ggplot object of Startrac diversity metrics
}
\description{
This function utilizes the Startrac approach derived from
\href{https://pubmed.ncbi.nlm.nih.gov/30479382/}{PMID: 30479382}.
Required to run the function, the "type" variable needs to include the
difference in where the cells were derived. The output of this function
will produce 3 indices: \strong{expa} (clonal expansion), \strong{migra}
(cross-tissue migration), and \strong{trans} (state transition). In order
to understand the underlying analyses of the outputs please
read and cite the linked manuscript.
}
\examples{
#Getting the combined contigs
combined <- combineTCR(contig_list, 
                        samples = c("P17B", "P17L", "P18B", "P18L", 
                                    "P19B","P19L", "P20B", "P20L"))

#Getting a sample of a Seurat object
scRep_example  <- get(data("scRep_example"))
scRep_example  <- combineExpression(combined, scRep_example)
scRep_example$Patient <- substring(scRep_example$orig.ident,1,3)
scRep_example$Type <- substring(scRep_example$orig.ident,4,4) 

#Using StartracDiversity()
StartracDiversity(scRep_example, 
                  type = "Type", 
                  group.by = "Patient")

}
\author{
Liangtao Zheng
}
\concept{SC_Functions}
