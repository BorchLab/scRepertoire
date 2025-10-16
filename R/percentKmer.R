#' Analyze K-mer Motif Composition
#' 
#' @description This function calculates and visualizes the frequency of k-mer 
#' motifs for either nucleotide (nt) or amino acid (aa) sequences. It produces 
#' a heatmap showing the relative composition of the most variable motifs across
#'  samples or groups.
#'
#' @details The function first calculates k-mer frequencies for each sample/group.
#' By default, it then identifies the 30 most variable motifs based on the Median
#' Absolute Deviation (MAD) across all samples and displays their frequencies 
#' in a heatmap.
#'
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' # Using percentKmer()
#' percentKmer(combined, 
#'             chain = "TRB", 
#'             motif.length = 3)
#' 
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()], or [combineExpression()]
#' @param chain The TCR/BCR chain to use. Accepted values: `TRA`, `TRB`, `TRG`,
#'  `TRD`, `IGH`, or `IGL` (for both light chains).
#' @param cloneCall Defines the clonal sequence grouping. Accepted values 
#' are: `nt` (CDR3 nucleotide sequence) or `aa` (CDR3 amino acid sequence).
#' @param group.by A column header in the metadata or lists to group the analysis 
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed as 
#' by list element or active identity in the case of single-cell objects.
#' @param order.by A character vector defining the desired order of elements 
#' of the `group.by` variable. Alternatively, use `alphanumeric` to sort groups 
#' automatically.
#' @param motif.length The length of the kmer to analyze
#' @param min.depth Minimum count a motif must reach to be retained in the 
#' output (`>= 1`). **Default:** `3`.
#' @param top.motifs Return the n most variable motifs as a function of 
#' median absolute deviation
#' @param exportTable If `TRUE`, returns a data frame or matrix of the results 
#' instead of a plot.
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals]
#' @param ... Additional arguments passed to the ggplot theme
#' 
#' @importFrom stats mad
#' @importFrom immApex calculateMotif
#' 
#' @export
#' 
#' @concept Summarize_Repertoire
#' 
#' @return A ggplot object displaying a heatmap of motif percentages.
#' If `exportTable = TRUE`, a matrix of the raw data is returned.
percentKmer <- function(input.data, 
                        chain = "TRB", 
                        cloneCall = "aa",
                        group.by = NULL, 
                        order.by = NULL,
                        motif.length = 3,
                        min.depth = 3,
                        top.motifs = 30,
                        exportTable = FALSE, 
                        palette = "inferno",
                        ...) {
  
  if(!cloneCall %in% c("aa", "nt")) {
    stop("Please select either nucleotide (nt) or amino acid (aa) sequences for cloneCall")
  }
  motifs.to.save <- NULL
  sco <- .is.seurat.or.se.object(input.data)
  input.data <- .dataWrangle(input.data, 
                             group.by, 
                             .theCall(input.data, cloneCall, 
                                      check.df = FALSE, silent = TRUE), 
                             chain)
  cloneCall <- .theCall(input.data, cloneCall)
  if(!is.null(group.by) && !sco) {
    input.data <- .groupList(input.data, group.by)
  }
  
  # Generating all motif counts
  motif.list <- lapply(input.data, function(x) {
      immApex::calculateMotif(x[[cloneCall]], 
                          motif.lengths = motif.length, 
                          min.depth = min.depth)
  })
  
  # Collecting motifs into a matrix
  unique.motifs <- unique(unlist(lapply(motif.list, `[`, , 1)))
  if(any(grepl(";", unique.motifs))) {
    unique.motifs <- unique.motifs[-grep(";", unique.motifs)] #remove seq deliminator
  }
  mat <- matrix(0, ncol = length(unique.motifs), nrow = length(input.data))
  colnames(mat) <- unique.motifs
  rownames(mat) <- names(input.data)
  
  # Populating matrix
  for (i in seq_along(motif.list)) {
    sample_motifs <- motif.list[[i]]
    col_indices <- match(sample_motifs[,1], colnames(mat))
    mat[i, na.omit(col_indices)] <- sample_motifs[!is.na(col_indices), 2]
    mat[i, ] <- mat[i, ] / sum(mat[i,])
  }
  
  # Filtering for top.motifs
  if(!is.null(top.motifs)) {
    mads <- apply(mat, 2, mad)
    motifs.to.save <- names(sort(mads, decreasing = TRUE))[seq_len(top.motifs)]
    mat <- mat[, colnames(mat) %in% motifs.to.save, drop = FALSE]
  }
  
  # Export table if asked for
  if (exportTable) {
    return(mat)
  }

  # Getting mat into a ggplot-compliant form
  mat_melt <- expand.grid(Var1 = rownames(mat), Var2 = colnames(mat), KEEP.OUT.ATTRS = FALSE)
  mat_melt$value <- as.vector(mat)
  if (!is.null(motifs.to.save)) {
    mat_melt$Var2 <- factor(mat_melt$Var2, levels = rev(motifs.to.save))
  }
  
  # Ordering group is order.by is not NULL
  if(!is.null(order.by)) {
    mat_melt <- .orderingFunction(vector = order.by,
                                  group.by = "Var1", 
                                  mat_melt)
  }
  
  #Plotting
  plot <- ggplot(mat_melt, aes(x=.data[["Var2"]], y = .data[["Var1"]], fill=.data[["value"]])) +
            geom_tile(lwd = 0.1, color = "black") + 
            scale_fill_gradientn(name = "Proportion", colors = .colorizer(palette,21)) +
            .themeRepertoire(...) + 
            coord_flip() + 
            theme(axis.title = element_blank())
  return(plot)
}
