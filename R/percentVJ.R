#' Quantifying the V and J gene usage across clones
#'
#' This function the proportion V and J genes used by 
#' grouping variables for an indicated **chain** to
#' produce a matrix of VJ gene pairings.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' percentVJ(combined, chain = "TRB")
#' 
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()], or [combineExpression()].
#' @param chain "TRA", "TRB", "TRG", "TRG", "IGH", "IGL"
#' @param group.by The variable to use for grouping
#' @param order.by A vector of specific plotting order or "alphanumeric"
#' to plot groups in order
#' @param exportTable Returns the data frame used for forming the graph
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals].
#' @export
#' @concept Summarize_Repertoire
#' @return ggplot of percentage of V and J gene pairings as a heatmap
#' 
percentVJ <- function(input.data,
                      chain = "TRB",
                      group.by = NULL, 
                      order.by = NULL,
                      exportTable = FALSE, 
                      palette = "inferno") {
  
  sco <- is_seurat_object(input.data) | is_se_object(input.data)
  input.data <- .data.wrangle(input.data, group.by, "CTgene", chain)
  if(!is.null(group.by) & !sco) {
    input.data <- .groupList(input.data, group.by)
  }
  
  if(chain %in% c("TRA", "TRG", "IGL")) {
    positions <- c(1,2)
  } else {
    positions <- c(1,3)
  }
  
  #Getting indicated genes across list
  gene_counts <- lapply(input.data, function(x) {
      tmp <- unlist(str_split(x[,"CTgene"], ";"))
      tmp <- str_split(tmp, "[.]", simplify = TRUE)
      strings <- paste0(tmp[,positions[1]], ";", tmp[,positions[2]])
      strings <- strings[-which(strings == "NA;")]
  })
  #Need total unique genes
  gene.dictionary <- unique(unlist(gene_counts))
  coordinates <- strsplit(gene.dictionary, ";")
  V_values <- unique(unlist(lapply(coordinates, function(coord) coord[1])))
  J_values <- unique(unlist(lapply(coordinates, function(coord) coord[2])))
  
  #Summarizing the gene usage across the list
  summary <- lapply(gene_counts, function(x) {
                 percentages <- unlist(prop.table(table(x)))
                 genes.to.add <- gene.dictionary [which(gene.dictionary  %!in% names(percentages))]
                 if(length(genes.to.add) >= 1) {
                   percentages.to.add <- rep(0, length(genes.to.add))
                   names(percentages.to.add) <- genes.to.add
                   percentages <- c(percentages, percentages.to.add)
                 }
                 percentages <- percentages[match(str_sort(names(percentages), numeric = TRUE), names(percentages))]
  })
  if (exportTable == TRUE) { 
    summary.matrix <- do.call(rbind,summary)
    return(summary.matrix) 
  }
  mat <- lapply(summary, function(x) {
    # Create an empty matrix
    result_matrix <- matrix(0, nrow = length(V_values), ncol = length(unique(J_values)))
    rownames(result_matrix) <- V_values
    colnames(result_matrix) <- J_values
    
    for (i in seq_len(length(x))) {
      coordinates <- unlist(strsplit(names(x)[i], ";"))
        result_matrix[coordinates[1], coordinates[2]] <- x[i]
    }
    result_matrix
  })
  
  #Melting matrix and Visualizing
  mat_melt <- melt(mat)
  if(!is.null(order.by)) {
    mat_melt <- .ordering.function(vector = order.by,
                                   group.by = "L1", 
                                   mat_melt)
  }
  
  plot <- ggplot(mat_melt, aes(y=Var1, x = Var2, fill=round(value*100,2))) +
    geom_tile(lwd= 0.25, color = "black") +
    scale_fill_gradientn(name = "Percentage", colors = .colorizer(palette,21)) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          axis.title = element_blank())
  if(length(unique(mat_melt$L1)) > 1) {
    plot <- plot + facet_wrap(~L1) 
  }
  return(plot)
}
