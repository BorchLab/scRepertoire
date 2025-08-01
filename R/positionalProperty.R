#' Plot Positional Physicochemical Property Analysis
#'
#' This function analyzes the physicochemical properties of amino acids at each
#' position along the CDR3 sequence. It calculates the mean property value and
#' the 95% confidence interval for each position across one or more groups,
#' visualizing the results as a line plot with a confidence ribbon.
#' 
#' @details
#' The function uses one of several established physicochemical property scales
#' to convert amino acid sequences into numerical vectors. More information for
#' the individual methods can be found at the following citations:
#'
#' **atchleyFactors:** [citation](https://pubmed.ncbi.nlm.nih.gov/15851683/)
#'
#' **crucianiProperties:** [citation](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/abs/10.1002/cem.856)
#'
#' **FASGAI:** [citation](https://pubmed.ncbi.nlm.nih.gov/18318694/)
#'
#' **kideraFactors:** [citation](https://link.springer.com/article/10.1007/BF01025492)
#'
#' **MSWHIM:** [citation](https://pubs.acs.org/doi/10.1021/ci980211b)
#'
#' **ProtFP:** [citation](https://pubmed.ncbi.nlm.nih.gov/24059694/)
#'
#' **stScales:** [citation](https://pubmed.ncbi.nlm.nih.gov/19373543/)
#'
#' **tScales:** [citation](https://www.sciencedirect.com/science/article/abs/pii/S0022286006006314)
#'
#' **VHSE:** [citation](https://pubmed.ncbi.nlm.nih.gov/15895431/)
#'
#' **zScales:** [citation](https://pubmed.ncbi.nlm.nih.gov/9651153/)
#' 
#'
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' # Using positionalProperty()
#' positionalProperty(combined, 
#'                    chain = "TRB",
#'                    method = "atchleyFactors", 
#'                    aa.length = 20)
#'                    
#' @param input.data The product of [combineTCR()], [combineBCR()], or 
#' [combineExpression()]
#' @param chain The TCR/BCR chain to use. Use Accepted values: `TRA`, `TRB`, 
#' `TRG`, `TRD`, `IGH`, or `IGL` (for both light chains)..
#' @param group.by A column header in the metadata or lists to group the analysis 
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed 
#' by list element or active identity in the case of single-cell objects.
#' @param order.by A character vector defining the desired order of elements 
#' of the `group.by` variable. Alternatively, use `alphanumeric` to sort groups 
#' automatically.
#' @param aa.length The maximum length of the CDR3 amino acid sequence. 
#' @param method Character string (one of the supported names) 
#' Defaults to `"atchleyFactors"`, but includes: `"crucianiProperties"`, 
#' `"FASGAI"`, `"kideraFactors"`, `"MSWHIM"`, `"ProtFP"`, `"stScales"`, 
#' `"tScales"`, `"VHSE"`, `"zScales"`
#' @param exportTable If `TRUE`, returns a data frame or matrix of the results 
#' instead of a plot.
#' @param palette Colors to use in visualization - input any [hcl.pals][grDevices::hcl.pals]
#' @param ... Additional arguments passed to the ggplot theme
#' @importFrom stats qt sd
#' @importFrom utils getFromNamespace
#' @export
#' @concept Summarize_Repertoire
#' @return A ggplot object displaying property by amino acid position.
#' If `exportTable = TRUE`, a matrix of the raw data is returned.
#' @author Florian Bach, Nick Borcherding

positionalProperty <- function(input.data, 
                               chain = "TRB", 
                               group.by = NULL, 
                               order.by = NULL,
                               aa.length = 20,
                               method = "atchleyFactors",
                               exportTable = FALSE, 
                               palette = "inferno",
                               ...)  {
  factors <- c("atchleyFactors", "crucianiProperties", "FASGAI", "kideraFactors", 
  "MSWHIM", "ProtFP", "stScales", "tScales", "VHSE", "zScales")
  if (!method %in% c(factors)) {
    stop("Please select a compatible method: ", paste0(factors, collapse = ", "))
  }
  
  amino.acids <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", 
                   "M", "F", "P", "S", "T", "W", "Y", "V")
  
  sco <- .is.seurat.or.se.object(input.data)
  input.data <- .dataWrangle(input.data, 
                             group.by, 
                             .theCall(input.data, "CTaa", 
                                      check.df = FALSE, silent = TRUE), 
                             chain)
  cloneCall <- .theCall(input.data, "CTaa")
  
  .aa.property.matrix <- getFromNamespace(".aa.property.matrix", "immApex")
  .padded.strings <- getFromNamespace(".padded.strings", "immApex")
  
  if(!is.null(group.by) & !sco) {
    input.data <- .groupList(input.data, group.by)
  }
  
  # Load the selected property matrix
  S <- .aa.property.matrix(method)
  S <- S[, amino.acids, drop = FALSE]
  
  # Main calculation loop over groups
  results_list <- lapply(names(input.data), function(current_group_name) {
    x <- input.data[[current_group_name]]
    input.sequences <- x[["CTaa"]]
    
    # Clean sequence vector
    if (any(grepl(";", input.sequences))) {
      input.sequences <- unlist(strsplit(input.sequences, ";"))
    }
    input.sequences <- na.omit(input.sequences)
    input.sequences <- input.sequences[nchar(input.sequences) > 0]
    
    if (length(input.sequences) == 0) return(NULL)
    
    n_seq <- length(input.sequences)
    k <- nrow(S)
    
    # Vectorized profile creation
    padded_seqs <- .padded.strings(input.sequences, 
                                   max.length =  aa.length, 
                                   pad = ".")
    char_matrix <- do.call(rbind, strsplit(padded_seqs, ""))
    index_matrix <- matrix(match(char_matrix, colnames(S)), nrow = n_seq)
    profile_array <- array(NA_real_, dim = c(k, aa.length, n_seq), dimnames = list(rownames(S), 1:aa.length, NULL))
    
    for (i in 1:k) {
      prop_values_matrix <- matrix(S[i, index_matrix], nrow = n_seq)
      profile_array[i, , ] <- t(prop_values_matrix)
    }
    
    # Summarization
    mean_vals <- apply(profile_array, c(1, 2), mean, na.rm = TRUE)
    sd_vals   <- apply(profile_array, c(1, 2), sd, na.rm = TRUE)
    n_obs     <- apply(profile_array, c(1, 2), function(val) sum(!is.na(val)))
    
    # 95% CI Calculation
    conf.level <- 0.95
    t_critical <- qt(1 - (1 - conf.level) / 2, df = pmax(1, n_obs - 1))
    margin_of_error <- t_critical * (sd_vals / sqrt(n_obs))
    
    # Construct summary data.frame
    summary_df <- as.data.frame.table(mean_vals, responseName = "mean")
    names(summary_df) <- c("property", "position", "mean")
    
    summary_df$sd <- as.vector(sd_vals)
    summary_df$n <- as.vector(n_obs)
    summary_df$ci_lower <- as.vector(mean_vals - margin_of_error)
    summary_df$ci_upper <- as.vector(mean_vals + margin_of_error)
    summary_df$group <- current_group_name
    
    return(summary_df)
  })
  
  # Combine results and handle ordering
  mat <- do.call(rbind, results_list[!sapply(results_list, is.null)])
  if (is.null(mat) || nrow(mat) == 0) {
    message("No sequences found to plot after filtering.")
    return(NULL)
  }
  #Removing any CI for calculations based on < 2 residues
  cols_to_modify <- c("sd", "ci_lower", "ci_upper")
  mat[mat[["n"]] <= 2, cols_to_modify] <- 0
  
  if (!is.null(order.by)) {
    mat <- .orderingFunction(mat, order.by, "group")
  }
  mat$position <- as.integer(mat$position)
  
  if (exportTable) {
    return(mat)
  }
  
  colors <- .colorizer(palette, length(unique(mat[["group"]])))
  # Visualization
  plot <- ggplot(mat, aes(x = as.factor(position), y = mean, group = group, color = group)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = group), alpha = 0.3, linetype = 0) +
    geom_line() +
    scale_color_manual(name = "Group", values = colors) +
    scale_fill_manual(name = "Group", values = colors) +
    labs(x = "Amino Acid Position", y = "Mean Property Value") +
    facet_grid(property~., scales = "free_y") +
    guides(fill = "none", 
           color = guide_legend(override.aes = list(fill = NA))) + 
    .themeRepertoire(...) + 
    theme(
      strip.background = element_rect(fill="white", color = "grey"),
      strip.text = element_text(face = "bold"),
      legend.key = element_rect(fill = NA)
    )
  
  return(plot)
}

