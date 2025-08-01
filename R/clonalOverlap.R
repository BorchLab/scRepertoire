#' Examining the clonal overlap between groups or samples
#'
#' This functions allows for the calculation and visualizations of 
#' various overlap metrics for clones. The methods include overlap 
#' coefficient (`overlap`), Morisita's overlap index 
#' (`morisita`), Jaccard index (`jaccard`), cosine 
#' similarity (`cosine`) or the exact number of clonal 
#' overlap (`raw`).
#' 
#' @details
#' The formulas for the indices are as follows:
#' 
#' **Overlap Coefficient:**
#' \deqn{overlap = \frac{\sum \min(a, b)}{\min(\sum a, \sum b)}}  
#' 
#' **Raw Count Overlap:**
#' \deqn{raw = \sum \min(a, b)}
#' 
#' **Morisita Index:**
#' \deqn{morisita = \frac{\sum a b}{(\sum a)(\sum b)}}  
#' 
#' **Jaccard Index:**
#' \deqn{jaccard = \frac{\sum \min(a, b)}{\sum a + \sum b - \sum \min(a, b)}}  
#' 
#' **Cosine Similarity:**
#' \deqn{cosine = \frac{\sum a b}{\sqrt{(\sum a^2)(\sum b^2)}}}  
#' 
#' Where:  
#' \itemize{  
#'   \item{\eqn{a} and \eqn{b} are the abundances of species \eqn{i} in groups A and B, respectively.}
#' }

#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' # Using clonalOverlap()
#' clonalOverlap(combined, 
#'               cloneCall = "aa", 
#'               method = "jaccard")
#'
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()], or [combineExpression()]
#' @param cloneCall Defines the clonal sequence grouping. Accepted values 
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino 
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains 
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL` 
#' (for both light chains), `both`.
#' @param method The method to calculate the `overlap`, `morisita`, 
#' `jaccard`, `cosine` indices or `raw` for the base numbers
#' @param group.by A column header in the metadata or lists to group the analysis 
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed 
#' by list element or active identity in the case of single-cell objects.
#' @param order.by A character vector defining the desired order of elements 
#' of the `group.by` variable. Alternatively, use `alphanumeric` to sort groups 
#' automatically.
#' @param exportTable If `TRUE`, returns a data frame or matrix of the results 
#' instead of a plot.
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals]
#' @param ... Additional arguments passed to the ggplot theme
#' 
#' @importFrom stats quantile
#' @export
#' @concept Visualizing_Clones
#' @return A ggplot object visualizing clonal overlap or a data.frame if
#'`exportTable = TRUE`.
clonalOverlap <- function(input.data, 
                          cloneCall = "strict", 
                          method = c("overlap", "morisita", "jaccard", "cosine", "raw"),
                          chain = "both", 
                          group.by = NULL,
                          order.by = NULL,
                          exportTable = FALSE,
                          palette = "inferno",
                          ...){
    method <- match.arg(method)
    if(method == "morisita") {
      return_type <- "freq"
    } else {
      return_type <- "unique"
    }
    input.data <- .dataWrangle(input.data, 
                               group.by, 
                               .theCall(input.data, cloneCall, 
                                        check.df = FALSE, silent = TRUE), 
                               chain)
    if(!is.null(order.by)) {
      if(length(order.by) == 1 && order.by == "alphanumeric") {
        input.data <- input.data[.alphanumericalSort(names(input.data))]
        
      } else {
        input.data <- input.data[order.by]
      }
    }
    
    cloneCall <- .theCall(input.data, cloneCall)
    
    sco <- .is.seurat.or.se.object(input.data)
    if(!is.null(group.by) & !sco) {
      input.data <- .groupList(input.data, group.by)
    }

    num_samples <- length(input.data[])
    names_samples <- names(input.data)
    length <- seq_len(num_samples)
    
    #Selecting Index Function
    indexFunc <- switch(method,
                        "morisita" = .morisitaCalc,
                        "jaccard"  = .jaccardCalc,
                        "raw"      = .rawCalc,
                        "overlap"  = .overlapCalc,
                        "cosine"  =  .cosineCalc,
                        stop("Invalid method provided"))
    
    #Calculating Index 
    coef_matrix <- data.frame(matrix(NA, num_samples, num_samples))
    coef_matrix <- .calculateIndex(input.data, 
                                   length, 
                                   cloneCall, 
                                   coef_matrix, 
                                   indexFunc, 
                                   return_type)
    
    #Data manipulation
    colnames(coef_matrix) <- names_samples
    rownames(coef_matrix) <- names_samples

    if (exportTable == TRUE) { 
      return(coef_matrix) 
    }
    
    coef_matrix_as_matrix <- as.matrix(coef_matrix)
    
    mat_melt <- data.frame(
      Var1 = rep(rownames(coef_matrix_as_matrix), ncol(coef_matrix_as_matrix)),
      Var2 = rep(colnames(coef_matrix_as_matrix), each = nrow(coef_matrix_as_matrix)),
      value = as.vector(coef_matrix_as_matrix)
    )
    
    mean_value <- mean(na.omit(mat_melt[,"value"]))
    
    plot <- ggplot(mat_melt, aes(x=.data[["Var1"]], y=.data[["Var2"]], fill=.data[["value"]])) +
                geom_tile() + 
                geom_tile(data = mat_melt[!is.na(mat_melt[,"value"]),], 
                          fill = NA, 
                          lwd= 0.25, 
                          color = "black") +
                labs(fill = .toCapitilize(method)) +
                geom_text(aes(label = round(value, digits = 3), 
                              color = ifelse(value <= mean_value,
                                             "white", "black")), 
                          na.rm = TRUE) +
                scale_fill_gradientn(colors = .colorizer(palette, 7), na.value = "white") +
                scale_color_identity() +
                .themeRepertoire() + 
                theme(axis.title = element_blank())
    return(plot) 
}

# Helper function to prepare data
.prepareDataFrame <- function(df, 
                              cloneCall, 
                              return_type = "unique") {
  if (return_type == "unique") {
    return(unique(df[, cloneCall]))
  } else if (return_type == "freq") {
    temp_df <- data.frame(table(df[, cloneCall]))
    colnames(temp_df) <- c(cloneCall, 'Count')
    temp_df[, 2] <- as.numeric(temp_df[, 2])
    return(temp_df)
  }
}

# Helper function for common loop and conditional structure
.calculateIndex <- function(df, 
                            length, 
                            cloneCall, 
                            coef_matrix, 
                            indexFunc, 
                            return_type = "unique") {
  for (i in seq_along(length)) {
    df_i <- .prepareDataFrame(df[[i]], cloneCall, return_type)
    for (j in seq_along(length)) {
      if (i >= j) { next }
      df_j <- .prepareDataFrame(df[[j]], cloneCall, return_type)
      coef_matrix[i, j] <- indexFunc(df_i, df_j)
    }
  }
  return(coef_matrix)
}

# Morisita Index calculation function
.morisitaCalc <- function(df_i, df_j) {
  merged <- merge(df_i, df_j, by = names(df_i)[1], all = TRUE)
  merged[is.na(merged)] <- 0
  
  X <- sum(merged[, 2])
  Y <- sum(merged[, 3])
  
  num <- 2 * sum(merged[, 2] * merged[, 3])
  den <- ((sum(df_i[, 2]^2) / (X^2)) + (sum(df_j[, 2]^2) / (Y^2))) * X * Y
  
  return(num / den)
}

# Jaccard Index calculation function
.jaccardCalc <- function(df_i, df_j) {
  overlap <- length(intersect(df_i, df_j))
  return(overlap / (length(df_i) + length(df_j) - overlap))
}

# Raw Index calculation function
.rawCalc <- function(df_i, df_j) {
  return(length(intersect(df_i, df_j)))
}

# Overlap Index calculation function
.overlapCalc <- function(df_i, df_j) {
  overlap <- length(intersect(df_i, df_j))
  return(overlap / min(length(df_i), length(df_j)))
}

# Overlap Index calculation function
.cosineCalc <- function(df_i, df_j) {
  all_species <- unique(c(df_i, df_j))
  vector_location1 <- as.integer(all_species %in% df_i)
  vector_location2 <- as.integer(all_species %in% df_j)
  
  return(sum(vector_location1 * vector_location2) / 
           (sqrt(sum(vector_location1^2)) * sqrt(sum(vector_location2^2))))
}
