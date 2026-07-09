#' Plot Clonal Homeostasis of the Repertoire
#'
#' This function calculates the space occupied by clone proportions.
#' The grouping of these clones is based on the parameter `clone.size`,
#' at default, `clone.size` will group the clones into bins of Rare = 0
#' to 0.0001, Small = 0.0001 to 0.001, etc. To adjust the proportions,
#' change the number or labeling of the clone.size parameter. If a matrix
#' output for the data is preferred, set `export.table` = TRUE.
#'
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list,
#'                         samples = c("P17B", "P17L", "P18B", "P18L",
#'                                     "P19B","P19L", "P20B", "P20L"))
#' clonalHomeostasis(combined, clone.call = "gene")
#'
#' @param input.data The product of [combineTCR()], [combineBCR()], or
#' [combineExpression()].
#' @param clone.size The cut points of the proportions.
#' @param clone.call Defines the clonal sequence grouping. Accepted values
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL`,
#' `IGK`, `Light` (for both light chains), or `both` (for TRA/B and Heavy/Light).
#' @param group.by A column header in the metadata or lists to group the analysis
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed
#' by list element or active identity in the case of single-cell objects.
#' @param order.by A character vector defining the desired order of elements
#' of the `group.by` variable. Alternatively, use `alphanumeric` to sort groups
#' automatically.
#' @param export.table If `TRUE`, returns a data frame or matrix of the results
#' instead of a plot.
#' @param palette Colors to use in visualization - input any
#' [hcl.pals][grDevices::hcl.pals].
#' @param cloneSize \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `clone.size` instead.
#' @param cloneCall \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `clone.call` instead.
#' @param exportTable \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `export.table` instead.
#' @param ... Additional arguments passed to the ggplot theme
#'
#' @export
#' @importFrom utils stack
#' @importFrom stats reshape
#' @concept Visualizing_Clones
#' @return A ggplot object visualizing clonal homeostasis, or a data.frame if
#' `export.table = TRUE`.
clonalHomeostasis <- function(input.data,
                              clone.size = NULL,
                              clone.call = NULL,
                              chain = "both",
                              group.by = NULL,
                              order.by = NULL,
                              export.table = NULL,
                              palette = "inferno",
                              # Deprecated arguments
                              cloneSize = NULL,
                              cloneCall = NULL,
                              exportTable = NULL,
                              ...) {

    # Handle deprecated arguments
    clone.size <- .deprecate_arg(cloneSize, clone.size, "cloneSize", "clone.size",
                                 "clonalHomeostasis",
                                 default = c(Rare = .0001, Small = .001, Medium = .01, Large = .1, Hyperexpanded = 1))
    clone.call <- .deprecate_arg(cloneCall, clone.call, "cloneCall", "clone.call",
                                 "clonalHomeostasis", default = "strict")
    export.table <- .deprecate_arg(exportTable, export.table, "exportTable", "export.table",
                                   "clonalHomeostasis", default = FALSE)

    clone.size <- c(None = 0, clone.size)
    input.data <- .dataWrangle(input.data,
                               group.by,
                               .theCall(input.data, clone.call,
                                        check.df = FALSE, silent = TRUE),
                               chain)
    clone.call <- .theCall(input.data, clone.call)
    sco <- .is.seurat.or.se.object(input.data)
    if(!is.null(group.by) & !sco) {
      input.data <- .groupList(input.data, group.by)
    }
    
    #Generating data matrix to store value
    mat <- matrix(0, length(input.data), length(clone.size) - 1,
                dimnames = list(names(input.data),
                names(clone.size)[-1]))

    #Assigning the clonal grouping
    input.data <- lapply(input.data, '[[', clone.call)
    input.data <- lapply(input.data, na.omit)
    fun <- function(x) { table(x)/length(x) }
    input.data <- lapply(input.data, fun)
    for (i in 2:length(clone.size)) {
        mat[,i-1] <- vapply(input.data, function (x)
                              sum(x[x > clone.size[i-1] & x <= clone.size[i]]), FUN.VALUE = numeric(1))
        colnames(mat)[i-1] <- paste0(names(clone.size[i]), ' (',
                                    clone.size[i-1], ' < X <= ',
                                    clone.size[i], ')') }
    if (export.table) {
      return(mat)
    }
    
    #Plotting
    mat_df <- as.data.frame(mat) 
    mat_df$Var1 <- rownames(mat_df)
    
    varying_cols <- names(mat_df)[grep("<", names(mat_df))]
    
    mat_melt <- reshape(mat_df,
                        varying = varying_cols,
                        v.names = "value",
                        timevar = "category",
                        times = varying_cols,
                        idvar = "Var1",
                        direction = "long")

    # Order the clonal groups by their proportion bins (largest to smallest)
    # rather than alphabetically. varying_cols follows the clone.size order
    # (ascending), so reversing gives Hyperexpanded -> Rare.
    mat_melt$category <- factor(mat_melt$category, levels = rev(varying_cols))

    if(!is.null(order.by)) {
      mat_melt <- .orderingFunction(vector = order.by,
                                    group.by = "Var1", 
                                    data.frame = mat_melt)
    }
    
    col <- length(unique(mat_melt$category))
    plot <- ggplot(mat_melt, aes(x=as.factor(.data[["Var1"]]), y=.data[["value"]], fill=.data[["category"]])) +
        geom_bar(stat = "identity", position="fill", 
                    color = "black", lwd= 0.25) +
        scale_fill_manual(values = rev(.colorizer(palette,col))) +
        labs(x = "Samples", 
             y = "Relative Abundance", 
             fill = "Clonal Group") +
        .themeRepertoire(...)
    return(plot)
}
