#' Demonstrate the relative abundance of clones by group or sample
#'
#' Displays the number of clones at specific frequencies by sample 
#' or group. Visualization can either be a line graph (
#' **scale** = FALSE) using calculated numbers or density 
#' plot (**scale** = TRUE). Multiple sequencing runs can 
#' be group together using the group parameter. If a matrix 
#' output for the data is preferred, set 
#' **exportTable** = TRUE.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' clonalAbundance(combined, 
#'                 cloneCall = "gene", 
#'                 scale = FALSE)
#'
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()], or [combineExpression()].
#' @param cloneCall How to call the clone - VDJC gene (**gene**), 
#' CDR3 nucleotide (**nt**), CDR3 amino acid (**aa**),
#' VDJC gene + CDR3 nucleotide (**strict**) or a custom variable 
#' in the data. 
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The variable to use for grouping
#' @param order.by A vector of specific plotting order or "alphanumeric"
#' to plot groups in order
#' @param scale Converts the graphs into density plots in order to show 
#' relative distributions.
#' @param exportTable Returns the data frame used for forming the graph
#' to the visualization.
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals].
#' @importFrom ggplot2 ggplot
#' @importFrom dplyr distinct
#' @export
#' @concept Visualizing_Clones
#' @author Nick Borcherding, Justin Reimertz
#' @return ggplot of the total or relative abundance of clones 
#' across quanta
clonalAbundance <- function(input.data, 
                            cloneCall = "strict", 
                            chain = "both", 
                            scale=FALSE, 
                            group.by = NULL, 
                            order.by = NULL,
                            exportTable = FALSE,
                            palette = "inferno") {
  Con.df <- NULL
  xlab <- "Abundance"
  input.data <- .data.wrangle(input.data, 
                              group.by, 
                              .theCall(input.data, cloneCall, check.df = FALSE), 
                              chain)
  cloneCall <- .theCall(input.data, cloneCall)
  
  names <- names(input.data)
  if (!is.null(group.by)) {
  	Con.df <- bind_rows(lapply(seq_along(input.data), function(x) {
  		input.data %>%
  			.parseContigs(x, names, cloneCall) %>%
  			mutate("{group.by}" := input.data[[x]][[group.by]])
  	})) %>% 
  		dplyr::distinct() %>%
  	  as.data.frame()
  	col <- length(unique(Con.df[[group.by]]))
    fill <- group.by
    if(!is.null(order.by)) {
        Con.df <- .ordering.function(vector = order.by,
                                     group.by = group.by, 
                                     data.frame =  Con.df)
    }
    if (scale == TRUE) { 
        ylab <- "Density of Clones"
        plot <- ggplot(Con.df, aes(x=Abundance, 
                                   fill=Con.df[,group.by])) +
                      geom_density(aes(y=after_stat(scaled)), 
                                   alpha=0.5, 
                                   lwd=0.25, 
                                   color="black", 
                                   bw=0.5)  +
                      scale_fill_manual(values = .colorizer(palette,col)) +
                      labs(fill = fill)
    } else { 
        ylab <- "Number of Clones"
        plot <- ggplot(Con.df, aes(x=Abundance, group.by = values, 
                               color = Con.df[,group.by])) +
                        geom_line(stat="count") +
                        scale_color_manual(values = .colorizer(palette,col)) +
                        labs(color = fill)
    }
  } else {
    for (i in seq_along(input.data)) {
      data1 <- .parseContigs(input.data, i, names, cloneCall) 
      Con.df<- rbind.data.frame(Con.df, data1) %>%
                  dplyr::distinct() %>%
                  as.data.frame()
    }
    Con.df <- data.frame(Con.df)
    if(!is.null(order.by)) {
      Con.df <- .ordering.function(vector = order.by,
                                   group.by = "values", 
                                   data.frame = Con.df)
    }
    
    col <- length(unique(Con.df$values))
    fill <- "Samples"
    if (scale == TRUE) { 
      ylab <- "Density of Clones"
      plot <- ggplot(Con.df, aes(Abundance, fill=values)) +
                      geom_density(aes(y=after_stat(scaled)), 
                                   alpha=0.5, 
                                   lwd=0.25, 
                                   color="black", 
                                   bw=0.5) +
                      scale_fill_manual(values = .colorizer(palette,col)) +
                      labs(fill = fill)
    } else { 
      ylab <- "Number of Clones"
      plot <- ggplot(Con.df, aes(x=Abundance, group = values, 
                               color = values)) +
                      geom_line(stat="count") +
                      scale_color_manual(values = .colorizer(palette,col)) +
                      labs(color = fill)
  } }
  if (exportTable == TRUE) { 
    return(Con.df) 
  }
  plot <- plot + scale_x_log10() + ylab(ylab) + xlab(xlab) +
    theme_classic()
  return(plot) 
}
