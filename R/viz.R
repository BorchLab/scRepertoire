#' Quantify the unique clonotypes in the filtered contigs.
#'
#' This function takes the output from combineTCR(), combineBCR(), or 
#' expression2List() and quantifies unique clonotypes. The unique clonotypes 
#' can be either reported as a raw output or scaled to the total number of 
#' clonotypes recovered using the scale parameter. Multiple sequencing 
#' runs can be group together using the group parameter. If a matrix output 
#' for the data is preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' quantContig(combined, cloneCall="strict", scale = TRUE)
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The column header used for grouping.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param scale Converts the graphs into percentage of unique clonotypes.
#' @param exportTable Returns the data frame used for forming the graph
#' @import ggplot2
#' @export
#' @return ggplot of the total or relative unique clonotypes
quantContig <- function(df, cloneCall = "strict", chain = "both", 
                        scale=FALSE, group.by = NULL, 
                        split.by = NULL,order=TRUE,
                        exportTable = FALSE) {
    if (length(group.by) > 1) { stop("Only one item in the group.by variable can 
                                    be listed.") }
    df <- list.input.return(df, split.by)
    cloneCall <- theCall(cloneCall)
    df <- checkBlanks(df, cloneCall)
    if (chain != "both") {
      for(i in seq_along(df)) {
        df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
      }
    }
    if (!is.null(group.by)) {
        x <- group.by
        labs <- group.by
        Con.df <- data.frame(matrix(NA, length(df), 4))
        colnames(Con.df) <- c("contigs","values", "total", group.by)
        for (i in seq_along(df)) {
            Con.df[i,1] <- length(unique(df[[i]][,cloneCall]))
            Con.df[i,2] <- names(df)[i]
            Con.df[i,3] <- length(df[[i]][,cloneCall])
            location <- which(colnames(df[[i]]) == group.by)
            Con.df[i,4] <- df[[i]][1,location] }
        col <- length(unique(Con.df[,group.by]))
        if (scale == TRUE) { y <- "scaled"
            Con.df$scaled <- Con.df$contigs/Con.df$total*100
            ylab <- "Percent of Unique Clonotype"

        } else { y <- "contigs"
            x <- group.by
            ylab <- "Unique Clonotypes"}
    } else {
        x <- "values"
        labs <- "Samples"
        Con.df <- data.frame(matrix(NA, length(df), 3))
        colnames(Con.df) <- c("contigs","values", "total")
        for (i in seq_along(df)) {
            Con.df[i,1] <- length(unique(df[[i]][,cloneCall]))
            Con.df[i,2] <- names(df)[i]
            Con.df[i,3] <- length(df[[i]][,cloneCall]) }
        col <- length(unique(Con.df$values))
        if (scale == TRUE) { y <- "scaled"
            Con.df$scaled <- Con.df$contigs/Con.df$total*100
            ylab <- "Percent of Unique Clonotype"
        } else { y <- "contigs"
            ylab <- "Unique Clonotypes" } }
    if (exportTable == TRUE) { return(Con.df) }
    # set the x-axis order the same as input
    if (order == TRUE){
        order_x <- Con.df[,x] %>% fct_inorder()}
    else {
        order_x <- Con.df[,x]}
    plot <- ggplot(aes(x=order_x, y=Con.df[,y],
            fill=order_x), data = Con.df) +
        stat_summary(geom = "errorbar", fun.data = mean_se, 
            position = "dodge", width=.5) + labs(fill = labs) +
        stat_summary(fun=mean, geom="bar", color="black", lwd=0.25)+
        theme_classic() + xlab("Samples") + ylab(ylab) +
        scale_fill_manual(values = colorblind_vector(col))
    return(plot) }


#' Demonstrate the relative abundance of clonotypes by group or sample.
#'
#' This function takes the output of combineTCR(), combineBCR(), or 
#' expression2List() and displays the number of clonotypes at specific 
#' frequencies by sample or group. Visualization can either be a line 
#' graph using calculated numbers or if scale = TRUE, the output will 
#' be a density plot. Multiple sequencing runs can be group together 
#' using the group parameter. If a matrix output for the data is 
#' preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' abundanceContig(combined, cloneCall = "gene", scale = FALSE)
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The column header for which you would like to analyze the data.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param scale Converts the graphs into density plots in order to show 
#' relative distributions.
#' @param exportTable Returns the data frame used for forming the graph
#' to the visualization.
#' @importFrom ggplot2 ggplot
#' @export
#' @return ggplot of the total or relative abundance of clonotypes 
#' across quanta
abundanceContig <- function(df, cloneCall = "strict", chain = "both", 
                            scale=FALSE, group.by = NULL, 
                            split.by = NULL, exportTable = FALSE) {
    df <- list.input.return(df,split.by)
    Con.df <- NULL
    xlab <- "Abundance"
    cloneCall <- theCall(cloneCall)
    df <- checkBlanks(df, cloneCall)
    names <- names(df)
    if (!is.null(group.by)) {
        for (i in seq_along(df)) {
            if (chain != "both") {
              df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
            }
            data1 <- parseContigs(df, i, names, cloneCall)
            label <- df[[i]][1,group.by]
            data1[,paste(group.by)] <- label
            Con.df<- rbind.data.frame(Con.df, data1) }
        Con.df <- data.frame(Con.df)
        col <- length(unique(Con.df[,group.by]))
        fill <- group.by
        if (scale == TRUE) { ylab <- "Density of Clonotypes"
            plot <- ggplot(Con.df, aes(x=Abundance, fill=Con.df[,group.by])) +
                geom_density(aes(y=..scaled..), alpha=0.5, 
                    lwd=0.25, color="black", bw=0.5)  +
                scale_fill_manual(values = colorblind_vector(col)) +
                labs(fill = fill)
        } else { ylab <- "Number of Clonotypes"
            plot <- ggplot(Con.df, aes(x=Abundance, group.by = values, 
                    color = Con.df[,group.by])) +
                geom_line(stat="count") +
                scale_color_manual(values = colorblind_vector(col)) +
                labs(color = fill)}
    } else{
        for (i in seq_along(df)) {
          if (chain != "both") {
            df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
          }
            data1 <- parseContigs(df, i, names, cloneCall)
            Con.df<- rbind.data.frame(Con.df, data1) }
        col <- length(unique(Con.df$values))
        fill <- "Samples"
        if (scale == TRUE) { ylab <- "Density of Clonotypes"
            plot <- ggplot(Con.df, aes(Abundance, fill=values)) +
                geom_density(aes(y=..scaled..), alpha=0.5, lwd=0.25, 
                    color="black", bw=0.5) +
                scale_fill_manual(values = colorblind_vector(col)) +
                labs(fill = fill)
        } else { ylab <- "Number of Clonotypes"
            plot <- ggplot(Con.df, aes(x=Abundance, group = values, 
                    color = values)) +
                geom_line(stat="count") +
                scale_color_manual(values = colorblind_vector(col)) +
                labs(color = fill)} }
    if (exportTable == TRUE) { return(Con.df)}
    plot <- plot + scale_x_log10() + ylab(ylab) + xlab(xlab) +
        theme_classic()
return(plot) }

#' Demonstrate the distribution of lengths filtered contigs.
#'
#' This function takes the output of combineTCR(), combineBCR(), or
#' expression2List() and displays either the nucleotide (nt) or amino 
#' acid (aa) sequence length. The sequence length visualized can be 
#' selected using the chains parameter, either the combined clonotype 
#' (both chains) or across all single chains. Visualization can either 
#' be a histogram or if scale = TRUE, the output will be a density plot. 
#' Multiple sequencing runs can be group together using the 
#' group parameter. If a matrix output for the data is preferred, set 
#' exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' lengthContig(combined, cloneCall="aa", chain = "both")
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - CDR3 nucleotide (nt), 
#' CDR3 amino acid (aa).
#' @param group.by The group header for which you would like to analyze 
#' the data.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param scale Converts the graphs into density plots in order to show 
#' relative distributions.
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param exportTable Returns the data frame used for forming the graph.
#' @importFrom stringr str_split
#' @importFrom ggplot2 ggplot
#' @export
#' @return ggplot of the discrete or relative length distributions of 
#' clonotype sequences
lengthContig <- function(df, 
                         cloneCall = "aa", 
                         group.by = NULL, 
                         split.by = NULL, 
                         scale = FALSE, 
                         chain = "both", 
                         exportTable = FALSE) {
    df <- list.input.return(df, split.by)
    if(cloneCall == "nt") { cloneCall <- "CTnt"
        ylab <- "CDR3 (NT)"
    } else if (cloneCall == "aa") { cloneCall <- "CTaa" 
        ylab <- "CDR3 (AA)"
    } else { stop("Please make a selection of the type of
                CDR3 sequence to analyze by using `cloneCall`") }
    AB <- grep("TRA", df[[1]]$TCR1)
    IG <- grep("IGH", df[[1]]$TCR1)
    if (length(AB) > 1 && length(IG) == 0) {
      c1 <- "TRA"
      c2 <- "TRB"
    } else if (length(AB) > 1) {
      c1 <- "IGH"
      c2 <- "IGL"
    } else {
      c1 <- "TRD"
      c2 <- "TRG"
    }
    xlab <- "Length"
    Con.df <- NULL
    Con.df <- lengthDF(df, cloneCall, chain, group.by, c1, c2)
    names <- names(df)
    if (!is.null(group.by)) { 
        fill = group.by
        col <- length(unique(Con.df[,group.by]))
        if (scale == TRUE) { yplus <- "Percent of "
            plot <- ggplot(Con.df, aes(fill=Con.df[,group.by],
                length,(..count..)/sum(..count..)*100)) + 
                geom_density(aes(y=..scaled..),alpha=.5,lwd=.25,color="black")
        } else { yplus <- "Number of "
            plot<-ggplot(Con.df,aes(as.factor(length),fill=Con.df[,group.by]))+
                geom_bar(position = position_dodge2(preserve = "single"), 
                color="black", lwd=0.25, width=0.9)  +
                scale_x_discrete(breaks = round(seq(min(Con.df$length), 
                max(Con.df$length), by = 5),10)) }
    } else if (is.null(group.by)){ 
        fill <- "Samples"
        col <- length(unique(Con.df$values))
        if (scale == TRUE) { yplus <- "Percent of "
            plot <- ggplot(Con.df, aes(length, (..count..)/sum(..count..)*100, 
                fill=values)) + geom_density(aes(y=..scaled..), alpha=0.5, 
                lwd=0.25, color="black")
    }  else { yplus <- "Number of "
        plot <- ggplot(Con.df, aes(as.factor(length), fill=values)) +
            geom_bar(position = position_dodge2(preserve = "single"), 
            color="black", lwd=0.25) +
            scale_x_discrete(breaks = round(seq(min(Con.df$length), 
            max(Con.df$length), by = 5),10))} }
    plot <- plot + scale_fill_manual(values = colorblind_vector(col)) +
        labs(fill = fill) + ylab(paste(yplus, ylab, sep="")) +
        xlab(xlab) + theme_classic()
    if (exportTable == TRUE) { return(Con.df) }
    return(plot)}

#' Demonstrate the difference in clonal proportion between clonotypes
#'
#' This function produces an alluvial or area graph of the proportion of 
#' the indicated clonotypes for all or selected samples. Clonotypes can be 
#' selected using the clonotypes parameter with the specific sequence of 
#' interest or using the number parameter with the top n clonotypes by 
#' proportion to be visualized. If multiple clonotypes have the same proportion 
#' and are within the selection by the number parameter, all the clonotypes 
#' will be visualized. In this instance, if less clonotypes are desired, 
#' reduce the number parameter.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' compareClonotypes(combined, numbers = 10, 
#' samples = c("PX_P", "PX_T"), cloneCall="aa")
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param samples The specific samples to isolate for visualization.
#' @param clonotypes The specific sequences of interest.
#' @param numbers The top number clonotype sequences per group
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param graph The type of graph produced, either "alluvial" or "area".
#' @param exportTable Returns the data frame used for forming the graph.
#' @import ggplot2
#'
#' @export
#' @return ggplot of the proportion of total sequencing read of 
#' selecting clonotypes
compareClonotypes <- function(df, cloneCall = "strict", chain = "both", samples = NULL, 
                        clonotypes = NULL, numbers = NULL, split.by = NULL, 
                        graph = "alluvial", exportTable = FALSE){
    df <- list.input.return(df, split.by)
    cloneCall <- theCall(cloneCall)
    df <- checkBlanks(df, cloneCall)
    if (!is.null(numbers) & !is.null(clonotypes)) {
        stop("Make sure your inputs are either numbers or clonotype sequences.")
    }
    Con.df <- NULL
    for (i in seq_along(df)) {
      if (chain != "both") {
        df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
      }
        tbl <- as.data.frame(table(df[[i]][,cloneCall]))
        tbl[,2] <- tbl[,2]/sum(tbl[,2])
        colnames(tbl) <- c("Clonotypes", "Proportion")
        tbl$Sample <- names(df[i])
        Con.df <- rbind.data.frame(Con.df, tbl)
    }
    if (!is.null(samples)) {
        Con.df <- Con.df[Con.df$Sample %in% samples,] }
    if (!is.null(clonotypes)) {
        Con.df <- Con.df[Con.df$Clonotypes %in% clonotypes,] }
    if (!is.null(numbers)) {
        top <- Con.df %>%
          group_by(Con.df[,3]) %>%
          top_n(n = numbers, wt = Proportion)
        Con.df <- Con.df[Con.df$Clonotypes %in% top$Clonotypes,] }
    if (nrow(Con.df) < length(unique(Con.df$Sample))) {
        stop("Reasses the filtering strategies here, there is not 
            enough clonotypes to examine.") }
    if (exportTable == TRUE) { return(Con.df)}
    
    plot <- ggplot(Con.df, aes(x = Sample, fill = Clonotypes, group = Clonotypes,
                    stratum = Clonotypes, alluvium = Clonotypes, 
                    y = Proportion, label = Clonotypes)) +
                theme_classic() +
                theme(axis.title.x = element_blank())
    if (graph == "alluvial") {
        plot = plot +  geom_stratum() + geom_flow(stat = "alluvium")
    } else if (graph == "area") {
        plot = plot +
            geom_area(aes(group = Clonotypes), color = "black") }
    return(plot)
}

#' Scatter plot comparing the expansion of two samples
#'
#' This function produces a scatter plot directly comparing the specific clonotypes
#' between two samples. The clonotypes will be categorized by counts into singlets or multiplets, 
#' either exclusive or shared between the selected samples. Visualization inspired 
#' by the work of \href{https://pubmed.ncbi.nlm.nih.gov/32103181/}{Wu, T, et al 2020}. 
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' scatterClonotype(combined, x.axis = "PY_P", y.axis = "PY_T",
#' graph = "proportion")
#' 
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param x.axis name of the list element to appear on the x.axis
#' @param y.axis name of the list element to appear on the y.axis
#' @param dot.size either total or the name of the list element to 
#' use for size of dots
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param graph graph either proportion or raw clonotype count
#' @param exportTable Returns the data frame used for forming the graph.
#' 
#' @import ggplot2
#' 
#' @export
#' @return ggplot of the relative clonotype numbers

scatterClonotype <- function(df, cloneCall ="strict", 
                             x.axis = NULL, y.axis = NULL,
                             chain = "both",
                             dot.size = "total", 
                             split.by = NULL,
                             graph = "proportion", 
                             exportTable = FALSE) {
  df <- list.input.return(df, split.by)
  cloneCall <- theCall(cloneCall)
  axes <- which(names(df) %in% c(x.axis, y.axis, dot.size))
  if (chain != "both") {
    for (i in axes) {
      df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
    }
  }
  x.df <- as.data.frame(table(df[[x.axis]][,cloneCall]))
  colnames(x.df)[2] <- x.axis
  y.df <- as.data.frame(table(df[[y.axis]][,cloneCall]))
  colnames(y.df)[2] <- y.axis
  combined.df <- merge(x.df, y.df, by = "Var1", all = TRUE)
  if (dot.size != "total") {
    if (dot.size %!in% colnames(combined.df)) {
      size.df <- as.data.frame(table(df[[dot.size]][,cloneCall]))
      colnames(size.df)[2] <- dot.size
      combined.df <- merge(combined.df, size.df, by = "Var1", all = TRUE) }
    combined.df[is.na(combined.df)] <- 0
    combined.df[,paste0("size", ".fraction")] <- combined.df[,dot.size]/sum(combined.df[,dot.size])
    labeling <- unique(c(x.axis, y.axis, dot.size))
  } else {
    combined.df[is.na(combined.df)] <- 0
    labeling <- unique(c(x.axis, y.axis)) }
  combined.df[,"class"] <- NA
  combined.df[,"sum"] <- rowSums(combined.df[,labeling])
  for (i in seq_along(labeling)) {
    if (length(labeling) > 2) {
      combined.df[,"class"] <- ifelse(combined.df[,labeling[i]] == 1 & rowSums(combined.df[,labeling[which(labeling != labeling[i])]]) == 0, 
                                      paste0(labeling[i], ".singlet"), combined.df[,"class"])
      combined.df[,"class"] <- ifelse(combined.df[,labeling[i]] > 1 & rowSums(combined.df[,labeling[which(labeling != labeling[i])]]) == 0, 
                                      paste0(labeling[i], ".multiplet"), combined.df[,"class"])
    } else if (length(labeling) == 2) {
      combined.df[,"class"] <- ifelse(combined.df[,labeling[i]] == 1 & combined.df[,labeling[which(labeling != labeling[i])]] == 0, 
                                      paste0(labeling[i], ".singlet"), combined.df[,"class"])
      combined.df[,"class"] <- ifelse(combined.df[,labeling[i]] > 1 & combined.df[,labeling[which(labeling != labeling[i])]] == 0, 
                                      paste0(labeling[i], ".multiplet"), combined.df[,"class"])
  }}
  combined.df[,"class"] <- ifelse(combined.df[,y.axis] >= 1 & combined.df[,x.axis] >= 1, paste0("dual.expanded"), combined.df[,"class"])
  combined.df[,paste0(x.axis, ".fraction")] <- combined.df[,x.axis]/sum(combined.df[,x.axis])
  combined.df[,paste0(y.axis, ".fraction")] <- combined.df[,y.axis]/sum(combined.df[,y.axis])
  if (graph == "proportion") {
    x <- combined.df[,paste0(x.axis, ".fraction")]
    y <- combined.df[,paste0(y.axis, ".fraction")]
  } else if (graph == "count") {
    x <- combined.df[,x.axis]
    y <- combined.df[,y.axis] }
  if (dot.size != "total") {
    size <- combined.df[,dot.size]
  } else { size <- combined.df[,"sum"] }
  if (exportTable == TRUE) { return(combined.df) }
  plot <- ggplot(combined.df, aes(x=x, y = y, color = class)) + 
    theme_classic() + 
    scale_color_manual(values = colorblind_vector(length(unique(combined.df$class)))) + 
    xlab(x.axis) + ylab(y.axis) + labs(size = "Total n")
  if (graph == "proportion") {
    plot <- plot + geom_abline(slope = 1, intercept = 0, alpha = 0.4, lty=2)  + 
      scale_y_sqrt() + scale_x_sqrt() 
  } else if (graph == "count") {
    plot <- plot + ylim(0, max(x,y)) + xlim(0, max(x,y)) 
  }
  plot <- plot + geom_jitter(aes(size = size))
  return(plot)  
}

#' Hierarchical clustering of clonotypes on clonotype size and 
#' Jensen-Shannon divergence
#'
#' This function produces a hierarchical clustering of clonotypes by sample 
#' using the Jensen-Shannon distance and discrete gamma-GPD spliced threshold 
#' model in 
#' \href{https://bioconductor.org/packages/devel/bioc/html/powerTCR.html}{powerTCR R package}.
#' Please read and cite PMID: 30485278 if using the function for analyses. 
#' 
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' clonesizeDistribution(combined, cloneCall = "strict", method="ward.D2")
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene),
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param threshold Numerical vector containing the thresholds 
#' the grid search was performed over.
#' @param method The clustering paramater for the dendrogram.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param exportTable Returns the data frame used for forming the graph.
#' @importFrom  dplyr bind_rows
#' @importFrom ggplot2 ggplot
#' @importFrom powerTCR fdiscgammagpd get_distances
#' @export
#' @return ggplot dendrogram of the clone size distribution

clonesizeDistribution <- function(df,  cloneCall ="strict", 
                                  chain = "both", 
                                  method = "ward.D2", 
                                  threshold = 1, 
                                  split.by = NULL, 
                                  exportTable = FALSE) {
        df <- list.input.return(df, split.by)
        cloneCall <- theCall(cloneCall)
        df <- checkBlanks(df, cloneCall)
        data <- bind_rows(df)
        unique_df <- unique(data[,cloneCall])
        Con.df <- data.frame(matrix(NA, length(unique_df), length(df)))
        Con.df <- data.frame(unique_df, Con.df, stringsAsFactors = FALSE)
        colnames(Con.df)[1] <- "clonotype"
        for (i in seq_along(df)) {
          if (chain != "both") {
            df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
          }
            data <- df[[i]]
            data <- data.frame(table(data[,cloneCall]), 
                        stringsAsFactors = FALSE)
            colnames(data) <- c(cloneCall, "Freq")
            for (y in seq_along(unique_df)){
                    clonotype.y <- Con.df$clonotype[y]
                    location.y <- which(clonotype.y == data[,cloneCall])
                    Con.df[y,i+1] <- data[location.y[1],"Freq"] }
        }
        colnames(Con.df)[2:(length(df)+1)] <- names(df)
        Con.df[is.na(Con.df)] <- 0
        list <- list()
        for (i in seq_along(df)) {
            list[[i]] <- Con.df[,i+1]
            list[[i]] <- suppressWarnings(fdiscgammagpd(list[[i]], useq = threshold))
            }
        names(list) <- names(df)
        grid <- 0:10000
        distances <- get_distances(list, grid, modelType="Spliced")
        hclust <- hclust(as.dist(distances), method = method)
        hcd <- as.dendrogram(hclust)
        plot <- plot(hcd)
        if (exportTable == TRUE) { return(distances) }
        return(plot)
}

#This is the basic color palette for the package
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", "#7301A8FF", "#9C179EFF", 
              "#BD3786FF", "#D8576BFF","#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

#Making lodes to function in alluvial plots
#' @importFrom ggalluvial to_lodes_form
makingLodes <- function(meta2, color, alpha, facet, set.axes) {
    if (!is.null(color) & !is.null(alpha) & !is.null(facet)) {
        lodes <- to_lodes_form(meta2,key="x",value="stratum",id="alluvium",
        axes=set.axes,diffuse=c(as.name(color),as.name(alpha),as.name(facet)))
    } else  if (!is.null(color) & !is.null(alpha) & is.null(facet)) {
        lodes <- to_lodes_form(meta2,key="x",value="stratum",id="alluvium",
        axes = set.axes, diffuse = c(as.name(color), as.name(alpha)))
    } else if (!is.null(color) & is.null(alpha) & !is.null(facet)) {
        lodes <- to_lodes_form(meta2,key="x",value="stratum",id ="alluvium",
        axes =set.axes, diffuse = c(as.name(color), as.name(facet)))
    } else if (is.null(color) & is.null(alpha) & !is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", 
        id="alluvium",axes=set.axes,diffuse=c(as.name(alpha),as.name(facet)))
    } else if (is.null(color) & is.null(alpha) & !is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", 
        id = "alluvium", axes = set.axes, diffuse = c(as.name(facet)))
    } else if (!is.null(color) & is.null(alpha) & is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", 
        id = "alluvium", axes = set.axes, diffuse = c(as.name(color)))
    } else if (is.null(color) & !is.null(alpha) & is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", 
        id = "alluvium", axes = set.axes, diffuse = c(as.name(alpha)))
    } else { lodes <- to_lodes_form(meta2, key = "x", value = "stratum", 
            id = "alluvium", axes = set.axes)}
    return(lodes) }


#' Visualizing the distribution of gene usage
#' 
#' This function will allow for the visualizing the distribution 
#' of the any VDJ and C gene of the TCR or BCR using heatmap or
#' bar chart. This function requires assumes two chains were used in 
#' defining clonotype, if not, it will default to the only chain 
#' present regardless of the chain parameter.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' vizGenes(combined, gene = "V", chain = "TRB", plot = "bar", scale = TRUE)
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), 
#' or combineExpression().
#' @param gene Which part of the immune receptor to visualize - V, D, J, C
#' @param chain indicate the specific chain should be used - 
#' e.g. "TRA", "TRG", "IGH", "IGL" (no both option here)
#' @param plot The type of plot to return - heatmap or bar 
#' @param y.axis Variable to separate the y-axis, can be both categorical or other gene 
#' gene segments such as V, D, J, or C.
#' @param order Categorical variable to organize the x-axis, either "gene" or "variance"
#' @param scale Converts the proportion of total genes 
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param exportTable Returns the data frame used for forming the graph.
#' @import ggplot2
#' @importFrom stringr str_split
#' @importFrom stats sd
#' @export
#' @return ggplot bar diagram of vgene counts

vizGenes <- function(df, 
                      gene = "V",
                      chain = "TRA", 
                      plot = "heatmap",  
                      y.axis = "sample", 
                      order = "gene",
                      scale = TRUE, 
                      split.by = NULL,
                      exportTable = FALSE) {
    df <- list.input.return(df, split.by = split.by)
    for(i in seq_along(df)) {
      df[[i]] <- off.the.chain(df[[i]], chain, "CTgene")
    }
    df <- bound.input.return(df)
    if (y.axis %!in% colnames(df) | is.null(y.axis)) {
      if (y.axis %!in% c("V", "D", "J", "C")) {
        y.axis <- "element.names"
      } else {
        df <- select.gene(df, chain, y.axis)
        colnames(df)[ncol(df)] <- y.axis
      }
    }
    df <- select.gene(df, chain, gene)
    df <- subset(df, !is.na(df[,ncol(df)])) #remove NA values
    df <- subset(df, df[,ncol(df)] != "NA") #remove values that are character "NA"
    df <- subset(df, df[,ncol(df)] != "") #remove rows with non genes
    #df <- table(df[,ncol(df)], df[,y.axis])
    
    if (!is.null(y.axis) && y.axis != "element.names") {
      df <- df %>%
        group_by(df[,ncol(df)], df[,y.axis], element.names) %>%
        count() 
    } else {
      df <- df %>%
        group_by(df[,ncol(df)], element.names) %>%
        count() 
    }
    df <- df %>%
      group_by(element.names) %>%
      mutate(sum = sum(n))
    col.lab <- "Total n"
    if (scale == TRUE) {
     df[,"n"] <- df[,"n"]/df[,"sum"]
     col.lab <- "Scaled Values"
    } 
    colnames(df)[1:2] <- c("Var1", "Var2")
    df <- df %>%
      group_by(Var1, Var2) %>%
      mutate(varcount = sum(n), 
             sd = sd(n, na.rm = TRUE),
             mean = mean(n))
    if (order == "variance") {
      varOrder <- order(df$varcount, decreasing = TRUE)
      df$Var1 <- factor(df$Var1, levels = unique(df$Var1[varOrder]))
    }
    if (plot == "bar") {
      df2 <- unique(df[,c("Var1", "Var2", "sd", "mean")])
      plot <- ggplot(df2, aes(x=Var1, y = mean)) + 
          geom_bar(stat = "identity") + 
          geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                        position=position_dodge(.9)) + 
          theme_classic() + 
          theme(axis.title.x = element_blank(), #remove titles
                axis.title.y = element_blank(), 
                  axis.ticks.x = element_blank(), #removes ticks
                  axis.text.x = element_text(angle = 90, 
                                      vjust = 0.5, hjust=1, size=rel(0.5))) + 
        facet_grid(Var2~.)
    } else if (plot == "heatmap") {
     
      plot <- ggplot(df, aes(x=Var1, y = Var2, fill = n)) + 
              geom_tile() + 
              theme_classic() + 
              labs(fill = col.lab) + 
              theme(axis.title.x = element_blank(), #remove titles
              axis.ticks.x = element_blank(), #removes ticks
              axis.text.x = element_text(angle = 90, 
                                       vjust = 0.5, hjust=1, size=rel(0.5)), 
              axis.title.y = element_blank(), 
              axis.text.y = element_text(size=rel(0.5))) + 
              scale_fill_gradientn(colors = rev(colorblind_vector(5)))
    }
    if (exportTable == TRUE) { return(df) }
    return(plot)
}
