#' Calculate Startrac-based Diversity Indices 
#' 
#' @description This function utilizes the STARTRAC approach to calculate T cell 
#' diversity metrics based on the work of Zhang et al. (2018, Nature) 
#' [PMID: 30479382](https://pubmed.ncbi.nlm.nih.gov/30479382/). It can compute 
#' three distinct indices: clonal expansion (`expa`), cross-tissue migration 
#' (`migr`), and state transition (`tran`). 
#' 
#' @details
#' The function requires a `type` variable in the metadata, which specifies the
#' tissue origin or any other categorical variable for migration analysis.
#'
#' **Indices:**
#' \itemize{
#'   \item{\strong{expa (Clonal Expansion):}} Measures the extent of clonal 
#'         proliferation within a T cell cluster. It is calculated as 
#'         `1 - normalized Shannon entropy`. A higher value indicates greater 
#'         expansion of a few clones.
#'   \item{\strong{migr (Cross-Tissue Migration):}} Quantifies the movement of 
#'         clonal T cells across different tissues (as defined by the `type`
#'         parameter). It is based on the entropy of a clonotype's distribution 
#'         across tissues.
#'   \item{\strong{tran (State Transition):}} Measures the developmental 
#'         transition of clonal T cells between different functional clusters. 
#'         It is based on the entropy  of a clonotype's distribution across 
#'         clusters.
#' }
#'
#' **Pairwise Analysis:**
#' The `pairwise` parameter enables the calculation of migration or transition
#' between specific pairs of tissues or clusters, respectively.
#' \itemize{
#'   \item{For migration (`index = "migr"`), set `pairwise` to the `type` column
#'         (e.g., `pairwise = "Type"`).}
#'   \item{For transition (`index = "tran"`), set `pairwise` to `"cluster"`.}
#' }
#'
#' The function loops over every **unordered** pair of clusters (or types),
#' subsets the cells to that pair, and recomputes the index inside the resulting
#' 2-category subspace. The exported table has one row per pair per anchor
#' category:
#' \itemize{
#'   \item{\strong{group:}} the `group.by` level the row was computed in.
#'   \item{\strong{cluster:}} the category the score is anchored to, i.e. the
#'         one supplying the cell weights.
#'   \item{\strong{partner:}} for `index = "tran"`, the other member of the
#'         pair.
#'   \item{\strong{comparison:}} the pair itself, always written in sorted
#'         order so the same pair carries the same label in every group.
#'   \item{\strong{value:}} the cell-weighted mean clonal entropy over the pair.
#'         With two categories it is bounded 0-1, where 0 means the clones
#'         making up the anchor are exclusive to it and 1 means every clone is
#'         split evenly with the partner.
#' }
#'
#' The pair label carries **no direction**: `"1 vs 3"` is not "from 1 to 3". Two
#' rows of the same pair can differ because each is weighted by a different
#' anchor, not because a transition was measured one way or the other. STARTRAC
#' indices are a static description of clonal sharing and cannot, on their own,
#' order that sharing in time.
#'
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' # Getting a sample of a Seurat object
#' scRep_example  <- get(data("scRep_example"))
#' scRep_example  <- combineExpression(combined, scRep_example)
#' scRep_example$Patient <- substring(scRep_example$orig.ident,1,3)
#' scRep_example$Type <- substring(scRep_example$orig.ident,4,4) 
#' 
#' # Calculate a single index (expansion)
#' StartracDiversity(scRep_example, 
#'                   type = "Type", 
#'                   group.by = "Patient",
#'                   index = "expa")
#'                   
#' # Calculate pairwise transition 
#' StartracDiversity(scRep_example, 
#'                   type = "Type", 
#'                   group.by = "Patient",
#'                   index = "tran",
#'                   pairwise = "cluster") 
#'
#' @param sc.data The single-cell object after [combineExpression()].
#' For SCE objects, the cluster variable must be in the meta data under
#' "cluster".
#' @param clone.call Defines the clonal sequence grouping. Accepted values
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL`,
#' `IGK`, `Light` (for both light chains), or `both` (for TRA/B and Heavy/Light).
#' @param index A character vector specifying which indices to calculate.
#' Options: "expa", "migr", "tran". Default is all three.
#' @param type The metadata variable that specifies tissue type for migration
#' analysis.
#' @param group.by A column header in the metadata or lists to group the analysis
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed as
#' by list element or active identity in the case of single-cell objects.
#' @param pairwise The metadata column to be used for pairwise comparisons.
#' Set to the `type` variable for pairwise migration or "cluster" for
#' pairwise transition.
#' @param export.table If `TRUE`, returns a data frame or matrix of the results
#' instead of a plot.
#' @param palette Colors to use in visualization - input any
#' [hcl.pals][grDevices::hcl.pals].
#' @param cloneCall \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `clone.call` instead.
#' @param exportTable \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `export.table` instead.
#' @param ... Additional arguments passed to the ggplot theme
#' @importFrom stats reshape
#' @export
#' @concept SC_Functions
#' @return A ggplot object visualizing STARTRAC diversity metrics or data.frame if
#'`exportTable = TRUE`.
#' @author Liangtao Zheng
StartracDiversity <- function(sc.data,
                              clone.call = NULL,
                              chain = "both",
                              index = c("expa", "migr", "tran"),
                              type = NULL,
                              group.by = NULL,
                              pairwise = NULL,
                              export.table = NULL,
                              palette = "inferno",
                              # Deprecated arguments
                              cloneCall = NULL,
                              exportTable = NULL,
                              ...) {

  # Handle deprecated arguments
  clone.call <- .deprecate_arg(cloneCall, clone.call, "cloneCall", "clone.call",
                               "StartracDiversity", default = "strict")
  export.table <- .deprecate_arg(exportTable, export.table, "exportTable", "export.table",
                                 "StartracDiversity", default = FALSE)

  if(!all(index %in% c("expa", "migr", "tran"))) {
    stop("Please select 'expa', 'migr', and/or 'tran' for index.")
  }
  if (!is.null(pairwise) && length(index) > 1) {
    stop("Pairwise analysis can only be performed for a single index ('migr' or 'tran').")
  }
  if (!is.null(pairwise) && !index %in% c("migr", "tran")) {
    stop("Pairwise analysis is only supported for 'migr' or 'tran' indices.")
  }
  
  
  # Prepare data
  df <- .grabMeta(sc.data)
  clone.call <- .theCall(df, clone.call)
  barcodes <- rownames(df)
  colnames(df)[ncol(df)] <- "cluster"
    
  if (is.null(group.by)) {
    if (!"orig.ident" %in% colnames(df)) {
      stop("Please select a group.by variable.")
    }
    group.by <- "orig.ident"
  }
  group.levels <- unique(df[,group.by])
  
  if (chain != "both") {
    df <- .offTheChain(df, chain, clone.call)
  }

  # Process clonotypes
  df <- df %>%
    group_by(across(all_of(c(group.by, clone.call)))) %>%
    dplyr::mutate(n = n()) %>%
    as.data.frame()

  rownames(df) <- barcodes
  remove.pos <- which(is.na(df[,clone.call]) | df[,clone.call] == "")
  if (length(remove.pos) > 0) {
    df <- df[-remove.pos,]
  }

  processed <- data.frame(
    Cell_Name = rownames(df),
    clone.id = df[,clone.call],
    patient = df[,group.by],
    cluster = df[,"cluster"],
    loc = df[,type],
    stringsAsFactors = FALSE
  )
  processed[processed == "NA"] <- NA
  processed <- na.omit(processed)
  
  # Calculate indices
  mat.list <- lapply(group.levels, function(level) {
    subset_data <- processed[processed$patient == level,]
    if (!is.null(pairwise)) {
      comparison_col <- if (index == "migr") {
        "loc"
      } else { # index == "tran"
        "cluster"
      }
      .calculatePairwiseIndices(subset_data, index, comparison_col)
    } else {
      .calculateIndices(subset_data, index)
    }
  })
  names(mat.list) <- as.character(group.levels)

  mat <- bind_rows(mat.list, .id = "group")
  if (nrow(mat) == 0) {
    warning("No data available for calculation. Returning NULL.")
    return(NULL)
  }

  if(!is.null(pairwise)) {
    mat$variable <- index[1]
    mat <- mat[!is.nan(mat$value),]
    col.order <- intersect(c("group", "cluster", "partner", "comparison",
                             "variable", "value"), colnames(mat))
    mat <- mat[, c(col.order, setdiff(colnames(mat), col.order)), drop = FALSE]
    rownames(mat) <- NULL
  }

  if (export.table) {
    return(mat)
  }
  # Plotting logic
  if (!is.null(pairwise)) {
    
    # "partner" is only present when the pairing was over clusters, which is
    # the one case where the anchor is itself a member of the pair
    if ("partner" %in% colnames(mat)) {
        # One facet per anchor cluster, one box per partner within it, so every
        # box is a single unambiguous pair summarized across the group.by levels
        cluster.levels <- .alphanumericalSort(c(mat$cluster, mat$partner))
        num_colors <- length(cluster.levels)
        mat$cluster <- factor(mat$cluster, levels = cluster.levels)
        mat$partner <- factor(mat$partner, levels = cluster.levels)
        plot <- ggplot(mat, aes(x = .data$partner, y = .data$value)) +
          geom_boxplot(aes(fill = .data$partner), outlier.alpha = 0, na.rm = TRUE) +
          labs(y = "Pairwise Index Score", x = "Partner Cluster") +
          facet_wrap(~ cluster, labeller = as_labeller(function(x) paste("Anchor:", x)))
    } else {
      col_name <- colnames(mat)[grepl("comparison", colnames(mat))]
      num_colors <- length(unique(mat[[col_name]]))
      plot <- ggplot(mat, aes(x = .data[[col_name]], y = .data$value)) +
        geom_boxplot(aes(fill = .data[[col_name]]), outlier.alpha = 0, na.rm = TRUE) +
        labs(y = "Pairwise Index Score", x = "Comparison")
    }
  } else {
    mat_melt <- reshape(mat,
                        varying = index,
                        v.names = "value",
                        timevar = "variable",
                        times = index,
                        direction = "long")
    values <- .alphanumericalSort(unique(mat_melt$cluster))
    mat_melt$cluster <- factor(mat_melt$cluster, levels = values)
    mat_melt$value <- as.numeric(mat_melt$value)
    num_colors <- length(unique(mat_melt$cluster))
    
    plot <- ggplot(mat_melt, aes(x = cluster, y = .data[["value"]])) +
      geom_boxplot(aes(fill = cluster), outlier.alpha = 0, na.rm = TRUE) +
      labs(y = "Index Score", 
           x = "Clusters") +
      theme(axis.title.x = element_blank())
  }
  
  plot <- plot + 
    .themeRepertoire(...) + 
    guides(fill = "none") +
    scale_fill_manual(values = .colorizer(palette, num_colors))
  
  if(length(index) > 1) {
    plot <- plot + facet_grid(variable ~ ., scales = "free_y") 
  } 
  
  return(plot)
}


# Helper function for standard index calculation
.calculateIndices <- function(processed, indices) {
  if (nrow(processed) == 0) return(NULL)
  clonotype.dist.cluster <- table(processed[,c("clone.id", "cluster")])
  clonotype.dist.loc <- table(processed[,c("clone.id", "loc")])
  
  # Return NULL if no clusters are found
  if (ncol(clonotype.dist.cluster) == 0) return(NULL)
  
  calIndex.matrix <- data.frame(cluster = colnames(clonotype.dist.cluster))
  
  if ("expa" %in% indices) {
    entropy_val <- .mcolEntropy(clonotype.dist.cluster)
    entropy_max <- log2(colSums(clonotype.dist.cluster > 0))
    expa <- 1 - (entropy_val / entropy_max)
    calIndex.matrix$expa <- expa
  }
  
  # Check if there are clonotypes to process for migr/tran
  if (nrow(clonotype.dist.cluster) > 0) {
    clonotype.data <- data.frame(clone.id = rownames(clonotype.dist.cluster))
    weights.mtx <- sweep(clonotype.dist.cluster, 2, colSums(clonotype.dist.cluster), "/")
    
    if ("migr" %in% indices && "loc" %in% colnames(processed) && nrow(clonotype.dist.loc) > 0) {
      # Ensure clonotypes for migration calculation exist in the cluster distribution
      shared_clones_migr <- intersect(rownames(clonotype.dist.cluster), rownames(clonotype.dist.loc))
      if(length(shared_clones_migr) > 0) {
        clonotype.data$migr <- .mrowEntropy(clonotype.dist.loc[shared_clones_migr,,drop=FALSE])
        migr_matrix <- t(weights.mtx[shared_clones_migr, , drop=FALSE]) %*% as.matrix(clonotype.data$migr)
        calIndex.matrix$migr <- migr_matrix[,1]
      } else {
        calIndex.matrix$migr <- NA
      }
    }
    
    if ("tran" %in% indices) {
      clonotype.data$tran <- .mrowEntropy(clonotype.dist.cluster)
      tran_matrix <- t(weights.mtx) %*% as.matrix(clonotype.data$tran)
      calIndex.matrix$tran <- tran_matrix[,1]
    }
  } else {
    # If no clonotypes, set indices to NA
    if ("migr" %in% indices) calIndex.matrix$migr <- NA
    if ("tran" %in% indices) calIndex.matrix$tran <- NA
  }
  
  for (col in names(calIndex.matrix)) {
    if (is.numeric(calIndex.matrix[[col]])) {
      calIndex.matrix[[col]][is.nan(calIndex.matrix[[col]]) | is.infinite(calIndex.matrix[[col]])] <- NA
    }
  }
  
  return(calIndex.matrix)
}

# Helper function for pairwise index calculation
#' @importFrom utils combn
.calculatePairwiseIndices <- function(processed, index, pairwise_col) {
  if (nrow(processed) < 2) return(NULL)
  
  # Sorting keeps the pair label identical across group.by levels - unique()
  # alone returns categories in order of first appearance, so the same
  # unordered pair would be written "1 vs 3" in one group and "3 vs 1" in the next
  unique_items <- .alphanumericalSort(as.character(processed[[pairwise_col]]))
  if (length(unique_items) < 2) return(NULL)

  pairs <- combn(unique_items, 2, simplify = FALSE)

  pairwise_results <- lapply(pairs, function(p) {
    pair_data <- processed[as.character(processed[[pairwise_col]]) %in% p,]
    pair_data <- droplevels(pair_data)

    if (index == "migr") {
      dist_table <- table(pair_data[,c("clone.id", "loc")])
      clonotype_dist_cluster <- table(pair_data[,c("clone.id", "cluster")])
    } else { # tran
      dist_table <- table(pair_data[,c("clone.id", "cluster")])
      clonotype_dist_cluster <- dist_table
    }
    
    if(nrow(dist_table) == 0) return(NULL)
    
    clonotype_data <- data.frame(clone.id = rownames(dist_table),
                                 value = .mrowEntropy(dist_table))
    
    weights_mtx <- sweep(clonotype_dist_cluster, 2, colSums(clonotype_dist_cluster), "/")
    
    # Ensure clone.ids match for matrix multiplication
    shared_clones <- intersect(rownames(weights_mtx), clonotype_data$clone.id)
    if (length(shared_clones) == 0) return(NULL)
    
    weights_mtx_filtered <- weights_mtx[shared_clones,, drop=FALSE]
    clonotype_data_filtered <- clonotype_data[clonotype_data$clone.id %in% shared_clones,]
    
    result_matrix <- t(weights_mtx_filtered) %*% as.matrix(clonotype_data_filtered$value)
    
    res <- data.frame(cluster = rownames(result_matrix),
                      value = result_matrix[,1],
                      stringsAsFactors = FALSE)
    res$comparison <- paste(p, collapse = " vs ")
    # For transition the anchor is itself a member of the pair, so the other
    # member can be named outright rather than parsed out of the label
    if (identical(pairwise_col, "cluster")) {
      res$partner <- ifelse(res$cluster == p[1], p[2], p[1])
    }
    rownames(res) <- NULL
    return(res)
  })
  
  bind_rows(pairwise_results)
}

# Entropy of each row of the input matrix
.mrowEntropy <- function(x) {
  freqs <- sweep(x, 1, rowSums(x), "/")
  H <- -rowSums(ifelse(freqs > 0, freqs * log2(freqs), 0))
  return(H)
}

# Entropy of each column of the input matrix
.mcolEntropy <- function(x) {
  freqs <- sweep(x, 2, colSums(x), "/")
  H <- -colSums(ifelse(freqs > 0, freqs * log2(freqs), 0))
  return(H)
}

