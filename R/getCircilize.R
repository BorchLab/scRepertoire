#' Generate Data Frame to Plot Chord Diagram
#'
#' This function will take the meta data from the product of
#' [combineExpression()] and generate a relational data frame to
#' be used for a chord diagram. Each chord will represent the number of
#' clones unique and shared across the multiple `group.by` variable.
#' If using the downstream circlize R package, please read and cite the
#' following [manuscript](https://pubmed.ncbi.nlm.nih.gov/24930139/).
#' If looking for more advanced ways for circular visualizations, there
#' is a great [cookbook](https://jokergoo.github.io/circlize_book/book/)
#' for the circlize package.
#'
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list,
#'                         samples = c("P17B", "P17L", "P18B", "P18L",
#'                                     "P19B","P19L", "P20B", "P20L"))
#'
#' # Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' scRep_example <- combineExpression(combined,
#'                                    scRep_example)
#'
#' # Getting data frame output for Circlize
#' circles <- getCirclize(scRep_example,
#'                        group.by = "seurat_clusters")
#'
#' # Multi-level grouping for hierarchical chord diagrams
#' scRep_example$Patient <- substring(scRep_example$orig.ident, 1, 3)
#' circles <- getCirclize(scRep_example,
#'                        group.by = c("Patient", "seurat_clusters"))
#'
#' # Get rich output with sector metadata
#' result <- getCirclize(scRep_example,
#'                       group.by = "seurat_clusters",
#'                       include.metadata = TRUE)
#'
#' @param sc.data The single-cell object after [combineExpression()].
#' @param clone.call Defines the clonal sequence grouping. Accepted values
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param group.by A column header (or vector of column headers for hierarchical
#' grouping) in the metadata to group the analysis by (e.g., "sample", "treatment").
#' If `NULL`, data will be analyzed by active identity. When multiple columns are
#' provided, they are combined with "_" separator for multi-level annotations.
#' @param method The method for calculating link values: `"unique"` (default) counts
#' unique shared clones, `"abundance"` sums clone frequencies, `"jaccard"` calculates
#' Jaccard similarity, `"overlap"` calculates overlap coefficient.
#' @param proportion Calculate the relationship by unique clones (`FALSE`, default)
#' or normalized by proportion (`TRUE`).
#' @param symmetric If `TRUE` (default), returns symmetric relationships. If `FALSE`,
#' returns directional flow showing proportion of source's clones found in destination.
#' @param include.self Include counting the clones within a single group.by comparison.
#' @param include.metadata If `TRUE`, returns a list with links data frame and
#' sector-level metadata including cell counts, clone counts, and expansion metrics.
#' @param min.shared Minimum number of shared clones to include a link (default 0).
#' @param top.links Keep only the top N links by value. If `NULL` (default), keep all.
#' @param filter.sectors Character vector of sectors to include. If `NULL`, include all.
#' @param palette Colors to use for sector color suggestions - input any
#' [hcl.pals][grDevices::hcl.pals].
#' @param cloneCall \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `clone.call` instead.
#'
#' @export
#' @concept SC_Functions
#' @return A data frame of shared clones between groups formatted for
#' [chordDiagram][circlize::chordDiagram]. If `include.metadata = TRUE`, returns
#' a list with `links` (the edge data frame), `sectors` (sector-level statistics),
#' and `colors` (suggested colors for each sector).
#' @author Dillon Corvino, Nick Borcherding
getCirclize <- function(sc.data,
                        clone.call = NULL,
                        group.by = NULL,
                        method = c("unique", "abundance", "jaccard", "overlap"),
                        proportion = FALSE,
                        symmetric = TRUE,
                        include.self = TRUE,
                        include.metadata = FALSE,
                        min.shared = 0,
                        top.links = NULL,
                        filter.sectors = NULL,
                        palette = "inferno",
                        # Deprecated arguments
                        cloneCall = NULL) {

  # Handle deprecated arguments
  clone.call <- .deprecate_arg(cloneCall, clone.call, "cloneCall", "clone.call",
                               "getCirclize", default = "strict")

  method <- match.arg(method)

  meta <- .grabMeta(sc.data)
  clone.call <- .theCall(meta, clone.call)

  if(is.null(group.by)) {
    group.by <- "ident"
  }


  # Handle multi-level grouping
  if(length(group.by) > 1) {
    # Create compound grouping variable
    meta$`.compound_group` <- apply(meta[, group.by, drop = FALSE], 1,
                                    function(x) paste(x, collapse = "_"))
    group.var <- ".compound_group"
    hierarchy.info <- meta[, c(group.by, ".compound_group"), drop = FALSE]
    hierarchy.info <- unique(hierarchy.info)
  } else {
    group.var <- group.by
    hierarchy.info <- NULL
  }

  # Filter sectors if specified
  if(!is.null(filter.sectors)) {
    meta <- meta[meta[[group.var]] %in% filter.sectors, ]
    if(nrow(meta) == 0) {
      stop("No cells remaining after filtering by filter.sectors")
    }
  }

  unique.sectors <- unique(meta[[group.var]])

  # Build group pairs based on symmetry
  if(symmetric) {
    group_pairs <- expand.grid(group1 = unique.sectors, group2 = unique.sectors,
                               stringsAsFactors = FALSE)
    group_pairs <- unique(t(apply(group_pairs, 1, function(x) sort(x))))
    group_pairs <- as.data.frame(group_pairs, stringsAsFactors = FALSE)
    colnames(group_pairs) <- c("from", "to")
  } else {
    # Asymmetric: all ordered pairs
    group_pairs <- expand.grid(from = unique.sectors, to = unique.sectors,
                               stringsAsFactors = FALSE)
  }

  if(!include.self) {
    group_pairs <- group_pairs[group_pairs[,1] != group_pairs[,2], ]
  }

  # Count clones across all identities
  clone.table <- .cloneCounter(meta, group.var, clone.call)
  clone.table[[clone.call]] <- as.character(clone.table[[clone.call]])

  group_pairs$value <- NA

  for(i in seq_len(nrow(group_pairs))) {
    pair1 <- group_pairs[i, 1]
    pair2 <- group_pairs[i, 2]

    clone1.data <- clone.table[clone.table[,1] == pair1 & clone.table[["n"]] > 0, ]
    clone2.data <- clone.table[clone.table[,1] == pair2 & clone.table[["n"]] > 0, ]

    clone1 <- clone1.data[[clone.call]]
    clone2 <- clone2.data[[clone.call]]

    common <- intersect(clone1, clone2)

    # Calculate value based on method
    if(method == "unique") {
      value <- length(common)

      if(symmetric && pair1 == pair2) {
        # For self-links in symmetric mode, subtract shared clones
        tmp <- clone.table[clone.table[[clone.call]] %in% common &
                           clone.table[,1] != pair1 &
                           clone.table[["n"]] > 0, ]
        shared.clones <- unique(tmp[[clone.call]])
        value <- value - length(shared.clones)
      }

    } else if(method == "abundance") {
      # Sum of minimum frequencies for shared clones
      if(length(common) == 0) {
        value <- 0
      } else {
        freq1 <- clone1.data[clone1.data[[clone.call]] %in% common, "n"]
        names(freq1) <- clone1.data[clone1.data[[clone.call]] %in% common, clone.call]
        freq2 <- clone2.data[clone2.data[[clone.call]] %in% common, "n"]
        names(freq2) <- clone2.data[clone2.data[[clone.call]] %in% common, clone.call]
        value <- sum(pmin(freq1[common], freq2[common]))
      }

    } else if(method == "jaccard") {
      union_clones <- union(clone1, clone2)
      if(length(union_clones) == 0) {
        value <- 0
      } else {
        value <- length(common) / length(union_clones)
      }

    } else if(method == "overlap") {
      min_size <- min(length(clone1), length(clone2))
      if(min_size == 0) {
        value <- 0
      } else {
        value <- length(common) / min_size
      }
    }

    # Apply proportion normalization if requested (for unique method)
    if(proportion && method == "unique") {
      denominator <- length(unique(clone.table[clone.table[,1] == pair2, clone.call]))
      if(denominator > 0) {
        value <- value / denominator
      } else {
        value <- 0
      }
    }

    # For asymmetric mode, calculate directional proportion
    if(!symmetric && pair1 != pair2) {
      n_source <- length(clone1)
      if(n_source > 0) {
        value <- length(common) / n_source
      } else {
        value <- 0
      }
    }

    group_pairs$value[i] <- value
  }

  # Filter by minimum shared value
  if(min.shared > 0) {
    group_pairs <- group_pairs[group_pairs$value >= min.shared, ]
  }

  # Keep only top links
  if(!is.null(top.links) && nrow(group_pairs) > top.links) {
    group_pairs <- group_pairs[order(group_pairs$value, decreasing = TRUE), ]
    group_pairs <- group_pairs[seq_len(top.links), ]
  }

  # If not including metadata, return simple data frame
 if(!include.metadata) {
    return(group_pairs)
  }

  # Build sector metadata
  sector.stats <- lapply(unique.sectors, function(s) {
    sector.meta <- meta[meta[[group.var]] == s, ]
    sector.clones <- clone.table[clone.table[,1] == s, ]

    n.cells <- nrow(sector.meta)
    n.clones <- nrow(sector.clones[sector.clones$n > 0, ])

    # Count shared clones (appearing in other sectors)
    sector.clone.ids <- sector.clones[sector.clones$n > 0, clone.call]
    other.clones <- clone.table[clone.table[,1] != s & clone.table$n > 0, clone.call]
    n.shared <- length(intersect(sector.clone.ids, other.clones))

    # Simple expansion metric: 1 - (unique clones / total cells)
    expansion <- ifelse(n.cells > 0, 1 - (n.clones / n.cells), 0)

    data.frame(
      sector = s,
      n.cells = n.cells,
      n.clones = n.clones,
      n.shared = n.shared,
      expansion = expansion,
      stringsAsFactors = FALSE
    )
  })
  sector.stats <- do.call(rbind, sector.stats)

  # Add hierarchy info if multi-level grouping
  if(!is.null(hierarchy.info)) {
    sector.stats <- merge(sector.stats, hierarchy.info,
                          by.x = "sector", by.y = ".compound_group",
                          all.x = TRUE)
  }

  # Generate colors
  n.sectors <- length(unique.sectors)
  colors <- .colorizer(palette, n.sectors)
  names(colors) <- unique.sectors

  result <- list(
    links = group_pairs,
    sectors = sector.stats,
    colors = colors
  )

  return(result)
}
