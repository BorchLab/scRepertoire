#' Cluster clones by sequence similarity
#'
#' This function clusters TCRs or BCRs based on the edit distance or alignment
#' score of their CDR3 sequences. It can operate on either nucleotide (`nt`)
#' or amino acid (`aa`) sequences and can optionally enforce that clones share
#' the same V and/or J genes. The output can be the input object with an added
#' metadata column for cluster IDs, a sparse adjacency matrix, or an `igraph`
#' graph object representing the cluster network.
#'
#' @details
#' The clustering process is as follows:
#' 1.  The function retrieves the relevant chain data from the input object.
#' 2.  It calculates the distance between all sequences within each group
#'     (or across the entire dataset if `group.by` is `NULL`).
#' 3.  An edge list is constructed, connecting sequences that meet the similarity
#'     `threshold`.
#' 4.  The `threshold` parameter behaves differently based on its value:
#'     - **`threshold` < 1 (e.g., 0.85):** Interpreted as a *normalized*
#'       distance. A higher value means greater similarity is required.
#'     - **`threshold` >= 1 (e.g., 2):** Interpreted as a maximum *raw* edit
#'       distance. A lower value means greater similarity is required.
#' 5.  **Distance Metrics:**
#'     - **Levenshtein/Hamming/Damerau:** Standard edit distance calculations.
#'       `"hamming"` is only defined for equal-length sequences; pairs of
#'       unequal length are assigned the maximum distance and so never cluster.
#'       CDR3 lengths vary, so the engine's equal-length notice is expected
#'       and most pairs are excluded by length alone. This mirrors the common
#'       BCR convention of partitioning by V gene, J gene, and junction length
#'       before scoring (use `use.V = TRUE`, `use.J = TRUE`). Use
#'       `"levenshtein"` if you want sequences of differing length to be
#'       comparable.
#'     - **Alignment (NW/SW):** If `dist.type` is "nw" (Needleman-Wunsch) or
#'       "sw" (Smith-Waterman), alignment scores are calculated using the
#'       specified substitution matrix (`dist.mat`). These scores are converted
#'       to a distance-like metric for clustering.
#' 6.  An `igraph` graph is built from the edge list.
#' 7.  A clustering algorithm is run on the graph (default: connected components).
#'
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list,
#'                        samples = c("P17B", "P17L", "P18B", "P18L",
#'                                    "P19B","P19L", "P20B", "P20L"))
#'
#' # Standard Levenshtein clustering (85% similarity)
#' sub_combined <- clonalCluster(combined[c(1,2)],
#'                               chain = "TRA",
#'                               sequence = "aa",
#'                               threshold = 0.85)
#'
#' # Alignment-based clustering using BLOSUM80
#' sub_combined_nw <- clonalCluster(combined[c(1,2)],
#'                                  chain = "TRA",
#'                                  dist.type = "nw",
#'                                  dist.mat = "BLOSUM80",
#'                                  threshold = 0.85)
#'
#' # Export the graph object instead
#' graph_obj <- clonalCluster(combined[c(1,2)],
#'                            chain = "TRA",
#'                            export.graph = TRUE)
#'
#' @param input.data The product of [combineTCR()],
#' [combineBCR()] or [combineExpression()].
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL`,
#' `IGK`, `Light` (for both light chains), or `both` (for TRA/B and Heavy/Light).
#' @param sequence Clustering based on either `aa` or `nt` sequences. This
#' selects `CTaa` or `CTnt` as the source of the CDR3 used for all distance
#' calculations.
#' @param threshold The similarity threshold. If < 1, treated as normalized
#' similarity (higher is stricter). If >= 1, treated as raw edit distance
#' (lower is stricter).
#' @param group.by A column header in the metadata or lists to group the analysis
#' by (e.g., "sample", "treatment"). If `NULL`, clusters will be calculated across
#' all sequences.
#' @param use.V If `TRUE`, sequences must share the same V gene to be
#' clustered together.
#' @param use.J If `TRUE`, sequences must share the same J gene to be
#' clustered together.
#' @param dist.type The distance metric to use. Options: `"levenshtein"` (default),
#' `"hamming"`, `"damerau"` (allows transpositions), `"nw"` (Needleman-Wunsch),
#' or `"sw"` (Smith-Waterman).
#' @param dist.mat The substitution matrix to use for alignment-based metrics
#' (`"nw"` or `"sw"`). Options: `"BLOSUM45"`, `"BLOSUM50"`, `"BLOSUM62"`,
#' `"BLOSUM80"` (default), `"BLOSUM100"`, `"PAM30"`, `"PAM40"`, `"PAM70"`, `"PAM120"`,
#' `"PAM250"`, or `"identity"`.
#' @param normalize Method for normalizing distances. Options: `"none"`,
#' `"maxlen"` (divide by max sequence length), or `"length"` (default, divide
#' by mean sequence length). If `threshold < 1`, this controls how the
#' similarity is calculated.
#' @param gap.open Penalty for opening a gap in alignment metrics (default: -10).
#' @param gap.extend Penalty for extending a gap in alignment metrics (default: -1).
#' @param cluster.method The clustering algorithm to use. Defaults to `"components"`,
#' which finds connected subgraphs.
#' @param cluster.prefix A character prefix to add to the cluster names (e.g.,
#' "cluster.").
#' @param export.graph If `TRUE`, the function returns an `igraph`
#' object of the sequence network.
#' @param export.adj.matrix If `TRUE`, the function returns a sparse
#' adjacency matrix (`dgCMatrix`) of the network.
#' @param dist_type \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `dist.type` instead.
#' @param dist_mat \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `dist.mat` instead.
#' @param gap_open \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `gap.open` instead.
#' @param gap_extend \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `gap.extend` instead.
#' @param exportGraph \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `export.graph` instead.
#' @param exportAdjMatrix \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `export.adj.matrix` instead.
#' @importFrom igraph graph_from_edgelist E E<- V V<- as_data_frame
#' as_adjacency_matrix membership set_vertex_attr
#' @importFrom dplyr left_join ungroup
#' @importFrom rlang %||%
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors DataFrame
#' @importFrom Matrix sparseMatrix
#' @importFrom stats reshape
#'
#' @export
#' @concept Visualizing_Clones
#' @return
#' Depending on the export parameters, one of the following:
#' \itemize{
#'   \item An amended `input.data` object with a new metadata column containing cluster IDs (default).
#'   \item An `igraph` object if `export.graph = TRUE`.
#'   \item A sparse `dgCMatrix` object if `export.adj.matrix = TRUE`.
#' }

clonalCluster <- function(input.data,
                          chain = "TRB",
                          sequence = "aa",
                          threshold = 0.85,
                          group.by = NULL,
                          dist.type = NULL,
                          dist.mat = NULL,
                          normalize = "length",
                          gap.open = NULL,
                          gap.extend = NULL,
                          cluster.method = "components",
                          cluster.prefix = "cluster.",
                          use.V = TRUE,
                          use.J = FALSE,
                          export.adj.matrix = NULL,
                          export.graph = NULL,
                          # Deprecated arguments
                          dist_type = NULL,
                          dist_mat = NULL,
                          gap_open = NULL,
                          gap_extend = NULL,
                          exportAdjMatrix = NULL,
                          exportGraph = NULL) {

  # Handle deprecated arguments
  dist.type <- .deprecate_arg(dist_type, dist.type, "dist_type", "dist.type",
                              "clonalCluster", default = "levenshtein")
  dist.mat <- .deprecate_arg(dist_mat, dist.mat, "dist_mat", "dist.mat",
                             "clonalCluster", default = "BLOSUM80")
  gap.open <- .deprecate_arg(gap_open, gap.open, "gap_open", "gap.open",
                             "clonalCluster", default = -10)
  gap.extend <- .deprecate_arg(gap_extend, gap.extend, "gap_extend", "gap.extend",
                               "clonalCluster", default = -1)
  export.graph <- .deprecate_arg(exportGraph, export.graph, "exportGraph", "export.graph",
                                 "clonalCluster", default = FALSE)
  export.adj.matrix <- .deprecate_arg(exportAdjMatrix, export.adj.matrix, "exportAdjMatrix", "export.adj.matrix",
                                      "clonalCluster", default = FALSE)
  
  if (export.graph && export.adj.matrix) {
    stop("Please set only one of `export.graph` or `export.adj.matrix` to TRUE.")
  }
  
  #Prepping any single-cell object
  is.list <- inherits(input.data, "list")
  if(!is.list && !.is.seurat.or.se.object(input.data)) {
    input.data <- .checkList(input.data)
  }
  if (chain == "both") {
    # "both" means the two chains of whichever receptor is present: TRA/TRB for
    # T cells, Heavy/Light for B cells. getIR() parses the V(D)J(C) block by
    # chain name, so handing it TCR names for BCR data shifts the gene columns:
    # the heavy slot reads the D gene as `j` and the light slot reads the C gene
    # as `j`, which silently corrupts `use.J` filtering.
    chains_to_get <- if (.detectReceptor(input.data) == "B") {
      c("Heavy", "Light")
    } else {
      c("TRA", "TRB")
    }
    chain_data <- lapply(chains_to_get, function(x) {
      getIR(input.data, chains = x, sequence.type = sequence, group.by = group.by)
    })
  } else {
    chain_data <- getIR(input.data, chains = .chainConverter(chain), sequence.type = sequence, group.by = group.by)
    chain_data <- list(chain_data)
  }
  
  # Flatten the list if nested and get all barcodes
  if (is.list(chain_data[[1]]) && !is.data.frame(chain_data[[1]])) {
    chain_data <- unlist(chain_data, recursive = FALSE)
  }
  all_barcodes <- unique(do.call(rbind, chain_data)[["barcode"]])
  
  # Connected-components clustering only needs connectivity, so the cheaper
  # "star" edge expansion is exact there; community methods need the full edge
  # multiplicity, so they use "clique".
  expand <- if (identical(cluster.method, "components")) "star" else "clique"

  # Apply the network function to each data frame and combine into one edge list.
  # The engine signals per call (e.g. the Hamming equal-length notice, a
  # message() in current immApex and a warning() in older ones), which becomes
  # one copy per chain per group.by level. Collect both kinds, then re-emit each
  # distinct text once so a real problem is not buried in repeats.
  held <- list(warning = character(0), message = character(0))
  result_list <- withCallingHandlers(
    lapply(chain_data, function(y) {
      y <- y[!is.na(y[,1]),]
      .buildNetwork(y, use.V, use.J, threshold,
                    dist.type, dist.mat, normalize,
                    gap.open, gap.extend, expand)
    }),
    warning = function(w) {
      held$warning <<- c(held$warning, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      held$message <<- c(held$message, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  for (txt in unique(held$warning)) warning(txt, call. = FALSE)
  for (txt in unique(held$message)) message(sub("\n$", "", txt))
  full_edge_list <- do.call(rbind, result_list)
  
  # if group.by is invoked
  if (!is.null(group.by) && !is.null(names(chain_data))) {
    row_counts <- sapply(result_list, NROW)
    full_edge_list$group <- rep(names(chain_data), times = row_counts)
  }
  
  # pseudocount for same sequences
  full_edge_list$dist[full_edge_list$dist == 0] <- 1e-6
  
  # Add filtering
  if(threshold >= 1) {
    full_edge_list <- full_edge_list[full_edge_list$dist >= 1 & full_edge_list$dist <= threshold,]
  } else {
    max_ned <- 1-threshold
    full_edge_list <- full_edge_list[full_edge_list$dist <= max_ned, ]
  }
  full_edge_list <- unique(full_edge_list) 
  
  if (nrow(full_edge_list) == 0) {
    warning("No clusters found, returning original data.")
    return(input.data)
  }
  
  # Create the graph object
  full_g <- igraph::graph_from_data_frame(full_edge_list, 
                                          directed = FALSE)
  igraph::E(full_g)$weight <- full_edge_list$dist
  if (!is.null(group.by)) {
    vertex_map_df <- data.frame(
      name = c(full_edge_list$from, full_edge_list$to),
      group = rep(full_edge_list[["group"]], 2))
    unique_vertex_map <- vertex_map_df[!duplicated(vertex_map_df$name), ]
    vertex_group_lookup <- setNames(unique_vertex_map$group, unique_vertex_map$name)
    graph_vertex_names <- igraph::V(full_g)$name
    igraph::V(full_g)$group <- vertex_group_lookup[graph_vertex_names]
  }
  
  # Perform clustering using the specified method
  clusters <- .clusterGraph(cluster.method, full_g)
  cluster_membership <- igraph::membership(clusters)
  
  # Rename clusters sequentially with the specified prefix
  unique_clusters <- unique(cluster_membership)
  cluster_names <- paste0(cluster.prefix, seq_along(unique_clusters))
  name_map <- setNames(cluster_names, unique_clusters)
  renamed_membership <- name_map[as.character(cluster_membership)]
  names(renamed_membership) <- names(cluster_membership)
  
  # Add vertex communities 
  igraph::V(full_g)$cluster <- renamed_membership
  
  # Returning Graph
  if(export.graph) {
    # Adding additional vertex information if graph being exported
    vertex_df <- igraph::as_data_frame(full_g, what = "vertices")
    colnames(vertex_df)[1] <- "barcode" 
    full_meta_long <- do.call(rbind, chain_data)
    
    if (nrow(full_meta_long) > 0) {
      meta_to_process <- unique(full_meta_long[, c("barcode", "cdr3_aa", "v", "j")])
      meta_indexed <- meta_to_process %>% 
        group_by(barcode) %>%
        mutate(chain_num = dplyr::row_number()) %>%
        ungroup() %>%
        as.data.frame()
      original_meta <- reshape(
        data = meta_indexed, 
        idvar = "barcode",
        timevar = "chain_num",
        direction = "wide",
        v.names = c("cdr3_aa", "v", "j"),
        sep = "" )
      # Loop to add more vertex info
      vertex_names <- V(full_g)$name
      match_indices <- match(vertex_names, original_meta$barcode)
      for (col in colnames(original_meta)) {
        if (col == "barcode") next 
        full_g <- set_vertex_attr(
          graph = full_g,
          name = col,
          value = original_meta[[col]][match_indices])
      }
    }
    return(full_g)
  }
  
  # Returning Adjacency Matrix
  if (export.adj.matrix) {
    adj_from_graph <- igraph::as_adjacency_matrix(full_g, 
                                                  sparse = TRUE,
                                                  attr = "weight")
    barcodes_in_graph <- rownames(adj_from_graph)
    
    adjacency_matrix <- Matrix::sparseMatrix(
      i = {}, j = {}, x = 0.0,
      dims = c(length(all_barcodes), length(all_barcodes)),
      dimnames = list(all_barcodes, all_barcodes)
    )
    adjacency_matrix[barcodes_in_graph, barcodes_in_graph] <- adj_from_graph
    return(adjacency_matrix)
  }
  
  # Attaching to input.data
  bound <- igraph::as_data_frame(full_g, what = "vertices")
  # Vertex attributes come back as name, [group], cluster. When group.by is set
  # the second column is the grouping variable, so select the cluster column by
  # name instead of by position and drop the helper group column; indexing
  # position 2 previously wrote the group labels into the .Cluster column.
  cluster.col <- ifelse(chain == "both", "Multi.Cluster", paste0(chain, ".Cluster"))
  bound <- bound[, c("name", "cluster"), drop = FALSE]
  colnames(bound) <- c("name", cluster.col)
  
  #Adding to potential single-cell object
  if(.is.seurat.or.se.object(input.data)) {
    PreMeta <- bound[,-1, drop = FALSE]
    if (.is.seurat.object(input.data)) { 
      col.name <- names(PreMeta) %||% colnames(PreMeta)
      input.data[[col.name]] <- PreMeta
    } else {
      combined_col_names <- unique(c(colnames(colData(sc.data)), colnames(PreMeta)))
      full_data <- merge(colData(sc.data), PreMeta[rownames, , drop = FALSE], by = "row.names", all.x = TRUE)
      rownames(full_data) <- full_data[, 1]
      full_data  <- full_data[, -1]
      colData(sc.data) <- DataFrame(full_data[, combined_col_names])
    }
  } else {
    colnames(bound)[1] <- "barcode"
    input.data <- lapply(input.data, function(df) {
      dplyr::left_join(df, bound, by = "barcode")
    })
  }
  return(input.data)
}  

# Handling clustering options
#' @importFrom igraph cluster_fast_greedy cluster_walktrap cluster_louvain 
#' cluster_leiden cluster_spinglass cluster_edge_betweenness components
.clusterGraph <- function(clustering.method, graph, ...) {
  if (!igraph::is_igraph(graph)) {
    stop("The 'graph' argument must be a valid 'igraph' object.", call. = FALSE)
  }
  
  if (!is.character(clustering.method) || length(clustering.method) != 1) {
    stop("The 'clustering.method' argument must be a single character string.", call. = FALSE)
  }
  result <- switch(
    clustering.method,
    "fast_greedy" = igraph::cluster_fast_greedy(graph, ...),
    "walktrap" = igraph::cluster_walktrap(graph, ...),
    "louvain" = igraph::cluster_louvain(graph, ...),
    "leiden" = igraph::cluster_leiden(graph, ...),
    "spinglass" = igraph::cluster_spinglass(graph, ...),
    "edge_betweenness" = igraph::cluster_edge_betweenness(graph, ...),
    "components" = igraph::components(graph, mode = "weak", ...),
    
    # Default case: Handle unsupported methods
    stop(
      paste0(
        "Unsupported clustering.method: '", clustering.method, "'. ",
        "Supported methods are: 'fast_greedy', 'walktrap', 'louvain', ",
        "'leiden', 'spinglass', 'edge_betweenness', and 'components'."
      ),
      call. = FALSE
    )
  )
  
  return(result)
}

#' @importFrom immApex buildNetwork
.buildNetwork <- function(df, use.V, use.J, threshold,
                          dist.type, dist.mat, normalize,
                          gap.open, gap.extend, expand = "clique") {
  # A normalized `threshold` (e.g. 0.85) means "keep >= 85% similar", i.e. keep
  # normalized distance <= 1 - threshold. buildNetwork interprets a <1 threshold
  # as the max normalized distance, so pass the tight bound (1 - threshold)
  # instead of the loose 0.85. This lets the engine prune early (and the
  # downstream `dist <= 1 - threshold` filter below stays a no-op) rather than
  # returning a huge edge list that is then discarded. Raw thresholds (>= 1) are
  # passed through unchanged.
  bn_threshold <- if (threshold < 1) 1 - threshold else threshold

  # `cdr3_aa` is getIR()'s fixed output name for the CDR3 column, not a claim
  # about the alphabet: it carries nucleotides when the caller asked for
  # `sequence = "nt"` (getIR pulls CTnt instead of CTaa). There is no cdr3_nt
  # column to point at here.
  edge_list <- buildNetwork(df,
                            seq_col   = "cdr3_aa",
                            v_col     = "v",
                            j_col     = "j",
                            filter.v  = use.V,
                            filter.j  = use.J,
                            ids       = df[["barcode"]],
                            threshold = bn_threshold,
                            dist_type = dist.type,
                            dist_mat  = dist.mat,
                            normalize = normalize,
                            gap_open  = gap.open,
                            gap_extend= gap.extend,
                            expand    = expand)

  return(edge_list)
}