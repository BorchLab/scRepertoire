# Adding Global Variables
# data('v_gene','j_gene', 'c_gene', 'd_gene')
# note that currently the Rcpp internals have hardcoded column names so
# if some breaking change here is made, the Rcpp code will need to be updated,
# or functions need to be adjusted to intake expected column names that
# uses these variables
utils::globalVariables(c("v_gene", "j_gene", "c_gene", "d_gene", "chain"))

heavy_lines <- c("IGH", "cdr3_aa1", "cdr3_nt1", "vgene1")
light_lines <- c("IGLC", "cdr3_aa2", "cdr3_nt2", "vgene2")
l_lines <- c("IGLct", "cdr3", "cdr3_nt", "v_gene")
k_lines <- c("IGKct", "cdr3", "cdr3_nt", "v_gene")
h_lines <- c("IGHct", "cdr3", "cdr3_nt", "v_gene")
tcr1_lines <- c("TCR1", "cdr3_aa1", "cdr3_nt1")
tcr2_lines <- c("TCR2", "cdr3_aa2", "cdr3_nt2")
data1_lines <- c("TCR1", "cdr3", "cdr3_nt")
data2_lines <- c("TCR2", "cdr3", "cdr3_nt")
CT_lines <- c("CTgene", "CTnt", "CTaa", "CTstrict")

utils::globalVariables(c(
    "heavy_lines", "light_lines", "l_lines", "k_lines", "h_lines", "tcr1_lines",
    "tcr2_lines", "data1_lines", "data2_lines", "CT_lines"
))

#' @title Combine T Cell Receptor Contig Data
#'
#' @description This function consolidates a list of TCR sequencing results to
#' the level of  the individual cell barcodes. Using the `samples` and
#' `ID` parameters, the function will add the strings as prefixes to
#' prevent issues with repeated  barcodes. The resulting new barcodes will
#' need to match the Seurat or SCE object in order to use,
#' [combineExpression()]. Several levels of filtering exist -
#' `remove.na`, `remove.multi`, or `filter.multi` are parameters
#' that control how the function deals with barcodes with multiple chains
#' recovered.
#'
#' @examples
#' combined <- combineTCR(contig_list,
#'                         samples = c("P17B", "P17L", "P18B", "P18L",
#'                                     "P19B","P19L", "P20B", "P20L"))
#'
#' @param input.data List of filtered contig annotations or
#' outputs from [loadContigs()].
#' @param samples The labels of samples (recommended).
#' @param ID The additional sample labeling (optional).
#' @param remove.na This will remove any chain without values.
#' @param remove.multi This will remove barcodes with greater than 2 chains.
#' @param filter.multi This option will allow for the selection of the 2
#' corresponding chains with the highest expression for a single barcode.
#' @param filter.nonproductive This option will allow for the removal of
#' nonproductive chains if the variable exists in the contig data. Default
#' is set to TRUE to remove nonproductive contigs.
#' @param retain.sequences Optionally retain full-length sequences as additional
#' per-chain columns without altering `CTaa`/`CTnt`/`CTgene`/`CTstrict`. `FALSE`
#' (default) keeps the legacy output. `TRUE` retains the standardized
#' `sequence` and `sequence_aa` columns (when present from [loadContigs()]) as
#' `sequence_nt1`/`sequence_nt2` and `sequence_aa1`/`sequence_aa2`. A character
#' vector retains specific columns, e.g.
#' `c("sequence_alignment", "germline_alignment")` for downstream dowser use.
#' @param removeNA \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `remove.na` instead.
#' @param removeMulti \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `remove.multi` instead.
#' @param filterMulti \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `filter.multi` instead.
#' @param filterNonproductive \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `filter.nonproductive` instead.
#'
#' @export
#' @concept Loading_and_Processing_Contigs
#' @return List of clones for individual cell barcodes
#'
combineTCR <- function(input.data,
                       samples = NULL,
                       ID = NULL,
                       remove.na = NULL,
                       remove.multi = NULL,
                       filter.multi = NULL,
                       filter.nonproductive = NULL,
                       retain.sequences = FALSE,
                       # Deprecated arguments
                       removeNA = NULL,
                       removeMulti = NULL,
                       filterMulti = NULL,
                       filterNonproductive = NULL) {

    # Handle deprecated arguments
    remove.na <- .deprecate_arg(removeNA, remove.na, "removeNA", "remove.na",
                                "combineTCR", default = FALSE)
    remove.multi <- .deprecate_arg(removeMulti, remove.multi, "removeMulti", "remove.multi",
                                   "combineTCR", default = FALSE)
    filter.multi <- .deprecate_arg(filterMulti, filter.multi, "filterMulti", "filter.multi",
                                   "combineTCR", default = FALSE)
    filter.nonproductive <- .deprecate_arg(filterNonproductive, filter.nonproductive,
                                           "filterNonproductive", "filter.nonproductive",
                                           "combineTCR", default = TRUE)

    input.data <- .checkList(input.data)
    input.data <- .checkContigs(input.data)
    out <- NULL
    final <- NULL
    for (i in seq_along(input.data)) {
        if(c("chain") %in% colnames(input.data[[i]])) {
          input.data[[i]] <- subset(input.data[[i]], chain != "Multi")
        }
        if(c("productive") %in% colnames(input.data[[i]]) & filter.nonproductive) {
          input.data[[i]] <- subset(input.data[[i]], productive %in% c(TRUE, "TRUE", "True", "true"))
        }
        input.data[[i]]$sample <- samples[i]
        input.data[[i]]$ID <- ID[i]
        if (filter.multi) {
          input.data[[i]] <- .filteringMulti(input.data[[i]])
        }
    }
    #Prevents error caused by list containing elements with 0 rows
    blank.rows <- which(unlist(lapply(input.data, nrow)) == 0)
    if(length(blank.rows) > 0) {
      input.data <- input.data[-blank.rows]
      if(!is.null(samples)) {
        samples <- samples[-blank.rows]
      }
      if(!is.null(ID)) {
        ID <- ID[-blank.rows]
      }
    }
    if (!is.null(samples)) {
      out <- .modifyBarcodes(input.data, samples, ID)
    } else {
      out <- input.data
    }
    for (i in seq_along(out)) {
        data2 <- .makeGenes(cellType = "T", out[[i]])
        retain_cols <- .resolveRetainCols(retain.sequences, colnames(data2))
        Con.df <- .constructConDfAndParseTCR(data2, retain_cols)
        Con.df <- .assignCT(cellType = "T", Con.df)
        Con.df[Con.df == "NA_NA" | Con.df == "NA;NA_NA;NA"] <- NA
        Con.df <- .renameSeqCols(Con.df, retain_cols)
        data3 <- merge(data2[,-which(names(data2) %in% c("TCR1","TCR2"))],
            Con.df, by = "barcode")

        columns_to_include <- c("barcode")
        # Conditionally add columns based on user input
        if (!is.null(samples)) {
          columns_to_include <- c(columns_to_include, "sample")
        }
        if (!is.null(ID)) {
          columns_to_include <- c(columns_to_include, "ID")
        }

        # Add TCR and CT lines which are presumably always needed
        columns_to_include <- c(columns_to_include, tcr1_lines, tcr2_lines, CT_lines)
        # Optionally retain full-length sequence columns (additive, never CT*)
        if (length(retain_cols) > 0) {
          columns_to_include <- c(columns_to_include, .seqOutCols(retain_cols))
        }

        # Subset the data frame based on the dynamically built list of columns
        data3 <- data3[, columns_to_include]

        final[[i]] <- data3
    }
    name_vector <- character(length(samples))
    for (i in seq_along(samples)) {
        if (!is.null(samples) && !is.null(ID)) {
            curr <- paste(samples[i], "_", ID[i], sep="")
        } else if (!is.null(samples) & is.null(ID)) {
            curr <- paste(samples[i], sep="")
        }
        name_vector[i] <- curr
    }
    names(final) <- name_vector
    for (i in seq_along(final)){
      final[[i]]<-final[[i]][!duplicated(final[[i]]$barcode),]
      final[[i]]<-final[[i]][rowSums(is.na(final[[i]])) < 10, ]
      final[[i]][final[[i]] == "NA"] <- NA
    }
    if (remove.na) {
      final <- .removingNA(final)
    }
    if (remove.multi) {
      final <- .removingMulti(final)
    }
    #Adding list element names to output if samples NULL
    if(is.null(samples)) {
      names(final) <- paste0("S", seq_len(length(final)))
    }
    final
}

#' Combine B Cell Receptor Contig Data
#'
#' This function consolidates a list of BCR sequencing results to the level
#' of the individual cell barcodes. Using the samples and ID parameters,
#' the function will add the strings as prefixes to prevent issues with
#' repeated barcodes. The resulting new barcodes will need to match the
#' Seurat or SCE object in order to use, [combineExpression()]. Unlike
#' [combineTCR()], combineBCR produces a column `CTstrict` based on the
#' edit distance clustering from [clonalCluster()]. The `CTstrict` column
#' is formatted as `Heavy_Light` (underscore-separated) for downstream
#' compatibility. Connected clones are labeled with `cluster.X`, while
#' unconnected clones (singlets) are labeled with the V gene and CDR3
#' sequence (e.g., `IGHV3-64.CAKSYS..._IGKV3-15.CQQYSN...`).
#'
#' @examples
#' # Data derived from the 10x Genomics intratumoral NSCLC B cells
#' BCR <- read.csv("https://www.borch.dev/uploads/contigs/b_contigs.csv")
#' combined <- combineBCR(BCR,
#'                        samples = "Patient1",
#'                        threshold = 0.85)
#'
#' @param input.data List of filtered contig annotations or outputs from
#' [loadContigs()].
#' @param samples A character vector of sample labels. Must be the same length
#' as the input list.
#' @param ID An optional character vector for additional sample identifiers.
#' @param call.related.clones Logical. If `TRUE`, uses `clonalCluster()` to
#' identify related clones based on sequence similarity. If `FALSE`, defines
#' clones by the exact V-gene and CDR3 amino acid sequence.
#' @param group.by The column header used for to group clones.
#' If (`NULL``), clusters will be calculated across samples.
#' @param threshold The similarity threshold passed to `clonalCluster()` if
#' `call.related.clones = TRUE`. See `?clonalCluster` for details.
#' @param chain The chain to use for clustering when `call.related.clones = TRUE`.
#' Passed to `clonalCluster()`. Default is `"IGH"`.
#' @param sequence The sequence type (`"nt"` or `"aa"`) to use for clustering.
#' Passed to `clonalCluster()`, which reads `CTnt` or `CTaa` accordingly.
#' Default is `"nt"`.
#' @param dist.type The distance metric to use. Options: `"levenshtein"` (default),
#' `"hamming"`, `"damerau"`, `"nw"` (Needleman-Wunsch), or `"sw"` (Smith-Waterman).
#' `"hamming"` only compares sequences of equal length; unequal-length pairs are
#' given the maximum distance and never cluster, which emits a notice. See
#' [clonalCluster()] for details.
#' @param dist.mat The substitution matrix to use for alignment-based metrics
#' (`"nw"` or `"sw"`). Options include `"BLOSUM62"`, `"PAM30"`, etc.
#' @param normalize Method for normalizing distances. Options: `"none"` (default),
#' `"maxlen"`, or `"length"`.
#' @param gap.open Penalty for opening a gap in alignment metrics (default: -10).
#' @param gap.extend Penalty for extending a gap in alignment metrics (default: -1).
#' @param use.V Logical. If `TRUE`, sequences must share the same V gene to be
#' clustered together.
#' @param use.J Logical. If `TRUE`, sequences must share the same J gene to be
#' clustered together.
#' @param cluster.method The clustering algorithm to use. Defaults to `"components"`,
#' which finds connected subgraphs.
#' @param remove.multi Logical. If `TRUE`, removes cells that have more than
#' one distinct heavy or light chain after processing.
#' @param filter.multi Logical. If `TRUE`, filters multi-chain cells to retain
#' only the most abundant IGH and IGL/IGK chains.
#' @param remove.na This will remove any chain without values.
#' @param filter.nonproductive Logical. If `TRUE`, removes non-productive contigs
#' from the analysis.
#' @param retain.sequences Optionally retain full-length sequences as additional
#' per-chain columns without altering `CTaa`/`CTnt`/`CTgene`/`CTstrict`. `FALSE`
#' (default) keeps the legacy output. `TRUE` retains the standardized
#' `sequence` and `sequence_aa` columns (when present from [loadContigs()]) as
#' `sequence_nt1`/`sequence_nt2` and `sequence_aa1`/`sequence_aa2`. A character
#' vector retains specific columns, e.g.
#' `c("sequence_alignment", "germline_alignment")` for downstream dowser use.
#' @param removeNA \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `remove.na` instead.
#' @param removeMulti \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `remove.multi` instead.
#' @param filterMulti \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `filter.multi` instead.
#' @param filterNonproductive \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `filter.nonproductive` instead.
#' @param dist_type \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `dist.type` instead.
#' @param dist_mat \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `dist.mat` instead.
#' @param gap_open \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `gap.open` instead.
#' @param gap_extend \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `gap.extend` instead.
#'
#' @export
#' @concept Loading_and_Processing_Contigs
#' @return A list of data frames, where each data frame represents a sample.
#' Each row corresponds to a unique cell barcode, with columns detailing the
#' BCR chains and the assigned clone ID.
combineBCR <- function(input.data,
                       samples = NULL,
                       ID = NULL,
                       chain = "IGH",
                       sequence = "nt",
                       dist.type = NULL,
                       dist.mat = NULL,
                       normalize = "length",
                       gap.open = NULL,
                       gap.extend = NULL,
                       call.related.clones = TRUE,
                       group.by = NULL,
                       threshold = 0.85,
                       cluster.method = "components",
                       use.V = TRUE,
                       use.J = TRUE,
                       remove.na = NULL,
                       remove.multi = NULL,
                       filter.multi = NULL,
                       filter.nonproductive = NULL,
                       retain.sequences = FALSE,
                       # Deprecated arguments
                       removeNA = NULL,
                       removeMulti = NULL,
                       filterMulti = NULL,
                       filterNonproductive = NULL,
                       dist_type = NULL,
                       dist_mat = NULL,
                       gap_open = NULL,
                       gap_extend = NULL) {

  # Handle deprecated arguments
  remove.na <- .deprecate_arg(removeNA, remove.na, "removeNA", "remove.na",
                              "combineBCR", default = FALSE)
  remove.multi <- .deprecate_arg(removeMulti, remove.multi, "removeMulti", "remove.multi",
                                 "combineBCR", default = FALSE)
  filter.multi <- .deprecate_arg(filterMulti, filter.multi, "filterMulti", "filter.multi",
                                 "combineBCR", default = TRUE)
  filter.nonproductive <- .deprecate_arg(filterNonproductive, filter.nonproductive,
                                         "filterNonproductive", "filter.nonproductive",
                                         "combineBCR", default = TRUE)
  dist.type <- .deprecate_arg(dist_type, dist.type, "dist_type", "dist.type",
                              "combineBCR", default = "levenshtein")
  dist.mat <- .deprecate_arg(dist_mat, dist.mat, "dist_mat", "dist.mat",
                             "combineBCR", default = "BLOSUM80")
  gap.open <- .deprecate_arg(gap_open, gap.open, "gap_open", "gap.open",
                             "combineBCR", default = -10)
  gap.extend <- .deprecate_arg(gap_extend, gap.extend, "gap_extend", "gap.extend",
                               "combineBCR", default = -1)

  # Resolve optional full-length sequence retention against the input columns.
  # Empty when retain.sequences = FALSE, keeping the default path untouched.
  # Use the intersection of columns across ALL samples so a column missing from
  # any one sample is never selected (a per-sample mismatch otherwise crashes
  # .parseBCR with "undefined columns selected" or silently drops the column).
  common_cols <- Reduce(intersect, lapply(.checkList(input.data), colnames))
  retain_cols <- .resolveRetainCols(retain.sequences, common_cols)
  seq_bases <- .seqBases(retain_cols)

  # Initial Contig Processing and Filtering
  processed_list <- input.data %>%
    .checkList() %>%
    .checkContigs() %>%
    unname() %>%
    purrr::imap(function(x, i) {
      x <- subset(x, chain %in% c("IGH", "IGK", "IGL"))
      if (!is.null(ID)) x$ID <- ID[i]
      if (filter.nonproductive && "productive" %in% colnames(x)) {
        x <- subset(x, tolower(productive) == "true")
      }
      if (filter.multi) {
        # Keep IGH / IGK / IGL info in save_chain
        x$save_chain <- x$chain
        # Collapse IGK and IGL chains
        x$chain <- ifelse(x$chain == "IGH", "IGH", "IGLC")
        x <- .filteringMulti(x)
        # Get back IGK / IGL distinction
        x$chain <- x$save_chain
        x$save_chain <- NULL
      }
      x
    }) %>%
    # Add sample/ID prefixes
    (function(x) {
      if (!is.null(samples)) {
        .modifyBarcodes(x, samples, ID)
      } else { 
        x
      }
    }) %>%
    # Reshape data to one row per barcode with columns for each chain
    lapply(function(x) {
      data2 <- data.frame(x)
      data2 <- .makeGenes(cellType = "B", data2)
      unique_df <- unique(data2$barcode)
      # Widen the per-cell target/source vectors when retaining sequences so the
      # full-length columns ride the same .parseBCR pass as the gene/CDR3 data.
      heavy_t <- c(heavy_lines, paste0(seq_bases, "1"))
      light_t <- c(light_lines, paste0(seq_bases, "2"))
      h_s <- c(h_lines, retain_cols)
      k_s <- c(k_lines, retain_cols)
      l_s <- c(l_lines, retain_cols)
      Con.df <- data.frame(matrix(NA, length(unique_df), 1 + length(heavy_t) + length(light_t)))
      colnames(Con.df) <- c("barcode", heavy_t, light_t)
      Con.df$barcode <- unique_df
      Con.df <- .parseBCR(Con.df, unique_df, data2, heavy_t, light_t, h_s, k_s, l_s)
      Con.df <- .assignCT(cellType = "B", Con.df)
      if(!is.null(group.by)) { #retain group.by variable for clustering
        barcode_groups <- data2[!duplicated(data2$barcode), c("barcode", group.by)]
        Con.df[[group.by]] <- barcode_groups[[group.by]][match(Con.df$barcode, barcode_groups$barcode)]
      }
      Con.df %>% 
        mutate(length1 = nchar(cdr3_nt1)) %>%
        mutate(length2 = nchar(cdr3_nt2))
    })
  
  # Getting CTstrict based on clusters
  if (call.related.clones) {
    clusters <- clonalCluster(processed_list,
                              sequence = sequence,
                              chain = chain,
                              threshold = threshold,
                              group.by = group.by,
                              use.V = use.V,
                              use.J = use.J,
                              cluster.method = cluster.method,
                              dist.type = dist.type,
                              dist.mat = dist.mat,
                              normalize = normalize,
                              gap.open = gap.open,
                              gap.extend = gap.extend)
  }
  
  # Defining element names for the final output
  list_names <- if (!is.null(samples)) {
    if (is.null(ID)) samples else paste0(samples, "_", ID)
  } else {
    paste0("S", seq_along(processed_list))
  }
  
  final_list <- purrr::map2(processed_list, seq_along(processed_list), function(df, i) {
    # Assigning CTstrict
    if (call.related.clones) {
      # Get the cluster column from clonalCluster output by name. clonalCluster
      # returns the input untouched when no edges pass the threshold, so the
      # column can legitimately be absent; positional indexing would silently
      # pick up whatever column happened to be last.
      cluster.col <- ifelse(chain == "both", "Multi.Cluster",
                            paste0(chain, ".Cluster"))
      cluster_col <- if (cluster.col %in% colnames(clusters[[i]])) {
        clusters[[i]][[cluster.col]]
      } else {
        rep(NA_character_, nrow(df))
      }
      
      # ========== CTstrict FORMATTING LOGIC ==========
      seq_col <- ifelse(sequence == "aa", "cdr3_aa", "cdr3_nt")
      heavy_seq_col <- paste0(seq_col, "1")
      light_seq_col <- paste0(seq_col, "2")
      
      # Create unique identifiers for each chain (vgene.sequence format)
      heavy_unique <- ifelse(
        !is.na(df[, "vgene1"]) & !is.na(df[, heavy_seq_col]),
        paste0(df[, "vgene1"], ".", df[, heavy_seq_col]),
        "NA"
      )
      
      light_unique <- ifelse(
        !is.na(df[, "vgene2"]) & !is.na(df[, light_seq_col]),
        paste0(df[, "vgene2"], ".", df[, light_seq_col]),
        "NA"
      )
      
      if (chain == "both") {
        # Both chains in same network - use cluster ID for both parts
        heavy_part <- ifelse(is.na(cluster_col), heavy_unique, cluster_col)
        light_part <- ifelse(is.na(cluster_col), light_unique, cluster_col)
        
      } else if (chain == "IGH") {
        # Only heavy chain clustered
        heavy_part <- ifelse(is.na(cluster_col), heavy_unique, cluster_col)
        light_part <- light_unique  # Light chain always uses unique ID
        
      } else if (chain %in% c("IGL", "IGK", "Light")) {
        # Only light chain clustered
        heavy_part <- heavy_unique  # Heavy chain always uses unique ID
        light_part <- ifelse(is.na(cluster_col), light_unique, cluster_col)
        
      } else {
        # Fallback for any other chain specification
        heavy_part <- heavy_unique
        light_part <- light_unique
      }
      
      # Combine into CTstrict with underscore separator
      df[, "CTstrict"] <- paste0(heavy_part, "_", light_part)
      
    } else {
      df[, "CTstrict"] <- paste0(df[, "vgene1"], ".", df[, "cdr3_aa1"], "_",
                                 df[, "vgene2"], ".", df[, "cdr3_aa2"])
    }
    # Adding samples/ID if applicable
    if (!is.null(samples)) df$sample <- samples[i]
    if (!is.null(ID)) df$ID <- ID[i]
    
    # Cleaning up the "NA"
    df[df == "NA_NA" | df == "NA.NA_NA.NA" | df == "NA;NA_NA;NA" | df == "NA"] <- NA
    
    # Select, reorder, and filter final columns
    col_selection <- c("barcode", "sample", "ID",
                       heavy_lines[c(1, 2, 3)], light_lines[c(1, 2, 3)], CT_lines)
    # Optionally retain full-length sequence columns (additive, never CT*)
    if (length(retain_cols) > 0) {
      col_selection <- c(col_selection,
                         paste0(seq_bases, "1"), paste0(seq_bases, "2"))
    }
    col_selection <- col_selection[col_selection %in% names(df)] # Keep only existing cols
    df <- df[, col_selection]
    # Tidy ';'-joined retained sequences: drop "NA" placeholder tokens left by
    # multi-contig chains where a contig lacked a full-length sequence.
    if (length(retain_cols) > 0) {
      for (sc in c(paste0(seq_bases, "1"), paste0(seq_bases, "2"))) {
        if (sc %in% names(df)) df[[sc]] <- .cleanSeqJoin(df[[sc]])
      }
    }
    df <- df[!duplicated(df$barcode), ]
    df <- df[rowSums(is.na(df)) < (ncol(df) - 1), ] 
    return(df)
  })
  
  # Set the names of the final list
  names(final_list) <- list_names
  
  # Final Optional Filtering
  if (remove.na) final_list <- .removingNA(final_list)
  if (remove.multi) final_list <- .removingMulti(final_list)
  
  return(final_list)
}
