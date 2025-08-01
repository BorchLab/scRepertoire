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
#' `removeNA`, `removeMulti`, or `filterMulti` are parameters
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
#' @param removeNA This will remove any chain without values.
#' @param removeMulti This will remove barcodes with greater than 2 chains.
#' @param filterMulti This option will allow for the selection of the 2
#' corresponding chains with the highest expression for a single barcode.
#' @param filterNonproductive This option will allow for the removal of
#' nonproductive chains if the variable exists in the contig data. Default
#' is set to TRUE to remove nonproductive contigs.
#'
#' @export
#' @concept Loading_and_Processing_Contigs
#' @return List of clones for individual cell barcodes
#'
combineTCR <- function(input.data,
                       samples = NULL,
                       ID = NULL,
                       removeNA = FALSE,
                       removeMulti = FALSE,
                       filterMulti = FALSE,
                       filterNonproductive = TRUE) {

    input.data <- .checkList(input.data)
    input.data <- .checkContigs(input.data)
    out <- NULL
    final <- NULL
    for (i in seq_along(input.data)) {
        if(c("chain") %in% colnames(input.data[[i]])) {
          input.data[[i]] <- subset(input.data[[i]], chain != "Multi")
        }
        if(c("productive") %in% colnames(input.data[[i]]) & filterNonproductive) {
          input.data[[i]] <- subset(input.data[[i]], productive %in% c(TRUE, "TRUE", "True", "true"))
        }
        input.data[[i]]$sample <- samples[i]
        input.data[[i]]$ID <- ID[i]
        if (filterMulti) {
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
        Con.df <- .constructConDfAndParseTCR(data2)
        Con.df <- .assignCT(cellType = "T", Con.df)
        Con.df[Con.df == "NA_NA" | Con.df == "NA;NA_NA;NA"] <- NA
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
    if (removeNA) {
      final <- .removingNA(final)
    }
    if (removeMulti) {
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
#' Seurat or SCE object in order to use, [combineExpression()].
#' Unlike [combineTCR()], combineBCR produces a column
#' `CTstrict` of an index of nucleotide sequence and the
#' corresponding V and J genes using [clonalCluster()]. 
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
#' Passed to `clonalCluster()`. Default is `"both"`.
#' @param sequence The sequence type (`"nt"` or `"aa"`) to use for clustering.
#' Passed to `clonalCluster()`. Default is `"nt"`.
#' @param use.V Logical. If `TRUE`, sequences must share the same V gene to be
#' clustered together.
#' @param use.J Logical. If `TRUE`, sequences must share the same J gene to be
#' clustered together.
#' @param cluster.method The clustering algorithm to use. Defaults to `"components"`, 
#' which finds connected subgraphs.
#' @param removeMulti Logical. If `TRUE`, removes cells that have more than
#' one distinct heavy or light chain after processing.
#' @param filterMulti Logical. If `TRUE`, filters multi-chain cells to retain
#' only the most abundant IGH and IGL/IGK chains.
#' @param removeNA This will remove any chain without values.
#' @param filterNonproductive Logical. If `TRUE`, removes non-productive contigs
#' from the analysis.
#'
#' @export
#' @concept Loading_and_Processing_Contigs
#' @return A list of data frames, where each data frame represents a sample.
#' Each row corresponds to a unique cell barcode, with columns detailing the
#' BCR chains and the assigned clone ID.
combineBCR <- function(input.data,
                       samples = NULL,
                       ID = NULL,
                       chain = "both",
                       sequence = "nt",
                       call.related.clones = TRUE,
                       group.by = NULL,
                       threshold = 0.85,
                       cluster.method = "components",
                       use.V = TRUE,
                       use.J = TRUE,
                       removeNA = FALSE,
                       removeMulti = FALSE,
                       filterMulti = TRUE,
                       filterNonproductive = TRUE) {

    # Initial Contig Processing and Filtering 
    processed_list <- input.data %>%
        .checkList() %>%
        .checkContigs() %>%
        unname() %>%
        purrr::imap(function(x, i) {
            x <- subset(x, chain %in% c("IGH", "IGK", "IGL"))
            if (!is.null(ID)) x$ID <- ID[i]
            if (filterNonproductive && "productive" %in% colnames(x)) {
                x <- subset(x, tolower(productive) == "true")
            }
            if (filterMulti) {
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
            } else { # https://github.com/BorchLab/scRepertoire/pull/450
                x
            }
        }) %>%
        # Reshape data to one row per barcode with columns for each chain
        lapply(function(x) {
            data2 <- data.frame(x)
            data2 <- .makeGenes(cellType = "B", data2)
            unique_df <- unique(data2$barcode)
            Con.df <- data.frame(matrix(NA, length(unique_df), 9))
            colnames(Con.df) <- c("barcode", heavy_lines, light_lines)
            Con.df$barcode <- unique_df
            Con.df <- .parseBCR(Con.df, unique_df, data2)
            Con.df <- .assignCT(cellType = "B", Con.df)
            Con.df %>% mutate(length1 = nchar(cdr3_nt1)) %>%
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
                                cluster.method = cluster.method)
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
        cluster_col <- clusters[[i]][, ncol(clusters[[i]])]
        df[, "CTstrict"] <- cluster_col
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
      col_selection <- col_selection[col_selection %in% names(df)] # Keep only existing cols
      df <- df[, col_selection]
      df <- df[!duplicated(df$barcode), ]
      df <- df[rowSums(is.na(df)) < (ncol(df) - 1), ] 
      return(df)
    })
    
    # Set the names of the final list
    names(final_list) <- list_names
    
    # Final Optional Filtering
    if (removeNA) final_list <- .removingNA(final_list)
    if (removeMulti) final_list <- .removingMulti(final_list)
    
    return(final_list)
}
