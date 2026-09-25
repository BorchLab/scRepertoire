# test script for clonalCluster.R 

library(testthat)
library(igraph)
library(SingleCellExperiment) 
library(Matrix)

# --- Global Setup ---
# Using the user-provided setup for general structure tests
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))

set.seed(42)
combined <- lapply(combined, function(x) {
  x <- x[sample(nrow(x), nrow(x) * .25),]
  x
})

# BCR contigs for the receptor-specific tests below
BCR_SOURCE <- read.csv("https://www.borch.dev/uploads/contigs/b_contigs.csv")

# --- Existing Tests (Preserved) ---

test_that("Basic functionality and default output structure", {
  clustered_list <- clonalCluster(combined[1:2])
  expect_length(clustered_list, length(combined[1:2]))
  expect_true(all(sapply(clustered_list, is.data.frame)))
  expect_true("TRB.Cluster" %in% names(clustered_list[[1]]))
  
  # Check for data integrity
  expect_equal(nrow(clustered_list[[1]]), nrow(combined[[1]]))
  expect_false(all(is.na(clustered_list[[1]]$TRB.Cluster)))
})

test_that("export.graph  = TRUE returns a valid igraph object", {
  graph_obj <- clonalCluster(combined[1:2], export.graph = TRUE)
  expect_s3_class(graph_obj, "igraph")
  expect_gt(vcount(graph_obj), 0)
  expect_gt(ecount(graph_obj), 0)
  expect_true("cluster" %in% igraph::vertex_attr_names(graph_obj))
  expect_true("weight" %in% igraph::edge_attr_names(graph_obj))
})

test_that("export.adj.matrix = TRUE returns a valid sparse matrix", {
  adj_matrix <- clonalCluster(combined[3:4], export.adj.matrix  = TRUE)
  all_barcodes <- unique(do.call(rbind, combined[3:4])[["barcode"]])
  num_barcodes <- length(all_barcodes)
  expect_s4_class(adj_matrix, "dgCMatrix")
  expect_equal(dim(adj_matrix), c(num_barcodes, num_barcodes))
  expect_equal(rownames(adj_matrix), all_barcodes)
})

test_that("chain parameter works correctly", {
  clustered_tra <- clonalCluster(combined[5:6], chain = "TRA")
  expect_true("TRA.Cluster" %in% names(clustered_tra[[1]]))
  expect_false("TRB.Cluster" %in% names(clustered_tra[[1]]))
  clustered_both <- clonalCluster(combined[5:6], chain = "both")
  expect_true("Multi.Cluster" %in% names(clustered_both[[1]]))
})

test_that("group.by parameter functions without error", {
  clustered_grouped <- clonalCluster(combined[1:2], 
                                     group.by = "sample")
  expect_true("TRB.Cluster" %in% names(clustered_grouped[[1]]))
  expect_false(all(is.na(clustered_grouped[[1]]$TRB.Cluster)))
})

test_that("Different `cluster.method` options work", {
  louvain_graph <- clonalCluster(combined[5:6], 
                                 cluster.method = "louvain", 
                                 export.graph  = TRUE)
  expect_s3_class(louvain_graph, "igraph")
  expect_true("cluster" %in% igraph::vertex_attr_names(louvain_graph))
})

test_that("Input validation and error handling", {
  expect_error(
    clonalCluster(combined, export.graph  = TRUE, export.adj.matrix = TRUE),
    "Please set only one of `export.graph` or `export.adj.matrix` to TRUE."
  )
  expect_error(
    clonalCluster(combined, cluster.method = "invalid_method"),
    "Unsupported clustering.method"
  )
})


test_that("Alignment (NW/SW) and Matrix selection", {
  # Create a small controlled dataset to ensure alignment logic runs
  toy_data <- combined[1]
  
  # Test Needleman-Wunsch with BLOSUM62
  res_nw <- clonalCluster(toy_data, 
                          dist.type = "nw", 
                          dist.mat = "BLOSUM62", 
                          threshold = 0.8) 
  expect_true("TRB.Cluster" %in% names(res_nw[[1]]))
  
  # Test Damerau (Transposition)
  res_dam <- clonalCluster(toy_data, dist.type = "damerau")
  expect_true("TRB.Cluster" %in% names(res_dam[[1]]))
})

test_that("Normalization parameters function correctly", {
  # normalize = "maxlen"
  res_max <- clonalCluster(combined[1], normalize = "maxlen", threshold = 0.1)
  expect_true("TRB.Cluster" %in% names(res_max[[1]]))
  
  # normalize = "length" (mean length)
  res_len <- clonalCluster(combined[1], normalize = "length", threshold = 0.1)
  expect_true("TRB.Cluster" %in% names(res_len[[1]]))
})
# --- Receptor detection, group.by column naming, and notice handling ---------

test_that(".detectReceptor distinguishes BCR from TCR input", {
  expect_equal(.detectReceptor(combined), "T")

  bcr <- combineBCR(BCR_SOURCE, samples = "P1", call.related.clones = FALSE)
  expect_equal(.detectReceptor(bcr), "B")

  # No CTgene to inspect falls back to the TCR chain names
  expect_equal(.detectReceptor(list(data.frame(barcode = "a"))), "T")
})

test_that("chain = 'both' resolves Heavy/Light for BCR, not TRA/TRB", {
  bcr <- combineBCR(BCR_SOURCE, samples = "P1", call.related.clones = FALSE)

  # TRA/TRB names shift the V(D)J parse on BCR data: the heavy slot picks up
  # the D gene as `j` and the light slot picks up the C gene as `j`.
  heavy_wrong <- immApex::getIR(bcr, chains = "TRA", sequence.type = "nt")
  heavy_right <- immApex::getIR(bcr, chains = "Heavy", sequence.type = "nt")
  expect_true(any(grepl("^IGHD", na.omit(heavy_wrong$j))))
  expect_true(all(grepl("^IGHJ", na.omit(heavy_right$j))))

  light_wrong <- immApex::getIR(bcr, chains = "TRB", sequence.type = "nt")
  light_right <- immApex::getIR(bcr, chains = "Light", sequence.type = "nt")
  expect_true(any(grepl("C$", na.omit(light_wrong$j))))
  expect_true(all(grepl("^IG[KL]J", na.omit(light_right$j))))

  # The J filter must therefore act on real J genes end to end
  res <- clonalCluster(bcr, chain = "both", sequence = "nt",
                       threshold = 0.85, use.V = TRUE, use.J = TRUE)
  expect_true("Multi.Cluster" %in% names(res[[1]]))
  expect_false(all(is.na(res[[1]]$Multi.Cluster)))
})

test_that("group.by writes cluster IDs, not group labels, into the .Cluster column", {
  res <- clonalCluster(combined[1:2], chain = "TRB", sequence = "aa",
                       threshold = 0.85, group.by = "sample")

  expect_true("TRB.Cluster" %in% names(res[[1]]))
  # The stale positional rename leaked the grouping variable into this column
  expect_false("cluster" %in% names(res[[1]]))
  ids <- na.omit(unlist(lapply(res, `[[`, "TRB.Cluster")))
  expect_true(all(grepl("^cluster\\.", ids)))
  expect_false(any(ids %in% unique(combined[[1]]$sample)))

  # Grouping is respected: a cluster never spans two groups
  bound <- do.call(rbind, lapply(res, function(x) x[, c("sample", "TRB.Cluster")]))
  bound <- bound[!is.na(bound$TRB.Cluster), ]
  expect_equal(max(tapply(bound$sample, bound$TRB.Cluster,
                          function(x) length(unique(x)))), 1L)
})

test_that("repeated engine notices are collapsed to one per distinct message", {
  bcr <- combineBCR(BCR_SOURCE, samples = "P1", call.related.clones = FALSE)
  bcr[[1]]$grp <- rep(c("a", "b"), length.out = nrow(bcr[[1]]))

  msgs <- character(0)
  withCallingHandlers(
    clonalCluster(bcr, chain = "both", sequence = "nt", dist.type = "hamming",
                  threshold = 0.9, group.by = "grp", use.V = TRUE, use.J = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  hamming <- grep("Hamming", msgs, value = TRUE)
  # 2 chains x 2 groups = 4 raw notices, surfaced once
  expect_lte(length(hamming), 1L)
})

test_that("sequence = 'nt' clusters on nucleotides, not amino acids", {
  bcr <- combineBCR(BCR_SOURCE, samples = "P1", call.related.clones = FALSE)

  # getIR names the CDR3 column cdr3_aa regardless of alphabet; confirm the
  # nt request actually pulls CTnt
  nt <- immApex::getIR(bcr, chains = "Heavy", sequence.type = "nt")
  aa <- immApex::getIR(bcr, chains = "Heavy", sequence.type = "aa")
  expect_true(all(grepl("^[ACGTN]+$", na.omit(nt$cdr3_aa))))
  expect_false(all(grepl("^[ACGTN]+$", na.omit(aa$cdr3_aa))))

  res_nt <- clonalCluster(bcr, chain = "IGH", sequence = "nt", threshold = 0.85)
  res_aa <- clonalCluster(bcr, chain = "IGH", sequence = "aa", threshold = 0.85)
  expect_false(identical(res_nt[[1]]$IGH.Cluster, res_aa[[1]]$IGH.Cluster))
})
