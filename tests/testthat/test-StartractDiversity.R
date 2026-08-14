# test script for StartracDiversity.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list, samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"))
scRep_example <- combineExpression(combined, scRep_example)
scRep_example$Patient <- substring(scRep_example$orig.ident, 1, 3)
scRep_example$Type <- substring(scRep_example$orig.ident, 4, 4)

# Minimal frame in the shape .calculateIndices()/.calculatePairwiseIndices()
# receive, for exercising the guards that real data rarely reaches
toy_processed <- data.frame(
  Cell_Name = paste0("cell", seq_len(4)),
  clone.id = c("a", "a", "b", "b"),
  patient = "P1",
  cluster = c("1", "2", "1", "2"),
  loc = c("B", "L", "B", "L"),
  stringsAsFactors = FALSE
)

test_that("Input validation and error handling work correctly", {
  expect_error(
    StartracDiversity(scRep_example, index = c("migr", "tran"), pairwise = "Type"),
    "Pairwise analysis can only be performed for a single index"
  )

  expect_error(
    StartracDiversity(scRep_example, index = "expa", pairwise = "Type"),
    "Pairwise analysis is only supported for 'migr' or 'tran' indices."
  )

  expect_error(
    StartracDiversity(scRep_example, type = "Type", index = "clonality"),
    "Please select 'expa', 'migr', and/or 'tran' for index."
  )

  # group.by is required when the object carries no orig.ident to fall back on
  no_ident <- scRep_example
  no_ident$orig.ident <- NULL
  expect_error(
    StartracDiversity(no_ident, type = "Type"),
    "Please select a group.by variable."
  )

  # Nothing left to calculate after dropping cells without clones
  no_clones <- scRep_example
  no_clones$CTstrict <- NA
  expect_warning(
    result <- StartracDiversity(no_clones, type = "Type", group.by = "Patient",
                                export.table = TRUE),
    "No data available for calculation"
  )
  expect_null(result)
})


test_that("group.by defaults to orig.ident and chain can be restricted", {
  default_group <- StartracDiversity(scRep_example,
                                     type = "Type",
                                     export.table = TRUE)
  expect_setequal(unique(default_group$group), unique(scRep_example$orig.ident))

  # Single-chain requests route through .offTheChain()
  single_chain <- StartracDiversity(scRep_example,
                                    type = "Type",
                                    group.by = "Patient",
                                    chain = "TRB",
                                    export.table = TRUE)
  expect_s3_class(single_chain, "data.frame")
  expect_true(all(c("group", "cluster", "expa", "migr", "tran") %in% names(single_chain)))
  expect_setequal(unique(single_chain$group), unique(scRep_example$Patient))
})


test_that("Output format and structure are correct", {
  # Returns a ggplot object by default
  plot_output <- StartracDiversity(scRep_example, 
                                   type = "Type", 
                                   group.by = "Patient")
  expect_s3_class(plot_output, "ggplot")
  
  # Returns a data.frame when export.table = TRUE
  table_output <- StartracDiversity(scRep_example, 
                                    type = "Type", 
                                    group.by = "Patient", 
                                    export.table = TRUE)
  expect_s3_class(table_output, "data.frame")
  
  # Check standard output columns
  expect_true(all(c("group", "cluster", "expa", "migr", "tran") %in% names(table_output)))
})


test_that("Helper functions for entropy calculate correctly", {
  # Test matrix
  mat <- matrix(c(1, 1, 0, 4, 2, 0, 3, 0, 0), nrow = 3, byrow = TRUE)
  
  expected_row_entropy <- c(1, 0.9182958, 0)
  expect_equal(.mrowEntropy(mat), expected_row_entropy, tolerance = 1e-6)

  expected_col_entropy <- c(1.405639, 0.9182958, NA)
  expect_equal(.mcolEntropy(mat), expected_col_entropy, tolerance = 1e-2)
})


test_that("Standard index calculations are mathematically correct", {
  results <- StartracDiversity(scRep_example, 
                               type = "Type", 
                               group.by = "Patient", 
                               export.table = TRUE)
  

  expect_equal(results$expa[results$cluster == "1"][1], 0, tolerance = 1e-4)
  expect_equal(results$expa[results$cluster == "2"][1], 0, tolerance = 1e-4)
  expect_equal(results$expa[results$cluster == "3"][1], 0, tolerance = 1e-4)
  
  expect_equal(results$migr[results$cluster == "1"][1], 0, tolerance = 1e-4)
  expect_equal(results$migr[results$cluster == "2"][1], 0, tolerance = 1e-4)
  expect_equal(results$migr[results$cluster == "3"][1], 0, tolerance = 1e-4)
  
  expect_equal(results$tran[results$cluster == "1"][1], 0, tolerance = 1e-4)
  expect_equal(results$tran[results$cluster == "2"][1], 0, tolerance = 1e-4)
  expect_equal(results$tran[results$cluster == "3"][1], 0.119958, tolerance = 1e-4)
})

test_that("Pairwise calculations are correct", {
  # Test 1: Pairwise migration (index = "migr", pairwise = "Type")
  pairwise_migr_results <- StartracDiversity(scRep_example,
                                             type = "Type",
                                             group.by = "Patient",
                                             index = "migr",
                                             pairwise = "Type",
                                             export.table = TRUE)
  expect_s3_class(pairwise_migr_results, "data.frame")
  expect_true(all(c("group", "cluster", "value", "comparison") %in% names(pairwise_migr_results)))
  expect_true(all(grepl("\\w vs \\w", pairwise_migr_results$comparison)))
  expect_equal(unique(pairwise_migr_results$comparison), "B vs L")
  
  
  # Test 2: Pairwise transition (index = "tran", pairwise = "cluster")
  pairwise_tran_results <- StartracDiversity(scRep_example,
                                             type = "Type",
                                             group.by = "Patient",
                                             index = "tran",
                                             pairwise = "cluster",
                                             export.table = TRUE)
  expect_s3_class(pairwise_tran_results, "data.frame")
  expect_true(all(c("group", "cluster", "partner", "value", "comparison") %in% names(pairwise_tran_results)))
  expect_true(all(grepl("\\w vs \\w", pairwise_tran_results$comparison)))
  expect_true(is.numeric(pairwise_tran_results$value))

  # Pairs are unordered, so a given pair must carry one label in every group
  pair_key <- vapply(strsplit(pairwise_tran_results$comparison, " vs "),
                     function(x) paste(sort(x), collapse = "|"), character(1))
  expect_equal(length(unique(pairwise_tran_results$comparison)),
               length(unique(pair_key)))

  # The anchor and its partner are the two members of the labelled pair
  members <- strsplit(pairwise_tran_results$comparison, " vs ")
  expect_true(all(mapply(function(m, cl, pt) setequal(m, c(cl, pt)),
                         members,
                         pairwise_tran_results$cluster,
                         pairwise_tran_results$partner)))

  # Test 3: group carries the group.by level, not the list position
  expect_setequal(unique(pairwise_tran_results$group), unique(scRep_example$Patient))
  expect_setequal(unique(pairwise_migr_results$group), unique(scRep_example$Patient))

  standard_results <- StartracDiversity(scRep_example,
                                        type = "Type",
                                        group.by = "Patient",
                                        export.table = TRUE)
  expect_setequal(unique(standard_results$group), unique(scRep_example$Patient))
})


test_that("Pairwise plots are built for both pairing variables", {
  # Anchor/partner faceting, only reachable when pairing over clusters
  tran_plot <- StartracDiversity(scRep_example,
                                 type = "Type",
                                 group.by = "Patient",
                                 index = "tran",
                                 pairwise = "cluster")
  expect_s3_class(tran_plot, "ggplot")
  built <- ggplot2::ggplot_build(tran_plot)
  expect_equal(tran_plot$labels$x, "Partner Cluster")
  # One facet per anchor cluster, and every cluster appears as an anchor
  expect_true("cluster" %in% names(built$layout$layout))
  expect_setequal(as.character(built$layout$layout$cluster),
                  as.character(unique(tran_plot$data$cluster)))

  # Pairing over a non-cluster variable keeps the comparison on the x-axis
  migr_plot <- StartracDiversity(scRep_example,
                                 type = "Type",
                                 group.by = "Patient",
                                 index = "migr",
                                 pairwise = "Type")
  expect_s3_class(migr_plot, "ggplot")
  expect_equal(migr_plot$labels$x, "Comparison")
  # One box per comparison, summarizing across the group.by levels
  expect_equal(nrow(ggplot2::ggplot_build(migr_plot)$data[[1]]),
               length(unique(migr_plot$data$comparison)))
})


test_that(".calculateIndices() handles degenerate input", {
  # No cells at all
  expect_null(.calculateIndices(toy_processed[0, ], c("expa", "migr", "tran")))

  # Cells present but no cluster assignment, so there is nothing to index
  no_cluster <- toy_processed
  no_cluster$cluster <- NA
  expect_null(.calculateIndices(no_cluster, c("expa", "migr", "tran")))

  # Clusters present but no clonotypes: every index is undefined, not zero
  no_clones <- toy_processed
  no_clones$clone.id <- NA
  empty_result <- .calculateIndices(no_clones, c("expa", "migr", "tran"))
  expect_s3_class(empty_result, "data.frame")
  expect_equal(empty_result$cluster, c("1", "2"))
  expect_true(all(is.na(empty_result$expa)))
  expect_true(all(is.na(empty_result$migr)))
  expect_true(all(is.na(empty_result$tran)))

  # Two clones split evenly over two clusters and two locations: maximal
  # transition and migration, and no expansion within either cluster
  full_result <- .calculateIndices(toy_processed, c("expa", "migr", "tran"))
  expect_equal(full_result$tran, c(1, 1))
  expect_equal(full_result$migr, c(1, 1))
  expect_equal(full_result$expa, c(0, 0))
})


test_that(".calculatePairwiseIndices() handles degenerate input", {
  # Fewer than two cells, or only one category, leaves nothing to pair
  expect_null(.calculatePairwiseIndices(toy_processed[1, ], "tran", "cluster"))

  one_cluster <- toy_processed
  one_cluster$cluster <- "1"
  expect_null(.calculatePairwiseIndices(one_cluster, "tran", "cluster"))

  # Two clusters but no clonotypes to distribute across them
  no_clones <- toy_processed
  no_clones$clone.id <- NA
  expect_equal(nrow(.calculatePairwiseIndices(no_clones, "tran", "cluster")), 0)

  # Both members of the pair are anchored, and a clone split evenly across
  # the pair scores 1 from either side
  pair_result <- .calculatePairwiseIndices(toy_processed, "tran", "cluster")
  expect_equal(nrow(pair_result), 2)
  expect_equal(pair_result$comparison, c("1 vs 2", "1 vs 2"))
  expect_equal(pair_result$cluster, c("1", "2"))
  expect_equal(pair_result$partner, c("2", "1"))
  expect_equal(pair_result$value, c(1, 1))

  # Pairing over location leaves the anchor on clusters, so no partner column
  loc_result <- .calculatePairwiseIndices(toy_processed, "migr", "loc")
  expect_false("partner" %in% names(loc_result))
  expect_equal(unique(loc_result$comparison), "B vs L")

  # Order of appearance must not change the label
  shuffled <- toy_processed[c(2, 1, 4, 3), ]
  expect_equal(.calculatePairwiseIndices(shuffled, "tran", "cluster")$comparison,
               pair_result$comparison)
})

