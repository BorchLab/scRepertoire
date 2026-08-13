# test script for StartracDiversity.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list, samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"))
scRep_example <- combineExpression(combined, scRep_example)
scRep_example$Patient <- substring(scRep_example$orig.ident, 1, 3)
scRep_example$Type <- substring(scRep_example$orig.ident, 4, 4)

test_that("Input validation and error handling work correctly", {
  expect_error(
    StartracDiversity(scRep_example, index = c("migr", "tran"), pairwise = "Type"),
    "Pairwise analysis can only be performed for a single index"
  )
  
  expect_error(
    StartracDiversity(scRep_example, index = "expa", pairwise = "Type"),
    "Pairwise analysis is only supported for 'migr' or 'tran' indices."
  )
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


