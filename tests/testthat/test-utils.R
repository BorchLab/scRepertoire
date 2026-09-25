# test script for utils.R - testcases are NOT comprehensive!

# Mock data frame simulating contig information
mock_df_contigs <- data.frame(
  barcode = rep(c("bc1", "bc2", "bc3", "bc4"), each = 2),
  chain = rep(c("TRA", "TRB"), 4),
  reads = c(100, 90, 200, 180, 50, 40, 300, 290),
  CTgene = c("TRAV1_TRBV1", "TRAV1_TRBV1", "TRAV2_TRBV2", "TRAV2_TRBV2", 
             "TRAV3_TRBV3", "TRAV3_TRBV3", "TRAV4_TRBV4;TRAV5_TRBV5", "TRAV4_TRBV4"),
  CTnt = c("AAA_CCC", "AAA_CCC", "GGG_TTT", "GGG_TTT", 
           "CCC_AAA", "CCC_AAA", "GGG_TTT;AAA_CCC", "GGG_TTT"),
  CTaa = c("K_C", "K_C", "G_F", "G_F", "P_K", "P_K", "G_F;K_C", "G_F"),
  CTstrict = c("TRAV1;AAA_TRBV1;CCC", "TRAV1;AAA_TRBV1;CCC", "TRAV2;GGG_TRBV2;TTT", "TRAV2;GGG_TRBV2;TTT",
               "TRAV3;CCC_TRBV3;AAA", "TRAV3;CCC_TRBV3;AAA", "TRAV4;GGG;TRAV5_TRBV4;TTT", "TRAV4;GGG_TRBV4;TTT"),
  sample = rep(c("s1", "s2"), each = 4),
  cluster = rep(c("c1", "c1", "c2", "c2"), 2)
)

# Mock metadata for clone counting
mock_meta <- data.frame(
  barcode = paste0("cell", 1:10),
  group = rep(c("group1", "group2"), each = 5),
  CTnt = c("AAA", "AAA", "GGG", "CCC", NA, "AAA", "GGG", "GGG", "TTT", "TTT"),
  clone.size = c(3, 3, 2, 1, NA, 3, 2, 2, 2, 2)
)

# Mock list of data frames
mock_df_list <- list(
  sample1 = data.frame(
    barcode = c("s1_bc1", "s1_bc2"),
    CTnt = c("AAA", "GGG"),
    CTaa = c("K", "G")
  ),
  sample2 = data.frame(
    barcode = c("s2_bc1"),
    CTnt = c("TTT"),
    CTaa = c("F")
  ),
  sample3_all_na = data.frame(
    barcode = c("s3_bc1", "s3_bc2"),
    CTnt = c(NA, NA),
    CTaa = c(NA, NA)
  ),
  sample4_empty = data.frame()
)

# Mock Seurat/SCE objects - simplified S4 objects for testing
setClass("mock_seurat", representation(meta.data = "data.frame", active.ident = "factor"))
setClass("mock_sce", representation(colData = "DataFrame"))



test_that(".toCapitilize works correctly", {
  expect_equal(.toCapitilize("hello"), "Hello")
  expect_equal(.toCapitilize("WORLD"), "World")
  expect_equal(.toCapitilize("miXeD"), "Mixed")
  expect_equal(.toCapitilize(""), "")
  expect_equal(.toCapitilize(c("one", "two")), c("One", "Two"))
})

test_that(".orderingFunction works correctly", {
  df <- data.frame(group = c("g2", "g10", "g1"))
  
  # Test alphanumeric ordering
  df_alpha <- .orderingFunction(df, group.by = "group", vector = "alphanumeric")
  expect_equal(levels(df_alpha$group), c("g1", "g2", "g10"))
  
  # Test custom vector ordering
  custom_order <- c("g10", "g1", "g2")
  df_custom <- .orderingFunction(df, group.by = "group", vector = custom_order)
  expect_equal(levels(df_custom$group), custom_order)
})

test_that(".chainPositionParser works correctly", {
  expect_equal(.chainPositionParser("TRA"), 1)
  expect_equal(.chainPositionParser("trb"), 2)
  expect_equal(.chainPositionParser("IGH"), 1)
  expect_equal(.chainPositionParser("IGL"), 2)
  
  # Test error on invalid input
  expect_error(.chainPositionParser("INVALID"), 
               "'INVALID' is not a valid entry for the `chain` argument.")
})

test_that(".offTheChain works correctly", {
  df <- data.frame(
    barcode = c("bc1", "bc2", "bc3"),
    CTgene = c("TRA.V1_TRB.V1", "IGH.V2_IGL.V2", "TRA.V3_NA"),
    CTaa = c("K_C", "G_F", "P_NA")
  )
  
  # Test TRA chain extraction
  df_tra <- .offTheChain(df, chain = "TRA", cloneCall = "CTaa", check = FALSE)
  expect_equal(df_tra$CTaa, c("K", "G", "P"))
  
  # Test TRB chain extraction
  df_trb <- .offTheChain(df, chain = "TRB", cloneCall = "CTaa", check = FALSE)
  expect_equal(df_trb$CTaa, c("C", "F", NA))
})

test_that(".cloneCounter works correctly", {
  counts <- .cloneCounter(mock_meta, group.by = "group", cloneCall = "CTnt")
  
  # Check dimensions
  expect_equal(nrow(counts), 8)
  expect_equal(ncol(counts), 5)
  
  # Check calculations for a specific clone
  aaa_group1 <- subset(counts, group == "group1" & CTnt == "AAA")
  expect_equal(aaa_group1$n, 2)
  expect_equal(aaa_group1$group.sum, 4)
  expect_equal(aaa_group1$clone.sum, 3)
})

test_that(".colorizer works correctly", {
  # Test default palette
  colors <- .colorizer(n = 5)
  expect_type(colors, "character")
  expect_length(colors, 5)
  
  # Test different palette
  colors_viridis <- .colorizer(palette = "viridis", n = 10)
  expect_type(colors_viridis, "character")
  expect_length(colors_viridis, 10)
  expect_false(identical(colors, colors_viridis[1:5]))
})

test_that(".checkBlanks works correctly", {
  filtered_list <- .checkBlanks(mock_df_list, cloneCall = "CTnt")
  
  # Expecting sample3 (all NA) and sample4 (empty) to be removed
  expect_length(filtered_list, 2)
  expect_equal(names(filtered_list), c("sample1", "sample2"))
})

test_that(".groupList works correctly", {
  df_bound <- bind_rows(mock_df_list[1:2], .id = "sample")
  df_grouped <- .groupList(df_bound, group.by = "sample")
  
  expect_type(df_grouped, "list")
  expect_length(df_grouped, 2)
  expect_equal(names(df_grouped), c("sample1", "sample2"))
  expect_equal(nrow(df_grouped$sample1), 2)
})

test_that(".checkList works correctly", {
  # Test with a data frame
  df <- data.frame(a = 1)
  expect_type(.checkList(df), "list")
  expect_length(.checkList(df), 1)
  
  # Test with a list
  ls <- list(df)
  expect_identical(.checkList(ls), ls)
})

test_that(".checkContigs works correctly", {
  df_list <- list(data.frame(a = c(1, ""), b = c("", 2)))
  checked_list <- .checkContigs(df_list)
  
  expect_true(is.na(checked_list[[1]][2, "a"]))
  expect_true(is.na(checked_list[[1]][1, "b"]))
})


test_that(".modifyBarcodes works correctly", {
  df <- list(data.frame(barcode = "bc1"), data.frame(barcode = "bc2"))
  
  # Test with samples only
  modified <- .modifyBarcodes(df, samples = c("s1", "s2"), ID = NULL)
  expect_equal(modified[[1]]$barcode, "s1_bc1")
  expect_equal(modified[[2]]$barcode, "s2_bc2")
  
  # Test with samples and ID
  modified_id <- .modifyBarcodes(df, samples = c("s1", "s2"), ID = c("id1", "id2"))
  expect_equal(modified_id[[1]]$barcode, "s1_id1_bc1")
  expect_equal(modified_id[[2]]$barcode, "s2_id2_bc2")
})

test_that(".removingNA works correctly", {
  df <- list(data.frame(a = c(1, NA), b = c(3, 4)))
  cleaned <- .removingNA(df)
  expect_equal(nrow(cleaned[[1]]), 1)
  expect_equal(cleaned[[1]]$a, 1)
})

test_that(".removingMulti works correctly", {
  df <- list(data.frame(CTnt = c("A;B", "C", "D;E")))
  cleaned <- .removingMulti(df) 
  expect_equal(nrow(cleaned[[1]]), 1)
  expect_equal(cleaned[[1]]$CTnt, "C")
})

test_that(".filteringMulti works correctly", {
  df <- data.frame(
    barcode = c("bc1", "bc1", "bc1", "bc2", "bc2"),
    chain = c("TRA", "TRB", "TRG", "TRA", "TRB"),
    reads = c(100, 90, 80, 200, 180)
  )
  filtered <- .filteringMulti(df)
  
  expect_equal(nrow(filtered), 4)
  expect_false("TRG" %in% subset(filtered, barcode == "bc1")$chain)
  expect_equal(subset(filtered, barcode == "bc1" & chain == "TRA")$reads, 100)
})

test_that(".theCall and .convertClonecall work together", {
  

  # Test dictionary conversion
  expect_equal(.theCall(mock_df_contigs, "gene"), "CTgene")
  expect_equal(.theCall(mock_df_contigs, "nt"), "CTnt")
  expect_equal(.theCall(mock_df_contigs, "amino"), "CTaa")
  expect_equal(.theCall(mock_df_contigs, "strict"), "CTstrict")
  
  # Test custom variable with warning
  expect_message(.theCall(mock_df_contigs, "cluster"), "A custom variable cluster will be used to call clones")
})

test_that(".assignCT works correctly", {
  con_df_t <- data.frame(TCR1 = "TRAV1", TCR2 = "TRBV1", cdr3_nt1 = "AAA", cdr3_nt2 = "CCC", cdr3_aa1 = "K", cdr3_aa2 = "C")
  con_df_b <- data.frame(IGH = "IGHV1", IGLC = "IGLC1", cdr3_nt1 = "GGG", cdr3_nt2 = "TTT", cdr3_aa1 = "G", cdr3_aa2 = "F")
  
  # Test for T cells
  t_assigned <- .assignCT("T", con_df_t)
  expect_equal(t_assigned$CTgene, "TRAV1_TRBV1")
  expect_equal(t_assigned$CTnt, "AAA_CCC")
  expect_equal(t_assigned$CTaa, "K_C")
  expect_equal(t_assigned$CTstrict, "TRAV1;AAA_TRBV1;CCC")
  
  # Test for B cells
  b_assigned <- .assignCT("B", con_df_b)
  expect_equal(b_assigned$CTgene, "IGHV1_IGLC1")
  expect_equal(b_assigned$CTnt, "GGG_TTT")
  expect_equal(b_assigned$CTaa, "G_F")
})

test_that(".makeGenes works correctly", {
  data_t <- data.frame(chain = c("TRA", "TRB"), v_gene = c("vA", "vB"), j_gene = c("jA", "jB"), c_gene = c("cA", "cB"), d_gene = c(NA, "dB"))
  data_b <- data.frame(chain = c("IGH", "IGK"), v_gene = c("vH", "vK"), j_gene = c("jH", "jK"), c_gene = c("cH", "cK"), d_gene = c("dH", NA))
  
  # Test T-cell gene creation
  t_genes <- .makeGenes("T", data_t)
  expect_equal(t_genes$TCR1[1], "vA.jA.cA")
  expect_equal(t_genes$TCR2[2], "vB.dB.jB.cB")
  
  # Test B-cell gene creation
  b_genes <- .makeGenes("B", data_b)
  expect_equal(b_genes$IGHct[1], "vH.dH.jH.cH")
  expect_equal(b_genes$IGKct[2], "vK.jK.cK")
})

test_that(".lengthDF works correctly", {
  df_list <- list(sample1 = data.frame(CTnt = c("AAA_CCC", "GGG_TTT"), group = c("A", "B")))
  
  # Test with chain = "both"
  ldf_both <- .lengthDF(df_list, cloneCall = "CTnt", chain = "both", group = "group")
  expect_equal(names(ldf_both), c("length", "CT", "group", "values"))
  expect_equal(ldf_both$length, c(6, 6))
  
  # Test with specific chain
  ldf_chain <- .lengthDF(df_list, cloneCall = "CTnt", chain = "TRA", group = "group")
  expect_equal(names(ldf_chain), c("length", "CT", "values", "chain", "group"))
  expect_equal(ldf_chain$length, c(7, 7)) # Length of first part of the string
})




# --- .deprecate_arg ----------------------------------------------------------

test_that(".deprecate_arg returns the right value in each branch", {
  # Only the new argument
  expect_equal(.deprecate_arg(NULL, "new", "old_nm", "new_nm", "fn"), "new")
  # Neither argument falls back to the default
  expect_equal(.deprecate_arg(NULL, NULL, "old_nm", "new_nm", "fn",
                              default = "dflt"), "dflt")
  expect_null(.deprecate_arg(NULL, NULL, "old_nm", "new_nm", "fn"))
  # Only the old argument is honoured. Signalling depends on the call depth
  # lifecycle sees, so it is covered separately below through an exported
  # function rather than asserted on this direct call.
  expect_equal(
    suppressWarnings(.deprecate_arg("old", NULL, "old_nm", "new_nm", "fn")),
    "old"
  )
  # Both provided prefers the new one
  expect_equal(
    suppressWarnings(.deprecate_arg("old", "new", "old_nm", "new_nm", "fn")),
    "new"
  )
})

# A minimal combined object for the deprecation-routing test below
mock_combined_deprecate <- combineTCR(contig_list[1:2], samples = c("A", "B"))

test_that("deprecation warnings reach the calling user, not just the test harness", {
  # lifecycle decides whether to signal from `user_env`. Its default lands two
  # frames above the signaller, which is .deprecate_arg's caller -- a frame
  # inside this namespace -- so every deprecation was classified as indirect and
  # deprecate_soft() stayed silent for real users. TESTTHAT_PKG hid that during
  # testing, because from_testthat() also counts as direct. Unset it and call
  # from the global environment, which is what a user's script looks like.
  # "warning" defeats the once-per-session throttle without bypassing the
  # directness check, which is what this test is about.
  withr::local_envvar(TESTTHAT_PKG = NA)
  withr::local_options(lifecycle_verbosity = "warning")

  caller <- function(x) clonalBin(x, cloneCall = "gene")
  environment(caller) <- globalenv()

  expect_warning(caller(mock_combined_deprecate), regexp = "cloneCall.*deprecated")
})
