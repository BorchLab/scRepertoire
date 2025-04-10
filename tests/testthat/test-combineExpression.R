# test script for combineExpression.R - testcases are NOT comprehensive!

test_that("combineExpression works with seurat objects", {
	data("scRep_example")
	combined <- getCombined()
  	combined_test <- combineExpression(combined, scRep_example)
	
	#Seurat object test
	expect_length(combined_test@meta.data, 13)
	expect_equal(combined_test@meta.data[, 1:6], scRep_example@meta.data[, 1:6])
	expect_equal(
		combined_test@meta.data[, 7:13],
		getdata("seuratFunctions", "combineExpression_new_metadata")
	)
	
	#Single-cell experiment test
	#as.SingleCellExperiment will add an ident column 13 --> 14
	sce_test <- Seurat::as.SingleCellExperiment(scRep_example)
	sce_test <- combineExpression(combined, sce_test)
	expect_length(sce_test@colData, 14)
	expect_equal(as.data.frame(sce_test@colData[, 1:6]), scRep_example@meta.data[, 1:6])
	
	
	#Testing group.by parameter
	combined <- addVariable(combined, 
	                        variable.name = "Patient", 
	                        variables = c("P17", "P17", "P18", "P18", "P19", "P19", "P20", "P20"))
	combined_test.2 <- combineExpression(combined, scRep_example, group.by = "Patient")
	
	expect_equal(
	  combined_test.2@meta.data[, 7:13],
	  getdata("seuratFunctions", "combineExpression_grouping_metadata")
	)
})
