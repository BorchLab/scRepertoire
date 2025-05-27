getdata <- function(dir, name) {
	readRDS(paste("testdata/", dir, "/", name, ".rds", sep = ""))
}

# shortcut get the combined contigs that are used in most testcases
getCombined <- function() {
	getdata("combineContigs", "combined")
}
