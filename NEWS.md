# scRepertoire VERSION 2.9.2

## BUG FIXES

* `getCirclize()` self-links (the diagonal, `include.self = TRUE`) no longer always return 0. The private-clone subtraction counted zero-count clone-by-group rows from the underlying cross-tabulation as shared, cancelling the full count. It now restricts to clones actually observed (`n > 0`) in other groups, so the diagonal correctly reports each group's private clones. The same `n > 0` guard was applied to the `n.shared` statistic in the `include.metadata = TRUE` output.
* `clonalHomeostasis()` now orders the clonal groups by their proportion bins (`Hyperexpanded -> Large -> Medium -> Small -> Rare`) instead of alphabetically. The fill categories were being sorted as a character vector, which placed `Rare` before `Small`.

# scRepertoire VERSION 2.9.1

## NEW FEATURES

* Full-length sequence retention (resolves #535). `loadContigs()` now keeps a standardized `sequence`/`sequence_aa` column for every format (reconstructed from framework/CDR regions for 10x and MiXCR, copied natively for AIRR-family inputs, `NA` where unavailable), plus `sequence_alignment`/`germline_alignment` when present. This is additive and does not change existing columns.
* `combineTCR()` and `combineBCR()` gain `retain.sequences` (default `FALSE`). When set, full-length sequences are carried into the combined object as per-chain columns (`sequence_nt1`/`sequence_nt2`, `sequence_aa1`/`sequence_aa2`, and optionally `sequence_alignment*`/`germline*`) without altering `CTaa`/`CTnt`/`CTgene`/`CTstrict`. The retained columns flow through `combineExpression()` into single-cell metadata.
* `exportClones(format = "airr")` now populates `sequence`/`sequence_aa` (and alignment/germline) fields when those sequences were retained.
* New `exportDowser()` function (and `exportClones(format = "dowser")`) reshapes BCR/TCR data into an AIRR data frame ready for `dowser::formatClones()`, enabling B cell lineage / phylogenetic analysis in the Immcantation framework (resolves #578). `exportDowser()` errors with guidance when the required IMGT-gapped alignments are absent.

## UNDERLYING CHANGES

* `clonalCluster()` now passes the tight similarity threshold to `immApex::buildNetwork()` instead of the loose value that was refiltered afterward. The engine prunes early rather than returning a large edge list that is mostly discarded. Raw thresholds are unchanged, and cluster assignments match the previous behavior.
* `clonalCluster()` selects the edge expansion from `cluster.method`. Connected-components clustering uses the memory-efficient star expansion, which is exact for components. Community-detection methods use the full clique expansion for exact edge multiplicity.

## BUG FIXES

* `combineBCR(retain.sequences = ...)` now resolves the retained columns against the intersection of columns shared by all samples, rather than the first sample only. This prevents an `"undefined columns selected"` error (or silent dropping of a column) when samples carry heterogeneous fields, matching `combineTCR()`.
* Retained per-chain sequences from multi-contig chains no longer carry literal `"NA"` placeholder tokens (e.g. `"seq1;NA"`) when a contig lacks a full-length sequence; such tokens are dropped from the join and an all-placeholder value becomes a real `NA`.

## DEPENDENCIES

* Added `dowser` to `Suggests` for the new lineage export path.
* `clonalCluster()` requires `immApex (>= 1.7.1)` for the new `expand` argument and the faster network engine.

## DOCUMENTATION

* Documented full-length sequence retention and the `dowser` export path in the main vignette and the "Combining Contigs into Clones" article.
* Added a "Exporting scRepertoire Clones to dowser" section to the Immcantation article covering the `retain.sequences` to `exportDowser()` to `dowser::formatClones()` workflow.

# scRepertoire VERSION 2.7.3

## BUG FIXES

* Fixed `combineBCR()` assigning the same `group.by` value to all cells instead of per-barcode values.
* Fixed `clonalCluster()` failing when `group.by` produces a single group by correcting the condition for adding the group column to the edge list.
* Fixed test for CTstrict clustering pattern to match the new default `chain = "IGH"` behavior in `combineBCR()`.

# scRepertoire VERSION 2.7.2

## NEW FEATURES

* New `clonalBin()` function to bin clones by frequency or proportion without requiring a single-cell object. Adds `clonalFrequency`, `clonalProportion`, and `cloneSize` columns to the output of `combineTCR()`, `combineBCR()`, or `combineExpression()`. Supports custom bin thresholds, optional grouping by metadata variable, and chain filtering.
* New `vizCirclize()` function for quick chord diagram visualization without manual circlize code. Supports directional arrows, custom colors, and sector annotations.
* `getCirclize()` major enhancements:
  * Multi-level hierarchical grouping via `group.by` accepting a vector of columns
  * New `method` parameter: `"unique"`, `"abundance"`, `"jaccard"`, `"overlap"`
  * New `symmetric` parameter for directional flow analysis
  * New `include.metadata` parameter returning sector statistics
  * New filtering options: `min.shared`, `top.links`, `filter.sectors`
  * Built-in color palette generation with `palette` parameter
* `alluvialClones()` major enhancements:
  * New `top.clones`, `min.freq`, `highlight.clones`, and `highlight.color` parameters
  * Visual customization: `stratum.width`, `flow.alpha`, `show.labels`, `label.size`
  * New `order.strata` parameter for controlling level ordering within each stratum
  * Enhanced `export.table` output now includes `freq`, `prop`, and `rank` columns
* `combineBCR()` defaults the clustering call to "IGH" instead of "both"

## API CHANGES

* Soft-deprecated camelCase arguments across all exported functions in favor of dot.notation (e.g., `cloneCall` to `clone.call`, `exportTable` to `export.table`, `cloneSize` to `clone.size`, `filterNA` to `filter.na`, `addLabel` to `add.label`, `clonalSplit` to `clonal.split`). All deprecated arguments will continue to work with a deprecation warning until version 3.0.0.

## BUG FIXES

* Fixed `combineExpression()` failing with "undefined columns selected" when input data already contained `clonalFrequency`/`clonalProportion` columns from a prior `clonalBin()` call.

# scRepertoire VERSION 2.6.2

## UNDERLYING CHANGES

* Expanded functionality in `combineBCR()` and `clonalCluster()`:
  * New metrics beyond normalized Levenshtein edit distances
  * Support for raw and normalized-based calculations
  * Support for distance matrices to allow for alignment
* Added support for declaring `chains = "IGL"`, `"IGK"`, or `"Light"` to get all light chains in downstream quantification.

## BUG FIXES

* Fixed handling of multiple chains in `combineBCR()`, specifically in formatting CTstrict.

# scRepertoire VERSION 2.6.1

Update to match Bioconductor Release 3.22 on 2025/10/30.

## BUG FIXES

* Fixed `order.by` issue in `positionalProperty()`.
* Fixed individual chain call for `combineExpression()`.
* Fixed issue with removing kmer with ";" in `percentKmer()`.

# scRepertoire VERSION 2.5.7

## BUG FIXES

* Fixed chain handling for BCR genes in `percentGeneUsage()` and propagated to wrappers: `percentGenes()`, `percentVJ()`, and `vizGenes()`.

# scRepertoire VERSION 2.5.6

## BUG FIXES

* Fixed passing `group.by` for `combineBCR()`.

# scRepertoire VERSION 2.5.5

## UNDERLYING CHANGES

* Updated unit tests for ggplot2 v4.

# scRepertoire VERSION 2.5.3

## BUG FIXES

* Fixed `clonalProportion` calculation to use grouping properly during `combineExpression()`.
* Fixed immunarch support for `exportClones()` TRA/Light chain column handling.

# scRepertoire VERSION 2.5.2

## UNDERLYING CHANGES

* Added support for mouse genes in `quietBCRgenes()` and `quietTCRgenes()`.

# scRepertoire VERSION 2.5.1

## UNDERLYING CHANGES

* Introduced pairwise calculations to `StartracDiversity()`.
* Internal function conversion for `clonalSizeDistribution()` - removed cubature, truncdist, and VGAM from dependencies.
* Increased speed of `clonalSizeDistribution()`.

# scRepertoire VERSION 2.5.0

## UNDERLYING CHANGES

* Updated and improved code for `loadContigs()`.
* Consolidated support for discrete AIRR formats under the umbrella of AIRR.
* Added `"tcrpheno"` and `"immunarch"` to `exportClones()`.
* Converted `exportClones()` to base R to reduce dependencies.
* Added dandelionR and tcrpheno vignettes to pkgdown site.
* `percentAA()` refactored to minimize dependencies and use immApex `calculateFrequency()`.
* `positionalEntropy()` refactored to minimize dependencies and use immApex `calculateEntropy()`.
* `clonalDiversity()` refactored for performance - now calculates a single diversity metric at a time and includes new estimators like `"gini"`, `"d50"`, and supports hill numbers.
* `percentKmer()` refactored to use immApex `calculateMotif()` for both aa and nt sequences. No longer calculates all possible motifs, but only motifs present.
* `clonalCluster()` now allows for dual-chain clustering, V/J filtering, normalized or straight edit distance calculations, and return of clusters, igraph objects, or adjacency matrix.
* `combineBCR()` offers single/dual chain clustering, aa or nt sequences, adaptive filtering of V and J genes, and normalized or straight edit distance calculations.
* `percentGeneUsage()` is now the underlying function for `percentGenes()`, `percentVJ()`, and `vizGenes()` and allows for percent, proportion, and raw count quantification.
* Added common theme (internal `.themeRepertoire()`) to all plots and allow users to pass arguments to it.

## BUG FIXES

* `clonalCompare()` issue with plotting a 0 row data frame now errors with message.
* `clonalScatter()` `group.by`/axes call now works for non-single-cell objects.
* Fixed issue with NULL and "none" `group.by` in `combineExpression()`.
* Allowing multi groupings via `x.axis` and `group.by` in `clonalDiversity()`.

# scRepertoire VERSION 2.3.4

## UNDERLYING CHANGES

* Updated internal `.parseContigs()` to function with more complex groupings.
* Added `annotateInvariant()` functionality for mouse and human TCRs.
* Added `quietTCRgenes()`, `quietBCRgenes()`, `quietVDJgenes()`.
* Fixed issue with `clonalCompare()` assertthat statements.
* Started integration with immApex API package.

# scRepertoire VERSION 2.3.2

## BUG FIXES

* Fixed issue with denominator in `getCirclize()`.
* Fixed chain issue with `clonalCompare()` - expanded assertthat statement.

# scRepertoire VERSION 2.2.1

Update for Bioconductor version 3.20.

## NEW FEATURES

* Added support for BCRs for loading ParseBio sequences.
* Added `quietBCRgenes()`, `quietTCRgenes()`, and `quietVDJgenes()` for filtering out known TCR and/or BCR gene signatures.

## UNDERLYING CHANGES

* Added `Seurat` to the `Suggests` field in the DESCRIPTION file.

# scRepertoire VERSION 2.0.8

## NEW FEATURES

* Added `getContigDoublets()` experimental function to identify TCR and BCR doublets as a preprocessing step to `combineExpression()`.
* Added `proportion` argument to `clonalCompare()` so that when set to FALSE, the comparison will be based on frequency normalized by per-sample repertoire diversity.

## UNDERLYING CHANGES

* Fixed issue with single chain output for `clonalLength()`.
* Removed unnecessary code remnant in `clonalLength()`.
* Allow one sample to be plotted by `percentVJ()`.
* Fixed issue with `positionalProperty()` and `exportTable`.
* Fixed issue with `loadContigs()` edge case when TRUST4 data only has 1 row.
* Converted documentation to use markdown (roxygen2md).
* Imported `lifecycle`, `purrr`, `withr`.
* Fixed issue with `clonalCluster()` and `exportGraph = TRUE`.
* Improved performance of `combineBCR()` by a constant factor with C++.
* Restructured functions to `exportTable` before plotting.

# scRepertoire VERSION 2.0.7

## BUG FIXES

* Fixed issue with `group.by` in `clonalOverlap()`.
* Fixed issue with `group.by` in `clonalCompare()`.

# scRepertoire VERSION 2.0.6

## BUG FIXES

* Fixed issue with custom column headers for clones.

# scRepertoire VERSION 2.0.5

## UNDERLYING CHANGES

* Added type checks using assertthat.
* Updated conditional statements in `constructConDFAndparseTCR.cpp`.
* Fixed issue in `clonalQuant()` and factor-based `group.by` variable.

# scRepertoire VERSION 2.0.4

## UNDERLYING CHANGES

* `getCirclize()` refactored to prevent assumptions and added `include.self` argument.
* Added `.count.clones()` internal function for `getCirclize()` and `clonalNetwork()`.
* Added `order.by` parameter to visualizations to specifically call order of plotting using a vector or `"alphanumeric"`.
* Fixed issue with `clonalLength()` and NA handling.
* `clonalCompare()` now retains the original clonal info if using `relabel.clones`.
* Added Dandelion support to `loadContigs()` and testthat.
* Fixed issue with `positionalProperty()` assumption that clones will all have 20 amino acids.
* Fixed issue with `positionalProperty()` and removing non-amino acids.
* Fixed IGH/K/L mistaking gene issue in `vizGenes()`.
* Added error message for NULL results in `clonalCluster()` with `export.graph = TRUE`.
* Fixed issue with "full.clones" missing in `combineExpression()` when using 1 chain.

# scRepertoire VERSION 2.0.3

## UNDERLYING CHANGES

* Modified support for Omniscope format to allow for dual chains.
* Added ParseBio support to `loadContigs()` and testthat.
* Added support for productive variable to `loadContigs()` for BD, Omniscope, and Immcantation formats.
* Replaced numerical indexing with name indexing for `loadContigs()`.
* `combineBCR()` and `combineTCR()` now allow for unproductive contig inclusions with new `filterNonproductive` parameter.
* `combineBCR()` will now prompt user if `samples` is not included instead of erroring.
* Added base threshold by length for internal `.lvCompare()`.
* Ensured internal `.lvCompare()` only looks at first set of sequences in multi-sequence chain.
* Fixed bug in exporting graph for `clonalCluster()`.
* Fixed conflict in functions between igraph and dplyr packages.

# scRepertoire VERSION 2.0.2

## BUG FIXES

* `clonalOccupy()` rewrite counting and NA handling.

# scRepertoire VERSION 2.0.1

## UNDERLYING CHANGES

* `clonalOverlay()` arguments now `cutpoint` and use `cut.category` to select either `clonalProportion` or `clonalFrequency`.

# scRepertoire VERSION 2.0.0 (2024-01-10)

## NEW FEATURES

* Added `percentAA()`, `percentGenes()`, `percentVJ()`, `percentKmer()`, `exportClones()`, `positionalEntropy()`, `positionalProperty()`.
* Renamed functions: `compareClonotypes` to `clonalCompare()`, `clonotypeSizeDistribution` to `clonalSizeDistribution()`, `scatterClonotypes` to `clonalScatter()`, `quantContig` to `clonalQuant()`, `highlightClonotypes` to `highlightClones()`, `lengthContigs` to `clonalLength()`, `occupiedscRepertoire` to `clonalOccupy()`, `abundanceContig` to `clonalAbundance()`, `alluvialClonotypes` to `alluvialClones()`.
* Added features to `clonalCompare()` to allow for highlighting sequences and relabeling clonotypes.

## UNDERLYING CHANGES

* Removed internal `.quiet()` function.
* `.theCall()` now allows for a custom header/variable and checks the colnames.
* Replaced data arguments to be more descriptive: `df` to `input.data`, `dir` to `input`, `sc` to `sc.data`.
* Deep clean on documentation for increased consistency and explainability.
* `StartracDiversity()` metric re-implemented to remove startrac-class object intermediary.
* Implemented powerTCR locally to reduce dependencies and continue support.
* Universalized underlying function language and intermediate variables.
* License change to MIT.
* `group.by` and `split.by` consolidated into single `group.by` parameter.
* Added support for Immcantation pipeline, .json, Omniscope, and MiXCR formats for `loadContigs()`.
* Made GitHub.io website for support/vignettes/FAQ.
* Added testthat for all exported and internal functions.
* Fixed issue with `clonalQuant()` for instance of `scale = FALSE` and `group.by` being set.
* `clonalDiversity()` no longer automatically orders samples.
* Removed `order` parameter from `clonalQuant()`, `clonalLength()`, and `clonalAbundance()`.
* `x.axis` parameter in `clonalDiversity()` separated from `group.by` parameter.
* Filtering chains will not eliminate non-matching chains.

## DEPRECATED AND DEFUNCT

* Deprecated `stripBarcodes()`.
* Deprecated `expression2List()` (now only an internal function).
* Deprecated `checkContigs()`.

# scRepertoire VERSION 1.11.0

* Rebasing for Bioconductor version.

# scRepertoire VERSION 1.7.5

* Fixed `combineBCR()` to allow for non-related sequences.

# scRepertoire VERSION 1.7.4

## NEW FEATURES

* `checkContigs()` function to quantify the percentages of NA values by genes or sequences.
* `exportClones` to `clonalNetwork()` to isolate clones shared across identities.

## UNDERLYING CHANGES

* Fixed issue with `clonalDiversity()` and skipping boots.
* Fixed underlying assumptions with `clonalBias()`.
* Added reads variable to `parseAIRR`.
* Fixed handling of samples parameter in combine functions.
* Removed need to relevel the cloneType factor after `combineExpression()`.
* Set up `lapply()` for `combineBCR()` and `clusterTCR()` - no more pairwise distance matrix calculation.
* `loadContigs()` support for data.frames or lists of contigs.
* Added examples for `loadContigs()`.
* Removed requirement for T cell type designation - will combine A/B and G/D simultaneously.
* Updated `combineBCR()` to chunk nucleotide edit distance calculations by V gene and give option to skip edit distance calculation with `call.related.clones = FALSE`.
* Updated `clusterTCR()` to use `lvCompare()` and base edit distances on V gene usage.

# scRepertoire VERSION 1.7.3

## UNDERLYING CHANGES

* Fixed misspellings for parse functions.
* Optimized WAT3R and 10x `loadContigs()`.
* Removed `combineTRUST4` - superseded by `loadContigs()`.
* Added support of TRUST4 for `combineBCR()`.
* Added support for BD in `loadContigs()`.
* `loadContigs()` TRUST4 parsing allows for all NA values in a chain.
* `combineExpression()` `group.by = NULL` will now collapse the whole list.

# scRepertoire VERSION 1.7.2

## UNDERLYING CHANGES

* `clonalDiversity()` now has `skip.boots` to stop bootstrapping and downsampling.

# scRepertoire VERSION 1.7.0

## UNDERLYING CHANGES

* Rebumping the version change with new release of Bioconductor.
* Added mean call to the heatmap of `vizGenes()`.
* `filteringMulti` in `combineTCR` now checks to remove list elements with 0 cells.
* Removed `top_n()` call (deprecated), using `slice_max()` without ties.
* Added `arrange()` call during `parseTCR()` to organize the chains.
* Corrected the gd flip in combine functions.
* Removed viridis call in `clonalNetwork()` that was leading to errors.
* Matched syntax for strict clonotype in `combineBCR()`.
* Added `group.by` variable to all applicable visualizations.
* Added `return.boots` to `clonalDiversity()` to allow export of all bootstrapped values.

# scRepertoire VERSION 1.5.4

## BUG FIXES

* Modified `grabMeta()` internal function to no longer assume the active identity is clusters.
* `checkBlanks()` now checks if a blank was detected before trying to remove it.
* `clonalNetwork()` automatically resulted in default error message - bug now removed.
* `clonalBias()` now adds z-score of bias when matrix is exported. `exportTable` parameter is now fixed.

# scRepertoire VERSION 1.5.3

## UNDERLYING CHANGES

* Added `loadContigs()` for non-10X formatted single-cell data.
* Removed `combineTRUST4`, superseded by `loadContigs()`.
* `combineTCR()` now allows for > 3 recovered TCRs per barcode.
* Re-added filtering steps to `combineTCR()`, will detect if data is from 10X and automatically remove nonproductive or multi chains.
* Updated `parseTCR()` to include evaluation for gamma/delta chains.

# scRepertoire VERSION 1.5.2

## UNDERLYING CHANGES

* `highlightClones()` now returns the specific clones instead of clonotype 1, etc.
* `clonalCompare()` `numbers` parameter now for group-wide numbers and not overall top X numbers.
* Fixed issue with `clonalDiversity()` that caused errors when `group.by` parameter was used.
* Modified `parseBCR()` to reduce complexity and assume lambda >> kappa.
* Fixed `clonalCluster()` function broken with Seurat Objects.
* `checkContigs` now ensures data frames and that "" are converted into NAs.
* Modified `makeGenes()` internal function changing `na.omit` to `str_replace_na()` and separating the BCR calls by chain to prevent combination errors.

# scRepertoire VERSION 1.3.5

## UNDERLYING CHANGES

* Modified `parseBCR()` to check for contents of the chains. Resolved issue with placing light chain into heavy chain slots when 2 contigs are present.
* Updated `checkBlanks()` to include NA evaluation and placed the check in all viz functions.
* Added `clonalNetwork()` function.
* Modified diversity visualization to remove outliers and place graphs on a single line.
* Modified `clonalOverlay()` to use new internal `getCoord()` function like `clonalNetwork()`.
* Added `threshold` parameter to `clonalSizeDistribution()`.
* Added support for single-cell objects to `clonalCluster()`.

# scRepertoire VERSION 1.3.4

## UNDERLYING CHANGES

* Modified `clonalCluster()` and `combineBCR()` to speed up comparison and use less memory.
* `filteringMulti` now isolates the top contig by chain, then for barcodes with chains > 2, isolates the top expressing chains.
* Modified `makeGenes()` internal function to use `str_c()`.
* Added `threshold` parameter to `combineTRUST4` for B cell manipulation.
* Changed `combineTCR` function to prevent cell type mix up.
* `vizGenes()` can now look at other component genes of the receptor and `separate` parameter replaced by `y.axis`.
* Added `clonalBias()` function for inter-cluster comparison.
* Fixed `clonalCluster()` and `combineBCR()` assumption that you will have unrelated clones.

# scRepertoire VERSION 1.3.3

## UNDERLYING CHANGES

* `combineBCR()` auto naming function updated to actually name the list elements.
* Added `createHTOContigList()` function to create contig list of multiplexed experiments. Fixed issue with `groupBy` variable.
* Added Inverse Pielou metric to diversity call.
* Added `include.na` and `split.by` to `clonalOccupy()` and changed labeling depending on frequency vs proportion.
* Added support for single-cell objects for most visualizations, list organizing by single-cell object can be called using `split.by` variable.
* All `group` and `groupBy` parameters are now `group.by`.

# scRepertoire VERSION 1.3.2

## UNDERLYING CHANGES

* Added `dot.size` parameter to `clonalScatter()`.
* `filteringMulti` now subsets clonotypes with contigs >= 2, to prevent 2 of the same chains.
* Changed how coldata is added to SCE objects using merge instead of union.
* Can now add BCR and TCR simultaneously by making large list.
* Scatter plotting code allows for user to select `dot.size` as a variable on the x or y axis.
* Removed `regressClonotype` function - too many dependencies required.
* Added `chain` option to visualizations and `combineExpression()` to allow users to facilitate single chains - removed `chain` option from `combineTCR`/`combineBCR`/`combineTRUST4`.
* Added NA filter to `combineTCR`/`combineBCR`/`combineTRUST4` for cell barcodes with only NA values.
* Added NA filter to `expression2List()` for cells with NA clonotypes.
* Updated `vizGenes()` to order genes automatically by highest to lowest variance.
* Updated `vizGenes()` to pull the correct genes based on selection.
* Updated parse method for V/J/D placement in TRB/Heavy chains.
* Simplified `clonalDiversity()` to allow for more options in organizing plot and box plots.
* `combineExpression()` adds the `groupBy` variable to Frequency, allowing for multiple calculations to be saved in the meta data.

# scRepertoire VERSION 1.2.3

## UNDERLYING CHANGES

* Changed the access of the sample data to github.io repo.

# scRepertoire VERSION 1.2.2

## DEPRECATED AND DEFUNCT

* Removed Startrac-based functions to pass build on Bioconductor.
* Deprecated `StartracDiversity()`.

# scRepertoire VERSION 1.2.0

## UNDERLYING CHANGES

* Added support for `SingleCellExperiment` format.

## DEPRECATED AND DEFUNCT

* Deprecated `combineSeurat` in favor of `combineExpression()`.
* Deprecated `seurat2List` in favor of `expression2List()`.

# scRepertoire VERSION 1.1.4

## NEW FEATURES

* Added `proportion` to `combineExpression()` function.
* Added `clonalCluster()` and `clonalOverlay()` functions.

## UNDERLYING CHANGES

* Replaced `hammingCompare` with `lvCompare` to enable superior clonotype calling in `combineBCR()`.
* Added downsampling to diversity calculations.
* Fixed Clonal Overlap Coefficient issue - was comparing unique barcodes and not clonotypes.
* Added `checkBlanks()` function to remove list elements without clonotypes.
* Re-added Startrac metrics by stripping down the package.
* Heavily modified dependencies to reduce total number.

# scRepertoire VERSION 1.0.0

## UNDERLYING CHANGES

* Removed dependencies ggfittext and ggdendrogram.
* `clonalSizeDistribution()` now returns a `plot()` function.

# scRepertoire VERSION 0.99.0

* Initial Bioconductor submission.
* Added `getCirclize()`, `exportTable` to visualization functions.
* Added `screp_example` data to package.
* Added `vizGenes()` function and support for monocle in `combineExpression()`.
