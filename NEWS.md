# scRepertoire VERSION 2.6.1

Update to 2.6.0 to match Bioconductor Release 3.22 on 2025/10/30

## UNDERLYING CHANGES
* Temporarily-requiring Seurat version < 5.3.0 due to bug in DimPlot()

## BUG FIXES
* Fixed order.by issue in `positionalProperty()`
* Fixed individual chain call for `combineExpression()`
* Fixed issue with removing kmer with ";" in `percentKmer()`

# scRepertoire VERSION 2.5.7

## BUG FIXES
* Fixed chain handling for BCR genes in `percentGeneUsage()` and propagates to wrappers: `percentGenes`, `percentVJ` and `vizGenes`

# scRepertoire VERSION 2.5.6

## BUG FIXES
* Fixed passing group.by for `combineBCR()`

# scRepertoire VERSION 2.5.5

* Update unit tests for ggplot2 v4

# scRepertoire VERSION 2.5.3

## BUG FIXES
* Fixed clonalProportion calculation to use grouping properly during `combineExpression()`
* Fixed immunarch support for `exportClones()` - TRA/Light chain column handling

# scRepertoire VERSION 2.5.2

* Added support for mouse genes in `quietBCRgenes()` and `quietTCRgenes()`

# scRepertoire VERSION 2.5.1

* Introduced pairwise calculations to ```StartracDiversity()```
* Internal function conversion for ```clonalSizeDistribution()``` - remove 
cubature, truncdist and VGAM from dependencies.
* Increased speed of ```clonalSizeDistribution()```

# scRepertoire VERSION 2.5.0

## UNDERLYING CHANGES
* Update/Improve code for ```loadContigs()```
* Consolidated support for discrete AIRR formats under the umbrella of AIRR
* Added `"tcrpheno"` and `immunarch` to ```exportClones()```
* Converted ```exportClones()``` to base R to reduce dependencies
* Added dandelionR vignette to pkgdown site
* Added tcrpheno vignette to pkgdown site
* ```percentAA()``` refactored to minimize dependencies and use immApex 
```calculateFrequency()```
* ```positionalEntropy()``` refactored to minimize dependencies and use immApex
```calculateEntropy()```
* ```clonalDiversity()``` refactored for performance -  it now calculates a 
single diversity metric at a time and includes new  estimators like "gini", 
"d50", and supports hill numbers).
* ```percentKmer()``` refactored to use immApex ```calculateMotif``` for both aa and
nt sequences. No longer calculates all possible motifs, but only motifs present. 
* ```clonalCluster()``` now allows for dual-chain clustering, V/J filtering, 
normalized or straight edit distance calculations, and return of clusters, 
igraph objects or adjacency matrix
* ```combineBCR()``` offers single/dual chain clustering, aa or nt sequences, adaptive filtering of V and J genes and normalized or straight edit distance calculations
* ```percentGeneUsage()``` now is the underlying function for ```percentGenes()```, ```percentVJ()```, 
and ```vizGenes()``` and allows for percent, proportion and raw count quantification.
* Added common theme (internal ```.themeRepertoire()```) to all plots and allow users to pass arguments to it

## BUG FIXES
* ```clonalCompare()``` issue with plotting a 0 row data frame now errors with message
* ```clonalScatter()``` group.by/axes call now works for non-single-cell objects
* Fixed issue with NULL and "none" group.by in ```combineExpression()```
* Allowing multi groupings via x.axis and group.by in ```clonalDiversity()```

# scRepertoire VERSION 2.3.4

## UNDERLYING CHANGES
* Update internal ```.parseContigs``` to function with more complex groupings
* Add ```annotateInvariant()``` functionality for mouse and human TCRs
* Add ```quietTCRgenes()```, ```quietBCRgenes()```, ```quietVDJgenes()```
* Fixed issue with ```clonalCompare()``` assertthat statements
* Started integration with immApex API package

# scRepertoire VERSION 2.3.2

## UNDERLYING CHANGES
* Fixed issue with denominator in ```getCirclize()```
* Fixing chain issue with clonalCompare() - expanded assertthat statement

# scRepertoire VERSION 2.2.1

* Rebasing for the purposes of Bioconductor version 3.20

## NEW FEATURES
* Added support for BCRs for loading ParseBio sequences.
* Added `quietBCRgenes()` and `quietTCRgenes()` from Ibex and Trex and `quietVDJgenes()` as a convenience that runs both. The functions filter out known TCR and/or BCR gene signatures.

## UNDERLYING CHANGES
* Added `Seurat` to the `Suggests` field in the DESCRIPTION file.

# scRepertoire VERSION 2.0.8

## NEW FEATURES
* Added ```getContigDoublets()``` experimental function to identify TCR and BCR doublets as a preprocessing step to ```combineExpression()```
* Added **proportion** argument to ```clonalCompare()``` so that when set to FALSE, the comparison will be based on frequency normalized by per-sample repertoire diversity.

## UNDERLYING CHANGES
* Fixed issue with single chain output for ```clonalLength()```
* Removed unnecessary code remnant in ```clonalLength()```
* Allow one sample to be plotted by ```percentVJ()```
* Fixed issue with ```positionalProperty()``` and exportTable
* Fixed issue with ```loadContigs()``` edge case when TRUST4 data only has 1 row.
* convert documentation to use markdown (`roxygen2md`)
* import `lifecycle`, `purrr`, `withr`
* suppressed "using discrete variable for alpha is not recommended" warning in alluvialClones unit tests.
* Fixed issue with ```clonalCluster()``` and exportGraph = TRUE
* improve performance of ```combineBCR()``` by a constant factor with C++
* Restructured functions to exportTable before plotting

# scRepertoire VERSION 2.0.7

## UNDERLYING CHANGES
* Fixed issue with "group.by" in ```clonalOverlap()```
* Fixed issue with "group.by" in ```clonalCompare()```

# scRepertoire VERSION 2.0.6

## UNDERLYING CHANGES
* Fixed issue with custom column headers for clones

# scRepertoire VERSION 2.0.5

## UNDERLYING CHANGES
* added type checks using assertthat
* updated conditional statements in constructConDFAndparseTCR.cpp
* Fixed issue in ```clonalQuant()``` and factor-based **group.by** variable

# scRepertoire VERSION 2.0.4

## UNDERLYING CHANGES
* ```getCirclize()``` refactored to prevent assumptions and added **include.self** argument
* Added ```.count.clones()``` internal function for ```getCirclize()``` and ```clonalNetwork()```
* Added **order.by** parameter to visualizations to specifically call order of plotting using a vector or can use "alphanumeric" to plot things in order
* Fix issue with ```clonalLength()``` and NA handling
* ```clonalCompare()``` now retains the original clonal info if using **relabel.clones**
* Add Dandelion support in to ```loadContigs()``` and testthat
* Fixed issue with ```positionalProperty()``` assumption that the clones will all have 20 amino acids.
* Fixed issue with ```positionalProperty()``` and removing non-amino acids.
* Fixed IGH/K/L mistaking gene issue in ```vizGenes()```
* Add error message for NULL results in ```clonalCluster()``` with **export.graph = TRUE**
* Fixed issue with "full.clones" missing in ```combineExpression()``` when using 1 chain

# scRepertoire VERSION 2.0.3

## UNDERLYING CHANGES

* Modified support for Omniscope format to allow for dual chains
* Added ParseBio support in to ```loadContigs()``` and testthat
* Added support for productive variable to ```loadContigs()``` for BD, Omniscope, and Immcantation formats
* Replace numerical indexing with name indexing for ```loadContigs()```
* ```combineBCR()``` and ```combineTCR()``` no allow for unproductive contig inclusions with new **filterNonproductive** parameter.
* ```combineBCR()``` will now prompt user if **samples** is not included instead of erroring.
* Added base threshold by length for internal ```.lvCompare()```
* Ensured internal ```.lvCompare()``` only looks at first set of sequences in multi-sequence chain.
* Fixed bug in exporting graph for ```clonaCluster()```
* Fixed conflict in functions between igraph and dplyr packages

# scRepertoire VERSION 2.0.2

## UNDERLYING CHANGES

* ```clonalOccupy()``` rewrite counting and NA handling

# scRepertoire VERSION 2.0.1 

## UNDERLYING CHANGES

* ```clonalOverlay()``` arguments now cutpoint and use cut.category to select either clonalProportion or clonalFrequency

# scRepertoire VERSION 2.0.0 (2024-01-10)

## NEW FEATURES

* Added ```percentAA()```
* Added ```percentGenes()```
* Added ```percentVJ()```
* Added ```percentKmer()```
* Added ```exportClones()``` 
* Added ```positionalEntropy()``` 
* Added ```positionalProperty()``` 
* Changed compareClonotypes to ```clonalCompare()```
* Changed clonotypeSizeDistribution to ```clonalSizeDistribution()```
* Changed scatterClonotypes to ```clonalScatter()```
* Changed quantContig to ```clonalQuant()```
* Changed highlightClonotypes to ```highlightClones()```
* Changed lengthContigs to ```clonalLength()```
* Changed occupiedscRepertoire to ```clonalOccupy()```
* Changed abundanceContig to ```clonalAbundance()```
* Changed alluvialClonotypes to ```alluvialClones()```
* Added features to ```clonalCompare()``` to allow for highlighting sequences, relabeling clonotypes.

## UNDERLYING CHANGES

* Removed internal **.quiet()** function.
* **.theCall()** now allows for a custom header/variable and checks the colnames. 
* Replaced data arguments to be more descriptive: *df* is now *input.data*, *dir* is now *input*, and *sc* is now *sc.data*
* Deep clean on the documentation for each function for increased consistency and explainability
* ```StartracDiversity()``` metric re-implemented to remove startrac-class object intermediary
* Implemented powerTCR locally to reduce dependencies and continue support
* Universalized underlying function language and intermediate variables
* License change to MIT
* **group.by** and **split.by** have been consolidated into single **group.by** parameter
* Added support for Immcantation pipeline, .json, Omniscope, and MiXCR formats for ```loadContigs()```
* Made GitHub.io website for support/vignettes/FAQ
* Restructured NEWS Tracking
* Added testthat for all exported and internal functions
* Fixed issue with ```clonalQuant()``` for instance of **scale** = FALSE and **group.by** being set.
* ```clonalDiversity()``` no longer automatically orders samples.
* Remove **order** parameter from ```clonalQuant()```, ```clonalLength()```, and ```clonalAbundance()```
* **x.axis** parameter in ```clonalDiversity()``` separated from **group.by** parameter
* filtering chains will not eliminate none matching chains.

## DEPRECATED AND DEFUNCT

* Deprecate stripBarcodes()
* Deprecate expression2List() (now only an internal function).
* Deprecate checkContigs()

# scRepertoire VERSION 1.11.0

* Rebasing for the purposes of bioconductor version

# scRepertoire VERSION 1.7.5

* Fixed combineBCR() to allow for non-related sequences

# scRepertoire VERSION 1.7.4

## NEW FEATURES

* checkContigs() function to quantify the percentages of NA values by genes or sequences
* exportClones to clonalNetwork() to isolate clones shared across identities.

## UNDERLYING CHANGES

* Fix issue with clonalDiversity() and skipping boots
* Fixing underlying assumptions with clonalBias()
* Adding reads variable to parseAIRR
* Fixing handling of samples parameter in combine contain functions
* removed need to relevel the cloneType factor after combineExpression()
* set up lapply() for combineBCR() and clusterTCR() - no more pairwise distance matrix calculation
* loadContigs() support for data.frames or lists of contigs
* Added examples for loadContigs() to test function
* Removed requirement for T cell type designation - will combine A/B and G/D simultaneously
* Updated combineBCR() to chunk nucleotide edit distance calculations by V gene and give option to skip edit distance calculation with call.related.clones = FALSE
* Updated clusterTCR() to use lvcompare() function and base edit distances of V gene usage. 


# scRepertoire VERSION 1.7.3

* Fix misspellings for parse contains functions
* Optimize WAT3R and 10x loadContigs()
* Remove combineTRUST4 - superseded by loadContigs()
* Added support of TRUST4 for combineBCR()
* Added support for BD in loadContigs()
* loadContigs() TRUST4 parsing allows for all NA values in a chain
* combineExpression() group.by = NULL will now collapse the whole list. 


# scRepertoire VERSION 1.7.2

* ClonalDiversity() now has skip.boots to stop bootstrapping and downsampling


# scRepertoire VERSION 1.7.0

* Rebumping the version change with new release of Bioconductor
* Added mean call to the heatmap of vizGenes()
* To combineTCR, filteringMulti now checks to remove list elements with 0 cells.
* Removed top_n() call as it is now deprecated, using slice_max() without ties.
* Add arrange() call during parseTCR() to organize the chains
* Correct the gd flip in the combineContig and subsequent functions
* Removed viridis call in the clonalNetwork() function that was leading to errors
* Matched syntax for strict clonotype in combineBCR()
* Added group.by variable to all applicable visualizations
* Added return.boots to clonalDiversity(), allow for export of all bootstrapped values


# scRepertoire VERSION 1.5.4

* modified grabMeta() internal function to no longer assume the active Identity is clusters. 
* checkBlanks() now checks if a blank was detected before trying to remove it
* clonalNetwork() automatically resulted in default error message, bug now removed.
* clonotypeBias now adds z-score of bias when matrix is exported. exportTable parameter is now fixed.

# scRepertoire VERSION 1.5.3

* Added loadContigs for non-10X formatted single-cell data
* removed combineTRUST4, superseded by loadContigs
* combineTCR() now allows for > 3 recovered TCRs per barcode
* Readded the filtering steps to combineTCR(), will detect if data is from 10X and automatically remove nonproductive or multi chains. 
* Updated parseTCR() to include evaluation for gamma/delta chains. 

# scRepertoire VERSION 1.5.2

* Arbitrarily numbering system to match new bio conductor dev version
* highlightClonotypes() now returns the specific clones instead of clonotype 1, ...
* compareClonotypes numbers parameter now for group-wide numbers and not overall top X numbers
* Fixed issue with clonalDiversity that cause errors when group.by parameter was used. 
* modified parseBCR() to reduce complexity and assume lambda >> kappa
* fixed clusterTCR() function broken with Seurat Objects
* checkContigs no ensures data frames and that "" are converted into NAs
* modified makeGenes() internal function changing na.omit to str_replace_na() and separating the BCR calls by chain to prevent combination errors. 

# scRepertoire VERSION 1.3.5

* Modified parseBCR() to check for contents of the chains. Resolve issue with placing light chain into heavy chain slots when 2 contigs are present. 
* Updated checkBlanks to include NA evaluation and placed the check in all .viz functions
* Added clonalNetwork() function
* Modified diversity visualization to remove outliers and place graphs on a single line
* Modified clonalOverlay() to use new internal getCoord() function like clonalNetwork()
* Added threshold parameter to clonesizeDistribution()
* Added support for single-cell objects to clusterTCR()


# scRepertoire VERSION 1.3.4

* Modification in clusterTCR() and combineBCR() to speed up the comparison and use less memory
* FilteringMulti, now isolates the top contig by chain, then for barcodes with chains > 2, isolates the top expressing chains. This substantially increases the speed of the filtering step. 
* Modified makeGenes() internal function to use strings str_c()
* Added threshold parameter to combineTRUST4 for B cell manipulation 
* Changed combineTCR function to prevent cell type mix up and clarified in function documentation.
* vizGenes can now be used to look at other component genes of the receptor and "separate" parameter was replaced by "y.axis" parameter.
* Added clonotypeBias() function for inter-cluster comparison.
* Fixed clusterTCR() and combineBCR() assumption that you will have unrelated clones. 


# scRepertoire VERSION 1.3.3

* CombineBCR() auto naming function updated to actually name the list elements. 
* Added createHTOContigList() function to create contain list of multiplexed experiments. Fixed issue with groupBy variable
* Added Inv.Pielou matric to diversity call - this is essentially 1-shannon/ln(length). Due to the bootstrapping the length with be constant.
* Added include.na and split.by to occupiedscRepertoire and changed labeling depending on frequency vs proportion
* Added support for single-cell objects for most visualizations, list organizing by single-cell object can be called using split.by variable
* All group and groupBy parameters are now group.by.


# scRepertoire VERSION 1.3.2

* This is the new numbering scheme apologies - we are all up-to-date now and now cell ranger >= 5 will # work on bioconductor, so let's all just take that as a win.
* added dot.size parameter to scatterClonotype
* filteringMulti now subsets clonotypes with contains >=2, to prevent 2 of the same chains
* changed how coldata is added to SCE objects using merge instead of union
* Can now add BCR and TCR simultaneously by making large list
* scatter plotting code is not so ugly and allows for user to select dot.size as a variable on the x or y axis
* Removed regressClonotype function - too many dependencies required, adding an additional vignette to go through the process
* Added chain option to visualizations and combineExpression to allow users to facilitate single chains - removed chain option from combineTCR/BCR/TRUST4 (the combined object will have both chains no matter what)
* Added NA filter to combineTCR/BCR/TRUST4 for cell barcodes with only NA values
* Added NA filter to expression2List() for cells with NA clonotypes.
* Updated VizGene to order the genes automatically by highest to lowest variance
* Updated VizGene to pull the correct genes based on selection
* Updated parse method - old version had issue with place V-->J-->D in the TRB/Heavy chains 
* Simplified the clonalDiversity() to allow for more options in organizing plot and box plots. 
* CombineExpression() adds the groupBy variable to Frequency, allowing for multiple calculations to be saved in the meta data. 
* Default color scheme now uses viridis plasma, because it I am on transfusion medicine.

# scRepertoire VERSION 1.2.2

* added the combineTRUST4 function to parse contigs from TUST4 pipeline
* added the filter of contigs by chain in the combineTCR, combineBCR, and combineTRUST4 functions
* no longer require the ID in the combineTCR/BCR/TRUST4 functions
* added jaccard index for overlap analysis
* replaced vizVgene with vizGene - allowing users to look at any gene in the combinedContig object
* Fixed coloring scale on the overlap analysis
* Added regressClonotype function using harmony to remove the clonotype effect on feature space
* allowed occupiedRepertoire to use proportion.
* added scatterClonotype function to Viz.R

# scRepertoire VERSION 1.2.1

* number of changes to the parseTCR/BCR functions to limit assumptions
* Changed grabMeta to include assessment of colnames
* fixed lengthDF handling of single chains with multi chains stored - ;
* Added labels to alluvialClonotype and occupiedClonotype plotting

# scRepertoire VERSION 1.1.4

* replaced hammingCompare with lvCompare to enable superior clonotype calling in combineBCR function.
* added proportion to combineExpression() function so users no longer need to know absolute frequencies when combining the contiguous information. 
* added clusterTCR() and clonalOverlay() functions.
* added downsampling to the diversity calculations
* replaced hammingCompare with lvCompare to enable superior clonotype calling in combineBCR function.
* added proportion to combineExpression() function so users no longer need to know absolute frequencies when combining the contiguous information. 
* added clusterTCR() and clonalOverlay() functions.
* added downsampling to the diversity calculations
* Clonal Overlap Coefficient issue fixed, was comparing unique barcodes and not clonotypes
* Added function checkBlanks to remove list elements without clonotypes, this prevents errors for visualizations
* Re-added Startrac metrics by stripping down the package and adding it piecemeal
* Heavily modified dependencies to reduce total number 

# scRepertoire VERSION 1.0.0

* removed dependencies ggfittext and ggdendrogram
* clonesizeDistribution now returns a plot() function


# scRepertoire VERSION 0.99.18

* Updated author information in the vignette


# scRepertoire VERSION 0.99.17

* Updated NEWS formatting
* Edited DESCRIPTION to Single Cell Experiment R package
* Updated information in the vignette


# scRepertoire VERSION 0.99.16

* Added ```getCirclize()```


# scRepertoire VERSION 0.99.15

* Modified numerator for index function


# scRepertoire VERSION 0.99.14

* Removed bracket from indexing function

# scRepertoire VERSION 0.99.13

* Added exportTable to remaining viz functions
* Modified morisita index to correct error

# scRepertoire VERSION 0.99.12

* Reducing the size of the screp_example to fulfill < 5 mB requirement. Randomly samples 100 cells and removed RNA counts from Seurat object

# scRepertoire VERSION 0.99.11

* Updated compareClonotype to allow for clonotype comparisons

# scRepertoire VERSION 0.99.10

* Bioconductor did not detect the version update.

# scRepertoire VERSION 0.99.9

* Bioconductor had no love - changed the Seurat package to imports instead of required, see if that will address the compiling issue that results in a killed: 9 error.

# scRepertoire VERSION 0.99.8

* Passed checks on system, let's see how much bioconductor hates it

# scRepertoire VERSION 0.99.7

* But really this time, changed the colData import

# scRepertoire VERSION 0.99.6

* Changed colData import

# scRepertoire VERSION 0.99.5

* Added screp_example data to package
* Added visVgene function for visualizing the distribution of V genes in TCR
* Added support for monocle to combineExpression function
* Updated documentation for combineTCR() and combineBCR()
* Updated documentation to utilize SingleCellExperiment formats
* Updated Vignette to utilize SingleCellExperiment formats
* Added Author information to vignette
* Add intro and conclusion to vignette
* Removed html knitted vignette
* Removed descriptive code snippets

# scRepertoire VERSION 0.99.4

* Modified expression2List() to allow for variables across meta data

# scRepertoire VERSION 0.99.1

* Changed R (>= 3.6) to R (>= 4.0)


# scRepertoire VERSION 0.99.0

* Changed DESCRIPTION version to 0.99.0
* Removed file seurat_example.rda, accidentally committed
* Deleted git attributes
* reduced Seurat object size for alluvialClonotype in vignette
* Changed the alluvialClonotype assessment to account for only 1 condition


# scRepertoire VERSION 1.2.3

* Changed the access of the sample data to github.io repo:
readRDS(url("https://ncborcherding.github.io/vignettes/scRepertoire_example.rds"))


# scRepertoire VERSION 1.2.2

* Removed Startrac-based functions in order to pass build on Bioconductor.

DEPRECATED AND DEFUNCT

* Deprecate StartracDiversity()


# scRepertoire VERSION 1.2.0

SIGNIFICANT USER-VISIBLE CHANGES

* Added support for ```SingleCellExperiment``` format.


DEPRECATED AND DEFUNCT

* Deprecate combineSeurat in favor or combineExpression().
* Deprecate seurat2List in favor of expression2List().
