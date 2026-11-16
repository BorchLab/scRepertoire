# test script for clonalSizeDistribution.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list,
                      samples = c("P17B", "P17L", "P18B", "P18L",
                                  "P19B", "P19L", "P20B", "P20L"))
combined <- addVariable(combined,
                        variable.name = "Type",
                        variables = rep(c("B", "L"), 4))

test_that("Output formats (distance matrix vs. plot) are correct", {
  # Test that export.table = TRUE returns a matrix
  dist_matrix <- clonalSizeDistribution(combined[1:2], export.table = TRUE)
  expect_true(is.matrix(dist_matrix))
  
  # Test that the default behavior returns a ggplot object
  plot_output <- clonalSizeDistribution(combined[1:2])
  expect_s3_class(plot_output, "ggplot")
})

test_that("Distance matrix has correct properties", {
  dist_matrix <- clonalSizeDistribution(combined[1:3], export.table = TRUE)
  num_samples <- length(combined[1:3])
  
  # Matrix should be square and have correct dimensions
  expect_equal(nrow(dist_matrix), num_samples)
  expect_equal(ncol(dist_matrix), num_samples)
  
  # Should have correct row and column names
  expect_equal(rownames(dist_matrix), names(combined)[1:3])
  expect_equal(colnames(dist_matrix), names(combined)[1:3])
  
  # Diagonal should be zero
  expect_true(all(diag(dist_matrix) == 0))
  
  # Matrix should be symmetric
  expect_equal(dist_matrix, t(dist_matrix))
  
  # All values should be non-negative
  expect_true(all(dist_matrix >= 0))
})

test_that("`group.by` parameter correctly aggregates data", {
  dist_matrix_grouped <- clonalSizeDistribution(combined, group.by = "Type", export.table = TRUE)
  num_groups <- 2
  expect_equal(nrow(dist_matrix_grouped), num_groups)
  expect_equal(ncol(dist_matrix_grouped), num_groups)
  expect_true(all(rownames(dist_matrix_grouped) %in% c("B", "L")))
})

test_that("`method` parameter is passed to hclust correctly", {
  expect_silent(clonalSizeDistribution(combined[1:2], method = "ward.D2"))
  expect_silent(clonalSizeDistribution(combined[1:2], method = "complete"))
  expect_silent(clonalSizeDistribution(combined[1:2], method = "single"))
})

test_that("`clone.call` and `chain` parameters execute correctly", {
  expect_silent(
    res_aa <- clonalSizeDistribution(combined[1:2], clone.call = "aa", export.table = TRUE)
  )
  expect_true(is.matrix(res_aa))
  
  expect_silent(
    res_tra <- clonalSizeDistribution(combined[1:2], chain = "TRA", export.table = TRUE)
  )
  expect_true(is.matrix(res_tra))
})

test_that("Plotting output has the correct dendrogram components", {
  plot_output <- clonalSizeDistribution(combined[5:6])
  expect_s3_class(plot_output$layers[[1]]$geom, "GeomSegment")
  expect_s3_class(plot_output$layers[[2]]$geom, "GeomText")
  expect_s3_class(plot_output$layers[[3]]$geom, "GeomPoint")
  expect_s3_class(plot_output$coordinates, "CoordFlip")
  expect_true(plot_output$scales$scales[[1]]$trans$name == "reverse")
})

#################################################################
# Section 0: Testing Internal GPD Distribution Functions
#################################################################
# These pin the in-house GPD functions to closed form rather than to any
# external package, so they stand on their own now that evmix is gone.

test_that("Internal GPD functions reduce to the exponential when xi = 0", {
  u <- 2.5
  sigma <- 1.7
  x <- seq(u + 0.01, u + 20, length.out = 50)
  p <- c(0.001, 0.1, 0.5, 0.9, 0.999)

  expect_equal(.dgpd_internal(x, u = u, sigmau = sigma, xi = 0),
               dexp(x - u, rate = 1 / sigma))
  expect_equal(.pgpd_internal(x, u = u, sigmau = sigma, xi = 0),
               pexp(x - u, rate = 1 / sigma))
  expect_equal(.qgpd_internal(p, u = u, sigmau = sigma, xi = 0),
               u + qexp(p, rate = 1 / sigma))

  # The limiting branch is taken for any |xi| < 1e-6
  expect_equal(.pgpd_internal(x, u = u, sigmau = sigma, xi = 1e-9),
               pexp(x - u, rate = 1 / sigma))
})

test_that("Internal GPD functions match the closed-form GPD when xi != 0", {
  for (params in list(c(u = 0, sigmau = 1, xi = 0.3),
                      c(u = 5, sigmau = 2.5, xi = 0.9),
                      c(u = 100, sigmau = 0.4, xi = -0.25))) {
    u <- params[["u"]]
    sigma <- params[["sigmau"]]
    xi <- params[["xi"]]
    # Stay strictly inside the support, which is bounded above when xi < 0
    top <- if (xi < 0) u - sigma / xi else u + 30
    x <- seq(u + 1e-3, top - 1e-3, length.out = 50)

    for (phiu in c(1, 0.35)) {
      s <- 1 + xi * (x - u) / sigma
      expect_equal(.dgpd_internal(x, u, sigma, xi, phiu = phiu),
                   (phiu / sigma) * s^(-1 / xi - 1))
      expect_equal(.pgpd_internal(x, u, sigma, xi, phiu = phiu),
                   1 - phiu * s^(-1 / xi))

      # The quantile function inverts the distribution function
      p <- .pgpd_internal(x, u, sigma, xi, phiu = phiu)
      expect_equal(.qgpd_internal(p, u, sigma, xi, phiu = phiu), x)
    }
  }
})

test_that("Internal GPD functions handle the support boundaries", {
  u <- 3
  sigma <- 2

  # No mass at or below the threshold
  expect_equal(.dgpd_internal(c(u - 1, u), u = u, sigmau = sigma, xi = 0.2),
               c(0, 0))
  expect_equal(.pgpd_internal(c(u - 1, u), u = u, sigmau = sigma, xi = 0.2),
               c(0, 0))

  # A negative shape bounds the support above at u - sigmau/xi
  xi <- -0.5
  upper <- u - sigma / xi
  expect_equal(.dgpd_internal(upper + 1, u = u, sigmau = sigma, xi = xi), 0)
  expect_equal(.pgpd_internal(upper + 1, u = u, sigmau = sigma, xi = xi), 1)

  # Below the threshold a mixture (phiu < 1) leaves the bulk to the other
  # component, so the value is undefined rather than zero
  expect_true(is.na(.pgpd_internal(u - 1, u = u, sigmau = sigma, xi = 0.2,
                                   phiu = 0.4)))
  expect_true(is.na(.qgpd_internal(0.1, u = u, sigmau = sigma, xi = 0.2,
                                   phiu = 0.4)))

  # log = TRUE and lower.tail = FALSE agree with the untransformed values
  x <- c(4, 7, 12)
  expect_equal(.dgpd_internal(x, u, sigma, 0.2, log = TRUE),
               log(.dgpd_internal(x, u, sigma, 0.2)))
  expect_equal(.pgpd_internal(x, u, sigma, 0.2, lower.tail = FALSE),
               1 - .pgpd_internal(x, u, sigma, 0.2))

  # Arguments recycle to a common length
  expect_length(.dgpd_internal(x, u = c(1, 2, 3), sigmau = sigma, xi = 0.2), 3)
  expect_error(.dgpd_internal(x, u, sigmau = -1, xi = 0.2), "positive")
  expect_error(.qgpd_internal(1.5, u, sigma, 0.2), "\\[0,1\\]")
})

test_that(".rgpd_internal draws from the right distribution", {
  set.seed(42)
  u <- 10
  sigma <- 3
  xi <- 0.25
  r <- .rgpd_internal(2e4, u = u, sigmau = sigma, xi = xi)

  expect_length(r, 2e4)
  expect_true(all(r > u))
  # Empirical quantiles should track the theoretical ones
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(as.vector(quantile(r, probs)),
               .qgpd_internal(probs, u = u, sigmau = sigma, xi = xi),
               tolerance = 0.05)

  # Draws from the bulk fraction are not part of the GPD tail
  set.seed(42)
  mixed <- .rgpd_internal(1000, u = u, sigmau = sigma, xi = xi, phiu = 0.5)
  expect_true(mean(is.na(mixed)) > 0.4 && mean(is.na(mixed)) < 0.6)
})

test_that(".fgpd_internal recovers the generating parameters", {
  set.seed(7)
  u <- 0
  sigma <- 2
  xi <- 0.3
  x <- .rgpd_internal(5e4, u = u, sigmau = sigma, xi = xi)

  fit <- .fgpd_internal(x, u = u)

  expect_type(fit$mle, "double")
  expect_length(fit$mle, 2)
  expect_null(names(fit$mle))
  expect_equal(fit$mle[1], sigma, tolerance = 0.05)
  expect_equal(fit$mle[2], xi, tolerance = 0.1)
  expect_equal(fit$sigmau, fit$mle[1])
  expect_equal(fit$xi, fit$mle[2])
  expect_equal(fit$phiu, 1)
  expect_true(fit$conv)

  # The threshold is strict, so phiu is the observed exceedance fraction
  fit_thresh <- .fgpd_internal(x, u = quantile(x, 0.8))
  expect_equal(fit_thresh$phiu, 0.2, tolerance = 0.01)

  expect_error(.fgpd_internal(c(1, 2, 3), u = 10), "above threshold")
})

test_that(".nlgpd_internal is the negative GPD excess log-likelihood", {
  x <- c(6, 7, 8, 12, 20)
  u <- 5
  sigma <- 2
  xi <- 0.4
  yu <- (x - u) / sigma

  expected <- -(-length(x) * log(sigma) - (1 / xi + 1) * sum(log(1 + xi * yu)))
  expect_equal(.nlgpd_internal(c(sigma, xi), x, u = u), expected)

  # Equivalent to summing the log density of the excesses
  expect_equal(.nlgpd_internal(c(sigma, xi), x, u = u),
               -sum(.dgpd_internal(x, u, sigma, xi, log = TRUE)))

  # Invalid parameters give an infinite penalty, clamped when finitelik = TRUE
  expect_equal(.nlgpd_internal(c(sigma, -5), x, u = u), Inf)
  expect_equal(.nlgpd_internal(c(sigma, -5), x, u = u, finitelik = TRUE), 1e6)
})

#################################################################
# Section 2: Testing Core Modeling Functions
#################################################################

test_that(".fdiscgammagpd returns a correctly structured list", {
  # Generate sample data that looks like clonal counts
  set.seed(42)
  x <- rpois(100, lambda = 2) + 1 
  
  fit <- .fdiscgammagpd(x, useq = 5)
  
  expect_type(fit, "list")
  expect_named(fit, c("x", "shift", "init", "useq", "nllhuseq", "optim", "nllh", "mle", "fisherInformation"))
  expect_type(fit$mle, "double")
  expect_equal(length(fit$mle), 6)
  expect_named(fit$mle, c("phi", "shape", "rate", "thresh", "sigma", "xi"))
})

test_that("NLL functions return a single positive number", {
  # Mock parameters and data
  params <- c(1.5, 2.5) # Log-scale parameters for shape, rate
  data <- c(1, 2, 3, 4)
  
  nll_gamma <- .discgammanll(log(params), dat = data, thresh = 5, phiu = 0.1, shift = 0)
  expect_type(nll_gamma, "double")
  expect_length(nll_gamma, 1)
  expect_true(nll_gamma > 0)
  
  # Mock parameters for GPD
  gpd_params <- c(log(1), 0.1) # log(sigma), xi
  gpd_data <- c(6, 7, 8)
  
  nll_gpd <- .discgpdnll(gpd_params, dat = gpd_data, thresh = 5, phiu = 0.1)
  expect_type(nll_gpd, "double")
  expect_length(nll_gpd, 1)
  expect_true(nll_gpd > 0)
})

test_that("Discrete distribution functions work as expected", {
  # Test that the probability mass functions are non-negative
  ddg <- .ddiscgamma(x = 1:5, shape = 1, rate = 1, thresh = 6, phiu = 0.5, shift = 0)
  expect_true(all(ddg >= 0))
  
  ddgpd <- .ddiscgpd(x = 6:10, thresh = 5, sigma = 1, xi = 0.1, phiu = 0.5)
  expect_true(all(ddgpd >= 0))
  
  # Test the random number generators
  set.seed(42)
  rands <- .rdiscgamma(100, shape = 2, rate = 1, thresh = 10, shift = 0)
  expect_length(rands, 100)
  expect_true(all(rands < 10)) # All values should be below the threshold
  expect_true(all(rands == floor(rands))) # All values should be integers
})

#################################################################
# Section 3: Testing Distance Calculation Functions
#################################################################


test_that(".JS_dist returns a single non-negative number", {
  set.seed(1)
  fit1 <- .fdiscgammagpd(rpois(100, 4) + 1, useq = 10)
  fit2 <- .fdiscgammagpd(rpois(100, 5) + 1, useq = 10)
  
  jsd <- .JS_dist(fit1, fit2, grid = 0:1000)
  
  expect_type(jsd, "double")
  expect_length(jsd, 1)
  expect_gte(jsd, 0) 
})

#################################################################
# Section 4: Testing Threshold Grid Search Robustness
#################################################################

test_that("Samples with a thin tail no longer abort the fit", {
  # P18L and P20L have almost no clones above the default threshold, which
  # drove the optimizer into an invalid region and errored out of optim().
  # Every threshold then failed, which.min() returned integer(0), and the
  # subscript two lines later raised "attempt to select less than one element".
  for (k in c(4, 5)) {
    dist_matrix <- clonalSizeDistribution(combined[seq_len(k)],
                                          export.table = TRUE)
    expect_true(is.matrix(dist_matrix))
    expect_equal(dim(dist_matrix), c(k, k))
    expect_true(all(is.finite(dist_matrix)))
    expect_true(all(dist_matrix >= 0))
    expect_equal(dist_matrix, t(dist_matrix))
    expect_true(all(diag(dist_matrix) == 0))
  }
})

test_that("`threshold` accepts a vector and the grid search picks per sample", {
  dist_matrix <- clonalSizeDistribution(combined, threshold = 1:5,
                                        export.table = TRUE)
  expect_true(is.matrix(dist_matrix))
  expect_equal(dim(dist_matrix), c(8, 8))
  expect_true(all(is.finite(dist_matrix)))

  # Thresholds that cannot be fit are skipped rather than aborting the search
  set.seed(5)
  x <- rpois(500, 3) + 1
  fit <- suppressWarnings(.fdiscgammagpd(x, useq = 1:6))
  expect_length(fit$nllhuseq, 6)
  expect_true(any(is.finite(fit$nllhuseq)))
  expect_equal(fit$mle[["thresh"]],
               (1:6)[which.min(fit$nllhuseq)])
})

test_that("NLL functions penalize invalid parameters instead of erroring", {
  # optim() searches log(sigma) unconstrained, so exp() can underflow to zero
  # or overflow to infinity; that must read as an infinite penalty
  expect_equal(.discgpdnll(c(-800, 0.1), dat = c(6, 7, 8), thresh = 5,
                           phiu = 0.1), Inf)
  expect_equal(.discgpdnll(c(800, 0.1), dat = c(6, 7, 8), thresh = 5,
                           phiu = 0.1), Inf)
  expect_equal(.discgammanll(c(800, 800), dat = c(1, 2), thresh = 5,
                             phiu = 0.1), Inf)

  # A single excess leaves var() undefined, so the starting scale can be NA
  expect_equal(.nlgpd_internal(c(NA_real_, 0.1), c(6, 7), u = 5), Inf)
  expect_equal(.nlgpd_internal(c(NA_real_, 0.1), c(6, 7), u = 5,
                               finitelik = TRUE), 1e6)
})

test_that("An unfittable model reports why rather than a subscript error", {
  # No clone sits above the threshold, so neither component can be identified
  x <- c(rep(0L, 50), 1L, 2L)
  expect_error(suppressWarnings(.fdiscgammagpd(x, useq = 100)),
               "could not be fit at any threshold")

  # Zero-length input must not recycle up to a length-one NA, which would
  # score an unfittable threshold as a perfect fit
  expect_error(.gpd_recycle(numeric(0), 1, 1, 0, 1), "empty")
})
