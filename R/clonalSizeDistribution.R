#' Plot powerTCR Clustering Based on Clonal Size
#'
#' This function produces a hierarchical clustering of clones by sample 
#' using discrete gamma-GPD spliced threshold model. If using this 
#' model please read and cite powerTCR (more info available at 
#' [PMID: 30485278](https://pubmed.ncbi.nlm.nih.gov/30485278/)).
#' 
#' @details
#' The probability density function (pdf) for the **Generalized Pareto Distribution (GPD)** is given by:
#'  \deqn{f(x|\mu, \sigma, \xi) = \frac{1}{\sigma} \left( 1 + \xi \left( \frac{x - \mu}{\sigma} \right) \right)^{-\left( \frac{1}{\xi} + 1 \right)}}
#' 
#' Where:
#' \itemize{
#'   \item{\eqn{\mu} is a location parameter}
#'   \item{\eqn{\sigma > 0} is a scale parameter}
#'   \item{\eqn{\xi} is a shape parameter}
#'   \item{\eqn{x \ge \mu} if \eqn{\xi \ge 0} and \eqn{\mu \le x \le \mu - \sigma/\xi} if \eqn{\xi < 0}}
#' }
#'               
#' The probability density function (pdf) for the **Gamma Distribution** is given by:
#' \deqn{f(x|\alpha, \beta) = \frac{x^{\alpha-1} e^{-x/\beta}}{\beta^\alpha \Gamma(\alpha)}}
#' 
#' Where:
#' \itemize{
#'   \item{\eqn{\alpha > 0} is the shape parameter}
#'   \item{\eqn{\beta > 0} is the scale parameter}
#'   \item{\eqn{x \ge 0}}
#'   \item{\eqn{\Gamma(\alpha)} is the gamma function of \eqn{\alpha}}
#' }
#' 
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list,
#'                        samples = c("P17B", "P17L", "P18B", "P18L",
#'                                    "P19B","P19L", "P20B", "P20L"))
#' 
#' # Using clonalSizeDistribution()
#' clonalSizeDistribution(combined,
#'                        clone.call = "strict",
#'                        method="ward.D2")
#'
#' @param input.data The product of [combineTCR()],
#' [combineBCR()], or [combineExpression()].
#' @param clone.call Defines the clonal sequence grouping. Accepted values
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL`,
#' `IGK`, `Light` (for both light chains), or `both` (for TRA/B and Heavy/Light).
#' @param threshold Numerical vector containing the thresholds
#' the grid search was performed over.
#' @param method The clustering parameter for the dendrogram.
#' @param group.by A column header in the metadata or lists to group the analysis
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed as
#' by list element or active identity in the case of single-cell objects.
#' @param export.table If `TRUE`, returns a data frame or matrix of the results
#' instead of a plot.
#' @param palette Colors to use in visualization - input any
#' [hcl.pals][grDevices::hcl.pals].
#' @param cloneCall \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `clone.call` instead.
#' @param exportTable \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use `export.table` instead.
#' @param ... Additional arguments passed to the ggplot theme
#'
#' @importFrom ggdendro dendro_data segment label
#' @importFrom stats hclust optim pgamma as.dist as.formula qgamma runif xtabs
#' @export
#' @concept Visualizing_Clones
#' @return A ggplot object visualizing dendrogram of clonal size distribution
#'  or a data.frame if `exportTable = TRUE`.
#' @author Hillary Koch
#'
clonalSizeDistribution <- function(input.data,
                                   clone.call = NULL,
                                   chain = "both",
                                   method = "ward.D2",
                                   threshold = 1,
                                   group.by = NULL,
                                   export.table = NULL,
                                   palette = "inferno",
                                   # Deprecated arguments
                                   cloneCall = NULL,
                                   exportTable = NULL,
                                   ...) {

  # Handle deprecated arguments
  clone.call <- .deprecate_arg(cloneCall, clone.call, "cloneCall", "clone.call",
                               "clonalSizeDistribution", default = "strict")
  export.table <- .deprecate_arg(exportTable, export.table, "exportTable", "export.table",
                                 "clonalSizeDistribution", default = FALSE)

  x <- y <- xend <- yend <- NULL
  input.data <- .dataWrangle(input.data,
                             group.by,
                             .theCall(input.data, clone.call,
                                      check.df = FALSE, silent = TRUE),
                             chain)
  clone.call <- .theCall(input.data, clone.call)
  sco <- .is.seurat.or.se.object(input.data)
  if(!is.null(group.by) & !sco) {
    input.data <- .groupList(input.data, group.by)
  }
  data <- bind_rows(input.data)
  unique_df<- unique(data[,clone.call])

  # Create long-format summary table
  summary_df <- dplyr::bind_rows(input.data, .id = "sample") %>%
    dplyr::group_by(sample, .data[[clone.call]]) %>%
    dplyr::summarise(Freq = dplyr::n(), .groups = 'drop')
  wide_matrix <- xtabs(as.formula(paste("Freq ~", clone.call, "+ sample")),
                       data = summary_df)
  Con.df <- as.data.frame.matrix(wide_matrix)
  Con.df[[clone.call]] <- rownames(Con.df)
  rownames(Con.df) <- NULL
  Con.df <- Con.df[, c(clone.call, setdiff(names(Con.df), clone.call))]
  
  # Fit models
  list <- lapply(seq_len(ncol(Con.df))[-1], function(x) {
    withCallingHandlers(
      suppressWarnings(.fdiscgammagpd(Con.df[,x], useq = threshold)),
      error = function(e) {
        stop("clonalSizeDistribution() could not fit the clonal size model ",
             "for \"", colnames(Con.df)[x], "\": ", conditionMessage(e),
             call. = FALSE)
      })
  })
  
  names(list) <- names(input.data)
  grid <- 0:10000
  distances <- .get_distances(list, grid, modelType="Spliced")
  mat_melt <- dendro_data(hclust(as.dist(distances), method = method), 
                          type = "rectangle")
  
  #Plotting
  plot <- ggplot() + 
            geom_segment(data = segment(mat_melt), 
                         aes(x = x, 
                             y = y, 
                             xend = xend, 
                             yend = yend)) +
            geom_text(data = label(mat_melt), 
                              aes(x = x, y = -0.02, 
                                  label = label, 
                                  hjust = 0), 
                      size = 4) +
            geom_point(data = label(mat_melt), 
                              aes(x = x, 
                                  y = -0.01, 
                                  color = as.factor(label)), 
                       size = 2) + 
            coord_flip() +
            scale_y_reverse(expand = c(0.2, 0)) + 
            scale_color_manual(values = .colorizer(palette, nrow(label(mat_melt)))) + 
            .themeRepertoire(..., grid_lines = "none") + 
            guides(color = "none") + 
            theme(axis.title = element_blank(), 
                  axis.ticks.y = element_blank(), 
                  axis.text.y = element_blank()) 
  
  if (export.table) {
    return(distances)
  }
  return(plot)
}

#################################################################
# Section 0: Internal GPD Distribution Functions
#################################################################
# Density, distribution, quantile, random generation and maximum-likelihood
# fitting for the Generalized Pareto Distribution (GPD). These replace
# evmix::dgpd/pgpd/qgpd/rgpd/fgpd, dropped as a dependency ahead of evmix's
# CRAN deprecation.
#
# The parameterization follows evmix: threshold u, scale sigmau > 0, shape xi,
# and tail fraction phiu, where the GPD occupies the upper phiu of the
# distribution. Writing the scaled excess y = (x - u)/sigmau and s = 1 + xi*y,
#
#   F(x) = 1 - phiu * s^(-1/xi)                        for x > u
#   f(x) = (phiu/sigmau) * s^(-1/xi - 1)               for x > u, s > 0
#   x_p  = u + (sigmau/xi) * (((1-p)/phiu)^(-xi) - 1)   for p > 1 - phiu
#
# with the xi -> 0 (exponential) limits exp(-y), (phiu/sigmau)*exp(-y) and
# u - sigmau*log((1-p)/phiu). As in evmix the limiting branch is taken when
# abs(xi) < 1e-6, and the parameters are recycled to a common length.

# Recycle the GPD arguments to a common length
.gpd_recycle <- function(q, u, sigmau, xi, phiu) {
  if (length(q) == 0L) {
    # Recycling an empty vector against scalar parameters would silently
    # produce a length-one NA. Callers treat a failure here as "no data above
    # this threshold", which is the right reading
    stop("quantiles input is empty.")
  }
  if (!all(is.finite(sigmau)) || any(sigmau <= 0)) {
    stop("sigmau must be finite and positive.")
  }
  if (any(!is.finite(phiu)) || any(phiu < 0) || any(phiu > 1)) {
    stop("phiu must be in [0,1].")
  }
  n <- max(length(q), length(u), length(sigmau), length(xi), length(phiu))
  list(n      = n,
       q      = rep(q,      length.out = n),
       u      = rep(u,      length.out = n),
       sigmau = rep(sigmau, length.out = n),
       xi     = rep(xi,     length.out = n),
       phiu   = rep(phiu,   length.out = n))
}

# GPD Density
.dgpd_internal <- function(x, u = 0, sigmau = 1, xi = 0, phiu = 1,
                           log = FALSE) {
  pars <- .gpd_recycle(x, u, sigmau, xi, phiu)
  x <- pars$q
  u <- pars$u
  sigmau <- pars$sigmau
  xi <- pars$xi
  phiu <- pars$phiu
  x[is.infinite(x)] <- NA

  yu <- (x - u) / sigmau
  syu <- 1 + xi * yu
  # Strictly greater than the threshold, so the density at x == u is zero
  inSupport <- (yu > 0) & (syu > 0)

  d <- x
  d[which(!inSupport)] <- -Inf

  expBranch <- abs(xi) < 1e-6
  idx <- which(expBranch & inSupport)
  if (length(idx) > 0L) {
    d[idx] <- -log(sigmau[idx]) - yu[idx]
  }
  idx <- which(!expBranch & inSupport)
  if (length(idx) > 0L) {
    d[idx] <- -log(sigmau[idx]) - (1 / xi[idx] + 1) * log(syu[idx])
  }

  d <- d + log(phiu)
  if (!log) {
    d <- exp(d)
  }
  d
}

# GPD Distribution Function
.pgpd_internal <- function(q, u = 0, sigmau = 1, xi = 0, phiu = 1,
                           lower.tail = TRUE) {
  pars <- .gpd_recycle(q, u, sigmau, xi, phiu)
  q <- pars$q
  u <- pars$u
  sigmau <- pars$sigmau
  xi <- pars$xi
  phiu <- pars$phiu
  q[is.infinite(q)] <- NA

  yu <- pmax(q - u, 0) / sigmau
  syu <- pmax(1 + xi * yu, 0)
  above <- yu > 0

  p <- q
  below <- which(!above)
  # Below the threshold the GPD alone (phiu == 1) has no mass; in a mixture
  # the bulk model supplies it, so the value is left undefined
  p[below] <- ifelse(phiu[below] == 1, 0, NA)

  expBranch <- abs(xi) < 1e-6
  idx <- which(expBranch & above)
  if (length(idx) > 0L) {
    p[idx] <- 1 - exp(-yu[idx])
  }
  idx <- which(!expBranch & above)
  if (length(idx) > 0L) {
    p[idx] <- 1 - syu[idx]^(-1 / xi[idx])
  }

  p <- 1 - (1 - p) * phiu
  if (!lower.tail) {
    p <- 1 - p
  }
  p
}

# GPD Quantile Function
.qgpd_internal <- function(p, u = 0, sigmau = 1, xi = 0, phiu = 1,
                           lower.tail = TRUE) {
  if (any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("p must be in [0,1].")
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  pars <- .gpd_recycle(p, u, sigmau, xi, phiu)
  p <- pars$q
  u <- pars$u
  sigmau <- pars$sigmau
  xi <- pars$xi
  phiu <- pars$phiu

  above <- p > (1 - phiu)
  quant <- p
  below <- which(!above)
  quant[below] <- ifelse(phiu[below] == 1, u[below], NA)

  expBranch <- abs(xi) < 1e-6
  idx <- which(expBranch & above)
  if (length(idx) > 0L) {
    quant[idx] <- u[idx] - sigmau[idx] * log((1 - p[idx]) / phiu[idx])
  }
  idx <- which(!expBranch & above)
  if (length(idx) > 0L) {
    quant[idx] <- u[idx] + sigmau[idx] *
      (((1 - p[idx]) / phiu[idx])^(-xi[idx]) - 1) / xi[idx]
  }
  quant
}

# GPD Random Generation
#' @importFrom stats rexp
.rgpd_internal <- function(n = 1, u = 0, sigmau = 1, xi = 0, phiu = 1) {
  if (any(xi == 1)) {
    stop("shape cannot be 1.")
  }
  n <- max(n, length(u), length(sigmau), length(xi), length(phiu))
  pars <- .gpd_recycle(rep(NA_real_, n), u, sigmau, xi, phiu)
  u <- pars$u
  sigmau <- pars$sigmau
  xi <- pars$xi
  phiu <- pars$phiu

  # Draws below the tail fraction belong to the bulk and are left as NA
  inTail <- runif(n) < phiu
  r <- rep(NA_real_, n)

  expBranch <- abs(xi) < 1e-6
  idx <- which(expBranch & inTail)
  if (length(idx) > 0L) {
    # evmix::rgpd passes the threshold where rexp() expects a count; the
    # number of draws is what is wanted here
    r[idx] <- u[idx] + sigmau[idx] * rexp(length(idx))
  }
  idx <- which(!expBranch & inTail)
  if (length(idx) > 0L) {
    r[idx] <- u[idx] + sigmau[idx] *
      (runif(length(idx))^(-xi[idx]) - 1) / xi[idx]
  }
  r
}

# Negative Log-Likelihood of the GPD Excesses Above the Threshold
#
#   l(sigmau, xi) = -nu*log(sigmau) - (1/xi + 1)*sum(log(1 + xi*y))
#                   + nu*log(phiu)
#
# with the xi -> 0 form -nu*log(sigmau) - sum(y) + nu*log(phiu). The
# likelihood is -Inf outside the support or for invalid parameters, which
# `finitelik` clamps to a large finite value for the gradient-based optimizers
# that cannot step away from an infinity.
.nlgpd_internal <- function(pvector, x, u = 0, phiu = 1, finitelik = FALSE) {
  sigmau <- pvector[1]
  xi <- pvector[2]

  x[!is.finite(x)] <- NA
  xu <- x[which(x > u)]
  nu <- length(xu)

  # A lone excess leaves var() undefined, so the starting scale can arrive as
  # NA; report an infinite penalty rather than letting it reach the test below
  if (!is.finite(sigmau) || !is.finite(xi) || nu == 0L) {
    nllh <- Inf
    if (finitelik) nllh <- 1e6
    return(nllh)
  }

  yu <- (xu - u) / sigmau
  syu <- 1 + xi * yu

  if (min(syu) <= 0 || sigmau <= 0 || phiu <= 0 || phiu > 1) {
    loglik <- -Inf
  } else if (abs(xi) < 1e-6) {
    loglik <- -nu * log(sigmau) - sum(yu) + nu * log(phiu)
  } else {
    loglik <- -nu * log(sigmau) - (1 / xi + 1) * sum(log(syu)) +
      nu * log(phiu)
  }

  nllh <- -loglik
  if (finitelik && is.infinite(nllh)) {
    nllh <- sign(nllh) * 1e6
  }
  nllh
}

# Maximum-Likelihood Fit of the GPD Above a Threshold
#' @importFrom stats var optim
.fgpd_internal <- function(x, u = 0, phiu = NULL, pvector = NULL,
                           method = "BFGS", control = list(maxit = 10000),
                           finitelik = TRUE, ...) {
  x <- x[!is.infinite(x)]
  n <- length(x)
  if (method == "L-BFGS-B" || method == "BFGS") {
    finitelik <- TRUE
  }

  # Missing values count toward n, so they are treated as below threshold
  xu <- x[which(x > u)]
  nu <- length(xu)
  if (nu < 1) {
    stop("no elements of x are above threshold")
  }

  if (is.null(pvector)) {
    # Gumbel moment match for the scale, sigmau = sqrt(6*var(y))/pi, and a
    # fixed mildly heavy-tailed start for the shape
    pvector <- c(sqrt(6 * var(xu - u)) / pi, 0.1)
  }

  if (is.null(phiu)) {
    phiu <- nu / n
  } else if (phiu == 0) {
    stop("tail probability must be in (0, 1]")
  }

  nllh <- .nlgpd_internal(pvector, xu, u, phiu)
  if (is.infinite(nllh)) {
    pvector[2] <- 0.1
    nllh <- .nlgpd_internal(pvector, xu, u, phiu)
  }
  if (is.infinite(nllh)) {
    stop("initial parameter values are invalid")
  }

  fit <- optim(par = as.vector(pvector), fn = .nlgpd_internal, x = xu, u = u,
               phiu = phiu, finitelik = finitelik, control = control,
               method = method, hessian = TRUE, ...)

  conv <- TRUE
  if (fit$convergence != 0 || any(fit$par == pvector) ||
      abs(fit$value) >= 1e6) {
    conv <- FALSE
    warning("check convergence")
  }

  list(x = as.vector(x), init = as.vector(pvector), optim = fit, conv = conv,
       mle = fit$par, rate = phiu, nllh = fit$value, n = n, u = u,
       sigmau = fit$par[1], xi = fit$par[2], phiu = phiu)
}

#################################################################
# Section 1: Main Fitting Engine
#################################################################
# Fit a Spliced Gamma and GPD Model
#' @importFrom methods is
.fdiscgammagpd <- function(x, useq, shift = NULL, pvector=NULL,
                          std.err = TRUE, method = "Nelder-Mead", ...){
  if(!is(x, "numeric")){
    stop("x must be numeric.")
  }
  
  if(!is.null(shift)){
    if(!is(shift, "numeric")){
      stop("shift must be numeric.")
    }
    if(shift != round(shift)){
      stop("shift must be an integer.")
    }
  }
  
  if(!is(useq, "numeric")){
    stop("useq must be numeric.")
  }
  
  if(any(x != round(x))){
    stop("all elements in x must be integers.")
  }
  
  if(any(useq != round(useq))){
    stop("all elements in useq must be integers.")
  }
  
  if(!is.null(pvector) & !(length(pvector) == 5)){
    stop("pvector must contain 5 elements.")
  }
  
  if(!(is.logical(std.err))){
    stop("std.err must be TRUE or FALSE.")
  }
  
  if(is.null(shift)){
    shift <- min(x)
  }
  
  if (is.null(pvector)) {
    pvector <- rep(NA,5)
    s <- log(mean(x+0.5))-mean(log(x+0.5))
    k <- (3-s + sqrt((s-3)^2 + 24*s))/12/s
    pvector[1] <- k
    pvector[2] <- k/mean(x)
    pvector[3] <- as.vector(quantile(x, 0.9))
    
    xu <- x[x>=pvector[3]]
    initfgpd <- .fgpd_internal(xu, min(xu)-10^(-5))
    pvector[4] <- initfgpd$mle[1]
    pvector[5] <- initfgpd$mle[2]
  }
  
  bulk <- lapply(seq_along(useq),
                 function(idx,x,useq) x < useq[idx], x=x, useq=useq)
  tail <- lapply(seq_along(useq),
                 function(idx,x,useq) x >= useq[idx], x=x, useq=useq)
  phiu <- lapply(seq_along(useq),
                 function(idx,tail) mean(tail[[idx]]), tail=tail)
  
  gammfit <- list()
  gpdfit <- list()
  nllhu <- rep(NA, length(useq))
  for(i in seq_along(useq)){
    gammfit[[i]] <- tryCatch(expr = .fdiscgamma(pvector[1:2],x[bulk[[i]]],
                                               useq[i],
                                               phiu[[i]],
                                               shift,
                                               method="Nelder-Mead"),
                             error = function(err) NA)
    gpdfit[[i]] <- tryCatch(expr = .fdiscgpd(pvector[4:5],
                                            x[tail[[i]]],
                                            useq[i],
                                            phiu[[i]],
                                            method="Nelder-Mead"),
                            error = function(err) {
                              tryCatch(expr = {
                                pvec3 <- as.vector(quantile(x,1-phiu[[i]]))
                                xu <- x[x>=pvec3]
                                initfgpd.adj <-
                                  .fgpd_internal(x, min(xu)-10^(-5))
                                pvec4 <- initfgpd.adj$mle[1]
                                pvec5 <- initfgpd.adj$mle[2]
                                .fdiscgpd(c(pvec4,pvec5),
                                          x[tail[[i]]],
                                          useq[i],
                                          phiu[[i]],
                                          method="Nelder-Mead")
                              }, error = function(err2) NA)
                            })
    nllhu[i] <- tryCatch(expr = gammfit[[i]]$value + gpdfit[[i]]$value,
                         error = function(err) NA)
  }
  
  bestfit <- which.min(nllhu)
  if(length(bestfit) == 0){
    stop("the spliced gamma-GPD model could not be fit at any threshold in ",
         "`useq`. This happens when the counts above the threshold carry too ",
         "little information to identify a tail, for instance when nearly ",
         "every clone above it is the same size. Supplying a different ",
         "`threshold` may help.", call. = FALSE)
  }
  fit.out <- list(gammfit[[bestfit]], gpdfit[[bestfit]])
  names(fit.out) <- c("bulk","tail")
  mle <- c(mean(x >= useq[bestfit]),
           exp(fit.out$bulk$par),
           useq[bestfit],
           exp(fit.out$tail$par[1]),
           fit.out$tail$par[2])
  names(mle) <- c("phi","shape","rate","thresh","sigma","xi")
  if(std.err){
    H <- fit.out$bulk$hessian %>% rbind(matrix(rep(0,4),nrow = 2)) %>%
      cbind(rbind(matrix(rep(0,4),nrow = 2),fit.out$tail$hessian))
    fisherInf <- tryCatch(expr = solve(H), error = function(err) NA)
    out <- list(x = as.vector(x), shift = shift, init = as.vector(pvector),
                useq = useq, nllhuseq = nllhu,
                optim = fit.out, nllh = nllhu[bestfit],
                mle=mle, fisherInformation = fisherInf)
  } else{
    out <- list(x = as.vector(x), shift = shift, init = as.vector(pvector),
                useq = useq, nllhuseq = nllhu,
                optim = fit.out, nllh = nllhu[bestfit],
                mle=mle)
  }
  out
}

#################################################################
# Section 2: Core Modeling Functions
#################################################################

# Discrete Gamma Probability Mass Function
#' @importFrom methods is
.ddiscgamma <- function(x, shape, rate, thresh, phiu, shift = 0, log = FALSE){
  if(any(x != floor(x))){
    stop("x must be an integer")
  }
  
  out <- rep(0, length(x))
  
  up <- pgamma(x+1-shift, shape=shape, rate=rate)
  down <- pgamma(x-shift, shape=shape, rate=rate)
  
  if(!log){
    b <- pgamma(thresh-shift, shape=shape, rate=rate)
    out[x < thresh] <- ((1-phiu)*(up-down)/b)[x < thresh]
  } else{
    b <- pgamma(thresh-shift, shape=shape, rate=rate, log.p = TRUE)
    out[x < thresh] <- (log(1-phiu)+log(up-down) - b)[x < thresh]
  }
  out
}

.pdiscgamma <- function(q, shape, rate, thresh, phiu, shift = 0){
  probs <- .ddiscgamma(0:q, shape, rate, thresh, phiu, shift)
  sum(probs)
}

.qtruncgamma_internal <- function(p, shape, rate, a=0, b=Inf) {
  Fa <- pgamma(a, shape, rate)
  Fb <- pgamma(b, shape, rate)
  p_adjusted <- p * (Fb - Fa) + Fa
  qgamma(p_adjusted, shape, rate)
}

.qdiscgamma <- function(p, shape, rate, thresh, phiu, shift = 0){
  .qtruncgamma_internal(p/(1-phiu), 
                        shape=shape, 
                        rate=rate, 
                        b=thresh-shift) %>% 
    floor() + shift
}

.rtruncgamma_internal <- function(n, shape, rate, a = 0, b = Inf) {
  Fa <- pgamma(a, shape, rate)
  Fb <- pgamma(b, shape, rate)
  u <- runif(n, min = Fa, max = Fb)
  qgamma(u, shape, rate)
}

.rdiscgamma <- function(n, shape, rate, thresh, shift = 0){
  .rtruncgamma_internal(n, 
                        shape = shape, 
                        rate = rate, 
                        b = thresh - shift) %>% 
    floor() + shift
}

# Discrete GPD Probability Mass Function
.ddiscgpd <- function(x, thresh, sigma, xi, phiu, log = FALSE){
  up <- .pgpd_internal(x+1, u=thresh, sigmau=sigma, xi=xi)
  down <- .pgpd_internal(x, u=thresh, sigmau=sigma, xi=xi)
  
  if(!log){
    phiu*(up-down)
  } else{
    log(phiu) + log(up-down)
  }
}

.pdiscgpd <- function(q, thresh, sigma, xi, phiu){
  probs <- .ddiscgpd(thresh:q, thresh, sigma, xi, phiu)
  sum(probs)
}

.qdiscgpd <- function(p, thresh, sigma, xi, phiu){
  .qgpd_internal(p/phiu, u=thresh, sigmau = sigma, xi=xi) %>% floor
}

.rdiscgpd <- function(n, thresh, sigma, xi){
  .rgpd_internal(n, u=thresh, sigmau=sigma, xi=xi) %>% floor
}

# Spliced Gamma-GPD Probability Mass Function
.ddiscgammagpd <- function(x, fit=NULL, shape, rate, u, sigma, xi,
                           phiu=NULL, shift = 0, log = FALSE){
  if(!is.null(fit)){
    if(!all(c("x", "shift", "init", "useq", "nllhuseq", "nllh",
              "optim", "mle") %in% names(fit))){
      stop("\"fit\" is not of the correct structure. It must be one model
                 fit from fdiscgammagpd.")
    }
    phiu <- fit$mle['phi']
    shape <- fit$mle['shape']
    rate <- fit$mle['rate']
    u <- fit$mle['thresh']
    sigma <- fit$mle['sigma']
    xi <- fit$mle['xi']
    shift <- fit$shift
  }
  if(!is(x, "numeric")){
    stop("x must be numeric.")
  }
  
  if(!is(shift, "numeric")){
    stop("shift must be numeric.")
  }
  
  if(any(x != floor(x))){
    stop("x must be an integer")
  }
  
  if(shift != round(shift)){
    stop("shift must be an integer.")
  }
  
  if(any(c(shape, rate, sigma) <= 0)){
    stop("shape, rate, and sigma must all be positive.")
  }
  
  if(!is.null(phiu)){
    if(phiu < 0 | phiu > 1){
      stop("phiu must be in [0,1].")
    }
  }
  
  if(is.null(phiu)){
    phiu <- 1-.pdiscgamma(u-1, shape=shape, rate=rate,
                          thresh=Inf, phiu = 0, shift=shift)
  }
  
  out <- rep(NA, length(x))
  
  if(sum(x>=u) != 0){
    out[x>=u] <- .ddiscgpd(x[x>=u], u, sigma, xi, phiu, log=log)
  }
  
  if(sum(x<u) != 0){
    out[x<u] <- .ddiscgamma(x[x<u], shape, rate, u, phiu, shift, log=log)
  }
  out
}

# --- Negative Log-Likelihood (NLL) Functions ---

# NLL for Discrete Gamma
.discgammanll <- function(param, dat, thresh, phiu, shift=0){
  shape <- exp(param[1])
  rate <- exp(param[2])
  
  # As in .discgpdnll, keep an unconstrained search inside the valid region
  if(!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0){
    return(Inf)
  }
  
  if(any(dat > thresh-1)){ warning("data must be less than the threshold") }
  
  ll <- log(.ddiscgamma(dat, shape, rate, thresh, phiu, shift))
  
  sum(-ll)
}

# NLL for Discrete GPD
.discgpdnll <- function(param, dat, thresh, phiu){
  sigma <- exp(param[1])
  xi <- param[2]
  
  # Nelder-Mead searches log(sigma) unconstrained, so it can underflow sigma
  # to zero or overflow it to infinity. The likelihood is undefined there, so
  # report an infinite penalty the optimizer can step away from instead of
  # raising an error that would discard the whole fit.
  if(!is.finite(sigma) || sigma <= 0 || !is.finite(xi)){
    return(Inf)
  }
  
  ll <- log(.ddiscgpd(dat, thresh, sigma, xi, phiu))
  
  sum(-ll)
}

# --- Fitting Wrappers ---

# Wrapper for Gamma Optimization
.fdiscgamma <- function(param, dat, thresh, phiu, shift = 0, method, ...){
  opt <- optim(log(param), .discgammanll, dat=dat, thresh=thresh,
               phiu=phiu, shift=shift, method=method, hessian = TRUE, ...)
  opt
}

# Wrapper for GPD Optimization
.fdiscgpd <- function(param, dat, thresh, phiu, method, ...){
  opt <- optim(c(log(param[1]),param[2]), .discgpdnll, dat=dat, thresh=thresh,
               phiu=phiu, method=method, hessian = TRUE, ...)
  opt
}

#################################################################
# Section 3: Distance Calculation
#################################################################

# Calculate Pairwise Distances
#' @importFrom methods is
.get_distances <- function(fits, grid, modelType = "Spliced"){
  if(!is(grid, "numeric")){
    stop("grid must be numeric.")
  }
  if(any(grid != round(grid))){
    stop("all elements in grid must be integers.")
  }
  if(!(modelType %in% c("Spliced", "Desponds"))){
    stop("modelType must be either \"Spliced\" or \"Desponds\".")
  }
  
  distances <- matrix(rep(0, length(fits)^2), nrow = length(fits))
  if(!is.null(names(fits))){
    rownames(distances) <- colnames(distances) <- names(fits)
  }
  
  for(i in seq_len((length(fits)-1))){
    for(j in (i+1):length(fits)){
      distances[i,j] <- .JS_dist(fits[[i]],
                                 fits[[j]],
                                 grid,
                                 modelType = modelType)
    }
  }
  distances <- distances + t(distances)
  distances
}

.JS_dist <- function(fit1, fit2, grid, modelType = "Spliced"){
  if(!(modelType %in% c("Spliced", "Desponds"))){
    stop("modelType must be either \"Spliced\" or \"Desponds\".")
  }
  
  if(modelType == "Spliced"){
    if(!all(c("x", "init", "useq", "nllhuseq", "nllh",
              "optim", "mle") %in% names(fit1))){
      stop("\"fit1\" is not of the correct structure. It must be a model
                 fit from fdiscgammagpd.")
    }
    if(!all(c("x", "init", "useq", "nllhuseq", "nllh",
              "optim", "mle") %in% names(fit2))){
      stop("\"fit2\" is not of the correct structure. It must be a model
                 fit from fdiscgammagpd.")
    }
    shiftp <- fit1$shift
    shiftq <- fit2$shift
    phip <- fit1$mle['phi']
    phiq <- fit2$mle['phi']
    shapep <- fit1$mle['shape']
    shapeq <- fit2$mle['shape']
    ratep <- fit1$mle['rate']
    rateq <- fit2$mle['rate']
    threshp <- fit1$mle['thresh']
    threshq <- fit2$mle['thresh']
    sigmap <- fit1$mle['sigma']
    sigmaq <- fit2$mle['sigma']
    xip <- fit1$mle['xi']
    xiq <- fit2$mle['xi']
    
    out <- .JS_spliced(grid, shiftp, shiftq, phip, phiq, shapep, shapeq,
                       ratep, rateq, threshp, threshq, sigmap, sigmaq,
                       xip, xiq)
  } else if(modelType == "Desponds"){
    if(!all(c("min.KS", "Cmin", "powerlaw.exponent",
              "pareto.alpha") == names(fit1))){
      stop("\"fit1\" is not of the correct structure. It must be a model
                 fit from fdesponds.")
    }
    if(!all(c("min.KS", "Cmin", "powerlaw.exponent",
              "pareto.alpha") == names(fit2))){
      stop("\"fit2\" is not of the correct structure. It must be a model
                 fit from fdesponds.")
    }
    Cminp <- fit1['Cmin']
    Cminq <- fit2['Cmin']
    alphap <- fit1['pareto.alpha']
    alphaq <- fit2['pareto.alpha']
    out <- .JS_desponds(grid, Cminp, Cminq, alphap, alphaq)
  }
  out
}

# Jensen-Shannon Distance Calculation
#' @importFrom methods is
.JS_spliced <- function(grid, shiftp, shiftq, phip, phiq, shapep, shapeq, ratep,
                        rateq, threshp, threshq, sigmap, sigmaq, xip, xiq){
  if(!is(grid, "numeric")){
    stop("grid must be numeric.")
  }
  
  if(any(grid != round(grid))){
    stop("all elements in grid must be integers.")
  }
  
  if(any(!is(c(shiftp, shiftq, phip, phiq, shapep, shapeq,
               ratep, rateq, threshp, threshq,
               sigmap, sigmaq, xip, xiq), "numeric"))){
    stop("shiftp, shiftq, phip, phiq, shapep, shapeq, ratep, rateq,
              threshp, threshq, sigmap, sigmaq, xip, and xiq must be numeric.")
  }
  
  if(shiftp != round(shiftp) | shiftq != round(shiftq)){
    stop("shiftp and shiftq must be integers.")
  }
  
  if(any(c(shapep, shapeq, ratep, rateq, sigmap, sigmaq) <= 0)){
    stop("shapep, shapeq, ratep, rateq, sigmap, and sigmaq must be 
             greater than 0.")
  }
  
  if(any(c(phip, phiq) > 1) | any(c(phip, phiq) < 0)){
    stop("phip and phiq must be in [0,1].")
  }
  
  if(ratep <= 0 | rateq <= 0){
    stop("ratep and rateq must be greater than 0.")
  }
  
  if(shapep <= 0 | shapeq <= 0){
    stop("shapep and shapeq must be greater than 0.")
  }
  
  if(threshp != round(threshp) | threshq != round(threshq)){
    stop("threshp and threshq must be integers.")
  }
  
  K <- max(grid)
  
  P <- .ddiscgammagpd(min(grid):K, shape = shapep, rate = ratep,
                      u=threshp, sigma = sigmap,
                      xi = xip, phiu = phip, shift=shiftp,
                      log = FALSE)
  adjp <- which(P == 0)
  if(length(adjp) != 0){
    P[adjp] <- .dgpd_internal(adjp+0.5, u=threshp,
                              sigmau = sigmap, xi = xip, phiu = phip)
  }
  
  Q <- .ddiscgammagpd(min(grid):K, shape = shapeq, rate = rateq,
                      u=threshq, sigma = sigmaq,
                      xi = xiq, phiu = phiq, shift=shiftq,
                      log = FALSE)
  adjq <- which(Q == 0)
  if(length(adjq) != 0){
    Q[adjq] <- .dgpd_internal(adjq+0.5, u=threshq,
                              sigmau = sigmaq, xi = xiq, phiu = phiq)
  }
  
  M <- 0.5*(P+Q)
  pzero <- which(P == 0)
  qzero <- which(Q == 0)
  
  sum1 <- sum2 <- rep(NA, length(grid))
  sum1 <- P*(log(P) - log(M))
  sum2 <- Q*(log(Q) - log(M))
  
  if(length(intersect(pzero, qzero)) != 0){
    sum1[intersect(pzero, qzero)] <- 0
    sum2[intersect(pzero, qzero)] <- 0
  }
  if(length(setdiff(pzero, qzero)) != 0){
    sum1[setdiff(pzero, qzero)] <- 0
  }
  if(length(setdiff(qzero, pzero)) != 0){
    sum2[setdiff(qzero, pzero)] <- 0
  }
  
  out <- sqrt(0.5*(sum(sum1) + sum(sum2)))
  out
}

# Core JS Distance Logic
#' @importFrom stats integrate
#' @importFrom methods is
.JS_desponds <- function(grid, Cminp, Cminq, alphap, alphaq){
  if(!is(grid, "numeric")){
    stop("grid must be numeric.")
  }
  
  if(any(grid != round(grid))){
    stop("all elements in grid must be integers.")
  }
  
  if(any(!is(c(Cminp, Cminq, alphap, alphaq), "numeric"))){
    stop("Cminp, Cminq, alphap, and alphaq must be numeric.")
  }
  
  if(Cminp != round(Cminp) | Cminq != round(Cminq)){
    stop("Cminp and Cminq must be integers.")
  }
  
  if(alphap <= 0 | alphaq <= 0){
    stop("alphap and alphaq must be greater than 0.")
  }
  
  lower <- min(grid)
  upper <- max(grid)
  
  out <- integrate(.eval_desponds,
                   lower = lower, upper = upper,
                   Cminp = Cminp, Cminq = Cminq,
                   alphap = alphap, alphaq = alphaq)$value %>% sqrt
  out
}

.dpareto_internal <- function(x, scale, shape) {
  ifelse(x < scale, 0, (shape * (scale^shape)) / (x^(shape + 1)))
}

.eval_desponds <- function(t, Cminp, Cminq, alphap, alphaq){
  M <- 0.5*(.dpareto_internal(t, scale=Cminp, shape=alphap) +
              .dpareto_internal(t, scale=Cminq, shape=alphaq))
  
  one <- .dpareto_internal(t, scale=Cminp, shape=alphap)
  two <- .dpareto_internal(t, scale=Cminq, shape=alphaq)
  
  if(one == 0 & two == 0){
    out <- 0
  } else if(one == 0 & two != 0){
    out <- .dpareto_internal(t, scale=Cminq, shape=alphaq) *
      (log(.dpareto_internal(t, scale=Cminq, shape=alphaq))-log(M))
  } else if(one != 0 & two == 0){
    out <- .dpareto_internal(t, scale=Cminp, shape=alphap) *
      (log(.dpareto_internal(t, scale=Cminp, shape=alphap))-log(M))
  } else{
    out <- .dpareto_internal(t, scale=Cminp, shape=alphap) *
      (log(.dpareto_internal(t, scale=Cminp, shape=alphap))-log(M)) +
      .dpareto_internal(t, scale=Cminq, shape=alphaq) *
      (log(.dpareto_internal(t, scale=Cminq, shape=alphaq))-log(M))
  }
  out
}
