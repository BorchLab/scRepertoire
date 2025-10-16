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
#'                        cloneCall = "strict", 
#'                        method="ward.D2")
#'
#' @param input.data The product of [combineTCR()],
#' [combineBCR()], or [combineExpression()].
#' @param cloneCall Defines the clonal sequence grouping. Accepted values 
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino 
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains 
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL` 
#' (for both light chains), `both`.
#' @param threshold Numerical vector containing the thresholds
#' the grid search was performed over.
#' @param method The clustering parameter for the dendrogram.
#' @param group.by A column header in the metadata or lists to group the analysis 
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed as 
#' by list element or active identity in the case of single-cell objects.
#' @param exportTable If `TRUE`, returns a data frame or matrix of the results 
#' instead of a plot.
#' @param palette Colors to use in visualization - input any
#' [hcl.pals][grDevices::hcl.pals].
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
                                   cloneCall ="strict", 
                                   chain = "both", 
                                   method = "ward.D2", 
                                   threshold = 1, 
                                   group.by = NULL,
                                   exportTable = FALSE, 
                                   palette = "inferno",
                                   ...) {
  x <- y <- xend <- yend <- mpg_div_hp <- NULL
  input.data <- .dataWrangle(input.data, 
                             group.by, 
                             .theCall(input.data, cloneCall, 
                                      check.df = FALSE, silent = TRUE), 
                             chain)
  cloneCall <- .theCall(input.data, cloneCall)
  sco <- .is.seurat.or.se.object(input.data)
  if(!is.null(group.by) & !sco) {
    input.data <- .groupList(input.data, group.by)
  }
  data <- bind_rows(input.data)
  unique_df<- unique(data[,cloneCall])
  
  # Create long-format summary table
  summary_df <- dplyr::bind_rows(input.data, .id = "sample") %>%
    dplyr::group_by(sample, .data[[cloneCall]]) %>%
    dplyr::summarise(Freq = dplyr::n(), .groups = 'drop')
  wide_matrix <- xtabs(as.formula(paste("Freq ~", cloneCall, "+ sample")), 
                       data = summary_df)
  Con.df <- as.data.frame.matrix(wide_matrix)
  Con.df[[cloneCall]] <- rownames(Con.df)
  rownames(Con.df) <- NULL
  Con.df <- Con.df[, c(cloneCall, setdiff(names(Con.df), cloneCall))]
  
  # Fit models
  list <- lapply(seq_len(ncol(Con.df))[-1], function(x) {
    suppressWarnings(.fdiscgammagpd(Con.df[,x], useq = threshold))
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
  
  if (exportTable) { 
    return(distances) 
  }
  return(plot)
}

#################################################################
# Section 1: Main Fitting Engine
#################################################################
# Fit a Spliced Gamma and GPD Model
#' @importFrom evmix fgpd
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
    initfgpd <- fgpd(xu, min(xu)-10^(-5))
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
                              pvec3 <- as.vector(quantile(x,1-phiu[[i]]))
                              xu <- x[x>=pvec3]
                              initfgpd.adj <-
                                fgpd(x, min(xu)-10^(-5))
                              pvec4 <- initfgpd.adj$mle[1]
                              pvec5 <- initfgpd.adj$mle[2]
                              tryCatch(expr = .fdiscgpd(c(pvec4,pvec5),
                                                       x[tail[[i]]],
                                                       useq[i],
                                                       phiu[[i]],
                                                       method="Nelder-Mead"),
                                       error = function(err2) NA)
                            })
    nllhu[i] <- tryCatch(expr = gammfit[[i]]$value + gpdfit[[i]]$value,
                         error = function(err) NA)
  }
  
  bestfit <- which.min(nllhu)
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
#' @importFrom evmix pgpd
.ddiscgpd <- function(x, thresh, sigma, xi, phiu, log = FALSE){
  up <- pgpd(x+1, u=thresh, sigmau=sigma, xi=xi)
  down <- pgpd(x, u=thresh, sigmau=sigma, xi=xi)
  
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

#' @importFrom evmix qgpd
.qdiscgpd <- function(p, thresh, sigma, xi, phiu){
  qgpd(p/phiu, u=thresh, sigmau = sigma, xi=xi) %>% floor
}

#' @importFrom evmix rgpd
.rdiscgpd <- function(n, thresh, sigma, xi){
  rgpd(n, u=thresh, sigmau=sigma, xi=xi) %>% floor
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
  
  if(any(dat > thresh-1)){ warning("data must be less than the threshold") }
  
  ll <- log(.ddiscgamma(dat, shape, rate, thresh, phiu, shift))
  
  sum(-ll)
}

# NLL for Discrete GPD
.discgpdnll <- function(param, dat, thresh, phiu){
  sigma <- exp(param[1])
  xi <- param[2]
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
#' @importFrom evmix dgpd
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
    P[adjp] <- dgpd(adjp+0.5, u=threshp,
                    sigmau = sigmap, xi = xip, phiu = phip)
  }
  
  Q <- .ddiscgammagpd(min(grid):K, shape = shapeq, rate = rateq,
                      u=threshq, sigma = sigmaq,
                      xi = xiq, phiu = phiq, shift=shiftq,
                      log = FALSE)
  adjq <- which(Q == 0)
  if(length(adjq) != 0){
    Q[adjq] <- dgpd(adjq+0.5, u=threshq,
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
