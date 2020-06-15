#' @import graphics stats utils

to.Fx <- function(x) {
  # A faster version of ecdf(x)(x)
  rank(x, ties.method="max")/length(x)
}

#' Constructs a simultaneous credible band 
#' 
#' \code{sim.cred.band} returns a simultaneous band over a finite set of
#' covariate points given either a sample from the posterior of the
#' regression surface or a function \code{FUN(x, params)} and a sample from
#' the posterior of the parameters.
#' 
#' If \code{design} is \code{NULL} (default), it is taken to be the identity
#' matrix of dimension \code{ncol(params)}, so that the rows of params
#' are treated as draws from the posterior \code{FUN(x, params)}.
#'
#' The 'asymptotic' method assumes that the marginal posteriors of 
#' the \code{FUN(x, params)} are asymptotically normal and is usually
#' significantly faster and less memory-intensive than the 'quantile'
#' method, which makes no such assumption.
#' 
#' @param params A numeric matrix whose rows are draws from the posterior
#'  distribution of either the regression surface or the parameter vector.
#' @param design (Optional) A numeric matrix whose rows are covariate points
#' over which the band is to be constructed.
#' @param FUN (Optional) a function of the form \code{function(x, params)}
#' that takes a row of \code{design} and the entire \code{params} matrix and
#' returns a vector of the same length of \code{x} representing the regression
#' surface.
#' @param cred.level Numeric; the credible level.
#' @param method Either "asymptotic" (default) or "quantile"; see details.
#' @param sides One of "both" (default), "upper", or "lower".
#'              Which bounds should be constructed?
#' @param est.FUN The function used to produce estimates of the regression
#' surface. Default is \code{mean}.
#' @param var.FUN The function used to quantify the variability of the
#' regression surface posterior. Default is \code{sd}.
#' @param point.estimate If not null, replaces the mean and sets the reference 
#'                       around which the standard error is computed.
#'                       Useful for bootstrapping methods.
#'                       Treated as a row of the \code{params} matrix.
#' @param track A numeric vector of indices indicating which rows (default none)
#' of the design matrix should have the sample of the corresponding
#' \code{FUN(x, params)} returned.
#' @param verbose Logical (default \code{FALSE}); print progress?
#' 
#' @return An object of class \code{sim.cred.band}, which contains:
#' \describe{
#'   \item{\code{upper}}{A numeric vector of upper bounds.}
#'   \item{\code{lower}}{A numeric vector of lower bounds.}
#'   \item{\code{cred.level}}{As provided.}
#'   \item{\code{method}}{As provided.}
#'   \item{\code{sides}}{As provided.}
#'   \item{\code{est}}{Posterior estimate of the regression surface.}
#'   \item{\code{est.FUN}}{As provided.}
#'   \item{\code{var}}{Summary of posterior variability of the regression
#'                     surface.}
#'   \item{\code{var.FUN}}{As provided.}
#'   \item{\code{W}}{An estimate of the extremal errors.}
#'   \item{\code{W.crit}}{The critical quantile of W.}
#'   \item{\code{trace}}{The posterior samples of the regression surface
#'                       indicated by the \code{track} argument.}
#' }
#' 
#' @examples
#' ### Sample from regression surface posterior
#' reg.surf.sample <- matrix(rnorm(1000, mean=1:10), ncol=2, byrow=TRUE)
#' sim.cred.band(reg.surf.sample, cred.level=0.80)
#' 
#' ### Parametric case
#' design <- cbind(1, 1:10)
#' params <- matrix(rnorm(200, mean=1:2), ncol=2, byrow=TRUE)
#' sim.cred.band(params, design)
#' 
#' ### With custom function
#' params.sd <- cbind(1 / rgamma(100, 1), params)
#' FUN.sd <- function(x, params) { params[, -1] %*% t(x) / params[, 1] }
#' sim.cred.band(params.sd, design, FUN.sd)
#' 
#' @export
sim.cred.band <- function(params,
                          design=NULL,
                          FUN=function(x, params) { params %*% t(x) },
                          cred.level=0.95,
                          method=c('asymptotic', 'quantile'),
                          sides=c('both', 'upper', 'lower'),
                          est.FUN=mean,
                          var.FUN=sd,
                          point.estimate=NULL,
                          track=numeric(0),
                          verbose=FALSE) {
  
  # Validate and shape input
  params <- data.matrix(params)
  M <- nrow(params)
  
  if (is.null(design)) {
    N <- ncol(params)
    nonpar <- TRUE
  } else {
    design <- data.matrix(design)
    N <- nrow(design)
    nonpar <- FALSE
  }
  
  est <- var <- numeric(N)
  
  method <- method[1]
  stopifnot(method %in% c('asymptotic', 'quantile'))
  if (method == 'asymptotic') {
    m <- numeric(N)
    s <- numeric(N)
  }
  
  sides <- sides[1]
  stopifnot(sides %in% c('both', 'upper', 'lower'))
  
  sim.cred.band <- list(upper=rep(Inf, N),
                        lower=rep(-Inf, N),
                        cred.level=cred.level,
                        method=method,
                        sides=sides,
                        W.crit=NA,
                        W=NA,
                        est=rep(NA, N),
                        est.FUN=est.FUN,
                        var=rep(NA, N),
                        var.FUN=var.FUN,
                        trace=matrix(NA, nrow=M, ncol=length(track)))
  colnames(sim.cred.band$trace) <- track
  class(sim.cred.band) <- 'sim.cred.band'
  
  # Compute W
  W <- rep(-Inf, M)
  
  # This is iterative to avoid storing an NxM matrix (or two)
  for (i in 1:N) {
    
    if (verbose && (i %% 100 == 0)) {
      cat(i, "/", N, "\n")
    }
    
    if (nonpar) {
      fx <- params[, i]
    } else {
      fx <- FUN(design[i, , drop=FALSE], params)
    }
    
    if (i %in% track) {
      sim.cred.band$trace[, which(track == i)] <- fx
    }
    
    est[i] <- est.FUN(fx)
    var[i] <- var.FUN(fx)
    
    if (method == 'asymptotic') {
      if (is.null(point.estimate)) {
        m[i] <- mean(fx)
        s[i] <- sd(fx)
      } else {
        if (nonpar) {
          point.fx <- point.estimate[i]
        } else {
          point.fx <- FUN(design[i, , drop=FALSE], point.estimate)
        }
        m[i] <- point.fx
        s[i] <- sqrt(mean((fx - point.fx) ^ 2))
      }
      
      if (sides == "both") {
        z <- abs(fx - m[i]) / s[i]
      } else if (sides == "upper") {
        z <- (fx - m[i]) / s[i]
      } else if (sides == "lower") {
        z <- (m[i] - fx) / s[i]
      }
      
    } else {
      if (sides == "both") {
        Fx <- to.Fx(fx)
        Gx <- 1 - to.Fx(-fx)
        z <- pmax(1-Fx, Gx)
      } else if (sides == "upper") {
        z <- 1 - to.Fx(-fx)
      } else if (sides == "lower") {
        z <- 1 - to.Fx(fx)
      }
    }
    W <- pmax(W, z)
  }
  
  sim.cred.band$W <- W
  
  # Estimate W.crit
  W.crit <- quantile(W, cred.level, type=1)
  sim.cred.band$W.crit <- W.crit
  
  # Compute bounds
  if (method == 'asymptotic') {
    if (sides == "both" | sides == "upper") {
      sim.cred.band$upper <- m + W.crit * s
    }
    
    if (sides == "both" | sides == "lower") {
      sim.cred.band$lower <- m - W.crit * s
    }
    
  } else {
    if (nonpar) {
      if (sides %in% c('both', 'upper')) {
        sim.cred.band$upper <- apply(params, 2, function(fx) {
          -quantile(-fx, prob=1-W.crit, type=1)
        })
      }
      if (sides %in% c('both', 'lower')) {
        sim.cred.band$lower <- apply(params, 2, function(fx) {
          quantile(fx, prob=1-W.crit, type=1)
        })
      }
    
    } else {
      if (sides == 'both') {
        bounds <- apply(design, 1, function(x, params) {
          fx <- FUN(t(x), params)
          upper <- -quantile(-fx, prob=1-W.crit, type=1)
          lower <-  quantile( fx, prob=1-W.crit, type=1)
          c(lower, upper)
        }, params=params)
        sim.cred.band$lower <- bounds[1, ]
        sim.cred.band$upper <- bounds[2, ]
      } else if (sides == 'upper') {
        sim.cred.band$upper <- apply(design, 1, function(x, params) {
          -quantile(-fx, prob=1-W.crit, type=1)
        }, params=params)
      } else if (sides == 'lower') {
        sim.cred.band$lower <- apply(design, 1, function(x, params) {
           quantile( fx, prob=1-W.crit, type=1)
        }, params=params)
      }
    }
  }
  
  sim.cred.band$est <- est
  sim.cred.band$var <- var

  sim.cred.band
}

#' Constructs a credible subset pair
#' 
#' \code{credsubs} returns a credible subset pair over a finite set of
#' covariate points given either a sample from the posterior of the
#' regression surface or a function \code{FUN(x, params)} and a sample from
#' the posterior of the parameters.
#' 
#' If design is NULL (default), it is taken to be the identity
#' matrix of dimension ncol(params), so that the rows of params
#' are treated as draws from the posterior FUN(x, params).
#'
#' The 'asymptotic' method assumes that the marginal posteriors of 
#' the FUN(x, params) are asymptotically normal and is usually
#' significantly faster and less memory-intensive than the 'quantile'
#' method, which makes no such assumption.
#' 
#' @param params A numeric matrix whose rows are draws from the posterior
#'               distribution of either the regression surface or the
#'               parameter vector.
#' @param design (Optional) A numeric matrix whose rows are covariate points
#'               over which the band is to be constructed.
#' @param FUN (Optional) a function of the form \code{function(x, params)}
#'            that takes a row of \code{design} and the entire \code{params}
#'            matrix and returns a vector of the same length of \code{x}
#'            representing the regression surface.
#' @param cred.level Numeric; the credible level.
#' @param threshold Numeric; the value of \code{FUN} above which
#'                  a covariate is included in the target subset.
#' @param method Either "asymptotic" (default) or "quantile"; see details.
#' @param step.down Logical (default \code{TRUE}); should the step-down
#'                  procedure be used?
#' @param sides One of "both" (default), "exclusive", or "inclusive".
#'              Which bounds should be constructed?
#' @param est.FUN The function used to produce estimates of the regression
#'                surface. Default is \code{mean}.
#' @param var.FUN The function used to quantify the variability of the
#'                regression surface posterior. Default is \code{sd}.
#' @param point.estimate If not null, replaces the mean and sets the reference 
#'                       around which the standard error is computed.
#'                       Useful for bootstrapping methods.
#'                       Treated as a row of the \code{params} matrix.
#' @param track A numeric vector of indices indicating which rows (default none)
#'              of the design matrix should have the sample of the corresponding
#'              \code{FUN(x, params)} returned.
#' @param verbose Logical (default \code{FALSE}); print progress?
#' 
#' @return An object of class \code{credsubs}, which contains:
#' \describe{
#'   \item{\code{exclusive}}{A logical vector indicating membership in
#'                           the exclusive credible subset.}
#'   \item{\code{inclusive}}{A logical vector indicating membership in
#'                           the inclusive credible subset.}
#'   \item{\code{cred.level}}{As provided.}
#'   \item{\code{threshold}}{As provided.}
#'   \item{\code{method}}{As provided.}
#'   \item{\code{step.down}}{As provided.}
#'   \item{\code{sides}}{As provided.}
#'   \item{\code{est}}{Posterior estimate of the regression surface.}
#'   \item{\code{est.FUN}}{As provided.}
#'   \item{\code{var}}{Summary of posterior variability of the regression
#'                     surface.}
#'   \item{\code{var.FUN}}{As provided.}
#'   \item{\code{W}}{An estimate of the extremal errors.}
#'   \item{\code{W.crit}}{The critical quantile of W.}
#'   \item{\code{trace}}{The posterior samples of the regression surface
#'                       indicated by the \code{track} argument.}
#' }
#' 
#' @examples
#' ### Sample from regression surface posterior
#' reg.surf.sample <- matrix(rnorm(1000, mean=1:10), ncol=2, byrow=TRUE)
#' credsubs(reg.surf.sample, cred.level=0.80)
#' 
#' ### Parametric case
#' design <- cbind(1, 1:10)
#' params <- matrix(rnorm(200, mean=1:2), ncol=2, byrow=TRUE)
#' credsubs(params, design)
#' 
#' ### With custom function
#' params.sd <- cbind(1 / rgamma(100, 1), params)
#' FUN.sd <- function(x, params) { params[, -1] %*% t(x) / params[, 1] }
#' credsubs(params.sd, design, FUN.sd, threshold=1)
#' 
#' @export
credsubs <- function(params,
                     design=NULL,
                     FUN=function(x, params) { params %*% t(x) },
                     cred.level=0.95,
                     threshold=0,
                     method=c('asymptotic', 'quantile'),
                     step.down=TRUE,
                     sides=c('both', 'exclusive', 'inclusive'),
                     est.FUN=mean,
                     var.FUN=sd,
                     point.estimate=NULL,
                     track=numeric(0),
                     verbose=FALSE) {
  
  # Validate and shape input
  params <- data.matrix(params)
  M <- nrow(params)
  
  if (is.null(design)) {
    N <- ncol(params)
    nonpar <- TRUE
  } else {
    design <- data.matrix(design)
    N <- nrow(design)
    nonpar <- FALSE
  }
  
  if (verbose) {
    cat("Computing credible subgroups over", N, "points using",
        M, "posterior draws.\n")
  }
  
  method <- method[1]
  stopifnot(method %in% c('asymptotic', 'quantile'))
  if (method == 'asymptotic') {
    m <- numeric(N)
    s <- numeric(N)
  }
  
  sides <- sides[1]
  stopifnot(sides %in% c('both', 'exclusive', 'inclusive'))
  scb.sides <- ifelse(sides == 'both', 'both',
                      ifelse(sides == 'exclusive', 'lower',
                             'upper'))
  
  credsubs <- list(exclusive=rep(NA, N),
                   inclusive=rep(NA, N),
                   cred.level=cred.level,
                   threshold=threshold,
                   method=method,
                   step.down=step.down,
                   sides=sides,
                   est=rep(NA, N),
                   est.FUN=est.FUN,
                   var=rep(NA, N),
                   var.FUN=var.FUN,
                   W=rep(NA, M),
                   W.crit=NA,
                   trace=matrix(NA, nrow=M, ncol=length(track)))
  class(credsubs) <- 'credsubs'
  
  credsubs$exclusive <- rep(FALSE, N)
  credsubs$inclusive <- rep(TRUE, N)
  
  test.set <- 1:N
  reject.set <- numeric(0)
  
  repeat {
    test.set <- setdiff(test.set, reject.set)
    if (nonpar) {
      test.par <- test.set
    } else {
      test.par <- 1:ncol(params)
    }
    
    # Runs on first iteration only
    if (length(test.set) == N) {
      scb.track <- track
    } else {
      scb.track <- numeric(0)
    }
    
    sim.cred.band <- sim.cred.band(params=params[, test.par, drop=FALSE],
                                   design=design[test.set, , drop=FALSE],
                                   FUN=FUN,
                                   cred.level=cred.level,
                                   method=method,
                                   sides=scb.sides,
                                   est.FUN=est.FUN,
                                   var.FUN=var.FUN,
                                   point.estimate=point.estimate,
                                   track=scb.track,
                                   verbose=verbose)
    
    over.set  <- test.set[sim.cred.band$lower > threshold]
    under.set <- test.set[sim.cred.band$upper < threshold]
    
    credsubs$exclusive[over.set]  <- TRUE
    credsubs$inclusive[under.set] <- FALSE
    
    credsubs$est[test.set] <- sim.cred.band$est
    credsubs$var[test.set] <- sim.cred.band$var
    
    credsubs$W.crit <- sim.cred.band$W.crit
    credsubs$W <- sim.cred.band$W
    
    reject.set <- union(over.set, under.set)
    
    # Runs on first iteration only
    if (length(test.set) == N) {
      credsubs$trace <- sim.cred.band$trace
    }
    
    if (verbose) {
      cat(length(test.set), "hypotheses tested,",
          length(reject.set), "rejected.\n")
    }
    
    if (!step.down ||
        length(reject.set) == 0 ||
        length(test.set) == length(reject.set)) {
      break
    }
  }
  
  credsubs
}

#' Compute the maximum credible levels at which conclusions may be drawn
#' 
#' For each covariate point, \code{credsubs.level} computes the maximum
#' credible level at which a conclusion may be drawn at each point, and
#' what that conclusion is.
#' 
#' If design is NULL (default), it is taken to be the identity
#' matrix of dimension ncol(params), so that the rows of params
#' are treated as draws from the posterior FUN(x, params).
#'
#' The 'asymptotic' method assumes that the marginal posteriors of 
#' the FUN(x, params) are asymptotically normal and is usually
#' significantly faster and less memory-intensive than the 'quantile'
#' method, which makes no such assumption.
#' 
#' By default (\code{z.store = "ram"}), the maximum credible level computation
#' stores a potentially very large amount of intermediate computation results
#' in memory. If not enough memory is available, \code{z.store = "disk"}
#' uses the \code{ff} package to store the intermediate results on disk,
#' which can still be fairly quick if the storage is fast (e.g. a local SSD).
#' Alternatively, \code{z.store = "recompute"} discards the intermediate
#' results and recomputes whenever needed. This uses minimal memory, but
#' is usually the slowest option.
#' 
#' @param params A numeric matrix whose rows are draws from the posterior
#'               distribution of either the regression surface or the
#'               parameter vector.
#' @param design (Optional) A numeric matrix whose rows are covariate points
#'               over which the band is to be constructed.
#' @param FUN (Optional) a function of the form \code{function(x, params)}
#'            that takes a row of \code{design} and the entire \code{params}
#'            matrix and returns a vector of the same length of \code{x}
#'            representing the regression surface.
#' @param threshold Numeric; the value of \code{FUN} above which
#'                  a covariate is included in the target subset.
#' @param method Either "asymptotic" (default) or "quantile"; see details.
#' @param step.down Logical (default \code{TRUE}); should the step-down
#'                  procedure be used?
#' @param sides One of "both" (default), "exclusive", or "inclusive".
#'              Which bounds should be constructed?
#' @param est.FUN The function used to produce estimates of the regression
#'                surface. Default is \code{mean}.
#' @param var.FUN The function used to quantify the variability of the
#'                regression surface posterior. Default is \code{sd}.
#' @param point.estimate If not null, replaces the mean and sets the reference 
#'                       around which the standard error is computed.
#'                       Useful for bootstrapping methods.
#'                       Treated as a row of the \code{params} matrix.
#' @param track A numeric vector of indices indicating which rows (default none)
#'              of the design matrix should have the sample of the corresponding
#'              \code{FUN(x, params)} returned.
#' @param verbose Logical (default \code{FALSE}); print progress?
#' @param z.store How should certain intermediate computations be handled?
#'                See details.
#' 
#' @return An object of class \code{credsubs.level}, which contains:
#' \describe{
#'   \item{\code{level}}{A numeric vector indicating the maximum credible
#'                       level at which a conclusion may be drawn at each
#'                       covariate point.}
#'   \item{\code{sign}}{A numeric vector indicating the which credible subsets
#'                      of which each covariate point is a member at the
#'                      credible level indicated by \code{level}. Exclusive
#'                      and inclusive: 1, inclusive only: 0, neither: -1.}
#'   \item{\code{threshold}}{As provided.}
#'   \item{\code{method}}{As provided.}
#'   \item{\code{step.down}}{As provided.}
#'   \item{\code{sides}}{As provided.}
#'   \item{\code{est}}{Posterior estimate of the regression surface.}
#'   \item{\code{est.FUN}}{As provided.}
#'   \item{\code{var}}{Summary of posterior variability of the regression
#'                     surface.}
#'   \item{\code{var.FUN}}{As provided.}
#'   \item{\code{trace}}{The posterior samples of the regression surface
#'                       indicated by the \code{track} argument.}
#' }
#' 
#' @examples
#' ### Sample from regression surface posterior
#' reg.surf.sample <- matrix(rnorm(1000, mean=1:10), ncol=2, byrow=TRUE)
#' credsubs.level(reg.surf.sample)
#' 
#' ### Parametric case
#' design <- cbind(1, 1:10)
#' params <- matrix(rnorm(200, mean=1:2), ncol=2, byrow=TRUE)
#' credsubs(params, design)
#' 
#' ### With custom function
#' params.sd <- cbind(1 / rgamma(100, 1), params)
#' FUN.sd <- function(x, params) { params[, -1] %*% t(x) / params[, 1] }
#' credsubs(params.sd, design, FUN.sd, threshold=1)
#' 
#' @export
credsubs.level <- function(params, design=NULL,
                           FUN=function(x, params) { params %*% t(x) },
                           threshold=0,
                           method=c('asymptotic', 'quantile'),
                           step.down=TRUE,
                           sides=c('both', 'exclusive', 'inclusive'),
                           est.FUN=mean,
                           var.FUN=sd,
                           point.estimate=NULL,
                           track=numeric(0),
                           verbose=FALSE,
                           z.store=c("ram", "recompute", "disk")) {
  
  # Validate and shape input
  params <- data.matrix(params)
  M <- nrow(params)
  
  if (is.null(design)) {
    N <- ncol(params)
    nonpar <- TRUE
  } else {
    design <- data.matrix(design)
    N <- nrow(design)
    nonpar <- FALSE
  }
  
  if (verbose) {
    cat("Finding maximum credible level at", N, "points using",
        M, "posterior draws.\n")
  }
  
  method <- method[1]
  stopifnot(method %in% c('asymptotic', 'quantile'))
  if (method == 'asymptotic') {
    m <- numeric(N)
    s <- numeric(N)
  } else if (method == 'quantile') {
    Fxt <- Gxt <- numeric(N)
  }
  
  sides <- sides[1]
  stopifnot(sides %in% c('both', 'exclusive', 'inclusive'))
  
  credsubs.level <- list(level=rep(NA, N),
                         sign=rep(NA, N),
                         threshold=threshold,
                         method=method,
                         step.down=step.down,
                         sides=sides)
  class(credsubs.level) <- 'credsubs.level'
  
  if (z.store[1] == "ram") {
    z.store <- matrix(0, nrow=M, ncol=N)
    recompute.z <- FALSE
  } else if (z.store[1] == "disk") {
    if (requireNamespace("ff", quietly=TRUE)) {
      z.store <- ff::ff(0, dim=c(M, N))
      recompute.z <- FALSE
    } else {
      warning("Package `ff`` required to use function 'credsubs.level'",
              "with option z.store='disk' set. Proceeding using",
              "z.store='recompute' instead.")
      recompute.z <- TRUE
    }
  } else if (z.store[1] == "recompute") {
    recompute.z <- TRUE
  } else {
    warning("option z.store must be one of 'ram', 'recompute', or 'disk'.",
            "Given: ", z.store[1])
    return(sim.cred.band)
  }
  
  W <- rep(-Inf, M)
  max.i <- numeric(M)
  m <- s <- est <- numeric(N)
  q <- sgn <- rep(NA, N)
  test.set <- 1:N
  
  for (i in 1:N) {

    if (verbose & i %% 100 == 0) {
      cat("prep", i, "/", N, "\n")
    }
    
    if (nonpar) {
      fx <- params[, i]
    } else {
      fx <- FUN(design[i, , drop=FALSE], params)
    }
    
    if (method == 'asymptotic') {
      
      if (is.null(point.estimate)) {
        m[i] <- mean(fx)
        s[i] <- sd(fx)
      } else {
        if (nonpar) {
          point.fx <- point.estimate[i]
        } else {
          point.fx <- FUN(design[i, , drop=FALSE], point.estimate)
        }
        m[i] <- point.fx
        s[i] <- sqrt(mean((fx - point.fx) ^ 2))
      }
      
      est[i] <- m[i]
      if (sides == "both") {
        z <- abs(fx - m[i]) / s[i]
      } else if (sides == "inclusive") {
        z <- (fx - m[i]) / s[i]
      } else if (sides == "exclusive") {
        z <- (m[i] - fx) / s[i]
      }
    } else {
      est[i] <- median(fx)
      Fxt[i] <- mean(fx <= threshold)
      Gxt[i] <- mean(fx <  threshold)
      
      if (sides == "both") {
        Fx <- to.Fx(fx)
        Gx <- 1 - to.Fx(-fx)
        z <- pmax(1-Fx, Gx)
      } else if (sides == "inclusive") {
        z <- 1 - to.Fx(-fx)
      } else if (sides == "exclusive") {
        z <- 1 - to.Fx(fx)
      }
    }
    
    if (!recompute.z) {
      z.store[, i] <- z
    }
    
    max.i <- ifelse(z > W, i, max.i)
    W <- pmax(W, z)
  }
  
  sgn <- ifelse(est > threshold &
                  (sides == 'both' | sides == 'exclusive'),
                1,
                ifelse(est < threshold &
                         (sides == 'both' | sides == 'inclusive'),
                       -1,
                       0))
  
  q[sgn == 0] <- 1
  test.set <- which(sgn != 0)
  
  first.max.i <- max.i
  
  prev.q <- 0
  recompute.m <- numeric(0)
  
  while (length(test.set) > 0) {
    
    if (verbose & (N - length(test.set)) %% 10 == 0) {
      cat("compute", N - length(test.set), "/", N,
          "recompute", length(recompute.m), "\n")
    }
    
    # Update W values
    if (length(recompute.m) > 0) {
      if (recompute.z) {
        if (method == 'asymptotic') {
          rcm <- recompute.m
        } else {
          rcm <- 1:M
        }
        if (nonpar) {
          fx <- params[rcm, test.set, drop=FALSE]
        } else {
          fx <- FUN(design[test.set, , drop=FALSE],
                    params[rcm, , drop=FALSE])
        }
        
        if (method == 'asymptotic') {
          if (sides == "both") {
            Z <- t(abs(t(fx) - m[test.set]) / s[test.set])
          } else if (sides == "inclusive") {
            Z <- t((t(fx) - m[test.set]) / s[test.set])
          } else if (sides == "exclusive") {
            Z <- t((m[test.set] - t(fx)) / s[test.set])
          }
        } else { # method == quantile
          if (sides == "both") {
            Fx <- t(apply(fx[recompute.m, , drop=FALSE], 1,
                          function(x, m) {
                            colMeans(t(t(m) <= x))
                          }, m=fx))
            Gx <- 1 - t(apply(-fx[recompute.m, , drop=FALSE], 1,
                            function(x, m) {
                              colMeans(t(t(m) <= x))
                            }, m=-fx))
            Z <- pmax(1-Fx, Gx)
          } else if (sides == "inclusive") {
            Z <- 1 - t(apply(-fx[recompute.m, , drop=FALSE], 1,
                           function(x, m) {
                             colMeans(t(t(m) <= x))
                           }, m=-fx))
          } else if (sides == "exclusive") {
            Z <- 1 - t(apply(fx[recompute.m, , drop=FALSE], 1,
                             function(x, m) {
                               colMeans(t(t(m) <= x))
                             }, m=fx))
          }
        }
      } else {
        Z <- z.store[recompute.m, test.set, drop=FALSE]
      }
      
      max.i[recompute.m] <- apply(Z, 1, function(z) { test.set[which.max(z)] })
      W[recompute.m] <- apply(Z, 1, max)
    }
    
    # Update empirical distribution of W
    FW <- ecdf(W)
    
    # Compute adjusted p-values
    if (method == 'asymptotic') {
      if (sides == 'both' | sides == 'exclusive') {
        lower <- 1 - FW((m[test.set] - threshold) / s[test.set])
      } else {
        lower <- rep(1, length(test.set))
      }
      
      if (sides == 'both' | sides == 'inclusive') {
        upper <- 1 - FW((threshold - m[test.set]) / s[test.set])
      } else {
        upper <- rep(1, length(test.set))
      }
    } else {
      if (sides == 'both' | sides == 'exclusive') {
        lower <- 1 - FW(1 - Fxt[test.set])
      } else {
        lower <- rep(1, length(test.set))
      }
      
      if (sides == 'both' | sides == 'inclusive') {
        upper <- 1 - FW(Gxt[test.set])
      } else {
        upper <- rep(1, length(test.set))
      }
    }
    
    p <- pmin(lower, upper)
    
    if (step.down) {
      lowest.p.i <- test.set[which.min(p)]
      q[lowest.p.i] <- prev.q <- max(prev.q, p[which(test.set == lowest.p.i)])
    } else {
      q <- p
      break
    }
    
    # these are the MCMC draws of W that change
    # when the selected covariate point is removed
    
    recompute.m <- which(max.i == lowest.p.i)
    
    test.set <- test.set[which(test.set != lowest.p.i)]
    
    i <- N - length(test.set)
    if (verbose && i %% 100 == 0) {
      cat("step", i, "/", N, "(", length(recompute.m), "/", M, ")\n")
    }
    
  }
  
  credsubs.level$level <- 1 - q
  credsubs.level$sign <- sgn
  credsubs.level$threshold <- threshold
  credsubs.level$step.down <- step.down
  
  credsubs.level
}

#' Build a credible subset calculator
#' 
#' This function builds a \code{shiny} application in the specified directory
#' that gives the maximum credible level at an entered covariate point.
#' 
#' The calculator creates a subdirectory according to \code{name} in
#' the directory specified by \code{dir}, and places in it files
#' \code{server.R}, \code{ui.R}, and \code{config.RData}. This application
#' requires the \code{shiny} package to run, and can be executed by passing
#' the directory path to \code{run.shiny.calc()}. The produced application
#' directory may be moved from its original location.
#' 
#' @param credsubs.level An object of class \code{credsubs.level}.
#' @param cov.space A data frame whose rows are human-readable
#'                  covariate points corresponding to the entries
#'                  of \code{credsubs.level$level}.
#' @param name A character string indicating the name of the application.
#' @param dir The directory in which to place the application.
#' @param title A character string to be displayed as the application title.
#' @param instructions A character string to be displayed as instructions.
#'                     HTML allowed.
#' 
#' @export
build.shiny.calc <- function(credsubs.level,
                             cov.space,
                             name="calc",
                             dir=".",
                             title="Credible Subsets Calculator",
                             instructions="Select a covariate point.") {
  
  # Shape and validate input
  cov.space <- as.data.frame(cov.space)
  stopifnot(length(credsubs.level$level) == nrow(cov.space))
  
  full.name <- paste0(dir, "/", name)
  if (!dir.exists(full.name)) {
    dir.create(paste0(full.name))
  }
  
  shiny.dir <- system.file("shiny", "calc", package="credsubs")
  file.copy(paste0(shiny.dir, "/server.R"),
            paste0(full.name, "/server.R"),
            overwrite=TRUE)
  file.copy(paste0(shiny.dir, "/ui.R"),
            paste0(full.name, "/ui.R"),
            overwrite=TRUE)
  
  # Save file
  cat(paste0("Building calculator in ",
             full.name, "/", "\n"))
  cat(paste0("Execute run.shiny.calc('", full.name , "') to run.", "\n"))
  save(credsubs.level, cov.space, title, instructions,
       file=paste0(full.name, "/config.RData"))
}

#' Run a calculator
#' 
#' Runs the specified calculator using the \code{shiny} package.
#' Calculators must be built using \code{build.shiny.calc}.
#' 
#' The \code{app.dir} argument need not exactly match the value recommended
#' by \code{build.shiny.calc()}, as long as it points to the correct directory.
#' For example, \code{"./calc/"}, \code{"calc/"}, and \code{"calc"} are
#' all equivalent. If no value is supplied, an example is run.
#' 
#' @param app.dir A character string pointing to the application directory.
#' 
#' @export
run.shiny.calc <- function(app.dir=system.file("shiny",
                                               "alzheimers",
                                               package = "credsubs")) {
  if (requireNamespace("shiny", quietly=TRUE)) {
    shiny::runApp(app.dir, display.mode = "normal")
  } else {
    warning("Package `shiny` needed for this function to work.", " ",
            "Please install it.")
  }
}

#' @export
plot.sim.cred.band <- function(x, ...) {
  sim.cred.band <- x
  
  o <- order(sim.cred.band$est)
  plot(sim.cred.band$est[o],
       ylim=c(min(sim.cred.band$lower), max(sim.cred.band$upper)),
       xlab="Rank", ylab="Value",
       ...)
  points(sim.cred.band$lower[o], ...)
  points(sim.cred.band$upper[o], ...)
}

#' @export
plot.credsubs <- function(x, ...) {
  credsubs <- x
  
  # Validate method, estimates, and variabilities
  stopifnot(credsubs$method[1] == "asymptotic")
  stopifnot(identical(credsubs$est.FUN, mean))
  stopifnot(identical(credsubs$var.FUN, sd))
  
  plot(credsubs$est - credsubs$threshold, credsubs$var,
       col=ifelse(
         credsubs$est -
           credsubs$W.crit * credsubs$var > credsubs$threshold,
         "forestgreen", ifelse(
           credsubs$est +
             credsubs$W.crit * credsubs$var >= credsubs$threshold,
           "gold", "firebrick"
         )),
       xlab="Mean - threshold", ylab="Standard deviation", ...)
  
  abline(v=0, lty=2)
  abline(h=0, lty=2)
  abline(b=1/credsubs$W.crit, a=0)
  abline(b=-1/credsubs$W.crit, a=0)
}

#' @export
plot.credsubs.level <- function(x, ...) {
  credsubs.level <- x
  
  signed.level <- credsubs.level$level * credsubs.level$sign
  plot(sort(signed.level), xlab="Rank", ylab="Signed credible level", ...)
  abline(h=0, lty=2)
}

#' @export
print.sim.cred.band <- function(x, ...) {
  sim.cred.band <- x
  
  cat("Simultaneous credible band\n")
  cat("\n")
  cat("Credible level:", sim.cred.band$cred.level, "\n")
  cat("Method:", sim.cred.band$method, "\n")
  cat("Critical W:", sim.cred.band$W.crit, "\n")
  cat("\n")
  cat("Upper bounds:", head(sim.cred.band$upper),
      ifelse(length(sim.cred.band$upper) > 6, "...\n", "\n"))
  cat("Lower bounds:", head(sim.cred.band$lower),
      ifelse(length(sim.cred.band$lower) > 6, "...\n", "\n"))
}

#' @export
print.credsubs <- function(x, ...) {
  credsubs <- x
  
  cat("Credible subsets\n")
  cat("\n")
  cat("Credible level:", credsubs$cred.level, "\n")
  cat("Threshold:", credsubs$threshold, "\n")
  cat("Method:", credsubs$method, "\n")
  cat("Step-down:", credsubs$step.down, "\n")
  cat("Critical W:", credsubs$W.crit, "\n")
  cat("\n")
  cat("Exclusive credible subset:", head(credsubs$exclusive),
      ifelse(length(credsubs$exclusive) > 6, "...\n", "\n"))
  cat("Inclusive credible subset:", head(credsubs$inclusive),
      ifelse(length(credsubs$inclusive) > 6, "...\n", "\n"))
}

#' @export
print.credsubs.level <- function(x, ...) {
  credsubs.level <- x
  cat("Maximum credible levels\n")
  cat("\n")
  cat("Threshold:", credsubs.level$threshold, "\n")
  cat("Method:", credsubs.level$method, "\n")
  cat("Step-down:", credsubs.level$step.down, "\n")
  cat("\n")
  cat("Maximum levels:", head(credsubs.level$level),
      ifelse(length(credsubs.level$level) > 6, "...\n", "\n"))
  cat("Signs:", head(credsubs.level$sign),
      ifelse(length(credsubs.level$sign) > 6, "...\n", "\n"))
}
