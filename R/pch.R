#' Generate sample from the piecewise constant hazard model
#'
#' @param n sample size
#' @param cuts cuts defining the intervals in the piecewise constant hazard model
#' @param alpha values of the hazard (over each interval) in the piecewise constant hazard model
#' @return a matrix of dimension \code{K * J}
#' @examples
#' \dontrun{
#' data <- rpch(1000, 1, c(0.2, 0.9))
#' plot(ecdf(data))
#' }
#' @export
rpch <- function(n, cuts, alpha) {
  u <- runif(n)
  k <- length(alpha)
  if (length(cuts) != (k - 1))
    stop("error: length(cp) must be equal to length(alpha)-1 !")
  cuts0 <- c(0, cuts)
  if (k > 1) {
    thresh <- exp(-cumsum(alpha[-k] * diff(cuts0)))
    if (n <= 200) {
      seg <- apply(matrix(rep(thresh, n), byrow = T, nrow = n) > u, 1, sum) +
        1
    } else {
      seg <- rep(NA, n)
      for (i in 1:n) seg[i] <- sum(thresh > u[i]) + 1
    }
  } else {
    seg <- rep(1, n)
  }
  res <- cuts0[seg] - (log(u) + cumsum(c(0, alpha[-k] * diff(cuts0)))[seg]) / alpha[seg]
  return(res)
}

pch.cumhaz <- function(time, cuts, alpha) {
  if (length(cuts) == 0)
    return(time * alpha[1])
  k <- length(alpha)
  I <- k - apply(matrix(rep(time, k - 1), byrow = TRUE, nrow = k - 1) < cuts, 2,
                 sum)
  cuts0 <- c(0, cuts)
  return(alpha[I] * (time - cuts0[I]) + c(0, cumsum(alpha * c(diff(cuts0), 0))[-k])[I])
}


mle.pchsurv <- function(start, end, delta, cuts, weights = NULL) {
  if (is.null(weights))
    weights <- rep(1, length(end))
  k <- length(cuts) + 1
  cuts0 <- c(0, cuts)
  # The exhaustive statistics
  Jend <- cut(end, breaks = c(0, cuts, Inf), labels = 1:k, right = FALSE)
  A <- sapply(split(weights * delta, Jend), sum)
  aux <- sapply(split(weights, Jend), sum)[-1]
  Bend <- sapply(split(weights * (end - cuts0[Jend]), Jend), sum)
  Bend[-k] <- Bend[-k] + rev(cumsum(rev(aux))) * diff(cuts0)
  # The exhaustive statistics
  Jstart <- cut(start, breaks = c(0, cuts, Inf), labels = 1:k, right = FALSE)
  aux <- sapply(split(weights, Jstart), sum)[-1]
  Bstart <- sapply(split(weights * (start - cuts0[Jstart]), Jstart), sum)
  Bstart[-k] <- Bstart[-k] + rev(cumsum(rev(aux))) * diff(cuts0)
  B <- Bend - Bstart
  a <- log(A / B)
  hazard <- NULL
  if (length(cuts) > 0)
    hazard <- stepfun(cuts, exp(a))
  output <- (list(A = A, B = B, a = a, hazard = hazard, cuts = cuts, start = start,
                  end = end))
  class(output) <- "pchsurv"
  return(output)
}

loglik.pchsurv <- function(start, end, delta, cuts, a, weights = NULL) {
  if (is.null(weights)) {
    weights <- rep(1, length(end))
  }
  k <- length(cuts) + 1
  Jend <- cut(end, breaks = c(0, cuts, Inf), labels = 1:k, right = FALSE)
  Jstart <- cut(start, breaks = c(0, cuts, Inf), labels = 1:k, right = FALSE)
  cuts0 <- c(0, cuts)
  cumhazEnd <- exp(a)[Jend] * (end - cuts0[Jend]) + c(0, cumsum(exp(a) * c(diff(cuts0),
                                                                           0))[-k])[Jend]
  cumhazStart <- exp(a)[Jstart] * (start - cuts0[Jstart]) + c(0, cumsum(exp(a) *
                                                                          c(diff(cuts0), 0))[-k])[Jstart]
  return(list(loglik = weights * ((delta == 1) * a[Jend] - cumhazEnd + cumhazStart),
              alpha = exp(a)))
}
# Without left-truncation
mle.pchsurv <- function(time, delta, cuts, weights = NULL) {
  if (is.null(weights))
    weights <- rep(1, length(time))
  k <- length(cuts) + 1
  J <- cut(time, breaks = c(0, cuts, Inf), labels = 1:k, right = FALSE)
  cuts0 <- c(0, cuts)
  A <- sapply(split(weights * delta, J), sum)
  aux <- sapply(split(weights, J), sum)[-1]
  B <- sapply(split(weights * (time - cuts0[J]), J), sum)
  B[-k] <- B[-k] + rev(cumsum(rev(aux))) * diff(cuts0)
  a <- log(A/B)
  hazard <- NULL
  if (length(cuts) > 0)
    hazard <- stepfun(cuts, exp(a))
  return(list(a = a, hazard = hazard, Var = A/B^2))
}
loglik.pchsurv <- function(time, delta, cuts, a, weights = NULL) {
  if (is.null(weights))
    weights <- rep(1, length(time))
  k <- length(cuts) + 1
  J <- cut(time, breaks = c(0, cuts, Inf), labels = 1:k, right = FALSE)
  cuts0 <- c(0, cuts)
  cumhaz <- exp(a)[J] * (time - cuts0[J]) + c(0, cumsum(exp(a) * c(diff(cuts0),
                                                                   0))[-k])[J]
  return(sum(weights * ((delta == 1) * a[J] - cumhaz)))
}


plot.pchsurv <- function(fit, xlab = "Time", ylab = "Hazard", lwd = 2, main = "",
                         lty = 1, ...) {
  plot(c(0, fit$cuts, max(fit$end) + 10), exp(c(fit$a, fit$a[length(fit$a)])),
       type = "s", lwd = lwd, xlab = xlab, ylab = ylab, main = main, lty = lty,
       ...)
  # plot(end,hazard(end),type='S',lwd=lwd,xlab=xlab,ylab=ylab,main=main,...)
}

lines.pchsurv <- function(fit, lwd = 2, lty = 1, ...) {
  lines(c(0, fit$cuts, max(fit$end) + 10), c(fit$a, fit$a[length(fit$a)]), type = "s",
        lwd = lwd, lty = lty, ...)
}

logsumexp <- function(l) {
  i <- which.max(l)
  res <- l[i] + log1p(sum(exp(l[-i] - l[i])))
  if (is.nan(res))
    res <- -Inf
  return(res)
}
# eta = matrix(0,n,3); eta[which(diff(sort(birth))>0)+1,]<-0.5
FB <- function(le, eta = NULL, birth, graphical = TRUE) {
  n <- nrow(le)
  K <- ncol(le)
  if (is.null(eta)) {
    eta <- matrix(0, n, K)
    eta[which(diff(sort(birth)) > 0) + 1, ] <- 0.5
  }
  # eta[which(diff(sort(birth),)==0)+1,]<-0
  leta <- log(eta)
  l1meta <- log(1 - eta)
  lF <- lB <- matrix(NA, n, K)
  # forward
  lF[1, 1] <- le[1, 1]
  lF[1, -1] <- -Inf
  for (i in 2:n) {
    aux <- rbind(lF[i - 1, ] + l1meta[i, ], c(-Inf, lF[i - 1, -K]) + leta[i,
                                                                          ])
    lF[i, ] <- apply(aux, 2, logsumexp)
    for (k in 1:K) lF[i, k] <- logsumexp(aux[, k]) + le[i, k]
  }
  # backward
  lB[n, ] <- lB[n - 1, ] <- -Inf
  lB[n, K] <- 0
  lB[n - 1, K - 1] <- leta[n, K - 1] + le[n, K]
  lB[n - 1, K] <- l1meta[n, K] + le[n, K]
  for (i in seq(n - 1, 2, by = -1)) {
    aux <- rbind(l1meta[i, ] + le[i, ] + lB[i, ], c(leta[i, -1] + le[i, -1] +
                                                      lB[i, -1], -Inf))
    lB[i - 1, ] <- apply(aux, 2, logsumexp)
  }
  # verification
  diff(range(apply(lF + lB, 1, logsumexp)))
  # posterior distribution
  lpev <- lF[n, K]
  weights <- exp(lF + lB - lpev)
  # return results
  return(list(weights = weights, le = le, lF = lF, lB = lB, lpev = lpev, pev = exp(lpev),
              leta = leta))
}

# Without covariates
survseg <- function(start, end, delta, bp = NULL, cuts = NULL, itermax = 200, tol = 1e-10,
                    verbose = FALSE) {
  n <- length(delta)
  # only one segment, simple
  if (k == 1) {
    weights <- rep(1, length(start))
    # if (method == 'exp') exp_res = exp_loglik(start, end, delta) loglik = exp_res$loglik
    # alpha = exp_res$alpha if (method == 'pch')
    a <- mle.pchsurv(start, end, delta, cuts)$a
    pch_res <- loglik.pchsurv(start, end, delta, cuts, a)
    loglik <- pch_res$loglik
    alpha <- exp(a)
    return(list(start = start, end = end, delta = delta, k = k, n = n, loglik = sum(loglik),
                alpha = alpha))
  }
  # more than one, EM
  if (k > 1) {
    # initialize weights
    bp <- c(0, bp, n)
    p <- 0.7
    weights <- matrix(rep((1 - p)/(k - 1), k * n), ncol = k)
    alpha <- matrix(rep(0, k * (length(cuts) + 1)), ncol = k)
    for (j in 1:k) {
      weights[((1:n) > bp[j]) & ((1:n) <= bp[j + 1]), j] <- p
    }
    # main loop
    loglik <- 0
    le <- matrix(rep(NA, k * n), ncol = k)
    # alpha=rep(NA,k)
    seg0 <- FB(matrix(0, ncol = k, nrow = n), birth = birth)
    for (iter in 1:itermax) {
      # update log-evidence if (method=='exp') for (j in 1:k) {
      # exp_res=exp_loglik(start,end,delta,weights=weights[,j]) le[,j]=exp_res$loglik
      # alpha[j]=exp_res$alpha } if (method=='pch')
      for (j in 1:k) {
        a <- mle.pchsurv(start, end, delta, cuts, weights = weights[, j])$a
        pch_res <- loglik.pchsurv(start, end, delta, cuts, a = a, weights = weights[,
                                                                                    j])
        le[, j] <- pch_res$loglik
        alpha[, j] <- exp(a)
      }
      # forward/backward
      le[is.nan(le)] <- -Inf
      seg <- FB(le, birth = birth)
      weights <- seg$weights
      # update loglik
      oldloglik <- loglik
      loglik <- seg$lpev - seg0$lpev
      # test for convergence
      if ((abs(loglik - oldloglik)/abs(loglik)) < tol)
        break
      if (verbose)
        print(c(iter, loglik))
      # if (verbose) print(apply(weights,2,sum))
    }
    return(list(start = start, end = end, delta = delta, k = k, n = n, loglik = loglik,
                seg = seg, alpha = alpha))
  }
}


