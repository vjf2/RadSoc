#Customized netlogit function

netlogitM<-function(y, x, intercept=TRUE, mode="digraph", diag=FALSE,
         nullhyp=c("qap", "qapspp", "qapy", "qapx", "qapallx", 
                   "cugtie", "cugden", "cuguman", "classical"), test.statistic = 
           c("z-value","beta"), tol=1e-7, reps=1000, weights)
{
  gfit <- function(glist, mode, diag) {
    y <- gvectorize(glist[[1]], mode = mode, diag = diag, 
                    censor.as.na = TRUE)
    x <- vector()
    for (i in 2:length(glist)) x <- cbind(x, gvectorize(glist[[i]], 
                                                        mode = mode, diag = diag, censor.as.na = TRUE))
    if (!is.matrix(x)) 
      x <- matrix(x, ncol = 1)
    mis <- is.na(y) | apply(is.na(x), 1, any)
    glm.fit(x[!mis, ], y[!mis], weights=weights[!mis], family = binomial(), intercept = FALSE)
  }
  
  gfitlm <- function(glist, mode, diag, tol) {
    y <- gvectorize(glist[[1]], mode = mode, diag = diag, 
                    censor.as.na = TRUE)
    x <- vector()
    for (i in 2:length(glist)) x <- cbind(x, gvectorize(glist[[i]], 
                                                        mode = mode, diag = diag, censor.as.na = TRUE))
    if (!is.matrix(x)) 
      x <- matrix(x, ncol = 1)
    mis <- is.na(y) | apply(is.na(x), 1, any)
    list(qr(x[!mis, ], tol = tol), y[!mis])
  }
  y <- as.sociomatrix.sna(y)
  x <- as.sociomatrix.sna(x)
  if (is.list(y) || ((length(dim(y)) > 2) && (dim(y)[1] > 1))) 
    stop("y must be a single graph in netlogit.")
  if (length(dim(y)) > 2) 
    y <- y[1, , ]
  if (is.list(x) || (dim(x)[2] != dim(y)[2])) 
    stop("Homogeneous graph orders required in netlogit.")
  nx <- stackcount(x) + intercept
  n <- dim(y)[2]
  g <- list(y)
  if (intercept) 
    g[[2]] <- matrix(1, n, n)
  if (nx - intercept == 1){g[[2 + intercept]] <- x
  }else for (i in 1:(nx - intercept)) g[[i + 1 + intercept]] <- x[i, 
                                                                 , ]
  if (any(sapply(lapply(g, is.na), any))) 
    warning("Missing data supplied to netlogit; this may pose problems for certain null hypotheses.  Hope you know what you're doing....")
  
  fit.base <- gfit(g, mode = mode, diag = diag)
  fit <- list()
  fit$coefficients <- fit.base$coefficients
  fit$fitted.values <- fit.base$fitted.values
  fit$residuals <- fit.base$residuals
  fit$se <- sqrt(diag(chol2inv(fit.base$qr$qr)))
  tstat <- match.arg(test.statistic)
  fit$test.statistic <- tstat
  if (tstat == "beta") 
    fit$tstat <- fit$coefficients
  else if (tstat == "z-value") 
    fit$tstat <- fit$coefficients/fit$se
  fit$linear.predictors <- fit.base$linear.predictors
  fit$n <- length(fit.base$y)
  fit$df.model <- fit.base$rank
  fit$df.residual <- fit.base$df.residual
  fit$deviance <- fit.base$deviance
  fit$null.deviance <- fit.base$null.deviance
  fit$df.null <- fit.base$df.null
  fit$aic <- fit.base$aic
  fit$bic <- fit$deviance + fit$df.model * log(fit$n)
  fit$qr <- fit.base$qr
  fit$ctable <- table(as.numeric(fit$fitted.values >= 0.5), 
                      fit.base$y, dnn = c("Predicted", "Actual"))
  # if (NROW(fit$ctable) == 1) {
  #   if (rownames(fit$ctable) == "0") 
  #     fit$ctable <- rbind(fit$ctable, c(0, 0))
  #   else fit$ctable <- rbind(c(0, 0), fit$ctable)
  #   rownames(fit$ctable) <- c("0", "1")
  # }
  nullhyp <- match.arg(nullhyp)
  if ((nullhyp %in% c("qap", "qapspp")) && (nx == 1)) 
    nullhyp <- "qapy"
  if (nullhyp == "classical") {
    cvm <- chol2inv(fit$qr$qr)
    se <- sqrt(diag(cvm))
    tval <- fit$coefficients/se
    fit$dist <- NULL
    fit$pleeq <- pt(tval, fit$df.residual)
    fit$pgreq <- pt(tval, fit$df.residual, lower.tail = FALSE)
    fit$pgreqabs <- 2 * pt(abs(tval), fit$df.residual, lower.tail = FALSE)
  }
  else if (nullhyp %in% c("cugtie", "cugden", "cuguman")) {
    repdist <- matrix(0, reps, nx)
    for (i in 1:nx) {
      gr <- g
      for (j in 1:reps) {
        gr[[i + 1]] <- switch(nullhyp, cugtie = rgraph(n, 
                                                       mode = mode, diag = diag, replace = FALSE, 
                                                       tielist = g[[i + 1]]), cugden = rgraph(n, tprob = gden(g[[i + 
                                                                                                                   1]], mode = mode, diag = diag), mode = mode, 
                                                                                              diag = diag), cuguman = (function(dc, n) {
                                                                                                rguman(1, n, mut = dc[1], asym = dc[2], null = dc[3], 
                                                                                                       method = "exact")
                                                                                              })(dyad.census(g[[i + 1]]), n))
        repfit <- gfit(gr, mode = mode, diag = diag)
        if (tstat == "beta") 
          repdist[j, i] <- repfit$coef[i]
        else repdist[j, i] <- repfit$coef[i]/sqrt(diag(chol2inv(repfit$qr$qr)))[i]
      }
    }
    fit$dist <- repdist
    fit$pleeq <- apply(sweep(fit$dist, 2, fit$tstat, "<="), 
                       2, mean)
    fit$pgreq <- apply(sweep(fit$dist, 2, fit$tstat, ">="), 
                       2, mean)
    fit$pgreqabs <- apply(sweep(abs(fit$dist), 2, abs(fit$tstat), 
                                ">="), 2, mean)
  }
  else if (nullhyp == "qapy") {
    repdist <- matrix(0, reps, nx)
    gr <- g
    for (i in 1:reps) {
      gr[[1]] <- rmperm(g[[1]])
      repfit <- gfit(gr, mode = mode, diag = diag)
      if (tstat == "beta") 
        repdist[i, ] <- repfit$coef
      else repdist[i, ] <- repfit$coef/sqrt(diag(chol2inv(repfit$qr$qr)))
    }
    fit$dist <- repdist
    fit$pleeq <- apply(sweep(fit$dist, 2, fit$tstat, "<="), 
                       2, mean)
    fit$pgreq <- apply(sweep(fit$dist, 2, fit$tstat, ">="), 
                       2, mean)
    fit$pgreqabs <- apply(sweep(abs(fit$dist), 2, abs(fit$tstat), 
                                ">="), 2, mean)
  }
  else if (nullhyp == "qapx") {
    repdist <- matrix(0, reps, nx)
    for (i in 1:nx) {
      gr <- g
      for (j in 1:reps) {
        gr[[i + 1]] <- rmperm(gr[[i + 1]])
        repfit <- gfit(gr, mode = mode, diag = diag)
        if (tstat == "beta") 
          repdist[j, i] <- repfit$coef[i]
        else repdist[j, i] <- repfit$coef[i]/sqrt(diag(chol2inv(repfit$qr$qr)))[i]
      }
    }
    fit$dist <- repdist
    fit$pleeq <- apply(sweep(fit$dist, 2, fit$tstat, "<="), 
                       2, mean)
    fit$pgreq <- apply(sweep(fit$dist, 2, fit$tstat, ">="), 
                       2, mean)
    fit$pgreqabs <- apply(sweep(abs(fit$dist), 2, abs(fit$tstat), 
                                ">="), 2, mean)
  }
  else if (nullhyp == "qapallx") {
    repdist <- matrix(0, reps, nx)
    gr <- g
    for (i in 1:reps) {
      for (j in 1:nx) gr[[1 + j]] <- rmperm(g[[1 + j]])
      repfit <- gfit(gr, mode = mode, diag = diag)
      if (tstat == "beta") 
        repdist[i, ] <- repfit$coef
      else repdist[i, ] <- repfit$coef/sqrt(diag(chol2inv(repfit$qr$qr)))
    }
    fit$dist <- repdist
    fit$pleeq <- apply(sweep(fit$dist, 2, fit$tstat, "<="), 
                       2, mean)
    fit$pgreq <- apply(sweep(fit$dist, 2, fit$tstat, ">="), 
                       2, mean)
    fit$pgreqabs <- apply(sweep(abs(fit$dist), 2, abs(fit$tstat), 
                                ">="), 2, mean)
  }
  else if ((nullhyp == "qap") || (nullhyp == "qapspp")) {
    xsel <- matrix(TRUE, n, n)
    if (!diag) 
      diag(xsel) <- FALSE
    if (mode == "graph") 
      xsel[upper.tri(xsel)] <- FALSE
    repdist <- matrix(0, reps, nx)
    for (i in 1:nx) {
      xfit <- gfitlm(g[1 + c(i, (1:nx)[-i])], mode = mode, 
                     diag = diag, tol = tol)
      xres <- g[[1 + i]]
      xres[xsel] <- qr.resid(xfit[[1]], xfit[[2]])
      if (mode == "graph") 
        xres[upper.tri(xres)] <- t(xres)[upper.tri(xres)]
      for (j in 1:reps) {
        repfit <- gfit(c(g[-(1 + i)], list(rmperm(xres))), 
                       mode = mode, diag = diag)
        if (tstat == "beta") 
          repdist[j, i] <- repfit$coef[nx]
        else repdist[j, i] <- repfit$coef[nx]/sqrt(diag(chol2inv(repfit$qr$qr)))[nx]
      }
    }
    fit$dist <- repdist
    fit$pleeq <- apply(sweep(fit$dist, 2, fit$tstat, "<="), 
                       2, mean)
    fit$pgreq <- apply(sweep(fit$dist, 2, fit$tstat, ">="), 
                       2, mean)
    fit$pgreqabs <- apply(sweep(abs(fit$dist), 2, abs(fit$tstat), 
                                ">="), 2, mean)
  }
  fit$nullhyp <- nullhyp
  fit$names <- paste("x", 1:(nx - intercept), sep = "")
  if (intercept) 
    fit$names <- c("(intercept)", fit$names)
  fit$intercept <- intercept
  class(fit) <- "netlogit"
  fit
}
