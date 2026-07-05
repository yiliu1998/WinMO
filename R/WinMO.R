#' WinMO: Win estimands (WR/NB/WO/DOOR) with 1-3 ordinal endpoints under missingness
#'
#' Main entry point that dispatches to the corresponding IPW/AIPW influence-function
#' implementations for 1, 2, or 3 hierarchical ordinal endpoints.
#'
#' @details
#' Provide outcome column names via \code{Y}. If \code{length(Y)=1}, uses the 1-endpoint
#' implementation; if \code{length(Y)=2}, uses the 2-endpoint implementation; if
#' \code{length(Y)=3}, uses the 3-endpoint implementation.
#'
#' This function assumes the helper functions
#' \code{IPW_WinStat_1end_IF}, \code{AIPW_WinStat_1end_IF},
#' \code{IPW_WinStat_2end_IF}, \code{AIPW_WinStat_2end_IF},
#' \code{IPW_WinStat_3end_IF}, \code{AIPW_WinStat_3end_IF}
#' are available in the package namespace.
#'
#' @param data A \code{data.frame} containing \code{A}, \code{X}, and outcomes \code{Y}.
#' @param A Character scalar. Treatment column name (coded 0/1 or coercible to integer).
#' @param X Character vector. Covariate column names.
#' @param Y Character vector of length 1, 2, or 3 giving outcome column names in
#'   hierarchical order (primary, secondary, tertiary).
#' @param method Character. \code{"AIPW"} (default) or \code{"IPW"}.
#' @param eps Numeric. Truncation for probabilities to avoid division by zero.
#'
#' @return A named \code{list} with elements (at least)
#' \code{pW}, \code{pL}, \code{pT}, \code{WR}, \code{NB}, \code{WO}, \code{DOOR},
#' and \code{SE_IF}. Additional fields include \code{method}, \code{K}, \code{Y}, \code{call}.
#'
#' @examples
#' \dontrun{
#' # 1 endpoint
#' res1 <- WinMO(dat, A = "A", X = c("X1","X2"), Y = "Y1", method = "AIPW")
#'
#' # 2 endpoints
#' res2 <- WinMO(dat, A = "A", X = c("X1","X2"), Y = c("Y1","Y2"), method = "IPW")
#'
#' # 3 endpoints
#' res3 <- WinMO(dat, A = "A", X = c("X1","X2"), Y = c("Y1","Y2","Y3"))
#' }
#'
#' @export
WinMO <- function(data,
                  A,
                  X,
                  Y,
                  method = c("AIPW", "IPW"),
                  eps = 1e-6) {
  method <- match.arg(method)

  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!is.character(A) || length(A) != 1L) stop("`A` must be a single column name (character).")
  if (!A %in% names(data)) stop("`A` is not a column in `data`.")
  if (!is.character(X) || length(X) < 1L) stop("`X` must be a non-empty character vector.")
  if (!all(X %in% names(data))) stop("Some `X` columns are not in `data`.")
  if (!is.numeric(eps) || length(eps) != 1L || eps <= 0) stop("`eps` must be a positive scalar.")

  if (is.list(Y)) Y <- unlist(Y, use.names = FALSE)
  if (!is.character(Y)) stop("`Y` must be outcome column name(s) as character.")
  if (!(length(Y) %in% c(1L, 2L, 3L))) stop("`Y` must have length 1, 2, or 3.")
  if (!all(Y %in% names(data))) stop("Some `Y` columns are not in `data`.")

  K <- length(Y)

  out <- switch(
    paste0(method, "_", K),
    "IPW_1"  = IPW_WinStat_1end_IF(data = data, A = A, X = X, Y  = Y[1], eps = eps),
    "AIPW_1" = AIPW_WinStat_1end_IF(data = data, A = A, X = X, Y  = Y[1], eps = eps),

    "IPW_2"  = IPW_WinStat_2end_IF(data = data, A = A, X = X, Y1 = Y[1], Y2 = Y[2], eps = eps),
    "AIPW_2" = AIPW_WinStat_2end_IF(data = data, A = A, X = X, Y1 = Y[1], Y2 = Y[2], eps = eps),

    "IPW_3"  = IPW_WinStat_3end_IF(data = data, A = A, X = X, Y1 = Y[1], Y2 = Y[2], Y3 = Y[3], eps = eps),
    "AIPW_3" = AIPW_WinStat_3end_IF(data = data, A = A, X = X, Y1 = Y[1], Y2 = Y[2], Y3 = Y[3], eps = eps),

    stop("Unsupported combination of `method` and number of endpoints.")
  )

  out$call <- match.call()
  out$method <- method
  out$K <- K
  out$Y <- Y
  out
}


WinMO_logratio_SE <- function(pW, pL, pT, Psi, eps = 1e-6) {
  V <- crossprod(Psi / nrow(Psi))

  WR_hat <- if (pL > 0) pW / pL else NA_real_
  WO_hat <- (pW + 0.5 * pT) / (pL + 0.5 * pT)

  num_WR <- max(pW, eps)
  den_WR <- max(pL, eps)
  num_WO <- max(pW + 0.5 * pT, eps)
  den_WO <- max(pL + 0.5 * pT, eps)

  grad_log_WR <- c(1 / num_WR, -1 / den_WR, 0)
  grad_NB <- c(1, -1, 0)
  grad_log_WO <- c(1 / num_WO, -1 / den_WO, 0.5 / num_WO - 0.5 / den_WO)
  grad_DOOR <- c(1, 0, 0.5)

  G <- rbind(grad_log_WR, grad_NB, grad_log_WO, grad_DOOR)
  se_tmp <- sqrt(pmax(diag(G %*% V %*% t(G)), 0))

  SE_log <- c(WR = unname(se_tmp[1]), WO = unname(se_tmp[3]))
  SE <- c(WR = if (is.finite(WR_hat)) WR_hat * SE_log["WR"] else NA_real_,
          NB = unname(se_tmp[2]),
          WO = if (is.finite(WO_hat)) WO_hat * SE_log["WO"] else NA_real_,
          DOOR = unname(se_tmp[4]))
  list(SE = SE, SE_log = SE_log)
}

IPW_WinStat_1end_IF <- function(data, A, X, Y, eps = 1e-6) {
  n <- nrow(data)
  Av <- as.integer(data[[A]])
  Xmat <- cbind(1, as.matrix(data[, X, drop = FALSE]))
  p <- ncol(Xmat)
  Yv <- data[[Y]]
  Rv <- as.integer(!is.na(Yv))
  pa <- mean(Av == 1)
  lev <- sort(unique(Yv[!is.na(Yv)]))
  L <- length(lev)
  Yint <- match(Yv, lev)
  Yint[is.na(Yint)] <- 0

  fit_pi <- function(a) {
    idx <- which(Av == a)
    Xsub <- Xmat[idx, , drop = FALSE]
    Rsub <- Rv[idx]
    fit <- glm.fit(x = Xsub, y = Rsub, family = binomial())
    b <- fit$coefficients
    eta <- drop(Xmat %*% b)
    pihat <- pmin(pmax(plogis(eta), eps), 1 - eps)
    W <- as.numeric(plogis(drop(Xsub %*% b)) * (1 - plogis(drop(Xsub %*% b))))
    H <- crossprod(Xsub * sqrt(W))
    Hinv <- tryCatch(solve(H), error = function(e) MASS::ginv(H))
    list(pi = pihat, Hinv = Hinv)
  }

  f0 <- fit_pi(0); f1 <- fit_pi(1)
  pi0 <- f0$pi; pi1 <- f1$pi
  Hinv0 <- f0$Hinv; Hinv1 <- f1$Hinv

  w0 <- (Av == 0) / (1 - pa) * Rv / pi0
  w1 <- (Av == 1) / pa * Rv / pi1

  M0 <- M1 <- numeric(L)
  for (i in 1:L) {
    IYi <- as.integer(Yint == i)
    M0[i] <- mean(w0 * IYi)
    M1[i] <- mean(w1 * IYi)
  }

  pW <- pL <- 0
  for (i in 1:L) for (j in 1:L) {
    if (i > j) pW <- pW + M1[i] * M0[j]
    if (i < j) pL <- pL + M1[i] * M0[j]
  }
  pT <- max(0, 1 - pW - pL)

  psi_beta0 <- ((Av == 0) * (Rv - pi0)) * (Xmat %*% Hinv0)
  psi_beta1 <- ((Av == 1) * (Rv - pi1)) * (Xmat %*% Hinv1)

  Gamma0 <- Gamma1 <- vector("list", L)
  for (i in 1:L) {
    IYi <- as.integer(Yint == i)
    res0 <- IYi - M0[i]; res1 <- IYi - M1[i]
    Gamma0[[i]] <- colMeans((Av == 0) * (-w0 * (1 - pi0)) * res0 * Xmat)
    Gamma1[[i]] <- colMeans((Av == 1) * (-w1 * (1 - pi1)) * res1 * Xmat)
  }

  psi_M0 <- psi_M1 <- matrix(0, n, L)
  for (i in 1:L) {
    IYi <- as.integer(Yint == i)
    main0 <- w0 * (IYi - M0[i])
    main1 <- w1 * (IYi - M1[i])
    corr0 <- as.numeric(psi_beta0 %*% matrix(Gamma0[[i]], ncol = 1))
    corr1 <- as.numeric(psi_beta1 %*% matrix(Gamma1[[i]], ncol = 1))
    psi_M0[, i] <- main0 - corr0
    psi_M1[, i] <- main1 - corr1
  }

  psi_pW <- psi_pL <- numeric(n)
  for (i in 1:L) for (j in 1:L) {
    if (i > j) {
      psi_pW <- psi_pW + psi_M1[, i] * M0[j] + M1[i] * psi_M0[, j]
    } else if (i < j) {
      psi_pL <- psi_pL + psi_M1[, i] * M0[j] + M1[i] * psi_M0[, j]
    }
  }
  psi_pT <- -psi_pW - psi_pL
  Psi <- cbind(psi_pW, psi_pL, psi_pT)
  SE_out <- WinMO_logratio_SE(pW, pL, pT, Psi, eps)
  SE <- SE_out$SE
  SE_log <- SE_out$SE_log

  list(pW = pW, pL = pL, pT = pT,
       WR = if (pL > 0) pW / pL else NA_real_,
       NB = pW - pL,
       WO = (pW + 0.5 * pT) / (pL + 0.5 * pT),
       DOOR = pW + 0.5 * pT,
       SE_IF = SE,
       SE_log_IF = SE_log)
}

AIPW_WinStat_1end_IF <- function(data, A, X, Y, eps = 1e-6) {
  n <- nrow(data)
  Av <- as.integer(data[[A]])
  Xmat <- cbind(1, as.matrix(data[, X, drop = FALSE]))
  p <- ncol(Xmat)
  Yv <- data[[Y]]
  Rv <- as.integer(!is.na(Yv))
  pa <- mean(Av == 1)

  lev <- sort(unique(Yv[!is.na(Yv)]))
  L <- length(lev)
  Yint <- match(Yv, lev)
  Yint[is.na(Yint)] <- 0

  fit_pi <- function(a) {
    idx <- which(Av == a)
    Xsub <- Xmat[idx, , drop = FALSE]
    Rsub <- Rv[idx]
    fit <- glm.fit(x = Xsub, y = Rsub, family = binomial())
    b <- fit$coefficients
    eta <- drop(Xmat %*% b)
    pihat <- pmin(pmax(plogis(eta), eps), 1 - eps)
    list(pi = pihat, b = b)
  }
  f0 <- fit_pi(0); f1 <- fit_pi(1)
  pi0 <- f0$pi; pi1 <- f1$pi
  w0 <- (Av == 0) / (1 - pa) * Rv / pi0
  w1 <- (Av == 1) / pa * Rv / pi1

  fit_mu <- function(a) {
    idx <- which(Av == a & Rv == 1)
    Ysub <- factor(Yint[idx], levels = seq_len(L))
    Xsub <- as.data.frame(Xmat[idx, , drop = FALSE])

    if (L == 2) {
      fit <- glm(Ysub ~ ., data = Xsub, family = binomial())
      p <- predict(fit, newdata = as.data.frame(Xmat), type = "response")
      pred <- cbind(1 - p, p)
    } else {
      fit <- nnet::multinom(Ysub ~ ., data = Xsub, trace = FALSE)
      pred <- predict(fit, newdata = as.data.frame(Xmat), type = "probs")
      if (is.vector(pred)) pred <- matrix(pred, ncol = 1)
      pred <- matrix(pred, nrow = n)
      if (ncol(pred) < L) {
        full_pred <- matrix(0, nrow = n, ncol = L)
        colnames(full_pred) <- as.character(seq_len(L))
        seen_levels <- colnames(pred)
        full_pred[, seen_levels] <- pred
        pred <- full_pred
      }
    }
    pmin(pmax(pred, eps), 1 - eps)
  }
  mu0 <- fit_mu(0)
  mu1 <- fit_mu(1)

  psi_beta0 <- ((Av == 0) * (Rv - pi0)) * Xmat
  psi_beta1 <- ((Av == 1) * (Rv - pi1)) * Xmat

  psi_gamma0 <- matrix(0, n, p)
  psi_gamma1 <- matrix(0, n, p)
  for (i in 1:L) {
    IYi <- as.integer(Yint == i)
    mu_i0 <- mu0[, i]; mu_i1 <- mu1[, i]
    U0 <- (Av == 0) * Rv * (IYi - mu_i0) * Xmat
    U1 <- (Av == 1) * Rv * (IYi - mu_i1) * Xmat
    W0 <- (Av == 0) * Rv * mu_i0 * (1 - mu_i0)
    W1 <- (Av == 1) * Rv * mu_i1 * (1 - mu_i1)
    J0 <- crossprod(Xmat * sqrt(W0))
    J1 <- crossprod(Xmat * sqrt(W1))
    Jinv0 <- tryCatch(solve(J0), error = function(e) MASS::ginv(J0))
    Jinv1 <- tryCatch(solve(J1), error = function(e) MASS::ginv(J1))
    psi_gamma0 <- psi_gamma0 + U0 %*% Jinv0
    psi_gamma1 <- psi_gamma1 + U1 %*% Jinv1
  }

  Gamma0 <- Gamma1 <- vector("list", L)
  Delta0 <- Delta1 <- vector("list", L)
  for (i in 1:L) {
    IYi <- as.integer(Yint == i)
    Gamma0[[i]] <- colMeans((Av == 0) * (-w0 * (1 - pi0)) * (IYi - mu0[, i]) * Xmat)
    Gamma1[[i]] <- colMeans((Av == 1) * (-w1 * (1 - pi1)) * (IYi - mu1[, i]) * Xmat)
    Delta0[[i]] <- colMeans((Av == 0) * (1 - w0) * (-mu0[, i] * (1 - mu0[, i])) * Xmat)
    Delta1[[i]] <- colMeans((Av == 1) * (1 - w1) * (-mu1[, i] * (1 - mu1[, i])) * Xmat)
  }

  M0 <- M1 <- numeric(L)
  for (i in 1:L) {
    IYi <- as.integer(Yint == i)
    M0[i] <- mean(w0 * (IYi - mu0[, i]) + mu0[, i])
    M1[i] <- mean(w1 * (IYi - mu1[, i]) + mu1[, i])
  }

  psi_M0 <- psi_M1 <- matrix(0, n, L)
  for (i in 1:L) {
    IYi <- as.integer(Yint == i)
    main0 <- w0 * (IYi - mu0[, i]) + mu0[, i] - M0[i]
    main1 <- w1 * (IYi - mu1[, i]) + mu1[, i] - M1[i]
    corr0 <- as.numeric(psi_beta0 %*% Gamma0[[i]] + psi_gamma0 %*% Delta0[[i]])
    corr1 <- as.numeric(psi_beta1 %*% Gamma1[[i]] + psi_gamma1 %*% Delta1[[i]])
    psi_M0[, i] <- main0 - corr0
    psi_M1[, i] <- main1 - corr1
  }

  pW <- pL <- 0
  for (i in 1:L) for (j in 1:L) {
    if (i > j) pW <- pW + M1[i] * M0[j]
    if (i < j) pL <- pL + M1[i] * M0[j]
  }
  pT <- max(0, 1 - pW - pL)

  psi_pW <- psi_pL <- numeric(n)
  for (i in 1:L) for (j in 1:L) {
    if (i > j)
      psi_pW <- psi_pW + psi_M1[, i] * M0[j] + M1[i] * psi_M0[, j]
    else if (i < j)
      psi_pL <- psi_pL + psi_M1[, i] * M0[j] + M1[i] * psi_M0[, j]
  }
  psi_pT <- -psi_pW - psi_pL

  Psi <- cbind(psi_pW, psi_pL, psi_pT)
  SE_out <- WinMO_logratio_SE(pW, pL, pT, Psi, eps)
  SE <- SE_out$SE
  SE_log <- SE_out$SE_log

  list(pW = pW, pL = pL, pT = pT,
       WR = if (pL > 0) pW / pL else NA_real_,
       NB = pW - pL,
       WO = (pW + 0.5 * pT) / (pL + 0.5 * pT),
       DOOR = pW + 0.5 * pT,
       SE_IF = SE,
       SE_log_IF = SE_log)
}

IPW_WinStat_2end_IF <- function(data, A, X, Y1, Y2, eps = 1e-6) {
  n <- nrow(data)
  Av <- as.integer(data[[A]])
  Xmat <- cbind(1, as.matrix(data[, X, drop = FALSE]))
  p <- ncol(Xmat)
  pa <- mean(Av == 1)

  Y1v <- data[[Y1]]
  Y2v <- data[[Y2]]
  R1 <- as.integer(!is.na(Y1v))
  R2 <- as.integer(!is.na(Y2v))
  R12 <- as.integer(R1 * R2)

  lev1 <- sort(unique(Y1v[!is.na(Y1v)]))
  lev2 <- sort(unique(Y2v[!is.na(Y2v)]))
  L1 <- length(lev1)
  L2 <- length(lev2)
  Y1int <- match(Y1v, lev1)
  Y2int <- match(Y2v, lev2)
  Y1int[is.na(Y1int)] <- 0
  Y2int[is.na(Y2int)] <- 0

  fit_pi <- function(a, Rtilde) {
    idx <- which(Av == a)
    Xsub <- Xmat[idx, , drop = FALSE]
    Rsub <- Rtilde[idx]
    fit <- glm.fit(x = Xsub, y = Rsub, family = binomial())
    b <- fit$coefficients
    eta <- drop(Xmat %*% b)
    pihat <- pmin(pmax(plogis(eta), eps), 1 - eps)
    W <- as.numeric(plogis(drop(Xsub %*% b)) * (1 - plogis(drop(Xsub %*% b))))
    H <- crossprod(Xsub * sqrt(W))
    Hinv <- tryCatch(solve(H), error = function(e) MASS::ginv(H))
    list(pi = pihat, Hinv = Hinv)
  }

  f1_0 <- fit_pi(0, R1)
  f1_1 <- fit_pi(1, R1)
  f2_0 <- fit_pi(0, R12)
  f2_1 <- fit_pi(1, R12)

  pi1_0 <- f1_0$pi; pi1_1 <- f1_1$pi
  pi2_0 <- f2_0$pi; pi2_1 <- f2_1$pi
  Hinv1_0 <- f1_0$Hinv; Hinv1_1 <- f1_1$Hinv
  Hinv2_0 <- f2_0$Hinv; Hinv2_1 <- f2_1$Hinv

  w1_0 <- (Av == 0) / (1 - pa) * R1 / pi1_0
  w1_1 <- (Av == 1) / pa * R1 / pi1_1
  w2_0 <- (Av == 0) / (1 - pa) * R12 / pi2_0
  w2_1 <- (Av == 1) / pa * R12 / pi2_1

  M1_0 <- M1_1 <- numeric(L1)
  M2_0 <- M2_1 <- array(0, dim = c(L1, L2))
  for (i1 in 1:L1) {
    IY1 <- as.integer(Y1int == i1)
    M1_0[i1] <- mean(w1_0 * IY1)
    M1_1[i1] <- mean(w1_1 * IY1)
    for (i2 in 1:L2) {
      IY12 <- as.integer(Y1int == i1 & Y2int == i2)
      M2_0[i1, i2] <- mean(w2_0 * IY12)
      M2_1[i1, i2] <- mean(w2_1 * IY12)
    }
  }

  pW <- pL <- 0
  for (i1 in 1:L1) for (j1 in 1:L1) {
    if (i1 > j1) {
      pW <- pW + M1_1[i1] * M1_0[j1]
    } else if (i1 < j1) {
      pL <- pL + M1_1[i1] * M1_0[j1]
    } else {
      for (i2 in 1:L2) for (j2 in 1:L2) {
        if (i2 > j2) pW <- pW + M2_1[i1, i2] * M2_0[i1, j2]
        if (i2 < j2) pL <- pL + M2_1[i1, i2] * M2_0[i1, j2]
      }
    }
  }
  pT <- max(0, 1 - pW - pL)

  psi_beta1_0 <- ((Av == 0) * (R1 - pi1_0)) * (Xmat %*% Hinv1_0)
  psi_beta1_1 <- ((Av == 1) * (R1 - pi1_1)) * (Xmat %*% Hinv1_1)
  psi_beta2_0 <- ((Av == 0) * (R12 - pi2_0)) * (Xmat %*% Hinv2_0)
  psi_beta2_1 <- ((Av == 1) * (R12 - pi2_1)) * (Xmat %*% Hinv2_1)

  Gamma1_0 <- Gamma1_1 <- vector("list", L1)
  Gamma2_0 <- Gamma2_1 <- vector("list", L1 * L2)
  for (i1 in 1:L1) {
    IY1 <- as.integer(Y1int == i1)
    res0 <- IY1 - M1_0[i1]; res1 <- IY1 - M1_1[i1]
    Gamma1_0[[i1]] <- colMeans((Av == 0) * (-w1_0 * (1 - pi1_0)) * res0 * Xmat)
    Gamma1_1[[i1]] <- colMeans((Av == 1) * (-w1_1 * (1 - pi1_1)) * res1 * Xmat)
    for (i2 in 1:L2) {
      idx <- (i1 - 1) * L2 + i2
      IY12 <- as.integer(Y1int == i1 & Y2int == i2)
      res20 <- IY12 - M2_0[i1, i2]; res21 <- IY12 - M2_1[i1, i2]
      Gamma2_0[[idx]] <- colMeans((Av == 0) * (-w2_0 * (1 - pi2_0)) * res20 * Xmat)
      Gamma2_1[[idx]] <- colMeans((Av == 1) * (-w2_1 * (1 - pi2_1)) * res21 * Xmat)
    }
  }

  psi_M1_0 <- psi_M1_1 <- matrix(0, n, L1)
  psi_M2_0 <- psi_M2_1 <- array(0, dim = c(n, L1, L2))
  for (i1 in 1:L1) {
    IY1 <- as.integer(Y1int == i1)
    main0 <- w1_0 * (IY1 - M1_0[i1])
    main1 <- w1_1 * (IY1 - M1_1[i1])
    corr0 <- as.numeric(psi_beta1_0 %*% matrix(Gamma1_0[[i1]], ncol = 1))
    corr1 <- as.numeric(psi_beta1_1 %*% matrix(Gamma1_1[[i1]], ncol = 1))
    psi_M1_0[, i1] <- main0 - corr0
    psi_M1_1[, i1] <- main1 - corr1
    for (i2 in 1:L2) {
      idx <- (i1 - 1) * L2 + i2
      IY12 <- as.integer(Y1int == i1 & Y2int == i2)
      main20 <- w2_0 * (IY12 - M2_0[i1, i2])
      main21 <- w2_1 * (IY12 - M2_1[i1, i2])
      corr20 <- as.numeric(psi_beta2_0 %*% matrix(Gamma2_0[[idx]], ncol = 1))
      corr21 <- as.numeric(psi_beta2_1 %*% matrix(Gamma2_1[[idx]], ncol = 1))
      psi_M2_0[, i1, i2] <- main20 - corr20
      psi_M2_1[, i1, i2] <- main21 - corr21
    }
  }

  psi_pW <- psi_pL <- numeric(n)
  for (i1 in 1:L1) for (j1 in 1:L1) {
    if (i1 > j1) {
      psi_pW <- psi_pW + psi_M1_1[, i1] * M1_0[j1] + M1_1[i1] * psi_M1_0[, j1]
    } else if (i1 < j1) {
      psi_pL <- psi_pL + psi_M1_1[, i1] * M1_0[j1] + M1_1[i1] * psi_M1_0[, j1]
    } else {
      for (i2 in 1:L2) for (j2 in 1:L2) {
        if (i2 > j2) {
          psi_pW <- psi_pW + psi_M2_1[, i1, i2] * M2_0[i1, j2] + M2_1[i1, i2] * psi_M2_0[, i1, j2]
        } else if (i2 < j2) {
          psi_pL <- psi_pL + psi_M2_1[, i1, i2] * M2_0[i1, j2] + M2_1[i1, i2] * psi_M2_0[, i1, j2]
        }
      }
    }
  }

  psi_pT <- -psi_pW - psi_pL
  Psi <- cbind(psi_pW, psi_pL, psi_pT)
  SE_out <- WinMO_logratio_SE(pW, pL, pT, Psi, eps)
  SE <- SE_out$SE
  SE_log <- SE_out$SE_log

  list(pW = pW, pL = pL, pT = pT,
       WR = if (pL > 0) pW / pL else NA_real_,
       NB = pW - pL,
       WO = (pW + 0.5 * pT) / (pL + 0.5 * pT),
       DOOR = pW + 0.5 * pT,
       SE_IF = SE,
       SE_log_IF = SE_log)
}

AIPW_WinStat_2end_IF <- function(data, A, X, Y1, Y2, eps = 1e-6) {
  n <- nrow(data)
  Av <- as.integer(data[[A]])
  Xmat <- cbind(1, as.matrix(data[, X, drop = FALSE]))
  p <- ncol(Xmat)

  # outcomes
  Y1v <- data[[Y1]]; Y2v <- data[[Y2]]
  R1  <- as.integer(!is.na(Y1v))
  R12 <- as.integer(!is.na(Y1v) & !is.na(Y2v))
  pa <- mean(Av == 1)

  lev1 <- sort(unique(Y1v[!is.na(Y1v)]))
  lev2 <- sort(unique(Y2v[!is.na(Y2v)]))
  L1 <- length(lev1); L2 <- length(lev2)
  Y1int <- match(Y1v, lev1); Y2int <- match(Y2v, lev2)
  Y1int[is.na(Y1int)] <- 0; Y2int[is.na(Y2int)] <- 0

  ## fit pi for R1 and R12
  fit_pi <- function(a, Rv) {
    idx <- which(Av == a)
    Xsub <- Xmat[idx, , drop = FALSE]
    Rsub <- Rv[idx]
    fit <- glm.fit(x = Xsub, y = Rsub, family = binomial())
    b <- fit$coefficients
    eta <- drop(Xmat %*% b)
    pi_hat <- pmin(pmax(plogis(eta), eps), 1 - eps)
    list(pi = pi_hat, b = b)
  }
  f1_0 <- fit_pi(0, R1); f1_1 <- fit_pi(1, R1)
  f2_0 <- fit_pi(0, R12); f2_1 <- fit_pi(1, R12)
  pi1_0 <- f1_0$pi; pi1_1 <- f1_1$pi
  pi2_0 <- f2_0$pi; pi2_1 <- f2_1$pi

  w1_0 <- (Av == 0) / (1 - pa) * R1 / pi1_0
  w1_1 <- (Av == 1) / pa * R1 / pi1_1
  w2_0 <- (Av == 0) / (1 - pa) * R12 / pi2_0
  w2_1 <- (Av == 1) / pa * R12 / pi2_1

  ## outcome models
  fit_mu1 <- function(a) {
    idx <- which(Av == a & R1 == 1)
    Ysub <- factor(Y1int[idx], levels = seq_len(L1))
    Xsub <- as.data.frame(Xmat[idx, , drop = FALSE])
    fit <- nnet::multinom(Ysub ~ ., data = Xsub, trace = FALSE)
    pred <- predict(fit, newdata = as.data.frame(Xmat), type = "probs")
    if (is.vector(pred)) pred <- matrix(pred, ncol = 1)
    pred <- matrix(pred, nrow = n)
    if (ncol(pred) < L1) {
      full_pred <- matrix(0, n, L1)
      colnames(full_pred) <- as.character(seq_len(L1))
      full_pred[, colnames(pred)] <- pred
      pred <- full_pred
    }
    pmin(pmax(pred, eps), 1 - eps)
  }

  fit_mu2 <- function(a) {
    mu2 <- array(0, dim = c(n, L1, L2))
    for (i1 in 1:L1) {
      idx <- which(Av == a & R12 == 1 & Y1int == i1)
      if (length(idx) < 2) next
      Ysub <- factor(Y2int[idx], levels = seq_len(L2))
      Xsub <- as.data.frame(Xmat[idx, , drop = FALSE])
      fit <- nnet::multinom(Ysub ~ ., data = Xsub, trace = FALSE)
      pred <- predict(fit, newdata = as.data.frame(Xmat), type = "probs")
      if (is.vector(pred)) pred <- matrix(pred, ncol = 1)
      pred <- matrix(pred, nrow = n)
      if (ncol(pred) < L2) {
        full_pred <- matrix(0, n, L2)
        colnames(full_pred) <- as.character(seq_len(L2))
        full_pred[, colnames(pred)] <- pred
        pred <- full_pred
      }
      mu2[, i1, ] <- pmin(pmax(pred, eps), 1 - eps)
    }
    mu2
  }

  mu1_0 <- fit_mu1(0); mu1_1 <- fit_mu1(1)
  mu2_0 <- fit_mu2(0); mu2_1 <- fit_mu2(1)

  ## psi_beta, psi_gamma
  psi_beta1_0 <- ((Av == 0) * (R1 - pi1_0)) * Xmat
  psi_beta1_1 <- ((Av == 1) * (R1 - pi1_1)) * Xmat
  psi_beta2_0 <- ((Av == 0) * (R12 - pi2_0)) * Xmat
  psi_beta2_1 <- ((Av == 1) * (R12 - pi2_1)) * Xmat

  psi_gamma1_0 <- matrix(0, n, p)
  psi_gamma1_1 <- matrix(0, n, p)
  psi_gamma2_0 <- matrix(0, n, p)
  psi_gamma2_1 <- matrix(0, n, p)

  for (i1 in 1:L1) {
    IY1 <- as.integer(Y1int == i1)
    mu_i10 <- mu1_0[, i1]; mu_i11 <- mu1_1[, i1]
    U10 <- (Av == 0) * R1 * (IY1 - mu_i10) * Xmat
    U11 <- (Av == 1) * R1 * (IY1 - mu_i11) * Xmat
    W10 <- (Av == 0) * R1 * mu_i10 * (1 - mu_i10)
    W11 <- (Av == 1) * R1 * mu_i11 * (1 - mu_i11)
    J10 <- crossprod(Xmat * sqrt(W10))
    J11 <- crossprod(Xmat * sqrt(W11))
    psi_gamma1_0 <- psi_gamma1_0 + U10 %*% MASS::ginv(J10)
    psi_gamma1_1 <- psi_gamma1_1 + U11 %*% MASS::ginv(J11)
    for (i2 in 1:L2) {
      IY12 <- as.integer(Y1int == i1 & Y2int == i2)
      mu_i20 <- mu2_0[, i1, i2]; mu_i21 <- mu2_1[, i1, i2]
      U20 <- (Av == 0) * R12 * (IY12 - mu_i20) * Xmat
      U21 <- (Av == 1) * R12 * (IY12 - mu_i21) * Xmat
      W20 <- (Av == 0) * R12 * mu_i20 * (1 - mu_i20)
      W21 <- (Av == 1) * R12 * mu_i21 * (1 - mu_i21)
      J20 <- crossprod(Xmat * sqrt(W20))
      J21 <- crossprod(Xmat * sqrt(W21))
      psi_gamma2_0 <- psi_gamma2_0 + U20 %*% MASS::ginv(J20)
      psi_gamma2_1 <- psi_gamma2_1 + U21 %*% MASS::ginv(J21)
    }
  }

  ## compute M and psi_M
  M1_0 <- M1_1 <- numeric(L1)
  M2_0 <- M2_1 <- array(0, dim = c(L1, L2))
  for (i1 in 1:L1) {
    IY1 <- as.integer(Y1int == i1)
    M1_0[i1] <- mean(w1_0 * (IY1 - mu1_0[, i1]) + mu1_0[, i1])
    M1_1[i1] <- mean(w1_1 * (IY1 - mu1_1[, i1]) + mu1_1[, i1])
    for (i2 in 1:L2) {
      IY12 <- as.integer(Y1int == i1 & Y2int == i2)
      M2_0[i1, i2] <- mean(w2_0 * (IY12 - mu2_0[, i1, i2]) + mu2_0[, i1, i2])
      M2_1[i1, i2] <- mean(w2_1 * (IY12 - mu2_1[, i1, i2]) + mu2_1[, i1, i2])
    }
  }

  psi_M1_0 <- psi_M1_1 <- matrix(0, n, L1)
  psi_M2_0 <- psi_M2_1 <- array(0, dim = c(n, L1, L2))

  for (i1 in 1:L1) {
    IY1 <- as.integer(Y1int == i1)
    psi_M1_0[, i1] <- (w1_0 * (IY1 - mu1_0[, i1]) + mu1_0[, i1] - M1_0[i1]) -
      as.numeric(psi_beta1_0 %*% colMeans((Av == 0) * (-w1_0 * (1 - pi1_0)) * (IY1 - mu1_0[, i1]) * Xmat) +
                   psi_gamma1_0 %*% colMeans((Av == 0) * (1 - w1_0) * (-mu1_0[, i1] * (1 - mu1_0[, i1])) * Xmat))

    psi_M1_1[, i1] <- (w1_1 * (IY1 - mu1_1[, i1]) + mu1_1[, i1] - M1_1[i1]) -
      as.numeric(psi_beta1_1 %*% colMeans((Av == 1) * (-w1_1 * (1 - pi1_1)) * (IY1 - mu1_1[, i1]) * Xmat) +
                   psi_gamma1_1 %*% colMeans((Av == 1) * (1 - w1_1) * (-mu1_1[, i1] * (1 - mu1_1[, i1])) * Xmat))

    for (i2 in 1:L2) {
      IY12 <- as.integer(Y1int == i1 & Y2int == i2)
      psi_M2_0[, i1, i2] <- (w2_0 * (IY12 - mu2_0[, i1, i2]) + mu2_0[, i1, i2] - M2_0[i1, i2]) -
        as.numeric(psi_beta2_0 %*% colMeans((Av == 0) * (-w2_0 * (1 - pi2_0)) * (IY12 - mu2_0[, i1, i2]) * Xmat) +
                     psi_gamma2_0 %*% colMeans((Av == 0) * (1 - w2_0) * (-mu2_0[, i1, i2] * (1 - mu2_0[, i1, i2])) * Xmat))

      psi_M2_1[, i1, i2] <- (w2_1 * (IY12 - mu2_1[, i1, i2]) + mu2_1[, i1, i2] - M2_1[i1, i2]) -
        as.numeric(psi_beta2_1 %*% colMeans((Av == 1) * (-w2_1 * (1 - pi2_1)) * (IY12 - mu2_1[, i1, i2]) * Xmat) +
                     psi_gamma2_1 %*% colMeans((Av == 1) * (1 - w2_1) * (-mu2_1[, i1, i2] * (1 - mu2_1[, i1, i2])) * Xmat))
    }
  }

  ## win probability
  pW <- pL <- 0
  for (i1 in 1:L1) for (j1 in 1:L1) {
    if (i1 > j1) pW <- pW + M1_1[i1] * M1_0[j1]
    if (i1 < j1) pL <- pL + M1_1[i1] * M1_0[j1]
    if (i1 == j1) {
      for (i2 in 1:L2) for (j2 in 1:L2) {
        if (i2 > j2) pW <- pW + M2_1[i1, i2] * M2_0[j1, j2]
        if (i2 < j2) pL <- pL + M2_1[i1, i2] * M2_0[j1, j2]
      }
    }
  }
  pT <- max(0, 1 - pW - pL)

  ## IF for pW pL
  psi_pW <- psi_pL <- numeric(n)
  for (i1 in 1:L1) for (j1 in 1:L1) {
    if (i1 > j1) {
      psi_pW <- psi_pW + psi_M1_1[, i1] * M1_0[j1] + M1_1[i1] * psi_M1_0[, j1]
    } else if (i1 < j1) {
      psi_pL <- psi_pL + psi_M1_1[, i1] * M1_0[j1] + M1_1[i1] * psi_M1_0[, j1]
    } else {
      for (i2 in 1:L2) for (j2 in 1:L2) {
        if (i2 > j2)
          psi_pW <- psi_pW + psi_M2_1[, i1, i2] * M2_0[j1, j2] + M2_1[i1, i2] * psi_M2_0[, j1, j2]
        if (i2 < j2)
          psi_pL <- psi_pL + psi_M2_1[, i1, i2] * M2_0[j1, j2] + M2_1[i1, i2] * psi_M2_0[, j1, j2]
      }
    }
  }
  psi_pT <- -psi_pW - psi_pL

  Psi <- cbind(psi_pW, psi_pL, psi_pT)
  SE_out <- WinMO_logratio_SE(pW, pL, pT, Psi, eps)
  SE <- SE_out$SE
  SE_log <- SE_out$SE_log

  list(pW = pW, pL = pL, pT = pT,
       WR = if (pL > 0) pW / pL else NA_real_,
       NB = pW - pL,
       WO = (pW + 0.5*pT) / (pL + 0.5*pT),
       DOOR = pW + 0.5*pT,
       SE_IF = SE,
       SE_log_IF = SE_log)
}

IPW_WinStat_3end_IF <- function(data, A, X, Y1, Y2, Y3, eps = 1e-6) {
  n <- nrow(data)
  Av <- as.integer(data[[A]])
  Xmat <- cbind(1, as.matrix(data[, X, drop = FALSE]))
  p <- ncol(Xmat)
  pa <- mean(Av == 1)

  Y1v <- data[[Y1]]
  Y2v <- data[[Y2]]
  Y3v <- data[[Y3]]
  R1 <- as.integer(!is.na(Y1v))
  R2 <- as.integer(!is.na(Y2v))
  R3 <- as.integer(!is.na(Y3v))
  R12 <- as.integer(R1 * R2)
  R123 <- as.integer(R1 * R2 * R3)

  lev1 <- sort(unique(Y1v[!is.na(Y1v)]))
  lev2 <- sort(unique(Y2v[!is.na(Y2v)]))
  lev3 <- sort(unique(Y3v[!is.na(Y3v)]))
  L1 <- length(lev1)
  L2 <- length(lev2)
  L3 <- length(lev3)
  Y1int <- match(Y1v, lev1)
  Y2int <- match(Y2v, lev2)
  Y3int <- match(Y3v, lev3)
  Y1int[is.na(Y1int)] <- 0
  Y2int[is.na(Y2int)] <- 0
  Y3int[is.na(Y3int)] <- 0

  fit_pi <- function(a, Rtilde) {
    idx <- which(Av == a)
    Xsub <- Xmat[idx, , drop = FALSE]
    Rsub <- Rtilde[idx]
    fit <- glm.fit(x = Xsub, y = Rsub, family = binomial())
    b <- fit$coefficients
    eta <- drop(Xmat %*% b)
    pihat <- pmin(pmax(plogis(eta), eps), 1 - eps)
    W <- as.numeric(plogis(drop(Xsub %*% b)) * (1 - plogis(drop(Xsub %*% b))))
    H <- crossprod(Xsub * sqrt(W))
    Hinv <- tryCatch(solve(H), error = function(e) MASS::ginv(H))
    list(pi = pihat, Hinv = Hinv)
  }

  f1_0 <- fit_pi(0, R1)
  f1_1 <- fit_pi(1, R1)
  f2_0 <- fit_pi(0, R12)
  f2_1 <- fit_pi(1, R12)
  f3_0 <- fit_pi(0, R123)
  f3_1 <- fit_pi(1, R123)

  pi1_0 <- f1_0$pi; pi1_1 <- f1_1$pi
  pi2_0 <- f2_0$pi; pi2_1 <- f2_1$pi
  pi3_0 <- f3_0$pi; pi3_1 <- f3_1$pi

  w1_0 <- (Av == 0) / (1 - pa) * R1 / pi1_0
  w1_1 <- (Av == 1) / pa * R1 / pi1_1
  w2_0 <- (Av == 0) / (1 - pa) * R12 / pi2_0
  w2_1 <- (Av == 1) / pa * R12 / pi2_1
  w3_0 <- (Av == 0) / (1 - pa) * R123 / pi3_0
  w3_1 <- (Av == 1) / pa * R123 / pi3_1

  M1_0 <- M1_1 <- numeric(L1)
  M2_0 <- M2_1 <- array(0, dim = c(L1, L2))
  M3_0 <- M3_1 <- array(0, dim = c(L1, L2, L3))
  for (i1 in 1:L1) {
    IY1 <- as.integer(Y1int == i1)
    M1_0[i1] <- mean(w1_0 * IY1)
    M1_1[i1] <- mean(w1_1 * IY1)
    for (i2 in 1:L2) {
      IY12 <- as.integer(Y1int == i1 & Y2int == i2)
      M2_0[i1, i2] <- mean(w2_0 * IY12)
      M2_1[i1, i2] <- mean(w2_1 * IY12)
      for (i3 in 1:L3) {
        IY123 <- as.integer(Y1int == i1 & Y2int == i2 & Y3int == i3)
        M3_0[i1, i2, i3] <- mean(w3_0 * IY123)
        M3_1[i1, i2, i3] <- mean(w3_1 * IY123)
      }
    }
  }

  pW <- pL <- 0
  for (i1 in 1:L1) for (j1 in 1:L1) {
    if (i1 > j1) pW <- pW + M1_1[i1] * M1_0[j1]
    if (i1 < j1) pL <- pL + M1_1[i1] * M1_0[j1]
    if (i1 == j1) {
      for (i2 in 1:L2) for (j2 in 1:L2) {
        if (i2 > j2) pW <- pW + M2_1[i1, i2] * M2_0[i1, j2]
        if (i2 < j2) pL <- pL + M2_1[i1, i2] * M2_0[i1, j2]
        if (i2 == j2) {
          for (i3 in 1:L3) for (j3 in 1:L3) {
            if (i3 > j3) pW <- pW + M3_1[i1, i2, i3] * M3_0[i1, i2, j3]
            if (i3 < j3) pL <- pL + M3_1[i1, i2, i3] * M3_0[i1, i2, j3]
          }
        }
      }
    }
  }
  pT <- max(0, 1 - pW - pL)

  psi_beta1_0 <- ((Av == 0) * (R1 - pi1_0)) * (Xmat %*% f1_0$Hinv)
  psi_beta1_1 <- ((Av == 1) * (R1 - pi1_1)) * (Xmat %*% f1_1$Hinv)
  psi_beta2_0 <- ((Av == 0) * (R12 - pi2_0)) * (Xmat %*% f2_0$Hinv)
  psi_beta2_1 <- ((Av == 1) * (R12 - pi2_1)) * (Xmat %*% f2_1$Hinv)
  psi_beta3_0 <- ((Av == 0) * (R123 - pi3_0)) * (Xmat %*% f3_0$Hinv)
  psi_beta3_1 <- ((Av == 1) * (R123 - pi3_1)) * (Xmat %*% f3_1$Hinv)

  psi_pW <- psi_pL <- numeric(n)
  for (i1 in 1:L1) for (j1 in 1:L1) {
    if (i1 > j1) {
      psi_pW <- psi_pW + (Av == 1) * w1_1 * (Y1int == i1) * M1_0[j1] -
        (Av == 0) * w1_0 * (Y1int == j1) * M1_1[i1]
    } else if (i1 < j1) {
      psi_pL <- psi_pL + (Av == 1) * w1_1 * (Y1int == i1) * M1_0[j1] -
        (Av == 0) * w1_0 * (Y1int == j1) * M1_1[i1]
    } else {
      for (i2 in 1:L2) for (j2 in 1:L2) {
        if (i2 > j2) {
          psi_pW <- psi_pW + (Av == 1) * w2_1 * (Y1int == i1 & Y2int == i2) * M2_0[i1, j2] -
            (Av == 0) * w2_0 * (Y1int == i1 & Y2int == j2) * M2_1[i1, i2]
        } else if (i2 < j2) {
          psi_pL <- psi_pL + (Av == 1) * w2_1 * (Y1int == i1 & Y2int == i2) * M2_0[i1, j2] -
            (Av == 0) * w2_0 * (Y1int == i1 & Y2int == j2) * M2_1[i1, i2]
        } else {
          for (i3 in 1:L3) for (j3 in 1:L3) {
            if (i3 > j3) {
              psi_pW <- psi_pW + (Av == 1) * w3_1 * (Y1int == i1 & Y2int == i2 & Y3int == i3) * M3_0[i1, i2, j3] -
                (Av == 0) * w3_0 * (Y1int == i1 & Y2int == i2 & Y3int == j3) * M3_1[i1, i2, i3]
            } else if (i3 < j3) {
              psi_pL <- psi_pL + (Av == 1) * w3_1 * (Y1int == i1 & Y2int == i2 & Y3int == i3) * M3_0[i1, i2, j3] -
                (Av == 0) * w3_0 * (Y1int == i1 & Y2int == i2 & Y3int == j3) * M3_1[i1, i2, i3]
            }
          }
        }
      }
    }
  }

  psi_pT <- -psi_pW - psi_pL
  Psi <- cbind(psi_pW, psi_pL, psi_pT)
  SE_out <- WinMO_logratio_SE(pW, pL, pT, Psi, eps)
  SE <- SE_out$SE
  SE_log <- SE_out$SE_log

  list(pW = pW, pL = pL, pT = pT,
       WR = if (pL > 0) pW / pL else NA_real_,
       NB = pW - pL,
       WO = (pW + 0.5 * pT) / (pL + 0.5 * pT),
       DOOR = pW + 0.5 * pT,
       SE_IF = SE,
       SE_log_IF = SE_log)
}

AIPW_WinStat_3end_IF <- function(data, A, X, Y1, Y2, Y3, eps = 1e-6) {
  n <- nrow(data)
  Av <- as.integer(data[[A]])
  Xmat <- cbind(1, as.matrix(data[, X, drop = FALSE]))
  pa <- mean(Av == 1)

  Y1v <- data[[Y1]]; Y2v <- data[[Y2]]; Y3v <- data[[Y3]]
  R1 <- as.integer(!is.na(Y1v))
  R2 <- as.integer(!is.na(Y2v))
  R3 <- as.integer(!is.na(Y3v))
  R12 <- R1 * R2
  R123 <- R1 * R2 * R3

  lev1 <- sort(unique(Y1v[!is.na(Y1v)]))
  lev2 <- sort(unique(Y2v[!is.na(Y2v)]))
  lev3 <- sort(unique(Y3v[!is.na(Y3v)]))
  L1 <- length(lev1); L2 <- length(lev2); L3 <- length(lev3)
  Y1int <- match(Y1v, lev1); Y2int <- match(Y2v, lev2); Y3int <- match(Y3v, lev3)
  Y1int[is.na(Y1int)] <- 0; Y2int[is.na(Y2int)] <- 0; Y3int[is.na(Y3int)] <- 0

  fit_pi <- function(a, Rtilde) {
    idx <- which(Av == a)
    Xsub <- Xmat[idx, , drop = FALSE]
    Rsub <- Rtilde[idx]
    fit <- glm.fit(x = Xsub, y = Rsub, family = binomial())
    b <- fit$coefficients
    eta <- drop(Xmat %*% b)
    pi_hat <- pmin(pmax(plogis(eta), eps), 1 - eps)
    list(pi = pi_hat, b = b)
  }
  f1_0 <- fit_pi(0, R1); f1_1 <- fit_pi(1, R1)
  f2_0 <- fit_pi(0, R12); f2_1 <- fit_pi(1, R12)
  f3_0 <- fit_pi(0, R123); f3_1 <- fit_pi(1, R123)
  pi1_0 <- f1_0$pi; pi1_1 <- f1_1$pi
  pi2_0 <- f2_0$pi; pi2_1 <- f2_1$pi
  pi3_0 <- f3_0$pi; pi3_1 <- f3_1$pi

  w1_0 <- (Av == 0)/(1 - pa)*R1/pi1_0
  w1_1 <- (Av == 1)/pa*R1/pi1_1
  w2_0 <- (Av == 0)/(1 - pa)*R12/pi2_0
  w2_1 <- (Av == 1)/pa*R12/pi2_1
  w3_0 <- (Av == 0)/(1 - pa)*R123/pi3_0
  w3_1 <- (Av == 1)/pa*R123/pi3_1

  fit_mu1 <- function(a) {
    idx <- which(Av == a & R1 == 1)
    Ysub <- factor(Y1int[idx], levels = seq_len(L1))
    Xsub <- as.data.frame(Xmat[idx, , drop = FALSE])
    fit <- nnet::multinom(Ysub ~ ., data = Xsub, trace = FALSE)
    pred <- predict(fit, newdata = as.data.frame(Xmat), type = "probs")
    if (is.vector(pred)) pred <- matrix(pred, ncol = 1)
    pred <- matrix(pred, nrow = n)
    if (ncol(pred) < L1) {
      full_pred <- matrix(0, n, L1)
      colnames(full_pred) <- as.character(seq_len(L1))
      full_pred[, colnames(pred)] <- pred
      pred <- full_pred
    }
    pmin(pmax(pred, eps), 1 - eps)
  }

  fit_mu2 <- function(a) {
    idx <- which(Av == a & R12 == 1)
    Yjoint <- interaction(Y1int[idx], Y2int[idx], drop = TRUE)
    L12 <- L1 * L2
    Xsub <- as.data.frame(Xmat[idx, , drop = FALSE])
    fit <- nnet::multinom(Yjoint ~ ., data = Xsub, trace = FALSE)
    pred <- predict(fit, newdata = as.data.frame(Xmat), type = "probs")
    if (is.vector(pred)) pred <- matrix(pred, ncol = 1)
    pred <- matrix(pred, nrow = n)
    if (ncol(pred) < L12) {
      full_pred <- matrix(0, n, L12)
      colnames(full_pred) <- as.character(seq_len(L12))
      full_pred[, colnames(pred)] <- pred
      pred <- full_pred
    }
    pmin(pmax(pred, eps), 1 - eps)
  }

  fit_mu3 <- function(a) {
    idx <- which(Av == a & R123 == 1)
    Yjoint <- interaction(Y1int[idx], Y2int[idx], Y3int[idx], drop = TRUE)
    L123 <- L1 * L2 * L3
    Xsub <- as.data.frame(Xmat[idx, , drop = FALSE])
    fit <- nnet::multinom(Yjoint ~ ., data = Xsub, trace = FALSE)
    pred <- predict(fit, newdata = as.data.frame(Xmat), type = "probs")
    if (is.vector(pred)) pred <- matrix(pred, ncol = 1)
    pred <- matrix(pred, nrow = n)
    if (ncol(pred) < L123) {
      full_pred <- matrix(0, n, L123)
      colnames(full_pred) <- as.character(seq_len(L123))
      full_pred[, colnames(pred)] <- pred
      pred <- full_pred
    }
    pmin(pmax(pred, eps), 1 - eps)
  }

  mu1_0 <- fit_mu1(0); mu1_1 <- fit_mu1(1)
  mu2_0 <- fit_mu2(0); mu2_1 <- fit_mu2(1)
  mu3_0 <- fit_mu3(0); mu3_1 <- fit_mu3(1)

  M1_0 <- M1_1 <- numeric(L1)
  M2_0 <- M2_1 <- array(0, c(L1, L2))
  M3_0 <- M3_1 <- array(0, c(L1, L2, L3))

  for (i1 in 1:L1) {
    IY1 <- as.integer(Y1int == i1)
    M1_0[i1] <- mean(w1_0 * (IY1 - mu1_0[, i1]) + mu1_0[, i1])
    M1_1[i1] <- mean(w1_1 * (IY1 - mu1_1[, i1]) + mu1_1[, i1])
    for (i2 in 1:L2) {
      IY12 <- as.integer(Y1int == i1 & Y2int == i2)
      idx12 <- (i1 - 1) * L2 + i2
      M2_0[i1, i2] <- mean(w2_0 * (IY12 - mu2_0[, idx12]) + mu2_0[, idx12])
      M2_1[i1, i2] <- mean(w2_1 * (IY12 - mu2_1[, idx12]) + mu2_1[, idx12])
      for (i3 in 1:L3) {
        IY123 <- as.integer(Y1int == i1 & Y2int == i2 & Y3int == i3)
        idx123 <- (i1 - 1) * L2 * L3 + (i2 - 1) * L3 + i3
        M3_0[i1, i2, i3] <- mean(w3_0 * (IY123 - mu3_0[, idx123]) + mu3_0[, idx123])
        M3_1[i1, i2, i3] <- mean(w3_1 * (IY123 - mu3_1[, idx123]) + mu3_1[, idx123])
      }
    }
  }

  pW <- pL <- 0
  for (i1 in 1:L1) for (j1 in 1:L1) {
    if (i1 > j1) pW <- pW + M1_1[i1] * M1_0[j1]
    if (i1 < j1) pL <- pL + M1_1[i1] * M1_0[j1]
    if (i1 == j1) {
      for (i2 in 1:L2) for (j2 in 1:L2) {
        if (i2 > j2) pW <- pW + M2_1[i1, i2] * M2_0[i1, j2]
        if (i2 < j2) pL <- pL + M2_1[i1, i2] * M2_0[i1, j2]
        if (i2 == j2) {
          for (i3 in 1:L3) for (j3 in 1:L3) {
            if (i3 > j3) pW <- pW + M3_1[i1, i2, i3] * M3_0[i1, i2, j3]
            if (i3 < j3) pL <- pL + M3_1[i1, i2, i3] * M3_0[i1, i2, j3]
          }
        }
      }
    }
  }
  pT <- max(0, 1 - pW - pL)

  psi_pW <- psi_pL <- numeric(n)
  for (i1 in 1:L1) for (j1 in 1:L1) {
    if (i1 > j1) {
      psi_pW <- psi_pW + (w1_1 * ((Y1int == i1) - mu1_1[, i1]) + mu1_1[, i1] - M1_1[i1]) * M1_0[j1] +
        M1_1[i1] * (w1_0 * ((Y1int == j1) - mu1_0[, j1]) + mu1_0[, j1] - M1_0[j1])
    } else if (i1 < j1) {
      psi_pL <- psi_pL + (w1_1 * ((Y1int == i1) - mu1_1[, i1]) + mu1_1[, i1] - M1_1[i1]) * M1_0[j1] +
        M1_1[i1] * (w1_0 * ((Y1int == j1) - mu1_0[, j1]) + mu1_0[, j1] - M1_0[j1])
    } else {
      for (i2 in 1:L2) for (j2 in 1:L2) {
        if (i2 > j2) {
          psi_pW <- psi_pW + (w2_1 * ((Y1int == i1 & Y2int == i2) - mu2_1[, (i1 - 1) * L2 + i2]) +
                                mu2_1[, (i1 - 1) * L2 + i2] - M2_1[i1, i2]) * M2_0[i1, j2] +
            M2_1[i1, i2] * (w2_0 * ((Y1int == i1 & Y2int == j2) - mu2_0[, (i1 - 1) * L2 + j2]) +
                              mu2_0[, (i1 - 1) * L2 + j2] - M2_0[i1, j2])
        } else if (i2 < j2) {
          psi_pL <- psi_pL + (w2_1 * ((Y1int == i1 & Y2int == i2) - mu2_1[, (i1 - 1) * L2 + i2]) +
                                mu2_1[, (i1 - 1) * L2 + i2] - M2_1[i1, i2]) * M2_0[i1, j2] +
            M2_1[i1, i2] * (w2_0 * ((Y1int == i1 & Y2int == j2) - mu2_0[, (i1 - 1) * L2 + j2]) +
                              mu2_0[, (i1 - 1) * L2 + j2] - M2_0[i1, j2])
        } else {
          for (i3 in 1:L3) for (j3 in 1:L3) {
            if (i3 > j3) {
              psi_pW <- psi_pW + (w3_1 * ((Y1int == i1 & Y2int == i2 & Y3int == i3) - mu3_1[, (i1 - 1) * L2 * L3 + (i2 - 1) * L3 + i3]) +
                                    mu3_1[, (i1 - 1) * L2 * L3 + (i2 - 1) * L3 + i3] - M3_1[i1, i2, i3]) * M3_0[i1, i2, j3] +
                M3_1[i1, i2, i3] * (w3_0 * ((Y1int == i1 & Y2int == i2 & Y3int == j3) - mu3_0[, (i1 - 1) * L2 * L3 + (i2 - 1) * L3 + j3]) +
                                      mu3_0[, (i1 - 1) * L2 * L3 + (i2 - 1) * L3 + j3] - M3_0[i1, i2, j3])
            } else if (i3 < j3) {
              psi_pL <- psi_pL + (w3_1 * ((Y1int == i1 & Y2int == i2 & Y3int == i3) - mu3_1[, (i1 - 1) * L2 * L3 + (i2 - 1) * L3 + i3]) +
                                    mu3_1[, (i1 - 1) * L2 * L3 + (i2 - 1) * L3 + i3] - M3_1[i1, i2, i3]) * M3_0[i1, i2, j3] +
                M3_1[i1, i2, i3] * (w3_0 * ((Y1int == i1 & Y2int == i2 & Y3int == j3) - mu3_0[, (i1 - 1) * L2 * L3 + (i2 - 1) * L3 + j3]) +
                                      mu3_0[, (i1 - 1) * L2 * L3 + (i2 - 1) * L3 + j3] - M3_0[i1, i2, j3])
            }
          }
        }
      }
    }
  }
  psi_pT <- -psi_pW - psi_pL

  Psi <- cbind(psi_pW, psi_pL, psi_pT)
  SE_out <- WinMO_logratio_SE(pW, pL, pT, Psi, eps)
  SE <- SE_out$SE
  SE_log <- SE_out$SE_log

  list(pW = pW, pL = pL, pT = pT,
       WR = if (pL > 0) pW / pL else NA_real_,
       NB = pW - pL,
       WO = (pW + 0.5 * pT) / (pL + 0.5 * pT),
       DOOR = pW + 0.5 * pT,
       SE_IF = SE,
       SE_log_IF = SE_log)
}

