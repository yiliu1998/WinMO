# WinMO

`WinMO` is an R package that implements inverse probability weighting (IPW) and augmented IPW (AIPW) estimators for estimation and inference of commonly used win estimands (win ratio [WR], win odds [WO], net benefit [NB], and desirability of the ordinal outcome ranking [DOOR]) defined by hierarchical **ordinal** endpoints subject to **missing data** (e.g., due to noncompliance and dropout) in **randomized** clinical trials. The package supports analyses with up to three ordinal endpoints, with an arbitrary number of levels for each endpoint. For both IPW and AIPW estimators, we develope the influence-function-based asymptotic variance estimators for inference. 

## Installation
To install the latest version of the R package from GitHub, please run following commands in R:

```r
if (!require("devtools"))
install.packages("devtools")
devtools::install_github("yiliu1998/WinMO")
```

## Usage

Below is a minimal reproducible example that simulates a randomized clinical trial with two hierarchical ordinal endpoints (`Y1`, `Y2`), induces missingness depending on covariates, and then estimates win estimands using IPW and AIPW via `WinMO`.

```r
library(WinMO)
set.seed(20260224)

## 1) Generate toy RCT data
n <- 10000
X <- matrix(rnorm(n * 3), n, 3) # baseline covariates
colnames(X) <- paste0("X", 1:3)
A <- rbinom(n, 1, 0.5)  # randomized treatment

# latent utilities for two ordinal endpoints
U1 <- 0.6*A + X[,1] + 4*X[,2] + rnorm(n)
U2 <- 0.8*A + 1.6*X[,1] + 5.2*X[,3] + rnorm(n)

# map to ordinal categories (3-level Y1, 4-level Y2)
Y1 <- cut(U1, breaks = c(-Inf, -0.3, 0.8, Inf), labels = c("Poor","OK","Good"), right = TRUE)
Y2 <- cut(U2, breaks = c(-Inf, -0.7, 0.2, 1.1, Inf), labels = c("A","B","C","D"), right = TRUE)
Y1 <- as.character(Y1)
Y2 <- as.character(Y2)

## 2) Induce missingness depending only on covariates (non-monotone)
logit <- function(z) 1/(1+exp(-z))
pR1 <- logit(0.3 - 0.5*X[,1] + 0.4*X[,2] - 0.2*X[,3])
pR2 <- logit(0.2 + 0.3*X[,1] - 0.4*X[,2] + 0.5*X[,3])
R1 <- rbinom(n, 1, pR1)
R2 <- rbinom(n, 1, pR2)
Y1[R1 == 0] <- NA
Y2[R2 == 0] <- NA

# observed data frame
dat <- data.frame(A = A, X1 = X[,1], X2 = X[,2], X3 = X[,3], Y1 = Y1, Y2 = Y2)

## 3) Run WinMO function (two hierarchical endpoints: Y1 primary, Y2 secondary)

Xnames <- c("X1","X2","X3")
res_ipw  <- WinMO(dat, A = "A", X = Xnames, Y = c("Y1","Y2"), method = "IPW")
res_aipw <- WinMO(dat, A = "A", X = Xnames, Y = c("Y1","Y2"), method = "AIPW")

## 4）Summarize results
summarize_winmo <- function(res) {
  est <- c(WR = res$WR, WO = res$WO, NB = res$NB, DOOR = res$DOOR)
  se  <- as.numeric(res$SE_IF[c("WR","WO","NB","DOOR")])
  ci  <- t(mapply(wald_ci, est, se))
  out <- data.frame(
    estimand = names(est),
    estimate = as.numeric(est),
    se       = se,
    ci_lower = ci[, "lower"],
    ci_upper = ci[, "upper"],
    row.names = NULL
  )
  out
}
tab_ipw  <- summarize_winmo(res_ipw)
tab_aipw <- summarize_winmo(res_aipw)
```
View the results: 

```r
tab_ipw
```
This gives
```r
  estimand    estimate          SE       lower       upper
1       WR  0.91471376 0.036836174  0.84529192  0.98983704
2       NB -0.03738413 0.016865417 -0.07043974 -0.00432852
3       WO  0.92792616 0.031343528  0.86848335  0.99143749
4     DOOR  0.48130793 0.008432709  0.46478013  0.49783574

```

```r
tab_aipw
```
This gives
  
```r
  estimand    estimate          SE       lower      upper
1       WR  0.90308152 0.027042103  0.85160531  0.9576693
2       NB -0.04256247 0.012371842 -0.06681083 -0.0183141
3       WO  0.91835028 0.022764610  0.87479900  0.9640697
4     DOOR  0.47871877 0.006185921  0.46659458  0.4908430
```

## Reference

Yi Liu, Huiman Barnhart, Sean O'Brien, Yuliya Lokhnygina, and Roland A. Matsouaka (2026+). Estimation and Inference for Win Estimands with Multiple Ordinal Endpoints Subject to Missing Data. Unpublished pre-print. 

## Contact

Please email Yi Liu (yi.liu.biostat@gmail.com) if you have any questions or find any issues about this package. 

