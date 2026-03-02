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
# devtools::install_github("YOUR_GITHUB/WinMO")
library(WinMO)

set.seed(20260224)

## 1) Generate toy RCT data
n <- 10000
X <- matrix(rnorm(n * 3), n, 3)
colnames(X) <- paste0("X", 1:3)

A <- rbinom(n, 1, 0.5)  # randomized treatment

# latent utilities for two ordinal endpoints
U1 <- 0.6*A + X[,1] + 4*X[,2] + rnorm(n)
U2 <- 0.8*A + 1.6*X[,1] + 5.2*X[,3] + rnorm(n)

# map to ordinal categories (3-level Y1, 4-level Y2)
Y1 <- cut(U1, breaks = c(-Inf, -0.3, 0.8, Inf),
          labels = c("Poor","OK","Good"), right = TRUE)
Y2 <- cut(U2, breaks = c(-Inf, -0.7, 0.2, 1.1, Inf),
          labels = c("A","B","C","D"), right = TRUE)

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

dat <- data.frame(
  A  = A,
  X1 = X[,1], X2 = X[,2], X3 = X[,3],
  Y1 = Y1, Y2 = Y2
)

## 3) Run WinMO (two hierarchical endpoints: Y1 primary, Y2 secondary)
Xnames <- c("X1","X2","X3")

res_ipw  <- WinMO(dat, A = "A", X = Xnames, Y = c("Y1","Y2"), method = "IPW")
res_aipw <- WinMO(dat, A = "A", X = Xnames, Y = c("Y1","Y2"), method = "AIPW")
```

Using normal approximation, we can also get the 95% confidence intervals for each estimand: 

```r
wald_ci <- function(est, se) {
  z <- qnorm(0.975)
  c(lower = est - z * se, upper = est + z * se)
}

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

num_cols <- sapply(tab_ipw, is.numeric)
tab_ipw[num_cols]  <- lapply(tab_ipw[num_cols],  round, 3)
tab_aipw[num_cols] <- lapply(tab_aipw[num_cols], round, 3)
```

Finally, we can print and view the results of this toy example: 

```r
> cat("\n--- IPW ---\n");  print(tab_ipw)
--- IPW ---
  estimand estimate    se ci_lower ci_upper
1       WR    0.915 0.037    0.843    0.987
2       WO    0.928 0.031    0.866    0.989
3       NB   -0.037 0.017   -0.070   -0.004
4     DOOR    0.481 0.008    0.465    0.498

> cat("\n--- AIPW ---\n"); print(tab_aipw)
--- AIPW ---
  estimand estimate    se ci_lower ci_upper
1       WR    0.903 0.027    0.850    0.956
2       WO    0.918 0.023    0.874    0.963
3       NB   -0.043 0.012   -0.067   -0.018
4     DOOR    0.479 0.006    0.467    0.491
```

## Reference

Yi Liu, Huiman Barnhart, Sean O'Brien, Yuliya Lokhnygina, and Roland A. Matsouaka (2026+). Estimation and Inference for Win Estimands with Multiple Ordinal Endpoints Subject to Missing Data. Unpublished pre-print. 

## Contact

Please email Yi Liu (yi.liu.biostat@gmail.com) if you have any questions or find any issues about this package. 

