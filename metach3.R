#' # Tests of Statistical Significance of Combined Results
#' ## 1. Methods based on the uniform distribution
#' Tippett (1931) method
# function and data used
source('~/Desktop/green_book/com_test.R')
setwd("~/Desktop/green_book")
datsex <- read.csv("datb_ch2.csv", header=TRUE, skip=1)
p1 <- min(datsex$sig.level.p); p1
k <- nrow(datsex)
1 - (1 - .05)^(1/k)
#' reject H0 if the following statement is TRUE. Therefore, we reject the hypothesis of no sex difference in any study and conclude that a sex difference is exhibited in at least one study
p1 < 1 - (1 - .05)^(1/k)
com_test(datsex$sig.level.p)
#' Genearalization of Tippett's procedure: Wilkinson (1951)
#'
#' Compare the second smallest p value: r = 2; if TRUE, then reject the null hypothesis that no sex difference is exhibited in any study
orddat <- datsex[order(datsex$sig.level.p), ]
r <- 2
orddat$sig.level.p[r]
qbeta(.05, r, k - r + 1)
orddat$sig.level.p[r] < qbeta(.05, r, k - r + 1)
com_test(datsex$sig.level.p, type = "unif_W")
#' ## 2. The inverse chi-square method Fisher (1932)
#' reject H0 if p >= c
p <- -2*sum(log(datsex$sig.level.p));p
c <- qchisq(.05, 2*nrow(datsex), lower.tail = FALSE);c
p >= c
com_test(datsex$sig.level.p, type = "invF")
#' ## 3. The inverse normal method
zi <- qnorm(datsex$sig.level.p)
z <- sum(zi)/sqrt(nrow(datsex))
pnorm(z)
com_test(datsex$sig.level.p, type = "invn")
#' combined null hypothesis is rejected at the 5% level of significance
#'
#' ## 4. Logit method
Li <- sum(log(datsex$sig.level.p/(1-datsex$sig.level.p)))
Li
Li_star <- abs(Li)*sqrt((3/pi^2)*(5*nrow(datsex) + 4)/
                                (nrow(datsex)*(5*nrow(datsex) + 2)))
Li_star
pt(Li_star, df = (5*nrow(datsex) + 4), lower.tail=FALSE)
com_test(datsex$sig.level.p, type = "logitm")
#' reject the combined null hypothesis of no sex difference in any study.
