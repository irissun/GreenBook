# Green book Chapter 3
setwd("~/Desktop/green_book")
datsex <- read.csv("datb_ch2.csv", header=TRUE, skip=1)
com_test <- function(pi, k=NA, alpha = .05, r = 2, type = "unif_T"){
        if(is.na(k)) k <- length(pi)
        switch(type,
              unif_T = c("Method: Uniform Tippett",
                            stats = min(pi),
                            crival = 1 - (1-alpha)^(1/k),
                         ifelse(min(pi) < 1 - (1-alpha)^(1/k),
                                "reject H0", "cannot reject H0")),
              unif_W = c("Method: Uniform Wilkinson",
                            stats = pi[order(pi)][r],
                            crival = qbeta(alpha, r, k - r + 1),
                            ifelse(pi[order(pi)][r] < qbeta(alpha, r, k - r + 1),
                                   "reject H0", "cannot reject H0")),
              invF = c("Method: Inverse chi-squared Fisher",
                           stats = -2*sum(log(pi)),
                           crival = qchisq(alpha, df = 2*k, lower.tail = FALSE),
                           ifelse(-2*sum(log(pi)) >= qchisq(alpha, df = 2*k, lower.tail = FALSE),
                                  "reject H0", "cannot reject H0") ),
              invn = c("Method: Inverse normal method",
                           stats = sum(qnorm(pi))/sqrt(k),
                           crival = qnorm(alpha/2),
                           ifelse(sum(qnorm(pi))/sqrt(k) < qnorm(alpha/2),
                                  "reject H0", "cannot reject H0" )  ),
              logitm =c("Method: The Logit Method",
                           stats = sum(log(pi/(1-pi))),
                           l_st = sqrt(.3*(5*k + 4) /(k * (5*k + 2))) *
                                abs(sum(log(pi/(1-pi)))),
                           crival = qt(.05/2, 54, lower.tail = F),
                           ifelse(sqrt(.3*(5*k + 4) /(k * (5*k + 2))) *
                                          abs(sum(log(pi/(1-pi)))) >
                                          qt(.05/2, 54, lower.tail = F),
                                  "reject H0", "cannot reject H0")
                      )
              )
}

#com_test(datsex$sig.level.p)
#com_test(datsex$sig.level.p, type = "unif_W")
#com_test(datsex$sig.level.p, type = "invF")
#com_test(datsex$sig.level.p, type = "invn")
#com_test(datsex$sig.level.p, type = "logitm")

