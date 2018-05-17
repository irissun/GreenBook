#' ---
#' title: Ch 4 Vote Counting Methods
#' author: Iris Sun
#' date: October, 26 2017
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#'
#' Example in page 62
#' n1 = n2 = 30
setwd("~/Desktop/green_book")
dat <- read.csv("dat_tab7.csv", skip = 1)
u <- sum(ifelse(dat$g > 0, 1, 0))
k <- nrow(dat)
pc <- u/k; pc
#' ## confidence interval: use normal theory
#' 90-percent
cl1 <- pc + c(-1, 1) * qnorm(.95) *sqrt(pc *(1- pc)/k); cl1
#' when n = 30, the corresponded effect sizes interval is
n <- 30
qt(cl1[[1]], df = n*2 - 2)/sqrt(n/2)
qt(cl1[[2]], df = n*2 - 2)/sqrt(n/2)
#' ## confidence interval: use chi-square theory
calp <- qchisq(.90, 1)
b <- calp/k
cl2 <- (2*pc + b + c(-1, 1) *sqrt(b^2 + 4*b*pc*(1-pc)))/(2*(1+b)); cl2
#' when n = 30, the corresponded effect sizes interval is
n <- 30
qt(cl2[[1]], df = n*2 - 2)/sqrt(n/2)
qt(cl2[[2]], df = n*2 - 2)/sqrt(n/2)
#' Chi-square method produces a shorter confidence interval
#' We cannot state that one interval is superior. However, when confidence intervals differ in length, one would normally prefer the shorter interval.
