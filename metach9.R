#' ---
#' title: Statistical Methods for Meta-Analysis
#' author: Iris Sun
#' date: 2017-06-13
#' ---


#' # Ch9 Random Effect Models for Effect Sizes
#' Data from page 195 Table 2: Studies of the effects of open education on students'attitude toward
attdata <- data.frame(study = 1:11,
                   n1 = c(131, 40, 40, 90, 40, 79, 84, 78, 38, 38, 20),
                   n2 = c(138, 40, 40, 90, 40, 49, 45, 55, 110, 93, 23),
                   d = c(.158, -.254, .261, -.043, .649, .503, .458, .577,
                         .588, .392, -.055))
#' s^2(d): observed variability: Formula 7 in page 194
s2d <- sum((attdata$d - mean(attdata$d))^2)/(nrow(attdata)-1)
s2d
#' Conditional sampling variance: conVar
attdata$n <- (attdata$n1 + attdata$n2)
attdata$nt <- (attdata$n1 * attdata$n2)/(attdata$n)
#' How to compute the gamma
#getAnywhere(.cmicalc)
cmi <- function (mi)
{
        cmi <- ifelse(mi <= 1, NA, exp(lgamma(mi/2) - log(sqrt(mi/2)) -
                                               lgamma((mi - 1)/2)))
        return(cmi)
}
attdata$a <- (attdata$n - 2)*(cmi(attdata$n - 2))^2/(attdata$n - 4)
attdata$c1 <- 1/attdata$nt
attdata$c2 <- (attdata$a - 1)/attdata$a
# Formular 9 in page 194
attdata$conVar<- attdata$c1 + attdata$c2 * (attdata$d)^2
#' Variability in the underlying population paramters or between studies variance: formular 10 in page 194
betVar <- s2d - sum(attdata$conVar)/nrow(attdata)

#' k by k matrix B
k <- nrow(attdata)
B <- matrix(-1, k, k)
diag(B) <- (k*attdata$a - 1)
B <- 1/k * B
#' k by k diagonal matrix D
D <- matrix(0, k, k)
diag(D) <- attdata$conVar + betVar
diag(D)
e <- matrix(1, nrow = k, ncol=1)
Estdel <- .301 # calculated in the following section: dw
#' The variance of page 196 Formular 14
2*sum(diag(B) * diag(D)^2)^2 + 4 * (t(e)%*%(B)%*%(D)%*%(B)%*%(e)) * Estdel^2

#' ## Approximate method to compute the variance of effect size
attdata$vard <- (attdata$n1 + attdata$n2)/(attdata$n1 * attdata$n2) +
        attdata$d^2/(2*(attdata$n1 + attdata$n2))
#' weighted d is
dp <- sum(attdata$d/attdata$vard)/sum(1/attdata$vard)
#' The Homogeneity test
qtest <- sum((attdata$d - dp)^2/attdata$vard)
ifelse(qtest > qchisq(.95, (nrow(attdata)-1)),
                       "reject the homogeneity, random effect",
                       "fixed effect")
#'## Estimating the mean effect size
attdata$vi <- diag(D)
dw  <- sum(attdata$d/attdata$vi)/sum(1/attdata$vi)
#' variance of mean effect size
vardw <- 1/sum(1/attdata$vi)
#' Confidence interval
dw + c(-1, 1) * qnorm(.975) *sqrt(vardw)
#' ## Empirical Bayes Estimation for Random Effects Models
#' weighted d is `dw`
#' Variane of the weighted d is `vardw`
bayes_RE <- function(n1, n2, d, k , est){
        cmi <- function (mi)
        {
                cmi <- ifelse(mi <= 1, NA, exp(lgamma(mi/2) - log(sqrt(mi/2)) -
                                                       lgamma((mi - 1)/2)))
                return(cmi)
        }
        ntli <- (n1*n2)/(n1 + n2)
        ai <- (n1 + n2 -2) * (cmi(n1 + n2 -2))^2 / (n1 + n2 -4)
        vi <- 1/ntli + ((ai - 1)/ai)*d^2
        wi <- (1/(vi)^2)/sum(1/vi^2)
        est_Del <- sum(wi*d)
        s2d <- sum((d - mean(d))^2)/(k-1)
        est_varDel <- s2d - sum(vi)/k
        repDel <- list()
        repVarDel <- list()
        est_repDel <- list(est_Del)
        est_repVarDel <- list(est_varDel)
        i <- 1
        while(TRUE){
                repDel <- (d + vi*est_repDel[[i]]/est_repVarDel[[i]])/
                        (1 + vi/est_repVarDel[[i]])
                repVarDel <- 1/(vi + est_repVarDel[[i]])
                est_repDel[[i+1]] <- sum(repDel)/k
                est_repVarDel[[i+1]] <- 1/k *(sum(repDel^2 + repVarDel) -
                                                    (sum(repDel))^2/k)
                if(round(est_repDel[[i]],est) == round(est_repDel[[i+1]], est)) break
                   i <- i + 1
        }
        result <- list(est_repDel = est_repDel,
                       est_repVarDel = est_repVarDel)
        result <- do.call(cbind, result)
        return(result)
}
attRE <- bayes_RE(n1 = attdata$n1, n2 = attdata$n2, d = attdata$d, k = 11, est = 3)
attRE
