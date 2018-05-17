#' ---
#' title: Statistical Methods for Meta-Analysis
#' author: Iris Sun
#' date: 2017-06-13
#' ---

#' # Ch 5 Esimation of a Single Effect Size: Parametric and Non Parametric Methods
#' ## Effect size computation
es_com <- function(ex = ex,
                   con = con){
        ne <- length(ex)
        nc <- length(con)
        sde <- sd(ex)
        sdc <- sd(con)
        me <- mean(ex)
        mc <- mean(con)
        nst <- ne + nc -2
        ntl <- (ne * nc)/(ne + nc)
        nt <- ne + nc
        sp <- sqrt(((ne - 1) *sde^2 + (nc - 1) * sdc^2)/(nst))
        cmi <- function (mi){
        cmi <- ifelse(mi <= 1, NA,
                      exp(lgamma(mi/2) -
                                  log(sqrt(mi/2)) -
                                  lgamma((mi - 1)/2)))
                return(cmi)
        }
        ###------------###
        #Glass's g
        gg <- (me - mc)/(sdc)
        var_gg <- nst/((nst-2)*ntl) + gg^2*(nst/(nst-2) - 1/cmi(nst)^2)
        #Cohen's d
        dsp <-(me - mc)/sp
        var_dsp <- nst/((nst-2)*ntl) + dsp^2*(nst/(nst-2) - 1/cmi(nst)^2)
        #Hedges's g
        heg <- cmi(nst) * dsp
        var_heg <- (cmi(nst))^2 * var_dsp
        var_app <- 1/nst + heg^2/(2*nt)
        #MLE
        gml <- sqrt(nt/nst) * dsp
        var_gml <- (nt/nst) * var_dsp
        #Shrunken
        gs <- (nt - 4) * dsp/(nst *cmi(nst))
        var_gs <- ((nt-4)/(nst *cmi(nst)))^2 * var_dsp
        ###-----------result------------###
        result <- data.frame(gg = gg, var_gg = var_gg,
                             dsp = dsp, var_dsp = var_dsp,
                             heg = heg, var_heg = var_heg,
                             var_app = var_app,
                             gml = gml, var_gml = var_gml,
                             gs = gs, var_gs = var_gs)
        return(result)
}
ex1 <- data.frame(ex = rnorm(1000, 4, 2),
                  con = rnorm(1000, 1, 2))
ES <- es_com(ex = ex1$ex, con = ex1$con)
ES
#' ## Variance-stabilizing transformation
varsta_ci <- function(d, n1, n2){
        #page 88 by Hedges and Olkin (1983)
        a <- sqrt(4+ 2*(n1/n2) + 2*(n2/n2))
        hd <- sqrt(2) * asinh(d/a)
        etaCI <- hd + c(-1, 1) * qnorm(.975)/sqrt(n1 + n2)
        ahd <- function(x){
                a*sinh(x/sqrt(2))
        }
        ci1 <- ahd(etaCI)
        #page 89 by Kraemer (1983)
        v <- (n1 + n2)*(n1 + n2 -2)/(n1*n2)
        rho <- d/sqrt(d^2 + v)
        df <- n1 + n2 - 2
        uci <- c(-1, 1)*sqrt(qt(.975, df)^2/(df + qt(.975, df)^2))
        rhoci <- (uci - rho)/(uci*rho -1)
        ci2 <- rhoci*sqrt(v)/sqrt(1 - rhoci^2)
        ci <- list(ci1 = ci1,
                   ci2 = ci2)
        return(ci)
}
varsta_ci(d = .57, n1=10, n2=10)
varsta_ci(d = ES$heg, n1 = 1000, n2 = 1000)
#' ## Exact Confidence Intervals for Effect Sizes
exact_ci <- function(n1 = 10,
                     n2 = 10,
                     g = .6,
                     est = 3,
                     alpha = .05){
        df <- n1 + n2 - 2
        delta <- seq(-2, 4, by = .1)
        compute_delta <- function(x){
                j <- 0
                while(TRUE){
                        checki <- c()
                        for (i in 1:length(delta)){
                checki[i] <- qt(x, df = df, ncp= delta[i]*sqrt((n1*n2)/(n1+n2))) *
                        sqrt((n1 + n2)/(n1*n2))
                }
                        listi <- data.frame(checki = checki, delta = delta)
                        a <- subset(listi, checki > g)[1,]
                        b <- tail(subset(listi, checki < g),1)
                        if(round(a$delta, est) == round(b$delta, est)) break
                        j <- j+1
                        dif <- (b$delta - a$delta)/100
                        delta <- sort(seq(a$delta, b$delta, by = dif))
                }
                ci <- round(a$delta, est)
                return(ci)
        }
        ci <- c(compute_delta(alpha/2), compute_delta(1-alpha/2))
        return(ci)
}
exact_ci(n1 = 10, n2 = 10, g = .6, est = 3, alpha = .05)
#' ## Robust and Nonparametric estimators of effect size
#' For the robust estimators, the computation of a2, ... a_(n-1) is not clear.
#' see Sarhn & Greenberg, 1962, pp. 218-251: Order Statistics
#' KA: Kraemer and Andrews (KA)
#' requires that both the pretest x and posttest y scores be available for each individual in the experimental and control groups of an experiment
#' Data: Table 5 (p. 94): Systolic blood pressure data on pre and post experimental and control data from 40 hypertensives
bpdat <- data.frame(ex_x = c(134, 135, 135, 136, 145, 147, 148, 150,
                             151, 153, 153, 155, 156, 158, 162, 165,
                             167, 168, 179, 180),
                    ex_y = c(130, 131, 135, 136, 136, 138, 124, 126,
                             104, 142, 114, 166, 153, 169, 127, 130,
                             120, 121, 149, 150),
                    co_x = c(139, 140, 141, 143, 151, 152, 152, 153,
                             153, 154, 154, 159, 160, 160, 162, 163,
                             165, 169, 175, 176),
                    co_y = c(130, 131, 144, 146, 128, 156, 161, 162,
                             160, 131, 158, 166, 150, 186, 188, 153,
                             144, 147, 169, 170))
#' Functions to compute the nonparametric effect sizes
#' KA1 computed here is based on page 95. "The result obtained by Kraemer and Andrews using slightly different definitions is -1.85
comKA <- function(exPre, exPost, coPre, coPost, n){
        #----KA1
        ## control
        nc1 <- sum(ifelse(coPre < median(coPost), 1, 0))
        pc1 <- ifelse(nc1 == 0, 1/(n + 1), nc1/n)
        d_c1 <- qnorm(pc1)
        ## experimental
        ne1 <- sum(ifelse(exPre < median(exPost), 1, 0))
        pe1 <- ifelse(ne1 == 0, 1/(n + 1), ne1/n)
        d_e1 <- qnorm(pe1)
        KA1 <- d_e1 - d_c1
        #----KA2
        ## control
        nc2 <- sum(ifelse(coPost > median(coPre), 1, 0))
        pc2 <- ifelse(nc2 == 0, 1/(n + 1), nc2/n)
        d_c2 <- qnorm(pc2)
        ## experimental
        ne2 <- sum(ifelse(exPost > median(exPre), 1, 0))
        pe2 <- ifelse(ne2 == 0, 1/(n + 1), ne2/n)
        d_e2 <- qnorm(pe2)
        KA2 <- d_e2 - d_c2
        #---------------KA3
        ne3 <- sum(ifelse(exPost > exPre, 1, 0))
        ne3 <- ne3 + (sum(exPost == exPre)/2)
        pe3 <- ifelse(ne3 == 0, 1/(n + 1), ne3/n)
        d_e3 <- qnorm(pe3)
        nc3 <- sum(ifelse(coPost > coPre, 1, 0))
        nc3 <- nc3 + (sum(coPost == coPre)/2)
        pc3 <- ifelse(nc3 == 0, 1/(n + 1), nc3/n)
        d_c3 <- qnorm(pc3)
        KA3 <- d_e3 - d_c3
        #---------------KAr1
        nq1 <- sum(ifelse(coPost > (median(exPost)/median(exPre))*coPre, 1, 0))
        pq1 <- ifelse(nq1 == 0, 1/(n + 1),nq1/n)
        KAr1 <- qnorm(pq1)
        #---------------KAr1
        nq2 <- sum(ifelse(exPost > (median(coPost)/median(coPre))*exPre, 1, 0))
        pq2 <- ifelse(nq2 == 0, 1/(n + 1),nq2/n)
        KAr2 <- qnorm(pq2)
        #---------------KAind1
        nind1 <- sum(ifelse(coPost < median(exPost), 1, 0))
        pind1 <- ifelse(nind1 == 0, 1/(n+1), nind1/n)
        KAind1 <- qnorm(pind1)
        #---------------KAind2
        nind2 <- sum(ifelse(exPost > median(coPost), 1, 0))
        pind2 <- ifelse(nind2 == 0, 1/(n+1), nind2/n)
        KAind2 <- qnorm(pind2)
        #---------------results
        KA <- as.data.frame(cbind(KA1, KA2, KA3, KAr1, KAr2, KAind1, KAind2))
        return(KA)
}
comKA(bpdat$ex_x, bpdat$ex_y, bpdat$co_x, bpdat$co_y, n=nrow(bpdat))
#' if normal distribution KAr1 = KAind1; KAr2 = KAind2
#'
#' if standard deviation of experiment and control groups are same, then KAr1 = KAind1 = KAr2 = KAind2
#'
#' if in the pretest standard deviation of experimental and control groups are same, and so for the posttest. KA1 = KA2
