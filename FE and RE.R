#' Data used to illustrate
attdata1 <- data.frame(study = 1:11,
                      n1 = c(131, 40, 40, 90, 40, 79, 84, 78, 38, 38, 20),
                      n2 = c(138, 40, 40, 90, 40, 49, 45, 55, 110, 93, 23),
                      d = c(.158, -.254, .261, -.043, .649, .503, .458, .577,
                            .588, .392, -.055))
attdata1
#' Function used
cmi <- function (mi)
{
        cmi <- ifelse(mi <= 1, NA, exp(lgamma(mi/2) - log(sqrt(mi/2)) -
                                               lgamma((mi - 1)/2)))
        return(cmi)
}
#' if the sample size is equal
attdata <- attdata1
attdata$n1 <- 20:30
attdata$n2 <- 20:30
mi <- attdata$n1 + attdata$n2 - 2
cmit <- cmi(mi)
nt <- attdata$n1 * attdata$n2/(attdata$n1 + attdata$n2)
attdata$g <- cmit * attdata$d
attdata$var_he <- 1/nt + (1 - (mi - 2)/(mi*cmit^2))* attdata$g^2
avg <- sum((attdata$n1 + attdata$n2)*attdata$g)/sum(attdata$n1 + attdata$n2)
attdata$var_ho <- 1/attdata$n1 + 1/attdata$n2 + avg^2/(2*(attdata$n1 + attdata$n2))
attdata$wi <- 1/attdata$var_he
attdata$wi_ho <- 1/attdata$var_ho
#' Fixed effect model
g_he <- sum(attdata$wi* attdata$g)/sum(attdata$wi); g_he
v_he <- 1/sum(attdata$wi); v_he
g_ho <- sum(attdata$wi_ho* attdata$g)/sum(attdata$wi_ho); g_ho; mean(attdata$g)
v_ho <- 1/sum(attdata$wi_ho); v_ho
#' Random effects model using HE
q_he <- sum(attdata$wi*(attdata$g - g_he)^2); q_he
c_he <- sum(attdata$wi) - sum(attdata$wi^2)/sum(attdata$wi)
tau_he <- (q_he - (nrow(attdata) - 1))/c_he; tau_he
attdata$wi_her <- 1/(attdata$var_he + tau_he)
g_her <- sum(attdata$wi_her* attdata$g)/sum(attdata$wi_her); g_her
#' Random effects model using HE
q_ho <- sum(attdata$wi_ho*(attdata$g - g_ho)^2); q_ho
c_ho <- sum(attdata$wi_ho) - sum(attdata$wi_ho^2)/sum(attdata$wi_ho)
tau_ho <- (q_ho - (nrow(attdata) - 1))/c_ho; tau_ho
attdata$wi_hor <- 1/(attdata$var_ho + tau_ho)
g_hor <- sum(attdata$wi_hor* attdata$g)/sum(attdata$wi_hor); g_hor
#' if the sample size is not equal
attdata <- attdata1
mi <- attdata$n1 + attdata$n2 - 2
cmit <- cmi(mi)
nt <- attdata$n1 * attdata$n2/(attdata$n1 + attdata$n2)
attdata$g <- cmit * attdata$d
attdata$var_he <- 1/nt + (1 - (mi - 2)/(mi*cmit^2))* attdata$g^2
avg <- sum((attdata$n1 + attdata$n2)*attdata$g)/sum(attdata$n1 + attdata$n2)
attdata$var_ho <- 1/attdata$n1 + 1/attdata$n2 + avg^2/(2*(attdata$n1 + attdata$n2))
attdata$wi <- 1/attdata$var_he
attdata$wi_ho <- 1/attdata$var_ho
#' Fixed effect model
g_he <- sum(attdata$wi* attdata$g)/sum(attdata$wi); g_he
v_he <- 1/sum(attdata$wi); v_he
g_ho <- sum(attdata$wi_ho* attdata$g)/sum(attdata$wi_ho); g_ho; mean(attdata$g)
v_ho <- 1/sum(attdata$wi_ho); v_ho
#' Random effects model using HE
q_he <- sum(attdata$wi*(attdata$g - g_he)^2); q_he
c_he <- sum(attdata$wi) - sum(attdata$wi^2)/sum(attdata$wi)
tau_he <- (q_he - (nrow(attdata) - 1))/c_he; tau_he
attdata$wi_her <- 1/(attdata$var_he + tau_he)
g_her <- sum(attdata$wi_her* attdata$g)/sum(attdata$wi_her); g_her
#' Random effects model using HE
q_ho <- sum(attdata$wi_ho*(attdata$g - g_ho)^2); q_ho
c_ho <- sum(attdata$wi_ho) - sum(attdata$wi_ho^2)/sum(attdata$wi_ho)
tau_ho <- (q_ho - (nrow(attdata) - 1))/c_ho; tau_ho
attdata$wi_hor <- 1/(attdata$var_ho + tau_ho)
g_hor <- sum(attdata$wi_hor* attdata$g)/sum(attdata$wi_hor); g_hor
