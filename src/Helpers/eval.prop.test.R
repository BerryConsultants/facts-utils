#
# carry out a test of 2 proportions "n" times
# n - number of tests
# Nc the control arm sample size
# Nt the treatment arm sample size
# Rc the rate on control
# Rt the rate on treatment
# alpha - the alpha level for significance
#
# returns the proprotion of the n tests that are significant
#
# there are 3 versions, using respectively
# a t.test (as in FACTS)
# the R prop.test() that uses chi-squared test
# the R fisher.test()
#

eval.dich.t.test <- function(n, Nc, Nt, Rc, Rt, alpha) {
  sig <- 0
  C <- rbinom(n, Nc, Rc)
  T <- rbinom(n, Nt, Rt)
  
  for (i in 1:n) {
    if (T[i] == 0) T[i] <- 1
    if (C[i] == 0) C[i] <- 1
    res <- t.test(x=rep(c(0,1), c(Nt - T[i], T[i])), y=rep(c(0,1), c(Nc-C[i], C[i])), alternative ="greater", var.equal=TRUE)
    if (res$p.value < alpha) sig <- sig + 1
  }
  
  return (sig / n)
}

eval.prop.test <- function(n, Nc, Nt, Rc, Rt, alpha) {
  sig <- 0
  C <- rbinom(n, Nc, Rc)
  T <- rbinom(n, Nt, Rt)
  
  for (i in 1:n) {
    res <- prop.test(x=c(T[i], C[i]), n=c(Nt, Nc), alternative ="greater")
    if (res$p.value < alpha) sig <- sig + 1
  }
  
  return (sig / n)
}

eval.fisher.test <- function(n, Nc, Nt, Rc, Rt, alpha) {
  sig <- 0
  C <- rbinom(n, Nc, Rc)
  T <- rbinom(n, Nt, Rt)
  
  for (i in 1:n) {
    test.data <- matrix(c(T[i], C[i], Nt-T[i], Nc-C[i]), nrow=2)
    res <- fisher.test(test.data)#, alternative ="greater", simulate.p.value=TRUE)
    if (res$p.value < alpha) sig <- sig + 1
  }
  
  return (sig / n)
}

# to test for example the type-1 error control when testing treatment vs control
# with sample sizes 50 and 100, and a response rate of 0.1 (both groups) at alpha 0.02
# eval.fisher.test(100000, 100, 50, 0.1, 0.1, 0.02)
# to test the power if the response group has a response rate of 0.3
# eval.fisher.test(100000, 100, 50, 0.1, 0.3, 0.02)
#
# the tests are set up to test single sided with the alternate hypothesis beign that
# the treatment response is greater than the control response
