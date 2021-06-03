#
# model FACTS tte predictor VSRs
# 

dich.single <- function(rate, alpha, beta) {
  Y <- rexp(100000, rate)
  Y2 <- alpha + beta * Y
  Y3 <- exp(Y2) / (1 + exp(Y2))
  return(Y3)
}

dich <- function(CHR, HR, alpha, beta) {
  rates <- NULL
  for (h in HR) {
    rates <- append(rates, mean(dich.single(CHR * h, alpha, beta)))
  }
  return(rates)
}

hist.dich.single <- function(rate, alpha, beta) {
  Y3 <- dich.single(rate, alpha, beta)
  hist(Y3)
  return(mean(Y3))
}