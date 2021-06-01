repeated_lm <- function(Nsims=100, Ndata=40, seed_sd=1, x_sd=2, y_sd=2) {
  intercept <- vector("numeric", Nsims)
  coeff <- vector("numeric", Nsims)
  
  for (i in 1:Nsims) {
    seed <- rnorm(Ndata, 0, seed_sd)
    xs <- seed + rnorm(Ndata, 0, x_sd)
    ys <- seed + rnorm(Ndata, 0, y_sd)
    res <- lm(ys ~ xs, list(xs, ys))
    intercept[i] <- res$coefficients[1]
    coeff[i] <- res$coefficients[2]
  }
  return(list(intercept, coeff))
}

explore_lm <- function() {
  # first with the default values
  res <- repeated_lm()
  hist(res[[1]], main="Histogram of intercept, 40 data pts")
  hist(res[[2]], main="Histogram of slope, 40 data pts")
  
  # now increase the number of data points
  res <- repeated_lm(Ndata=1000)
  hist(res[[1]], main="Histogram of intercept, 1000 data pts")
  hist(res[[2]], main="Histogram of slope, 1000 data pts")
  
  # now increase the variability in the seed c/w the observation error
  res <- repeated_lm(seed_sd=10)
  hist(res[[1]], main="Histogram of intercept, 40 data pts, seed sd 10")
  hist(res[[2]], main="Histogram of slope, 40 data pts, seed sd 10")
}