#
# TTE predictive probability calculation
#

# fn - file name of a patients file to read in
#
# returns the patients file as a dataframe
#
read.pat.file <- function (fn) {
  pat <- read.csv(fn, SKIP=1, HEADER=TRUE)
  return(pat)
}

#
# calculates a log rank test p-value on a patients dataframe
# returns a single sided p-value
# the HR.lt1 flag is used to indicate whether the desired HR is gr or lt 1
#
log.rank.pv <- function(pat, HR.lt1 = TRUE) {
  res <- survdiff(Surv(Duration, Outcome) ~ Dose, data=pat)
  pval <- 1-pchisq(res$chisq,length(res$n)-1)
  HR = (res$obs[2]/res$exp[2])/(res$obs[1]/res$exp[1])
  pval <- pval / 2
  if (HR.lt1 & HR > 1) pval <- 1
  if (!HR.lt1 & HR < 1) pval <- 1
  return(pval)
}

#
# add simulated followup data to the pat dataframe
# 
# pat     the dataframe to be augmented
# lambda  the hazard rate for each dose
# fup.len the maximum followup time
# fup.individ whether the followup time is per subject or from the recruitment of the last subject
#

sim.fup.data <- function(pat, lambda, fup.len, fup.individ=FALSE){
  D <- max(pat$Dose)
  for (d in 1:D) {
    pat.ix <- which(pat$Dose == d & pat$Outcome == 0)
    obs <- rexp(length(pat.ix), lambda[d])
    pat[pat.ix,]$Duration <- pat[pat.ix,]$Duration + obs
  }
  
  if (fup.individ) {
    # patients where we observe an event in the follow up
    pat.ix <- which(pat$Duration < fup.len & pat$Outcome == 0)
    pat[pat.ix,]$Outcome <- 1
    
    # patients who are censored by the end of the follow up time
    pat.ix <- which(pat$Duration > fup.len)
    pat[pat.ix,]$Duration <- fup.len
  } else {
    # patients where we observe an event in the follow up
    end.acc <- max(pat$Date/7)
    pat.max <- (end.acc - pat$Date/7) + fup.len
    pat.ix <- which(pat$Duration < pat.max & pat$Outcome == 0)
    pat[pat.ix,]$Outcome <- 1
    
    # patients who are censored by the end of the follow up time
    pat.ix <- which(pat$Duration > pat.max)
    pat[pat.ix,]$Duration <- pat.max[pat.ix]
  }
  return(pat)
}

pred.prob <- function(pat, n, alpha, fup, lamb.mu, lamb.sd, hr.mu, hr.sd) {
  
  sig <- 0
  
  for (i in 1:n) {
    lamb <- rnorm(n=1, mean=lamb.mu, sd=lamb.sd)
    hr <- rnorm(n=1, mean=hr.mu, sd=hr.sd)
    newpat <- sim.fup.data(pat, c(lamb, lamb*hr), fup.len=fup)
    p <- log.rank.pv(newpat)
    if (p <alpha) sig <- sig + 1
  }
  
  return(sig/n)
}