#wksf <- read.csv("weeks00009.csv",skip=1)

getcols <- function(wksf, row, colnm) {
  
  nms <- apply(array(1:20), MARGIN=1, FUN=function(x){paste(colnm,x,sep=".")})
  ix <- match(nms, names(wksf), nomatch=0)
  ix <- ix[ix > 0]
  as.vector(wksf[row, ix],"numeric")
}



error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

plot.dose.response <- function(wksf, row){
  browser()
  fit <- getcols(wksf, row, "Mean.resp")
  fit.sd <- getcols(wksf, row, "SD.Mean.resp")
  raw <- getcols(wksf, row, "Mean.Raw.Response")
  raw.se <- getcols(wksf, row, "SE.Mean.Raw.Response")
  
  min.y <- min(fit - 1.96 * fit.sd, raw - 1.96 * raw.se)
  max.y <- max(fit + 1.96 * fit.sd, raw + 1.96 * raw.se)
  
  x <- 0:(length(raw) -1)
  
  plot(x, raw, type="p", col="red", ylim=c(min.y, max.y), xlab="Arms", ylab="Mean change from baseline",
       main="Change in pain score")
  arrows(x, raw+1.96*raw.se, x, raw-1.96*raw.se, angle=90, code=3, col="red", length=0.1)
  legend(4,1.9, legend=c("Raw data mean & 95% CI", "Fitted mean", "Fitted 95% CI"),
         col=c("red", "green", "green"), lwd=c(1,2,2), lty=c(1,1,2), pch=c(1,NA,NA))
  lines(x, fit, type="l", col="green", lwd=2)
  lines(x, fit+1.96*fit.sd, col="green", lty=2, lwd=2)
  lines(x, fit-1.96*fit.sd, col="green", lty=2, lwd=2)
  op <- par(bg="white")
  
  par(op)
  
}

hist.subj <- function(wksf, row) {
  browser()
  subj <- getcols(wksf, row, "Alloc")
  comp <- getcols(wksf, row, "Complete")
  pral <- getcols(wksf, row, "Pr.Alloc.")
  
  colnames <- apply(array(0:(length(subj)-1)), MARGIN=1, FUN=function(x){paste("Arm",x,sep=" ")})
  
  # rescal probability of alloc so 1 would match the highest "subj" bar
  pral <- pral *max(subj)
  
  data <- rbind(subj,comp,pral)
  
  barplot(data, beside=TRUE, names.arg=colnames, legend.text=c("Subjects Allocated", "Subjects Complete", "Pr(allocation)"),
          main="Subject Allocation Status")
}

