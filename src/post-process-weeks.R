#
# Function to post process weeks files in order to determine
# Early stopping thresholds with desired operating characteristics
#

extract.qoi.vals <- function(agg.wks.df, qoi.nm, final.qoi.nm, off.end.val=NA){
  
  # a data.frame with the scenario name and sim ID
  qoi.val.df <- unique(data.frame(Scen.ID=agg.wks.df$Scenario.ID, Sim=agg.wks.df$Sim))
  
  # the max interim number
  int.no <-unique(agg.wks.df$InterimNumber)
  max.int <- length(int.no) - 1 # ignore the final interim - "999"
  
  # now add the qoi value at each interim to the 
  for (i in 1:max.int) {
    int.df <- agg.wks.df[agg.wks.df$Int==i, c("Scenario.ID", "Sim", qoi.nm)]
    qoi.val.df <- merge(qoi.val.df, int.df, by.x=c("Scen.ID", "Sim"), by.y=c("Scenario.ID", "Sim"))
    colnames(qoi.val.df)[colnames(qoi.val.df)==qoi.nm] <- paste("Int", sep="", toString(i))
  }
  
  # now add the final eval qoi at the final analysis "interim"
  int.df <- agg.wks.df[agg.wks.df$Int==999, c("Scenario.ID", "Sim", final.qoi.nm)]
  qoi.val.df <- merge(qoi.val.df, int.df, by.x=c("Scen.ID", "Sim"), by.y=c("Scenario.ID", "Sim"))
  colnames(qoi.val.df)[colnames(qoi.val.df)==final.qoi.nm] <- "Int999"
  
  # just to make sure that all the data frames returned are in the same row order
  order.by.Scen.Sim = order(qoi.val.df$Scen.ID, qoi.val.df$Sim)
  qoi.val.df = qoi.val.df[order.by.Scen.Sim, ]
  
  return(qoi.val.df)
}

thresh.for.type1 <- function(succ.qoi, error.at.int, final.thresh) {
  
  # assume it would only be a success if it would be a success at the end
  # assume success is a p-value test
  # we are only concerned with early success of sims that would not have been successful at the end
  succ.qoi <- succ.qoi[succ.qoi$Int999 > final.thresh,]
  thresh <- rep(0.0, length(error.at.int))
  col.off = match("Int1",colnames(succ.qoi)) - 1
  
  for (i in 1:length(error.at.int)) {
    vals.at.int <- sort(succ.qoi[, col.off+i], decreasing=TRUE)
    # assume the test for success is that qoi.val is greater than threshold
    # error.at.int[i] is the fraction of errors we are happto to have at interim i
    ix <- max(1, floor(length(vals.at.int) * error.at.int[i]))
    thresh[i] <- vals.at.int[ix]
    succ.qoi <- succ.qoi[succ.qoi[,col.off+i] < thresh[i],]
  }
  return(thresh)
}


thresh.for.type2 <- function(futil.qoi, error.at.int, final.thresh) {
  
  # assume it would only be a success if it would be a success at the end
  # assume success is a p-value test
  # we are only concerned with early futility of sims that would not have been futile at the end
  # full.futil.qoi <- futil.qoi
  futil.qoi <- futil.qoi[futil.qoi$Int999 < final.thresh,]
  thresh <- rep(0.0, length(error.at.int))
  stopped <- rep(0.0, length(error.at.int))
  col.off = match("Int1",colnames(futil.qoi)) - 1
  
  for (i in 1:length(error.at.int)) {
    vals.at.int <- sort(futil.qoi[, col.off+i])
    # assume the test for futility is that qoi.val is less than threshold
    # error.at.int[i] is the fraction of errors we are happy to to have at interim i
    ix <- max(1, floor(length(vals.at.int) * error.at.int[i]))
    thresh[i] <- vals.at.int[ix]
    futil.qoi <- futil.qoi[futil.qoi[,col.off+i] >= thresh[i],]
    
    # stopped[i] <- sum(full.futil.qoi[,col.off+i] < thresh[i])
    # full.futil.qoi <- full.futil.qoi[full.futil.qoi[,col.off+i] >= thresh[i],]
  }
  return(thresh)
}

eval.thresh <- function(succ.qoi, futil.qoi, succ.thresh, futil.thresh, final.thresh) 
  {
  col.off = match("Int1",colnames(succ.qoi)) - 1
  futil.count <- rep(0.0, length(futil.thresh)+1)
  succ.count <- rep(0.0, length(succ.thresh)+1)
 
  for (i in 1:length(succ.thresh)) {
    # assume the test for futility is that qoi.val is less than threshold
    stop.for.futil <- futil.qoi[,col.off+i] < futil.thresh[i]
    futil.count[i] <- sum(stop.for.futil)
    
    # assume the test for success is that qoi.val is greater than threshold
    # also assume success & futilty thresholds don't cross!
    stop.for.succ <- succ.qoi[,col.off+i] > succ.thresh[i]
    succ.count[i] <- sum(stop.for.succ)
    
    stopped <- stop.for.futil | stop.for.succ
    
    # remove the stopped sims
    futil.qoi <- futil.qoi[ !stopped, ]
    succ.qoi <- succ.qoi[ !stopped,]
  }
  
  # and final success / futility
  futil.count[length(futil.count)] <- sum(futil.qoi$Int999 > final.thresh)
  succ.count[length(succ.count)] <- sum(futil.qoi$Int999 < final.thresh)
  
  return(list(succ.count, sum(succ.count), futil.count, sum(futil.count)))
}

#eval.thresh(my.qoi.val.succ[my.qoi.val.succ$Scen.ID==1,], 
#            my.qoi.val.fut[my.qoi.val.fut$Scen.ID==1,], 
#            succ.thresh = c(0.8, 0.75, 0.7, 0.6, 0.6), 
#            futil.thresh = c(0.05, 0.1, 0.15, 0.2, 0.3), 
#            final.thresh = 0.05)
