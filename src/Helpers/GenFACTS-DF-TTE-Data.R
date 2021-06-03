#
# Example script for generating an Dose Finding TTE external data file
#


# Uses the supplied model and parameter functions to generate
# endpoint data for each patient at each visit
#
RunEndpointModel = function(n.pats.per.arm, doses, cov.params, resp.params, 
                            pat.cov.func, pat.resp.func)
{
  patient.id = 1:(n.pats.per.arm * length(doses))   # Generate a list of patients
  
  # Allocate patients equally to each dose
  dose.allocation = rep(1:length(doses), times=n.pats.per.arm)
  
  # Get patient covariates
  pat.cov = pat.cov.func(patient.id, dose.allocation, cov.params)
  
  # Store patient responses in a data frame
  patient.responses = data.frame(id=patient.id, dose=dose.allocation)
  
  # Get responses from supplied model function
  responses = pat.resp.func(dose.allocation, pat.cov, resp.params)
    
  patient.responses[ , "TTE"] = responses
  
  return (patient.responses)
}


# Write the endpoint data to a file in a format suitable for importing into FACTS
#
# Ref: FACTS documentation, AIPF TTE User Guide p22
# 
#         #Patient ID, Treatment Index, Time to Event
#         1,                  1,        8.87
#         2,                  2,        9.34
#         3,                  1,        6.78
#         4,                  2,        10.23
#         5,                  1,        9.96
#         6,                  2,        5.6

#


GenFACTSSubjCts = function(file.dir, file.name, patient.responses)
{
  
  # the TTE format is much more straight forward to write than cts or dichotomous
  # no re-arranging of data is required
    
  # Write to file
  file.path = paste(file.dir, sep="/", file.name)
  write.table(patient.responses, file.path, sep = ",", col.names=FALSE, row.names=FALSE)
  
  # Tell user where file was saved
  cat("\n Simulated patient data written to:", file.path, "\n\n")
  
  return(patient.responses)
}

#
# Funcions for 'responder' model
#

#
# patient covariate function
# pat.ids - list of all patient ids
# pat.dose - list of dose index (1 = control, 2 = dose1, ...) for for each patient (should be same length at pat.ids)
# cov.params - list per dose, the probability for each dose of the subject being a responder
#
responder.cov.fn <- function(pat.ids, pat.dose, cov.params){
  # default probability of being a responder is 0
  cov <- rep(0, length(pat.ids))
  
  # if allocated a dose, probability of being a responder is specified
  for (i in 1:length(cov.params)){
    ix <- which(pat.dose==i)
    if (length(ix) > 0) {
      cov[ix] <- sample(0:1, size=length(ix), replace=TRUE, prob=c(1-cov.params[i], cov.params[i]))
    }
  }
  
  return(cov)
}

#
# patient dose response function

# resp.params : matrix of paramters for the response model, row per dose, 
#               two entries - mean time to event for (non-responders, responders)
#
# returns tte for each subj - one per row in pat.trt
#
responder.resp.fn <- function(pat.trt, pat.cov, resp.params){
  resp <- rep(0, length(pat.trt))
  
  trt <- unique(pat.trt)
  
  for (i in 1:length(trt)){
    ix <- which(pat.trt==i  & pat.cov==0)
    if (length(ix) > 0) {
      resp[ix] <- rexp(length(ix), 1 / resp.params[i,1])
    }
    
    ix <- which(pat.trt==i & pat.cov==1)
    if (length(ix) > 0) {
      resp[ix] <- rexp(length(ix), 1 / resp.params[i,2])
    }
  }
  
  return(resp)
}



#
# Examples of use
#

doses.4 <- c(1,2,3,4)
# entry per dose - Pr(responder)
cov.param.null <- c(0.1, 0.1, 0.1, 0.2)
cov.param.good <- c(0.1, 0.2, 0.25, 0.2)
cov.param.weak <- c(0.1, 0.15, 0.2, 0.2)

# row per dose - time to event non-responders, repsonders
resp.params.null <- matrix( c(20, 24, 20, 24, 20, 24, 24, 30), byrow=TRUE, nrow=4, ncol=2)
resp.params.D1.good <- matrix( c(20, 24, 24, 32, 22, 28, 24, 30), byrow=TRUE, nrow=4, ncol=2)
resp.params.D2.good <- matrix( c(20, 24, 22, 28, 24, 32, 24, 30), byrow=TRUE, nrow=4, ncol=2)
resp.params.D1.weak <- matrix( c(20, 24, 22, 28, 20, 24, 24, 30), byrow=TRUE, nrow=4, ncol=2)
resp.params.D2.weak <- matrix( c(20, 24, 20, 24, 22, 28, 24, 30), byrow=TRUE, nrow=4, ncol=2)

gen.responders <- function(dir.nm, file.nm, doses, cov.params, resp.params, n.pats.per.arm){

  pat.resp <- RunEndpointModel(n.pats.per.arm, doses,
                               cov.params, resp.params, 
                               responder.cov.fn, responder.resp.fn)
  
  resp<-GenFACTSSubjCts(dir.nm, file.nm, pat.resp)  
}

# Gen responder example:
#
# gen.responders("C:/FACTS/FACTS 3.5.4B test", 
# "df-tte-expl.dat", doses.4, cov.param.good, resp.params.D2.good,  1000)

# Gen directory of external files
#

gen.example.dir <- function(){
  gen.responders("C:/FACTS/FACTS 3.5.4B test/tte ext dir", "df-tte-D2.good.dat", 
                 doses.4, cov.param.good, resp.params.D2.good,  1000)
  gen.responders("C:/FACTS/FACTS 3.5.4B test/tte ext dir", "df-tte-D1.good.dat", 
                 doses.4, cov.param.good, resp.params.D1.good,  1000)
  gen.responders("C:/FACTS/FACTS 3.5.4B test/tte ext dir", "df-tte-null.good.dat", 
                 doses.4, cov.param.good, resp.params.null,  1000)
  
  gen.responders("C:/FACTS/FACTS 3.5.4B test/tte ext dir", "df-tte-D2.null.dat", 
                 doses.4, cov.param.null, resp.params.D2.good,  1000)
  gen.responders("C:/FACTS/FACTS 3.5.4B test/tte ext dir", "df-tte-D1.null.dat", 
                 doses.4, cov.param.null, resp.params.D1.good,  1000)
  gen.responders("C:/FACTS/FACTS 3.5.4B test/tte ext dir", "df-tte-null.null.dat", 
                 doses.4, cov.param.null, resp.params.null,  1000)  
}


# For trial execution mode
#
# gen.responders("C:/FACTS/FACTS 3.5.4B test/Analysis", 
# "df-tte-D2.good.dat", doses.e, cov.param.good, resp.params.D2.good, 30)
#
