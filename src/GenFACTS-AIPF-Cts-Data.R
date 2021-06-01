#
# Example script for generating an AIPF Cts external data file
#


# Uses the supplied model and parameter functions to generate
# endpoint data for each patient at each visit
#
RunEndpointModel = function(n.pats.per.arm, groups, visits, n.arms.per.group, cov.params, resp.params, 
                            baseline, add.resp.to.base, baseline.params, pat.base.func,
                            pat.cov.func, pat.resp.func)
{
  browser()
  
  patient.id = 1:(n.pats.per.arm * length(groups) * n.arms.per.group)   # Generate a list of patients
  
  # Allocate patients equally to each group
  # gen group indexs 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, etc.
  # then treatment indexes 1, 2, 1, 2, 1, 2 so equal number of subjects simulated within each group
  # the proportion randomised will be controlled by the simulator, this just ensures there be some to sample from
  group.allocation = rep(1:(length(groups)), each=n.pats.per.arm*n.arms.per.group)
  treat.allocation = rep((3-n.arms.per.group):2, length(groups) * n.pats.per.arm)
  
  if (baseline) pat.base = pat.base.func(patient.id, baseline.params)
    
  # Get patient covariates
  pat.cov = pat.cov.func(patient.id, group.allocation, treat.allocation, cov.params)
  
  # Store patient responses in a data frame
  patient.responses = data.frame(id=patient.id, group.id=group.allocation, trt=treat.allocation)
  
  if (baseline) patient.responses[, "Baseline"] = pat.base
  
  # Get responses for every visit
  for (v in 1:length(visits)) {
    
    # Get responses from supplied model function
    responses = pat.resp.func(pat.cov, pat.base, resp.params, visits[v], add.resp.to.base)
    
    # Label response columns by visit number
    visit.label = paste("Visit", sep="", toString(v))                     
    patient.responses[ ,visit.label] = responses
  }
  
  return (patient.responses)
}


# Write the endpoint data to a file in a format suitable for importing into FACTS
#
# Ref: FACTS documentation, AIPF TTE User Guide p22
# 
#         #Patient ID, Group Index, Treatment Index, Time to Event
#         1,          1,          1,        8.87
#         2,          1,          2,        9.34
#         3,          2,          1,        6.78
#         4,          2,          2,        10.23
#         5,          3,          1,        9.96
#         6,          3,          2,        5.6

#


GenFACTSSubjCts = function(file.dir, file.name, patient.responses, baseline)
{
  
  n.patients = length(patient.responses$id)
  
  # Assume that the first column in the data frame of patient responses is labelled
  # "Visit1" and that the last column contains the responses for the last visit
  first.column.index = match("Visit1",colnames(patient.responses))
  last.column.index = ncol(patient.responses)
  visit.column.names = names(patient.responses)[first.column.index:last.column.index]
  n.visits = length(visit.column.names)

  FACTS.data = data.frame(patient.id=patient.responses$id,
                          group.id=patient.responses$group.id,
                          arm.id.id=patient.responses$trt,
                          visit.id=1,
                          response=patient.responses$Visit1)
  
  # Set up the data frame for FACTS format
  if (baseline) {
    FACTS.data = rbind( data.frame(patient.id=patient.responses$id,
                            group.id=patient.responses$group.id,
                            arm.id.id=patient.responses$trt,
                            visit.id=0,
                            response=patient.responses$Baseline), FACTS.data)
  } 
  
  # Append responses for each successive visit
  if (n.visits > 1) for (visit.id in 2:n.visits)
  {
    visit.column.name = visit.column.names[visit.id]
    visit.data.to.add = data.frame(patient.id=patient.responses$id,
                                   group.id=patient.responses$group.id,
                                   arm.id.id=patient.responses$trt,
                                   visit.id=visit.id,
                                   response=patient.responses[ ,visit.column.name])
    FACTS.data = rbind(FACTS.data, visit.data.to.add)
  }
  
  # FACTS requires that all visits for each patient must be grouped together   
  order.by.patient = order(FACTS.data$patient.id)
  FACTS.data = FACTS.data[order.by.patient, ]
  
  # Write to file
  file.path = paste(file.dir, sep="/", file.name)
  write.table(FACTS.data, file.path, sep = ",", col.names=FALSE, row.names=FALSE)
  
  # Tell user where file was saved
  cat("\n Simulated patient data written to:", file.path, "\n\n")
  
  return(FACTS.data)
  
  
}


#
# patient baseline function
# pat.ids - list of all patient ids
# base.params - 4 element row: mean, sd of base, lower, upper
#
baseline.fnc <- function(pat.ids, base.params){
  
  base <- rnorm(length(pat.ids), base.params$mean, base.params$sd)
  
  # repeatedly re-allocate until none fall below 'lower'
  ix <- which(base < base.params$lower)
  while (length(ix) >0){
    base[ix] <- rnorm(length(ix), base.params$mean, base.params$sd)
    ix <- which(base < base.params$lower)
  }

  # repeatedly re-allocate until none fall above 'upper'
  ix <- which(base > base.params$upper)
  while (length(ix) >0){
    base[ix] <- rnorm(length(ix), base.params$mean, base.params$sd)
    ix <- which(base > base.params$upper)
  }
  
  return(base)
}


#
# Functions for 'responder' model
#

#
# patient covariate function
# pat.ids - list of all patient ids
# pat.group - list of group for each patient (should be same length at pat.ids)
# pat.treat - list of treatment flag (1 = control, 2 = treat) for for each patient (should be same length at pat.ids)
# cov.params - 2 element row per group, the probability for each group of the subject being a responder (if on control, if on treatment)
#
responder.cov.fnc <- function(pat.ids, pat.group, pat.treat, cov.params){
  # default probability of being a responder is 0
  cov <- rep(0, length(pat.ids))
  
  # if allocated a dose, probability of being a responder is specified
  for (i in 1:length(cov.params[,1])){
    ix <- which(pat.group==i & pat.treat == 1)
    if (length(ix) > 0) {
      cov[ix] <- sample(0:1, size=length(ix), replace=TRUE, prob=c(1-cov.params[i, 1], cov.params[i, 1]))
    }
    
    ix <- which(pat.group==i & pat.treat == 2)
    if (length(ix) > 0) {
      cov[ix] <- sample(0:1, size=length(ix), replace=TRUE, prob=c(1-cov.params[i, 2], cov.params[i, 2]))
    }
  }
  
  return(cov)
}

#
# patient dose response function

# resp.params : data frame of paramters for the response model, one row for all arms, 
#               difference between arms is simply in ppn of responders 
#
#               nonrpndr.mean    - for non responders
#               nonrpndr.sd  
#               rpndr.mean  - for respondiers
#               rpndr.sd
#
# visit.fraction - the fraction of the final response and variance observed at this visit
#
# returns response for each subj - one per row in pat.resp
#
responder.resp.fnc <- function(pat.cov, pat.base, resp.params, visit.fraction, add.resp.to.base){
  resp <- rep(0, length(pat.cov))
  
  vf <- visit.fraction            # just a short alias 
  sqrt.vf <- sqrt(visit.fraction) # 
  
  ix <- which(pat.cov==0)
  if (length(ix) > 0) {
    resp[ix] <- rnorm(length(ix), vf *resp.params[1,]$nonrpndr.mean, sqrt.vf * resp.params[1,]$nonrpndr.sd)
  }
    
  ix <- which(pat.cov==1)
  if (length(ix) > 0) {
    resp[ix] <- rnorm(length(ix), vf * resp.params[1,]$rpndr.mean, sqrt.vf * resp.params[1,]$rpndr.sd)
  }
  
  if (add.resp.to.base) resp <- resp + pat.base
  
  return(resp)
}



#
# Examples of use
#

groups.3 <- c(1,2,3)
# row per group - Pr(control responder), Pr(Treatment responder)
cov.param.expl <- rbind(c(0, 0.6), c(0.1, 0.4), c(0.15, 0.3))

visits <- c(0.5, 0.75, 0.875, 1)

# row per group - mean times to event
resp.params.expl <- data.frame(nonrpndr.mean = 2,
                          nonrpndr.sd = 2, 
                          rpndr.mean = 4,
                          rpndr.sd = 2)

base.expl <- data.frame( mean=7, sd=4, lower=4, upper=9)


gen.responders <- function(dir.nm, file.nm, groups, visits, n.arms.per.group, cov.params, resp.params, 
                           baseline, add.resp.to.base, baseline.params,
                           n.pats.per.arm) {
  browser()
  pat.resp <- RunEndpointModel(n.pats.per.arm, groups, visits, n.arms.per.group, 
                               cov.params, resp.params, baseline, add.resp.to.base, baseline.params, baseline.fnc, 
                               responder.cov.fnc, responder.resp.fnc)
  
  resp<-GenFACTSSubjCts(dir.nm, file.nm, pat.resp, baseline)  
}

# Gen responder example:
#
# gen.responders("C:/FACTS/FACTS 3.3.24 Test/FACTS 3.3.24 AIPF Cts Test Script Results", 
# "resp-cts.dat", groups.3, visits, 2, cov.param.expl, resp.params.expl,  
# baseline=FALSE, add.resp.to.base=FALSE, baseline.params=NULL, 100)
#
# gen.responders("C:/FACTS/FACTS 3.3.24 Test/FACTS 3.3.24 AIPF Cts Test Script Results", 
# "resp-cts-base-cfb.dat", groups.3, visits, 2, cov.param.expl, resp.params.expl,  
# baseline=TRUE, add.resp.to.base=FALSE, base.expl, 100)
#
# gen.responders("C:/FACTS/FACTS 3.3.24 Test/FACTS 3.3.24 AIPF Cts Test Script Results", 
# "resp-cts-base-raw.dat", groups.3, visits, 2, cov.param.expl, resp.params.expl,  
# baseline=TRUE, add.resp.to.base=TRUE, base.expl, 100)

# 
# Genresponder Alzheimers example 
#
alz.example <- function(){
  alz.groups <- c(1,2,3)
  alz.visits <- c(0.33, 0.66, 1)
  # row per group - Pr(control responder), Pr(Treatment responder)
  alz.cov.param <- rbind(c(0.2, 0.6), c(0.1, 0.3), c(0, 0.1))

  # row per group - mean times to event
  alz.resp.params <- data.frame(nonrpndr.mean = 0,
                                nonrpndr.sd = 5, 
                                rpndr.mean = 4,
                                rpndr.sd = 5)

  gen.responders("~/Documents/FACTS test/FACTS 5 testing/FACTS Standard Examples", 
                 "resp-cts.dat", alz.groups, alz.visits, n.arms.per.group=2,
                 alz.cov.param, alz.resp.params,
                 baseline=FALSE, add.resp.to.base=FALSE, 
                 baseline.params=NULL, n.pats.per.arm=100)
}
