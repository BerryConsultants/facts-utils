#
# Two generic funtions for creating an external virtual subject response file
# can be used for FACTS Core single endpoint dichotomous or continuous
# with visits and baseline.
# 
# RunEndpointModel creates a (short wide) dataframe with one row per subject
# functions to generate covariates and respones are passed in: pat.covariates.function and pat.response.function
# the function specific parameters are passed through as lists cov.params and resp.params
#
# GenFACTSSubjectsCts takes the (short wide) dataframe and creates a long thin one, with one row per subject per visit
# this is then written out to a csv file.
#

# 
# RunEndpointModel
# Uses the supplied model and parameter functions to generate
# endpoint data for each patient at each visit
#
# n.pats.per.arm - number of patients per arm to simulate
# num.visits - number of visits to simulate
# doses - vector of dose strengths
# cov.params - list of parameters for the covariates generation function
# resp.params - list of parameters for the response generation function
# pat.covariates.function - function that will generate all the "baseline" values for each subject, 
#                           which may include the final response
# pat.response.function - function that will generate all the responses for each subject at a specified visit
#                         given the subjects basline covariates and the dose allocated
#
RunEndpointModel = function(n.pats.per.arm, num.visits, doses, cov.params, resp.params, 
                            pat.covariates.function, pat.response.function, inc.baseline=FALSE)
{
  patient.id = 1:(n.pats.per.arm * length(doses))   # Generate a list of patients
  
  # Allocate patients equally to each dose
  dose.allocation = rep(1:length(doses), each=n.pats.per.arm)
 
  # Store patient responses in a data frame
  patient.responses = data.frame(id=patient.id, dose=dose.allocation)
     
  # Setup patient covariates, and baseline (if using)
  patient.responses = pat.covariates.function(patient.responses, cov.params, inc.baseline)
  
  # Get responses for every visit
  for (v in 1:num.visits) {
    
    # Get responses from supplied model function
    resp = pat.response.function(patient.responses, resp.params, v, num.visits, inc.baseline)
    
    # Label response columns by visit number
    visit.label = paste("Visit", sep="", toString(v))  
    patient.responses[ ,visit.label] = resp
  }
  
  return (patient.responses)
}


# Write the endpoint data to a file in a format suitable for importing into FACTS
#
# Ref: FACTS documentation, Dose finding specification.pdf, p7
# 
#         #Patient ID,Dose Index, Visit ID, Efficacy
#         1,          1,          1,        0.11
#         1,          1,          2,        0.17
#         2,          2,          1,        0.27
#         2,          2,          2,        0.27
#


GenFACTSSubjCts = function(file.name, patient.responses, inc.baseline=FALSE)
{
  # Assume that the first column in the data frame of patient responses is labelled
  # "Visit1" and that the last column contains the responses for the last visit
  first.column.index = match("Visit1",colnames(patient.responses))
  last.column.index = ncol(patient.responses)
  visit.column.names = names(patient.responses)[first.column.index:last.column.index]
  n.visits = length(visit.column.names)
  
  # Set up the data frame for FACTS format
  FACTS.data = data.frame(patient.id=patient.responses$id,
                          dose.id=patient.responses$dose,
                          visit.id=1,
                          response=patient.responses$Visit1)
  
  # Set up the data frame for FACTS format - include baseline
  if (inc.baseline) {
    FACTS.data = rbind( data.frame(patient.id=patient.responses$id,
                                   dose.id=patient.responses$dose,
                                   visit.id=0,
                                   response=patient.responses$Visit0), FACTS.data)
  } 
  
  # Append responses for each successive visit
  if (n.visits > 1) for (visit.id in 2:n.visits)
  {
    visit.column.name = visit.column.names[visit.id]
    visit.data.to.add = data.frame(patient.id=patient.responses$id,
                                   dose.id=patient.responses$dose,
                                   visit.id=visit.id,
                                   response=patient.responses[ ,visit.column.name])
    FACTS.data = rbind(FACTS.data, visit.data.to.add)
  }
  
  # FACTS requires that all visits for each patient must be grouped together   
  order.by.patient = order(FACTS.data$patient.id, FACTS.data$visit.id)
  FACTS.data = FACTS.data[order.by.patient, ]
  
  # Write to file
  write.table(FACTS.data, file.name, sep = ",", col.names=FALSE, row.names=FALSE)
  
  return(FACTS.data)
}



#
# Examples of use
#


# -----------------------------------
#
# Example 1
# Responder simulation - the patient population is made up of responders & non-repsonders
#
# -----------------------------------

#
# patient covariate function
# pat.ids - list of all patient ids
# pat.doses - list of dose for each patient (should be same length at pat.ids)
# cov.params - the probability of being a responder for each dose
#
responder.cov.fn <- function(pat.resp, cov.params, inc.baseline){
  # default probability of being a responder is 0
  cov <- rep(0, length(pat.resp[,1]))
  
  # if allocated a dose, probability of being a responder is specified
  for (i in 1:length(cov.params$pr.responder)){
    ix <- which(pat.resp$dose==i)
    if (length(ix) > 0) {
      cov[ix] <- sample(0:1, size=length(ix), replace=TRUE, 
                        prob=c(1-cov.params$pr.responder[i], cov.params$pr.responder[i]))
    }
  }
  
  pat.resp$cov <- cov
  
  if (inc.baseline){
    resp <- rnorm(length(pat.resp$id), cov.params$base.mean, cov.params$base.sd)
    pat.resp$Visit0 <- resp
  }
  
  return(pat.resp)
}

#
# patient dose response function using a "responder" model
# some subjects respond, others don't - there are separate N(mean, sd) distributions for
# responders and non-responders.
# The dose-reponse differs by the responder rate for each dose.
#
# pat.resp: data frame of patients, must have a column "cov" which is their flag indicating whether they are a reponder (1) or not (0)
# resp.params : list of paramters for the response model with values
#               resp.mean     - mean change from baseline for responders
#               resp.sd       - sd of change from baseline for responders
#               nonresp.mean  - mean change from baseline for non responders
#               nonresp.sd    - sd of change from baseline for responders
#               baseline      - baseline to add in to response, set to 0 if no baseline
#
# vist.time : ppn of time of observation from baseline to final obs, visit.time == 0 indicates baseline response
# dose.resp : dose dependent response (ignored in this model - all response is via the difference in proportion of responders)
#
# returns vector or repsonses - one per row in pat.resp
# function(pat.resp, resp.params, visit, max.visit, dose.resp, inc.baseline=FALSE){

responder.resp.fn <- function(pat.resp, resp.params, visit, max.visit, inc.baseline=FALSE){
  
  resp <- rep(0, length(pat.resp$id))
    
  # response of non repsonders
  ix <- which(pat.resp$cov==0)
  if (length(ix) > 0){
    resp[ix] <- rnorm(length(ix), resp.params$nonresp.mean * resp.params$nr.frac.of.resp[visit], 
                      resp.params$nonresp.sd * sqrt(resp.params$nr.frac.of.var[visit]))
  }
    
  #response of responders
  ix <- which(pat.resp$cov==1)
  if (length(ix) > 0){
    resp[ix] <- rnorm(length(ix), resp.params$resp.mean/resp.params$r.frac.of.resp[visit], 
                      resp.params$resp.sd/sqrt(resp.params$r.frac.of.var[visit]))
  }
    
  if (inc.baseline) { resp <- resp + pat.resp$Visit0 }
  
  return(resp)
}

#
# Example use
#

gen.responders <- function(file.nm, n.pats.per.arm, inc.baseline=TRUE){

  doses <- c(1,2,3,4,5,6,7,8)

  cov.params <- list( base.mean = 10,
                      base.sd = 4,
                      pr.responder=c(0.1, 0.1, 0.15, 0.2, 0.3, 0.5, 0.6, 0.3))
  
  resp.params <- list(resp.mean = 9,
                      resp.sd = 2, 
                      nonresp.mean = 5,
                      nonresp.sd = 2,
                      # 4 post baseline visits
                      nr.frac.of.resp = c(0.2, 0.5, 0.8, 1),
                      nr.frac.of.var = c(0.2, 0.5, 0.8, 1),
                      r.frac.of.resp = c(0.4, 0.6, 0.8, 1),
                      r.frac.of.var = c(0.4, 0.6, 0.8, 1))
                      
  pat.resp <- RunEndpointModel(n.pats.per.arm, num.visits=4, doses, cov.params, resp.params, 
                               responder.cov.fn, responder.resp.fn, inc.baseline)
  GenFACTSSubjCts(file.nm, pat.resp, inc.baseline)  
}

# gen.responders("resp-responder.dat", 10000)

# -----------------------------------
#
# Example 2
# Score is an integer over a limited range - e.g. a pain score.
# 
# -----------------------------------

#
# function to apply range to a vector of simulated responses
#

apply.range <- function(resp, resp.ub, resp.lb) {
  # apply range upper bound
  ix <- which(resp > resp.ub)
  if (length(ix) > 0) {
    resp[ix] <- resp.ub
  }
  
  # apply ranger lower bound
  ix <- which(resp < resp.lb)
  if (length(ix) > 0) {
    resp[ix] <- resp.lb
  }
  
  return(resp)
}

#
# function to calculate the patetients final respone
#

range.final.resp <- function(pat.resp, dose.resp, resp.sd) {
  resp <- rep(0, length(pat.resp$id))
  
  # get dose response 
  for (d in 1:length(dose.resp)){
    ix <- which(pat.resp$dose==d)
    if (length(ix) > 0){
      resp[ix] <- rnorm(length(ix), dose.resp[d], resp.sd)
    }
  }
  
  # add baseline ('cov') to change to get final absolute
  # in order to apply scale boundaries
  # basline will be subracted at the end if change from baseline required
  resp <- pat.resp$Visit0 + resp

  # we don't apply range or truncation or rounding at this point
  pat.resp$final <- resp
  
  return(pat.resp)
}


#
# patient covariate function
# pat.ids - list of all patient ids
# pat.doses - list of dose for each patient (should be same length at pat.ids)
# cov.params -  the mean & sd of baseline
#               the lower and upper bound cut-offs
#
range.cov.fn <- function(pat.resp, cov.params, inc.baseline){
  # create basline response
  base <- rnorm(length(pat.resp[,1]), cov.params$base.mean, cov.params$base.sd)
  
  # which subjects have a base line outside the limits
  ix <- which(base < cov.params$lowerb | base > cov.params$upperb)
  
  # keep replacing subject with baseline outside limits until none left
  while (length(ix) > 0) {
    base[ix] <- rnorm(length(ix), cov.params$base.mean, cov.params$base.sd)
    ix <- which(base < cov.params$lowerb | base > cov.params$upperb)
  }
  
  pat.resp[ ,"Visit0"] = base
  pat.resp <- range.final.resp(pat.resp, cov.params$dose.resp, cov.params$resp.sd)
  
  pat.resp$Visit0 <- round(pat.resp$Visit0)
  
  return(pat.resp)
}


#
# patient dose response function - using a restricted score
#
#
# pat.resp: data frame of patients,
# resp.params : data frame (single row) of paramters for the response model with values
#               resp.lowerb - upper bound of score, all final scores above this are capped at this
#               resp.upperb - lower bound of score, all final scores below this are capped at this
#               resp.sd     - sd of change from baseline for responders
#
# vist.time : time of observation from baseline - in this version, unused.
# dose.resp : the mean response (change from baseline) for each dose
#
# returns vector or repsonses - one per row in pat.resp
#
range.resp.fn <- function(pat.resp, resp.params, visit, max.visit, inc.baseline=FALSE){
  if (visit < max.visit) {
    resp <- pat.resp$final * resp.params$frac.of.resp[visit]
    resp <- rnorm(length(resp), resp, resp.params$resp.sd * sqrt(resp.params$frac.of.var[visit]))
  } else {
    resp <- pat.resp$final
  }
  
  resp <- apply.range(resp, resp.params$resp.upperb, resp.params$resp.lowerb)
  resp <- round(resp)
  
  if (!inc.baseline) {
    resp <- resp - pat.resp$Visit0
  }
  return(resp)
}


gen.rangeresp <- function(file.nm, n.pats.per.arm, inc.baseline=FALSE){
  doses <- c(1,2,3,4,5,6,7,8,9)
  
  cov.params <- list(base.mean = 7,
                     base.sd = 2,
                     lowerb = 5,
                     upperb = 9,
                     dose.resp = c(0, -0.03, -0.13, -0.39, -0.85, -01.31, -1.57, -1.67, -1.7),
                     resp.sd = 2.12)
  
  resp.params <- list(resp.sd=2.12,
                      resp.lowerb = 0,
                      resp.upperb = 10,
                      frac.of.resp=c(1,1,1),
                      frac.of.var=c(0.8,0.4,0))
  
  
  pat.resp <- RunEndpointModel(n.pats.per.arm, num.visits = 3, doses, cov.params, resp.params,
                               range.cov.fn, range.resp.fn, inc.baseline)
  
  GenFACTSSubjCts(file.nm, pat.resp, inc.baseline)  
}

# gen.rangeresp("resp-logistic.dat", 100)



