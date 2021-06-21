#
# Example script for generating an AIPF Cts or Dichotomous external data file
#


# Uses the supplied model and parameter functions to generate
# endpoint data for each patient at each visit
#
# n.pats.per.arm - number of patients per arm to simulate
# groups - the numbver of groups
# num.visits - number of visits to simulate
# n.arms.per.group - number of arms per group - currentlhy that must be either 1 or 2
# cov.params - list of parameters for the covariates generation function
# resp.params - list of parameters for the response generation function
# pat.covariates.function - function that will generate all the "baseline" values for each subject, 
#                           which may include the final response
# pat.response.function - function that will generate all the responses for each subject at a specified visit
#                         given the subjects basline covariates and the dose allocated
# inc.baseline - bool flag indicating if baseline is to be included
#
RunEndpointModel = function(n.pats.per.arm, num.groups, num.visits, n.arms.per.group, cov.params, resp.params, 
                            pat.cov.func, pat.resp.func, inc.baseline=FALSE)
{
  browser()
  
  patient.id = 1:(n.pats.per.arm * num.groups * n.arms.per.group)   # Generate a list of patients
  
  # Allocate patients equally to each group
  # gen group indexs 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, etc.
  # then treatment indexes 1, 2, 1, 2, 1, 2 so equal number of subjects simulated within each group
  # the proportion randomised will be controlled by the simulator, this just ensures there be some to sample from
  group.allocation = rep(1:num.groups, each=n.pats.per.arm*n.arms.per.group)
  treat.allocation = rep((3-n.arms.per.group):2, num.groups * n.pats.per.arm)
  
  # Store patient responses in a data frame
  patient.responses = data.frame(id=patient.id, group.id=group.allocation, trt=treat.allocation)
  
  # Get patient covariates
  patient.responses = pat.cov.func(patient.responses, cov.params, inc.baseline)
  
  if (inc.baseline){
    for (g in 1:num.groups){
      cat("Mean response at visit ", 0, " Group ", g, " Control:", 
          mean(patient.responses[patient.responses$group.id==g & patient.responses$trt==1, "Visit0"]), 
          " Treatment :", 
          mean(patient.responses[patient.responses$group.id==g & patient.responses$trt==2, "Visit0"]), 
          "\n")
    }
  }
  
  # Get responses for every visit
  for (v in 1:num.visits) {
    
    # Get responses from supplied model function
    resp = pat.resp.func(patient.responses, resp.params, v, num.visits, inc.baseline)
    
    # Label response columns by visit number
    visit.label = paste("Visit", sep="", toString(v))                     
    patient.responses[ ,visit.label] = resp
    
    for (g in 1:num.groups){
      cat("Mean response at visit ", v, " Group ", g, " Control:", 
          mean(patient.responses[patient.responses$group.id==g & patient.responses$trt==1, visit.label]), 
          " Treatment :",
          mean(patient.responses[patient.responses$group.id==g & patient.responses$trt==2, visit.label]), 
          "\n")
    }
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


GenFACTSSubj = function(file.name, patient.responses, baseline.inc=FALSE)
{
  
  browser()
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
  if (baseline.inc) {
    FACTS.data = rbind( data.frame(patient.id=patient.responses$id,
                            group.id=patient.responses$group.id,
                            arm.id.id=patient.responses$trt,
                            visit.id=0,
                            response=patient.responses$Visit0), FACTS.data)
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
  
  write.table(FACTS.data, file.name, sep = ",", col.names=FALSE, row.names=FALSE)
  
  # Tell user where file was saved
  cat("\n Simulated patient data written to:", file.name, "\n\n")
  
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
# Here responses are drawn from two distirbutions - the responders and non-responders
# Treatment effect differs by the proportion of subjects on that treatment that are responders.
# Thus a key parameter that is sampled for each subject is whether they are responders
# This is done in the covariates function so it is known for all subsequent draws fro that subject.
# As all other response simulation stems from that.
#
# -----------------------------------

#
# patient covariate function
# pat.resp - the patient responses data frame
# cov.params - the probability of being a responder for each dose
# inc.baseline - bool flag, if true a baselinbe respones ("visit0") is added
#
responder.cov.func <- function(pat.resp, cov.params, inc.baseline){
  # default probability of being a responder is 0
  cov <- rep(0, length(pat.resp$id))
  
  # if allocated a dose, probability of being a responder is specified
  for (i in 1:length(cov.params$pr.resp[,1])){
    ix <- which(pat.resp$group.id==i & pat.resp$trt == 1)
    if (length(ix) > 0) {
      cov[ix] <- sample(0:1, size=length(ix), replace=TRUE, prob=c(1-cov.params$pr.resp[i, 1], cov.params$pr.resp[i, 1]))
    }
    
    ix <- which(pat.resp$group.id==i & pat.resp$trt == 2)
    if (length(ix) > 0) {
      cov[ix] <- sample(0:1, size=length(ix), replace=TRUE, prob=c(1-cov.params$pr.resp[i, 2], cov.params$pr.resp[i, 2]))
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
# patient dose response function
#
# pat.resp: data frame of patients, must have a column "cov" which is their flag indicating whether they are a reponder (1) or not (0)
#
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
responder.resp.func <- function(pat.resp, resp.params, visit, max.visit, inc.baseline=FALSE){
  resp <- rep(0, length(pat.resp$id))
  
  vf <- resp.params$frac.of.resp[visit]            # just a short alias 
  sqrt.vf <- sqrt(vf) # 
  
  ix <- which(pat.resp$cov==0)
  if (length(ix) > 0) {
    resp[ix] <- rnorm(length(ix), vf *resp.params$nonrpndr.mean, sqrt.vf * resp.params$nonrpndr.sd)
  }
    
  ix <- which(pat.resp$cov==1)
  if (length(ix) > 0) {
    resp[ix] <- rnorm(length(ix), vf * resp.params$rpndr.mean, sqrt.vf * resp.params$rpndr.sd)
  }
  
  if (inc.baseline) resp <- resp + pat.resp$Visit0
  
  return(resp)
}



#
# Examples of use
#

gen.responders <- function(dir.nm, file.nm, n.pats.per.arm, num.groups, num.visits, n.arms.per.group, 
                           cov.params, resp.params, inc.baseline=FALSE) {
  browser()

  pat.resp <- RunEndpointModel(n.pats.per.arm, num.groups, num.visits, n.arms.per.group, 
                               cov.params, resp.params, 
                               responder.cov.func, responder.resp.func, inc.baseline)
  filename <- paste(dir.nm, file.nm, sep="\\")
  resp<-GenFACTSSubj(filename, pat.resp, inc.baseline)  
}

# Gen responder example:
#
# 
# Genresponder Alzheimers example 
#
alz.example <- function(){
  
  # row per group - Pr(control responder), Pr(Treatment responder)
  alz.cov.param <- list(pr.responder = rbind(c(0.2, 0.6), c(0.1, 0.3), c(0, 0.1)))

  # row per group - mean times to event
  alz.resp.params <- list(nonrpndr.mean = 0,
                          nonrpndr.sd = 5, 
                          rpndr.mean = 4,
                          rpndr.sd = 5,
                          frac.of.resp = c(0.33, 0.66, 1))

  gen.responders(".", "resp-cts.dat", 
                 n.pats.per.arm=100, num.groups=3, num.visits=3, n.arms.per.group=2,
                 alz.cov.param, alz.resp.params, inc.baseline=FALSE)
}

example2 <-function() {

	# row per group - Pr(control responder), Pr(Treatment responder)
	x2.cov.param <- list(pr.responder = rbind(c(0, 0.6), c(0.1, 0.4), c(0.15, 0.3)),
	                     base.mean=7, base.sd=4)

	# row per group - mean times to event
	x2.resp.params <- list(nonrpndr.mean = 2,
                          nonrpndr.sd = 2, 
                          rpndr.mean = 6,
                          rpndr.sd = 2,
                          frac.of.resp = c(0.5, 0.75, 0.875, 1))
	
	gen.responders(".", "ex2.dat", 
	               n.pats.per.arm=100, num.groups=3, num.visits=4, n.arms.per.group=2,
                 x2.cov.param, x2.resp.params, inc.baseline=TRUE )
}

#
# Example functions for Dichotomous model
# 

#
# We wnat to know the patient's final response for use in the visit model
# So the final response is simulate at the outset and stored as a patient covariate
#
dichot.cov.func <- function(pat.resp, cov.params, inc.baseline){
  cov <- rep(0, length(pat.resp$id))
  
  for (g in min(pat.resp$group.id):max(pat.resp$group.id)) {
    for (t in min(pat.resp$trt):max(pat.resp$trt)) {
      # ix - indexes of subjects in group g that are have treatment t
      ix <- which(pat.resp$group.id == g & pat.resp$trt == t)
      if (length(ix) > 0) {
          cov[ix] <- rbinom(n=length(ix), size=1, prob=cov.params[g,t])
      }
    }
  }
  
  pat.resp$cov <- cov
  return(pat.resp)
}

dichot.resp.func <- function(pat.resp, visit.params, visit, max.visit, inc.baseline=FALSE){
  # default probability of being a responder is 0
  visit.resp <- rep(0, length(pat.resp$id))
  
  # probability of being a responder at a visit is specified in terms of control or treatment 
  # and whether a final responder or not. We assume smae probabilities for all groups

  # control non-responders
  ix <- which(pat.resp$trt == 1 & pat.resp$cov == 0)
  if (length(ix) > 0) {
      visit.resp[ix] <- dichot.gen.visit.resp(length(ix), visit.params[visit, 1])
  }
    
  # control responders
  ix <- which(pat.resp$trt == 1 & pat.resp$cov == 1)
  if (length(ix) > 0) {
    visit.resp[ix] <- dichot.gen.visit.resp(length(ix), visit.params[visit, 2])
  }
    
  # treatment non-responders
  ix <- which(pat.resp$trt == 2 & pat.resp$cov == 0)
  if (length(ix) > 0) {
    visit.resp[ix] <- dichot.gen.visit.resp(length(ix), visit.params[visit, 3])
  }
    
  # treatment responders
  ix <- which(pat.resp$trt == 2 & pat.resp$cov == 1)
  if (length(ix) > 0) {
    visit.resp[ix] <- dichot.gen.visit.resp(length(ix), visit.params[visit, 4])
  }
  
  return(visit.resp)
}

# return n (0 or 1)s, with probability p that its a 1
dichot.gen.visit.resp <- function(n, p){
  return (sample(0:1, size=n, replace=TRUE, prob=c(1-p, p)))
}

# Dichotomous example 
#
dichot.example <- function(){
  # row per group - Pr(control responder), Pr(Treatment responder)
  dichot.resp.params <- rbind(c(0.2, 0.6), c(0.1, 0.3), c(0, 0.1))
  
  # row per visit - Pr(control non-reponder, controlrespdonder, treat non-responder, treat responder)
  dichot.visit.params <- rbind(c(0.2, 0.4, 0.3, 0.5), c(0.2, 0.5, 0.3, 0.6), c(0, 1, 0, 1))
  
  
  pat.data <- RunEndpointModel(n.pats.per.arm=100, num.groups=3, num.visits=3, 
                               n.arms.per.group=2, dichot.resp.params, dichot.visit.params,
                               dichot.cov.func, dichot.resp.func)
  facts.data <- GenFACTSSubj("resp-dichot.dat", pat.data)
}

