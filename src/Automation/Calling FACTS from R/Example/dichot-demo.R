## Set up Folders and Paths

# This is the directory where the parameter file and patient data must be located
# It will be where the MCMCM files are written

setwd("Z:/FACTS test/FACTS 6 Training/FACTS R interface/Example")


# This must be the location of the factR.R file
FactR.src = "../factR.R"

# This must be the location of the executable files
Exec.dir = "../WindowsExecutables"

# Load runFACTS
source(FactR.src)


# Test to check its working
# Copy an example patients file from the simulations results to this folder before running.
system.time(runFACTS(engine='dichot', data.file = 'patients00001.csv', param.file = 'bin1_e.param',
               mcmc.file.num = 1, rng.seed = 1, exec.path = Exec.dir))





# generates a data frame that can be used to drive a FACTS analysis
# dichotomous endpoint
# no visits
#
# n.per.arm: int, the number of subjects to be simulated for each arm
# rates: int[], the response rate to be simulated for each arm
#               the length if rates defines the number of arms
#
# returns a dataframe with n.per.arm * length(rates) simulated subjects
#

genBinaryData <- function(nPerArm, rates) {
  patientID <- 1:(nPerArm * length(rates)) # Generate a list of patients
  region <- rep(1, nPerArm * length(rates)) # all patients come from region 1
  date <- 1:(nPerArm * length(rates)) # Generate a list of enrolment dates - here simply one per day
  doseAlloc <- rep(1:length(rates), each=nPerArm) # Allocate patients equally to each dose
  lastVisit <- rep(1, nPerArm * length(rates)) # all patients have last visit data
  dropout <- rep(0, nPerArm * length(rates)) # no patients  have dropped out
  baseline <- rep(-9999, nPerArm * length(rates)) # not simulating baseline
  visit1 <- rep(0, nPerArm * length(rates)) # create the outcome vector
  
  # get responses for each dose
  for (d in 1:length(rates)){
    # get indices of patients on dose d
    ix <- which(doseAlloc==d)
    
    # assign them a final response based on the rate to simulate for dose d
    if (length(ix) > 0) {
      visit1[ix] <- sample(c(0,1), size=length(ix), replace=TRUE,
                           prob=c(1-rates[d], rates[d]))
    }
  } 
  
  dat <- data.frame(SubjectID=patientID, Region=region, Date=date,
                    Dose=doseAlloc, LastVisit=lastVisit,
                    Dropout=dropout, Baseline=baseline, Visit1=visit1, 
                    row.names=NULL)
  return(dat)
}

########### Toy Example Trial Sim ##########
### Constants

DATAFILE = "patients.csv"
MCMCFILE = "mcmc00000.csv"

#
# function to simulate an example data set with dichotomous endpoint
#
# nSims - the number of sims to run
# nBurnin - the number of MCMC smaples to discard 
#           (the number of MCMC samples is specified in the parameter file)
# details - a boolean. If TRUE the function returns a data frame with
#           the results of each individual simulation,
#           otherwise just the win proportion and probabilities of bein > control
#

runSims <- function(nSims=10, nBurnin=1000, rates=c(0.1, 0.1, 0.125, 0.15, 0.2, 0.25),
                    details=FALSE){
  winPpn = 0
  pr.gt.ctl.sum <- rep(0, length(rates)-1)
  
  if (details) {
    perSim <- data.frame(Sim=1)
  }
  
  for(sim in 1:nSims) {
  
    dat = genBinaryData(nPerArm=50, rates=rates)
    write.csv(dat,DATAFILE, row.names = FALSE)
    
    if (details) {
      perSim[sim, "Sim"] <- sim
      
      # record true rates and observed rates
      for ( d in 1:length(rates)) {
        perSim[sim, paste("sim.rate.", d, sep="")] <- rates[d]
      }
      for ( d in 1:length(rates)) {
        perSim[sim, paste("obs.rate.", d, sep="")] <- mean(dat[dat[,"Dose"]==d, "Visit1"])
      }
    }
    
    cat("run FACTS: ", sim, "\n")
    
    ret = runFACTS(engine='dichot', data.file = DATAFILE, param.file = 'bin1_e.param',
                 mcmc.file.num = 0, rng.seed = sim, exec.path = Exec.dir)
  
    dat = read.csv(MCMCFILE, skip=1)
    # discard burnin rows and just estimates of rate - the "Pi" columns
    dat = dat[(nBurnin+1):nrow(dat),grep("Pi", names(dat))]
    
    if (details) {
      # record est rate
      for ( d in 1:length(rates)) {
        perSim[sim, paste("est.rate.", d, sep="")] <- 
          mean(dat[,paste("Pi.", d, sep="")])
      }
    }
  
    # success if the first dose is not in the top 2 .. i.e. the resposnse on any 2 doses is > control
    success = apply(dat,1, 
                  FUN = function(x) {ifelse(length(x) - which(order(x)==1) >= 2, 1, 0)})
    
    if (details){
      perSim[sim, "Pr.Success"] <- mean(success)
      perSim[sim, "Success.flag"] <- ifelse(mean(success) > 0.9, 1,0)
    }
    
    winPpn = winPpn + ifelse(mean(success) > 0.9, 1,0)
  
    # example: calc pr(theta_d > theta_ctl)
    gt.ctl.flag <- apply(dat, 1, FUN=function(x){x[2:length(x)] > x[1]})
    pr.gt.ctl <- apply(gt.ctl.flag,1,sum)
    pr.gt.ctl <- pr.gt.ctl / length(gt.ctl.flag[1,])
    pr.gt.ctl.sum <- pr.gt.ctl.sum + pr.gt.ctl
    
    if (details) {
      for (d in 1:length(pr.gt.ctl)) {
        perSim[sim, paste("Pr.pi.", d+1, ">pi_ctl", sep="")] <- pr.gt.ctl[d]
      }
    }
  }
  
  cat("win proportion: ", winPpn/nSims, "\n")
  
  if (details)
    return (list(winPpn/nSims, pr.gt.ctl.sum/nSims, perSim))
  else
    return (list(winPpn/nSims, pr.gt.ctl.sum/nSims))
}



# example: calc pr(max)
# max.flag <- apply(dat, 1, FUN=function(x){x[2:length(x)] == max(x[2:length(x)])})
# pr.max <- apply(max.flag,1,sum)
# pr.max <- pr.max / length(max.flag[1,])

# example: calc pr(theta_d > theta_ctl)
# gt.ctl.flag <- apply(dat, 1, FUN=function(x){x[2:length(x)] > x[1]})
# pr.gt.ctl <- apply(gt.ctl.flag,1,sum)
# pr.gt.ctl <- pr.gt.ctl / length(gt.ctl.flag[1,])               