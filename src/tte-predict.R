#
# function to generate TTE predictor and TTE endpoint using the 
# FACTS predictor | endpoint model
#
# Z is subjects predictor time, Y is subjects endpoint time:
# Z ~ Exp(lambda_z * e ^ (Beta * Y))
#
# hr - control hazard rate
# HR - hazard ratio for comparison arm
# lz - lambda_z parameter for model
# beta - beta parameter for model
# fup - followup time
#
pred.tte.hr <- function(hr, HR, beta, lz, fup) {
  ctl <- rexp(100000, rate=hr)
  trt <- rexp(100000, rate=HR*hr)
  ctl_p <- sapply(ctl, FUN=function(x) rexp(1, rate=lz*exp(beta*x)))
  # ctl_pc predictor times censored by final event
  ctl_pc <- apply(cbind(ctl_p, ctl), MARGIN=1, FUN=min)
  #ctl_pcc predictor times censored by final event and max follow up time
  ctl_pcc <- apply(cbind(ctl_p, ctl, fup), MARGIN=1, FUN=min)
  
  trt_p <- sapply(trt, FUN=function(x) rexp(1, rate=lz*exp(beta*x)))
  trt_pc <- apply(cbind(trt_p, trt), MARGIN=1, FUN=min)
  trt_pcc <- apply(cbind(trt_p, trt, fup), MARGIN=1, FUN=min)
  
  # return:
  # median predictor time on control, median censored predictor time on control
  # HR = ratio of uncensored means
  # HR = ratio of observed events / total exposure
  return(list(median(ctl_p), median(ctl_pcc), mean(ctl_pc)/mean(trt_pc), 
              (sum(trt_pc<fup)/sum(trt_pcc))/(sum(ctl_pc<fup)/sum(ctl_pcc))))
}

#
# function to generate TTE predictor and TTE endpoint using the 
# FACTS endpoint model | predictor
#
# p_hr - predictor control hazard rate
# p_HR - predictor hazard ration for comparison arm
# f_hr - post predictor hazard rate for both arms
# p_z - probability that post predictor time is zero
# fup - follow up time 
#
tte.pred.hr <- function(p_hr, p_HR, f_hr, p_z, fup) {
  pfs1 <- rexp(100000, rate=p_hr)
  pfs2 <- rexp(100000, rate=p_HR*p_hr)
  os1 <- ifelse(rbinom(100000, size=1, prob=p_z) > 0, 0, rexp(100000, rate=f_hr))
  os2 <- ifelse(rbinom(100000, size=1, prob=p_z) > 0, 0, rexp(100000, rate=f_hr))
  ov1=pfs1+os1
  ov2=pfs2+os2
  
  #ov1_c, ov2_c: times censored by max follow up time
  ov1_c <- apply(cbind(ov1, fup), MARGIN=1, FUN=min)
  ov2_c <- apply(cbind(ov2, fup), MARGIN=1, FUN=min)
  
  # return:
  # median time on control, median censored time on control 
  # HR = ratio of uncensored means
  # HR = ratio of observed events / total exposure
  return(list(median(ov1), median(ov1_c), mean(ov1) / mean(ov2),
              (sum(ov2_c<fup)/sum(ov2_c))/(sum(ov1_c<fup)/sum(ov1_c))))
}

#
# function to return observed HR of two exp populations
# taking into account the censoring of finite followup
#
tte.hr <- function(hr, HR, fup) {
  ctl <- rexp(100000, rate=hr)
  trt <- rexp(100000, rate=HR*hr)

  #ctl_c, trt_c: times censored by max follow up time
  ctl_c <- apply(cbind(ctl, fup), MARGIN=1, FUN=min)
  trt_c <- apply(cbind(trt, fup), MARGIN=1, FUN=min)
  
  # return:
  # median time on control
  # HR = ratio of uncensored means
  # HR = ratio of observed events / total exposure
  return(list(median(ctl), median(ctl_c), mean(ctl)/mean(trt), 
              (sum(trt<fup)/sum(trt_c))/(sum(ctl<fup)/sum(ctl_c))))
}
