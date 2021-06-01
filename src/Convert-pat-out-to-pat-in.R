#
# Convert FACTS patient output files to FACTS patient input files for FACTS DF ME
#
# in.fn - name of input file
# out.fn - name of output file
# last.day - cut off data after this date
# visits - a vector of visit times (in days)
# v.ix.translate - needed, hopefully on temporarily, patients respose outputs are by 'endpoint' visit not
#   absolute visit. This matrix (nrows=number of repsones, ncols=overall number of visits), allows translation from
#   visit in the patient out file to the visit in the patient in file.
# 

FACTS.pat.convert.ME <- function(in.fn, out.fn, last.day, visits, v.ix.translate){
  browser()
  
  pat.in.df <- read.csv(in.fn, header=TRUE, skip=1)
  pat.out.created <- FALSE
  
  rnames <- c("Response.1", "Response.2", "Response.3", "Response.4")
  for (r.ix in 1:length(rnames)){
    for (v.ix in 1:length(visits)){
      resp.col.name <- paste(rnames[r.ix], v.ix, sep=".")
      column.index = match(resp.col.name, colnames(pat.in.df), nomatch=-1)
      if (column.index != -1) {
        # create visit records for response r.ix visit v.ix
        if (pat.out.created){
          pat.out.df <- rbind(pat.out.df, 
                              data.frame(Subject=pat.in.df$X..Subject, Dose=pat.in.df$Dose, 
                                         Endpoint=r.ix, Visit=v.ix.translate[r.ix, v.ix], 
                                         Response=pat.in.df[,column.index],
                                         Date=pat.in.df$Date + visits[v.ix]))
        } else {
          pat.out.df <- data.frame(Subject=pat.in.df$X..Subject, Dose=pat.in.df$Dose, 
                                   Endpoint=r.ix, Visit=v.ix.translate[r.ix, v.ix], 
                                   Response=pat.in.df[,column.index],
                                   Date=pat.in.df$Date + visits[v.ix])
          pat.out.created <- TRUE
        }
      }                         
    }  
  }
  browser()
  
  
  # remove entries beyond the suppled 'last.day' parameter and then drop Date field
  out.ix <- which(pat.out.df$Date < last.day)
  pat.out.df <- pat.out.df[out.ix, c("Subject", "Dose", "Endpoint", "Visit", "Response")]
  
  # sort records into Subject, Visit order
  sorted.ix <- order(pat.out.df$Subject, pat.out.df$Visit, pat.out.df$Endpoint)
  pat.out.df <- pat.out.df[sorted.ix,]
  write.csv(pat.out.df, out.fn, row.names=FALSE)
  
  return(pat.out.df)
}

#
# Convert FACTS patient output files to FACTS patient input files for DF TTE including predictors
#
# in.fn - name of input file
# out.fn - name of output file
# 
# Needs further work for TTE with TTE predictor

FACTS.pat.convert.TTE <- function(in.fn, out.fn, last.day){
  browser()
  
  pat.in.df <- read.csv(in.fn, header=TRUE, skip=1)
  
  pat.out.df <- data.frame(Subject=pat.in.df$X..Subject, Region=1, Dose=pat.in.df$Dose, 
                           Duration=pat.in.df$Duration, Outcome=pat.in.df$Outcome,
                           Predictor=pat.in.df$Predictor, Status=pat.in.df$Pred.Outcome,
                           Date=pat.in.df$Date, Date.Obs=pat.in.df$Date+pat.in.df$Duration, Dropout=pat.in.df$Dropout)
  
  # remove entries beyond the supplied 'last.day' parameter, censor events after 'last.day' and then drop Date field
  out.ix <- which(pat.out.df$Date < last.day)
  pat.out.df <- pat.out.df[out.ix,]
  
  censor.ix <- which(pat.out.df$Date.Obs > last.day)
  pat.out.df[censor.ix,]$Duration <- last.day - pat.out.df[censor.ix,]$Date
  pat.out.df[censor.ix,]$Outcome <- 0
  
  pat.out.df <- pat.out.df[out.ix, c("Subject", "Region", "Date", "Dose", "Duration", "Outcome", "Predictor", "Status", "Dropout")]

  write.csv(pat.out.df, out.fn, row.names=FALSE)
  
  return(pat.out.df)
}


#
# Convert FACTS patient output files to FACTS patient input files for AIPF Cts
#
# in.fn - name of input file
# out.fn - name of output file
# 
# 

FACTS.pat.convert.AIPF <- function(in.fn, out.fn, last.day, visits){
  browser()
  pat.in.df <- read.csv(in.fn, header=TRUE, skip=1)
  pat.out.created <- FALSE
  
  column.index = match("Baseline", colnames(pat.in.df), nomatch=-1)
  if (column.index != -1) {
    # create visit records for Baseline
    pat.out.df <- data.frame(Subject=pat.in.df$X.Subject, Group=pat.in.df$Group, 
                             Arm=pat.in.df$Arm, Visit=0, 
                             Response=pat.in.df[,column.index],
                             Last=pat.in.df$LastVisit.,
                             Date=pat.in.df$Date)
    pat.out.created <- TRUE
  }      
  
  
  for (v.ix in 1:length(visits)){
    resp.col.name <- paste("Visit", v.ix, sep=".")
    column.index = match(resp.col.name, colnames(pat.in.df), nomatch=-1)
    if (column.index != -1) {
      # create visit records for response r.ix visit v.ix
      
      if (pat.out.created){
        pat.out.df <- rbind(pat.out.df, 
                            data.frame(Subject=pat.in.df$X.Subject, Group=pat.in.df$Group, 
                                       Arm=pat.in.df$Arm, Visit=v.ix, 
                                       Response=pat.in.df[,column.index],
                                       Last=pat.in.df$LastVisit.,
                                       Date=pat.in.df$Date + visits[v.ix]))
      } else {
        pat.out.df <- data.frame(Subject=pat.in.df$X.Subject, Group=pat.in.df$Group, 
                                   Arm=pat.in.df$Arm, Visit=v.ix, 
                                   Response=pat.in.df[,column.index],
                                   Last=pat.in.df$LastVisit.,
                                   Date=pat.in.df$Date + visits[v.ix])
          pat.out.created <- TRUE
      }
    }                           
  }
  browser()
  

  keep.ix <- which(pat.out.df$Visit <= pat.out.df$Last+1)
  pat.out.df <- pat.out.df[keep.ix,]
  drop.ix <- which(pat.out.df$Visit == pat.out.df$Last+1)
  
  if (length(drop.ix > 0)) {  
    pat.out.df[drop.ix,]$Visit <- -1
  }
  
  # remove entries beyond the suppled 'last.day' parameter and then drop Date & last field
  out.ix <- which(pat.out.df$Date < last.day)
  pat.out.df <- pat.out.df[out.ix, c("Subject", "Group", "Arm", "Visit", "Response")]
  
  # sort records into Subject, Visit order
  sorted.ix <- order(pat.out.df$Subject, pat.out.df$Visit)
  pat.out.df <- pat.out.df[sorted.ix,]
  write.csv(pat.out.df, out.fn, row.names=FALSE)
  
  return(pat.out.df)
}

#
# Convert FACTS patient output files to FACTS patient input files for AIPF TTE including predictors
#
# in.fn - name of input file
# out.fn - name of output file
# 
# Needs further work for TTE with TTE predictor

FACTS.pat.convert.AIPF.TTE <- function(in.fn, out.fn, last.day){
  browser()
  
  pat.in.df <- read.csv(in.fn, header=TRUE, skip=1)
  
  pat.out.df <- data.frame(Subject=pat.in.df$X.Subject, Group=pat.in.df$Group, Arm=pat.in.df$Arm, 
                           Duration=pat.in.df$Duration, Outcome=pat.in.df$Outcome,
                           Date=pat.in.df$Date, Date.Obs=pat.in.df$Date+pat.in.df$Duration)
  
  # remove entries beyond the supplied 'last.day' parameter, censor events after 'last.day' and then drop Date field
  out.ix <- which(pat.out.df$Date < last.day)
  pat.out.df <- pat.out.df[out.ix,]
  
  censor.ix <- which(pat.out.df$Date.Obs > last.day)
  pat.out.df[censor.ix,]$Duration <- last.day - pat.out.df[censor.ix,]$Date
  pat.out.df[censor.ix,]$Outcome <- 0
  
  pat.out.df <- pat.out.df[out.ix, c("Subject", "Date", "Group", "Arm", "Duration", "Outcome")]
  
  write.csv(pat.out.df, out.fn, row.names=FALSE)
  
  return(pat.out.df)
}

tte.f.in <- "C:\\FACTS\\FACTS 3.6.2 testing\\test 4 arm cts predictor_results\\Accrual 1_Dropout 1_Response 1_Control Hazard 1_Predictor 1\\patients00001.csv"
tte.f.out <- "C:\\FACTS\\FACTS 3.6.2 testing\\test 4 arm cts predictor_results\\patients00001_converted.csv"
#FACTS.pat.convert.TTE(tte.f.in, tte.f.out, 300)

f.in <- "C:\\FACTS\\FACTS 3.6.2 testing\\Pain - 2 endpoint + LM example_results\\Accrual 1_Dropout 1_Linear\\patients00001.csv"
f.out <- "C:\\FACTS\\FACTS 3.6.2 testing\\Pain - 2 endpoint + LM example_results\\patients00001_converted.csv"
#FACTS.pat.convert.ME(f.in, f.out, 129, c(7, 14, 21, 28, 35, 42), c(c(1,2,3,4,5,6), c(2,0,0,0,0,0,0), c(0,0,0,0,0,0), C(0,0,0,0,0,0)))

#aipf.f.in <- "C:\\FACTS\\FACTS 3.6.2 testing\\AIPF Cts Alzheimers + baseline + nes_results\\Accrual 1_Dropout 1_Baseline 1_Nugget 2_Linear\\patients00001.csv"
aipf.f.in <- "C:\\FACTS\\FACTS 3.6.2 testing\\patients00001.csv"
aipf.f.out <- "C:\\FACTS\\FACTS 3.6.2 testing\\AIPF Cts Alzheimers + baseline + nes_results\\patients00001_converted.csv"

aipf.d.f.in <- "C:\\FACTS\\FACTS 3.6.2 testing\\AIPF Dich Cancer Trial RM - es_results\\Accrual 1_Dropout 1_MixedModerate\\patients00001.csv"
aipf.d.f.out <- "C:\\FACTS\\FACTS 3.6.2 testing\\AIPF Dich Cancer Trial RM - es_results\\patients00001_converted.csv"
#aipf.d.f.in <- "C:\\FACTS\\FACTS 3.6.2 testing\\AIPF Dich Cancer Trial LM+Borrow_results\\Accrual 1_Dropout 1_MixedModerate_Longitudinal 1\\patients00001.csv"
#aipf.d.f.out <- "C:\\FACTS\\FACTS 3.6.2 testing\\AIPF Dich Cancer Trial LM+Borrow_results\\patients00001_converted.csv"
#FACTS.pat.convert.AIPF(aipf.f.in, aipf.f.out, 400, c(13, 24, 52))
#FACTS.pat.convert.AIPF(aipf.d.f.in, aipf.d.f.out, 400, c(14, 28, 42, 56, 70, 84))
aipf.tte.f.in <- "C:\\FACTS\\FACTS 3.6.2 testing\\AIPF TTE Genomic Biomarker Ex1_results\\Accrual 1_Dropout 1_Median8_StrongWWX\\patients00001.csv"
aipf.tte.f.out <- "C:\\FACTS\\FACTS 3.6.2 testing\\AIPF TTE Genomic Biomarker Ex1_results\\patients00001_converted.csv"
#FACTS.pat.convert.AIPF.TTE(aipf.tte.f.in, aipf.tte.f.out, 600)