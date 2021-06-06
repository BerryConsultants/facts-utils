runFACTS = function(engine = c("contin","dichot","ME","TTE",
                               "ed_contin","ed_dichot","ed_tte"),
                    data.file = "patients.dat",
                    param.file = "nuk1_e.param",
                    mcmc.file.num = 0,
                    rng.seed,
                    exec.path = getwd()) {

  engine = match.arg(engine)
  if(.Platform$OS.type != "windows") {
    exec = switch(engine,
                  contin = "contin.x",
                  dichot = "dichot.x",
                  ME = "multep.x",
                  TTE = "tte.x",
                  ed_contin = "aipf_contin.x",
                  ed_dichot = "aipf_dichot.x",
                  ed_tte = "aipf_tte.x",
                  stop("Invalid FACTS engine name!"))
  }
  else {
    exec = switch(engine,
                  contin = "contin.exe",
                  dichot = "dichot.exe",
                  ME = "multep.exe",
                  TTE = "tte.exe",
                  ed_contin = "aipf_contin.exe",
                  ed_dichot = "aipf_dichot.exe",
                  ed_tte = "aipf_tte.exe",
                  stop("Invalid FACTS engine name!"))
  }

  cmd = paste(exec,
              "--analysis-mode",
              "f", param.file,
              "--analysis-data", data.file,
              "--mcmc-num", mcmc.file.num)
  if(!missing(rng.seed))
    cmd = paste(cmd,"--seed", rng.seed)
  cmd = paste(exec.path,cmd,sep="/")
  status = system(cmd)
  msg = (status==0)
}

