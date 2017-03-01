
#' @title Births Module for up to 5 race groups heterosexuals and MSM.
#'
#' @description Module function for births or entries into the sexually active
#'              population.
#'
#' @inheritParams aging_msm
#'
#' @details
#' New population members are added based on expected numbers of entries among
#' all five race/immigrant groups and two sexes, stochastically determined with draws from Poisson
#' distributions. THe proportion of men who are MSM are determined by \code{msm.frac}. 
#' THe proportion of men who are MSMF are determined by \code{msmf.frac}. 
#' For each new entry, a set of attributes is added for that node,
#' and the nodes are added onto the network objects. Only attributes that are
#' a part of the network model formula are updated as vertex attributes on the
#' network objects.
#'
#' @return
#' This function updates the \code{attr} list with new attributes for each new
#' population member, and the \code{nw} objects with new vertices.
#'
#' @keywords module msm het
#' @export
#'
births_shamp <- function(dat, at){

  ## Variables

  # Parameters
  b.B.f.rate <- dat$param$b.B.f.rate
  b.BI.f.rate <- dat$param$b.BI.f.rate
  b.H.f.rate <- dat$param$b.H.f.rate
  b.HI.f.rate <- dat$param$b.HI.f.rate 
  b.W.f.rate <- dat$param$b.W.f.rate
  
  b.B.m.rate <- dat$param$b.B.m.rate
  b.BI.m.rate <- dat$param$b.BI.m.rate
  b.H.m.rate <- dat$param$b.H.m.rate
  b.HI.m.rate <- dat$param$b.HI.m.rate 
  b.W.m.rate <- dat$param$b.W.m.rate

  b.method <- dat$param$b.method


  ## Process
  if (b.method == "fixed") {
    numB.f <- dat$epi$num.B.f[1]
    numBI.f <- dat$epi$num.BI.f[1]
    numH.f <- dat$epi$num.H.f[1]
    numHI.f <- dat$epi$num.HI.f[1]
    numW.f <- dat$epi$num.W.f[1]
    
    numB.m <- dat$epi$num.B.m[1]
    numBI.m <- dat$epi$num.BI.m[1]
    numH.m <- dat$epi$num.H.m[1]
    numHI.m <- dat$epi$num.HI.m[1]
    numW.m <- dat$epi$num.W.m[1]
    

  }
  if (b.method == "varying") {
    numB.f <- dat$epi$num.B.f[at - 1]
    numBI.f <- dat$epi$num.BI.f[at - 1]
    numH.f <- dat$epi$num.H.f[at - 1]
    numHI.f <- dat$epi$num.HI.f[at - 1]
    numW.f <- dat$epi$num.W.f[at - 1]
    
    numB.m <- dat$epi$num.B.m[at - 1]
    numBI.m <- dat$epi$num.BI.m[at - 1]
    numH.m <- dat$epi$num.H.m[at - 1]
    numHI.m <- dat$epi$num.HI.m[at - 1]
    numW.m <- dat$epi$num.W.m[at - 1]
  }

  nBirths.B.f <- rpois(1, b.B.f.rate * numB.f)
  nBirths.BI.f <- rpois(1, b.BI.f.rate * numBI.f)
  nBirths.H.f <- rpois(1, b.H.f.rate * numH.f)
  nBirths.HI.f <- rpois(1, b.HI.f.rate * numHI.f)
  nBirths.W.f <- rpois(1, b.W.f.rate * numW.f)

  nBirths.B.m <- rpois(1, b.B.m.rate * numB.m)
  nBirths.BI.m <- rpois(1, b.BI.m.rate * numBI.m)
  nBirths.H.m <- rpois(1, b.H.m.rate * numH.m)
  nBirths.HI.m <- rpois(1, b.HI.m.rate * numHI.m)
  nBirths.W.m <- rpois(1, b.W.m.rate * numW.m)
  
  nBirths <- nBirths.B.f + nBirths.BI.f +  nBirths.H.f +  nBirths.HI.f +
  nBirths.W.f +  nBirths.B.m +  nBirths.BI.m +  nBirths.H.m +
  nBirths.HI.m +  nBirths.W.m


  ## Update Attr

  if (nBirths > 0) {
    dat <- setBirthAttr_shamp(dat, at, nBirths.B.f, nBirths.BI.f, nBirths.H.f, nBirths.HI.f,
                              nBirths.W.f, nBirths.B.m, nBirths.BI.m, nBirths.H.m,
                              nBirths.HI.m, nBirths.W.m)
  }


  # Update Networks
  if (nBirths > 0) {
    for (i in 1:3) {
      dat$el[[i]] <- add_vertices(dat$el[[i]], nBirths)
    }
  }


  ## Output
  dat$epi$nBirths[at] <- nBirths

  return(dat)
}


setBirthAttr_shamp <- function(dat, at, nBirths.B.f, nBirths.BI.f, nBirths.H.f, nBirths.HI.f,
                               nBirths.W.f, nBirths.B.m, nBirths.BI.m, nBirths.H.m,
                               nBirths.HI.m, nBirths.W.m) {

  nBirths <- nBirths.B.f + nBirths.BI.f +  nBirths.H.f +  nBirths.HI.f +
    nBirths.W.f +  nBirths.B.m +  nBirths.BI.m +  nBirths.H.m +
    nBirths.HI.m +  nBirths.W.m

  # Set all attributes NA by default
  dat$attr <- lapply(dat$attr, {
    function(x)
      c(x, rep(NA, nBirths))
  })
  newIds <- which(is.na(dat$attr$active))

  # Demographic
  dat$attr$active[newIds] <- rep(1, nBirths)
  dat$attr$uid[newIds] <- dat$temp$max.uid + (1:nBirths)
  dat$temp$max.uid <- dat$temp$max.uid + nBirths
  dat$attr$sex.ident[newIds] <- rep("NA", nBirths)
  dat$attr$immig.loc[newIds] <- rep(0, nBirths)

  dat$attr$arrival.time[newIds] <- rep(at, nBirths)

  sex.race <- sample(rep(c("B.f", "BI.f", "H.f", "HI.f", "W.f","B.m", "BI.m", "H.m", "HI.m", "W.m"), 
                    c(nBirths.B.f, nBirths.BI.f, nBirths.H.f, nBirths.HI.f, nBirths.W.f,
                      nBirths.B.m, nBirths.BI.m, nBirths.H.m, nBirths.HI.m, nBirths.W.m)))
  newF <- which(sex.race == "B.f" | sex.race == "BI.f" | sex.race == "H.f" | sex.race == "HI.f" | sex.race == "W.f")
  newM <- which(sex.race == "B.m" | sex.race == "BI.m" | sex.race == "H.m" | sex.race == "HI.m" | sex.race == "W.m")
  sex<-rep(NA,length(newIds))
  sex[newF]<-"F"
  sex[newM]<-"M"
  dat$attr$sex[newIds] <- sex
  
  sex.ident<-c(newF,newM)
  sex.ident[newF]<-"f"
  sex.ident[newM]<-sample(c("msm","msmf","msf"),length(newM),
                                    prob=c(dat$param$msm.frac,dat$param$msmf.frac,(1-(dat$param$msm.frac+dat$param$msmf.frac))),
                                           replace=TRUE)
  dat$attr$sex.ident[newIds] <- sex.ident
 
  newB <- which(sex.race == "B.f" | sex.race == "B.m")
  newBI <- which(sex.race == "BI.f" | sex.race == "BI.m")
  newH <- which(sex.race == "H.f" | sex.race == "H.m")
  newHI <- which(sex.race == "HI.f" | sex.race == "HI.m")
  newW <- which(sex.race == "W.f" | sex.race == "W.m")
  new.msf <- which(sex.ident == "msf")
  new.msm <- which(sex.ident == "msm")
  new.msmf <- which(sex.ident == "msmf")
  race<-rep(NA,length(newIds))
  race[newB]<-"B"  
  race[newBI]<-"BI"
  race[newH]<-"H"
  race[newHI]<-"HI"
  race[newW]<-"W"
  
  dat$attr$race[newIds] <- race

  dat$attr$age[newIds] <- rep(dat$param$birth.age, nBirths)
  dat$attr$sqrt.age[newIds] <- sqrt(dat$attr$age[newIds])

 

  dat$attr$status[newIds] <- rep(0, nBirths)
  dat$attr$aq.class[newIds] <- rep(NA, nBirths)

 
  # dat$attr$inst.ai.class[newIds] <- sample(1:dat$param$num.inst.ai.classes,
  #                                         nBirths, replace = TRUE)
 
   #### females f. 
  
  if (length(newF) > 0) {
  dat$attr$tt.traj[newIds[intersect(newB,newF)]] <- sample(c(1, 2, 3, 4),
                                           nBirths.B.f, replace = TRUE,
                                           prob = dat$param$tt.traj.B.f.prob)
  dat$attr$tt.traj[newIds[intersect(newBI,newF)]] <- sample(c(1, 2, 3, 4),
                                           nBirths.BI.f, replace = TRUE,
                                           prob = dat$param$tt.traj.BI.f.prob)
  dat$attr$tt.traj[newIds[intersect(newH,newF)]] <- sample(c(1, 2, 3, 4),
                                           nBirths.H.f, replace = TRUE,
                                           prob = dat$param$tt.traj.H.f.prob)
  dat$attr$tt.traj[newIds[intersect(newHI,newF)]] <- sample(c(1, 2, 3, 4),
                                           nBirths.HI.f, replace = TRUE,
                                           prob = dat$param$tt.traj.HI.f.prob)
  dat$attr$tt.traj[newIds[intersect(newW,newF)]] <- sample(c(1, 2, 3, 4),
                                           nBirths.W.f, replace = TRUE,
                                           prob = dat$param$tt.traj.W.f.prob)
  }
  
  #### Males msf. 
  if (length(new.msf) > 0) {
  dat$attr$tt.traj[newIds[intersect(newB,new.msf)]] <- sample(c(1, 2, 3, 4),
                                                  length(intersect(newB,new.msf)), replace = TRUE,
                                                 prob = dat$param$tt.traj.B.msf.prob)
  dat$attr$tt.traj[newIds[intersect(newBI,new.msf)]] <- sample(c(1, 2, 3, 4),
                                                  length(intersect(newBI,new.msf)), replace = TRUE,
                                                  prob = dat$param$tt.traj.BI.msf.prob)
  dat$attr$tt.traj[newIds[intersect(newH,new.msf)]] <- sample(c(1, 2, 3, 4),
                                                  length(intersect(newH,new.msf)), replace = TRUE,
                                                 prob = dat$param$tt.traj.H.msf.prob)
  dat$attr$tt.traj[newIds[intersect(newHI,new.msf)]] <- sample(c(1, 2, 3, 4),
                                                  length(intersect(newHI,new.msf)), replace = TRUE,
                                                  prob = dat$param$tt.traj.HI.msf.prob)
  dat$attr$tt.traj[newIds[intersect(newW,new.msf)]] <- sample(c(1, 2, 3, 4),
                                                  length(intersect(newW,new.msf)), replace = TRUE,
                                                 prob = dat$param$tt.traj.W.msf.prob)
  }
  
  #### Males msm. 
  if(length(new.msm) > 0){
  dat$attr$tt.traj[newIds[intersect(newB,new.msm)]] <- sample(c(1, 2, 3, 4),
                                                      length(intersect(newB,new.msm)), replace = TRUE,
                                                      prob = dat$param$tt.traj.B.msm.prob)
  dat$attr$tt.traj[newIds[intersect(newBI,new.msm)]] <- sample(c(1, 2, 3, 4),
                                                      length(intersect(newBI,new.msm)), replace = TRUE,
                                                       prob = dat$param$tt.traj.BI.msm.prob)
  dat$attr$tt.traj[newIds[intersect(newH,new.msm)]] <- sample(c(1, 2, 3, 4),
                                                      length(intersect(newH,new.msm)), replace = TRUE,
                                                      prob = dat$param$tt.traj.H.msm.prob)
  dat$attr$tt.traj[newIds[intersect(newHI,new.msm)]] <- sample(c(1, 2, 3, 4),
                                                      length(intersect(newHI,new.msm)), replace = TRUE,
                                                       prob = dat$param$tt.traj.HI.msm.prob)
  dat$attr$tt.traj[newIds[intersect(newW,new.msm)]] <- sample(c(1, 2, 3, 4),
                                                      length(intersect(newW,new.msm)), replace = TRUE,
                                                      prob = dat$param$tt.traj.W.msm.prob)
  }
  
  #### Males msmf. 
  if(length(new.msmf) > 0){
    dat$attr$tt.traj[newIds[intersect(newB,new.msmf)]] <- sample(c(1, 2, 3, 4),
                                                            length(intersect(newB,new.msmf)), replace = TRUE,
                                                            prob = dat$param$tt.traj.B.msmf.prob)
    dat$attr$tt.traj[newIds[intersect(newBI,new.msmf)]] <- sample(c(1, 2, 3, 4),
                                                             length(intersect(newBI,new.msmf)), replace = TRUE,
                                                             prob = dat$param$tt.traj.BI.msmf.prob)
    dat$attr$tt.traj[newIds[intersect(newH,new.msmf)]] <- sample(c(1, 2, 3, 4),
                                                            length(intersect(newH,new.msmf)), replace = TRUE,
                                                            prob = dat$param$tt.traj.H.msmf.prob)
    dat$attr$tt.traj[newIds[intersect(newHI,new.msmf)]] <- sample(c(1, 2, 3, 4),
                                                             length(intersect(newHI,new.msmf)), replace = TRUE,
                                                             prob = dat$param$tt.traj.HI.msmf.prob)
    dat$attr$tt.traj[newIds[intersect(newW,new.msmf)]] <- sample(c(1, 2, 3, 4),
                                                            length(intersect(newW,new.msmf)), replace = TRUE,
                                                            prob = dat$param$tt.traj.W.msmf.prob)
  }
  
  
  # Circumcision
  dat$attr$circ[newIds[intersect(newB,newM)]] <- rbinom(nBirths.B.m, 1, dat$param$circ.B.prob)
  dat$attr$circ[newIds[intersect(newBI,newM)]] <- rbinom(nBirths.BI.m, 1, dat$param$circ.BI.prob)
  dat$attr$circ[newIds[intersect(newH,newM)]] <- rbinom(nBirths.H.m, 1, dat$param$circ.H.prob)
  dat$attr$circ[newIds[intersect(newHI,newM)]] <- rbinom(nBirths.HI.m, 1, dat$param$circ.HI.prob)
  dat$attr$circ[newIds[intersect(newW,newM)]] <- rbinom(nBirths.W.m, 1, dat$param$circ.W.prob)
  
  dat$attr$circ[newIds[newF]] <- 0
  
  # Role
  dat$attr$role.class[newIds[newF]] <- NA
  dat$attr$role.class[newIds[newM]] <- NA
  
  if (length(new.msm) > 0) {

  
  dat$attr$role.class[newIds[intersect(newB,new.msm)]] <- sample(c("I", "R", "V"),
                                              length(intersect(newB,new.msm)), replace = TRUE,
                                              prob = dat$param$role.B.msm.prob)
  dat$attr$role.class[newIds[intersect(newBI,new.msm)]] <- sample(c("I", "R", "V"),
                                              length(intersect(newBI,new.msm)), replace = TRUE,
                                              prob = dat$param$role.BI.msm.prob)
  dat$attr$role.class[newIds[intersect(newH,new.msm)]] <- sample(c("I", "R", "V"),
                                              length(intersect(newH,new.msm)), replace = TRUE,
                                              prob = dat$param$role.H.msm.prob)
  dat$attr$role.class[newIds[intersect(newHI,new.msm)]] <- sample(c("I", "R", "V"),
                                              length(intersect(newHI,new.msm)), replace = TRUE,
                                               prob = dat$param$role.HI.msm.prob)
  dat$attr$role.class[newIds[intersect(newW,new.msm)]] <- sample(c("I", "R", "V"),
                                              length(intersect(newW,new.msm)), replace = TRUE,
                                              prob = dat$param$role.W.msm.prob)
  
  }
  
  if (length(new.msmf) > 0) {
  dat$attr$role.class[newIds[intersect(newB,new.msmf)]] <- sample(c("I", "R", "V"),
                                                             length(intersect(newB,new.msmf)), replace = TRUE,
                                                             prob = dat$param$role.B.msmf.prob)
  dat$attr$role.class[newIds[intersect(newBI,new.msmf)]] <- sample(c("I", "R", "V"),
                                                              length(intersect(newBI,new.msmf)), replace = TRUE,
                                                              prob = dat$param$role.BI.msmf.prob)
  dat$attr$role.class[newIds[intersect(newH,new.msmf)]] <- sample(c("I", "R", "V"),
                                                             length(intersect(newH,new.msmf)), replace = TRUE,
                                                             prob = dat$param$role.H.msmf.prob)
  dat$attr$role.class[newIds[intersect(newHI,new.msmf)]] <- sample(c("I", "R", "V"),
                                                              length(intersect(newHI,new.msmf)), replace = TRUE,
                                                              prob = dat$param$role.HI.msmf.prob)
  dat$attr$role.class[newIds[intersect(newW,new.msmf)]] <- sample(c("I", "R", "V"),
                                                             length(intersect(newW,new.msmf)), replace = TRUE,
                                                             prob = dat$param$role.W.msmf.prob)
  }

  ins.quot <- rep(NA,nBirths)
 
  ins.quot[dat$attr$role.class[newIds] == "I"]  <- 1
  ins.quot[dat$attr$role.class[newIds] == "R"]  <- 0
 
  ins.quot.v <- runif(sum(dat$attr$role.class[newIds] == "V",na.rm = TRUE))
  x <- which(dat$attr$role.class[newIds] == "V")
  ins.quot[x]  <-ins.quot.v
                                 
  dat$attr$ins.quot[newIds] <- ins.quot


  # CCR5
  ccr5.B.f.prob <- dat$param$ccr5.B.f.prob
  ccr5.BI.f.prob <- dat$param$ccr5.BI.f.prob
  ccr5.H.f.prob <- dat$param$ccr5.H.f.prob
  ccr5.HI.f.prob <- dat$param$ccr5.HI.f.prob
  ccr5.W.f.prob <- dat$param$ccr5.W.f.prob
  
  
  
  dat$attr$ccr5[newIds[intersect(newB,newF)]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.B.f, replace = TRUE,
                                        prob = c(1 - sum(ccr5.B.f.prob),
                                                 ccr5.B.f.prob[2], ccr5.B.f.prob[1]))
  dat$attr$ccr5[newIds[intersect(newBI,newF)]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.BI.f, replace = TRUE,
                                        prob = c(1 - sum(ccr5.BI.f.prob),
                                                 ccr5.BI.f.prob[2], ccr5.BI.f.prob[1]))
  dat$attr$ccr5[newIds[intersect(newH,newF)]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.H.f, replace = TRUE,
                                        prob = c(1 - sum(ccr5.H.f.prob),
                                                 ccr5.H.f.prob[2], ccr5.H.f.prob[1]))
  dat$attr$ccr5[newIds[intersect(newHI,newF)]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.HI.f, replace = TRUE,
                                        prob = c(1 - sum(ccr5.HI.f.prob),
                                                 ccr5.HI.f.prob[2], ccr5.HI.f.prob[1]))
  dat$attr$ccr5[newIds[intersect(newW,newF)]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.W.f, replace = TRUE,
                                        prob = c(1 - sum(ccr5.W.f.prob),
                                                 ccr5.W.f.prob[2], ccr5.W.f.prob[1]))

  ccr5.B.m.prob <- dat$param$ccr5.B.m.prob
  ccr5.BI.m.prob <- dat$param$ccr5.BI.m.prob
  ccr5.H.m.prob <- dat$param$ccr5.H.m.prob
  ccr5.HI.m.prob <- dat$param$ccr5.HI.m.prob
  ccr5.W.m.prob <- dat$param$ccr5.W.m.prob
  
  
  
  dat$attr$ccr5[newIds[intersect(newB,newM)]] <- sample(c("WW", "DW", "DD"),
                                              nBirths.B.m, replace = TRUE,
                                              prob = c(1 - sum(ccr5.B.m.prob),
                                                       ccr5.B.m.prob[2], ccr5.B.m.prob[1]))
  dat$attr$ccr5[newIds[intersect(newBI,newM)]] <- sample(c("WW", "DW", "DD"),
                                               nBirths.BI.m, replace = TRUE,
                                               prob = c(1 - sum(ccr5.BI.m.prob),
                                                        ccr5.BI.m.prob[2], ccr5.BI.m.prob[1]))
  dat$attr$ccr5[newIds[intersect(newH,newM)]] <- sample(c("WW", "DW", "DD"),
                                              nBirths.H.m, replace = TRUE,
                                              prob = c(1 - sum(ccr5.H.m.prob),
                                                       ccr5.H.m.prob[2], ccr5.H.m.prob[1]))
  dat$attr$ccr5[newIds[intersect(newHI,newM)]] <- sample(c("WW", "DW", "DD"),
                                               nBirths.HI.m, replace = TRUE,
                                               prob = c(1 - sum(ccr5.HI.m.prob),
                                                        ccr5.HI.m.prob[2], ccr5.HI.m.prob[1]))
  dat$attr$ccr5[newIds[intersect(newW,newM)]] <- sample(c("WW", "DW", "DD"),
                                              nBirths.W.m, replace = TRUE,
                                              prob = c(1 - sum(ccr5.W.m.prob),
                                                       ccr5.W.m.prob[2], ccr5.W.m.prob[1]))

  # Degree
  dat$attr$deg.main[newIds] <- 0
  dat$attr$deg.pers[newIds] <- 0

  # One-off risk group
  dat$attr$riskg[newIds] <- sample(1:5, nBirths, TRUE)

  # UAI group
  p1.msm <- dat$param$cond.pers.always.prob.msm
  p2.msm <- dat$param$cond.inst.always.prob.msm
  rho.msm <- dat$param$cond.always.prob.corr.msm
  uai.always.msm <- bindata::rmvbin(nBirths, c(p1.msm, p2.msm), bincorr = (1 - rho.msm) * diag(2) + rho.msm)
  dat$attr$cond.always.pers.msm[newIds] <- uai.always.msm[, 1]
  dat$attr$cond.always.inst.msm[newIds] <- uai.always.msm[, 2]
  
  # UVI group
  p1.het <- dat$param$cond.pers.always.prob.het
  p2.het <- dat$param$cond.inst.always.prob.het
  rho.het <- dat$param$cond.always.prob.corr.het
  uvi.always.het <- bindata::rmvbin(nBirths, c(p1.het, p2.het), bincorr = (1 - rho.het) * diag(2) + rho.het)
  dat$attr$cond.always.pers.het[newIds] <- uvi.always.het[, 1]
  dat$attr$cond.always.inst.het[newIds] <- uvi.always.het[, 2]

  # PrEP
  dat$attr$prepStat[newIds] <- 0

  return(dat)
}

