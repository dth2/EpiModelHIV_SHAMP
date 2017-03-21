

#' @title Transmission Module for Heat Bath Exposure of MSMF to an MSM heatbath
#'
#' @description Stochastically simulates disease transmission from the heat bath given the current
#'              state of the active HIV negative population.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This function takes all active HIV negative MSMF nodes and calculates a
#' transmission probability for each individual based on age specific exposure to the heatbath
#' After transmission events, individual-level attributes for the infected
#' persons are updated and summary statistics for incidence calculated.
#'
#' The per-act transmission probability depends on the following elements:
#' an age and risk catagory specific probability of aquiring HIV from the heat bath per time step
#' The input values \code{dat$param$msm.aq.prob.'R'} were R is a race specific probability taking values (B,BI,H,HI,W).
#' 
#'
#' @return
#' For each new infection, the disease status, infection time, and related
#' HIV attributes are updated for the infected node. Summary statistics for
#' disease incidence overall, and by race and age groups are calculated and
#' stored on \code{dat$epi}.
#'
#' @keywords module
#' @export
#'

heatbath.adol <- function(dat, at){
  
  
  ## Infected from the heat bath
  #Attributes
  age<-floor(dat$attr$age)
  #riskg<-dat$attr$riskg
  race<-dat$attr$race
  active<-dat$attr$active
  status<-dat$attr$status
  prepStat<-dat$attr$prepStat
  prepClass<-dat$attr$prepClass
  sex.ident<-dat$attr$sex.ident
  

  #Parameters
  prep.hr <- dat$param$prep.class.hr
  prep.class.prob <- dat$param$prep.class.prob
  
  msm.aq.prob.B <- dat$param$msm.aq.prob.B
  msm.aq.prob.BI <- dat$param$msm.aq.prob.BI
  msm.aq.prob.H <- dat$param$msm.aq.prob.H
  msm.aq.prob.HI <- dat$param$msm.aq.prob.HI
  msm.aq.prob.W <- dat$param$msm.aq.prob.W


  #Select active msmf by race
  ids.B<-which(active==1 & sex.ident=="msmf" & race =="B")
  ids.BI<-which(active==1 & sex.ident=="msmf" & race =="BI")
  ids.H<-which(active==1 & sex.ident=="msmf" & race =="H")
  ids.HI<-which(active==1 & sex.ident=="msmf" & race =="HI")
  ids.W<-which(active==1 & sex.ident=="msmf" & race =="W")
  
  #Determine who has Contact and who becomes infected 
  ##HERE HERE HERE
  age.b<-age[ids.b]
  status.b<-status[ids.b]
  riskg.b<-riskg[ids.b]
  prepStat.b<-prepStat[ids.b]
  prepClass.b<-prepClass[ids.b]

  
  age.w<-age[ids.w]
  status.w<-status[ids.w]
  riskg.w<-riskg[ids.w]
  prepStat.w<-prepStat[ids.w]
  prepClass.w<-prepClass[ids.w]

  
  
  infect.new.b<-rep(0,length(ids.b))
  infect.new.w<-rep(0,length(ids.w))
  
  infect.new.mat.in.b<-cbind(ids.b,status.b,riskg.b,age.b,prepStat.b,prepClass.b,AI.b,infect.new.b)
  infect.new.mat.out.b<-NULL
  age.vec.b<-sort(unique(age.b))
  riskg.vec.b<-sort(unique(riskg.b))
  
  infect.new.mat.in.w<-cbind(ids.w,status.w,riskg.w,age.w,prepStat.w,prepClass.w,AI.w,infect.new.w)
  infect.new.mat.out.w<-NULL
  age.vec.w<-sort(unique(age.w))
  riskg.vec.w<-sort(unique(riskg.w))
  
#Calculate the number of Contacts.
  for(i in 1:length(riskg.vec.b)){
      for(j in 1:length(age.vec.b)){
        
      infect.new.mat.temp.b<-subset(infect.new.mat.in.b,infect.new.mat.in.b[,"riskg.b"]==riskg.vec.b[i] & infect.new.mat.in.b[,"age.b"]==age.vec.b[j])
      infect.new.mat.temp.b[,"AI.b"]<-rpois(length(infect.new.mat.temp.b[,"AI.b"]),heat.b[i,j])
      infect.new.mat.out.b<-rbind(infect.new.mat.out.b,infect.new.mat.temp.b)

      }
    }
  
  for(i in 1:length(riskg.vec.w)){
    for(j in 1:length(age.vec.w)){
      
      infect.new.mat.temp.w<-subset(infect.new.mat.in.w,infect.new.mat.in.w[,"riskg.w"]==riskg.vec.w[i] & infect.new.mat.in.w[,"age.w"]==age.vec.w[j])
      infect.new.mat.temp.w[,"AI.w"]<-rpois(length(infect.new.mat.temp.w[,"AI.w"]),heat.w[i,j])
      infect.new.mat.out.w<-rbind(infect.new.mat.out.w,infect.new.mat.temp.w)
      
    }
  }
  
 
    #Select those that had AI.
    AI.mat.out.b<-subset(infect.new.mat.out.b,infect.new.mat.out.b[,"AI.b"]>=1)
    AI.b<-AI.mat.out.b[,"ids.b"]
    AI.mat.out.b<-as.data.frame(AI.mat.out.b)
    
    AI.mat.out.w<-subset(infect.new.mat.out.w,infect.new.mat.out.w[,"AI.w"]>=1)
    AI.w<-AI.mat.out.w[,"ids.w"]
    AI.mat.out.w<-as.data.frame(AI.mat.out.w)

    #Calculate the transmissions
    #Include PrEP
    #Select those that had AI and are negative
    infect.new.mat.out.b<-subset(infect.new.mat.out.b,infect.new.mat.out.b[,"AI.b"]>=1 & infect.new.mat.out.b[,"status.b"]==0)
    infected.b<-NULL
    infect.new.mat.out.b<-as.data.frame(infect.new.mat.out.b)
    
    infect.new.mat.out.w<-subset(infect.new.mat.out.w,infect.new.mat.out.w[,"AI.w"]>=1 & infect.new.mat.out.w[,"status.w"]==0)
    infected.w<-NULL
    infect.new.mat.out.w<-as.data.frame(infect.new.mat.out.w)
    
    

    if (length(infect.new.mat.out.b$infect.new.b)>0){
      
    infect.new.mat.out.n.b<-subset(infect.new.mat.out.b, infect.new.mat.out.b[,"prepStat.b"] == 1 & infect.new.mat.out.b[,"prepClass.b"] == "n")
    infect.new.mat.out.n.b$infect.new.b<-sapply(infect.new.mat.out.n.b$ids.b, function(x) 
      max(rbinom(infect.new.mat.out.n.b$AI.b[infect.new.mat.out.n.b$ids.b==x],1,heat.trans.b*(1 - pce[1])),0)) 
    
    infect.new.mat.out.l.b<-subset(infect.new.mat.out.b, infect.new.mat.out.b[,"prepStat.b"] == 1 & infect.new.mat.out.b[,"prepClass.b"] == "l")
    infect.new.mat.out.l.b$infect.new.b<-sapply(infect.new.mat.out.l.b$ids.b, function(x) 
      max(rbinom(infect.new.mat.out.l.b$AI.b[infect.new.mat.out.l.b$ids.b==x],1,heat.trans.b*(1 - pce[2])),0))    
      
    infect.new.mat.out.m.b<-subset(infect.new.mat.out.b, infect.new.mat.out.b[,"prepStat.b"] == 1 & infect.new.mat.out.b[,"prepClass.b"] == "m")
    infect.new.mat.out.m.b$infect.new.b<-sapply(infect.new.mat.out.m.b$ids.b, function(x) 
      max(rbinom(infect.new.mat.out.m.b$AI.b[infect.new.mat.out.m.b$ids.b==x],1,heat.trans.b*(1 - pce[3])),0))   
      
    infect.new.mat.out.h.b<-subset(infect.new.mat.out.b, infect.new.mat.out.b[,"prepStat.b"] == 1 & infect.new.mat.out.b[,"prepClass.b"] == "h")
    infect.new.mat.out.h.b$infect.new.b<-sapply(infect.new.mat.out.h.b$ids.b, function(x) 
      max(rbinom(infect.new.mat.out.h.b$AI.b[infect.new.mat.out.h.b$ids.b==x],1,heat.trans.b*(1 - pce[4])),0))   
    
    infect.new.mat.out.none.b<-subset(infect.new.mat.out.b,infect.new.mat.out.b[,"prepStat.b"] == 0) 
    infect.new.mat.out.none.b$infect.new.b<-sapply(infect.new.mat.out.none.b$ids.b, function(x) 
      max(rbinom(infect.new.mat.out.none.b$AI.b[infect.new.mat.out.none.b$ids.b==x],1,heat.trans.b),0))   
    
    
    infect.new.mat.out.b<-rbind(infect.new.mat.out.l.b,infect.new.mat.out.m.b,infect.new.mat.out.h.b,infect.new.mat.out.none.b)
    infected.b<-as.vector(infect.new.mat.out.b$ids.b[infect.new.mat.out.b$infect.new.b==1])
    }
    
    if (length(infect.new.mat.out.w$infect.new.w)>0){
      
      infect.new.mat.out.n.w<-subset(infect.new.mat.out.w, infect.new.mat.out.w[,"prepStat.w"] == 1 & infect.new.mat.out.w[,"prepClass.w"] == "n")
      infect.new.mat.out.n.w$infect.new.w<-sapply(infect.new.mat.out.n.w$ids.w, function(x) 
        max(rbinom(infect.new.mat.out.n.w$AI.w[infect.new.mat.out.n.w$ids.w==x],1,heat.trans.w*(1 - pce[1])),0)) 
      
      infect.new.mat.out.l.w<-subset(infect.new.mat.out.w, infect.new.mat.out.w[,"prepStat.w"] == 1 & infect.new.mat.out.w[,"prepClass.w"] == "l")
      infect.new.mat.out.l.w$infect.new.w<-sapply(infect.new.mat.out.l.w$ids.w, function(x) 
        max(rbinom(infect.new.mat.out.l.w$AI.w[infect.new.mat.out.l.w$ids.w==x],1,heat.trans.w*(1 - pce[2])),0))    
      
      infect.new.mat.out.m.w<-subset(infect.new.mat.out.w, infect.new.mat.out.w[,"prepStat.w"] == 1 & infect.new.mat.out.w[,"prepClass.w"] == "m")
      infect.new.mat.out.m.w$infect.new.w<-sapply(infect.new.mat.out.m.w$ids.w, function(x) 
        max(rbinom(infect.new.mat.out.m.w$AI.w[infect.new.mat.out.m.w$ids.w==x],1,heat.trans.w*(1 - pce[3])),0))   
      
      infect.new.mat.out.h.w<-subset(infect.new.mat.out.w, infect.new.mat.out.w[,"prepStat.w"] == 1 & infect.new.mat.out.w[,"prepClass.w"] == "h")
      infect.new.mat.out.h.w$infect.new.w<-sapply(infect.new.mat.out.h.w$ids.w, function(x) 
        max(rbinom(infect.new.mat.out.h.w$AI.w[infect.new.mat.out.h.w$ids.w==x],1,heat.trans.w*(1 - pce[4])),0))   
      
      infect.new.mat.out.none.w<-subset(infect.new.mat.out.w,infect.new.mat.out.w[,"prepStat.w"] == 0) 
      infect.new.mat.out.none.w$infect.new.w<-sapply(infect.new.mat.out.none.w$ids.w, function(x) 
        max(rbinom(infect.new.mat.out.none.w$AI.w[infect.new.mat.out.none.w$ids.w==x],1,heat.trans.w),0))   
      
      
      infect.new.mat.out.w<-rbind(infect.new.mat.out.l.w,infect.new.mat.out.m.w,infect.new.mat.out.h.w,infect.new.mat.out.none.w)
      infected.w<-as.vector(infect.new.mat.out.w$ids.w[infect.new.mat.out.w$infect.new.w==1])
    }
    
      
  # Update attributes
    AI<-as.numeric(AI.b)
 if (length(AI) >= 1){
     dat$attr$AI.time[AI]<-ifelse(dat$attr$everAI[AI]==0,at,dat$attr$AI.time[AI])
     dat$attr$everAI[AI]<-1
     dat$attr$AI.adult.count[AI]<-dat$attr$AI.adult.count[AI] + as.numeric(AI.mat.out.b$AI.b[AI.mat.out.b$ids.b==AI])
     dat$attr$AI.adult.count.t[AI]<- as.numeric(AI.mat.out.b$AI.b[AI.mat.out.b$ids.b==AI])
}
     
     AI<-as.numeric(AI.w)
 if (length(AI) >= 1){
       dat$attr$AI.time[AI]<-ifelse(dat$attr$everAI[AI]==0,at,dat$attr$AI.time[AI])
       dat$attr$everAI[AI]<-1
       dat$attr$AI.adult.count[AI]<-dat$attr$AI.adult.count[AI] + as.numeric(AI.mat.out.w$AI.w[AI.mat.out.w$ids.w==AI])    
       dat$attr$AI.adult.count.t[AI]<- as.numeric(AI.mat.out.w$AI.w[AI.mat.out.w$ids.w==AI])
 }
    infected<-c(as.numeric(infected.b),as.numeric(infected.w))
 if (length(infected) >= 1){
    dat$attr$status[infected] <- 1
    dat$attr$inf.time[infected] <- at
    dat$attr$vl[infected] <- 0
    dat$attr$stage[infected] <- "AR"
    dat$attr$stage.time[infected] <- 0
    dat$attr$diag.status[infected] <- 0
    dat$attr$tx.status[infected] <- 0
      
    
    dat$attr$infector[infected] <- "heat"
    dat$attr$inf.role[infected] <- "heat"
    dat$attr$inf.type[infected] <- "heat"
    dat$attr$inf.diag[infected] <- "heat" 
    dat$attr$inf.tx[infected] <- "heat"
    dat$attr$inf.stage[infected] <- "heat"
    
    dat$attr$cum.time.on.tx[infected] <- 0
    dat$attr$cum.time.off.tx[infected] <- 0
    
    
    # Summary Output
    dat$epi$incid[at] <- dat$epi$incid[at] + length(infected)
    dat$epi$incid.B[at] <- dat$epi$incid.B[at] + sum(race[infected] == "B")
    dat$epi$incid.W[at] <- dat$epi$incid.W[at] + sum(race[infected] == "W")
    dat$epi$incid.heat[at] <- length(infected)
     }
    
    if (length(infected) < 1){ dat$epi$incid.heat[at] <-0}
  
  
    
  return(dat)
}





