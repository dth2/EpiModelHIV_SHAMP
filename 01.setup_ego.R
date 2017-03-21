##Running network simulations with EPIModelHIV and ergm.ego.

#Load required packages.
library(EpiModelHIV)
library(EpiModelHPC)
#devtools::install_github("statnet/ergm.ego-private", auth_token = )
#devtools::install_github("statnet/tergmLite")
suppressPackageStartupMessages(library(ergm.ego))

#Load the ego and alter data frames and look at them
#load("~/SHAMP/egonet/data/nsfg.egodata.rda")
#str(nsfg.egodata)

#change the age range to 18-55.
#nsfg.egodata$egos$age<-sample(18:55, size=length(nsfg.egodata$egos$age), replace = TRUE, prob = NULL)
##Fix race and sex attibutes.

#RACE.
#nsfg.egodata$egos$race<-ifelse(nsfg.egodata$egos$race=="Non-Hispanic Black","B",
#                               ifelse(nsfg.egodata$egos$race=="Non-Hispanic Black Imigrant","BI",
#                               ifelse(nsfg.egodata$egos$race=="Hispanic","H",
#                              ifelse(nsfg.egodata$egos$race=="Non-Hispanic Other","HI",
#                               ifelse(nsfg.egodata$egos$race=="Non-Hispanic White","W",nsfg.egodata$egos$race)))))

#nsfg.egodata$altersMain$race<-ifelse(nsfg.egodata$altersMain$race=="Non-Hispanic Black","B",
#                               ifelse(nsfg.egodata$altersMain$race=="Non-Hispanic Black Imigrant","BI",
#                               ifelse(nsfg.egodata$altersMain$race=="Hispanic","H",
#                               ifelse(nsfg.egodata$altersMain$race=="Non-Hispanic Other","HI",
#                               ifelse(nsfg.egodata$altersMain$race=="Non-Hispanic White","W",nsfg.egodata$altersMain$race)))))

#nsfg.egodata$altersCasual$race<-ifelse(nsfg.egodata$altersCasual$race=="Non-Hispanic Black","B",
#                               ifelse(nsfg.egodata$altersCasual$race=="Non-Hispanic Black Imigrant","BI",
#                               ifelse(nsfg.egodata$altersCasual$race=="Hispanic","H",
#                               ifelse(nsfg.egodata$altersCasual$race=="Non-Hispanic Other","HI",
#                               ifelse(nsfg.egodata$altersCasual$race=="Non-Hispanic White","W",nsfg.egodata$altersCasual$race)))))

#nsfg.egodata$altersOT$race<-ifelse(nsfg.egodata$altersOT$race=="Non-Hispanic Black","B",
#                               ifelse(nsfg.egodata$altersOT$race=="Non-Hispanic Black Imigrant","BI",
#                               ifelse(nsfg.egodata$altersOT$race=="Hispanic","H",
#                               ifelse(nsfg.egodata$altersOT$race=="Non-Hispanic Other","HI",
#                               ifelse(nsfg.egodata$altersOT$race=="Non-Hispanic White","W",nsfg.egodata$altersOT$race)))))
#SEX.
#nsfg.egodata$egos$sex<-ifelse(nsfg.egodata$egos$sex=="male","M",
#                               ifelse(nsfg.egodata$egos$sex=="female","F",
#                                      nsfg.egodata$egos$sex))
#nsfg.egodata$altersMain$sex<-ifelse(nsfg.egodata$altersMain$sex=="male","M",
#                                     ifelse(nsfg.egodata$altersMain$sex=="female","F",
#                                            nsfg.egodata$altersMain$sex))
#nsfg.egodata$altersCasual$sex<-ifelse(nsfg.egodata$altersCasual$sex=="male","M",
#                                       ifelse(nsfg.egodata$altersCasual$sex=="female","F"
#                                              ,nsfg.egodata$altersCasual$sex))

#nsfg.egodata$altersOT$sex<-ifelse(nsfg.egodata$altersOT$sex=="male","M",
#                                   ifelse(nsfg.egodata$altersOT$sex=="female","F",
#                                    nsfg.egodata$altersOT$sex))


##Add in variable not yet added: sqrt.age, immig, immig.loc, sex.ident, role.class.
#nsfg.egodata$egos$sqrt.age<-sqrt(nsfg.egodata$egos$age)
#nsfg.egodata$egos$immig<-rep("n",length(nsfg.egodata$egos$age))
#nsfg.egodata$egos$immig.loc<-rep(0,length(nsfg.egodata$egos$age))
#nsfg.egodata$egos$sex.ident[nsfg.egodata$egos$sex=="M"] <- "msf"
#nsfg.egodata$egos$sex.ident[nsfg.egodata$egos$sex=="F"] <- "f"

#nsfg.egodata$egos$role.class[nsfg.egodata$egos$sex=="M"] <- "I"
#nsfg.egodata$egos$role.class[nsfg.egodata$egos$sex=="F"] <- "R"

##Make the three ergm.ego objects.
#ego.obj_m<-as.egodata(nsfg.egodata$egos,alters=nsfg.egodata$altersMain,egoIDcol="ego")
#ego.obj_c<-as.egodata(nsfg.egodata$egos,alters=nsfg.egodata$altersCasual,egoIDcol="ego")
#ego.obj_i<-as.egodata(nsfg.egodata$egos,alters=nsfg.egodata$altersOT,egoIDcol="ego")

##Make sure the variables for casual and main are character.

##CHANGE NAMES TO deg.pers and deg.main.
#ego.obj_m$egos$deg.pers<-as.numeric(ego.obj_m$egos$casual)
#ego.obj_m$egos$deg.main<-as.numeric(ego.obj_m$egos$main)
#ego.obj_c$egos$deg.pers<-as.numeric(ego.obj_c$egos$casual)
#ego.obj_c$egos$deg.main<-as.numeric(ego.obj_c$egos$main)
#ego.obj_i$egos$deg.pers<-as.numeric(ego.obj_i$egos$casual)
#ego.obj_i$egos$deg.main<-as.numeric(ego.obj_i$egos$main)



##For the time being we can only use attributes in both the ego and alter file.
##This fits.
#fit.m<-ergm.ego(ego.obj_m ~edges + nodematch("sex", diff=TRUE) + nodefactor("race",base=1)
#                +nodefactor("deg.pers",base=1),
#                control=control.ergm.ego(ppopsize=50000, stats.est="bootstrap"))

#fit.c<-ergm.ego(ego.obj_c ~edges + nodematch("sex", diff=TRUE) +nodefactor("race",base=1)
#                +nodefactor("deg.main",base=1),
#                control=control.ergm.ego(ppopsize=50000))

#fit.i<-ergm.ego(ego.obj_i ~edges + nodematch("sex", diff=TRUE) + nodefactor("race",base=1)
#                + nodefactor("deg.pers",base=1),
#                control=control.ergm.ego(ppopsize=50000))


#Create additional required elements for est
# Time unit for simulation, relative to 1 day
time.unit <- 7
method<-1

#Simulation size.
sim.size<-40412

#Target Demograqphics
## May need these to target birthrates  (TABLE EGO RACE BY SEX  )
#num.B.f<-
#num.BI.f
#num.H.f
#num.HI.f
#num.W.f
#num.B.f
#num.BI.f
#num.H.f
#num.HI.f
#num.W.f

# Mean durations
rates_m <- .049505
durs_m  <- 1/rates_m
diss_m = ~offset(edges)

rates_c <- .03
durs_c  <- 1/rates_c
diss_c = ~offset(edges)

#Mortality
ages <- 18:55
age.unit <- 52

asmr.B.f <- c(rep(0, 17),
            1-(1-c(rep(0.00159, 7),
                   rep(0.00225, 10),
                   rep(0.00348, 21)))^(1 / age.unit),
            1)
asmr.BI.f <- c(rep(0, 17),
             1-(1-c(rep(0.00159, 7),
                    rep(0.00225, 10),
                    rep(0.00348, 21)))^(1 / age.unit),
             1)

asmr.H.f <- c(rep(0, 17),
            1-(1-c(rep(0.00103, 7),
                   rep(0.00133, 10),
                   rep(0.00214, 21)))^(1/age.unit),
            1)
asmr.HI.f <- c(rep(0, 17),
             1-(1-c(rep(0.00103, 7),
                    rep(0.00133, 10),
                    rep(0.00214, 21)))^(1/age.unit),
             1)
asmr.W.f <- c(rep(0, 17),
            1-(1-c(rep(0.00103, 7),
                   rep(0.00133, 10),
                   rep(0.00214, 21)))^(1/age.unit),
            1)


asmr.B.m <- c(rep(0, 17),
              1-(1-c(rep(0.00159, 7),
                     rep(0.00225, 10),
                     rep(0.00348, 21)))^(1 / age.unit),
              1)
asmr.BI.m <- c(rep(0, 17),
               1-(1-c(rep(0.00159, 7),
                      rep(0.00225, 10),
                      rep(0.00348, 21)))^(1 / age.unit),
               1)

asmr.H.m <- c(rep(0, 17),
              1-(1-c(rep(0.00103, 7),
                     rep(0.00133, 10),
                     rep(0.00214, 21)))^(1/age.unit),
              1)
asmr.HI.m <- c(rep(0, 17),
               1-(1-c(rep(0.00103, 7),
                      rep(0.00133, 10),
                      rep(0.00214, 21)))^(1/age.unit),
               1)
asmr.W.m <- c(rep(0, 17),
              1-(1-c(rep(0.00103, 7),
                     rep(0.00133, 10),
                     rep(0.9, 21)))^(1/age.unit),
              1)
exp.mort <- (mean(asmr.B.f[ages]) + mean(asmr.BI.f[ages]) + mean(asmr.H.f[ages])
             + mean(asmr.HI.f[ages]) + mean(asmr.W.f[ages]) + mean(asmr.B.m[ages]) 
             + mean(asmr.BI.m[ages]) + mean(asmr.H.m[ages])
             + mean(asmr.HI.m[ages]) + mean(asmr.W.m[ages]) ) / 10

coef.diss_m <- dissolution_coefs(dissolution = diss_m,
                               duration = durs_m / time.unit,
                               d.rate = exp.mort)

coef.diss_c <- dissolution_coefs(dissolution = diss_c,
                               duration = durs_c / time.unit,
                               d.rate = 0)

coef.diss_i <-dissolution_coefs(~offset(edges), 1)

#target.stats_m<-as.numeric(fit.m$target.stats[2:length(fit.m$target.stats)])
#target.stats_c<-as.numeric(fit.c$target.stats[2:length(fit.c$target.stats)])
#target.stats_i<-as.numeric(fit.i$target.stats[2:length(fit.i$target.stats)])

#target.stats.names_m<- names(fit.m$target.stats[2:length(fit.m$target.stats)])
#target.stats.names_c<- names(fit.c$target.stats[2:length(fit.c$target.stats)])
#target.stats.names_i<- names(fit.i$target.stats[2:length(fit.i$target.stats)])

#coef.form.crude_m<-fit.m$coef[2:length(fit.m$coef)]
#coef.form_m<-coef.form.crude_m
#coef.form_m[1]<- coef.form_m[1]- coef.diss_m$coef.adj
#constraints_m <- ~.

#coef.form.crude_c<-fit.c$coef[2:length(fit.c$coef)]
#coef.form_c<-coef.form.crude_c
#coef.form_c[1]<- coef.form_c[1]- coef.diss_c$coef.adj
#constraints_c <- ~.

#coef.form.crude_i<-fit.i$coef[2:length(fit.i$coef)]
#coef.form_i<-coef.form.crude_i
#coef.form_i[1]<- coef.form_i[1]
#constraints_i <- ~.

nw<-network.initialize(sim.size, directed = FALSE, hyper = FALSE, loops = FALSE,
                       multiple = FALSE, bipartite = FALSE)

#fit.m <- list(fit= fit.m, formation=fit.m$formula, target.stats= target.stats_m,
#            target.stats.names= target.stats.names_m, coef.form = coef.form_m,
#            coef.form.crude= coef.form.crude_m, coef.diss=coef.diss_m, constraints= constraints_m,
#            edapprox=TRUE)

#fit.c <- list(fit= fit.c, formation=fit.c$formula, target.stats= target.stats_c,
#              target.stats.names= target.stats.names_c, coef.form = coef.form_c,
#              coef.form.crude= coef.form.crude_c, coef.diss=coef.diss_c, constraints= constraints_c,
#              edapprox=TRUE)

#fit.i <- list(fit= fit.i, formation=fit.i$formula, target.stats= target.stats_i,
#              target.stats.names= target.stats.names_i, coef.form = coef.form_i,
#              coef.form.crude= coef.form.crude_i, coef.diss=coef.diss_i, constraints= constraints_i,
#              edapprox=TRUE)

param <- param_shamp()
init <- init_shamp()
control <- control_shamp(nsteps = 10)
est <- list(fit.m, fit.c, fit.i)
#save(est, file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fit.rda")

load(file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fit.rda")

netsim(est, param, init, control)
#sim <- netsim(est, param, init, control)
