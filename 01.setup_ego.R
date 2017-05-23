##Running network simulations with EPIModelHIV and ergm.ego.

#Load required packages.
library(EpiModelHIV)
library(EpiModelHPC)
#devtools::install_github("statnet/ergm.ego-private", auth_token = )
#devtools::install_github("statnet/tergmLite")
suppressPackageStartupMessages(library(ergm.ego))

#Load the ego and alter data frames and look at them
load("~/SHAMP/egonet/data/nsfg.egodata.rda")
str(nsfg.egodata)
data.params<-list ()

#change the age range to 18-59.
nsfg.egodata$egos$age<-sample(18:59, size=length(nsfg.egodata$egos$age), replace = TRUE, prob = NULL)
nsfg.egodata$egos$sqrt.age<-sqrt(nsfg.egodata$egos$age)

##TEMP get rid of race "O".
nsfg.egodata$egos$race<-ifelse(nsfg.egodata$egos$race=="O","W",nsfg.egodata$egos$race)
nsfg.egodata$altersMain$race<-ifelse(nsfg.egodata$altersMain$race=="O","W",nsfg.egodata$altersMain$race)
nsfg.egodata$altersCasual$race<-ifelse(nsfg.egodata$altersCasual$race=="O","W",nsfg.egodata$altersCasual$race)
nsfg.egodata$altersOT$race<-ifelse(nsfg.egodata$altersOT$race=="O","W",nsfg.egodata$altersOT$race)

##CHANGE NAMES TO deg.pers and deg.main.
nsfg.egodata$egos$deg.pers<-as.numeric(nsfg.egodata$egos$casual)
nsfg.egodata$egos$deg.main<-as.numeric(nsfg.egodata$egos$main)


new_data<-input_shamp(nsfg.egodata, data.params, immigration=TRUE, msm.msmf=TRUE)
data.params<-as.list(new_data[1])




##Make the three ergm.ego objects.
ego.obj_m<-as.egodata(new_data[[2]]$egos,alters=new_data[[2]]$altersMain,egoIDcol="ego")
ego.obj_c<-as.egodata(new_data[[2]]$egos,alters=new_data[[2]]$altersCasual,egoIDcol="ego")
ego.obj_i<-as.egodata(new_data[[2]]$egos,alters=new_data[[2]]$altersOT,egoIDcol="ego")




##This fits.
fit.m<-ergm.ego(ego.obj_m ~edges + nodematch("sex", diff=TRUE) + nodefactor("race",base=1)
                +nodefactor("deg.pers",base=1),
                control=control.ergm.ego(ppopsize=50000, stats.est="bootstrap"),
                set.control.ergm = control.ergm(MCMC.interval=5000,
                                                MCMC.samplesize=5000,
                                                MPLE.max.dyad.types = 1e10,
                                                init.method = "MPLE",
                                                MCMLE.maxit = 200,
                                                SAN.maxit=100,
                                                SAN.burnin.times = 500))


fit.c<-ergm.ego(ego.obj_c ~edges + nodematch("sex", diff=TRUE) + nodefactor("race",base=1)
                +nodefactor("deg.main",base=1),
                control=control.ergm.ego(ppopsize=50000, stats.est="bootstrap"))

fit.i<-ergm.ego(ego.obj_i ~edges + nodematch("sex", diff=TRUE) + nodefactor("race",base=1)
                + nodefactor("deg.pers",base=1),
                control=control.ergm.ego(ppopsize=50000, stats.est="bootstrap"))


#Create additional required elements for est
# Time unit for simulation, relative to 1 day
time.unit <- 7
method<-1

#Simulation size.
sim.size<-38362


# Mean durations
diss_m = ~offset(edges) 

diss_c = ~offset(edges)

#Mortality
ages <- 18:59
age.unit <- 52

asmr.B.f <- c(rep(0, 17),
            1-(1-c(rep(0.002, 7),
                   rep(0.001, 10),
                   rep(0.003, 25)))^(1 / age.unit),
            1)
asmr.BI.f <- c(rep(0, 17),
             1-(1-c(rep(0.002, 7),
                    rep(0.001, 10),
                    rep(0.003, 25)))^(1 / age.unit),
             1)

asmr.H.f <- c(rep(0, 17),
            1-(1-c(rep(0.002, 7),
                   rep(0.001, 10),
                   rep(0.003, 25)))^(1/age.unit),
            1)
asmr.HI.f <- c(rep(0, 17),
             1-(1-c(rep(0.002, 7),
                    rep(0.001, 10),
                    rep(0.003, 25)))^(1/age.unit),
             1)
asmr.W.f <- c(rep(0, 17),
            1-(1-c(rep(0.002, 7),
                   rep(0.001, 10),
                   rep(0.003, 25)))^(1/age.unit),
            1)


asmr.B.m <- c(rep(0, 17),
              1-(1-c(rep(0.002, 7),
                     rep(0.001, 10),
                     rep(0.003, 25)))^(1 / age.unit),
              1)
asmr.BI.m <- c(rep(0, 17),
               1-(1-c(rep(0.002, 7),
                      rep(0.001, 10),
                      rep(0.003, 25)))^(1 / age.unit),
               1)

asmr.H.m <- c(rep(0, 17),
              1-(1-c(rep(0.002, 7),
                     rep(0.001, 10),
                     rep(0.003, 25)))^(1/age.unit),
              1)
asmr.HI.m <- c(rep(0, 17),
               1-(1-c(rep(0.002, 7),
                      rep(0.001, 10),
                      rep(0.003, 25)))^(1/age.unit),
               1)
asmr.W.m <- c(rep(0, 17),
              1-(1-c(rep(0.002, 7),
                     rep(0.001, 10),
                     rep(0.003, 25)))^(1/age.unit),
              1)
exp.mort <- (mean(asmr.B.f[ages]) + mean(asmr.BI.f[ages]) + mean(asmr.H.f[ages])
             + mean(asmr.HI.f[ages]) + mean(asmr.W.f[ages]) + mean(asmr.B.m[ages]) 
             + mean(asmr.BI.m[ages]) + mean(asmr.H.m[ages])
             + mean(asmr.HI.m[ages]) + mean(asmr.W.m[ages]) ) / 10

coef.diss_m <- dissolution_coefs(dissolution = diss_m,
                               duration = data.params[[1]]$durs_m / time.unit,
                               d.rate = exp.mort)

coef.diss_c <- dissolution_coefs(dissolution = diss_c,
                               duration = data.params[[1]]$durs_c / time.unit,
                               d.rate = exp.mort)


coef.diss_i <-dissolution_coefs(~offset(edges), 1)

target.stats_m<-as.numeric(fit.m$target.stats[2:length(fit.m$target.stats)])
target.stats_c<-as.numeric(fit.c$target.stats[2:length(fit.c$target.stats)])
target.stats_i<-as.numeric(fit.i$target.stats[2:length(fit.i$target.stats)])

target.stats.names_m<- names(fit.m$target.stats[2:length(fit.m$target.stats)])
target.stats.names_c<- names(fit.c$target.stats[2:length(fit.c$target.stats)])
target.stats.names_i<- names(fit.i$target.stats[2:length(fit.i$target.stats)])

coef.form.crude_m<-fit.m$coef[2:length(fit.m$coef)]
coef.form_m<-coef.form.crude_m
coef.form_m[1]<- coef.form_m[1]- coef.diss_m$coef.adj
constraints_m <- ~.

coef.form.crude_c<-fit.c$coef[2:length(fit.c$coef)]
coef.form_c<-coef.form.crude_c
coef.form_c[1]<- coef.form_c[1]- coef.diss_c$coef.adj
constraints_c <- ~.

coef.form.crude_i<-fit.i$coef[2:length(fit.i$coef)]
coef.form_i<-coef.form.crude_i
coef.form_i[1]<- coef.form_i[1]
constraints_i <- ~.

nw<-network.initialize(sim.size, directed = FALSE, hyper = FALSE, loops = FALSE,
                       multiple = FALSE, bipartite = FALSE)

fit.m <- list(fit= fit.m, formation=fit.m$formula, target.stats= target.stats_m,
            target.stats.names= target.stats.names_m, coef.form = coef.form_m,
            coef.form.crude= coef.form.crude_m, coef.diss=coef.diss_m, constraints= constraints_m,
            edapprox=TRUE)

fit.c <- list(fit= fit.c, formation=fit.c$formula, target.stats= target.stats_c,
              target.stats.names= target.stats.names_c, coef.form = coef.form_c,
              coef.form.crude= coef.form.crude_c, coef.diss=coef.diss_c, constraints= constraints_c,
              edapprox=TRUE)

fit.i <- list(fit= fit.i, formation=fit.i$formula, target.stats= target.stats_i,
              target.stats.names= target.stats.names_i, coef.form = coef.form_i,
              coef.form.crude= coef.form.crude_i, coef.diss=coef.diss_i, constraints= constraints_i,
              edapprox=TRUE)

param <- param_shamp(data.params)
init <- init_shamp()
control <- control_shamp(nsteps = 5)
est <- list(fit.m, fit.c, fit.i)
save(est, file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fit.rda")

load(file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fit.rda")

sim<-netsim(est, param, init, control)
save(sim, file = "~/EpiModelHIV_shamp_modeling/scenarios/sim.rda")
save(demog4, file = "~/EpiModelHIV_shamp_modeling/scenarios/demogs.rda")
save(c(demog4,demog52,demog104,demog1040), file = "~/EpiModelHIV_shamp_modeling/scenarios/demogs.rda")


#sim <- netsim(est, param, init, control)
