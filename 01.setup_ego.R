##Running network simulations with EPIModelHIV and ergm.ego.

#Load required packages.
library(EpiModelHIV)
library(EpiModelHPC)
#devtools::install_github("statnet/ergm.ego-private",ref="3.7-compat", auth_token ="")
#devtools::install_github("statnet/tergmLite")
#devtools::install_github("statnet/ergm.ego")
library(ergm.ego)
library(tergmLite)

#Load the ego and alter data frames and look at them
load("~/SHAMP/egonet/data/all.egodata.rda")
str(all.egodata)


data.params<-list ()


new_data<-input_shamp(all.egodata, data.params, immigration=TRUE, msm.msmf=FALSE)
data.params<-as.list(new_data[1])


##Make the three ergm.ego objects.
ego.obj_c<-as.egodata(new_data[[2]]$egos,alters=new_data[[2]]$altersCohab,egoIDcol="ego", egoWt=new_data[[2]]$egos$weight)
ego.obj_p<-as.egodata(new_data[[2]]$egos,alters=new_data[[2]]$altersPers,egoIDcol="ego", egoWt=new_data[[2]]$egos$weight)
ego.obj_i<-as.egodata(new_data[[2]]$egos,alters=new_data[[2]]$altersOT,egoIDcol="ego", egoWt=new_data[[2]]$egos$weight)



###cohab partnership network.

fit.c<-ergm.ego(ego.obj_c ~edges + 
                    nodefactor("race.sex.pers",base=9) + 
                    nodefactor("race",base=5) +
                    nodematch("race",diff=TRUE) +
                    nodefactor("agecat", base=1) + 
                    absdiff("sqrt.age.adj") + 
                    offset(nodematch("sex", diff=FALSE)),
                  offset.coef = -Inf,
                  control=control.ergm.ego(ppopsize=50000, stats.est="asymptotic",
                                           ergm.control = control.ergm(MCMC.interval=6000,
                                                                       MCMC.samplesize=6000,
                                                                       MCMC.burnin = 6000,
                                                                       MPLE.max.dyad.types = 1e5,
                                                                       init.method = "MPLE",
                                                                       MCMLE.maxit = 250)))

summary(fit.c)


####Casual partnership network.


fit.p<-ergm.ego(ego.obj_p ~edges + 
                    nodefactor("race.sex.cohab",base=9) + 
                    nodefactor("race",base=5) + 
                    nodematch("race",diff=TRUE) +
                    nodefactor("agecat", base=1) + 
                    absdiff("sqrt.age.adj") + 
                    offset(nodematch("sex", diff=FALSE)),
                  offset.coef = -Inf,
                  control=control.ergm.ego(ppopsize=50000, stats.est="asymptotic",
                                           ergm.control = control.ergm(MCMC.interval=5000,
                                                                       MCMC.samplesize=5000,
                                                                       MCMC.burnin = 5000,
                                                                       MPLE.max.dyad.types = 1e5,
                                                                       init.method = "MPLE",
                                                                       MCMLE.maxit = 250)))

summary(fit.p)


######One time partnerships.

summary(ego.obj_i ~edges + 
          nodefactor("race", base=0) +
          nodematch("race",diff=TRUE) +
          nodefactor("agecat", base=0) + 
          nodefactor("deg.cohab.c",base=0) +
          nodefactor("deg.pers.c",base=0))

fit.i<-ergm.ego(ego.obj_i ~edges + 
                    nodefactor("race", base=5) +
                    nodematch("race",diff=TRUE) +
                    nodefactor("agecat", base=c(3,4)) + 
                    nodefactor("deg.cohab.c",base=1) +
                    nodefactor("deg.pers.c",base=1) +
                    offset(nodematch("sex", diff=FALSE)),
                  offset.coef = -Inf,
                  verbose=TRUE,
                  control=control.ergm.ego(ppopsize=50000, stats.est="asymptotic",
                                           ergm.control = control.ergm(SAN.maxit=50,
                                                                       SAN.burnin.times=100,
                                                                       MCMC.interval=7000,
                                                                       MCMC.samplesize=7000,
                                                                       MCMC.burnin = 7000,
                                                                       MPLE.max.dyad.types = 1e4,
                                                                       init.method = "MPLE",
                                                                       MCMLE.maxit = 250)))

summary(fit.i)



modelfits <- list(fit.c, fit.p, fit.i)
save(modelfits, file = "~/EpiModelHIV_shamp_modeling/scenarios/est/modelfits.rda")



#Create additional required elements for est
# Time unit for simulation, relative to 1 day
time.unit <- 7
method<-1

#Simulation size.
sim.size<-36368


# Mean durations
diss_c = ~offset(edges) 

diss_p = ~offset(edges)

#Mortality
ages <- 18:60
age.unit <- 52

asmr.B.f <- c(rep(0, 17),
            1-(1-c(rep(0.000405376, 12),
                   rep(0.000661066, 10),
                   rep(0.001378053, 10),
                   rep(0.003065837, 11)))^(1 / age.unit),
            1)
asmr.BI.f <- c(rep(0, 17),
             1-(1-c(rep(0.000405376, 12),
                    rep(0.000661066, 10),
                    rep(0.001378053, 10),
                    rep(0.003065837, 11)))^(1 / age.unit),
             1)

asmr.H.f <- c(rep(0, 17),
            1-(1-c(rep(0.000405376, 12),
                   rep(0.000661066, 10),
                   rep(0.001378053, 10),
                   rep(0.003065837, 11)))^(1/age.unit),
            1)
asmr.HI.f <- c(rep(0, 17),
             1-(1-c(rep(0.000405376, 12),
                    rep(0.000661066, 10),
                    rep(0.001378053, 10),
                    rep(0.003065837, 11)))^(1/age.unit),
             1)
asmr.W.f <- c(rep(0, 17),
            1-(1-c(rep(0.000405376, 12),
                   rep(0.000661066, 10),
                   rep(0.001378053, 10),
                   rep(0.003065837, 11)))^(1/age.unit),
            1)


asmr.B.m <- c(rep(0, 17),
              1-(1-c(rep(0.000853417, 12),
                     rep(0.001084014, 10),
                     rep(0.001982864, 10),
                     rep(0.005400669, 11)))^(1 / age.unit),
              1)
asmr.BI.m <- c(rep(0, 17),
               1-(1-c(rep(0.000853417, 12),
                      rep(0.001084014, 10),
                      rep(0.001982864, 10),
                      rep(0.005400669, 11)))^(1 / age.unit),
               1)

asmr.H.m <- c(rep(0, 17),
              1-(1-c(rep(0.000853417, 12),
                     rep(0.001084014, 10),
                     rep(0.001982864, 10),
                     rep(0.005400669, 11)))^(1/age.unit),
              1)
asmr.HI.m <- c(rep(0, 17),
               1-(1-c(rep(0.000853417, 12),
                      rep(0.001084014, 10),
                      rep(0.001982864, 10),
                      rep(0.005400669, 11)))^(1/age.unit),
               1)
asmr.W.m <- c(rep(0, 17),
              1-(1-c(rep(0.000853417, 12),
                     rep(0.001084014, 10),
                     rep(0.001982864, 10),
                     rep(0.005400669, 11)))^(1/age.unit),
              1)




exp.mort <- (mean(asmr.B.f[ages]) + mean(asmr.BI.f[ages]) + mean(asmr.H.f[ages])
             + mean(asmr.HI.f[ages]) + mean(asmr.W.f[ages]) + mean(asmr.B.m[ages]) 
             + mean(asmr.BI.m[ages]) + mean(asmr.H.m[ages])
             + mean(asmr.HI.m[ages]) + mean(asmr.W.m[ages]) ) / 10

coef.diss_c <- dissolution_coefs(dissolution = diss_c,
                               duration = data.params[[1]]$durs_c / time.unit,
                               d.rate = exp.mort)

coef.diss_p <- dissolution_coefs(dissolution = diss_p,
                               duration = data.params[[1]]$durs_p / time.unit,
                               d.rate = exp.mort)


coef.diss_i <-dissolution_coefs(~offset(edges), 1)

target.stats_c<-as.numeric(fit.c$target.stats[2:length(fit.c$target.stats)])
target.stats_p<-as.numeric(fit.p$target.stats[2:length(fit.p$target.stats)])
target.stats_i<-as.numeric(fit.i$target.stats[2:length(fit.i$target.stats)])

target.stats.names_c<- names(fit.c$target.stats[2:length(fit.c$target.stats)])
target.stats.names_p<- names(fit.p$target.stats[2:length(fit.p$target.stats)])
target.stats.names_i<- names(fit.i$target.stats[2:length(fit.i$target.stats)])


coef.form.crude_c<-fit.c$coef[2:length(fit.c$coef)]
coef.form_c<-coef.form.crude_c
coef.form_c[1]<-coef.form_c[1]+fit.c$coef[1]
coef.form_c[1]<- coef.form_c[1]- coef.diss_c$coef.adj
constraints_c <- ~.

coef.form.crude_p<-fit.p$coef[2:length(fit.p$coef)]
coef.form_p<-coef.form.crude_p
coef.form_p[1]<-coef.form_p[1]+fit.p$coef[1]
coef.form_p[1]<- coef.form_p[1]- coef.diss_p$coef.adj
constraints_p <- ~.

coef.form.crude_i<-fit.i$coef[2:length(fit.i$coef)]
coef.form_i<-coef.form.crude_i
coef.form_i[1]<-coef.form_i[1]+fit.i$coef[1]
coef.form_i[1]<- coef.form_i[1]
coef.form_i[1]<-coef.form_i[1]+ -log(52)
constraints_i <- ~.



nw<-network.initialize(sim.size, directed = FALSE, hyper = FALSE, loops = FALSE,
                       multiple = FALSE, bipartite = FALSE)

fit.c.s <- list(fit= fit.c, formation=fit.c$formula, target.stats= target.stats_c,
            target.stats.names= target.stats.names_c, coef.form = coef.form_c,
            coef.form.crude= coef.form.crude_c, coef.diss=coef.diss_c, constraints= constraints_c,
            edapprox=TRUE)

fit.p.s <- list(fit= fit.p, formation=fit.p$formula, target.stats= target.stats_p,
              target.stats.names= target.stats.names_p, coef.form = coef.form_p,
              coef.form.crude= coef.form.crude_p, coef.diss=coef.diss_p, constraints= constraints_p,
              edapprox=TRUE)

fit.i.s <- list(fit= fit.i, formation=fit.i$formula, target.stats= target.stats_i,
              target.stats.names= target.stats.names_i, coef.form = coef.form_i,
              coef.form.crude= coef.form.crude_i, coef.diss=coef.diss_i, constraints= constraints_i,
              edapprox=TRUE)

param <- param_shamp(data.params)
init <- init_shamp()
control <- control_shamp(nsteps = 5)
est <- list(fit.c.s, fit.p.s, fit.i.s)
#save(est, file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fit.rda")

save(est, file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fitsmall.rda")
save(data.params, file = "~/EpiModelHIV_shamp_modeling/scenarios/est/data.params.rda")

#load(file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fit.rda")

load(file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fitsmall.rda")

sim<-netsim(est, param, init, control)

save(sim, file = "~/EpiModelHIV_shamp_modeling/scenarios/sim.rda")
demog.table<-as.data.frame(rbind(demog4,demog52,demog104,demog208,demog312,demog416,demog520))
save(demog.table, file = "~/EpiModelHIV_shamp_modeling/scenarios/demogs.rda")



library(xlsx) #load the package
write.xlsx(x = demog.table, file = "~/EpiModelHIV_shamp_modeling/scenarios/demogs.xlsx",
           sheetName = "SHAMP Demog", row.names = FALSE)

netsim(est, param, init, control)
