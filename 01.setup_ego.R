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
load("~/SHAMP/egonet/data/nsfg.egodata.rda")
str(nsfg.egodata)


data.params<-list ()

#change the age range to 18-59.
nsfg.egodata$egos$age<-sample(18:59, size=length(nsfg.egodata$egos$age), replace = TRUE, prob = NULL)
nsfg.egodata$egos$sqrt.age<-sqrt(nsfg.egodata$egos$age)



new_data<-input_shamp(nsfg.egodata, data.params, immigration=TRUE, msm.msmf=FALSE)
data.params<-as.list(new_data[1])


##Make the three ergm.ego objects.
ego.obj_m<-as.egodata(new_data[[2]]$egos,alters=new_data[[2]]$altersMain,egoIDcol="ego")
ego.obj_c<-as.egodata(new_data[[2]]$egos,alters=new_data[[2]]$altersCasual,egoIDcol="ego")
ego.obj_i<-as.egodata(new_data[[2]]$egos,alters=new_data[[2]]$altersOT,egoIDcol="ego")



##Fit the models.


fit.m.ego<-ergm.ego(ego.obj_m ~edges + nodematch("sex", diff=FALSE) + nodefactor("race",base=5)
                    +nodefactor("deg.pers.c",base=1),
                    control=control.ergm.ego(ppopsize=100000, stats.est="asymptotic",
                    ergm.control = control.ergm(MCMC.interval=5000,
                                                    MCMC.samplesize=5000,
                                                    MCMC.burnin = 5000,
                                                    MPLE.max.dyad.types = 1e6,
                                                    init.method = "MPLE",
                                                    MCMLE.maxit = 50)))


fit.c.ego<-ergm.ego(ego.obj_c ~edges + nodematch("sex", diff=FALSE) + nodefactor("race",base=5)
                    +nodefactor("deg.main.c",base=1),
                    control=control.ergm.ego(ppopsize=100000, stats.est="asymptotic",
                    ergm.control = control.ergm(MCMC.interval=5000,
                                                MCMC.samplesize=5000,
                                                MCMC.burnin = 5000,
                                                MPLE.max.dyad.types = 1e6,
                                                init.method = "MPLE",
                                                MCMLE.maxit = 50)))


fit.i.ego<-ergm.ego(ego.obj_i ~edges + nodematch("sex", diff=FALSE) + nodefactor("race",base=5)
                + nodefactor("deg.main.c",base=1) + nodefactor("deg.pers.c",base=1),
                control=control.ergm.ego(ppopsize=100000, stats.est="asymptotic",
                set.control.ergm = control.ergm(MCMC.interval=5000,
                                                MCMC.burnin = 5000,
                                                MCMC.samplesize=5000,
                                                MPLE.max.dyad.types = 1e6,
                                                init.method = "MPLE",
                                                MCMLE.maxit = 100)))





#Create additional required elements for est
# Time unit for simulation, relative to 1 day
time.unit <- 7
method<-1

#Simulation size.
sim.size<-114435


# Mean durations
diss_m = ~offset(edges) 

diss_c = ~offset(edges)

#Mortality
ages <- 18:59
age.unit <- 52

asmr.B.f <- c(rep(0, 17),
            1-(1-c(rep(0.000405376, 12),
                   rep(0.000661066, 10),
                   rep(0.001378053, 10),
                   rep(0.003065837, 10)))^(1 / age.unit),
            1)
asmr.BI.f <- c(rep(0, 17),
             1-(1-c(rep(0.000405376, 12),
                    rep(0.000661066, 10),
                    rep(0.001378053, 10),
                    rep(0.003065837, 10)))^(1 / age.unit),
             1)

asmr.H.f <- c(rep(0, 17),
            1-(1-c(rep(0.000405376, 12),
                   rep(0.000661066, 10),
                   rep(0.001378053, 10),
                   rep(0.003065837, 10)))^(1/age.unit),
            1)
asmr.HI.f <- c(rep(0, 17),
             1-(1-c(rep(0.000405376, 12),
                    rep(0.000661066, 10),
                    rep(0.001378053, 10),
                    rep(0.003065837, 10)))^(1/age.unit),
             1)
asmr.W.f <- c(rep(0, 17),
            1-(1-c(rep(0.000405376, 12),
                   rep(0.000661066, 10),
                   rep(0.001378053, 10),
                   rep(0.003065837, 10)))^(1/age.unit),
            1)


asmr.B.m <- c(rep(0, 17),
              1-(1-c(rep(0.000853417, 12),
                     rep(0.001084014, 10),
                     rep(0.001982864, 10),
                     rep(0.005400669, 10)))^(1 / age.unit),
              1)
asmr.BI.m <- c(rep(0, 17),
               1-(1-c(rep(0.000853417, 12),
                      rep(0.001084014, 10),
                      rep(0.001982864, 10),
                      rep(0.005400669, 10)))^(1 / age.unit),
               1)

asmr.H.m <- c(rep(0, 17),
              1-(1-c(rep(0.000853417, 12),
                     rep(0.001084014, 10),
                     rep(0.001982864, 10),
                     rep(0.005400669, 10)))^(1/age.unit),
              1)
asmr.HI.m <- c(rep(0, 17),
               1-(1-c(rep(0.000853417, 12),
                      rep(0.001084014, 10),
                      rep(0.001982864, 10),
                      rep(0.005400669, 10)))^(1/age.unit),
               1)
asmr.W.m <- c(rep(0, 17),
              1-(1-c(rep(0.000853417, 12),
                     rep(0.001084014, 10),
                     rep(0.001982864, 10),
                     rep(0.005400669, 10)))^(1/age.unit),
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

target.stats_m<-as.numeric(fit.m.ego$target.stats[2:length(fit.m.ego$target.stats)])
target.stats_c<-as.numeric(fit.c.ego$target.stats[2:length(fit.c.ego$target.stats)])
target.stats_i<-as.numeric(fit.i.ego$target.stats[2:length(fit.i.ego$target.stats)])

target.stats.names_m<- names(fit.m.ego$target.stats[2:length(fit.m.ego$target.stats)])
target.stats.names_c<- names(fit.c.ego$target.stats[2:length(fit.c.ego$target.stats)])
target.stats.names_i<- names(fit.i.ego$target.stats[2:length(fit.i.ego$target.stats)])


coef.form.crude_m<-fit.m.ego$coef[2:length(fit.m.ego$coef)]
coef.form_m<-coef.form.crude_m
coef.form_m[1]<-coef.form_m[1]+fit.m.ego$coef[1]
coef.form_m[1]<- coef.form_m[1]- coef.diss_m$coef.adj
constraints_m <- ~.

coef.form.crude_c<-fit.c.ego$coef[2:length(fit.c.ego$coef)]
coef.form_c<-coef.form.crude_c
coef.form_c[1]<-coef.form_c[1]+fit.c.ego$coef[1]
coef.form_c[1]<- coef.form_c[1]- coef.diss_c$coef.adj
constraints_c <- ~.

coef.form.crude_i<-fit.i.ego$coef[2:length(fit.i.ego$coef)]
coef.form_i<-coef.form.crude_i
coef.form_i[1]<-coef.form_i[1]+fit.i.ego$coef[1]
coef.form_i[1]<- coef.form_i[1]
constraints_i <- ~.



nw<-network.initialize(sim.size, directed = FALSE, hyper = FALSE, loops = FALSE,
                       multiple = FALSE, bipartite = FALSE)

fit.m <- list(fit= fit.m.ego, formation=fit.m.ego$formula, target.stats= target.stats_m,
            target.stats.names= target.stats.names_m, coef.form = coef.form_m,
            coef.form.crude= coef.form.crude_m, coef.diss=coef.diss_m, constraints= constraints_m,
            edapprox=TRUE)

fit.c <- list(fit= fit.c.ego, formation=fit.c.ego$formula, target.stats= target.stats_c,
              target.stats.names= target.stats.names_c, coef.form = coef.form_c,
              coef.form.crude= coef.form.crude_c, coef.diss=coef.diss_c, constraints= constraints_c,
              edapprox=TRUE)

fit.i <- list(fit= fit.i.ego, formation=fit.i.ego$formula, target.stats= target.stats_i,
              target.stats.names= target.stats.names_i, coef.form = coef.form_i,
              coef.form.crude= coef.form.crude_i, coef.diss=coef.diss_i, constraints= constraints_i,
              edapprox=TRUE)

param <- param_shamp(data.params)
init <- init_shamp()
control <- control_shamp(nsteps = 520)
est <- list(fit.m, fit.c, fit.i)
save(est, file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fit.rda")

load(file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fit.rda")

sim<-netsim(est, param, init, control)

save(sim, file = "~/EpiModelHIV_shamp_modeling/scenarios/sim.rda")
demog.table<-as.data.frame(rbind(demog4,demog52,demog104,demog208,demog312,demog416,demog520))
save(demog.table, file = "~/EpiModelHIV_shamp_modeling/scenarios/demogs.rda")



library(xlsx) #load the package
write.xlsx(x = demog.table, file = "~/EpiModelHIV_shamp_modeling/scenarios/demogs.xlsx",
           sheetName = "SHAMP Demog", row.names = FALSE)

netsim(est, param, init, control)
