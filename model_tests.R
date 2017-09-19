##Testing model fits.

#Load required packages.
library(EpiModelHIV)
library(ergm.ego)
library(parallel)
np = detectCores()



load(file = "~/EpiModelHIV_SHAMP2/data/ego.obj_c.rda")
load(file = "~/EpiModelHIV_SHAMP2/data/ego.obj_p.rda")
load(file = "~/EpiModelHIV_SHAMP2/data/ego.obj_i.rda")



###cohab partnership network.

summary(ego.obj_c ~edges + 
                  nodefactor("race.sex",base=0) +
                  nodematch("race",diff=TRUE) +
                  nodefactor("agecat", base=1) + 
                  absdiff("sqrt.age.adj") + 
                  nodefactor("deg.pers.c",base=1))

            
fit.c<-ergm.ego(ego.obj_c ~edges + 
                    nodefactor("race",base=5) +
                    nodematch("race",diff=TRUE) +
                    nodefactor("agecat", base=1) + 
                    absdiff("sqrt.age.adj") + 
                    nodefactor("deg.pers.c",base=1) + 
                    offset(nodematch("sex", diff=FALSE)) + 
                    offset(concurrent()),
                    offset.coef = c(-Inf, -Inf),
                    constraints=~bd(maxout=1),
                    control=control.ergm.ego(ppopsize=50000, stats.est="asymptotic",
                                           ergm.control = control.ergm(MCMC.interval=7500,
                                                                       MCMC.samplesize=7500,
                                                                       MCMC.burnin = 7500,
                                                                       MPLE.max.dyad.types = 1e7,
                                                                       init.method = "zeros",
                                                                       MCMLE.maxit = 400,
                                                                       parallel = np, 
                                                                       parallel.type="PSOCK")))




plot(gof(fit.c, GOF="degree"))
plot(gof(fit.c, GOF="model"))


test<-simulate(fit.c)
test <- as.egodata(test)
degreedist.egodata(test, by="race.sex", prob=T)
summary(test ~ degree(0:5, by="race.sex"))

pdf(file = "XXX", height = 6, width = 12, pointsize = 16)
#par(mfrow = c(1,2), mar = c(4,3,2.5,1), mgp = c(2,.5,0))
degreedist.egodata(test, by="race.sex")
dev.off()

####Casual partnership network.

summary(ego.obj_p ~edges + 
          nodefactor("race",base=5) + 
          nodematch("race",diff=TRUE) +
          nodefactor("agecat", base=1) + 
          absdiff("sqrt.age.adj") + 
          nodefactor("deg.cohab.c",base=1) +
          concurrent(by="sex") + 
          concurrent(by="race"))

fit.p<-ergm.ego(ego.obj_p ~edges + 
                    nodefactor("race",base=5) + 
                    nodematch("race",diff=TRUE) +
                    nodefactor("agecat", base=1) + 
                    absdiff("sqrt.age.adj") + 
                    nodefactor("deg.cohab.c",base=1) +
                    concurrent(by="sex") + 
                    concurrent(by="race") +
                    offset(nodematch("sex", diff=FALSE)),
                  offset.coef = -Inf,
                  constraints=~bd(maxout=3),
                  control=control.ergm.ego(ppopsize=50000, stats.est="asymptotic",
                                           ergm.control = control.ergm(MCMC.interval=7000,
                                                                       MCMC.samplesize=7000,
                                                                       MCMC.burnin = 7000,
                                                                       MPLE.max.dyad.types = 1e7,
                                                                       init.method = "zeros",
                                                                       MCMLE.maxit = 350,
                                                                       parallel = np, 
                                                                       parallel.type="PSOCK")))

plot(gof(fit.p, GOF="degree"))
plot(gof(fit.p, GOF="model"))


test<-simulate(fit.p)
test <- as.egodata(test)
degreedist.egodata(test, by="race.sex", prob=T)
summary(test ~ degree(0:5, by="race.sex"))

pdf(file = "XXX", height = 6, width = 12, pointsize = 16)
#par(mfrow = c(1,2), mar = c(4,3,2.5,1), mgp = c(2,.5,0))
degreedist.egodata(test, by="race.sex")
dev.off()


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
                                                                       MCMLE.maxit = 350,
                                                                       parallel = np, 
                                                                       parallel.type="PSOCK")))

summary(fit.i)



modelfits <- list(fit.c, fit.p, fit.i)
save(modelfits, file = "~/EpiModelHIV_shamp_modeling/scenarios/est/modelfits.rda")
