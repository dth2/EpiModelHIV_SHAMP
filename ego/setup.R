##Running network simulations with EPIModelHIV and ergm.ego.

#Load required packages.
suppressPackageStartupMessages(library(EpiModelHIV))
suppressPackageStartupMessages(library(EpiModelHPC))
suppressPackageStartupMessages(library(ergm.ego))

##Read in ego and alter data for ergm.ego.
data.egos<-read.csv("ego_temp.csv",sep=",")
data.egos$race<-as.character(data.egos$race)
data.egos$sex<-as.character(data.egos$sex)
data.alters<-read.csv("alter_temp.csv",sep=",")
data.alters$race<-as.character(data.alters$race)
data.alters$sex<-as.character(data.alters$sex)

##Make the ergm.ego object.
ego.obj<-as.egodata(data.egos,alters=data.alters,egoIDcol="egoid")

##Create the network model using ergm.ego.

model<-ergm.ego(ego.obj ~edges + degree(1)+nodefactor("sex",base=1)+nodemix("race"), ppopsize=10000)


#Create additional required elements for est
# Time unit for simulation, relative to 1 day
time.unit <- 7

# Mean durations
rates <- .049505
durs  <- 1/rates
diss = ~offset(edges)

#Mortality
ages <- 18:39
age.unit <- 52
asmr.B <- c(rep(0, 17),
            1-(1-c(rep(0.00159, 7),
                   rep(0.00225, 10),
                   rep(0.00348, 5)))^(1 / age.unit),
            1)
asmr.BI <- c(rep(0, 17),
             1-(1-c(rep(0.00159, 7),
                    rep(0.00225, 10),
                    rep(0.00348, 5)))^(1 / age.unit),
             1)

asmr.H <- c(rep(0, 17),
            1-(1-c(rep(0.00103, 7),
                   rep(0.00133, 10),
                   rep(0.00214, 5)))^(1/age.unit),
            1)
asmr.HI <- c(rep(0, 17),
             1-(1-c(rep(0.00103, 7),
                    rep(0.00133, 10),
                    rep(0.00214, 5)))^(1/age.unit),
             1)
asmr.W <- c(rep(0, 17),
            1-(1-c(rep(0.00103, 7),
                   rep(0.00133, 10),
                   rep(0.00214, 5)))^(1/age.unit),
            1)

exp.mort <- (mean(asmr.B[ages]) + mean(asmr.BI[ages]) + mean(asmr.H[ages])
             + mean(asmr.HI[ages]) + mean(asmr.W[ages])) / 5

coef.diss <- dissolution_coefs(dissolution = diss,
                               duration = durs / time.unit,
                               d.rate = exp.mort)


target.stats<-as.numeric(model$target.stats[2:length(model$target.stats)])
                         
target.stats.names<- names(model$target.stats[2:length(model$target.stats)])
coef.form.crude<-model$coef[2:length(model$coef)]
coef.form<-coef.form.crude
coef.form[1]<- coef.form[1]- coef.diss$coef.adj
constraints <- ~.
                                                    
est <- list(fit= model, formation=model$formula, target.stats= target.stats, 
            target.stats.names= target.stats.names, coef.form = coef.form, 
            coef.form.crude= coef.form.crude, coef.diss, constraints= constraints, 
            edapprox=TRUE)
                                                    
param <- param_het()
init <- init_het(i.prev.male = 0.25, i.prev.feml = 0.25)
control <- control_het(nsteps = 2)

sim <- netsim(est, param, init, control)
