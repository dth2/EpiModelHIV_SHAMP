##Running network simulations with EPIModelHIV and ergm.ego.

#Load required packages.
library(EpiModelHIV)
library(EpiModelHPC)
#devtools::install_github("statnet/ergm.ego-private",ref="3.7-compat", auth_token ="")
#devtools::install_github("statnet/tergmLite")
#devtools::install_github("statnet/ergm.ego")
library(ergm.ego)
library(tergmLite)



##Load data
load(file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fitsmall.rda")
load(file = "~/EpiModelHIV_shamp_modeling/scenarios/est/data.params.rda")
source(file= "~/EpiModelHIV_SHAMP2/tests/age_mort.R")

############################################################################
##No infection 
param <- param_shamp(data.params,
                     URAI.prob = 0,
                     UIAI.prob = 0,
                     URVI.prob = 0,
                     UIVI.prob = 0)

init <- init_shamp()
control <- control_shamp(nsteps = 1040)

test_death_noinfection<-netsim(est, param, init, control)

##Outside of package
save(test_death_noinfection, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_death_noinfection.rda")

##Run the Check
sum(test_death_noinfection$epi$dth.gen,na.rm=TRUE)
sum(test_death_noinfection$epi$dth.dis,na.rm=TRUE)
sum(test_death_noinfection$epi$dth.age,na.rm=TRUE)


######################################################################################
##Max infection 
param <- param_shamp(data.params,
                     URAI.prob = .05,
                     UIAI.prob = .05,
                     URVI.prob = .05,
                     UIVI.prob = .05)

init <- init_shamp()
control <- control_shamp(nsteps = 1040)

test_death_maxinfection<-netsim(est, param, init, control)

##Outside of package
save(test_death_maxinfection, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_death_maxinfection.rda")

##Run the Check
sum(test_death_maxinfection$epi$dth.gen,na.rm=TRUE)
sum(test_death_maxinfection$epi$dth.dis,na.rm=TRUE)
sum(test_death_maxinfection$epi$dth.age,na.rm=TRUE)

time<-1:1040
test_death_maxinfection_plot<-plot(time,test_death_maxinfection$epi$i.prev[,1])

lines(test_death_maxinfection$epi$i.prev.B[,1],col="orange")
lines(test_death_maxinfection$epi$i.prev.BI[,1],col="yellow")
lines(test_death_maxinfection$epi$i.prev.H[,1],col="blue")
lines(test_death_maxinfection$epi$i.prev.HI[,1],col="lightblue")
lines(test_death_maxinfection$epi$i.prev.W[,1],col="green")
lines(test_death_maxinfection$epi$i.prev.m[,1],col="red")
lines(test_death_maxinfection$epi$i.prev.f[,1],col="pink")
lines(test_death_maxinfection$epi$i.prev.msmf[,1],col="yellow")
save(test_death_maxinfection, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_death_maxinfection_plot.rda")


######################################################################################
##No background mortality infection 


asmr.B.f <- asmr.BI.f <- asmr.H.f <- asmr.HI.f <- asmr.W.f <- c(rep(0,59),1) 
asmr.B.m <- asmr.BI.m <- asmr.H.m <- asmr.HI.m <- asmr.W.m <- c(rep(0,59),1)
  
param <- param_shamp(data.params)
init <- init_shamp()
control <- control_shamp(nsteps = 1040)

test_death_nobackmort<-netsim(est, param, init, control)

##Outside of package
save(test_death_nobackmort, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_death_nobackmort")

##Run the Check
sum(test_death_nobackmort$epi$dth.gen,na.rm=TRUE)
sum(test_death_nobackmort$epi$dth.dis,na.rm=TRUE)
sum(test_death_nobackmort$epi$dth.age,na.rm=TRUE)



##############################################################################################

#If you want to output to an excel
#library(xlsx) #load the package
#write.xlsx(x = demog.table, file = "~/EpiModelHIV_shamp_modeling/scenarios/demogs.xlsx",
#           sheetName = "SHAMP Demog", row.names = FALSE)


