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
##Birth counts Vs Death counts
param <- param_shamp(data.params)

init <- init_shamp()
control <- control_shamp(nsteps = 1040)

test_births<-netsim(est, param, init, control)

##Outside of package
save(test_births, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_births.rda")

##Run the Check
sum(test_births$epi$dth.gen,na.rm=TRUE)
sum(test_births$epi$dth.dis,na.rm=TRUE)
sum(test_births$epi$dth.age,na.rm=TRUE)
sum(test_births$epi$nBirths,na.rm=TRUE)

(sum(test_births$epi$dth.gen,na.rm=TRUE)+sum(test_births$epi$dth.age,na.rm=TRUE))-(sum(test_births$epi$nBirths,na.rm=TRUE))

##############################################################################################

#If you want to output to an excel
#library(xlsx) #load the package
#write.xlsx(x = demog.table, file = "~/EpiModelHIV_shamp_modeling/scenarios/demogs.xlsx",
#           sheetName = "SHAMP Demog", row.names = FALSE)


