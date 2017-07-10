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
param <- param_shamp(data.params)

init <- init_shamp()
control <- control_shamp(nsteps = 1040)

test_baseline<-netsim(est, param, init, control)

##Outside of package
save(test_baseline, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_baseline.rda")

##Run the Check
time<-1:1040
test_baseline_plot<-plot(time,test_baseline$epi$i.prev[,1])

lines(test_baseline$epi$i.prev.B[,1],col="orange")
lines(test_baseline$epi$i.prev.BI[,1],col="yellow")
lines(test_baseline$epi$i.prev.H[,1],col="blue")
lines(test_baseline$epi$i.prev.HI[,1],col="lightblue")
lines(test_baseline$epi$i.prev.W[,1],col="green")
lines(test_baseline$epi$i.prev.m[,1],col="red")
lines(test_baseline$epi$i.prev.f[,1],col="pink")
lines(test_baseline$epi$i.prev.msmf[,1],col="yellow")
save(test_baseline_plot, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_baseline_plot.rda")


#If you want to output to an excel
library(xlsx) #load the package
write.xlsx(x = demog.table, file = "~/EpiModelHIV_shamp_modeling/scenarios/demogs.xlsx",
           sheetName = "SHAMP Demog", row.names = FALSE)


