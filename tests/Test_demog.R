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

param <- param_shamp(data.params)
init <- init_shamp()
control <- control_shamp(nsteps = 1560)


test_demog<-netsim(est, param, init, control)


##Outside of package
save(test_demog, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_demog.rda")

##Run the Check


time<-1:1560
test_demog_plot<-plot(time,test_demog$epi$num.W.f[,1],col="black", ylim=c(0, 15000))

lines(test_demog$epi$num.BI.f[,1],col="orange")
lines(test_demog$epi$num.H.f[,1],col="yellow")
lines(test_demog$epi$num.HI.f[,1],col="blue")
lines(test_demog$epi$num.B.f[,1],col="lightblue")

lines(test_demog$epi$num.B.m[,1],col="green")
lines(test_demog$epi$num.BI.m[,1],col="red")
lines(test_demog$epi$num.H.m[,1],col="pink")
lines(test_demog$epi$num.HI.m[,1],col="grey")
lines(test_demog$epi$num.W.m[,1],col="purple")

save(test_demog_plot, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_demog_plot.rda")

test_demog_table<-as.data.frame(rbind(demog4,demog52,demog104,demog208,demog312,demog416,demog520,demog1040))

save(test_demog_table, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_demog_table.rda")


#library(xlsx) #load the package
#write.xlsx(x = test_demog_table, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_demog_table.xlsx",
#           sheetName = "SHAMP Demog", row.names = FALSE)

