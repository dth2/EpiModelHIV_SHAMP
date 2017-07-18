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
##No testing 
param <- param_shamp(data.params,
                     tt.traj.B.f.prob = c(1, 0.000, 0, 0),
                     tt.traj.BI.f.prob = c(1, 0.000, 0, 0),
                     tt.traj.H.f.prob = c(1, 0.000, 0, 0),
                     tt.traj.HI.f.prob = c(1, 0.000, 0, 0),
                     tt.traj.W.f.prob = c(1, 0.000, 0, 0),
                     tt.traj.B.msf.prob = c(1, 0.000, 0, 0),
                     tt.traj.BI.msf.prob = c(1, 0.000, 0, 0),
                     tt.traj.H.msf.prob = c(1, 0.000, 0, 0),
                     tt.traj.HI.msf.prob = c(1, 0.000, 0, 0),
                     tt.traj.W.msf.prob = c(1, 0.000, 0, 0),
                     tt.traj.B.msm.prob = c(1, 0, 0, 0),
                     tt.traj.BI.msm.prob = c(1, 0, 0, 0),
                     tt.traj.H.msm.prob = c(1, 0, 0, 0),
                     tt.traj.HI.msm.prob = c(1, 0, 0, 0),
                     tt.traj.W.msm.prob = c(1, 0, 0, 0),
                     tt.traj.B.msmf.prob = c(1, 0.000, 0, 0),
                     tt.traj.BI.msmf.prob = c(1, 0.000, 0, 0),
                     tt.traj.H.msmf.prob = c(1, 0.000, 0, 0),
                     tt.traj.HI.msmf.prob = c(1, 0.000, 0, 0),
                     tt.traj.W.msmf.prob = c(1, 0.000, 0, 0))

init <- init_shamp()
control <- control_shamp(nsteps = 1040)

test_testing_nevertest<-netsim(est, param, init, control)

##Outside of package
save(test_testing_nevertest, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_testing_nevertest.rda")

##Run the Check

table(test_testing_nevertest$attr[[1]]$evertest)
table(test_testing_nevertest$attr[[1]]$last.neg.test)
table(test_testing_nevertest$attr[[1]]$diag.time)
table(test_testing_nevertest$attr[[1]]$diag.status)


