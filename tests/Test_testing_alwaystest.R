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
##Testing dayly no window 
param <- param_shamp(data.params,
                     tt.traj.B.f.prob = c(0, 0.000, .5, .5),
                     tt.traj.BI.f.prob = c(0, 0.000, .5, .5),
                     tt.traj.H.f.prob = c(0, 0.000, .5, .5),
                     tt.traj.HI.f.prob = c(0, 0.000, .5, .5),
                     tt.traj.W.f.prob = c(0, 0.000, .5, .5),
                     tt.traj.B.msf.prob = c(0, 0.000, .5, .5),
                     tt.traj.BI.msf.prob = c(0, 0.000, .5, .5),
                     tt.traj.H.msf.prob = c(0, 0.000, .5, .5),
                     tt.traj.HI.msf.prob = c(0, 0.000, .5, .5),
                     tt.traj.W.msf.prob = c(0, 0.000, .5, .5),
                     tt.traj.B.msm.prob = c(0, 0.000, .5, .5),
                     tt.traj.BI.msm.prob = c(0, 0.000, .5, .5),
                     tt.traj.H.msm.prob = c(0, 0.000, .5, .5),
                     tt.traj.HI.msm.prob = c(0, 0.000, .5, .5),
                     tt.traj.W.msm.prob = c(0, 0.000, .5, .5),
                     tt.traj.B.msmf.prob = c(0, 0.000, .5, .5),
                     tt.traj.BI.msmf.prob = c(0, 0.000, .5, .5),
                     tt.traj.H.msmf.prob = c(0, 0.000, .5, .5),
                     tt.traj.HI.msmf.prob = c(0, 0.000, .5, .5),
                     tt.traj.W.msmf.prob = c(0, 0.000, .5, .5),
                     last.neg.test.B.f.int = 2,
                     last.neg.test.BI.f.int = 2,
                     last.neg.test.H.f.int = 2,
                     last.neg.test.HI.f.int = 2,
                     last.neg.test.W.f.int = 2,
                     last.neg.test.B.msf.int = 2,
                     last.neg.test.BI.msf.int = 2,
                     last.neg.test.H.msf.int = 2,
                     last.neg.test.HI.msf.int = 2,
                     last.neg.test.W.msf.int = 2,
                     last.neg.test.B.msm.int = 2,
                     last.neg.test.BI.msm.int = 2,
                     last.neg.test.H.msm.int = 2,
                     last.neg.test.HI.msm.int = 2,
                     last.neg.test.W.msm.int = 2,
                     last.neg.test.B.msmf.int = 2,
                     last.neg.test.BI.msmf.int = 2,
                     last.neg.test.H.msmf.int = 2,
                     last.neg.test.HI.msmf.int = 2,
                     last.neg.test.W.msmf.int = 2,
                     mean.test.B.f.int = 2,
                     mean.test.BI.f.int = 2,
                     mean.test.H.f.int = 2,
                     mean.test.HI.f.int = 2,
                     mean.test.W.f.int = 2,
                     mean.test.B.msf.int = 2,
                     mean.test.BI.msf.int = 2,
                     mean.test.H.msf.int = 2,
                     mean.test.HI.msf.int = 2,
                     mean.test.W.msf.int = 2,
                     mean.test.B.msm.int = 2,
                     mean.test.BI.msm.int = 2,
                     mean.test.H.msm.int = 2,
                     mean.test.HI.msm.int = 2,
                     mean.test.W.msm.int = 2,
                     mean.test.B.msmf.int = 2,
                     mean.test.BI.msmf.int = 2,
                     mean.test.H.msmf.int = 2,
                     mean.test.HI.msmf.int = 2,
                     mean.test.W.msmf.int = 2,
                     testing.pattern = "memoryless",
                     test.window.int = 1,
                     URVI.prob = .01,
                     UIVI.prob = .01)

init <- init_shamp()
control <- control_shamp(nsteps = 1040)

test_testing_daily<-netsim(est, param, init, control)

##Outside of package
save(test_testing_daily, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_testing_daily.rda")

##Run the Check

table(test_testing_nevertest$attr[[1]]$evertest)
table(test_testing_nevertest$attr[[1]]$last.neg.test)
table(test_testing_nevertest$attr[[1]]$diag.time)
table(test_testing_nevertest$attr[[1]]$diag.status)


