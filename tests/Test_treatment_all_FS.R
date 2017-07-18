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
##Treatment emedidiate on pos test.  All treaters.
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
                     tx.init.B.f.prob = 1,
                     tx.init.BI.f.prob = 1,
                     tx.init.H.f.prob = 1,
                     tx.init.HI.f.prob = 1,
                     tx.init.W.f.prob = 1,
                     tx.init.B.msf.prob = 1,
                     tx.init.BI.msf.prob = 1,
                     tx.init.H.msf.prob = 1,
                     tx.init.HI.msf.prob = 1,
                     tx.init.W.msf.prob = 1,
                     tx.init.B.msm.prob = 1,
                     tx.init.BI.msm.prob = 1,
                     tx.init.H.msm.prob = 1,
                     tx.init.HI.msm.prob = 1,
                     tx.init.W.msm.prob = 1,
                     tx.init.B.msmf.prob = 1,
                     tx.init.BI.msmf.prob = 1,
                     tx.init.H.msmf.prob = 1,
                     tx.init.HI.msmf.prob = 1,
                     tx.init.W.msmf.prob = 1,
                     URVI.prob = .01,
                     UIVI.prob = .01)

init <- init_shamp()
control <- control_shamp(nsteps = 1040)

test_treatment_all_FS<-netsim(est, param, init, control)

##Outside of package
save(test_treatment_all_FS, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_treatment_all_FS.rda")

##Run the Check


table(test_treatment_all_FS$attr[[1]]$diag.status,test_treatment_all_FS$attr[[1]]$tx.status)
table(test_treatment_all_FS$attr[[1]]$status)

