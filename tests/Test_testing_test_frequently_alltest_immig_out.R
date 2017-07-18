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
##Testing frequently no window immig leave
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
                     last.neg.test.B.f.int = 200,
                     last.neg.test.BI.f.int = 200,
                     last.neg.test.H.f.int = 200,
                     last.neg.test.HI.f.int = 200,
                     last.neg.test.W.f.int = 200,
                     last.neg.test.B.msf.int = 200,
                     last.neg.test.BI.msf.int = 200,
                     last.neg.test.H.msf.int = 200,
                     last.neg.test.HI.msf.int = 200,
                     last.neg.test.W.msf.int = 200,
                     last.neg.test.B.msm.int = 200,
                     last.neg.test.BI.msm.int = 200,
                     last.neg.test.H.msm.int = 200,
                     last.neg.test.HI.msm.int = 200,
                     last.neg.test.W.msm.int = 200,
                     last.neg.test.B.msmf.int = 200,
                     last.neg.test.BI.msmf.int = 200,
                     last.neg.test.H.msmf.int = 200,
                     last.neg.test.HI.msmf.int = 200,
                     last.neg.test.W.msmf.int = 200,
                     mean.test.B.f.int = 200,
                     mean.test.BI.f.int = 200,
                     mean.test.H.f.int = 200,
                     mean.test.HI.f.int = 200,
                     mean.test.W.f.int = 200,
                     mean.test.B.msf.int = 200,
                     mean.test.BI.msf.int = 200,
                     mean.test.H.msf.int = 200,
                     mean.test.HI.msf.int = 200,
                     mean.test.W.msf.int = 200,
                     mean.test.B.msm.int = 200,
                     mean.test.BI.msm.int = 200,
                     mean.test.H.msm.int = 200,
                     mean.test.HI.msm.int = 200,
                     mean.test.W.msm.int = 200,
                     mean.test.B.msmf.int = 200,
                     mean.test.BI.msmf.int = 200,
                     mean.test.H.msmf.int = 200,
                     mean.test.HI.msmf.int = 200,
                     mean.test.W.msmf.int = 200,
                     testing.pattern = "memoryless",
                     test.window.int = 1,
                     URVI.prob = .05,
                     UIVI.prob = .05,
                     immig.depart.BI.f=1,
                     immig.depart.HI.f=1,
                     immig.depart.BI.m=1,
                     immig.depart.HI.m=1,
                     immig.return.BI.f=.000001,
                     immig.return.HI.f=.000001,
                     immig.return.BI.m=.000001,
                     immig.return.HI.m=.000001)

init <- init_shamp()
control <- control_shamp(nsteps = 2040)

test_testing_frequently_alltest_immig_out<-netsim(est, param, init, control)


##Outside of package
save(test_testing_frequently_alltest_immig_out, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_testing_frequently_alltest_immig_out.rda")

##Run the Check

table(test_testing_frequently_alltest_immig_out$attr[[1]]$evertest,test_testing_frequently_alltest_immig_out$attr[[1]]$race)
table(test_testing_frequently_alltest_immig_out$attr[[1]]$diag.status,test_testing_frequently_alltest_immig_out$attr[[1]]$race)
table(test_testing_frequently_alltest_immig_out$attr[[1]]$status,test_testing_frequently_alltest_immig_out$attr[[1]]$race)
table(test_testing_frequently_alltest_immig_out$attr[[1]]$immig.loc,test_testing_frequently_alltest_immig_out$attr[[1]]$race)
