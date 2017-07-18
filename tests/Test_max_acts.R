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
##Max acts 
param <- param_shamp(data.params,


                     base.vi.main.B.rate = 100,
                     base.vi.main.BI.rate = 100,
                     base.vi.main.H.rate = 100,
                     base.vi.main.HI.rate = 100,
                     base.vi.main.W.rate = 100,
                     base.vi.pers.B.rate = 100,
                     base.vi.pers.BI.rate = 100,
                     base.vi.pers.H.rate = 100,
                     base.vi.pers.HI.rate = 100,
                     base.vi.pers.W.rate = 100,
                     
                     msm.aq.prob.B=0,
                     msm.aq.prob.BI=0,
                     msm.aq.prob.H=0,
                     msm.aq.prob.HI=0,
                     msm.aq.prob.W=0,
                     immig.aq.prob.BI.f=0,
                     immig.aq.prob.HI.f=0,
                     immig.aq.prob.BI.m=0,
                     immig.aq.prob.HI.m=0)

init <- init_shamp()
control <- control_shamp(nsteps = 1040)

test_max_acts<-netsim(est, param, init, control)

##Outside of package
save(test_max_acts, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_max_acts.rda")

##Run the Check



table(test_max_acts$attr[[1]]$inf.class)
table(test_max_acts$attr[[1]]$status,test_max_acts$attr[[1]]$race)

test_max_acts$epi$i.prev[,1][1040]
