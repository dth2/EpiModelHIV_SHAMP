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


                     base.vi.main.W.rate = 0,
                     base.vi.pers.W.rate = 0,
                     URVI.prob = .02,
                     UIVI.prob = .02)

init <- init_shamp()
control <- control_shamp(nsteps = 1040)

test_acts_min_W_acts<-netsim(est, param, init, control)

##Outside of package
save(test_acts_min_W_acts, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_acts_min_W_acts.rda")

##Run the Check



table(test_acts_min_W_acts$attr[[1]]$inf.class)
table(test_acts_min_W_acts$attr[[1]]$status,test_acts_min_W_acts$attr[[1]]$race)

test_acts_min_W_acts$epi$i.prev[,1][1040]
