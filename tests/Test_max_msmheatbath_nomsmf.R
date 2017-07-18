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

                     msm.aq.prob.B=.1,
                     msm.aq.prob.BI=.1,
                     msm.aq.prob.H=.1,
                     msm.aq.prob.HI=.1,
                     msm.aq.prob.W=.1,
                     msmf.frac=0)

init <- init_shamp()
control <- control_shamp(nsteps = 1040)

test_max_msmheatbath_nomsmf<-netsim(est, param, init, control)

##Outside of package
save(test_max_msmheatbath_nomsmf, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_max_msmheatbath_nomsmf.rda")

##Run the Check



table(test_max_msmheatbath_nomsmf$attr[[1]]$inf.class)
table(test_max_msmheatbath_nomsmf$attr[[1]]$status,test_max_msmheatbath_nomsmf$attr[[1]]$race)
table(test_max_msmheatbath_nomsmf$attr[[1]]$status,test_max_msmheatbath_nomsmf$attr[[1]]$sex.ident)
table(floor(test_max_msmheatbath_nomsmf$attr[[1]]$age),test_max_msmheatbath_nomsmf$attr[[1]]$sex.ident)
test_max_msmheatbath_nomsmf$epi$num.msmf[,1][1040]
test_max_msmheatbath_nomsmf$epi$i.prev.msmf[,1][1040]
