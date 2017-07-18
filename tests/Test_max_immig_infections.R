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
                     immig.depart.BI.f=1/(52*4),
                     immig.depart.HI.f=1/(52*4),
                     immig.depart.BI.m=1/(52*4),
                     immig.depart.HI.m=1/(52*4),
                     immig.return.BI.f=1/4,
                     immig.return.HI.f=1/4,
                     immig.return.BI.m=1/4,
                     immig.return.HI.m=1/4,
                     immig.aq.prob.BI.f=1,
                     immig.aq.prob.HI.f=1,
                     immig.aq.prob.BI.m=1,
                     immig.aq.prob.HI.m=1,  
                     
                     msm.aq.prob.B=0,
                     msm.aq.prob.BI=0,
                     msm.aq.prob.H=0,
                     msm.aq.prob.HI=0,
                     msm.aq.prob.W=0)

init <- init_shamp()
control <- control_shamp(nsteps = 1040)

test_max_immig_infections<-netsim(est, param, init, control)

##Outside of package
save(test_max_immig_infections, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_max_immig_infections.rda")

##Run the Check


sum(test_max_immig_infections$epi$incid.L,na.rm = TRUE)
sum(test_max_immig_infections$epi$incid.H,na.rm = TRUE)
sum(test_max_immig_infections$epi$incid.FA,na.rm = TRUE)

table(test_max_immig_infections$attr[[1]]$inf.class)

