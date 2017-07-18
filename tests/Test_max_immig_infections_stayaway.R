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
                     immig.return.BI.f=0,
                     immig.return.HI.f=0,
                     immig.return.BI.m=0,
                     immig.return.HI.m=0,
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

test_max_immig_infections_stayaway<-netsim(est, param, init, control)

##Outside of package
save(test_max_immig_infections_stayaway, file = "~/EpiModelHIV_shamp_modeling/scenarios/test_max_immig_infections_stayaway.rda")

##Run the Check


sum(test_max_immig_infections_stayaway$epi$incid.L,na.rm = TRUE)
sum(test_max_immig_infections_stayaway$epi$incid.H,na.rm = TRUE)
sum(test_max_immig_infections_stayaway$epi$incid.FA,na.rm = TRUE)

table(test_max_immig_infections_stayaway$attr[[1]]$inf.class)

test_max_immig_infections_stayaway$epi$i.prev.W[,1][1040]
test_max_immig_infections_stayaway$epi$i.prev.B[,1][1040]
test_max_immig_infections_stayaway$epi$i.prev.BI[,1][1040]
test_max_immig_infections_stayaway$epi$i.prev.H[,1][1040]
test_max_immig_infections_stayaway$epi$i.prev.HI[,1][1040]