library(EpiModel)
library(data.table)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(plyr)
library(pracma)
library(stats)
library(future)
library(doParallel)
rm(list=ls())

#You'll need to run the covid_sims_func.R file first before running the code below. 

# function to set-up and run the baseline simulations
simulate <- function(# control.icm params
  type = "SEIQHRF", 
  nsteps = 366, 
  nsims = 8,
  ncores = 4,
  prog.rand = FALSE,
  rec.rand = FALSE,
  fat.rand = TRUE,
  quar.rand = FALSE,
  hosp.rand = FALSE,
  disch.rand = TRUE,
  infection.FUN = infection.seiqhrf.icm,
  recovery.FUN = progress.seiqhrf.icm,
  departures.FUN = departures.seiqhrf.icm,
  arrivals.FUN = arrivals.icm,
  get_prev.FUN = get_prev.seiqhrf.icm,
  # init.icm params
  s.num = 9997,
  e.num=0,
  i.num = 3,
  q.num=0,
  h.num=0,
  r.num = 0,
  f.num = 0,
  # param.icm params
  inf.prob.e = 0.02, 
  act.rate.e = 10,
  inf.prob.i = 0.05, 
  act.rate.i = 10,
  inf.prob.q = 0.02, 
  act.rate.q = 2.5,                    
  quar.rate = 1/30, 
  hosp.rate = 1/100,
  disch.rate = 1/15,
  prog.rate = 1/10,
  prog.dist.scale = 5,
  prog.dist.shape = 1.5,
  rec.rate = 1/20,
  rec.dist.scale = 35,
  rec.dist.shape = 1.5,
  fat.rate.base = 1/50,
  hosp.cap = 40,
  fat.rate.overcap = 1/25,
  fat.tcoeff = 0.5,
  vital = TRUE,
  a.rate = (10.5/365)/1000, 
  a.prop.e = 0.01,
  a.prop.i = 0.001,
  a.prop.q = 0.01,
  ds.rate = (7/365)/1000, 
  de.rate = (7/365)/1000, 
  di.rate = (7/365)/1000,
  dq.rate = (7/365)/1000,
  dh.rate = (20/365)/1000,
  dr.rate = (7/365)/1000,
  out="mean"
) {
  
  control <- control.icm(type = type, 
                         nsteps = nsteps, 
                         nsims = nsims,
                         ncores = ncores,
                         prog.rand = prog.rand,
                         rec.rand = rec.rand,
                         infection.FUN = infection.FUN,
                         recovery.FUN = recovery.FUN,
                         arrivals.FUN = arrivals.FUN,
                         departures.FUN = departures.FUN,
                         get_prev.FUN = get_prev.FUN)
  
  init <- init.icm(s.num = s.num,
                   e.num = e.num,
                   i.num = i.num,
                   q.num = q.num,
                   h.num = h.num,
                   r.num = r.num,
                   f.num = f.num)
  
  param <-  param.icm(inf.prob.e = inf.prob.e, 
                      act.rate.e = act.rate.e,
                      inf.prob.i = inf.prob.i, 
                      act.rate.i = act.rate.i,
                      inf.prob.q = inf.prob.q, 
                      act.rate.q = act.rate.q,                    
                      quar.rate = quar.rate,
                      hosp.rate = hosp.rate,
                      disch.rate = disch.rate,
                      prog.rate = prog.rate,
                      prog.dist.scale = prog.dist.scale,
                      prog.dist.shape = prog.dist.shape,
                      rec.rate = rec.rate,
                      rec.dist.scale = rec.dist.scale,
                      rec.dist.shape = rec.dist.shape,
                      fat.rate.base = fat.rate.base,
                      hosp.cap = hosp.cap,
                      fat.rate.overcap = fat.rate.overcap,
                      fat.tcoeff = fat.tcoeff,
                      vital = vital,
                      a.rate = a.rate, 
                      a.prop.e = a.prop.e,
                      a.prop.i = a.prop.i,
                      a.prop.q = a.prop.q,
                      ds.rate = ds.rate, 
                      de.rate = de.rate, 
                      di.rate = di.rate,
                      dq.rate = dq.rate,
                      dh.rate = dh.rate,
                      dr.rate = dr.rate)
  
  sim <- icm.seiqhrf(param, init, control)
  sim_df <- as.data.frame(sim, out=out)
  
  return(list(sim=sim, df=sim_df))
}

tic()
baseline_sim <- simulate(ncores = 4, type = "SIS")
toc()

infectious_hygiene_delayed_ramp <- function(t) {
  ifelse(t < 30, 0.05, ifelse(t <= 60, 0.05 - (t - 30) * (0.05 - 
                                                            0.02)/30, 0.02))
}

infectious_hygiene_delayed_ramp_sim <- simulate(inf.prob.i = infectious_hygiene_delayed_ramp(1:366), type = "SIS")


#Simulations
#looping through 20 companies, random size between 10 and 1000
set.seed(1000)
infected_diff = list()
for (i in 1:100){
  inf_diff_sub <- tibble(numemp=NA, infect_difference_novac_vac=NA)
  numemp <- sample(10:1000, 1)
  infectious_noempvac <- simulate(inf.prob.i = 0.05, s.num = round(.3*numemp),  r.num =  round(.54*numemp), i.num = round(.06*numemp), act.rate.e= round(.1*numemp),  type = "SIR")
  infectious_empvac <- simulate(inf.prob.i = 0.05, s.num = .5*(round(.3*numemp)),  r.num = round(.54*numemp) + .5*(round(.3*numemp)), i.num = round(.06*numemp),act.rate.e= round(.1*numemp),  type = "SIR")

  baseline_vac <- infectious_empvac$df %>% 
  dplyr::select(time, s.num,  i.num,  num)
  baseline_novac <- infectious_noempvac$df %>% 
    dplyr::select(time, s.num,  i.num,  num)
  
  inf_diff_sub$infect_difference_novac_vac = mean(baseline_novac$i.num) -  mean(baseline_vac$i.num)
  inf_diff_sub$numemp = numemp
  infected_diff[[i]] = inf_diff_sub
}  

#Appending output:
inf_output = do.call("rbind", infected_diff)
colnames(inf_output) = c("num_employees", "inf_diff")

mean(inf_output$inf_diff)
#1.66277

#Small businesses have 100 or fewer employees.  They employ 48.8% of CA employees, and constitute 98.8% of CA businesses.
#1% are medium businesses, with 100-500 employees.
#.2% are large, with over 500 employees.


#Small firms
set.seed(1000)
infected_diff_small = list()
for (i in 1:80){
  inf_diff_sub <- tibble(numemp=NA, infect_difference_novac_vac=NA)
  numemp <- sample(2:100, 1)
  infectious_noempvac <- simulate(inf.prob.i = 0.05, s.num = round(.3*numemp),  r.num =  round(.54*numemp), i.num = round(.06*numemp), act.rate.e= round(.1*numemp),  type = "SIR")
  infectious_empvac <- simulate(inf.prob.i = 0.05, s.num = .5*(round(.3*numemp)),  r.num = round(.54*numemp) + .5*(round(.3*numemp)), i.num = round(.06*numemp),act.rate.e= round(.1*numemp),  type = "SIR")
  
  baseline_vac <- infectious_empvac$df %>% 
    dplyr::select(time, s.num,  i.num,  num)
  baseline_novac <- infectious_noempvac$df %>% 
    dplyr::select(time, s.num,  i.num,  num)
  
  inf_diff_sub$infect_difference_novac_vac = mean(baseline_novac$i.num) -  mean(baseline_vac$i.num)
  inf_diff_sub$numemp = numemp
  infected_diff_small[[i]] = inf_diff_sub
}  

#Appending output:
inf_output_small = do.call("rbind", infected_diff_small)
colnames(inf_output_small) = c("num_employees", "inf_diff")
mean(inf_output_small$inf_diff)
#0.1231088


#Medium firms
set.seed(1000)
infected_diff_medium = list()
for (i in 1:80){
  inf_diff_sub <- tibble(numemp=NA, infect_difference_novac_vac=NA)
  numemp <- sample(101:500, 1)
  infectious_noempvac <- simulate(inf.prob.i = 0.05, s.num = round(.3*numemp),  r.num =  round(.54*numemp), i.num = round(.06*numemp), act.rate.e= round(.1*numemp),  type = "SIR")
  infectious_empvac <- simulate(inf.prob.i = 0.05, s.num = .5*(round(.3*numemp)),  r.num = round(.54*numemp) + .5*(round(.3*numemp)), i.num = round(.06*numemp),act.rate.e= round(.1*numemp),  type = "SIR")
  
  baseline_vac <- infectious_empvac$df %>% 
    dplyr::select(time, s.num,  i.num,  num)
  baseline_novac <- infectious_noempvac$df %>% 
    dplyr::select(time, s.num,  i.num,  num)
  
  inf_diff_sub$infect_difference_novac_vac = mean(baseline_novac$i.num) -  mean(baseline_vac$i.num)
  inf_diff_sub$numemp = numemp
  infected_diff_medium[[i]] = inf_diff_sub
}  

#Appending output:
inf_output_medium = do.call("rbind", infected_diff_medium)
colnames(inf_output_medium) = c("num_employees", "inf_diff")
mean(inf_output_medium$inf_diff)
# 0.9291837



#Large firms
set.seed(1000)
infected_diff_large = list()
for (i in 1:80){
  inf_diff_sub <- tibble(numemp=NA, infect_difference_novac_vac=NA)
  numemp <- sample(501:5000, 1)
  infectious_noempvac <- simulate(inf.prob.i = 0.05, s.num = round(.3*numemp),  r.num =  round(.54*numemp), i.num = round(.06*numemp), act.rate.e= round(.1*numemp),  type = "SIR")
  infectious_empvac <- simulate(inf.prob.i = 0.05, s.num = .5*(round(.3*numemp)),  r.num = round(.54*numemp) + .5*(round(.3*numemp)), i.num = round(.06*numemp),act.rate.e= round(.1*numemp),  type = "SIR")
  
  baseline_vac <- infectious_empvac$df %>% 
    dplyr::select(time, s.num,  i.num,  num)
  baseline_novac <- infectious_noempvac$df %>% 
    dplyr::select(time, s.num,  i.num,  num)
  
  inf_diff_sub$infect_difference_novac_vac = mean(baseline_novac$i.num) -  mean(baseline_vac$i.num)
  inf_diff_sub$numemp = numemp
  infected_diff_large[[i]] = inf_diff_sub
}  

#Appending output:
inf_output_large = do.call("rbind", infected_diff_large)
colnames(inf_output_large) = c("num_employees", "inf_diff")
mean(inf_output_large$inf_diff)
#8.510959

#Function
infection_diff = function(minemp, maxemp, nsims){
  infected_diff = list()
  for (i in 1:nsims){
    inf_diff_sub <- tibble(numemp=NA, infect_difference_novac_vac=NA)
    numemp <- sample(minemp:maxemp, 1)
    infectious_noempvac <- simulate(inf.prob.i = 0.05, s.num = round(.3*numemp),  r.num =  round(.54*numemp), i.num = round(.06*numemp), act.rate.e= round(.1*numemp),  type = "SIR")
    infectious_empvac <- simulate(inf.prob.i = 0.05, s.num = .5*(round(.3*numemp)),  r.num = round(.54*numemp) + .5*(round(.3*numemp)), i.num = round(.06*numemp),act.rate.e= round(.1*numemp),  type = "SIR")
    
    baseline_vac <- infectious_empvac$df %>% 
      dplyr::select(time, s.num,  i.num,  num)
    baseline_novac <- infectious_noempvac$df %>% 
      dplyr::select(time, s.num,  i.num,  num)
    
    inf_diff_sub$infect_difference_novac_vac = mean(baseline_novac$i.num) -  mean(baseline_vac$i.num)
    inf_diff_sub$numemp = numemp
    infected_diff[[i]] = inf_diff_sub
  }
  #Appending output:
  inf_output = do.call("rbind", infected_diff)
  colnames(inf_output) = c("num_employees", "inf_diff")
  mean(inf_output_large$inf_diff)
  return(inf_output)
}

infection_diff_small = infection_diff(2, 100, 100)
infection_diff_medium = infection_diff(101, 500, 100)
infection_diff_large = infection_diff(501, 5000, 100)

all_sims = rbind(infection_diff_small, infection_diff_medium, infection_diff_large)
mean(all_sims$inf_diff)
#3.367658

sum(all_sims$num_employees)
# 319422



infection_diff_carate = function(minemp, maxemp, nsims){
  infected_diff = list()
  for (i in 1:nsims){
    inf_diff_sub <- tibble(numemp=NA, infect_difference_novac_vac=NA)
    numemp <- sample(minemp:maxemp, 1)
    infectious_noempvac <- simulate(inf.prob.i = 0.05, s.num = round(.2*numemp),  r.num =  round(.68*numemp), i.num = round(.04*numemp), act.rate.e= round(.08*numemp),  type = "SIR")
    infectious_empvac <- simulate(inf.prob.i = 0.05, s.num = .5*(round(.2*numemp)),  r.num = round(.68*numemp) + .5*(round(.2*numemp)), i.num = round(.04*numemp),act.rate.e= round(.08*numemp),  type = "SIR")
    
    baseline_vac <- infectious_empvac$df %>% 
      dplyr::select(time, s.num,  i.num,  num)
    baseline_novac <- infectious_noempvac$df %>% 
      dplyr::select(time, s.num,  i.num,  num)
    
    inf_diff_sub$infect_difference_novac_vac = mean(baseline_novac$i.num) -  mean(baseline_vac$i.num)
    inf_diff_sub$numemp = numemp
    infected_diff[[i]] = inf_diff_sub
  }
  #Appending output:
  inf_output = do.call("rbind", infected_diff)
  colnames(inf_output) = c("num_employees", "inf_diff")
  mean(inf_output_large$inf_diff)
  return(inf_output)
}

infection_diff_small_ca = infection_diff_carate(2, 100, 100)
infection_diff_medium_ca = infection_diff_carate(101, 500, 100)
infection_diff_large_ca = infection_diff_carate(501, 5000, 100)

mean(infection_diff_small_ca$inf_diff)
#0.1120458
mean(infection_diff_medium_ca$inf_diff)
# 0.6285451
mean(infection_diff_large_ca$inf_diff)
#5.677565


infection_diff_lowrate = function(minemp, maxemp, nsims){
  infected_diff = list()
  for (i in 1:nsims){
    inf_diff_sub <- tibble(numemp=NA, infect_difference_novac_vac=NA)
    numemp <- sample(minemp:maxemp, 1)
    infectious_noempvac <- simulate(inf.prob.i = 0.05, s.num = round(.54*numemp),  r.num =  round(.30*numemp), i.num = round(.06*numemp), act.rate.e= round(.1*numemp),  type = "SIR")
    infectious_empvac <- simulate(inf.prob.i = 0.05, s.num = .5*(round(.54*numemp)),  r.num = round(.30*numemp) + .5*(round(.54*numemp)), i.num = round(.06*numemp),act.rate.e= round(.1*numemp),  type = "SIR")
    
    baseline_vac <- infectious_empvac$df %>% 
      dplyr::select(time, s.num,  i.num,  num)
    baseline_novac <- infectious_noempvac$df %>% 
      dplyr::select(time, s.num,  i.num,  num)
    
    inf_diff_sub$infect_difference_novac_vac = mean(baseline_novac$i.num) -  mean(baseline_vac$i.num)
    inf_diff_sub$numemp = numemp
    infected_diff[[i]] = inf_diff_sub
  }
  #Appending output:
  inf_output = do.call("rbind", infected_diff)
  colnames(inf_output) = c("num_employees", "inf_diff")
  mean(inf_output_large$inf_diff)
  return(inf_output)
}

infection_diff_small_low = infection_diff_lowrate(2, 100, 100)
infection_diff_medium_low = infection_diff_lowrate(101, 500, 100)
infection_diff_large_low = infection_diff_lowrate(501, 5000, 100)

mean(infection_diff_small_low$inf_diff)
# 0.2636749
mean(infection_diff_medium_low$inf_diff)
# 1.614126
mean(infection_diff_large_low$inf_diff)
#15.8606
