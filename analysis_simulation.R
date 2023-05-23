############## DATA ANALYSIS #################################
library(tidyverse)
library(dtplyr)
library(reshape2)
library(foreach)
library(doParallel)
library(tictoc)

setwd(...) #Set your location of all the files here

# load functions for DPA 
source('helping_functions.R')

#load data
load('mydata.reg.rda')

# Real effects
load('realeffects.rda')
load("real_aa1_aastar0.rda")
load("real_aa0_aastar0.rda")
load("real_aa1_aastar1.rda")
colorval= "#b433de" #for plotting purposes later


################################################################
## What variables are we looking at?

summary(mydata.reg)

## To use our functions the raw data must be in given format:
# patient - numeric unique with values 1,...,n
# A - binary numeric with values 0 and 1
# Z - numeric
# E - binary numeric with values 0 and 1 (1 indicating event)
# study_exit - numeric, time to death time (max value 13)
# M - binary numeric with values 0 and 1 (1 indicating mediator event)
# M_time - numeric, time to mediator event time (max value 13) -- if the individual died (or got censored) before getting the mediator, this is then the death (or censoring) time


#################################################################
timebyy <- 0.5
mydata <- mydata.reg

# Fitting the model for the mediator
M_fits <- mydata %>% 
  pool_mediator_data(timeby=timebyy) %>% 
  pooled_logreg() %>% 
  pool_event_data(data=mydata,
                  timeby=timebyy)


# Fitting the model for the survival outcome
additive_fits <- my.additive.new(regformula = "~ A + M_rs + Z",
                                 event="E",
                                 startt="int_begin",
                                 stopt="int_end",
                                 dataset=M_fits$pooled_event_data)
additive_fits$coeff[,"M_rs"] <- ifelse(is.na(additive_fits$coeff[,"M_rs"]), 0, additive_fits$coeff[,"M_rs"]) #turn to 0 where indir effect is undefined


# Manipulate the results a bit and get them ready for estimating the survival functions
data_long <- M_fits$pooled_event_data
logreg_fits <- M_fits$logreg_coefs
obstimes.imput <- sort(unique(data_long$int_end[data_long$E==1] )) #get event times

additive_coefs <- data.frame(time = obstimes.imput,
                             mu_s = additive_fits$coeff[,"(Intercept)"],
                             alpha_s = additive_fits$coeff[,"A"],
                             beta_s = additive_fits$coeff[,"M_rs"],
                             rho_s = additive_fits$coeff[,"Z"]) %>% 
  mutate(cum_mu_s = cumsum(mu_s),
         cum_alpha_s = cumsum(alpha_s),
         cum_beta_s = cumsum(beta_s),
         cum_rho_s = cumsum(rho_s))








##############################################################
# Estimating the survival functions

res_a1_astar1 <- surv_fn(aastar=1,
                         aa=1,
                         ZZ=67,
                         logreg_fits = logreg_fits,
                         additive_coefs = additive_coefs,
                         data_long = data_long)
plot(res_a1_astar1$time, res_a1_astar1$surv_prob, type="l", ylim=c(0,1))

res_a1_astar0 <- surv_fn(aastar=0,
                         aa=1,
                         ZZ=67,
                         logreg_fits = logreg_fits,
                         additive_coefs = additive_coefs,
                         data_long = data_long)
plot(res_a1_astar0$time, res_a1_astar0$surv_prob, type="l", ylim=c(0,1))

res_a0_astar0 <- surv_fn(aastar=0,
                         aa=0,
                         ZZ=67,
                         logreg_fits = logreg_fits,
                         additive_coefs = additive_coefs,
                         data_long = data_long)
plot(res_a0_astar0$time, res_a0_astar0$surv_prob, type="l", ylim=c(0,1))







# Saving all estimated and real effects, getting them ready for plotting
relative_effects <- rbind(
  data.frame(
    time = res_a1_astar0$time, 
    eff = res_a1_astar0$surv_prob / res_a0_astar0$surv_prob,
    type = "Direct effect",
    scale = "Relative survival"
  ),
  data.frame(
    time = res_a1_astar1$time, 
    eff = res_a1_astar1$surv_prob / res_a1_astar0$surv_prob,
    type = "Indirect effect",
    scale = "Relative survival"
  ),
  data.frame(
    time = res_a1_astar1$time, 
    eff = res_a1_astar1$surv_prob / res_a0_astar0$surv_prob,
    type = "Total effect",
    scale = "Relative survival"
  )
)

relative_effects_real <- rbind(
  data.frame(
    time = real_aa1_aastar0$time, 
    eff = real_aa1_aastar0$surv_prob / real_aa0_aastar0$surv_prob,
    type = "Direct effect",
    scale = "Relative survival"
  ),
  data.frame(
    time = real_aa1_aastar1$time, 
    eff = real_aa1_aastar1$surv_prob / real_aa1_aastar0$surv_prob,
    type = "Indirect effect",
    scale = "Relative survival"
  ),
  data.frame(
    time = real_aa1_aastar1$time, 
    eff = real_aa1_aastar1$surv_prob / real_aa0_aastar0$surv_prob,
    type = "Total effect",
    scale = "Relative survival"
  )
)

difference_effects <- rbind(
  data.frame(
    time = res_a1_astar0$time, 
    eff = res_a1_astar0$surv_prob - res_a0_astar0$surv_prob,
    type = "Direct effect",
    scale = "Survival difference"
  ),
  data.frame(
    time = res_a1_astar1$time, 
    eff = res_a1_astar1$surv_prob - res_a1_astar0$surv_prob,
    type = "Indirect effect",
    scale = "Survival difference"
  ),
  data.frame(
    time = res_a1_astar1$time, 
    eff = res_a1_astar1$surv_prob - res_a0_astar0$surv_prob,
    type = "Total effect",
    scale = "Survival difference"
  )
)

difference_effects_real <- rbind(
  data.frame(
    time = real_aa1_aastar0$time, 
    eff = real_aa1_aastar0$surv_prob - real_aa0_aastar0$surv_prob,
    type = "Direct effect",
    scale = "Relative survival"
  ),
  data.frame(
    time = real_aa1_aastar1$time, 
    eff = real_aa1_aastar1$surv_prob - real_aa1_aastar0$surv_prob,
    type = "Indirect effect",
    scale = "Relative survival"
  ),
  data.frame(
    time = real_aa1_aastar1$time, 
    eff = real_aa1_aastar1$surv_prob - real_aa0_aastar0$surv_prob,
    type = "Total effect",
    scale = "Relative survival"
  )
)









################################################################################
# Non-parametric bootstrap
# NB! This code is quite slow, so if one wishes to skip it and see the plots in the end without it,
# one can follow the instructions in the plotting commands in the end of this file

N <- 500 #number of bootstrap samples
timebyy <- 0.5

# Let's make this faster
cl <- makeCluster(2)
registerDoParallel(cl)

tic()
BS_samples <- foreach(iter=1:N,
                      .combine="rbind",
                      .packages=c("dplyr", "reshape2", "dtplyr", "foreach", "stringr")) %dopar% {
                        
                        
                        mydata <- mydata.reg[sample(nrow(mydata.reg),replace=T),] %>% 
                          mutate(patient = 1:nrow(mydata.reg))
                        
                        # Fitting the model for the mediator
                        M_fits <- mydata %>% 
                          pool_mediator_data(timeby=timebyy) %>% 
                          pooled_logreg() %>% 
                          pool_event_data(data=mydata,
                                          timeby=timebyy)
                        
                        
                        # Fitting the model for the survival outcome
                        additive_fits <- my.additive.new(regformula = "~ A + M_rs + Z",
                                                         event="E",
                                                         startt="int_begin",
                                                         stopt="int_end",
                                                         dataset=M_fits$pooled_event_data)
                        additive_fits$coeff[,"M_rs"] <- ifelse(is.na(additive_fits$coeff[,"M_rs"]),
                                                               0,
                                                               additive_fits$coeff[,"M_rs"]) #turn to 0 where indir effect is undefined
                        
                        
                        # Manipulate the results a bit and get them ready for estimating the survival functions
                        data_long <- M_fits$pooled_event_data
                        logreg_fits <- M_fits$logreg_coefs
                        obstimes.imput <- sort(unique(data_long$int_end[data_long$E==1] )) #get event times
                        
                        additive_coefs <- data.frame(time = obstimes.imput,
                                                     mu_s = additive_fits$coeff[,"(Intercept)"],
                                                     alpha_s = additive_fits$coeff[,"A"],
                                                     beta_s = additive_fits$coeff[,"M_rs"],
                                                     rho_s = additive_fits$coeff[,"Z"]) %>% 
                          mutate(cum_mu_s = cumsum(mu_s),
                                 cum_alpha_s = cumsum(alpha_s),
                                 cum_beta_s = cumsum(beta_s),
                                 cum_rho_s = cumsum(rho_s))
                        
                        
                        # Estimate the surv fn
                        res_a1_astar1 <- surv_fn(aastar=1,
                                                 aa=1,
                                                 ZZ=67,
                                                 logreg_fits = logreg_fits,
                                                 additive_coefs = additive_coefs,
                                                 data_long = data_long)
                        
                        res_a1_astar0 <- surv_fn(aastar=0,
                                                 aa=1,
                                                 ZZ=67,
                                                 logreg_fits = logreg_fits,
                                                 additive_coefs = additive_coefs,
                                                 data_long = data_long)
                        
                        res_a0_astar0 <- surv_fn(aastar=0,
                                                 aa=0,
                                                 ZZ=67,
                                                 logreg_fits = logreg_fits,
                                                 additive_coefs = additive_coefs,
                                                 data_long = data_long)
                        
                        
                        # Return results
                        data.frame(
                          time = res_a1_astar1$time,
                          
                          dir_rel = res_a1_astar0$surv_prob / res_a0_astar0$surv_prob,
                          indir_rel = res_a1_astar1$surv_prob / res_a1_astar0$surv_prob,
                          total_rel = res_a1_astar1$surv_prob / res_a0_astar0$surv_prob,
                          
                          dir_diff = res_a1_astar0$surv_prob - res_a0_astar0$surv_prob,
                          indir_diff = res_a1_astar1$surv_prob - res_a1_astar0$surv_prob,
                          total_diff = res_a1_astar1$surv_prob - res_a0_astar0$surv_prob
                        )
                      }


toc()
stopCluster(cl)


BS_quantiles <- BS_samples %>% 
  group_by(time) %>% 
  summarise(
    dir_rel_lower = quantile(dir_rel, 0.025),
    dir_rel_upper = quantile(dir_rel, 0.975),
    
    indir_rel_lower = quantile(indir_rel, 0.025),
    indir_rel_upper = quantile(indir_rel, 0.975),
    
    total_rel_lower = quantile(total_rel, 0.025),
    total_rel_upper = quantile(total_rel, 0.975),
    
    dir_diff_lower = quantile(dir_diff, 0.025),
    dir_diff_upper = quantile(dir_diff, 0.975),
    
    indir_diff_lower = quantile(indir_diff, 0.025),
    indir_diff_upper = quantile(indir_diff, 0.975),
    
    total_diff_lower = quantile(total_diff, 0.025),
    total_diff_upper = quantile(total_diff, 0.975)
  ) %>% 
  ungroup()


BS_quantiles_forggplot <-rbind(
  data.frame(
    time = BS_quantiles$time,
    lower = BS_quantiles$dir_rel_lower,
    upper = BS_quantiles$dir_rel_upper,
    type = "Direct effect",
    scale = "rel"
  ),
  
  data.frame(
    time = BS_quantiles$time,
    lower = BS_quantiles$indir_rel_lower,
    upper = BS_quantiles$indir_rel_upper,
    type = "Indirect effect",
    scale = "rel"
  ),
  
  data.frame(
    time = BS_quantiles$time,
    lower = BS_quantiles$total_rel_lower,
    upper = BS_quantiles$total_rel_upper,
    type = "Total effect",
    scale = "rel"
  ),
  
  data.frame(
    time = BS_quantiles$time,
    lower = BS_quantiles$dir_diff_lower,
    upper = BS_quantiles$dir_diff_upper,
    type = "Direct effect",
    scale = "diff"
  ),
  
  data.frame(
    time = BS_quantiles$time,
    lower = BS_quantiles$indir_diff_lower,
    upper = BS_quantiles$indir_diff_upper,
    type = "Indirect effect",
    scale = "diff"
  ),
  
  data.frame(
    time = BS_quantiles$time,
    lower = BS_quantiles$total_diff_lower,
    upper = BS_quantiles$total_diff_upper,
    type = "Total effect",
    scale = "diff"
  )
) %>% 
  mutate(eff=1)









##################################################################################
# Plotting the results
relative_effects %>% 
  ggplot(aes(time, eff)) +
  geom_step(size=0.5)+
  geom_ribbon(data=BS_quantiles_forggplot %>% filter(scale=="rel"), #If one wishes to skip bootstrapping, one can comment this and the next row out
              aes(ymin=lower, ymax=upper), alpha=0.3)+
  geom_line(data=relative_effects_real, aes(time, eff), color=colorval)+
  facet_grid( ~ type)+
  theme_bw()+
  theme(strip.background = element_rect("white"),
        text = element_text(family = "serif", size=13))+
  ylab("Relative survival")+
  xlab("Years")


difference_effects %>% 
  ggplot(aes(time, eff)) +
  geom_step(size=0.5)+
  geom_ribbon(data=BS_quantiles_forggplot %>% filter(scale=="diff"), #If one wishes to skip bootstrapping, one can comment this and the next row out
              aes(ymin=lower, ymax=upper), alpha=0.3)+
  geom_line(data=difference_effects_real, aes(time, eff), color=colorval)+
  facet_grid( ~ type)+
  theme_bw()+
  theme(strip.background = element_rect("white"),
        text = element_text(family = "serif", size=13))+
  ylab("Survival difference")+
  xlab("Years")









