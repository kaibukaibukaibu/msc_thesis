######################################################################################
############################## Helping functions #####################################
######################################################################################

## To use these functions the raw data must be in given format:
# patient - numeric unique with values 1,...,n
# A - binary numeric with values 0 and 1
# Z - numeric
# E - binary numeric with values 0 and 1
# study_exit - numeric
# M - binary numeric with values 0 and 1
# M_time - numeric



# Turn data into long format for the pooled logistic regression model
pool_mediator_data <- function(data,
                               timestart=0,
                               timeend=13,
                               timeby=.5){
  
  timepoints <- seq(timestart, timeend, timeby)
  
  data %>% 
    slice(rep(1:n(), length(timepoints)-1)) %>% #rep as many times as nr of intervals
    mutate(int_index = rep(1:(length(timepoints)-1), each=nrow(data))) %>% 
    mutate(int_begin = rep(timepoints[1:(length(timepoints)-1)], each=nrow(data)),
           int_end = rep(timepoints[2:(length(timepoints))], each=nrow(data))) %>% 
    mutate(M = ifelse(int_begin > M_time | M_time >= int_end, 0, M)) %>% #if no M yet, turn to 0
    filter(M_time >= int_begin) #choose rows up to M=1
           #!(int_begin <= M_time & M_time < int_end & M==0)) #if censoring happens, consider its time to be the last interval break before it
}




# Pooled logistic regression model
# Returns estimated values for E[M_t | M_{t-1}=0, A=a, Z=z] (so for each time point and treatment levels)
pooled_logreg <- function(pooled_data){
  
  fit <- glm(M ~ as.factor(int_index) + A + Z, data=pooled_data, family=binomial)
  
  newdata <- data.frame(A = rep(c(0,1), each=length(unique(pooled_data$int_index))),
                        int_index = rep(unique(pooled_data$int_index), times=length(unique(pooled_data$A))),
                        int_begin = rep(unique(pooled_data$int_begin), times=length(unique(pooled_data$A))),
                        int_end = rep(unique(pooled_data$int_end), times=length(unique(pooled_data$A))))
  
  nint <- length(unique(pooled_data$int_index))
  baseline_data <- pooled_data %>% 
    group_by(patient) %>% 
    filter(int_index == min(int_index)) %>% 
    select(patient, Z, A) %>% 
    ungroup() %>% 
    arrange(patient)
  
  newdata2 <- data.frame(patient = rep(unique(pooled_data$patient), each=nint),
                         int_index = rep(unique(pooled_data$int_index), times=nrow(baseline_data)),
                         int_begin = rep(unique(pooled_data$int_begin), times=nrow(baseline_data)),
                         int_end = rep(unique(pooled_data$int_end), times=nrow(baseline_data)),
                         #Z = rep(median(baseline_data$Z), each=nrow(baseline_data)*nint), 
                         Z = rep(baseline_data$Z, each=nint), 
                         A = rep(baseline_data$A, each=nint))
  
  p_j2 <- c(predict(fit, newdata2, type = "response")) #probabilities for each patient for each time interval (even tho they might die before)
  
  coefs <- fit$coefficients

  res <- newdata2 %>% 
    bind_cols(p_j2)
  
  names(res)[(ncol(newdata2)+1)] <- "p"
  
  res <- res %>%
    arrange(int_index) %>% 
    select(patient, int_index, int_begin, int_end, A, Z, p)
  
  list(res = res,
       logreg_coefs = fit)
}






# Turn data into long format for the event
pool_event_data <- function(M_probs,
                            data,
                            timestart=0,
                            timeend=13,
                            timeby=0.5){
  
  logreg_coefs <- M_probs$logreg_coefs
  M_probs <- M_probs$res
  
  timepoints <- seq(timestart, timeend, timeby)
  
  t <- data %>% 
    slice(rep(1:n(), length(timepoints)-1)) %>% #rep as many times as nr of intervals
    mutate(int_index = rep(1:(length(timepoints)-1), each=nrow(data))) %>% 
    mutate(int_begin = rep(timepoints[1:(length(timepoints)-1)], each=nrow(data)),
           int_end = rep(timepoints[2:(length(timepoints))], each=nrow(data))) %>% 
    mutate(E = ifelse(int_begin > study_exit | study_exit >= int_end, 0, E)) %>% #if no E yet, turn to 0
    filter(study_exit >= int_begin) %>%  #choose rows up to E=1
    group_by(patient) %>% 
    mutate(int_end = ifelse(int_end==max(int_end), study_exit, int_end)) %>% 
    ungroup() %>% 
    
    group_by(patient) %>% 
    mutate(M = ifelse((int_begin <= M_time & M_time < int_end) | 
                        (M_time <= int_begin & M_time <= int_end), 1, 0)) %>% 

    ungroup() %>% 
    arrange(patient, int_index)
  

  t$M_rs <- if_else( #M_{r(s)} -- mediator value at the beginning of the interval
    t$patient == lag(t$patient),
    lag(t$M, default = 0), 0
  )
  
  
  t <- t %>% 
    mutate(M_rs = ifelse(is.na(M_rs), 0, M_rs)) 
  
  
  
  list(pooled_event_data = t,
       logreg_coefs = logreg_coefs,
       M_probs = M_probs)
}

  



# Additive hazards model (similar as in the code provided in the article by Aalen et al. (2019))
# This function returns mu_t, alpha_t, beta_t and rho_t
my.additive.new = function(regformula,event,startt, stopt, dataset){
  
  obestime_A=dataset[[stopt]][dataset[[event]]==1] #select all stopt times where the event happened
  obstimes = sort(unique(obestime_A)) #order all event times uniquely (no repeated times)
  maxxer = max(obstimes)
  
  b.matrix = foreach(i = obstimes,
                     .combine = "rbind") %do% { #Loop doing a linear regression at each death time point:

    print(paste0(round(which(obstimes==i)/length(obstimes)*100,2), "%"))
                       
    # select the subset for subjects under risk at time i: 
    sta<-dataset[[startt]] #start times of time intervals
    sto<-dataset[[stopt]] #stop times of time intervals
    obstimei = subset(dataset, sta<i & sto>=i) #take all observations where start < i <= stop, so all people under risk at that specific time
    
    E = as.numeric(obstimei[[stopt]] == i & obstimei[[event]]==1) #selecting all the events that happened at time i (in the end of the intervals)
    
    obstimei$event <- ifelse(obstimei[[stopt]]==i & obstimei[[event]]==1, 1, 0)
    
    fitter <- lm(paste("event", regformula), data = obstimei)
    
    coef(fitter)
  }

  #Output:
  list(coeff = b.matrix)
}








# Survival function by specific a and a* values (they should be 0/1 numeric values)
# Returns surv fn estimates for all event (death) times
surv_fn <- function(aastar,
                    aa,
                    ZZ,
                    timestart = 0,
                    timeend = 13,
                    timebyy = 0.5,
                    logreg_fits,
                    additive_coefs=additive_coefs,
                    data_long){
  
  timepoints <- seq(timestart, timeend, by=timebyy)
  event_times <- sort(unique(data_long$int_end[data_long$E==1]))
  
  # Getting estimated pooled logreg probabilities
  probs_aastar_all <- data_long %>% 
    filter(E == 1) %>% 
    mutate(A = aastar,
           Z= ZZ) %>% 
    select(A, Z, int_index)
  probs_aastar_all$p_pred_aastar <- c(predict(logreg_fits, probs_aastar_all, type = "response"))
  probs_aastar_all <- probs_aastar_all %>% 
    mutate(A = aa)
  probs_aastar_all$p_pred_aa <- c(predict(logreg_fits, probs_aastar_all, type = "response"))
  probs_aastar_all <- probs_aastar_all %>% 
    unique() %>% 
    select(!c(Z, A)) %>% 
    #rbind(c(1,0,0)) %>% #since last measured mediator is M_0=0
    arrange(int_index)
  
  # All estimated additive hazards model coefficients
  all_coefs <- additive_coefs %>% 
    mutate(a = aa,
           Z = ZZ) %>% 
    group_by(time) %>% 
    mutate(t_M = max(timepoints[timepoints <= time]),
           t_Mplus1 = t_M+timebyy,
           int_index = t_Mplus1 / timebyy) %>%
    ungroup() 
  
  
  foreach(tgrp = as.numeric(unique(all_coefs$int_index)),
          .combine = "rbind"#,
          #.packages = c("dplyr", "reshape2", "stringr")
  ) %do% {
    t <- all_coefs[all_coefs$int_index==tgrp,"time"]
    
    #Select estimated probabilities for 0->1 transitions of the mediators
    probs_aastar <- probs_aastar_all %>% 
      filter(int_index * timebyy <= min(t) + timebyy) %>% 
      select(!p_pred_aa)
    
    probs_aa <- probs_aastar_all %>% 
      filter(int_index * timebyy <= min(t)) %>% 
      select(!p_pred_aastar)
    
    
    ####### Component 1 -----------------------------------------------
    all_in_interval <- all_coefs %>% 
      filter(int_index == tgrp)
    
    comp1 <- exp(-(all_in_interval$cum_mu_s + all_in_interval$cum_alpha_s * all_in_interval$a + all_in_interval$cum_rho_s * all_in_interval$Z))
    
    
    ####### Component 2-3 ---------------------------------------------
    # If in the first interval, there is only one measured mediator (M_0=0 for all), so combination probability can only be 1
    if(tgrp == unique(all_coefs$int_index)[1]){
      all_in_interval %>% 
        arrange(time) %>% 
        mutate(surv_prob = comp1) %>% 
        select(time, surv_prob)
      
    } else { # If not in the first interval
      meds <- matrix(1,
                     nrow=nrow(probs_aastar)-1,
                     ncol=nrow(probs_aastar)-1)
      meds <- rbind(meds,
                    rep(0, times=nrow(probs_aastar)-1)) #add scenario where M=1 never happens
      
      if(nrow(meds)==2){
        mediator_combinations <- as.data.frame(meds) %>% 
          mutate(scenario = 1:n()) %>%  #add scenario id
          relocate(scenario) %>% 
          melt(id="scenario") %>% 
          rename(t_M = variable,
                 M = value) %>% 
          mutate(t_M = as.numeric(str_split_i(t_M, "V", 2))) %>% 
          arrange(scenario, t_M) %>% 
          mutate(t_Mminus1 = t_M-1,
                 Mminus1 = ifelse(t_M==1, 0, M),
                 int_index = t_M) %>% 
          relocate(scenario, t_Mminus1, t_M, Mminus1, M) %>% 
          left_join(probs_aastar, by=c("int_index"="int_index")) #probs of 0->1 for mediator
      } else {
        
        mediator_combinations <- data.frame(upper.tri(meds, diag=T) * meds) %>% 
          mutate(scenario = 1:n()) %>%  #add scenario id
          relocate(scenario) %>% 
          melt(id="scenario") %>% 
          rename(t_M = variable,
                 M = value) %>% 
          mutate(t_M = as.numeric(str_split_i(t_M, "X", 2))) %>% 
          arrange(scenario, t_M) %>% 
          mutate(t_Mminus1 = t_M-1,
                 Mminus1 = ifelse(t_M==1, 0, M),
                 int_index = t_M) %>% 
          relocate(scenario, t_Mminus1, t_M, Mminus1, M) %>% 
          left_join(probs_aastar, by=c("int_index"="int_index")) #probs of 0->1 for mediator
      }
      
      mediator_combinations$Mminus1 <- if_else(
        mediator_combinations$scenario == lag(mediator_combinations$scenario),
        lag(mediator_combinations$M, default = 0), 0
      )
      mediator_combinations$Mminus1 <- if_else(is.na(mediator_combinations$Mminus1), 0, mediator_combinations$Mminus1)
      
      mediator_combinations <- mediator_combinations %>% 
        mutate(
          prob_Mtowhatever = case_when(
            Mminus1==0 & M==0 ~ 1-p_pred_aastar,
            Mminus1==0 & M==1 ~ p_pred_aastar,
            Mminus1==1 ~ 1)) %>% 
        group_by(scenario) %>% 
        mutate(scenario_prob = prod(prob_Mtowhatever)) %>% 
        ungroup() %>% 
        arrange(scenario) %>% 
        select(!c(prob_Mtowhatever, prob_Mtowhatever)) %>% 
        mutate(int_index = int_index+1)
      
      mediator_combinations <- mediator_combinations %>%   
        rbind(data.frame(
          scenario = 1:nrow(probs_aastar),
          t_Mminus1 = rep(NA, times=nrow(probs_aastar)),
          t_M = rep(0, times=nrow(probs_aastar)),
          Mminus1 = rep(0, times=nrow(probs_aastar)),
          M = rep(0, times=nrow(probs_aastar)),
          int_index = rep(1, times=nrow(probs_aastar)),
          p_pred_aastar = rep(1, times=nrow(probs_aastar)),
          scenario_prob = unique(mediator_combinations$scenario_prob)
        )) %>% 
        arrange(scenario, int_index) %>% 
        mutate(int_index =as.factor(int_index),
               scenario = as.factor(scenario))
      
      
      #sum(unique(mediator_combinations$scenario_prob)) #this should be 1
      
      all_upto_lastininterval <- all_coefs %>% 
        filter(time <= max(t))
      
      k=all_upto_lastininterval %>% 
        slice(rep(1:n(), nrow(all_in_interval))) %>% 
        mutate(target_time = rep(all_in_interval$time, each=nrow(all_upto_lastininterval))) %>% 
        slice(rep(1:n(), length(unique(mediator_combinations$scenario)))) %>% 
        mutate(scenario = rep(unique(mediator_combinations$scenario), each=nrow(all_in_interval)*nrow(all_upto_lastininterval))) %>% 
        filter(time <= target_time) %>%
        mutate(int_index =as.factor(int_index),
               scenario = as.factor(scenario)) %>% 
        left_join(mediator_combinations, by=c("scenario"="scenario", "int_index"="int_index")) %>% 
        mutate(comp2 = -beta_s * M) %>% 
        group_by(scenario, target_time) %>% 
        summarise(comp23 = exp(sum(comp2)),
                  scenario_prob = unique(scenario_prob),
                  time = unique(target_time)) %>% 
        ungroup() %>% 
        mutate(comp23 = comp23 * scenario_prob) %>% 
        group_by(time) %>% 
        summarise(comp23 = sum(comp23),
                  time = unique(time)) %>% 
        ungroup() %>% 
        mutate(comp1 = comp1,
               surv_prob = comp1 * comp23) %>% 
        select(time, surv_prob)
      
      
      k
    }
  }
}

##############################################################
##############################################################