library(tidyverse)
library(foreach)
set.seed(11)

source_location <- "I:/MSc/Aalen simulations/Aalen_et_al_time_dependent_mediators_in_survial_analysis/to_git/"
setwd(source_location)


# number of patients
n <- 20000
delta <- 1/52 # 52 is the number of weeks in a year!
endofstudy <- 13
time<-seq(delta, endofstudy, by=delta) #define timepoints
time_indexes <- 1:length(time)
pr <- function(lincomb) exp(lincomb) / (1+exp(lincomb))


# Baseline covariate (age)
Z <- rnorm(n=n, mean=67, sd=3)
hist(Z)

# treatment A and M0
A_0 <- rep(-66.5, times=n)
A <- rbinom(prob=pr(A_0 + 1*Z + rnorm(n)), n=n, size=1) #trt depending on Z
table(A)
M0 <- rep(0, n) #we assume that mediator had not happened at study start


# Effect of treatment 0 on M
eta_1 <- rep(-7.2, times=length(time))
plot(eta_1)


# effect of treatment 1 on M (plus baseline)
gamma <- rep(-1, times=length(time))
plot(gamma)


# effect of time
eta_i <- spline(c(0,4,7,10,13), c(0, 0.013, 0.05, 0.1, 0.13), xout=time)
eta_i <- eta_i$y
plot(eta_i)
eta_i.mat <- matrix(rep(eta_i, times=n), byrow=T, nrow=n)


# effect of baseline covariate Z on M
teeta <- rep(0.003, times=length(time))
plot(teeta)
teeta.mat <- matrix(rep(teeta, times=n), byrow=T, nrow=n)
Z.mat <- matrix(rep(Z, times=length(time)), byrow=F, nrow=n)



# Generate M matrix and noise
M_linkomb_nonoise <- eta_1 + eta_i.mat + A %*% t(gamma) + teeta.mat*Z.mat
# so for each timepoint i for each individual we have eta_1 + eta_i*I(i!=1) + gamma*A + teeta*Z here

Sigma <- matrix(0, length(time), length(time)) #some noise
diag(Sigma) <- 0.05
mu <- rep(0, length(time))
noise <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)



# Probabilities of mediator happening for each time point (if it has not happened yet)
M.mat_pr <- pr(M_linkomb_nonoise + noise)
M.mat <- M_linkomb_nonoise

#first column, based only on M0
pr_0to1 <- M.mat_pr[,1]
M.mat[,1] <- sapply(pr_0to1, rbinom, n=1, size=1)

for(j in 2:ncol(M.mat)){
  for(i in 1:nrow(M.mat)){
    if(M.mat[i,(j-1)]==1){ #if mediator happened at a previous time point
      M.mat_pr[i,j] <- 1 #prob of M(t)=1 if M(t-1)=1 is 1
      M.mat[i,j] <- 1 #mediator value stays 1
    } else{ #if mediator no happen at previous time point
      M.mat[i,j] <- rbinom(1, 1, M.mat_pr[i,j]) #generate event based on transition probability
    }
  }
}

test <- cbind(A, M.mat)
table(test[,1], test[,2]) #checking how many mediators happened by M1
table(test[,1], test[,length(time)]) #checking how many eventually got the mediator
prop.table(table(test[,1], test[,length(time)]), margin = 1)


# combine into a data set
patient <- seq(1, n)
sim.dat <- as.data.frame(cbind(patient, A, M0, M.mat, Z))
names(sim.dat)[4:(ncol(sim.dat) - 1)] <- paste0("M", time_indexes)


# Effect of A on hazard
alpha_s <- spline(c(0,2,4,8,13), c(-0.0173, -0.0135, -0.0125, -0.0115, -0.0106), xout=time)
plot(alpha_s)
alpha_s <- alpha_s$y


# Matrix with effect of A on hazard 
A.mat <- sim.dat$A %*% t(alpha_s)


# Effect of M on hazard
beta_s <- spline(c(0,2,10,13), c(0.08, 0.1, 0.15, 0.15)/40, xout=time)
plot(beta_s) 
beta_s <- beta_s$y
beta_s.mat <- matrix(rep(beta_s, times=n), byrow=T, nrow=n)


# Matrix with effects of M on hazard - elementwise product 
Meff.mat <- M.mat * beta_s.mat


# Effect of Z on hazard
rho_s <- spline(c(0,9,11,13), c(0.2, 0.23, 0.22999, 0.23)/(50*67), xout=time)
plot(rho_s) 
rho_s <- rho_s$y
rho_s.mat <- matrix(rep(rho_s, times=n), byrow=T, nrow=n)
Zeff.mat <- rho_s.mat * Z.mat


# Death probabilities for each tiny interval
hmat <- A.mat + Meff.mat + Zeff.mat #combine dir and indir and baseline cov
beta0 <- abs(min(hmat)) #baseline, to make sure no probability is negative
hmat <- hmat + beta0 
pmat <- hmat*delta # an approximated probability of the probability of an event within an interval of size delta so we end up with event probabilities for small intervals


p.mat.org <- cbind(1:n, pmat)
sim.dat$Time <- rep(time[length(time)], n) #repeat each time point for each individual
p.mat <- p.mat.org


# Generate death events by probabilites
for(t in time_indexes) { #for each time point
  t.i <- time[t]
  p.vec <- p.mat[, (t+1)] #+1 because of the id indexes (vector of all observations at week t)
  E <- rep(0, length(p.vec))
  id <- rep(0, length(p.vec))
  
  for(j in 1:length(p.vec)) { #for each patient
    id[j] <- p.mat[j, 1]
    E[j] <- sample(c(0,1), 1, prob=c(1-p.vec[j], p.vec[j])) #event samples -- for each observation at each time interval
  }
  
  if(any(E==1)) { #if the event happens to any of the patients
    id.event <- id[E==1] #mark down the patient id-s where events happen
    sim.dat[sim.dat$patient %in% id.event, "Time"] <- time[t] #check for each patient if their id is in the event id vector, and then overwrite the times (if event happens) -- if event doesn't happen, this stays 13
    p.mat <- p.mat[E==0, ] #select patients for who the event didn't happen yet
  }
  #checking the process
  print(t)
  print(dim(p.mat))
  print(sum(E))
}

sim.dat$E <- 1*(sim.dat$Time<endofstudy) #set the event indicator to 1 if it happened before the end of study


# Censoring
ff <- runif(n, 0, 28)
sim.dat$C.Time <- ff
sim.dat$C <- 1*(sim.dat$C.Time < sim.dat$Time) #if censored before event time (or study end time), then set to 1
sim.dat[sim.dat$C==1, "E"] <- 0 #set event indicator to 0 for those which got censored
sim.dat[sim.dat$C==1, "Time"] <- sim.dat[sim.dat$C==1, "C.Time"]  #set time to event to time to censoring time for those that got censored


# Reshaping data to suited format
mydata <- reshape(sim.dat, paste0("M", time_indexes),
                  times = time-delta, # measurement times from first mediator measurement
                  direction="long", v.names="M")


# Ordering
mydata <- mydata[ order(mydata[,"patient"], mydata[,"time"]), ]
row.names(mydata) <- 1:nrow(mydata)
mydata <- mydata[, c("patient","A","M0", "M","Time","C", "E", "Z", "time")]
names(mydata)[9] <- "start"
mydata$stopt <- rep(time, n)

m.split <- split(mydata, mydata$patient) #divide data to n groups based on patient id-s

hax <- function(df) { #for each patient:
  df$E <- 1*(df$Time <= df$stopt & sum(df$E) == length(time) ) #if event time is before interval end and sum of event indicators is 8 (because of 8 intervals), then E=1
  df$C <- 1*(df$Time <= df$stopt & sum(df$C) == length(time) )#if event time is before interval end and sum of censor indicators is 8 (because of 8 intervals), then C=1
  df$stopt <- pmin(c(df$stopt), df$Time) #set interval stop to be event/censor time or the initial interval stop time
  keep <- 1:min(c(which(df$E==1), length(time), which(df$C==1))) #select all the intervals before event/censoring (or all the intervals)
  df <- df[keep, ]
}

m.split.hax <- lapply(m.split, hax)
mydata.hax <- do.call("rbind", m.split.hax)


mydata.hax[mydata.hax$E==1, "stopt"] <- mydata.hax[mydata.hax$E==1, "stopt"] + 
  runif(length(mydata.hax[mydata.hax$E==1, "stopt"]), 0, 0.0005) #add some noise for the interval stop times for those where the event happened

sim.dat.long <- mydata.hax[, c("patient","A","M0", "M","C","E", "Z", "start","stopt")]
mydata.fin <- sim.dat.long



# transform data to regular format
lastrow <- mydata.fin %>% 
  group_by(patient) %>% 
  filter(stopt == max(stopt)) %>% #select the last time row for each observation
  rename(study_exit = stopt) %>% 
  select(patient, A, E, study_exit, Z) %>% 
  ungroup()

M_status <- reshape2::melt(M.mat) %>% 
  mutate(M_time = max(time)/length(time)*Var2) %>% 
  rename(patient = Var1,
         M = value) %>% 
  group_by(patient) %>% 
  mutate(M_status = max(M)) %>% 
  filter(M_status == M) %>%
  filter(M_time == min(M_time)) %>% 
  mutate(M_time = ifelse(M_status==1, M_time, endofstudy)) %>% # get mediator event time
  ungroup() %>% 
  select(!c(Var2, M_status))

mydata.reg <- merge(lastrow, M_status, by = "patient") %>% 
  mutate(M = ifelse(study_exit < M_time, 0, M),
         M_time = ifelse(study_exit <= M_time, study_exit, M_time))
save(mydata.reg, file= "mydata.reg.rda")

# check data quality by some frequency tables
hist(mydata.reg$M_time)
hist(mydata.reg$study_exit)
prop.table(table(mydata.reg$M, mydata.reg$E),1)
prop.table(table(mydata.reg$A, mydata.reg$M),1)
prop.table(table(mydata.reg$A, mydata.reg$E),1)

# censored people
mydata.reg_cens <- mydata.reg %>% 
  filter(E==0)
hist(mydata.reg_cens$M_time, probability=T, breaks = 20)
hist(mydata.reg_cens$study_exit, probability=T, breaks = 20)

# dead people
mydata.reg_dead <- mydata.reg %>% 
  filter(E==1)
hist(mydata.reg_dead$M_time, probability=T, breaks = 20)
hist(mydata.reg_dead$study_exit, probability=T, breaks = 20)

#censoring percentage
round(nrow(mydata.reg_cens)/nrow(mydata.reg)*100,2)


#censoring before study end
mydata.reg_cens_befend <- mydata.reg_cens %>% 
  filter(study_exit < 13)

#censoring percentage
round(nrow(mydata.reg_cens_befend)/nrow(mydata.reg)*100,2)











# FIND REAL EFFECTS --------------------------------------------------------------------
realeffects <- data.frame(time,
                          mu_s = beta0,
                          alpha_s = alpha_s,
                          beta_s = beta_s,
                          rho_s = rho_s,

                          eta_1 = eta_1,
                          eta_s = eta_i,
                          gamma = gamma,
                          teeta = teeta,
                          delta)

save(realeffects, file="realeffects.rda")




# Calculating the survival function -----------------------------------

surv_fn_gene <- function(aastar,
                         aa,
                         ZZ,
                         timestart = 0,
                         timeend = 13,
                         timebyy = delta,
                         realeffects){

  timepoints <- c(timestart, realeffects$time)

  # Getting mediator transition 0->1 at time t probabilities
  probs_aastar_all <- realeffects %>%
    arrange(time) %>%
    mutate(aa = aa,
           aastar = aastar,

           p_pred_aa = exp(eta_1 + eta_s + gamma * aa + teeta * ZZ) / (1+exp(eta_1 + eta_s + gamma * aa + teeta * ZZ)),
           p_pred_aastar = exp(eta_1 + eta_s + gamma * aastar + teeta * ZZ) / (1+exp(eta_1 + eta_s + gamma * aastar + teeta * ZZ)),

           int_index = 1:n()) %>%
    select(int_index, p_pred_aa, p_pred_aastar) %>%
    filter(time <= timeend)


  # All estimated additive hazards model coefficients
  all_coefs <- realeffects %>%
    mutate(a = aa,
           Z = ZZ,
           cum_mu_s = cumsum(mu_s),
           cum_alpha_s = cumsum(alpha_s),
           cum_beta_s = cumsum(beta_s),
           cum_rho_s = cumsum(rho_s),
           int_index = 1:n()) %>%
    group_by(time) %>%
    mutate(t_M = max(timepoints[timepoints <= time]),
           t_Mplus1 = t_M+timebyy) %>%
    ungroup()


  foreach(tgrp = unique(all_coefs$int_index),
          .combine = "rbind"#,
          #.packages = c("dplyr", "reshape2", "stringr")
  ) %do% {
    print(tgrp)
    t <- all_coefs[all_coefs$int_index==tgrp,"time"]

    #Select estimated probabilities for 0->1 transitions
    probs_aastar <- probs_aastar_all %>%
      #filter(int_index <= tgrp+1) %>%
      filter(int_index <= tgrp+1) %>%
      select(!p_pred_aa)

    probs_aa <- probs_aastar_all %>%
      filter(int_index <= tgrp) %>%
      select(!p_pred_aastar)


    ####### Component 1 ---------------------------------------------------------------------------------------
    all_in_interval <- all_coefs %>%
      filter(min(t) <= time & time <= max(t))

    comp1 <- exp(-(all_in_interval$cum_mu_s + all_in_interval$cum_alpha_s * all_in_interval$a + all_in_interval$cum_rho_s * all_in_interval$Z)*all_in_interval$delta)


    ####### Component 2-3 (the sum over all mediator combinations) ---------------------------------------------
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
    
    if(nrow(meds)==2){
      mediator_combinations <- data.frame(upper.tri(meds, diag=T) * meds) %>% 
        mutate(scenario = 1:n()) %>%  #add scenario id
        relocate(scenario) %>% 
        reshape2::melt(id="scenario") %>% 
        rename(t_M = variable,
               M = value) %>% 
        mutate(t_M = as.numeric(str_split_i(t_M, "X", 2))) %>% 
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
        reshape2::melt(id="scenario") %>%
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
        mutate(prob_Mto1 = case_when(
          Mminus1 == 0 ~ p_pred_aastar,
          Mminus1 == 1 ~ 1),

          prob_Mtowhatever = case_when(
            Mminus1 == 1 ~ 1,
            Mminus1 == 0 & M == 0 ~ 1-p_pred_aastar,
            Mminus1 == 0 & M == 1 ~ p_pred_aastar),

          prob_Mtowhatever_Mbeforet = ifelse(int_index==max(int_index), 1, prob_Mtowhatever)) %>%
        group_by(scenario) %>%
        mutate(scenario_prob = prod(prob_Mtowhatever_Mbeforet)) %>%
        ungroup() %>%
        arrange(scenario) %>%
        select(!c(prob_Mtowhatever, prob_Mtowhatever_Mbeforet))

      #sum(unique(mediator_combinations$scenario_prob)) #this should be 1

      all_upto_lastininterval <- all_coefs %>%
        filter(time <= max(t))

      all_upto_lastininterval %>%
        slice(rep(1:n(), nrow(all_in_interval))) %>%
        mutate(target_time = rep(all_in_interval$time, each=nrow(all_upto_lastininterval))) %>%
        slice(rep(1:n(), length(unique(mediator_combinations$scenario)))) %>%
        mutate(scenario = rep(unique(mediator_combinations$scenario), each=nrow(all_in_interval)*nrow(all_upto_lastininterval))) %>%
        filter(time <= target_time) %>%
        left_join(mediator_combinations, by=c("scenario"="scenario", "int_index"="int_index")) %>%
        mutate(comp2 = -beta_s * Mminus1 * timebyy) %>% #should be prob0to1 here if we plant it on the probability
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



    }
  }
}


# Checking the survival curves
real_aa1_aastar0 <- surv_fn_gene(aastar=0,
                                 aa=1,
                                 ZZ=67,
                                 realeffects = realeffects)
save(real_aa1_aastar0, file="real_aa1_aastar0.rda")
plot(real_aa1_aastar0$time, real_aa1_aastar0$surv_prob, type="l", ylim=c(0,1))

real_aa1_aastar1 <- surv_fn_gene(aastar=1,
                                 aa=1,
                                 ZZ=67,
                                 realeffects = realeffects)
save(real_aa1_aastar1, file= "real_aa1_aastar1.rda")
plot(real_aa1_aastar1$time, real_aa1_aastar1$surv_prob, type="l", ylim=c(0,1))

real_aa0_aastar0 <- surv_fn_gene(aastar=0,
                                 aa=0,
                                 ZZ=67,
                                 realeffects = realeffects)
save(real_aa0_aastar0, file="real_aa0_aastar0.rda")
plot(real_aa0_aastar0$time, real_aa0_aastar0$surv_prob, type="l", ylim=c(0,1))





