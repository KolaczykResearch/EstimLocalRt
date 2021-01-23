library(EpiEstim)
library(incidence)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(data.table)
library(grid)
library(gridExtra)
library(lemon)
library(foreach)
library(doParallel)
library(HDInterval)
library(nimble)
source("model_bayesian.R")
source("Rt_misc.R")

# set a directory for input
hdir <- "../../Results/Simulation_Epidemics/"

# set a directory for output
sdir <- "../../Results/Rt_Estimation/"
 
ntask <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if ( is.na(ntask) ) ntask <- 1


# set identification error rates for simulated  epidemics in Hong Kong and Victoria, Australia
if(ntask <= 6)
{
  int_name <- "hongkong"
  date_seq <- seq(as.Date("2020-01-18"), as.Date("2020-04-26"), by = "day")
} else {
  int_name <- "australia"
  date_seq <- seq(as.Date("2020-01-19"), as.Date("2020-04-15"), by = "day")
}

if(ntask %in% c(1,7) ) {
  zeta_0_true <- 2
  xi_0_true <- 18
  zeta_1_true <- 2
  xi_1_true <- 8
} else if(ntask %in% c(2,8) ){
  zeta_0_true <- 2
  xi_0_true <- 18
  zeta_1_true <- 4
  xi_1_true <- 8
} else if(ntask %in% c(3,9) ){
  zeta_0_true <- 2
  xi_0_true <- 18
  zeta_1_true <- 8
  xi_1_true <- 8
} else if(ntask %in% c(4,10) ){
  zeta_0_true <- 2
  xi_0_true <- 8
  zeta_1_true <- 2
  xi_1_true <- 18
} else if(ntask %in% c(5,11) ){
  zeta_1_true <- 2
  xi_1_true <- 18
  zeta_0_true <- 4
  xi_0_true <- 8
} else if(ntask %in% c(6,12) ){
  zeta_1_true <- 2
  xi_1_true <- 18
  zeta_0_true <- 8
  xi_0_true <- 8
}

# set seed
set.seed(1)
 
# number of simualtions
n_sims <- 1000

# number of cores for parallel computing
numCores <- 28
cl <- makeCluster(numCores[1]-1)
registerDoParallel(cl)
 
# 95% credible intervals
low_quantile <- 0.025
high_quantile <- 1- low_quantile
ci <- 1- 2*low_quantile 



# true generation interval (GI)
name_csv <- paste0("infection_counts_",int_name,".csv") 
dir_csv <- paste0(hdir,name_csv)
gi_distr <-  fread(dir_csv, header = T) %>% select(GI) %>% subset(GI>0) %>% group_by(GI)%>%
  summarise(freq=n() )  %>% mutate(prob=freq/sum(freq)) %>% select(GI,prob) %>% 
  complete(GI= 0:(length(date_seq)-1),fill = list(prob = 0))%>% select(prob) %>% unlist
  
 
# mu_local, mu_import
sim_results_mean <-  fread(dir_csv, header = T) %>% select(sim_num,dates,exogenous)%>%
  mutate(exogenous=case_when(exogenous ==1 ~ "imported",exogenous== 0 ~ "local"))    %>%
  group_by(sim_num,dates,exogenous) %>% summarise(count=n())   %>%
  ungroup %>% complete(sim_num= 0:(n_sims-1),dates= as.character(date_seq), exogenous,fill = list(count = 0))%>%
  mutate(dates=as.Date(dates))%>%group_by(dates,exogenous) %>% summarise(means=mean(count))%>% 
  spread(exogenous,  means) 

# N_local N_import
name_csv <- paste0("diagnosed_counts_",int_name,".csv") 
dir_csv <- paste0(hdir,name_csv)
sim_results_wide <- fread(dir_csv, header = T) %>% select(sim_num,dates,exogenous)%>%
  mutate(exogenous=case_when(exogenous ==1 ~ "imported",exogenous== 0 ~ "local"))    %>%
  group_by(sim_num,dates,exogenous) %>% summarise(count=n())   %>%
  ungroup %>% complete(sim_num= 0:(n_sims-1),dates= as.character(date_seq), exogenous,fill = list(count = 0))%>%
  mutate(dates=as.Date(dates)) %>%spread(exogenous, count)

# tilde N_local, tilde N_import
sim_results_obs <- sim_results_wide %>% mutate(local_obs = rbinom(length(local), local, 1-rbeta(length(local),zeta_1_true,xi_1_true))+
                                                   rbinom(length(imported), imported, rbeta(length(imported),zeta_0_true,xi_0_true) ),
                                                 imported_obs=local+imported-local_obs ) %>%select(sim_num,dates,local_obs,imported_obs) %>%
  rename(local=local_obs,imported=imported_obs)

# true instantaneous reproduction number (GI, infected)
res_true <- getR_true(sim_results_mean, gi_distr)
plot_true <- data.frame(dates=date_seq ,means=res_true,Type="True(GI, infected)", country=int_name)


# true serial interval (SI)
name_csv <- paste0("symptomatic_counts_",int_name,".csv") 
dir_csv <- paste0(hdir,name_csv)
si_distr <-  fread(dir_csv, header = T) %>% select(SI) %>% subset(SI>0) %>% group_by(SI)%>%
  summarise(freq=n() )  %>% mutate(prob=freq/sum(freq)) %>% select(SI,prob) %>% 
  complete(SI= 0:(length(date_seq)-1),fill = list(prob = 0))%>% select(prob) %>% unlist
  
# 'true' instantaneous reproduction number (SI, diagnosed)
name_csv <- paste0("diagnosed_counts_",int_name,".csv") 
dir_csv <- paste0(hdir,name_csv)
sim_results_mean_diagnosed <-  fread(dir_csv, header = T) %>% select(sim_num,dates,exogenous)%>%
  mutate(exogenous=case_when(exogenous ==1 ~ "imported",exogenous== 0 ~ "local"))    %>%
  group_by(sim_num,dates,exogenous) %>% summarise(count=n())   %>%
  ungroup %>% complete(sim_num= 0:(n_sims-1),dates= as.character(date_seq), exogenous,fill = list(count = 0))%>%
  mutate(dates=as.Date(dates))%>%group_by(dates,exogenous) %>% summarise(means=mean(count))%>% 
  spread(exogenous,  means) 
res_true_diagnosed <- getR_true(sim_results_mean_diagnosed, si_distr)
plot_true_diagnosed <- data.frame(dates=date_seq ,means=res_true_diagnosed,Type="True(SI, diagnosed)", country=int_name)
  

# estimate instantaneous reproduction number
est_res <- res_est <- plot_est <- list() 

# EpiEstim (no sliding)
tau <- 0
date_start <- 2 
mean_prior <- std_prior <- 1 
t_start <- date_start:(length(date_seq)-tau )
t_end <- t_start +tau 

name_csv <- paste0("symptomatic_counts_",int_name,".csv") 
dir_csv <- paste0(hdir,name_csv)
si_distr_csv <-  fread(dir_csv, header = T) 

est_res[[1]]  <- foreach (i=1:n_sims,.combine=rbind,.packages=c("EpiEstim","dplyr","tidyr")) %dopar% { 
  sim_num <- i-1
  si_distr <- si_distr_csv %>% subset(sim_num==(i-1)) %>% select(SI) %>% subset(SI>0) %>% group_by(SI)%>%
    summarise(freq=n() )  %>% mutate(prob=freq/sum(freq)) %>% select(SI,prob) %>% 
    complete(SI= 0:(length(date_seq)-1),fill = list(prob = 0))%>% select(prob) %>% unlist
  dates <- as.character(date_seq[t_end])
  incid <- sim_results_obs %>% subset(sim_num==(i-1)) %>% select(-sim_num)
  Rt_res <- getR_posterior(incid, si_distr, mean_prior=mean_prior, std_prior=std_prior,t_start, t_end,type="local")
  shape <- Rt_res[[1]]
  scale <- Rt_res[[2]]
  cbind(sim_num,dates,shape,scale)
}

# EpiEstim (weekly sliding)
tau <- 6
t_start <- date_start:(length(date_seq)-tau )
t_end <- t_start +tau 

est_res[[2]]  <- foreach (i=1:n_sims,.combine=rbind,.packages=c("EpiEstim","dplyr","tidyr")) %dopar% { 
  sim_num <- i-1
  si_distr <- si_distr_csv %>% subset(sim_num==(i-1)) %>% select(SI) %>% subset(SI>0) %>% group_by(SI)%>%
    summarise(freq=n() )  %>% mutate(prob=freq/sum(freq)) %>% select(SI,prob) %>% 
    complete(SI= 0:(length(date_seq)-1),fill = list(prob = 0))%>% select(prob) %>% unlist
  dates <- as.character(date_seq[t_end])
  incid <- sim_results_obs %>% subset(sim_num==(i-1)) %>% select(-sim_num)
  Rt_res <- getR_posterior(incid, si_distr, mean_prior=mean_prior, std_prior=std_prior,t_start, t_end,type="local")
  shape <- Rt_res[[1]]
  scale <- Rt_res[[2]]
  cbind(sim_num,dates,shape,scale)
}


# Bayesian model (MCMC)
niter <- 10000
nburnin <- 1000
bayes.mod.params <-  c("R_local")
n_days <- length(date_seq)

est_res[[3]]  <- foreach (i=1:n_sims,.combine=rbind,.packages=c("nimble","dplyr","tidyr","HDInterval"),
                     .export=c("dsumofbin", "rsumofbin","ddeg","rdeg")) %dopar% { 
                       source("model_bayesian.R")
                       sim_num <- i-1
                       si_distr <- si_distr_csv %>% subset(sim_num==(i-1)) %>% select(SI) %>% subset(SI>0) %>% group_by(SI)%>%
                         summarise(freq=n() )  %>% mutate(prob=freq/sum(freq)) %>% select(SI,prob) %>% 
                         complete(SI= 0:(length(date_seq)-1),fill = list(prob = 0))%>% select(prob) %>% unlist
                       dates <- as.character(date_seq)
                       incid <- sim_results_obs %>% subset(sim_num==(i-1)) %>% select(-sim_num)
                       ni <- sim_results_wide %>% subset(sim_num==(i-1)) %>% select(-sim_num)%>% summarise(imported=sum(imported),local=sum(local))
                       beta_est <- BetaMixture::BM_Fit(c(rbeta(ni$imported,zeta_0_true,xi_0_true),
                                                         rbeta(ni$local,xi_1_true,zeta_1_true) ), 2, 0.01, 1)
                       zeta_0_est <- beta_est$Alpha[1]
                       xi_0_est <- beta_est$Beta[1]
                       zeta_1_est<- beta_est$Beta[2]
                       xi_1_est <-beta_est$Alpha[2]
                       sim.data  <- list(I_local_obs=incid$local,I_imported_obs=incid$imported )
                       sim.constants <- list(si_distr=si_distr,  n_days = n_days,zeta_0=zeta_0_est,xi_0=xi_0_est, zeta_1=zeta_1_est,xi_1=xi_1_est)
                       alpha_0_rand <- zeta_0_est/(zeta_0_est+xi_0_est)
                       alpha_1_rand <- zeta_1_est/(zeta_1_est+xi_1_est)
                       I_local_true_rand <- ((1-alpha_0_rand)*incid$local-alpha_0_rand*incid$imported)/(1-alpha_0_rand-alpha_1_rand)
                       I_local_true_rand <- ifelse(I_local_true_rand<0, 0,I_local_true_rand) 
                       I_local_true_rand <- round(ifelse(I_local_true_rand>incid$local + incid$imported,  incid$local + incid$imported,I_local_true_rand) )
                       I_imported_true_rand <- incid$local + incid$imported - I_local_true_rand
                       initsFunction <- function() list(alpha_0=rbeta(1,zeta_0_est,xi_0_est) ,alpha_1=rbeta(1,zeta_1_est,xi_1_est),mu_imported=rpois(sim.constants$n_days,1),mu_local=rpois(1,1), 
                                                        R_local=rgamma(sim.constants$n_days,1,1),
                                                        I_local_true= I_local_true_rand,
                                                        I_imported_true = I_imported_true_rand )
                       initsList <- initsFunction()
                       Rmodel <- nimbleModel(code=bayes.mod, constants =sim.constants,data = sim.data, inits = initsList )
                       mcmc.out <- nimbleMCMC(model =  Rmodel, inits = initsList ,
                                              nchains = 1, niter = niter,nburnin = nburnin,
                                              summary = TRUE, samples= TRUE, monitors = bayes.mod.params, setSeed = TRUE)
                       mcmc_R_local <- mcmc.out$summary[paste0("R_local[",1:n_days,"]"),]
                       tmp <- sapply(1:n_days, function (j) hdi(mcmc.out$samples[,paste0("R_local[",j,"]")],ci)) 
                       ave <- mcmc_R_local[,1]
                       low <- tmp[1,]
                       high <- tmp[2,]
                       cbind(sim_num,dates,ave,low,high)
                     }




# compute 95% credible intervals for Rt local estimation
type_seq <- c("Misidentified_EpiEstim(no sliding)", "Misidentified_EpiEstim(weekly sliding)","Misidentified_Bayesian")

for(i in 1:2)
{
  res_est[[i]] <- as.data.frame(est_res[[i]] ) %>% mutate(sim_num=as.integer(as.character(sim_num)),
                                                          dates=as.Date(as.character(dates)),shape=as.numeric(as.character(shape)), scale=as.numeric(as.character(scale)))
  tmp_ci <- sapply(seq_len(length(res_est[[i]]$shape)), function(j) hdi(qgamma, ci, shape=res_est[[i]]$shape[j], scale=res_est[[i]]$scale[j]))
  res_est[[i]] <- res_est[[i]] %>% mutate(ave=shape*scale,low=tmp_ci[1,] ,high=tmp_ci[2,],value= 1-pgamma(1,shape=shape,scale=scale))

  plot_est[[i]]  <- res_est[[i]] %>% select(sim_num,dates,ave,low,high) %>% group_by(dates) %>%
    summarise(means=mean(ave,na.rm=TRUE),Low=mean(low,na.rm=TRUE),High=mean(high,na.rm=TRUE)) %>% mutate(Type=type_seq[i],country=int_name,
       zeta_0_true=zeta_0_true,xi_0_true=xi_0_true,zeta_1_true=zeta_1_true,xi_1_true=xi_1_true)
  
}

i <- 3
res_est[[i]] <- as.data.frame(est_res[[i]]) %>% mutate(sim_num=as.integer(as.character(sim_num)),
                                             dates=as.Date(as.character(dates)),ave=as.numeric(as.character(ave)), low=as.numeric(as.character(low)),high=as.numeric(as.character(high)))
                                             
plot_est[[i]]  <- res_est[[i]] %>% select(sim_num,dates,ave,low,high)%>%group_by(dates) %>%
  summarise(means=mean(ave,na.rm=TRUE),Low=mean(low,na.rm=TRUE),High=mean(high,na.rm=TRUE)) %>% mutate(Type=type_seq[i], country=int_name,
   zeta_0_true=zeta_0_true,xi_0_true=xi_0_true,zeta_1_true=zeta_1_true,xi_1_true=xi_1_true)


plot_est[[length(type_seq)+1]] <- plot_true
plot_est[[length(type_seq)+2]] <- plot_true_diagnosed
data_long <- Reduce(function(x, y) merge(x, y, all=TRUE), plot_est)  
 
write.csv(data_long,file = paste0(sdir,"simulation_",int_name,"_",ntask,".csv"),row.names = FALSE )

stopCluster(cl)
