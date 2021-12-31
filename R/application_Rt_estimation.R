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


getBetaPara <- function(mu,sigma)
{
  k=(1-mu)/mu
  a=mu*((1-mu)*mu/sigma^2-1)
  b=a*k
  return(c(a,b))
}

# set a directory for output
sdir <- "../Results/Rt_Estimation/"

ntask <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if ( is.na(ntask) ) ntask <- 1

# set a directory for input
if(ntask == 1)
{
  int_name <- "hongkong"
  hdir <- "../Data/HongKong/"
} else {
  int_name <- "australia"
  hdir <- "../Data/Australia/"
}

 
name_csv <- "infection_counts_obs.csv"
dir_csv <- paste0(hdir,name_csv)
infection_counts_obs  <-  fread(dir_csv, header = T) %>% mutate(dates=as.Date(dates))
date_seq <- infection_counts_obs$dates

# set seed 
set.seed(1)

# serial interval 
si_distr <- discr_si(0:(length(date_seq)-1),2.23/0.37,sqrt(2.23/0.37^2) )

# 95% credible intervals
low_quantile <- 0.025
high_quantile <- 1- low_quantile
ci <- 1- 2*low_quantile 

est_res <- res_est <- plot_est <- list() 

# no identification error: EpiEstim (no sliding)
tau <- 0
date_start <- 2
t_start <- date_start:(length(date_seq)-tau )
t_end <- t_start +tau 
mean_prior <- std_prior <- 1 

dates <- as.character(date_seq[t_end])
incid <- infection_counts_obs
Rt_res <- getR_posterior(incid, si_distr, mean_prior=mean_prior, std_prior=std_prior,t_start, t_end,type="local")
shape <- Rt_res[[1]]
scale <- Rt_res[[2]]
est_res[[1]]  <- cbind(dates,shape,scale)

# error 1: alpha0 ~ .1, alpha1 ~ .3
mean_0_est <- .1
mean_1_est <- .3
sd_0_est <- sd_1_est <- .02
para_0 <- getBetaPara(mean_0_est,sd_0_est)
para_1 <- getBetaPara(mean_1_est,sd_1_est)
zeta_0_est <- para_0[1]
xi_0_est <- para_0[2]
zeta_1_est <- para_1[1]
xi_1_est <- para_1[2]
 

niter <- 10000
nburnin <- 1000
bayes.mod.params <-  c("R_local")
n_days <- length(date_seq)

dates <- as.character(date_seq)
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
est_res[[2]]  <- cbind(dates,ave,low,high)


# error 2: alpha0 ~ .3, alpha1 ~ .1
mean_0_est <- .3
mean_1_est <- .1
sd_0_est <- sd_1_est <- .02
para_0 <- getBetaPara(mean_0_est,sd_0_est)
para_1 <- getBetaPara(mean_1_est,sd_1_est)
zeta_0_est <- para_0[1]
xi_0_est <- para_0[2]
zeta_1_est <- para_1[1]
xi_1_est <- para_1[2]

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
est_res[[3]]  <- cbind(dates,ave,low,high)

# compute 95% credible intervals
type_seq <- c("No misidentification", "Error 1", "Error 2")

i=1
res_est[[i]] <- as.data.frame(est_res[[i]] ) %>% mutate(dates=as.Date(as.character(dates)),shape=as.numeric(as.character(shape)), scale=as.numeric(as.character(scale)))
tmp_ci <- sapply(seq_len(length(res_est[[i]]$shape)), function(j) hdi(qgamma, ci, shape=res_est[[i]]$shape[j], scale=res_est[[i]]$scale[j]))
res_est[[i]] <- res_est[[i]] %>% mutate(ave=shape*scale,low=tmp_ci[1,] ,high=tmp_ci[2,],value= 1-pgamma(1,shape=shape,scale=scale))
plot_est[[i]]  <- res_est[[i]] %>% select(dates,ave,low,high) %>% group_by(dates) %>%
    summarise(means=mean(ave,na.rm=TRUE),Low=mean(low,na.rm=TRUE),High=mean(high,na.rm=TRUE)) %>% mutate(Type=type_seq[i],country=int_name)

for(i in 2:3){
  res_est[[i]] <- as.data.frame(est_res[[i]]) %>% mutate(dates=as.Date(as.character(dates)),ave=as.numeric(as.character(ave)), low=as.numeric(as.character(low)),high=as.numeric(as.character(high)))
  
  plot_est[[i]]  <- res_est[[i]] %>% select(dates,ave,low,high)%>%group_by(dates) %>%
    summarise(means=mean(ave,na.rm=TRUE),Low=mean(low,na.rm=TRUE),High=mean(high,na.rm=TRUE)) %>% mutate(Type=type_seq[i], country=int_name)
}


data_long <- Reduce(function(x, y) merge(x, y, all=TRUE), plot_est)  

write.csv(data_long,file = paste0(sdir,"application_",int_name,".csv"),row.names = FALSE )



