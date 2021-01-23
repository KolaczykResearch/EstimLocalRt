## Part of these functions are adapt from R package EpiEstim , which is described in the following paper:
## A New Framework and Software to Estimate Time-Varying Reproduction Numbers During Epidemics 
## Anne Cori, Neil M. Ferguson, Christophe Fraser and Simon Cauchemez American Journal of Epidemiology 2013.

#' --------------------------------------------------------------------
#' Get overall the overall infectivity due to previously infected individuals.
#' 
#' @param incid - data frame with incid$local contains the incidence of cases due to local transmission and incid$imported contains the incidence of imported cases
#' @param si_distr - vector of probabilities giving the discrete distribution of the serial interval.
#' @return lambda - a vector which contains the overall infectivity at each time step
#'

getInfectivity <- function(incid, si_distr)
{
  T <- nrow(incid)
  lambda <- vector()
  lambda[1] <- NA
  for (t in seq(2, T))
  {
    lambda[t] <- sum(si_distr[seq_len(t)]*rowSums(incid[seq(t, 1), c("local", "imported")]), na.rm = TRUE)
  } 
  return(lambda)                                         
}

#' --------------------------------------------------------------------
#' Calculates the parameters of the Gamma posterior distribution for local/imported instantaneous reproduction number from the discrete SI distribution
#' 
#' @param incid - data frame with incid$local contains the incidence of cases due to local transmission and incid$imported contains the incidence of imported cases
#' @param si_distr - discrete distribution of the serial interval
#' @param mean_prior a positive number giving the mean of the common prior
#' distribution for all reproduction numbers
#' @param std_prior - a positive number giving the standard deviation of the
#' common prior distribution for all reproduction numbers
#' @param t_start - vector of positive integers giving the starting times of each
#' window over which the reproduction number will be estimated 
#' @param t_end - vector of positive integers giving the ending times of each
#' window over which the reproduction number will be estimated
#' @return posterior_para - list with posterior shape and scale parameters
#'
getR_posterior <- function(incid, si_distr, mean_prior=5, std_prior=5,
                           t_start, t_end,type="local") {
  a_prior <- (mean_prior / std_prior)^2
  b_prior <- std_prior^2 / mean_prior
  nb_time_periods <- length(t_start)
  lambda <- getInfectivity(incid, si_distr)
  a_posterior <- sapply(seq_len(nb_time_periods), function(t) 
    a_prior + sum(incid[seq(t_start[t], t_end[t]), type]))
  
  b_posterior <- sapply(seq_len(nb_time_periods), function(t) 
    1 / (1 / b_prior + sum(lambda[seq(t_start[t], t_end[t])], na.rm = TRUE)))
  posterior_para <- list(a_posterior, b_posterior)
  return(posterior_para)
}

#' --------------------------------------------------------------------
#' Calculates the true instantaneous reproduction number  
#' 
#' @param incid - data frame with incid$local contains the incidence of cases due to local transmission and incid$imported contains the incidence of imported cases
#' @param gi_distr - vector of probabilities giving the discrete distribution of the generation interval.
#' @return Rt_true - vector of true instantaneous reproduction number

getR_true <- function(incid, gi_distr)
{
  Rt_true <- incid$local/getInfectivity(incid, gi_distr)
  return(Rt_true)
}

