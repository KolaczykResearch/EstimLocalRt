# user-defined distribution: sum of two binomial distribution
dsumofbin <- nimbleFunction(
  run = function(x = double(0), size = double(1), prob = double(1),
                 log = integer(0, default = 0)) {
    returnType(double())
    Prob <- 0.0
    for (i in 0:x) {
      tmp <- dbinom(i, size[1], prob[1],0)*dbinom(x-i, size[2], prob[2],0)
      Prob <- Prob+tmp
    }
    if(log) return(log(Prob))
    else return(Prob)
  })

rsumofbin <- nimbleFunction(
  run = function(n = integer(0), size = double(1), prob = double(1)) {
    returnType(double())
    if(n != 1) print("rsumofbin only allows n = 1; using n = 1.")
    dev <- rbinom(1, size[1], prob[1])+rbinom(1, size[2], prob[2])  
    return(dev)
  })
  
# user-defined distribution: degenerate distribution
ddeg <- nimbleFunction(
  run = function(x = double(0), point = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    Prob <- dconstraint(1,x==point,log)
    return(Prob) 
  })

rdeg <- nimbleFunction(
  run = function(n = integer(0), point = double(0)) {
    returnType(double(0))
    if(n != 1) print("rdeg only allows n = 1; using n = 1.")
    dev <- runif(1,point,point)
    return(dev)
  })

# model
bayes.mod <- nimbleCode( {
  # bayesian model
  I_local_true[1] ~  dpois(mu_local)
  lambda[1] <- 0
  for (t in 2:n_days)
  {
    for (k in 2:t) {
      tmp_matrix[t-1,k-1] <- si_distr[k]*(I_imported_true[t-k+1]+I_local_true[t-k+1])
    }
    lambda[t] <- sum(tmp_matrix[t-1, 1:(t-1)])
    I_local_true[t] ~ dpois(lambda[t]*R_local[t])
  }
  
  for (i in 1:n_days) {
    I_imported_true[i] ~ dpois(mu_imported[i])
    vec1[i,1] <-  I_local_true[i]
    vec1[i,2] <-  I_imported_true[i]
    vec2[i,1] <- 1-alpha_1
    vec2[i,2] <- alpha_0
    I_local_obs[i] ~ dsumofbin(vec1[i,1:2],vec2[i,1:2] )
    I_imported_obs[i]  ~ ddeg(I_imported_true[i]+I_local_true[i]-I_local_obs[i])
  }
  # prior distritbuion 
  alpha_0 ~ dbeta(zeta_0, xi_0)
  alpha_1 ~ dbeta(zeta_1, xi_1)
  mu_local ~ dgamma(1, 1)
  for(i in 1:n_days)
  {
    R_local[i] ~ dgamma(1, 1)
    mu_imported[i] ~ dgamma(1,1)
  }
})
