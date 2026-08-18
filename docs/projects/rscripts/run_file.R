# Daniel Wang
# MCMC Algorithm
# 4/22/25

library(mvtnorm)
library(invgamma)

# _____________________________________________________________________________
# PDF of the gamma distribution with new parameterization

log.gamma.pdf <- function(y, v, mu) {
  
  expr <- (1/gamma(v)) * (v/mu)^v * y^(v-1) * exp(-v*y/mu)
  
  return(log(expr))
}
# _____________________________________________________________________________




# _____________________________________________________________________________
# Function to compute log-likelihood of the data (Y~exponential)

log_density_exp <- function( par, prior_par, par_index, y, X) {
  
  n = length(y)
  beta = matrix(par[par_index$beta], ncol=1)
  
  log_fn <- sum(dgamma(c(y), shape = 1, scale = exp(X%*%beta), log = T))
  
  beta.mean = prior_par$prior_mean[par_index$beta]
  beta.sd = diag(prior_par$prior_sd[par_index$beta])
  
  log_pd_beta = dmvnorm(x=c(beta), mean=beta.mean, sigma=beta.sd, log=T)
  
  return(log_fn + log_pd_beta)
}
# _____________________________________________________________________________




# _____________________________________________________________________________
# Function to compute log-likelihood of the data (Y~gamma)

log_density_gamma <- function( par, prior_par, par_index, y, X) {
  
  n = length(y)
  beta = matrix(par[par_index$beta], ncol=1)
  v = par[par_index$v]
  
  log_fn <- sum(log.gamma.pdf(c(y), v, exp(X %*% beta)))
  
  beta.mean = prior_par$prior_mean[par_index$beta]
  beta.sd = diag(prior_par$prior_sd[par_index$beta])
  v_shape = prior_par$v_shape
  v_scale = prior_par$v_scale
  
  log_pd_beta = dmvnorm(x=c(beta), mean=beta.mean, sigma=beta.sd, log=T)
  
  log_pd_v = dgamma(x=v, shape = v_shape, scale = v_scale, log = T)
  
  return(log_fn + log_pd_beta + log_pd_v)
}
# _____________________________________________________________________________




# _____________________________________________________________________________
## MCMC routine using MH algorithm
metro <- function(y, X, log_d, init_par, prior_par, par_index, steps, burnin) {
  
  par = init_par
  n=length(y)
  n_par = length(par)
  chain = matrix( 0, steps, n_par)
  
  # referencing different parameters depending on which distribution the response is fitted under
  group <- par_index
  n_group = length(group)
  
  pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(group[[j]]))
  pscale = rep(.0001, n_group)
  accept = rep(0, n_group)
  
  # Evaluate the log_post of the initial par
  log_density <- log_d
  log_post_prev = log_density(par, prior_par, par_index, y, X)
  if(!is.finite(log_post_prev)){
    print("Infinite log-posterior; choose better initial parameters")
    break
  }
  
  # ___________________________________________________________________________
  
  # Begin Algorithm
  chain[1,] = par
  for(i in 2:steps) {
    for (j in 1:n_group) {
      
      # Propose an update
      ind_j = group[[j]]
      proposal = par
      if (j==1) {
        proposal[ind_j]=rmvnorm(n=1, mean=par[ind_j], sigma=pcov[[j]]*pscale[j])
      } else {
        proposal[ind_j]=rgamma(n=1, shape=par[ind_j]/exp(pcov[[j]]*pscale[j]), scale=exp(pcov[[j]]*pscale[j]))
      }
      
      # Compute the log density for the proposal
      log_post = log_density(proposal, prior_par, par_index, y, X)
      
      # Only propose valid parameters during the burnin period
      if(i < burnin){
        while(!is.finite(log_post)){
          print('bad proposal')
          proposal = par
          if (j==1) {
            proposal[ind_j]=rmvnorm(n=1, mean=par[ind_j], sigma=pcov[[j]]*pscale[j])
          } else {
            proposal[ind_j]=rgamma(n=1, shape=par[ind_j]/exp(pcov[[j]]*pscale[j]), scale=exp(pcov[[j]]*pscale[j]))
          }
          log_post = log_density( proposal, prior_par, par_index, y, X)
        }		
      }
      
      # Evaluate the Metropolis-Hastings ratio
      if( log_post - log_post_prev > log(runif(1,0,1)) ){
        log_post_prev = log_post
        par[ind_j] = proposal[ind_j]
        accept[j] = accept[j] +1
      }
      chain[i,ind_j] = par[ind_j]
      
      # Proposal tuning scheme ------------------------------------------------
      if(i < burnin){
        # During the burnin period, update the proposal covariance in each step 
        # to capture the relationships within the parameters vectors for each 
        # transition.  This helps with mixing.
        if(i == 100)  pscale[j] = 1
        
        if(length(ind_j) > 1){
          if(100 <= i & i <= 2000){  
            temp_chain = chain[1:i,ind_j]
            pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
            
          } else if(2000 < i){  
            temp_chain = chain[(i-2000):i,ind_j]
            pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
          }
          if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )
        }
        
        # Tune the proposal covariance for each transition to achieve 
        # reasonable acceptance ratios.
        if(i %% 30 == 0){ 
          if(i %% 480 == 0){  
            accept[j] = 0  
            
          } else if( accept[j] / (i %% 480) < .4 ){ 
            pscale[j] = (.75^2)*pscale[j] 
            
          } else if( accept[j] / (i %% 480) > .5 ){ 
            pscale[j] = (1.25^2)*pscale[j]
          } 
        }
      }
    }
    # Restart the acceptance ratio at burnin.
    if(i == burnin)  accept[j] = 0
    
    if(i%%100==0)  cat('--->',i,'\n')
  }
  
  print(accept/(steps-burnin))
  cat("pscale = ",pscale,"\n")
  cat("pcov = ","\n")
  print(pcov)
  return(list( chain=chain[burnin:steps,], accept=accept/(steps-burnin),
               pscale=pscale))
}
