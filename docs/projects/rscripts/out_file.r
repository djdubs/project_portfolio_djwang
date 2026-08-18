# Daniel Wang
# MCMC Computations
# 4/22/25

setwd("projects")
source("rscripts/run_file.r")


# Dataset subsetted for pitchers in 2024 and batting average has been converted to a percentage  
pitching_2015_2024 <- read.csv("data/stats_pitching.csv", header = T)
pitch <- pitching_2015_2024[pitching_2015_2024$year==2024,]
pitch$batting_avg.per <- pitch$batting_avg*100

# Histogram of ERA
png("results/mcmc_p1.png", width = 8, height = 5, units = "in", res=300)
par(mar = c(5.1,4.1,1,1))
hist(pitch$p_era, freq = T, main = "Histogram of ERA", xlab = "ERA")
par(mar = c(5.1,4.1,4.1,2.1))
dev.off()

# _____________________________________________________________________________
## Exponential Model
# Implementing Bayesian computational algorithm on the dataset
set.seed(5)
X = cbind(1, pitch$batting_avg.per, pitch$exit_velocity_avg)
n = nrow(pitch)
p = ncol(X)

y = pitch$p_era

par_index = list( beta=1:p)

init_par = c( rep( 0, p))
prior_par = list( prior_mean=rep(0, p), prior_sd=c(1.23, 0.01, 0.01))
steps = 20000
burnin = 5000

metro.exp.out <- metro(y, X, log_density_exp, init_par, prior_par, par_index, steps, burnin)

# _____________________________________________________________________________
# Estimating the parameters
chain.exp.out <- metro.exp.out$chain
beta.hat.a <- colMeans(chain.exp.out)
var1.hat.a <- var(chain.exp.out[,2])
var2.hat.a <- var(chain.exp.out[,3])


# _____________________________________________________________________________
# Simulation Study with synthetic data
seeds <- seq(46,55)
metro.exp.sim <- list()

p = 3
beta <- beta.hat.a
par_true.a = c( beta)

n = nrow(pitch)

par_index = list( beta=1:p)

init_par = c( rep( 0, p))
prior_par = list( prior_mean=rep(0, p), prior_sd=c(1.23, 0.01, 0.01))
steps = 20000
burnin = 5000

for (i in seeds) {
  set.seed(i)
  X.syn <- cbind(1, rnorm(n, mean = 25, sd = 4), rnorm(n, mean = 88, sd = 1.5))
  
  y.syn <- rgamma(n, shape = 1, scale = exp(X.syn %*% beta))
  
  metro.exp.sim[[i-45]] <- metro(y.syn, X.syn, log_density_exp, init_par, prior_par, par_index, steps, burnin)
}

exp.chains <- lapply(metro.exp.sim, function(x) x$chain)

exp.est <- lapply(exp.chains, colMeans)
exp.sim.totals <- do.call(rbind, exp.chains)


# _____________________________________________________________________________
# MCMC trace plots and histograms
n_post = 15000
index_post = (steps - burnin - n_post + 1):(steps - burnin)
chain.exp.sim = metro.exp.sim[[1]]$chain[index_post,]
param.names <- c("Intercept", "Batting Avg.", "EV Avg.", "V")


par_mean = par_median = upper = lower = rep( NA, ncol(chain.exp.sim))
# pdf("tplots_expmod.pdf")
par(mfrow=c(3, 2))
for(r in 1:ncol(chain.exp.sim)){
  
  plot( NULL, xlab=NA, ylab=NA, xlim=c(1,length(index_post)), 
        main=paste(param.names[r]), ylim=range(chain.exp.sim[,r]) )
  
  lines( chain.exp.sim[,r], type='p', col='black')
  
  par_mean[r] = round( mean(chain.exp.sim[,r]), 4)
  par_median[r] = round( median(chain.exp.sim[,r]), 4)
  upper[r] = round( quantile( chain.exp.sim[,r], prob=.975), 4)
  lower[r] = round( quantile( chain.exp.sim[,r], prob=.025), 4)
  
  hist( chain.exp.sim[,r], ylab=NA, main=NA, freq=F,
        breaks=sqrt(nrow(chain.exp.sim)),
        xlab=paste0('Mean = ',toString(par_mean[r]),
                    ' Median = ',toString(par_median[r])))
  abline( v=upper[r], col='red', lwd=2, lty=2)
  abline( v=lower[r], col='purple', lwd=2, lty=2)
  
  abline( v=par_true.a[r], col='green', lwd=2, lty=2)
}
# dev.off()

# _____________________________________________________________________________
## Gamma Model
# Implementing Bayesian computational algorithm on the dataset
set.seed(20)
X = cbind(1, pitch$batting_avg.per, pitch$exit_velocity_avg)
n = nrow(pitch)
p = ncol(X)

y = pitch$p_era

par_index = list( beta=1:p, v=p+1)

init_par = c( rep( 0, p), exp(0))
prior_par = list( prior_mean=rep(0, p), prior_sd=c(1.23, 0.01, 0.01),
                  v_shape=2, v_scale=2)
steps = 20000
burnin = 5000

metro.gam.out <- metro(y, X, log_density_gamma, init_par, prior_par, par_index, steps, burnin)



# _____________________________________________________________________________
# Estimating the parameters
chain.gamma.out <- metro.gam.out$chain
beta.hat.b <- colMeans(chain.gamma.out)
var1.hat.b <- var(chain.gamma.out[,2])
var2.hat.b <- var(chain.gamma.out[,3])
var3.hat.b <- var(chain.gamma.out[,4])
params_list <- list(beta.hat.b, var1.hat.b, var2.hat.b, var3.hat.b)
save(params_list, file="results/params_list.RData")

# _____________________________________________________________________________
# Simulation Study with synthetic data
seeds <- seq(31, 40)
metro.gam.sim <- list()

p = 3
beta <- beta.hat.b[1:p]
v <- beta.hat.b[p+1]
par_true.b = c(beta, v)

n = nrow(pitch)

par_index = list( beta=1:p, v=p+1)

init_par = c( rep( 0, p), exp(0))
prior_par = list( prior_mean=rep(0, p), prior_sd=c(1.23, 0.01, 0.01),
                  v_shape=2, v_scale=2)
steps = 20000
burnin = 5000

for (i in seeds) {
  set.seed(i)
  
  X.syn <- cbind(1, rnorm(n, mean = 25, sd = 4), rnorm(n, mean = 88, sd = 1.5))
  
  y.syn <- rgamma(n, shape = v, scale = exp(X.syn %*% beta)/v)
  
  metro.gam.sim[[i-30]] <- metro(y.syn, X.syn, log_density_gamma, init_par, prior_par, par_index, steps, burnin)
}

gam.chains <- lapply(metro.gam.sim, function(x) x$chain)

gam.est <- lapply(gam.chains, colMeans)
gam.sim.totals <- do.call(rbind, gam.chains)
sim_res <- list(est = gam.est, totals = gam.sim.totals)
save(sim_res, file="results/sim_res.RData")

# _____________________________________________________________________________
# MCMC trace plots and histograms
n_post = 15000
index_post = (steps - burnin - n_post + 1):(steps - burnin)
chain.gam.sim = metro.gam.sim[[1]]$chain[index_post,]
param.names <- c("Intercept", "Batting Avg.", "EV Avg.", "Shape V")


par_mean = par_median = upper = lower = rep( NA, ncol(chain.gam.sim))
# pdf("tplots_gammod.pdf")
par(mfrow=c(2, 2))
for(r in 1:ncol(chain.gam.sim)){
  
  plot( NULL, xlab=NA, ylab=NA, xlim=c(1,length(index_post)), 
        main=paste(param.names[r]), ylim=range(chain.gam.sim[,r]) )
  
  lines( chain.gam.sim[,r], type='p', col='black')
  
  par_mean[r] = round( mean(chain.gam.sim[,r]), 4)
  par_median[r] = round( median(chain.gam.sim[,r]), 4)
  upper[r] = round( quantile( chain.gam.sim[,r], prob=.975), 4)
  lower[r] = round( quantile( chain.gam.sim[,r], prob=.025), 4)
  
  hist( chain.gam.sim[,r], ylab=NA, main=NA, freq=F,
        breaks=sqrt(nrow(chain.gam.sim)),
        xlab=paste0('Mean = ',toString(par_mean[r]),
                    ' Median = ',toString(par_median[r])))
  abline( v=upper[r], col='red', lwd=2, lty=2)
  abline( v=lower[r], col='purple', lwd=2, lty=2)
  
  abline( v=par_true.b[r], col='green', lwd=2, lty=2)
}
# dev.off()