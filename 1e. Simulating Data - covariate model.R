
# ------------------------ Simulation: Covariate Example -----------------------

library(jagsUI)

## Sampling design mirroring the average samples taken at each pool in the 
## real-world carp example. 

# ----------- Sampling Design ---------------
n.pool = 16
n.year = 4
n.site = 4
n.visit = 4
n.sample = 4
# -------------------------------------------

## Data generating parameters now vary within each pool:

# ------- Data Generating Parameters --------

# Large-scale data:
psi1 <- 0.5

bphi.0 <- 0
bphi.x <- -2

bgamma.0 <- 0
bgamma.x <- 2

x <- runif(n = n.pool, min = -2, max = 2)    
mu.phi <- bphi.0 + bphi.x * x      
mu.gamma <- bgamma.0 + bgamma.x * x
                               
phi <- plogis(mu.phi)
gamma <- plogis(mu.gamma)

## Small-scale data:
omega1 <- 0.5

brho.0 <- 0
brho.w <- -2

blambda.0 <- 0
blambda.w <- 2

w <- matrix(runif(n = n.pool * n.site, min = -2, max = 2), 
            nrow = n.site, ncol = n.pool)
mu.rho <- brho.0 + brho.w * w
mu.lambda <- blambda.0 + blambda.w * w 

rho <- plogis(mu.rho)
lambda <- plogis(mu.lambda)

## Visit-scale data:
theta1 <- 0.5 

beps.0 <- 0
beps.w <- -2

bdelta.0 <- 0
bdelta.w <- 2

mu.eps <- beps.0 + beps.w * w
mu.delta <- bdelta.0 + bdelta.w * w

eps <- plogis(mu.eps)
delta <- plogis(mu.delta)

## Detection Probability:
p <- 0.5
# --------------------------------------------

## Parameters to monitor:

# --------------------------------------------

## 1. Full Model:
## Large-Scale:
psi1_mean <- numeric(sims)
psi1_lcl <- numeric(sims)
psi1_ucl <- numeric(sims)

bphi0_mean <- numeric(sims)
bphi0_lcl <- numeric(sims)
bphi0_ucl <- numeric(sims)

bphix_mean <- numeric(sims)
bphix_lcl <- numeric(sims)
bphix_ucl <- numeric(sims)

bgamma0_mean <- numeric(sims)
bgamma0_lcl <- numeric(sims)
bgamma0_ucl <- numeric(sims)

bgammax_mean <- numeric(sims)
bgammax_lcl <- numeric(sims)
bgammax_ucl <- numeric(sims)

## Small-Scale:
omega1_mean <- numeric(sims)
omega1_lcl <- numeric(sims)
omega1_ucl <- numeric(sims)

brho0_mean <- numeric(sims)
brho0_lcl <- numeric(sims)
brho0_ucl <- numeric(sims)

brhow_mean <- numeric(sims)
brhow_lcl <- numeric(sims)
brhow_ucl <- numeric(sims)

blambda0_mean <- numeric(sims)
blambda0_lcl <- numeric(sims)
blambda0_ucl <- numeric(sims)

blambdaw_mean <- numeric(sims)
blambdaw_lcl <- numeric(sims)
blambdaw_ucl <- numeric(sims)

## Intra-Annual Scale:
theta1_mean <- numeric(sims)
theta1_lcl <- numeric(sims)
theta1_ucl <- numeric(sims)

beps0_mean <- numeric(sims)
beps0_lcl <- numeric(sims)
beps0_ucl <- numeric(sims)

bepsw_mean <- numeric(sims)
bepsw_lcl <- numeric(sims)
bepsw_ucl <- numeric(sims)

bdelta0_mean <- numeric(sims)
bdelta0_lcl <- numeric(sims)
bdelta0_ucl <- numeric(sims)

bdeltaw_mean <- numeric(sims)
bdeltaw_lcl <- numeric(sims)
bdeltaw_ucl <- numeric(sims)

## Detection:
p_mean <- numeric(sims)
p_lcl <- numeric(sims)
p_ucl <- numeric(sims)

## MCMC Convergence Information:
Rhat_sim <- list()


## Looping over each simulation:

sims <- 1000

# ------------ Simulation Study --------------

for(sim in 1:sims){
  
  ## Latent variables (z,s,a):
  
  ## Latent variable z[i,t] -
  ## true presence/absence of pool i during year t:
  
  z <- matrix(0, nrow = n.pool, ncol = n.year)
  
  z[, 1] <- rbinom(n.pool, size = 1, prob = psi1)
  
  for (i in 1:n.pool) {
    for (t in 2:n.year) {
      
      z[i, t] <- rbinom(1, size = 1, prob =
                             ((z[i, t-1] * phi[i]) +
                             ((1 - z[i, t-1]) * gamma[i])))
    }
  }
  
  
  ## Latent variable s[i,j,t] -
  ## true presence/absence of site j in pool i during year t:
  
  s <- array(0, dim = c(n.pool, n.site, n.year))
  
  for (i in 1:n.pool) {
    for (j in 1:n.site) {
      
      s[i, j, 1] <- rbinom(1, size = 1, prob = omega1 * z[i,1])
      
      for (t in 2:n.year) {
        
        s[i, j, t] <- rbinom(1, size = 1, prob =
                                  ((s[i,j,t-1] * rho[j,i]) +
                                  (( 1 - s[i,j,t-1]) * lambda[j,i])) * z[i,t])
      }
    }
  }
  
  ## Latent variable a[j,v,i,t] -
  ## true presence/absence of site j in pool i during visit v within year t:
  
  a <- array(0, dim = c(n.site, n.visit, n.pool, n.year))
  
  for (i in 1:n.pool) {
    for (j in 1:n.site) {
      for (t in 1:n.year) {
        
        a[j,1,i,t] <- rbinom(1, size = 1, prob = theta1 * s[i,j,t])
        
        for(v in 2:n.visit){
          
          a[j,v,i,t] <- rbinom(1, size = 1, prob =
                                 ((a[j,v-1,i,t] * eps[j,i]) +
                                    ((1 - a[j,v-1,i,t]) * delta[j,i])) * s[i,j,t])
        }
      }
    }
  }
  
  
  ## Simulating observations:
  
  y <- array(0, dim = c(n.visit, n.sample, n.site, n.pool, n.year))
  
  for (k in 1:n.sample) {
    for (v in 1:n.visit) {
      for (j in 1:n.site) {
        for (i in 1:n.pool) {
          for (t in 1:n.year) {
            
            y[v,k,j,i,t] <- rbinom(1, 1, prob = p * a[j,v,i,t])
            
          }
        }
      }
    }
  }
  
  ## Gathering the data for the JAGS model: 
  
  sim_jags <- list(n.year = n.year,
                   n.pool = n.pool,
                   n.site = n.site,
                   n.sample = n.sample,
                   n.visit = n.visit,
                   y = y,
                   x = x,
                   w = w)
  
  ## Generating the inits values:
  
  zst <- matrix(1, nrow = n.pool, ncol = n.year) 
  sst <- array(1, dim = c(n.pool, n.site, n.year))   
  ast <- array(1, dim = c(n.site, n.visit, n.pool, n.year))  
  
  inits <- function() {list(z = zst, s = sst, a = ast)}
  
  ## Specify the parameters to report:
  
  params <- c("psi1", 
              "bphi.0", "bphi.x",
              "bgamma.0", "bgamma.x",
              "omega1", 
              "brho.0", "brho.w",
              "blambda.0","blambda.w",
              "theta1", 
              "beps.0", "beps.w",
              "bdelta.0", "bdelta.w",
              "p")  
  
  ## Running the model:
  
  sim_results <- autojags(sim_jags,              
                          inits,                
                          params,                
                          'Simulation_Code/Sim_modelcovar.txt',          
                          n.chains = 3, 
                          n.thin = 10,           
                          n.burnin = 30000,
                          parallel = T,
                          save.all.iter = TRUE)
  
  ## Saving the results:
  
  psi1_mean[sim] <- sim_results$mean$psi1
  psi1_lcl[sim] <- sim_results$q2.5$psi1
  psi1_ucl[sim] <- sim_results$q97.5$psi1
  
  bphi0_mean[sim] <- sim_results$mean$bphi.0
  bphi0_lcl[sim] <- sim_results$q2.5$bphi.0
  bphi0_ucl[sim] <- sim_results$q97.5$bphi.0
  
  bphix_mean[sim] <- sim_results$mean$bphi.x
  bphix_lcl[sim] <- sim_results$q2.5$bphi.x
  bphix_ucl[sim] <- sim_results$q97.5$bphi.x
  
  bgamma0_mean[sim] <- sim_results$mean$bgamma.0
  bgamma0_lcl[sim] <- sim_results$q2.5$bgamma.0
  bgamma0_ucl[sim] <- sim_results$q97.5$bgamma.0
  
  bgammax_mean[sim] <- sim_results$mean$bgamma.x
  bgammax_lcl[sim] <- sim_results$q2.5$bgamma.x
  bgammax_ucl[sim] <- sim_results$q97.5$bgamma.x
  
  omega1_mean[sim] <- sim_results$mean$omega1
  omega1_lcl[sim] <- sim_results$q2.5$omega1
  omega1_ucl[sim] <- sim_results$q97.5$omega1
  
  brho0_mean[sim] <- sim_results$mean$brho.0
  brho0_lcl[sim] <- sim_results$q2.5$brho.0
  brho0_ucl[sim] <- sim_results$q97.5$brho.0
  
  brhow_mean[sim] <- sim_results$mean$brho.w
  brhow_lcl[sim] <- sim_results$q2.5$brho.w
  brhow_ucl[sim] <- sim_results$q97.5$brho.w
  
  blambda0_mean[sim] <- sim_results$mean$blambda.0
  blambda0_lcl[sim] <- sim_results$q2.5$blambda.0
  blambda0_ucl[sim] <- sim_results$q97.5$blambda.0
  
  blambdaw_mean[sim] <- sim_results$mean$blambda.w
  blambdaw_lcl[sim] <- sim_results$q2.5$blambda.w
  blambdaw_ucl[sim] <- sim_results$q97.5$blambda.w
  
  theta1_mean[sim] <- sim_results$mean$theta1
  theta1_lcl[sim] <- sim_results$q2.5$theta1
  theta1_ucl[sim] <- sim_results$q97.5$theta1
  
  beps0_mean[sim] <- sim_results$mean$beps.0
  beps0_lcl[sim] <- sim_results$q2.5$beps.0
  beps0_ucl[sim] <- sim_results$q97.5$beps.0
  
  bepsw_mean[sim] <- sim_results$mean$beps.w
  bepsw_lcl[sim] <- sim_results$q2.5$beps.w
  bepsw_ucl[sim] <- sim_results$q97.5$beps.w
  
  bdelta0_mean[sim] <- sim_results$mean$bdelta.0
  bdelta0_lcl[sim] <- sim_results$q2.5$bdelta.0
  bdelta0_ucl[sim] <- sim_results$q97.5$bdelta.0
  
  bdeltaw_mean[sim] <- sim_results$mean$bdelta.w
  bdeltaw_lcl[sim] <- sim_results$q2.5$bdelta.w
  bdeltaw_ucl[sim] <- sim_results$q97.5$bdelta.w
  
  p_mean[sim] <- sim_results$mean$p
  p_lcl[sim] <- sim_results$q2.5$p
  p_ucl[sim] <- sim_results$q97.5$p 

  Rhat_sim[[sim]] <- sim_results$Rhat
  
  
}


results <- data.frame(psi1_lcl, psi1_mean, psi1_ucl,
                      
                      bphi0_lcl, bphi0_mean, bphi0_ucl,
                      bphix_lcl, bphix_mean, bphix_ucl,
                      
                      bgamma0_lcl, bgamma0_mean, bgamma0_ucl,
                      bgammax_lcl, bgammax_mean, bgammax_ucl,
                      
                      omega1_lcl, omega1_mean, omega1_ucl,
                      
                      brho0_lcl, brho0_mean, brho0_ucl,
                      brhow_lcl, brhow_mean, brhow_ucl,
                      
                      blambda0_lcl, blambda0_mean, blambda0_ucl,
                      blambdaw_lcl, blambdaw_mean, blambdaw_ucl,
                      
                      theta1_lcl, theta1_mean, theta1_ucl,
                      
                      beps0_lcl, beps0_mean, beps0_ucl,
                      bepsw_lcl, bepsw_mean, bepsw_ucl,
                      
                      bdelta0_lcl, bdelta0_mean, bdelta0_ucl,
                      bdeltaw_lcl, bdeltaw_mean, bdeltaw_ucl,
                      
                      p_lcl, p_mean, p_ucl)
