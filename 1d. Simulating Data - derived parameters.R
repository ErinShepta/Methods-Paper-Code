
# ------------ Simulation: sampling design, data-generating parameters ---------

library(jagsUI)

## Specify the different sampling designs:
## n.pool (N), n.year (T), n.site (J), n.visit (V), n.sample (K)

## Specify the different data-generating parameters:
## initial occupancy - psi1, omega1, theta1
## persistence - phi, rho, eps
## colonization - gamma, lambda, delta
## detection - p

samp.des <- data.frame(n.pool = 16,
                       
                       n.year = 4,
                       
                       n.site = 4,
                       
                       n.visit = 4,
                       
                       n.sample = 4,
                       
                       psi1 = 0.5,
                       omega1 = 0.5,
                       theta1 = 0.5,
                       
                       phi = c(0.5,0.2, 0.8, 0.8),
                       rho = c(0.5,0.2, 0.8, 0.2),
                       eps = c(0.5,0.2, 0.8, 0.2),
                       
                       gamma = c(0.5, 0.8, 0.2, 0.2),
                       lambda = c(0.5, 0.8, 0.2, 0.8),
                       delta = c(0.5, 0.8, 0.2, 0.8),
                       
                       p = 0.5)

## Number of simulations:

sims <- 1000

Rhat_final <- list()

## looping over every sampling design:

for(samp in 1:nrow(samp.des)){
  
  ## Specify the design from the data.frame (samp.des):
  
  n.pool <- samp.des$n.pool[samp]
  n.year <- samp.des$n.year[samp]
  n.site <- samp.des$n.site[samp]
  n.visit <- samp.des$n.visit[samp]
  n.sample <- samp.des$n.sample[samp]
  
  ## Initial Occupancy Probabilities:
  
  psi1 <- samp.des$psi1[samp]
  omega1 <- samp.des$omega1[samp]
  theta1 <- samp.des$theta1[samp]
  
  ## Persistence Probabilities:
  
  phi <- samp.des$phi[samp]
  rho <- samp.des$rho[samp]
  eps <- samp.des$eps[samp]
  
  ## Colonization Probabilities:
  
  gamma <- samp.des$gamma[samp]
  lambda <- samp.des$lambda[samp]
  delta <- samp.des$delta[samp]
  
  ## Detection Probability:
  
  p <- samp.des$p[samp]
  
  
  ## Parameters to be saved for each simulation (results):
  
  psi1_mean <- numeric(sims)
  psi1_lcl <- numeric(sims)
  psi1_ucl <- numeric(sims)
  
  phi_mean <- numeric(sims)
  phi_lcl <- numeric(sims)
  phi_ucl <- numeric(sims)
  
  gamma_mean <- numeric(sims)
  gamma_lcl <- numeric(sims)
  gamma_ucl <- numeric(sims)
  
  omega1_mean <- numeric(sims)
  omega1_lcl <- numeric(sims)
  omega1_ucl <- numeric(sims)
  
  rho_mean <- numeric(sims)
  rho_lcl <- numeric(sims)
  rho_ucl <- numeric(sims)
  
  lambda_mean <- numeric(sims)
  lambda_lcl <- numeric(sims)
  lambda_ucl <- numeric(sims)
  
  theta1_mean <- numeric(sims)
  theta1_lcl <- numeric(sims)
  theta1_ucl <- numeric(sims)
  
  eps_mean <- numeric(sims)
  eps_lcl <- numeric(sims)
  eps_ucl <- numeric(sims)
  
  delta_mean <- numeric(sims)
  delta_lcl <- numeric(sims)
  delta_ucl <- numeric(sims)
  
  p_mean <- numeric(sims)
  p_lcl <- numeric(sims)
  p_ucl <- numeric(sims)
  
  ## Derived Parameters: 
  turn1_mean <- numeric(sims)
  turn1_lcl <- numeric(sims)
  turn1_ucl <- numeric(sims)
  
  turn2_mean <- numeric(sims)
  turn2_lcl <- numeric(sims)
  turn2_ucl <- numeric(sims)
  
  turn3_mean <- numeric(sims)
  turn3_lcl <- numeric(sims)
  turn3_ucl <- numeric(sims)
  
  turn1_s_mean <- numeric(sims)
  turn1_s_lcl <- numeric(sims)
  turn1_s_ucl <- numeric(sims)
  
  turn2_s_mean <- numeric(sims)
  turn2_s_lcl <- numeric(sims)
  turn2_s_ucl <- numeric(sims)
  
  turn3_s_mean <- numeric(sims)
  turn3_s_lcl <- numeric(sims)
  turn3_s_ucl <- numeric(sims)
  
  turn1_a_mean <- numeric(sims)
  turn1_a_lcl <- numeric(sims)
  turn1_a_ucl <- numeric(sims)
  
  turn2_a_mean <- numeric(sims)
  turn2_a_lcl <- numeric(sims)
  turn2_a_ucl <- numeric(sims)
  
  turn3_a_mean <- numeric(sims)
  turn3_a_lcl <- numeric(sims)
  turn3_a_ucl <- numeric(sims)
  
  ## Convergence Diagnostics:
  
  Rhat_sim <- list()
  
  ## Looping over every simulation:
  
  for(sim in 1:sims){
    
    ## Latent variables (z,s,a):
    
    ## Latent variable z[i,t] -
    ## true presence/absence of pool i during year t:
    
    z <- matrix(0, nrow = n.pool, ncol = n.year)
    
    for (i in 1:n.pool) {
      
      z[, 1] <- rbinom(n.pool, size = 1, prob = psi1)
      
      for (t in 2:n.year) {
        
        z[i, t] <- rbinom(1, size = 1, prob =
                            ((z[i, t-1] * phi) +
                               ((1 - z[i, t-1]) * gamma)))
      }
    }
    
    
    ## Latent variable s[i,j,t] -
    ## true presence/absence of site j in pool i during year t:
    
    s <- array(0, dim = c(n.site, n.pool, n.year))
    
    for (i in 1:n.pool) {
      for (j in 1:n.site) {
        
        s[j, i, 1] <- rbinom(1, size = 1, prob = omega1 * z[i,1])
        
        for (t in 2:n.year) {
          
          s[j, i, t] <- rbinom(1, size = 1, prob =
                                 ((s[j,i,t-1] * rho) +
                                    (( 1 - s[j,i,t-1]) * lambda)) * z[i,t])
        }
      }
    }
    
    ## Latent variable a[j,v,i,t] -
    ## true presence/absence of site j in pool i during visit v within year t:
    
    a <- array(0, dim = c(n.site, n.visit, n.pool, n.year))
    
    for (i in 1:n.pool) {
      for (j in 1:n.site) {
        for (t in 1:n.year) {
          
          a[j,1,i,t] <- rbinom(1, size = 1, prob = theta1 * s[j,i,t])
          
          for(v in 2:n.visit){
            
            a[j,v,i,t] <- rbinom(1, size = 1, prob =
                                   ((a[j,v-1,i,t] * eps) +
                                      ((1 - a[j,v-1,i,t]) * delta)) * s[j, i, t])
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
                     y = y)
    
    
    ## Generating the inits values:
    
    zst <- matrix(1, nrow = n.pool, ncol = n.year) 
    sst <- array(1, dim = c(n.site, n.pool, n.year))   
    ast <- array(1, dim = c(n.site, n.visit, n.pool, n.year))  
    
    inits <- function() {list(z = zst,
                              s = sst,
                              a = ast)}
    
    
    ## Specify the parameters to report:
    
    params <- c("psi1",
                "phi",
                "gamma",
                "omega1",
                "rho",
                "lambda",
                "theta1",
                "eps",
                "delta",
                "p",
                "turnover",
                "turnover_s",
                "turnover_a")                   
    
    ## Running the model:
    
    sim_results <- autojags(sim_jags,              
                            inits,                
                            params,                
                            'Simulation_Code/Sim_modelder.txt',          
                            n.chains = 3, 
                            n.thin = 2,           
                            n.burnin = 30000,
                            parallel = T,
                            save.all.iter = TRUE)
    
    ## Saving the results:
    
    psi1_mean[sim] <- sim_results$mean$psi1
    psi1_lcl[sim] <- sim_results$q2.5$psi1
    psi1_ucl[sim] <- sim_results$q97.5$psi1
    
    phi_mean[sim] <- sim_results$mean$phi
    phi_lcl[sim] <- sim_results$q2.5$phi
    phi_ucl[sim] <- sim_results$q97.5$phi
    
    gamma_mean[sim] <- sim_results$mean$gamma
    gamma_lcl[sim] <- sim_results$q2.5$gamma
    gamma_ucl[sim] <- sim_results$q97.5$gamma
    
    omega1_mean[sim] <- sim_results$mean$omega1
    omega1_lcl[sim] <- sim_results$q2.5$omega1
    omega1_ucl[sim] <- sim_results$q97.5$omega1
    
    rho_mean[sim] <- sim_results$mean$rho
    rho_lcl[sim] <- sim_results$q2.5$rho
    rho_ucl[sim] <- sim_results$q97.5$rho
    
    lambda_mean[sim] <- sim_results$mean$lambda
    lambda_lcl[sim] <- sim_results$q2.5$lambda
    lambda_ucl[sim] <- sim_results$q97.5$lambda
    
    theta1_mean[sim] <- sim_results$mean$theta1
    theta1_lcl[sim] <- sim_results$q2.5$theta1
    theta1_ucl[sim] <- sim_results$q97.5$theta1
    
    eps_mean[sim] <- sim_results$mean$eps
    eps_lcl[sim] <- sim_results$q2.5$eps
    eps_ucl[sim] <- sim_results$q97.5$eps
    
    delta_mean[sim] <- sim_results$mean$delta
    delta_lcl[sim] <- sim_results$q2.5$delta
    delta_ucl[sim] <- sim_results$q97.5$delta
    
    p_mean[sim] <- sim_results$mean$p
    p_lcl[sim] <- sim_results$q2.5$p
    p_ucl[sim] <- sim_results$q97.5$p 
    
    psi_eq_mean[sim] <- sim_results$mean$psi_eq
    psi_eq_lcl[sim] <- sim_results$q2.5$psi_eq
    psi_eq_ucl[sim] <- sim_results$q97.5$psi_eq
    
    omega_eq_mean[sim] <- sim_results$mean$omega_eq
    omega_eq_lcl[sim] <- sim_results$q2.5$omega_eq
    omega_eq_ucl[sim] <- sim_results$q97.5$omega_eq
    
    turn1_mean[sim] <- sim_results$mean$turnover[2]
    turn1_lcl[sim] <- sim_results$q2.5$turnover[2]
    turn1_ucl[sim] <- sim_results$q97.5$turnover[2]
    turn2_mean[sim] <- sim_results$mean$turnover[3]
    turn2_lcl[sim] <- sim_results$q2.5$turnover[3]
    turn2_ucl[sim] <- sim_results$q97.5$turnover[3]
    turn3_mean[sim] <- sim_results$mean$turnover[4]
    turn3_lcl[sim] <- sim_results$q2.5$turnover[4]
    turn3_ucl[sim] <- sim_results$q97.5$turnover[4]
    
    turn1_s_mean[sim] <- sim_results$mean$turnover_s[2]
    turn1_s_lcl[sim] <- sim_results$q2.5$turnover_s[2]
    turn1_s_ucl[sim] <- sim_results$q97.5$turnover_s[2]
    turn2_s_mean[sim] <- sim_results$mean$turnover_s[3]
    turn2_s_lcl[sim] <- sim_results$q2.5$turnover_s[3]
    turn2_s_ucl[sim] <- sim_results$q97.5$turnover_s[3]
    turn3_s_mean[sim] <- sim_results$mean$turnover_s[4]
    turn3_s_lcl[sim] <- sim_results$q2.5$turnover_s[4]
    turn3_s_ucl[sim] <- sim_results$q97.5$turnover_s[4]
    
    turn1_a_mean[sim] <- sim_results$mean$turnover_a[2]
    turn1_a_lcl[sim] <- sim_results$q2.5$turnover_a[2]
    turn1_a_ucl[sim] <- sim_results$q97.5$turnover_a[2]
    turn2_a_mean[sim] <- sim_results$mean$turnover_a[3]
    turn2_a_lcl[sim] <- sim_results$q2.5$turnover_a[3]
    turn2_a_ucl[sim] <- sim_results$q97.5$turnover_a[3]
    turn3_a_mean[sim] <- sim_results$mean$turnover_a[4]
    turn3_a_lcl[sim] <- sim_results$q2.5$turnover_a[4]
    turn3_a_ucl[sim] <- sim_results$q97.5$turnover_a[4]
    
    Rhat_sim[[sim]] <- sim_results$Rhat
    
  }
  
  ## gathering results into a data.frame:
  
  results <- data.frame(psi1_lcl, psi1_mean, psi1_ucl,
                        phi_lcl, phi_mean, phi_ucl,
                        gamma_lcl,gamma_mean, gamma_ucl,
                        omega1_lcl, omega1_mean, omega1_ucl,
                        rho_lcl, rho_mean, rho_ucl,
                        lambda_lcl, lambda_mean, lambda_ucl,
                        theta1_lcl, theta1_mean, theta1_ucl,
                        eps_lcl, eps_mean, eps_ucl,
                        delta_lcl, delta_mean, delta_ucl,
                        p_lcl, p_mean, p_ucl,
                        
                        psi_eq_lcl, psi_eq_mean, psi_eq_ucl,
                        omega_eq_lcl, omega_eq_mean, omega_eq_ucl,
                        
                        turn1_lcl, turn1_mean, turn1_ucl,
                        turn2_lcl, turn2_mean, turn2_ucl,
                        turn3_lcl, turn3_mean, turn3_ucl,
                        
                        turn1_s_lcl, turn1_s_mean, turn1_s_ucl,
                        turn2_s_lcl, turn2_s_mean, turn2_s_ucl,
                        turn3_s_lcl, turn3_s_mean, turn3_s_ucl,
                        
                        turn1_a_lcl, turn1_a_mean, turn1_a_ucl,
                        turn2_a_lcl, turn2_a_mean, turn2_a_ucl,
                        turn3_a_lcl, turn3_a_mean, turn3_a_ucl)
  
  ## Adding a column to specify the sampling design:
  
  results$samp <- samp
  
  ## If results do not exist, results = RESULTS
  ## else add new 'results' onto RESULTS data.frame
  ## ensure ALL results will be saved into one final data.frame:
  
  if(!exists('RESULTS')){
    RESULTS <- results
  }else{
    RESULTS <- rbind(RESULTS,results)
  }
  
  ## Convergence diagnostics:
  
  Rhat_final[[samp]] <- Rhat_sim
  
}





