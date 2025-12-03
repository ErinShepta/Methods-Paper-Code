# --------------------- Sampling Design Model (Intercept Only) -----------------

## Code for a dynamic multi-scale occupancy model looking at persistence and 
## colonization probabilities across 3 scales: annual large-scale, annual 
## small-scale, intra-annual small-scale. This is an intercept only model used
## in simulation studies to test the impact of sampling design and data 
## generating parameters.

sink("Simulation_Code/Sim_model.txt")

cat("model{

# ------------------------------------------------------------------------------

 # -------- Priors ---------

   # 1. Annual Unit[i,t] Priors:
      psi1 ~ dbeta(1, 1)        
      phi ~ dbeta(1, 1)         
      gamma ~ dbeta(1, 1)      
   # 2. Annual Subunit[i,j,t] Priors:
      omega1 ~ dbeta(1, 1)       
      rho ~ dbeta(1, 1)          
      lambda ~ dbeta(1, 1)        
   # 3. Intra-Annual SUbunit [j,v,i,t] Priors:
      theta1 ~ dbeta(1, 1)       
      eps ~ dbeta(1, 1)           
      delta ~ dbeta(1, 1)       
   # 4. Detection [k,v,j,i,t] Priors
      p ~ dbeta(1, 1)            
  # --------------------------


# ----------------------- 1. Annual Large-Scale Model --------------------------

  # Initial Occupancy: 
  for(i in 1:n.pool) { z[i,1] ~ dbern(psi1) }
  
  ## Colonization and Persistence:
  for(i in 1:n.pool){
   for(t in 2:n.year) {
      z[i,t] ~ dbern((z[i, t-1] * phi) + ((1 - z[i, t-1]) * gamma))
      
    }
  }

 # ---------------------- 2. Annual Small-Scale Model --------------------------
 
  # Initial Occupancy:
    for(i in 1:n.pool){
     for(j in 1:n.site){
       s[i,j,1] ~ dbern(omega1 * z[i,1]) 
    }
 }

  # Colonization and Persistence:
    for(i in 1:n.pool){
     for(j in 1:n.site){
       for(t in 2:n.year){
       s[i,j,t] ~ dbern((((s[i,j,t-1] * rho) + 
                        ((1 - s[i,j,t-1]) * lambda))) *
                        z[i,t])
                        
       
        }
      }
    }

 # -------------------- 3. Intra-Annual Small-Scale Model ----------------------

  # Initial Occupancy:
    for(t in 1:n.year) {
      for(i in 1:n.pool) {
       for(j in 1:n.site) {
      a[j,1,i,t] ~ dbern(theta1 * s[i,j,t])
     }
   }
 }
 
 # Colonization and Persistence:
 for(t in 1:n.year) {
  for(i in 1:n.pool) {
    for(j in 1:n.site) {
      for(v in 2:n.visit) {
        a[j,v,i,t] ~ dbern(((a[j,v-1,i,t] * eps) + 
                            ((1 - a[j,v-1,i,t]) * delta))  * 
                            s[i,j,t])
        }
      }
    }
  }

 # -------------------------- 4. Detection Probability -------------------------

 for(t in 1:n.year) {
   for(i in 1:n.pool) {
     for(j in 1:n.site) {
        for(v in 1:n.visit) {
          for(k in 1:n.sample) {
          y[v,k,j,i,t] ~ dbern(p * a[j,v,i,t])
           }
         }
       }
     }
   }

  } ", fill=TRUE)
sink()

# ------------------------------------------------------------------------------


# --------------------- Sampling Design Model (Derived Values) -----------------

## Code for a dynamic multi-scale occupancy model looking at persistence and 
## colonization probabilities across 3 scales: annual large-scale, annual 
## small-scale, intra-annual small-scale. This is an intercept only model used
## in simulation studies. Here we are deriving occupancy parameters (turnover,
## turnover_s, turnover_a) for three scales. 

sink("Simulation_Code/Sim_modelder.txt")

cat("model{

# ------------------------------------------------------------------------------

 # -------- Priors ---------

   # 1. Annual Unit[i,t] Priors:
      psi1 ~ dbeta(1, 1)        
      phi ~ dbeta(1, 1)         
      gamma ~ dbeta(1, 1)      
   # 2. Annual Subunit[i,j,t] Priors:
      omega1 ~ dbeta(1, 1)       
      rho ~ dbeta(1, 1)          
      lambda ~ dbeta(1, 1)        
   # 3. Intra-Annual SUbunit [j,v,i,t] Priors:
      theta1 ~ dbeta(1, 1)       
      eps ~ dbeta(1, 1)           
      delta ~ dbeta(1, 1)       
   # 4. Detection [k,v,j,i,t] Priors
      p ~ dbeta(1, 1)            
  # --------------------------


# ----------------------- 1. Annual Large-Scale Model --------------------------

  # Initial Occupancy: 
  for(i in 1:n.pool) { z[i,1] ~ dbern(psi1) }
  
  ## Colonization and Persistence:
  for(i in 1:n.pool){
   for(t in 2:n.year) {
      z[i,t] ~ dbern((z[i, t-1] * phi) + ((1 - z[i, t-1]) * gamma))
      
    }
  }
  
  for(t in 2:n.year){

    turnover[t] <- sum((1 - z[,t-1]) * z[,t]) / sum(z[,t])
  
  }



 # ---------------------- 2. Annual Small-Scale Model --------------------------
 
  # Initial Occupancy:
    for(i in 1:n.pool){
     for(j in 1:n.site){
       s[j,i,1] ~ dbern(omega1 * z[i,1]) 
    }
 }

  # Colonization and Persistence:
    for(i in 1:n.pool){
     for(j in 1:n.site){
       for(t in 2:n.year){
       s[j,i,t] ~ dbern((((s[j,i,t-1] * rho) + 
                        ((1 - s[j,i,t-1]) * lambda))) *
                        z[i,t])
                        
       
        }
      }
    }


    for(t in 2:n.year){

      turnover_s[t] <- sum((1 - s[,,t-1]) * s[,,t]) / sum(s[,,t])

      }

 # -------------------- 3. Intra-Annual Small-Scale Model ----------------------

  # Initial Occupancy:
    for(t in 1:n.year) {
      for(i in 1:n.pool) {
       for(j in 1:n.site) {
      a[j,1,i,t] ~ dbern(theta1 * s[j,i,t])
     }
   }
 }
 
 # Colonization and Persistence:
 for(t in 1:n.year) {
  for(i in 1:n.pool) {
    for(j in 1:n.site) {
      for(v in 2:n.visit) {
        a[j,v,i,t] ~ dbern(((a[j,v-1,i,t] * eps) + 
                            ((1 - a[j,v-1,i,t]) * delta))  * 
                            s[j,i,t])
        }
      }
    }
  }
 
for(t in 1:n.year){
  
  turnover_a[t] <- sum( (1 - a[,2:n.visit,,t]) * a[,2:n.visit,,t] ) /
                    sum( a[,2:n.visit,,t] )
  
  }
      


 # -------------------------- 4. Detection Probability -------------------------

 for(t in 1:n.year) {
   for(i in 1:n.pool) {
     for(j in 1:n.site) {
        for(v in 1:n.visit) {
          for(k in 1:n.sample) {
          y[v,k,j,i,t] ~ dbern(p * a[j,v,i,t])
           }
         }
       }
     }
   }

  } ", fill=TRUE)
sink()

# ------------------------------------------------------------------------------


# ---------------- Full Model (Heterogeneous Sampling Design) ------------------

## Code for a dynamic multi-scale occupancy model looking at persistence and 
## colonization probabilities across 3 scales: annual large-scale, annual 
## small-scale, intra-annual small-scale. Data generating scenario assumes that 
## sampling during each visit is uneven.  

sink("Simulation_Code/Sim_model5.txt")

cat("model{

# ------------------------------------------------------------------------------

 # -------- Priors ---------

   # 1. Annual Unit[i,t] Priors:
      psi1 ~ dbeta(1, 1)        
      phi ~ dbeta(1, 1)         
      gamma ~ dbeta(1, 1)      
   # 2. Annual Subunit[i,j,t] Priors:
      omega1 ~ dbeta(1, 1)     
      rho ~ dbeta(1, 1)          
      lambda ~ dbeta(1, 1) 
   # 3. Intra-Annual SUbunit [j,v,i,t] Priors:
      theta1 ~ dbeta(1, 1) 
      eps ~ dbeta(1, 1)           
      delta ~ dbeta(1, 1) 
   # 4. Detection [k,v,j,i,t] Priors
      p ~ dbeta(1, 1)            
  # --------------------------


# ----------------------- 1. Annual Large-Scale Model --------------------------

  # Initial Occupancy: 
  for(i in 1:n.pool) { z[i,1] ~ dbern(psi1) } 
  
  ## Colonization and Persistence:
  for(i in 1:n.pool){
   for(t in 2:n.year) {
      z[i,t] ~ dbern((z[i, t-1] * phi) + ((1 - z[i, t-1]) * gamma))

    }
  }
 # ---------------------- 2. Annual Small-Scale Model --------------------------
 
  # Initial Occupancy:
    for(i in 1:n.pool){
     for(j in 1:n.site){
     
       s[j,i,1] ~ dbern(omega1 * z[i,1]) 
       
    }
 }

  # Colonization and Persistence:
    for(i in 1:n.pool){
     for(j in 1:n.site){
       for(t in 2:n.year){
       s[j,i,t] ~ dbern((((s[j,i,t-1] * rho) + 
                        ((1 - s[j,i,t-1]) * lambda))) *
                        z[i,t])
      
      }
    }
  }
 # -------------------- 3. Intra-Annual Small-Scale Model ----------------------

  # Initial Occupancy:
    for(t in 1:n.year) {
      for(i in 1:n.pool) {
       for(j in 1:n.site) {
      a[j,1,i,t] ~ dbern(theta1 * s[j,i,t])
     }
   }
 }
 
 # Colonization and Persistence:
 for(t in 1:n.year) {
  for(i in 1:n.pool) {
    for(j in 1:n.site) {
      for(v in 2:n.visit) {
        a[j,v,i,t] ~ dbern(((a[j,v-1,i,t] * eps) + 
                            ((1 - a[j,v-1,i,t]) * delta))  * 
                            s[j,i,t])
        }
      }
    }
  }

 # -------------------------- 4. Detection Probability -------------------------

 for(t in 1:n.year) {
   for(i in 1:n.pool) {
     for(j in 1:n.site) {
        for(v in 1:n.visit) {
          for(k in 1:n.sample) {
          y[v,k,j,i,t] ~ dbern(p * a[j,v,i,t])
           }
         }
       }
     }
   }

  } ", fill=TRUE)
sink()

# ------------------------------------------------------------------------------


# -------------------------- Full Model with Covariates ------------------------

## Code for a dynamic multi-scale occupancy model looking at persistence and 
## colonization probabilities across 3 scales: annual large-scale, annual 
## small-scale, intra-annual small-scale. This is an intercept only model used
## in simulation studies to test the impact of sampling design and data 
## generating parameters.

sink("Simulation_Code/Sim_modelcovar.txt")

cat("model{

# ------------------------------------------------------------------------------

 # -------- Priors ---------

   # 1. Annual Unit[i,t] Priors:
      psi1 ~ dbeta(1, 1)    
      
      bphi.0 ~ dnorm(0,0.1)
      bphi.x ~ dnorm(0,0.1)
      bgamma.0 ~ dnorm(0,0.1)
      bgamma.x ~ dnorm(0,0.1)
      
   # 2. Annual Subunit[i,j,t] Priors:
      omega1 ~ dbeta(1, 1)   
      
      brho.0 ~ dnorm(0,0.1)
      brho.w ~ dnorm(0,0.1)
      blambda.0 ~ dnorm(0,0.1)
      blambda.w ~ dnorm(0,0.1)
      
   # 3. Intra-Annual SUbunit [j,v,i,t] Priors:
      theta1 ~ dbeta(1, 1)  
      
      beps.0 ~ dnorm(0,0.1)
      beps.w ~ dnorm(0,0.1)
      bdelta.0 ~ dnorm(0,0.1)
      bdelta.w ~ dnorm(0,0.1)
      
   # 4. Detection [k,v,j,i,t] Priors
      p ~ dbeta(1, 1)            
  # --------------------------


# ----------------------- 1. Annual Large-Scale Model --------------------------

  # Initial Occupancy: 
  for(i in 1:n.pool) { z[i,1] ~ dbern(psi1) }
  
  ## Colonization and Persistence:
  
  for(i in 1:n.pool){
  
      logit(phi[i]) <- (bphi.0 + bphi.x * x[i])
      logit(gamma[i]) <- (bgamma.0 + bgamma.x * x[i]) 
      
      }
  
  
  for(i in 1:n.pool){
   for(t in 2:n.year) {
   
      z[i,t] ~ dbern((z[i, t - 1] * phi[i]) + 
                    ((1 - z[i, t - 1]) * gamma[i]))
      
    }
  }

 # ---------------------- 2. Annual Small-Scale Model --------------------------
 
  # Initial Occupancy:
    for(i in 1:n.pool){
     for(j in 1:n.site){
       s[i,j,1] ~ dbern(omega1 * z[i,1]) 
    }
 }

  # Colonization and Persistence:
  
   for(i in 1:n.pool){
     for(j in 1:n.site){
       
       logit(rho[j,i]) <- (brho.0 + brho.w * w[j,i])
       logit(lambda[j,i]) <- (blambda.0 + blambda.w * w[j,i])
       
        }
      }
      
    
    for(i in 1:n.pool){
     for(j in 1:n.site){
       for(t in 2:n.year){
       
       s[i,j,t] ~ dbern((((s[i,j,t-1] * rho[j,i]) + 
                        ((1 - s[i,j,t-1]) * lambda[j,i]))) *
                        z[i,t])
                        
       
        }
      }
    }

 # -------------------- 3. Intra-Annual Small-Scale Model ----------------------

  # Initial Occupancy:
    for(t in 1:n.year) {
      for(i in 1:n.pool) {
       for(j in 1:n.site) {
      a[j,1,i,t] ~ dbern(theta1 * s[i,j,t])
     }
   }
 }
 
 # Colonization and Persistence:
 

  for(i in 1:n.pool) {
    for(j in 1:n.site) {
      
       logit(eps[j,i]) <- (beps.0 + beps.w * w[j,i])
       logit(delta[j,i]) <- (bdelta.0 + bdelta.w * w[j,i])

        }
      }
    
  
 for(t in 1:n.year) {
  for(i in 1:n.pool) {
    for(j in 1:n.site) {
      for(v in 2:n.visit) {
      
        a[j,v,i,t] ~ dbern(((a[j,v-1,i,t] * eps[j,i]) + 
                            ((1 - a[j,v-1,i,t]) * delta[j,i]))  * 
                            s[i,j,t])
        }
      }
    }
  }

 # -------------------------- 4. Detection Probability -------------------------

 for(t in 1:n.year) {
   for(i in 1:n.pool) {
     for(j in 1:n.site) {
        for(v in 1:n.visit) {
          for(k in 1:n.sample) {
          y[v,k,j,i,t] ~ dbern(p * a[j,v,i,t])
           }
         }
       }
     }
   }

  } ", fill=TRUE)
sink()

# ------------------------------------------------------------------------------