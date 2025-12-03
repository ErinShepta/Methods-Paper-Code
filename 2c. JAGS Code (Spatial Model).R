
# --------------------------- Spatial Variable Only ----------------------------

sink("Spatial_model.txt")

cat("model{
  
# ----------------------------------- Priors -----------------------------------
  
  # ------------------- Large-Scale ------------------------
  int.psi1 ~ dnorm(0, 1)             
   
  int.phi ~ dnorm(0, 1)          
  bphi.X ~ dnorm(0, 1)                ## Neighboring pool occupancy(z[n[i],t-1])

  int.gamma ~ dnorm(0, 1)            
  bgamma.X ~ dnorm(0, 1)              ## Neighboring pool occupancy(z[n[i],t-1])
  # ---------------------------------------------------------
  
  # --------------------- Small-Scale -----------------------
  int.omega1 ~ dnorm(0, 1)         
  
  int.rho ~ dnorm(0, 1)             
  brho.X ~ dnorm(0, 1)                ## Neighboring pool occupancy(z[n[i],t-1])
   
  int.lambda ~ dnorm(0, 1)          
  blambda.X ~ dnorm(0, 1)             ## Neighboring pool occupancy(z[n[i],t-1])
  # ----------------------------------------------------------
  
  # ----------------------- Visit Scale ----------------------
  int.theta1 ~ dnorm(0, 1) 
  
  int.eps ~ dnorm(0, 1)             
  beps.X ~ dnorm(0, 1)                ## Neighboring pool occupancy(z[n[i],t-1])
  
  int.delta ~ dnorm(0, 1)           
  bdelta.X ~ dnorm(0, 1)              ## Neighboring pool occupancy(z[n[i],t-1])
  # ----------------------------------------------------------
  
  # -------------------- Detection Scale ---------------------
  int.p ~ dnorm(0, 1)                 ## Reference: Electrofishing (EF)
  bp.X ~ dnorm(0, 1)                  ## Dozer Trawl (DT)
  bp.Y ~ dnorm(0, 1)                  ## Gillnet (GILL)
  # ----------------------------------------------------------
  
# ------------------------------------------------------------------------------

  
# ----------------------- Large-Scale Parameters -------------------------------  
  
  # -------------------- Initial Occupancy ------------------
  for(i in 1:n.pool){                    
      
      logit(psi1[i]) <- (int.psi1)
      z[i, 1] ~ dbern(psi1[i])
  }    
  # ---------------------------------------------------------
  
  # ---------- Persistence and Colonization -----------------
  for(i in 1:n.pool){
    for(t in 2:n.year){
      
      logit(phi[i,t-1]) <- (int.phi + bphi.X * z[n[i], t-1])
      logit(gamma[i,t-1]) <- (int.gamma + bgamma.X * z[n[i], t-1])
      
      z[i, t] ~ dbern(z[i, t - 1] * phi[i, t-1] +
                      (1 - z[i, t - 1]) * gamma[i, t-1])

    }
  }
  
  for(t in 1:n.year){ z[17,t] <- 1 }
    
  for(t in 2:n.year){

    avg_gamma[t] <- sum(gamma[,t-1]) / n.pool
    avg_phi[t] <- sum(phi[,t-1]) / n.pool
    
    psi_eq[t] <- (avg_gamma[t] / (avg_gamma[t] + (1 - avg_phi[t]) + 1.0E-10)) 
    turnover[t] <- sum((1 - z[,t-1]) * z[,t]) / (sum(z[,t]))
    
      }
  
  for(f in 1:n.front){
    for(t in 2:n.year){

     turnover_f[f,t] <- sum((1 - z[f_start[f]:f_end[f],t-1]) * 
                                 z[f_start[f]:f_end[f],t]) / 
                       (sum(z[f_start[f]:f_end[f],t]) + 1.0E-10)
    
    }
  }

  # ----------------------------------------------------------

# ------------------------------------------------------------------------------


# ----------------------- Small-Scale Parameters -------------------------------    

  # -------------------- Initial Occupancy ------------------   
  for(i in 1:n.pool){
    for(j in 1:n.sites[i]){
    
      logit(omega1[j,i]) <- (int.omega1)
      
      s[j,i,1] ~ dbern(omega1[j,i] * z[i,1])
      
    }
  }
  # ---------------------------------------------------------
  
  # ---------- Persistence and Colonization -----------------
  for(i in 1:n.pool){
    for(j in 1:n.sites[i]){
      for(t in 2:n.year){
        
        logit(rho[j,i,t-1]) <- (int.rho + 
                                brho.X * z[n[i],t-1])
        
        logit(lambda[j,i,t-1]) <- (int.lambda + blambda.X * z[n[i],t-1])
        
        s[j,i,t] ~ dbern((s[j,i,t-1] * rho[j,i,t-1] +
                         (1 - s[j,i,t-1]) * lambda[j,i,t-1]) * z[i,t])

      }
    }
  }
  # --------------------------------------------------------
  
  # ----------------- Derived Parameters -------------------
  for(i in 1:n.pool){
    for(t in 2:n.year){

    avg_lambda[i,t] <- sum(lambda[1:n.sites[i],i,t-1]) / n.sites[i]
    avg_rho[i,t]    <- sum(rho[1:n.sites[i],i,t-1]) / n.sites[i]
    
    omega_eq[i,t] <- (avg_lambda[i,t] / 
                      (avg_lambda[i,t] + (1 - avg_rho[i,t]) + 1.0E-10)) 
    
    turnover_s[i,t] <- sum((1 - s[1:n.sites[i],i,t-1]) * s[1:n.sites[i],i,t]) / 
                       (sum(s[1:n.sites[i],i,t]) + 1.0E-10)
    
      }
    }
  # ------------------------------------------------------------
  
# ------------------------------------------------------------------------------


# ----------------------- Small-Scale Parameters ------------------------------- 
  
  # -------------------- Initial Occupancy ------------------ 
    for(i in 1:n.pool){
      for(t in 1:n.year){
        for(j in 1:n.sites[i]){
        
        logit(theta1[i,j,t]) <- (int.theta1)
        a[j,1,i,t] ~ dbern(theta1[i,j,t] * s[j,i,t])
        
         }
       }
     }
  # --------------------------------------------------------

  # ---------- Persistence and Colonization -----------------
    for(i in 1:n.pool){
      for(t in 1:n.year){
        for(j in 1:n.sites[i]){
          for(v in 2:n.visit){
        
          logit(eps[j,v-1,i,t]) <- (int.eps + beps.X * z[n[i],t])
                                    
          logit(delta[j,v-1,i,t]) <- (int.delta + bdelta.X * z[n[i],t])
        
          a[j,v,i,t] ~ dbern((a[j,v-1,i,t] * eps[j,v-1,i,t] +
                             (1 - a[j,v-1,i,t]) * delta[j,v-1,i,t]) * s[j,i,t])
           }
         }
       }
     }
  # --------------------------------------------------------

  # ----------------- Derived Parameters -------------------
    for(i in 1:n.pool){
      for(t in 1:n.year){
        for(v in 2:n.visit){

          turnover_a[v,i,t] <- sum((1 - a[1:n.sites[i],v-1,i,t]) * 
                                        a[1:n.sites[i],v,i,t]) / 
                               (sum(a[1:n.sites[i],v,i,t]) + 1.0E-10)
          }
        }
      }
  # --------------------------------------------------------
  
  
# ----------------------- Detection-Level Parameters ---------------------------

  # -------------------- Detection Probability ------------------
  for(i in 1:N){
    
    logit(p[i]) <- (int.p +
                    bp.X * DT[i] +
                    bp.Y * GILL[i])
    
    y[i] ~ dbern(p[i] * a[Site_Num[i], Visit_Num[i], Pool_Num[i], Year_Num[i]])
    
  }      
  # ------------------------------------------------------------
  
    }  
    " , fill=TRUE)
sink()










# --------------------- With Front and Spatial Variable ------------------------

sink("Categorical_model.txt")

cat("model{
  
# ----------------------------------- Priors -----------------------------------
  
  # ------------------- Large-Scale ------------------------
  int.psi1 ~ dnorm(0, 1)             
   
  int.phi ~ dnorm(0, 1)          
  bphi.X ~ dnorm(0, 1)                ## Neighboring pool occupancy(z[n[i],t-1])
  bphi.Y ~ dnorm(0, 1)                ## IN[i] - Invasion Front
  bphi.Z ~ dnorm(0, 1)                ## PR[i] - Presence Front
  bphi.W ~ dnorm(0, 1)                ## UN[i] - Uninvaded Front

  int.gamma ~ dnorm(0, 1)            
  bgamma.X ~ dnorm(0, 1)              ## Neighboring pool occupancy(z[n[i],t-1])
  bgamma.Y ~ dnorm(0, 1)              ## IN[i] - Invasion Front
  bgamma.Z ~ dnorm(0, 1)              ## PR[i] - Presence Front
  bgamma.W ~ dnorm(0, 1)              ## UN[i] - Uninvaded Front 
  # ---------------------------------------------------------
  
  # --------------------- Small-Scale -----------------------
  int.omega1 ~ dnorm(0, 1)         
  
  int.rho ~ dnorm(0, 1)             
  brho.X ~ dnorm(0, 1)                ## Neighboring pool occupancy(z[n[i],t-1])
  brho.Y ~ dnorm(0, 1)                ## IN[i] - Invasion Front
  brho.Z ~ dnorm(0, 1)                ## PR[i] - Presence Front
  brho.W ~ dnorm(0, 1)                ## UN[i] - Uninvaded Front 
   
  int.lambda ~ dnorm(0, 1)          
  blambda.X ~ dnorm(0, 1)             ## Neighboring pool occupancy(z[n[i],t-1])
  blambda.Y ~ dnorm(0 ,1)             ## IN[i] - Invasion Front 
  blambda.Z ~ dnorm(0 ,1)             ## PR[i] - Presence Front
  blambda.W ~ dnorm(0, 1)             ## UN[i] - Uninvaded Front 
  # ----------------------------------------------------------
  
  # ----------------------- Visit Scale ----------------------
  int.theta1 ~ dnorm(0, 1) 
  
  int.eps ~ dnorm(0, 1)             
  beps.X ~ dnorm(0, 1)              ## Neighboring pool occupancy(z[n[i],t-1])
  beps.Y ~ dnorm(0, 1)              ## IN[i] - Invasion Front 
  beps.Z ~ dnorm(0, 1)              ## PR[i] - Presence Front 
  beps.W ~ dnorm(0, 1)              ## UN[i] - Uninvaded Front 
  
  int.delta ~ dnorm(0, 1)           
  bdelta.X ~ dnorm(0, 1)            ## Neighboring pool occupancy(z[n[i],t-1])
  bdelta.Y ~ dnorm(0, 1)            ## IN[i] - Invasion Front 
  bdelta.Z ~ dnorm(0, 1)            ## PR[i] - Presence Front 
  bdelta.W ~ dnorm(0, 1)            ## UN[i] - Uninvaded Front 
  # ----------------------------------------------------------
  
  # -------------------- Detection Scale ---------------------
  int.p ~ dnorm(0, 1)               ## Reference: Electrofishing (EF)
  bp.X ~ dnorm(0, 1)                ## Dozer Trawl (DT)
  bp.Y ~ dnorm(0, 1)                ## Gillnet (GILL)
  # ----------------------------------------------------------
  
# ------------------------------------------------------------------------------

  
# ----------------------- Large-Scale Parameters -------------------------------  
  
  # -------------------- Initial Occupancy ------------------
  for(i in 1:n.pool){                    
      
      logit(psi1[i]) <- (int.psi1)
      
      z[i, 1] ~ dbern(psi1[i])
  }    
  # ---------------------------------------------------------
  
  # ---------- Persistence and Colonization -----------------
  for(i in 1:n.pool){
    for(t in 2:n.year){
      
      logit(phi[i,t-1]) <- (int.phi + 
                            bphi.X * z[n[i], t-1] + 
                            bphi.Y * IN[i] + 
                            bphi.Z * PR[i] +
                            bphi.W * UN[i])
                            
      logit(gamma[i,t-1]) <- (int.gamma + 
                              bgamma.X * z[n[i], t-1] + 
                              bgamma.Y * IN[i] + 
                              bgamma.Z * PR[i] +
                              bgamma.W * UN[i])
      
      z[i, t] ~ dbern(z[i, t - 1] * phi[i, t-1] +
                        (1 - z[i, t - 1]) * gamma[i, t-1])

    }
  }
  # ----------------------------------------------------------
  
  
    for(t in 2:n.year){

    avg_gamma[t] <- sum(gamma[,t-1]) / n.pool
    avg_phi[t] <- sum(phi[,t-1]) / n.pool
    
    psi_eq[t] <- (avg_gamma[t] / (avg_gamma[t] + (1 - avg_phi[t]) + 1.0E-10)) 
    turnover[t] <- sum((1 - z[,t-1]) * z[,t]) / (sum(z[,t]))
    
      }
  
  for(f in 1:n.front){
    for(t in 2:n.year){

     turnover_f[f,t] <- sum((1 - z[f_start[f]:f_end[f],t-1]) * 
                                 z[f_start[f]:f_end[f],t]) / 
                       (sum(z[f_start[f]:f_end[f],t]) + 1.0E-10)
    
    }
  } 
# ------------------------------------------------------------------------------


# ----------------------- Small-Scale Parameters -------------------------------    

  # -------------------- Initial Occupancy ------------------   
  for(i in 1:n.pool){
    for(j in 1:n.sites[i]){
    
      logit(omega1[j,i]) <- (int.omega1)
      
      s[j,i,1] ~ dbern(omega1[j,i] * z[i,1])
      
    }
  }
  # ---------------------------------------------------------
  
  # ---------- Persistence and Colonization -----------------
  for(i in 1:n.pool){
    for(j in 1:n.sites[i]){
      for(t in 2:n.year){
        
        logit(rho[j,i,t-1]) <- (int.rho + 
                                brho.X * z[n[i],t-1] + 
                                brho.Y * IN[i] + 
                                brho.Z * PR[i] +
                                brho.W * UN[i])
        
        logit(lambda[j,i,t-1]) <- (int.lambda + 
                                   blambda.X * z[n[i],t-1] + 
                                   blambda.Y * IN[i] + 
                                   blambda.Z * PR[i] +
                                   blambda.W * UN[i])
        
        s[j,i,t] ~ dbern((s[j,i,t-1] * rho[j,i,t-1] +
                         (1 - s[j,i,t-1]) * lambda[j,i,t-1]) * z[i,t])

      }
    }
  }
  # --------------------------------------------------------
  
  # ----------------- Derived Parameters -------------------
  for(i in 1:n.pool){
    for(t in 2:n.year){

    avg_lambda[i,t] <- sum(lambda[1:n.sites[i],i,t-1]) / n.sites[i]
    avg_rho[i,t]    <- sum(rho[1:n.sites[i],i,t-1]) / n.sites[i]
    
    omega_eq[i,t] <- (avg_lambda[i,t] / 
                      (avg_lambda[i,t] + (1 - avg_rho[i,t]) + 1.0E-10)) 
    
    turnover_s[i,t] <- sum((1 - s[1:n.sites[i],i,t-1]) * s[1:n.sites[i],i,t]) / 
                       (sum(s[1:n.sites[i],i,t]) + 1.0E-10)
    
      }
    }
  # ------------------------------------------------------------
  
# ------------------------------------------------------------------------------


# ----------------------- Small-Scale Parameters ------------------------------- 
  
  # -------------------- Initial Occupancy ------------------ 
    for(i in 1:n.pool){
      for(t in 1:n.year){
        for(j in 1:n.sites[i]){
        
        logit(theta1[i,j,t]) <- (int.theta1)
        a[j,1,i,t] ~ dbern(theta1[i,j,t] * s[j,i,t])
        
         }
       }
     }
  # --------------------------------------------------------

  # ---------- Persistence and Colonization -----------------
    for(i in 1:n.pool){
      for(t in 1:n.year){
        for(j in 1:n.sites[i]){
          for(v in 2:n.visit){
        
          logit(eps[j,v-1,i,t]) <- (int.eps + 
                                    beps.X * z[n[i],t]  + 
                                    beps.Y * IN[i] + 
                                    beps.Z * PR[i] +
                                    beps.W * UN[i])
                                    
          logit(delta[j,v-1,i,t]) <- (int.delta + 
                                      bdelta.X * z[n[i],t]  + 
                                      bdelta.Y * IN[i] + 
                                      bdelta.Z * PR[i] +
                                      bdelta.W * UN[i])
        
          a[j,v,i,t] ~ dbern((a[j,v-1,i,t] * eps[j,v-1,i,t] +
                             (1 - a[j,v-1,i,t]) * delta[j,v-1,i,t]) * s[j,i,t])
           }
         }
       }
     }
  # --------------------------------------------------------

  # ----------------- Derived Parameters -------------------
    for(i in 1:n.pool){
      for(t in 1:n.year){
        for(v in 2:n.visit){

          turnover_a[v,i,t] <- sum((1 - a[1:n.sites[i],v-1,i,t]) * 
                                        a[1:n.sites[i],v,i,t]) / 
                               (sum(a[1:n.sites[i],v,i,t]) + 1.0E-10)
          }
        }
      }
  
    for(i in 1:n.pool){
      for(t in 1:n.year){
      
       turnover_a[1,i,t] <- 0
       avg_turn_a[i,t] <- sum(turnover_a[,i,t]) / (n.visit - 1)
       
      }
    }  
  # --------------------------------------------------------
  
  
# ----------------------- Detection-Level Parameters ---------------------------

  # -------------------- Detection Probability ------------------
  for(i in 1:N){
    
    logit(p[i]) <- (int.p +
                    bp.X * DT[i] +
                    bp.Y * GILL[i])
    
    y[i] ~ dbern(p[i] * a[Site_Num[i], Visit_Num[i], Pool_Num[i], Year_Num[i]])
    
  }      
  # ------------------------------------------------------------
  
    }  
    " , fill=TRUE)
sink()


