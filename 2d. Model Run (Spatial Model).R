
library(jagsUI)

# --------------------------- Spatial Model Run --------------------------------

params_spat <- c('int.psi1', 'int.phi', 'int.gamma',
                 'int.omega1', 'int.rho', 'int.lambda',
                 'int.theta1', 'int.eps', 'int.delta', 'int.p',
               
                 'bphi.X', 'bgamma.X',
                 'brho.X', 'blambda.X',
                 'beps.X', 'bdelta.X',
                 'bp.X', 'bp.Y',  
               
                 'psi_eq', 'turnover', 'turnover_f',
                 'avg_lambda','avg_rho', 'turnover_s', 
                 'turnover_a')

set.seed(123)

Model_spat <- jags(data = data,
                   inits = inits,
                   params_spat,
                   'Spatial_model.txt',
                   n.chains = 3,
                   n.iter = 2500000,
                   n.burnin = 30000,
                   n.adapt = 1000,
                   n.thin = 30,
                   parallel = T)

Model_spat$Rhat

# Model_spat_update <- update(Model_spat,
#                             parameters.to.save = params_spat,
#                             n.iter = 500000,
#                             n.thin = 30)
# 
# Model_spat_update$Rhat


# saveRDS(Model_spat, file = "Carp_Occupancy/Updated_SpatialResults.rds")

# ------------------------------------------------------------------------------



# --------------------------- Categories and Spatial ---------------------------

params_cats <- c('int.psi1', 'int.phi', 'int.gamma',
                 'int.omega1', 'int.rho', 'int.lambda',
                 'int.theta1', 'int.eps', 'int.delta', 'int.p',
                 
                 'bphi.X', 'bphi.Z','bphi.Y','bphi.W',
                 
                 'bgamma.X','bgamma.Z','bgamma.Y','bgamma.W',
                 
                 'brho.X','brho.Z','brho.Y','brho.W',
                 
                 'blambda.X','blambda.Z','blambda.Y','blambda.W',
                 
                 'beps.X','beps.Z','beps.Y','beps.W',
                 
                 'bdelta.X','bdelta.Z','bdelta.Y','bdelta.W',
                 
                 'bp.X', 'bp.Y',  
                 
                 'psi_eq', 'turnover', 'turnover_f',
                 'omega_eq', 'turnover_s', 
                 'avg_turn_a')

set.seed(123)

Model_cats <- jags(data = data,
                   inits = inits,
                   params_cats,
                   'Spatial_model.txt',
                   n.chains = 3,
                   n.iter = 2000000,
                   n.burnin = 10000,
                   n.adapt = 1000,
                   n.thin = 30,
                   parallel = T)

Model_cats$Rhat

Model_cats_update <- update(Model_cats,
                            parameters.to.save = params_cats,
                            n.iter = 500000)

# saveRDS(Model_update, file = "Carp_Occupancy/Updated_Results.rds")
    