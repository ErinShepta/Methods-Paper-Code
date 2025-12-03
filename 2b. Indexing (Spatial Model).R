
# ------------------------------- Simple Indexing ------------------------------

n.pool <- max(carp$Pool_Num) 
n.year <- max(carp$Year_Num) 
n.visit <- max(carp$Visit_Num)
n.sites <- numeric(n.pool)


for(i in 1:n.pool){
  
  n.sites[i] <- nrow(carp %>% 
                     filter(Pool_Num == i) %>% 
                       select(Site_Num) %>%
                       distinct(Site_Num)) }

n.front <- 4
f_start <- c(1, 8, 10, 15)
f_end <- c(7, 9, 14, 16)
f_pool <- c(7, 2, 5, 2)


data <- list(n.pool = n.pool,
             n.sites = n.sites,
             n.visit = n.visit,
             n.year = n.year,
             n.front = n.front,
             N = nrow(carp),
             Visit_Num = carp$Visit_Num,
             Site_Num = carp$Site_Num,
             Pool_Num = carp$Pool_Num,
             Year_Num = carp$Year_Num,
             f_start = f_start,
             f_end = f_end,
             f_pool = f_pool,
             DT = carp$DT,
             GILL = carp$GILL,
             n = n,
             INV = INV,
             PR = PR,
             UN = UN,
             y = carp$Silver)

str(data)




z_init <- matrix(nrow = data$n.pool, ncol = data$n.year)


for(i in 1:data$n.pool){
  for(j in 1:data$n.year){
    
    z_init[i, j] <- 1
    
  }
}


s_init <- array(dim = c(max(data$Site_Num),
                        data$n.pool,
                        data$n.year))

for(i in 1:data$n.pool){
  for(j in 1:data$n.sites[i]){
    for(k in 1:data$n.year){
      
      s_init[j,i,k] <- 1
      
    }
  }
}

a_init <- array(dim = c(max(data$Site_Num),
                        data$n.visit,
                        data$n.pool,
                        data$n.year))

for(i in 1:data$n.pool){
  for(t in 1:data$n.year){
    for(j in 1:data$n.sites[i]){
      for(v in 1:data$n.visit){
        
        a_init[j,v,i,t] <- 1
        
      }
    }
  }
}

inits <- function(){
  list(
    z = z_init,
    s = s_init,
    a = a_init
  )
}

