library(dplyr)
library(tidyr)
library(lubridate)

# --------------------------- 1. Carp Observations (y) -------------------------

## Importing and formatting data for carp observations/detection covariates 

carp <- read.csv("Carp_Occupancy/adult_det.csv", 
                 header = T, 
                 na.strings = "") %>%
                 filter(is.na(Silver) == F)

## Main 3 (used for analysis)
## DT (dozer trawl)
## Gillnet
## EF (electrofishing)

## Reference - EF
carp$DT <- ifelse(carp$Gear == "DT", 1, 0)
carp$GILL <- ifelse(carp$Gear == "Gillnet", 1, 0)

## Identity of the downstream pool (place-holder for pool 1: '17'): 
n <- c(17,1,2,1,4,5,6,7,8,9,10,11,12,13,12,15)

## Categorical Fronts (Establishment = reference):
INV <- c(0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0)
PR <- c(0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0)
UN <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1)


