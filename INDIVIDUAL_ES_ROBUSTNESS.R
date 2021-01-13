######## Single ES Robustness Calc. #####
## Written by: Aislyn Keyes
## Date last edited: 4-23-2020

# load libraries (install if necessary)
library(dplyr)
library(tidyr)
library(data.table)
library(tidyverse)
library(stringr)

## source functions for analysis
source("NEWrobustness_functions.R")
## disable warnings 
options(warn=-1)



# load data
# load all three salt marshes for one sequence
# Will run the code for all systems and save it as one dataframe for that sequence
# to run on a different sequence, adjust the file names
EPB.dat <- read.csv("RobustEPB_DirectBiomassSMALLLARGEOutput_ES.csv")
EPB.dat <- EPB.dat[,-1] # get rid of extra column so it doesn't throw off the rest of the code

CSM.dat <- read.csv("RobustCSM_DirectBiomassSMALLLARGEOutput_ES.csv")
CSM.dat <- CSM.dat[,-1]

BSQ.dat <- read.csv("RobustBSQ_DirectBiomassSMALLLARGEOutput_ES.csv")
BSQ.dat <- BSQ.dat[,-1]

# NOTE: Before beginning, you need to adjust things to match the sequence
# adjustments needed at the end of each section in indiv.r.XXX 
# and in the write.csv at the very end. changes should reflect the sequence



#########################################
#### Section 1: EPB #####################
#########################################

# wave attenuation robustness
wave.auc.EPB <- EPB.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "350")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

wave.R.EPB <- wave.auc.EPB/tail(EPB.dat$prop_removed,1)

# shoreline stabilization robustness
shore.auc.EPB <- EPB.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "450")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

shore.R.EPB <- shore.auc.EPB/tail(EPB.dat$prop_removed,1)


# CO2 sequestration robustness
CO2.auc.EPB <- EPB.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "550")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

CO2.R.EPB <- CO2.auc.EPB/tail(EPB.dat$prop_removed,1)

# Water Filtration robustness
waterfilt.auc.EPB <- EPB.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "650")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

waterfilt.R.EPB <- waterfilt.auc.EPB/tail(EPB.dat$prop_removed,1)

# Fishery robustness
Fishery.auc.EPB <- EPB.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "750")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

Fishery.R.EPB <- Fishery.auc.EPB/tail(EPB.dat$prop_removed,1)

# Birdwatching robustness
bird.auc.EPB <- EPB.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "850")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

bird.R.EPB <- bird.auc.EPB/tail(EPB.dat$prop_removed,1)

# waterfowl hunting robustness
hunt.auc.EPB <- EPB.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "950")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

hunt.R.EPB <- hunt.auc.EPB/tail(EPB.dat$prop_removed,1)


indiv.r.EPB <- data.frame(
  sequence = "IndirectREVERSE",
  system = "EPB",
  service = c("wave attenuation","shoreline stabilization","carbon sequestration",
              "water filtration","fishery","birdwatching","waterfowl hunting"),
  robustness = 
    c(wave.R.EPB$auc,shore.R.EPB$auc,CO2.R.EPB$auc,waterfilt.R.EPB$auc,
      Fishery.R.EPB$auc,bird.R.EPB$auc,hunt.R.EPB$auc)
)



#########################################
#### Section 2: CSM #####################
#########################################

# wave attenuation robustness
wave.auc.CSM <- CSM.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "350")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

wave.R.CSM <- wave.auc.CSM/tail(CSM.dat$prop_removed,1)

# shoreline stabilization robustness
shore.auc.CSM <- CSM.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "450")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

shore.R.CSM <- shore.auc.CSM/tail(CSM.dat$prop_removed,1)


# CO2 sequestration robustness
CO2.auc.CSM <- CSM.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "550")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

CO2.R.CSM <- CO2.auc.CSM/tail(CSM.dat$prop_removed,1)

# Water Filtration robustness
waterfilt.auc.CSM <- CSM.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "650")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

waterfilt.R.CSM <- waterfilt.auc.CSM/tail(CSM.dat$prop_removed,1)

# Fishery robustness
Fishery.auc.CSM <- CSM.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "750")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

Fishery.R.CSM <- Fishery.auc.CSM/tail(CSM.dat$prop_removed,1)

# Birdwatching robustness
bird.auc.CSM <- CSM.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "850")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

bird.R.CSM <- bird.auc.CSM/tail(CSM.dat$prop_removed,1)

# waterfowl hunting robustness
hunt.auc.CSM <- CSM.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "950")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

hunt.R.CSM <- hunt.auc.CSM/tail(CSM.dat$prop_removed,1)


indiv.r.CSM <- data.frame(
  sequence = "IndirectREVERSE",
  system = "CSM",
  service = c("wave attenuation","shoreline stabilization","carbon sequestration",
              "water filtration","fishery","birdwatching","waterfowl hunting"),
  robustness = 
    c(wave.R.CSM$auc,shore.R.CSM$auc,CO2.R.CSM$auc,waterfilt.R.CSM$auc,
      Fishery.R.CSM$auc,bird.R.CSM$auc,hunt.R.CSM$auc)
)


#########################################
#### Section 3: BSQ #####################
#########################################

# wave attenuation robustness
wave.auc.BSQ <- BSQ.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "350")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

wave.R.BSQ <- wave.auc.BSQ/tail(BSQ.dat$prop_removed,1)

# shoreline stabilization robustness
shore.auc.BSQ <- BSQ.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "450")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

shore.R.BSQ <- shore.auc.BSQ/tail(BSQ.dat$prop_removed,1)


# CO2 sequestration robustness
CO2.auc.BSQ <- BSQ.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "550")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

CO2.R.BSQ <- CO2.auc.BSQ/tail(BSQ.dat$prop_removed,1)

# Water Filtration robustness
waterfilt.auc.BSQ <- BSQ.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "650")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

waterfilt.R.BSQ <- waterfilt.auc.BSQ/tail(BSQ.dat$prop_removed,1)

# Fishery robustness
Fishery.auc.BSQ <- BSQ.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "750")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

Fishery.R.BSQ <- Fishery.auc.BSQ/tail(BSQ.dat$prop_removed,1)

# Birdwatching robustness
bird.auc.BSQ <- BSQ.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "850")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

bird.R.BSQ <- bird.auc.BSQ/tail(BSQ.dat$prop_removed,1)

# waterfowl hunting robustness
hunt.auc.BSQ <- BSQ.dat %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "950")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / 1)) %>%
  auc_wrapper_given()

hunt.R.BSQ <- hunt.auc.BSQ/tail(BSQ.dat$prop_removed,1)


indiv.r.BSQ <- data.frame(
  sequence = "IndirectREVERSE",
  system = "BSQ",
  service = c("wave attenuation","shoreline stabilization","carbon sequestration",
              "water filtration","fishery","birdwatching","waterfowl hunting"),
  robustness = 
    c(wave.R.BSQ$auc,shore.R.BSQ$auc,CO2.R.BSQ$auc,waterfilt.R.BSQ$auc,
      Fishery.R.BSQ$auc,bird.R.BSQ$auc,hunt.R.BSQ$auc)
)



#########################################
#### Section 4: combine & save ##########
#########################################

all.dat <- rbind(indiv.r.CSM,indiv.r.BSQ,indiv.r.EPB)

write.csv(all.dat, "IndivES_Robust_DirectBiomass_LargeSmall.csv")
