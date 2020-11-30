#############################################
######## Weighted ES Robustness Individual #####
## Written by: Aislyn Keyes
## Last edited: 13 Nov 2020

# load libraries
library(dplyr)

# load functions
source("NEWrobustness_functions.R")
## disable warnings 
options(warn=-1)

## load data ----
# Food web data
# Need the links data to pull ESPs and the biomass data for prop lost.
e.bio <- read.csv("Final Code and Data - Check/EPBBiomass.Estimates.csv") 
e.links <- read.csv("Final.EPB.edges.filtered.csv")

c.bio <- read.csv("Final Code and Data - Check/CSMBiomass.Estimates.csv") 
c.links <- read.csv("Final.CSM.edges.filtered.csv")

b.bio <- read.csv("Final Code and Data - Check/BSQBiomass.Estimates.csv") 
b.links <- read.csv("Final.BSQ.edges.filtered.csv")

# Output data - run this code once per sequence
EPB.dat <- read.csv("EPB_MatOutput_SuppSppHL_ES_Xtarget.csv")
EPB.dat <- EPB.dat[,-1] # get rid of extra column so it doesn't throw off the rest of the code

CSM.dat <- read.csv("CSM_MatOutput_SuppSppHL_ES_Xtarget.csv")
CSM.dat <- CSM.dat[,-1]

BSQ.dat <- read.csv("BSQ_MatOutput_SuppSppHL_ES_Xtarget.csv")
BSQ.dat <- BSQ.dat[,-1]

## Process ESP data----
# ESPs and biomasses
# EPB
e.esp <- e.links %>% subset(Type=="ES") # pull all esps
e.esp <- e.esp %>% subset(ConsumerSpeciesID!="850") # remove birdwatching
colnames(e.esp)[colnames(e.esp)=="ResourceSpeciesID"] <- "SpeciesID"
e.esp <- merge(e.esp,e.bio,by="SpeciesID",all.x=T) # merge to get esp biomass
e.esp <- e.esp[,-c(3,4)]
e.wave <- e.esp %>% subset(ConsumerSpeciesID=="350")
e.wave$total.bio <- sum(e.wave$all.biomass)
e.shore <- e.esp %>% subset(ConsumerSpeciesID=="450")
e.shore$total.bio <- sum(e.shore$all.biomass)
e.carbon <- e.esp %>% subset(ConsumerSpeciesID=="550")
e.carbon$total.bio <- sum(e.carbon$all.biomass)
e.water <- e.esp %>% subset(ConsumerSpeciesID=="650")
e.water$total.bio <- sum(e.water$all.biomass)
e.fish <- e.esp %>% subset(ConsumerSpeciesID=="750")
e.fish$total.bio <- sum(e.fish$all.biomass)
e.hunt <- e.esp %>% subset(ConsumerSpeciesID=="950")
e.hunt$total.bio <- sum(e.hunt$all.biomass)

# CSM
c.esp <- c.links %>% subset(Type=="ES") # pull all esps
c.esp <- c.esp %>% subset(ConsumerSpeciesID!="850") # remove birdwatching
colnames(c.esp)[colnames(c.esp)=="ResourceSpeciesID"] <- "SpeciesID"
c.esp <- merge(c.esp,c.bio,by="SpeciesID",all.x=T) # merge to get esp biomass
c.esp <- c.esp[,-3]
c.wave <- c.esp %>% subset(ConsumerSpeciesID=="350")
c.wave$total.bio <- sum(c.wave$all.biomass)
c.shore <- c.esp %>% subset(ConsumerSpeciesID=="450")
c.shore$total.bio <- sum(c.shore$all.biomass)
c.carbon <- c.esp %>% subset(ConsumerSpeciesID=="550")
c.carbon$total.bio <- sum(c.carbon$all.biomass)
c.water <- c.esp %>% subset(ConsumerSpeciesID=="650")
c.water$total.bio <- sum(c.water$all.biomass)
c.fish <- c.esp %>% subset(ConsumerSpeciesID=="750")
c.fish$total.bio <- sum(c.fish$all.biomass)
c.hunt <- c.esp %>% subset(ConsumerSpeciesID=="950")
c.hunt$total.bio <- sum(c.hunt$all.biomass)


# BSQ
b.esp <- b.links %>% subset(Type=="ES") # pull all esps
b.esp <- b.esp %>% subset(ConsumerSpeciesID!="850") # remove birdwatching
colnames(b.esp)[colnames(b.esp)=="ResourceSpeciesID"] <- "SpeciesID"
b.esp <- merge(b.esp,b.bio,by="SpeciesID",all.x=T) # merge to get esp biomass
b.esp <- b.esp[,-3]
b.wave <- b.esp %>% subset(ConsumerSpeciesID=="350")
b.wave$total.bio <- sum(b.wave$all.biomass)
b.shore <- b.esp %>% subset(ConsumerSpeciesID=="450")
b.shore$total.bio <- sum(b.shore$all.biomass)
b.carbon <- b.esp %>% subset(ConsumerSpeciesID=="550")
b.carbon$total.bio <- sum(b.carbon$all.biomass)
b.water <- b.esp %>% subset(ConsumerSpeciesID=="650")
b.water$total.bio <- sum(b.water$all.biomass)
b.hunt <- b.esp %>% subset(ConsumerSpeciesID=="950")
b.hunt$total.bio <- sum(b.hunt$all.biomass)


## Process Output data & Calculate INDIVIDUAL ES robustness ----
# EPB
epb.dat <- EPB.dat[,c("basal_removed_each","num_removed_tot","prop_removed","name_lost")] # pull direct removed, num_removed and name_lost columns

epb.dat <- epb.dat %>% 
  mutate(SpeciesID = name_lost) # duplicate column so we can separate values by column

epb.dat <- separate_rows(epb.dat, SpeciesID) # seperate rows for indirect losses

## wave atten
epb.dat$wave.lost <- e.wave$all.biomass[match(epb.dat$SpeciesID, e.wave$SpeciesID)]
epb.dat[["wave.lost"]][is.na(epb.dat[["wave.lost"]])] <- 0
epb.dat$total.wave <- sum(epb.dat$wave.lost)
# cumulative sum of esp losses
epb.dat$cum_wave <- cumsum(epb.dat$wave.lost)
epb.dat$prop_wave_lost <- epb.dat$cum_wave/epb.dat$total.wave
epb.dat$prop_wave_remain <- 1-epb.dat$prop_wave_lost

robust_auc(x=epb.dat$prop_removed, y=epb.dat$prop_wave_remain) 

## shoreline stab
epb.dat$shore.lost <- e.shore$all.biomass[match(epb.dat$SpeciesID, e.shore$SpeciesID)]
epb.dat[["shore.lost"]][is.na(epb.dat[["shore.lost"]])] <- 0
epb.dat$total.shore <- sum(epb.dat$shore.lost)
# cumulative sum of esp losses
epb.dat$cum_shore <- cumsum(epb.dat$shore.lost)
epb.dat$prop_shore_lost <- epb.dat$cum_shore/epb.dat$total.shore
epb.dat$prop_shore_remain <- 1-epb.dat$prop_shore_lost

robust_auc(x=epb.dat$prop_removed, y=epb.dat$prop_shore_remain) 

## carbon
epb.dat$carbon.lost <- e.carbon$all.biomass[match(epb.dat$SpeciesID, e.carbon$SpeciesID)]
epb.dat[["carbon.lost"]][is.na(epb.dat[["carbon.lost"]])] <- 0
epb.dat$total.carbon <- sum(epb.dat$carbon.lost)
# cumulative sum of esp losses
epb.dat$cum_carbon <- cumsum(epb.dat$carbon.lost)
epb.dat$prop_carbon_lost <- epb.dat$cum_carbon/epb.dat$total.carbon
epb.dat$prop_carbon_remain <- 1-epb.dat$prop_carbon_lost

robust_auc(x=epb.dat$prop_removed, y=epb.dat$prop_carbon_remain) 

## water filtration
epb.dat$water.lost <- e.water$all.biomass[match(epb.dat$SpeciesID, e.water$SpeciesID)]
epb.dat[["water.lost"]][is.na(epb.dat[["water.lost"]])] <- 0
epb.dat$total.water <- sum(epb.dat$water.lost)
# cumulative sum of esp losses
epb.dat$cum_water <- cumsum(epb.dat$water.lost)
epb.dat$prop_water_lost <- epb.dat$cum_water/epb.dat$total.water
epb.dat$prop_water_remain <- 1-epb.dat$prop_water_lost

robust_auc(x=epb.dat$prop_removed, y=epb.dat$prop_water_remain) 

## fishery
epb.dat$fish.lost <- e.fish$all.biomass[match(epb.dat$SpeciesID, e.fish$SpeciesID)]
epb.dat[["fish.lost"]][is.na(epb.dat[["fish.lost"]])] <- 0
epb.dat$total.fish <- sum(epb.dat$fish.lost)
# cumulative sum of esp losses
epb.dat$cum_fish <- cumsum(epb.dat$fish.lost)
epb.dat$prop_fish_lost <- epb.dat$cum_fish/epb.dat$total.fish
epb.dat$prop_fish_remain <- 1-epb.dat$prop_fish_lost

robust_auc(x=epb.dat$prop_removed, y=epb.dat$prop_fish_remain) 

## waterfowl hunt
epb.dat$hunt.lost <- e.hunt$all.biomass[match(epb.dat$SpeciesID, e.hunt$SpeciesID)]
epb.dat[["hunt.lost"]][is.na(epb.dat[["hunt.lost"]])] <- 0
epb.dat$total.hunt <- sum(epb.dat$hunt.lost)
# cumulative sum of esp losses
epb.dat$cum_hunt <- cumsum(epb.dat$hunt.lost)
epb.dat$prop_hunt_lost <- epb.dat$cum_hunt/epb.dat$total.hunt
epb.dat$prop_hunt_remain <- 1-epb.dat$prop_hunt_lost

robust_auc(x=epb.dat$prop_removed, y=epb.dat$prop_hunt_remain) 

write.csv(epb.dat, "Weighted-Extension/Individual_ES_Outputs/EPB_SuppSppHL_individual_output.csv")

# CSM
csm.dat <- CSM.dat[,c("basal_removed_each","num_removed_tot","prop_removed","name_lost")] # pull direct removed, num_removed and name_lost columns

csm.dat <- csm.dat %>% 
  mutate(SpeciesID = name_lost) # duplicate column so we can separate values by column

csm.dat <- separate_rows(csm.dat, SpeciesID) # seperate rows for indirect losses

## wave atten
csm.dat$wave.lost <- e.wave$all.biomass[match(csm.dat$SpeciesID, e.wave$SpeciesID)]
csm.dat[["wave.lost"]][is.na(csm.dat[["wave.lost"]])] <- 0
csm.dat$total.wave <- sum(csm.dat$wave.lost)
# cumulative sum of esp losses
csm.dat$cum_wave <- cumsum(csm.dat$wave.lost)
csm.dat$prop_wave_lost <- csm.dat$cum_wave/csm.dat$total.wave
csm.dat$prop_wave_remain <- 1-csm.dat$prop_wave_lost

robust_auc(x=csm.dat$prop_removed, y=csm.dat$prop_wave_remain) 

## shoreline stab
csm.dat$shore.lost <- e.shore$all.biomass[match(csm.dat$SpeciesID, e.shore$SpeciesID)]
csm.dat[["shore.lost"]][is.na(csm.dat[["shore.lost"]])] <- 0
csm.dat$total.shore <- sum(csm.dat$shore.lost)
# cumulative sum of esp losses
csm.dat$cum_shore <- cumsum(csm.dat$shore.lost)
csm.dat$prop_shore_lost <- csm.dat$cum_shore/csm.dat$total.shore
csm.dat$prop_shore_remain <- 1-csm.dat$prop_shore_lost

robust_auc(x=csm.dat$prop_removed, y=csm.dat$prop_shore_remain) 

## carbon
csm.dat$carbon.lost <- e.carbon$all.biomass[match(csm.dat$SpeciesID, e.carbon$SpeciesID)]
csm.dat[["carbon.lost"]][is.na(csm.dat[["carbon.lost"]])] <- 0
csm.dat$total.carbon <- sum(csm.dat$carbon.lost)
# cumulative sum of esp losses
csm.dat$cum_carbon <- cumsum(csm.dat$carbon.lost)
csm.dat$prop_carbon_lost <- csm.dat$cum_carbon/csm.dat$total.carbon
csm.dat$prop_carbon_remain <- 1-csm.dat$prop_carbon_lost

robust_auc(x=csm.dat$prop_removed, y=csm.dat$prop_carbon_remain) 

## water filtration
csm.dat$water.lost <- e.water$all.biomass[match(csm.dat$SpeciesID, e.water$SpeciesID)]
csm.dat[["water.lost"]][is.na(csm.dat[["water.lost"]])] <- 0
csm.dat$total.water <- sum(csm.dat$water.lost)
# cumulative sum of esp losses
csm.dat$cum_water <- cumsum(csm.dat$water.lost)
csm.dat$prop_water_lost <- csm.dat$cum_water/csm.dat$total.water
csm.dat$prop_water_remain <- 1-csm.dat$prop_water_lost

robust_auc(x=csm.dat$prop_removed, y=csm.dat$prop_water_remain) 

## fishery
csm.dat$fish.lost <- e.fish$all.biomass[match(csm.dat$SpeciesID, e.fish$SpeciesID)]
csm.dat[["fish.lost"]][is.na(csm.dat[["fish.lost"]])] <- 0
csm.dat$total.fish <- sum(csm.dat$fish.lost)
# cumulative sum of esp losses
csm.dat$cum_fish <- cumsum(csm.dat$fish.lost)
csm.dat$prop_fish_lost <- csm.dat$cum_fish/csm.dat$total.fish
csm.dat$prop_fish_remain <- 1-csm.dat$prop_fish_lost

robust_auc(x=csm.dat$prop_removed, y=csm.dat$prop_fish_remain) 

## waterfowl hunt
csm.dat$hunt.lost <- e.hunt$all.biomass[match(csm.dat$SpeciesID, e.hunt$SpeciesID)]
csm.dat[["hunt.lost"]][is.na(csm.dat[["hunt.lost"]])] <- 0
csm.dat$total.hunt <- sum(csm.dat$hunt.lost)
# cumulative sum of esp losses
csm.dat$cum_hunt <- cumsum(csm.dat$hunt.lost)
csm.dat$prop_hunt_lost <- csm.dat$cum_hunt/csm.dat$total.hunt
csm.dat$prop_hunt_remain <- 1-csm.dat$prop_hunt_lost

robust_auc(x=csm.dat$prop_removed, y=csm.dat$prop_hunt_remain) 

write.csv(csm.dat, "Weighted-Extension/Individual_ES_Outputs/CSM_SuppSppHL_individual_output.csv")

# BSQ
bsq.dat <- BSQ.dat[,c("basal_removed_each","num_removed_tot","prop_removed","name_lost")] # pull direct removed, num_removed and name_lost columns

bsq.dat <- bsq.dat %>% 
  mutate(SpeciesID = name_lost) # duplicate column so we can separate values by column

bsq.dat <- separate_rows(bsq.dat, SpeciesID) # seperate rows for indirect losses

## wave atten
bsq.dat$wave.lost <- e.wave$all.biomass[match(bsq.dat$SpeciesID, e.wave$SpeciesID)]
bsq.dat[["wave.lost"]][is.na(bsq.dat[["wave.lost"]])] <- 0
bsq.dat$total.wave <- sum(bsq.dat$wave.lost)
# cumulative sum of esp losses
bsq.dat$cum_wave <- cumsum(bsq.dat$wave.lost)
bsq.dat$prop_wave_lost <- bsq.dat$cum_wave/bsq.dat$total.wave
bsq.dat$prop_wave_remain <- 1-bsq.dat$prop_wave_lost

robust_auc(x=bsq.dat$prop_removed, y=bsq.dat$prop_wave_remain) 

## shoreline stab
bsq.dat$shore.lost <- e.shore$all.biomass[match(bsq.dat$SpeciesID, e.shore$SpeciesID)]
bsq.dat[["shore.lost"]][is.na(bsq.dat[["shore.lost"]])] <- 0
bsq.dat$total.shore <- sum(bsq.dat$shore.lost)
# cumulative sum of esp losses
bsq.dat$cum_shore <- cumsum(bsq.dat$shore.lost)
bsq.dat$prop_shore_lost <- bsq.dat$cum_shore/bsq.dat$total.shore
bsq.dat$prop_shore_remain <- 1-bsq.dat$prop_shore_lost

robust_auc(x=bsq.dat$prop_removed, y=bsq.dat$prop_shore_remain) 

## carbon
bsq.dat$carbon.lost <- e.carbon$all.biomass[match(bsq.dat$SpeciesID, e.carbon$SpeciesID)]
bsq.dat[["carbon.lost"]][is.na(bsq.dat[["carbon.lost"]])] <- 0
bsq.dat$total.carbon <- sum(bsq.dat$carbon.lost)
# cumulative sum of esp losses
bsq.dat$cum_carbon <- cumsum(bsq.dat$carbon.lost)
bsq.dat$prop_carbon_lost <- bsq.dat$cum_carbon/bsq.dat$total.carbon
bsq.dat$prop_carbon_remain <- 1-bsq.dat$prop_carbon_lost

robust_auc(x=bsq.dat$prop_removed, y=bsq.dat$prop_carbon_remain) 

## water filtration
bsq.dat$water.lost <- e.water$all.biomass[match(bsq.dat$SpeciesID, e.water$SpeciesID)]
bsq.dat[["water.lost"]][is.na(bsq.dat[["water.lost"]])] <- 0
bsq.dat$total.water <- sum(bsq.dat$water.lost)
# cumulative sum of esp losses
bsq.dat$cum_water <- cumsum(bsq.dat$water.lost)
bsq.dat$prop_water_lost <- bsq.dat$cum_water/bsq.dat$total.water
bsq.dat$prop_water_remain <- 1-bsq.dat$prop_water_lost

robust_auc(x=bsq.dat$prop_removed, y=bsq.dat$prop_water_remain) 

## waterfowl hunt
bsq.dat$hunt.lost <- e.hunt$all.biomass[match(bsq.dat$SpeciesID, e.hunt$SpeciesID)]
bsq.dat[["hunt.lost"]][is.na(bsq.dat[["hunt.lost"]])] <- 0
bsq.dat$total.hunt <- sum(bsq.dat$hunt.lost)
# cumulative sum of esp losses
bsq.dat$cum_hunt <- cumsum(bsq.dat$hunt.lost)
bsq.dat$prop_hunt_lost <- bsq.dat$cum_hunt/bsq.dat$total.hunt
bsq.dat$prop_hunt_remain <- 1-bsq.dat$prop_hunt_lost

robust_auc(x=bsq.dat$prop_removed, y=bsq.dat$prop_hunt_remain) 

write.csv(bsq.dat, "Weighted-Extension/Individual_ES_Outputs/BSQ_SuppSppHL_individual_output.csv")
