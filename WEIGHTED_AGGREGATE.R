#############################################
######## Weighted ES Robustness Aggregate #####
## Written by: Aislyn Keyes
## Last edited: 13 Nov 2020

# load libraries
library(tidyverse)
library(igraph)
library(dplyr)
library(stringi)

# load functions
source("NEWrobustness_functions.R")
## disable warnings 
options(warn=-1)

## load data Part I ----
# Ecosystem service provider data
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
all.e <- e.esp[,-c(2:4)]
all.e <- distinct(all.e) # pool all esps across services
all.e$total.bio <- sum(all.e$all.biomass) #total ESP biomass

# CSM
c.esp <- c.links %>% subset(Type=="ES") # pull all esps
c.esp <- c.esp %>% subset(ConsumerSpeciesID!="850") # remove birdwatching
colnames(c.esp)[colnames(c.esp)=="ResourceSpeciesID"] <- "SpeciesID"
c.esp <- merge(c.esp,c.bio,by="SpeciesID",all.x=T) # merge to get esp biomass
all.c <- c.esp[,-c(2:3)]
all.c <- distinct(all.c) # pool all esps across services
all.c$total.bio <- sum(all.c$all.biomass) #total ESP biomass

# BSQ
b.esp <- b.links %>% subset(Type=="ES") # pull all esps
b.esp <- b.esp %>% subset(ConsumerSpeciesID!="850") # remove birdwatching
colnames(b.esp)[colnames(b.esp)=="ResourceSpeciesID"] <- "SpeciesID"
b.esp <- merge(b.esp,b.bio,by="SpeciesID",all.x=T) # merge to get esp biomass
all.b <- b.esp[,-c(2:3)]
all.b <- distinct(all.b) # pool all esps across services
all.b$total.bio <- sum(all.b$all.biomass) #total ESP biomass

## Process Output data & Calculate AGGREGATE ES robustness ----
# EPB
epb.dat <- EPB.dat[,c("basal_removed_each","num_removed_tot","prop_removed","name_lost")] # pull direct removed, num_removed and name_lost columns

epb.dat <- epb.dat %>% 
  mutate(SpeciesID = name_lost) # duplicate column so we can separate values by column

epb.dat <- separate_rows(epb.dat, SpeciesID) # seperate rows for indirect losses

# merge output data with biomass data for ESPs only
epb.dat$all.biomass <- e.esp$all.biomass[match(epb.dat$SpeciesID, e.esp$SpeciesID)]
epb.dat[["all.biomass"]][is.na(epb.dat[["all.biomass"]])] <- 0
epb.dat$total.biomass <- sum(epb.dat$all.biomass)
# cumulative sum of esp losses
epb.dat$cum_loss <- cumsum(epb.dat$all.biomass)
epb.dat$prop_ES_lost <- epb.dat$cum_loss/epb.dat$total.biomass
epb.dat$prop_ES_remain <- 1-epb.dat$prop_ES_lost

robust_auc(x=epb.dat$prop_removed, y=epb.dat$prop_ES_remain) 
write.csv(epb.dat, "Weighted-Extension/Aggregate_ES_Outputs/EPB_SuppSppHL_weighted_output.csv")

# CSM
csm.dat <- CSM.dat[,c("basal_removed_each","num_removed_tot","prop_removed","name_lost")] # pull direct removed, num_removed and name_lost columns

csm.dat <- csm.dat %>% 
  mutate(SpeciesID = name_lost) # duplicate column so we can separate values by column

csm.dat <- separate_rows(csm.dat, SpeciesID) # seperate rows for indirect losses

# merge output data with biomass data for ESPs only
csm.dat$all.biomass <- c.esp$all.biomass[match(csm.dat$SpeciesID, c.esp$SpeciesID)]
csm.dat[["all.biomass"]][is.na(csm.dat[["all.biomass"]])] <- 0
csm.dat$total.biomass <- sum(csm.dat$all.biomass)
# cumulative sum of esp losses
csm.dat$cum_loss <- cumsum(csm.dat$all.biomass)
csm.dat$prop_ES_lost <- csm.dat$cum_loss/csm.dat$total.biomass
csm.dat$prop_ES_remain <- 1-csm.dat$prop_ES_lost

robust_auc(x=csm.dat$prop_removed, y=csm.dat$prop_ES_remain) 
write.csv(csm.dat, "Weighted-Extension/Aggregate_ES_Outputs/CSM_SuppSppHL_weighted_output.csv")

# BSQ
bsq.dat <- BSQ.dat[,c("basal_removed_each","num_removed_tot","prop_removed","name_lost")] # pull direct removed, num_removed and name_lost columns

bsq.dat <- bsq.dat %>% 
  mutate(SpeciesID = name_lost) # duplicate column so we can separate values by column

bsq.dat <- separate_rows(bsq.dat, SpeciesID) # seperate rows for indirect losses

# merge output data with biomass data for ESPs only
bsq.dat$all.biomass <- b.esp$all.biomass[match(bsq.dat$SpeciesID, b.esp$SpeciesID)]
bsq.dat[["all.biomass"]][is.na(bsq.dat[["all.biomass"]])] <- 0
bsq.dat$total.biomass <- sum(bsq.dat$all.biomass)
# cumulative sum of esp losses
bsq.dat$cum_loss <- cumsum(bsq.dat$all.biomass)
bsq.dat$prop_ES_lost <- bsq.dat$cum_loss/bsq.dat$total.biomass
bsq.dat$prop_ES_remain <- 1-bsq.dat$prop_ES_lost

robust_auc(x=bsq.dat$prop_removed, y=bsq.dat$prop_ES_remain) 
write.csv(bsq.dat, "Weighted-Extension/Aggregate_ES_Outputs/BSQ_SuppSppHL_weighted_output.csv")


#######################################################
## RANDOM ----
# load data Part II
# Food web data
# CHANGE SALT MARSH DATA FILE
e.nodes <- read.csv("Final_CSM_Nodes_ES.csv") 
e.edges <- read.csv("Final.CSM.edges.filtered.csv")

# create networks and pull adj mat
# CHANGE 
e.nodes2 <- subset(e.nodes, e.nodes$NodeType != "ES") # exclude ES from node list 
e.edges2 <- subset(e.edges, e.edges$Type == "Feeding") # subset interactions to remove all ES links

# create two graph objects: 1. including ES and all interactions, 2. excluding ES and ES interactions
e.netES <- graph.data.frame(e.edges,
                           directed = T,
                           vertices = e.nodes) # create network object WITH services

e.netSPP <- graph.data.frame(e.edges2,
                            directed = T,
                            vertices = e.nodes2) # create network object WITHOUT services

#Remove cannibalism 
e.netSPP <- simplify(e.netSPP, remove.loops = T)
e.netES <- simplify(e.netES, remove.loops=T)

# convert graph objects to adjacency matrices to be used in source functions
e.matES <- get.adjacency(e.netES, sparse = FALSE, attr = NULL) 
e.matSPP <- get.adjacency(e.netSPP, sparse=FALSE, attr = NULL)  


#### Start ALL randomization and robustness ----
# CSM
mat_basal_given <- colnames(e.matSPP)

# create empty data frame to fill with robustness values
rand <- data.frame(matrix(ncol=2,nrow=1000))
x <- c("Randomization","R_ES")
colnames(rand) <- x
rand$Randomization <- seq(1:1000)

# create empty list and fill it with 1000 randomization of sp list
n <- c(1:1000)
sequences <- vector(mode="list",length=1000)
for(i in n){
  sequences[[i]] <- sample(mat_basal_given)
}

mat_output_list_ES <- vector(mode="list",length=1000)
for(i in n) {
  mat_output_list_ES[[i]] <- collapse_wrap_given(N = e.matES, 
                                                 basal = sequences[[i]])
}

for( i in n ){
  mat_output_list_ES[[1]] <- mat_output_list_ES[[i]][,c(1,7,8)] # pull relevent columns
  mat_output_list_ES[[i]] <- mat_output_list_ES[[i]] %>% 
    mutate(SpeciesID = name_lost) # duplicate column so we can separate values by column
  mat_output_list_ES[[i]] <- separate_rows(mat_output_list_ES[[i]], SpeciesID) # seperate rows for indirect losses
  
  # merge output data with biomass data for ESPs only
  mat_output_list_ES[[i]]$all.biomass <- e.esp$all.biomass[match(mat_output_list_ES[[i]]$SpeciesID, all.e$SpeciesID)]
  mat_output_list_ES[[i]][["all.biomass"]][is.na(mat_output_list_ES[[i]][["all.biomass"]])] <- 0
  mat_output_list_ES[[i]]$total.biomass <- sum(mat_output_list_ES[[i]]$all.biomass)
  # cumulative sum of esp losses
  mat_output_list_ES[[i]]$cum_loss <- cumsum( mat_output_list_ES[[i]]$all.biomass)
  mat_output_list_ES[[i]]$prop_ES_lost <-  mat_output_list_ES[[i]]$cum_loss/ mat_output_list_ES[[i]]$total.biomass
  mat_output_list_ES[[i]]$prop_ES_remain <- 1- mat_output_list_ES[[i]]$prop_ES_lost
  
  #calculate robustness
  rand$ES_AUC[i] <- robust_auc(x=mat_output_list_ES[[i]]$prop_removed, y=mat_output_list_ES[[i]]$prop_ES_remain) 
  }


# unlist AUC then calculate robustness and store
rand$R_ES <- unlist(rand$ES_AUC)
rand <- rand[,-4]

min(rand$R_ES)
mean(rand$R_ES)
max(rand$R_ES)
write.csv(rand,"Weighted-Extension/Aggregate_ES_Outputs/CSM_Random_Weighted_Values.csv")


#### Start DIRECT randomization and robustness ----
# CSM
mat_basal_given <- all.c$SpeciesID # change all.X to match system

# create empty data frame to fill with robustness values
rand <- data.frame(matrix(ncol=2,nrow=100))
x <- c("Randomization","R_ES")
colnames(rand) <- x
rand$Randomization <- seq(1:100)

# create empty list and fill it with 1000 randomization of sp list
n <- c(1:100)
sequences <- vector(mode="list",length=100)
for(i in n){
  sequences[[i]] <- sample(mat_basal_given)
}

# create new empty list and fill it with the 100 randomized vectors + extra spp that we'll remove later. Need full web to run the robustness! 
sequences.total <- vector(mode="list",length=100)
for(i in n){
  sequences.total[[i]] <- c(sequences[[i]][which(mat_basal_given %in%
                                                   as.numeric(colnames(e.matSPP)))],# only the given basal species that ARE in the matrix
                            as.numeric(colnames(e.matSPP)[which(!colnames(e.matSPP) %in% mat_basal_given)]))# and also the rest of the species in the matrix that aren't in the basal list
}


mat_output_list_ES <- vector(mode="list",length=100)
for(i in n) {
  mat_output_list_ES[[i]] <- collapse_wrap_given(N = e.matES, 
                                                 basal=sequences.total[[i]])
}

mat_results_list_ES <- vector(mode="list",length=100)
for(i in n){
  mat_results_list_ES[[i]] <- subset(mat_output_list_ES[[i]], basal_removed_each %in% mat_basal_given | prop_remain == 1)
}

target_remove <- length(all.c$SpeciesID) #CHANGE SYSTEM
for( i in seq_along(mat_results_list_ES)){
  mat_results_list_ES[[i]]$prop_removed <- mat_results_list_ES[[i]]$num_removed_tot / target_remove
}

for( i in n ){
  mat_results_list_ES[[1]] <- mat_results_list_ES[[i]][,c(1,7,8)] # pull relevent columns
  mat_results_list_ES[[i]] <- mat_results_list_ES[[i]] %>% 
    mutate(SpeciesID = name_lost) # duplicate column so we can separate values by column
  mat_results_list_ES[[i]] <- separate_rows(mat_results_list_ES[[i]], SpeciesID) # seperate rows for indirect losses
  
  # merge results data with biomass data for ESPs only
  mat_results_list_ES[[i]]$all.biomass <- e.esp$all.biomass[match(mat_results_list_ES[[i]]$SpeciesID, all.e$SpeciesID)]
  mat_results_list_ES[[i]][["all.biomass"]][is.na(mat_results_list_ES[[i]][["all.biomass"]])] <- 0
  mat_results_list_ES[[i]]$total.biomass <- sum(mat_results_list_ES[[i]]$all.biomass)
  # cumulative sum of esp losses
  mat_results_list_ES[[i]]$cum_loss <- cumsum( mat_results_list_ES[[i]]$all.biomass)
  mat_results_list_ES[[i]]$prop_ES_lost <-  mat_results_list_ES[[i]]$cum_loss/ mat_results_list_ES[[i]]$total.biomass
  mat_results_list_ES[[i]]$prop_ES_remain <- 1- mat_results_list_ES[[i]]$prop_ES_lost
  
  #calculate robustness
  rand$ES_AUC[i] <- robust_auc(x=mat_results_list_ES[[i]]$prop_removed, y=mat_results_list_ES[[i]]$prop_ES_remain) 
}

# unlist AUC then calculate robustness and store
rand$R_ES <- unlist(rand$ES_AUC)
rand <- rand[,-4]

mean(rand$R_ES)
write.csv(rand,"Weighted-Extension/Aggregate_ES_Outputs/CSM_RandomESP_Weighted_Values.csv")


#### Start INDIRECT randomization and robustness
# CSM
indirect.spp <- read.csv("CSM_IndirectAll.csv", header=T) 
mat_basal_given <- indirect.spp$SpeciesID

# create empty data frame to fill with robustness values
rand <- data.frame(matrix(ncol=2,nrow=100))
x <- c("Randomization","R_ES")
colnames(rand) <- x
rand$Randomization <- seq(1:100)

# create empty list and fill it with 1000 randomization of sp list
n <- c(1:100)
sequences <- vector(mode="list",length=100)
for(i in n){
  sequences[[i]] <- sample(mat_basal_given)
}

# create new empty list and fill it with the 100 randomized vectors + extra spp that we'll remove later. Need full web to run the robustness! 
sequences.total <- vector(mode="list",length=100)
for(i in n){
  sequences.total[[i]] <- c(sequences[[i]][which(mat_basal_given %in%
                                                   as.numeric(colnames(e.matSPP)))],# only the given basal species that ARE in the matrix
                            as.numeric(colnames(e.matSPP)[which(!colnames(e.matSPP) %in% mat_basal_given)]))# and also the rest of the species in the matrix that aren't in the basal list
}


mat_output_list_ES <- vector(mode="list",length=100)
for(i in n) {
  mat_output_list_ES[[i]] <- collapse_wrap_given(N = e.matES, 
                                                 basal=sequences.total[[i]])
}

mat_results_list_ES <- vector(mode="list",length=100)
for(i in n){
  mat_results_list_ES[[i]] <- subset(mat_output_list_ES[[i]], basal_removed_each %in% mat_basal_given | prop_remain == 1)
}

target_remove <- length(indirect.spp$SpeciesID)
for( i in seq_along(mat_results_list_ES)){
  mat_results_list_ES[[i]]$prop_removed <- mat_results_list_ES[[i]]$num_removed_tot / target_remove
}

for( i in n ){
  mat_results_list_ES[[1]] <- mat_results_list_ES[[i]][,c(1,7,8)] # pull relevent columns
  mat_results_list_ES[[i]] <- mat_results_list_ES[[i]] %>% 
    mutate(SpeciesID = name_lost) # duplicate column so we can separate values by column
  mat_results_list_ES[[i]] <- separate_rows(mat_results_list_ES[[i]], SpeciesID) # seperate rows for indirect losses
  
  # merge results data with biomass data for ESPs only
  mat_results_list_ES[[i]]$all.biomass <- e.esp$all.biomass[match(mat_results_list_ES[[i]]$SpeciesID, all.e$SpeciesID)]
  mat_results_list_ES[[i]][["all.biomass"]][is.na(mat_results_list_ES[[i]][["all.biomass"]])] <- 0
  mat_results_list_ES[[i]]$total.biomass <- sum(mat_results_list_ES[[i]]$all.biomass)
  # cumulative sum of esp losses
  mat_results_list_ES[[i]]$cum_loss <- cumsum( mat_results_list_ES[[i]]$all.biomass)
  mat_results_list_ES[[i]]$prop_ES_lost <-  mat_results_list_ES[[i]]$cum_loss/ mat_results_list_ES[[i]]$total.biomass
  mat_results_list_ES[[i]]$prop_ES_remain <- 1- mat_results_list_ES[[i]]$prop_ES_lost
  
  #calculate robustness
  rand$ES_AUC[i] <- robust_auc(x=mat_results_list_ES[[i]]$prop_removed, y=mat_results_list_ES[[i]]$prop_ES_remain) 
}

# unlist AUC then calculate robustness and store
rand$R_ES <- unlist(rand$ES_AUC)
rand <- rand[,-4]

mean(rand$R_ES)
write.csv(rand,"Weighted-Extension/Aggregate_ES_Outputs/CSM_RandomSuppSPP_Weighted_Values.csv")



