##########################################
######## Final robustness calculations 
### Written by: Aislyn Keyes using functions written by Allison Barner
### Last Edited: 30 November 2020


# load packages
library(tidyverse)
library(igraph)
library(dplyr)
library(stringi)

## AB source functions for analysis----
source("NEWrobustness_functions.R")
## disable warnings 
options(warn=-1)


## load data ----
# For each of these, update to match the file name of the node, edge lists and sequence files
nodes.ES <- read.csv("Final_CSM_Nodes_ES.csv") #node list
edges.ES <- read.csv("Final.CSM.edges.filtered.csv") # edgelist
threats <- read.csv("CSM.RealisticSeq_AllThreats.csv", header=T) # threatened species sequence
rare.seq <- read.csv("RareCSM_seq_FINAL.csv", header=T) # rarity sequence
indirect.spp <- read.csv("CSM_IndirectAll.csv", header=T) # supporting species sequence
node.info <- read.csv("Final_CSM_Nodes_ES.csv", header=T) # will be used for abundance/biomass data
biomass.dat <- read.csv("CSMBiomass.Estimates.csv") # biomass data, 
# estimated using equations from Hechinger et al. (2011) A common scaling rule for abundance, energetics and prod.
# of parasitic and free-living species
service_nodes <- c(350,450,550,650,750,850,950) # for tracking in ES robustness calculation 

## food web data processing ----
nodes.spp <- subset(nodes.ES, nodes.ES$NodeType != "ES") # exclude ES from node list 
edges.spp <- subset(edges.ES, edges.ES$Type == "Feeding") # subset interactions to remove all ES links

# create two graph objects: 1. including ES and all interactions, 2. excluding ES and ES interactions
net.ES <- graph.data.frame(edges.ES,
                           directed = T,
                           vertices = nodes.ES) # create network object WITH services

net.spp <- graph.data.frame(edges.spp,
                            directed = T,
                            vertices = nodes.spp) # create network object WITHOUT services

# Remove cannibalism 
net.spp <- simplify(net.spp, remove.loops = T)
net.ES <- simplify(net.ES, remove.loops=T)

# convert graph objects to adjacency matrices to be used in source functions by AB
mat.ES <- get.adjacency(net.ES, sparse = FALSE, attr = NULL) 
mat.spp <- get.adjacency(net.spp, sparse=FALSE, attr = NULL) 


## sequence data processing ----
# organize the sequence data

# RARITY
rare.seq <- rare.seq[,-1] #remove column 1 to keep just SpeciesID

# INDIRECT ESP
indirect.seq <- indirect.spp[,1] # select column 1 for SpeciesID

# THREATENED
threats <- subset(threats, threats$Mortality == 1) # select species that have known mortality
threats <- as.data.frame(threats[,1]) # keep SpeciesID column
colnames(threats)[colnames(threats)=="threats[, 1]"] <- "SpeciesID" # rename to match node.info
threats <- distinct(threats) # remove duplicate entries
attach(node.info)
node.info2 <- node.info[,c("SpeciesID","Abundance.no..ha.")] # select relevent columns
detach(node.info)
threat.spp <- merge(threats,node.info2,by="SpeciesID") # merge 2 df by Species ID
attach(threat.spp)
threat.spp <- threat.spp[order(Abundance.no..ha.),] # sort by abundance
threat.seq <- threat.spp[,1] # select just the SpeciesID

# DIRECT ESP
direct <- subset(edges.ES, edges.ES$Type == "ES") # select interactions between species and ES to get direct ESPs
direct <- as.data.frame(direct[,1]) # # select only the "ResourceSpeciesID" to get the ESPs
colnames(direct)[colnames(direct)=="direct[, 1]"] <- "SpeciesID" # rename to match node.info
direct <- distinct(direct) # remove duplicate entries

direct.spp <- merge(direct,biomass.dat,by="SpeciesID", all.x=T) # merge 2 df by Species ID
attach(direct.spp)
direct.spp.bio <- direct.spp[order(all.biomass),] # sort by biomass
direct.bio.seq <- direct.spp.bio[,1] # select spp id column for biomass sequence
direct.rand.seq <- direct.spp$SpeciesID


## nonbasal and basal species for recalculating the x and y axis
# y axis of food web robustness should only consider species that can be secondarily lost (ie, no basal species)
resources <- data.frame(SpeciesID = nodes.spp$SpeciesID, InDegree = igraph::degree(net.spp, mode="in")) # calculate number of resources for each species
susc.spp <- resources[resources$InDegree>0,] # subset to the nonbasal, susceptible species
basal.spp <- resources[resources$InDegree==0,] # subset to the basal, NOT susceptible species
tot.susc <- nrow(susc.spp) # denominator for all sequences on food web y axis


###################################################################################

## Begin extinction sequences

## Most to least connected ----
## set sequence
mat_basal <- colnames(mat.spp) # list of species to remove, order is set in the collapes_wrap function for this sequence

## exclude ES, food web only
mat_output <- collapse_wrap(N = mat.spp, basal = mat_basal, n_rand = 1)

mat_results <- subset(mat_output, mat_output$id == "high") # select only the output for most-least connected sequence

# ADJUST Y AXIS
mat_results$nodes_susc <- tot.susc
mat_results$nontarget_lost_each <- stri_count_fixed(mat_results$name_lost,";") # number species LOST each step
mat_results$cum_nontarget_lost <- cumsum(mat_results$nontarget_lost_each)

mat_results$Y <- (mat_results$nodes_susc - mat_results$cum_nontarget_lost) / mat_results$nodes_susc

# ADJUST X AXIS
target_remove <- nrow(nodes.spp)
mat_results$prop_target_removed <- mat_results$num_removed_tot / target_remove

# plot to visualize
root <- ggplot(mat_results, aes(x=prop_target_removed, y=Y))
(root + geom_point(shape=19) 
  + geom_line() + theme_bw(base_size = 14) 
  + xlim(0,1) + ylim(0,1)
  + xlab("Proportion of target species removed")
  + ylab("Proportion of susceptible species remaining")
)

# calc robustness (AUC)
robust_auc(x = mat_results$prop_target_removed, y = mat_results$Y) 

write.csv(mat_results, "CSM_MatOutput_MLC_Xtarget.csv") #save output if desired

## include ES
mat_output_ES <- collapse_wrap(N = mat.ES, basal = mat_basal, n_rand = 1)

mat_results_ES <- subset(mat_output_ES, mat_output_ES$degree_type == "high") # select only the output for most-least connected sequence


# make a column to track number of services lost at each step
mat_results_ES$service_lost <- ifelse(str_count(string = mat_results_ES$name_lost, pattern = "350|450|550|650|750|850|950"), str_count(mat_results_ES$name_lost,".50") , NA)

# track number services remain to use in propES_remain column
mat_results_ES[["service_lost"]][is.na(mat_results_ES[["service_lost"]])] <- 0
mat_results_ES$cum_nontarget_lost <- cumsum(mat_results_ES$service_lost)
mat_results_ES$ES_remain <- c(7-cumsum(mat_results_ES$service_lost))
mat_results_ES$propES_remain <- c(mat_results_ES$ES_remain/7)

# ADJUST X AXIS
target_remove <- nrow(nodes.spp)
mat_results_ES$prop_removed <- mat_results_ES$num_removed_tot / target_remove

# revise output to get ES AUC 
# this tracks the secondary loss of ecosystem service nodes instead of species
mat_results_ES %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "350|450|550|650|750|850|950")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / length(service_nodes))) %>%
  auc_wrapper()

root <- ggplot(mat_results_ES, aes(x=prop_removed, y=propES_remain))
(root + geom_point(shape=19) 
  + geom_line() + theme_bw(base_size = 14) 
  + xlim(0,1) + ylim(0,1)
  + xlab("Proportion of target species removed")
  + ylab("Proportion of ecosystem services remaining")
)

write.csv(mat_results_ES, "CSM_MatOutput_MLC_ES_Xtarget.csv")


## Rarity --------------------------------------------------------------------------
## set sequence
mat_basal_given <- rare.seq # list of species to remove, order is set in the collapes_wrap function for this sequence

## exclude ES, food web only
mat_results <- collapse_wrap_given(N = mat.spp, basal = mat_basal_given) 

# ADJUST Y AXIS
mat_results$nodes_susc <- tot.susc
mat_results$nontarget_lost_each <- stri_count_fixed(mat_results$name_lost,";") # number species LOST each step
mat_results$cum_nontarget_lost <- cumsum(mat_results$nontarget_lost_each)

mat_results$Y <- (mat_results$nodes_susc - mat_results$cum_nontarget_lost) / mat_results$nodes_susc


# ADJUST X AXIS
target_remove <- nrow(nodes.spp)
mat_results$prop_target_removed <- mat_results$num_removed_tot / target_remove


root <- ggplot(mat_results, aes(x=prop_target_removed, y=Y))
(root + geom_point(shape=19) 
  + geom_line() + theme_bw(base_size = 14) 
  + xlim(0,1) + ylim(0,1)
  + xlab("Proportion of target species removed")
  + ylab("Proportion of susceptible species remaining")
)

# calc robustness (AUC)
robust_auc(x = mat_results$prop_target_removed, y = mat_results$Y) 

write.csv(mat_results, "CSM_MatOutput_Rarity_Xtarget.csv")

## include ES
mat_results_ES <- collapse_wrap_given(N = mat.ES, basal = mat_basal_given) # run ext seq like normal


# make a column to track number of services lost at each step
mat_results_ES$service_lost <- ifelse(str_count(string = mat_results_ES$name_lost, pattern = "350|450|550|650|750|850|950"), str_count(mat_results_ES$name_lost,".50") , NA)

# track number services remain to use in propES_remain column
mat_results_ES[["service_lost"]][is.na(mat_results_ES[["service_lost"]])] <- 0
mat_results_ES$cum_nontarget_lost <- cumsum(mat_results_ES$service_lost)
mat_results_ES$ES_remain <- c(7-cumsum(mat_results_ES$service_lost))
mat_results_ES$propES_remain <- c(mat_results_ES$ES_remain/7)

# ADJUST X AXIS
target_remove <- nrow(nodes.spp)
mat_results_ES$prop_removed <- mat_results_ES$num_removed_tot / target_remove

# revise output to get ES AUC
mat_results_ES %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "350|450|550|650|750|850|950")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / length(service_nodes))) %>%
  auc_wrapper_given()

root <- ggplot(mat_results_ES, aes(x=prop_removed, y=propES_remain))
(root + geom_point(shape=19) 
  + geom_line() + theme_bw(base_size = 14) 
  + xlim(0,1) + ylim(0,1)
  + xlab("Proportion of target species removed")
  + ylab("Proportion of ecosystem services remaining")
)

write.csv(mat_results_ES, "CSM_MatOutput_Rarity_ES_Xtarget.csv")


## Threatened --------------------------------------------------------------------------
## set sequence
mat_basal_given <- threat.seq 
mat_basal_given_total <- c(mat_basal_given[which(mat_basal_given %in%
                                                   as.numeric(colnames(mat.spp)))],# only the given basal species that ARE in the matrix
                           as.numeric(colnames(mat.spp)[which(!colnames(mat.spp) %in% mat_basal_given)]))# and also the rest of the species in the matrix that aren't in the basal list

## exclude ES, food web only
mat_results <- collapse_wrap_given(N = mat.spp, basal = mat_basal_given_total)


# ADJUST Y AXIS
mat_results$nodes_susc <- tot.susc
mat_results$nontarget_lost_each <- stri_count_fixed(mat_results$name_lost,";") # number species LOST each step
mat_results$cum_nontarget_lost <- cumsum(mat_results$nontarget_lost_each)

mat_results$Y <- (mat_results$nodes_susc - mat_results$cum_nontarget_lost) / mat_results$nodes_susc


# ADJUST X AXIS
target_remove <- length(mat_basal_given)
mat_results$prop_target_removed <- mat_results$num_removed_tot / target_remove

# we want to get the robustness of the system for ONLY the sequence we're interested in, not the full removal
mat_results <- subset(mat_results, basal_removed_each %in% mat_basal_given | prop_remain == 1)


root <- ggplot(mat_results, aes(x=prop_target_removed, y=Y))
(root + geom_point(shape=19) 
  + geom_line() + theme_bw(base_size = 14) 
  + xlim(0,1) + ylim(0,1)
  + xlab("Proportion of target species removed")
  + ylab("Proportion of susceptible species remaining")
)

# calc robustness (AUC)
robust_auc(x = mat_results$prop_target_removed, y = mat_results$Y) 

write.csv(mat_results, "CSM_MatOutput_Threat_Xtarget.csv")

## include ES
mat_results_ES <- collapse_wrap_given(N = mat.ES, basal = mat_basal_given_total) # run ext seq like normal

# we want to get the robustness of the system for ONLY the sequence we're interested in, not the full removal
mat_results_ES <- subset(mat_results_ES, basal_removed_each %in% mat_basal_given | prop_remain == 1)

# make a column to track number of services lost at each step
mat_results_ES$service_lost <- ifelse(str_count(string = mat_results_ES$name_lost, pattern = "350|450|550|650|750|850|950"), str_count(mat_results_ES$name_lost,".50") , NA)

# track number services remain to use in propES_remain column
mat_results_ES[["service_lost"]][is.na(mat_results_ES[["service_lost"]])] <- 0
mat_results_ES$cum_nontarget_lost <- cumsum(mat_results_ES$service_lost)
mat_results_ES$ES_remain <- c(7-cumsum(mat_results_ES$service_lost))
mat_results_ES$propES_remain <- c(mat_results_ES$ES_remain/7)

# ADJUST X AXIS
target_remove <- length(mat_basal_given)
mat_results_ES$prop_removed <- mat_results_ES$num_removed_tot / target_remove

# revise output to get ES AUC
mat_results_ES %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "350|450|550|650|750|850|950")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / length(service_nodes))) %>%
  auc_wrapper_given()

root <- ggplot(mat_results_ES, aes(x=prop_removed, y=propES_remain))
(root + geom_point(shape=19) 
  + geom_line() + theme_bw(base_size = 14) 
  + xlim(0,1) + ylim(0,1)
  + xlab("Proportion of target species removed")
  + ylab("Proportion of ecosystem services remaining")
)

write.csv(mat_results_ES, "CSM_MatOutput_Threat_ES_Xtarget.csv")


## LS ESP --------------------------------------------------------------------------
## set sequence
mat_basal_given <-  rev(direct.bio.seq)
mat_basal_given_total <- c(mat_basal_given[which(mat_basal_given %in%
                                                   as.numeric(colnames(mat.spp)))],# only the given basal species that ARE in the matrix
                           as.numeric(colnames(mat.spp)[which(!colnames(mat.spp) %in% mat_basal_given)]))# and also the rest of the species in the matrix that aren't in the basal list

## exclude ES, food web only
mat_results <- collapse_wrap_given(N = mat.spp, basal = mat_basal_given_total)


# ADJUST Y AXIS
mat_results$nodes_susc <- tot.susc
mat_results$nontarget_lost_each <- stri_count_fixed(mat_results$name_lost,";") # number species LOST each step
mat_results$cum_nontarget_lost <- cumsum(mat_results$nontarget_lost_each)

mat_results$Y <- (mat_results$nodes_susc - mat_results$cum_nontarget_lost) / mat_results$nodes_susc


# ADJUST X AXIS
target_remove <- length(mat_basal_given)
mat_results$prop_target_removed <- mat_results$num_removed_tot / target_remove



# we want to get the robustness of the system for ONLY the sequence we're interested in, not the full removal
mat_results <- subset(mat_results, basal_removed_each %in% mat_basal_given | prop_remain == 1)

write.csv(mat_results, "CSM_MatOutput_ESPls_Xtarget.csv")


root <- ggplot(mat_results, aes(x=prop_target_removed, y=Y))
(root + geom_point(shape=19) 
  + geom_line() + theme_bw(base_size = 14) 
  + xlim(0,1) + ylim(0,1)
  + xlab("Proportion of target species removed")
  + ylab("Proportion of susceptible species remaining")
)

# calc robustness (AUC)
robust_auc(x = mat_results$prop_target_removed, y = mat_results$Y) 

## include ES
mat_results_ES <- collapse_wrap_given(N = mat.ES, basal = mat_basal_given_total) # run ext seq like normal

# we want to get the robustness of the system for ONLY the sequence we're interested in, not the full removal
mat_results_ES <- subset(mat_results_ES, basal_removed_each %in% mat_basal_given | prop_remain == 1)

# make a column to track number of services lost at each step
mat_results_ES$service_lost <- ifelse(str_count(string = mat_results_ES$name_lost, pattern = "350|450|550|650|750|850|950"), str_count(mat_results_ES$name_lost,".50") , NA)

# track number services remain to use in propES_remain column
mat_results_ES[["service_lost"]][is.na(mat_results_ES[["service_lost"]])] <- 0
mat_results_ES$cum_nontarget_lost <- cumsum(mat_results_ES$service_lost)
mat_results_ES$ES_remain <- c(7-cumsum(mat_results_ES$service_lost))
mat_results_ES$propES_remain <- c(mat_results_ES$ES_remain/7)

# ADJUST X AXIS
target_remove <- length(mat_basal_given)
mat_results_ES$prop_removed <- mat_results_ES$num_removed_tot / target_remove

# revise output to get ES AUC
mat_results_ES %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "350|450|550|650|750|850|950")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / length(service_nodes))) %>%
  auc_wrapper_given()

root <- ggplot(mat_results_ES, aes(x=prop_removed, y=propES_remain))
(root + geom_point(shape=19) 
  + geom_line() + theme_bw(base_size = 14) 
  + xlim(0,1) + ylim(0,1)
  + xlab("Proportion of target species removed")
  + ylab("Proportion of ecosystem services remaining")
)

write.csv(mat_results_ES, "CSM_MatOutput_ESPls_ES_Xtarget.csv")

## SL ESP --------------------------------------------------------------------------
## set sequence
mat_basal_given <-  direct.bio.seq
mat_basal_given_total <- c(mat_basal_given[which(mat_basal_given %in%
                                                   as.numeric(colnames(mat.spp)))],# only the given basal species that ARE in the matrix
                           as.numeric(colnames(mat.spp)[which(!colnames(mat.spp) %in% mat_basal_given)]))# and also the rest of the species in the matrix that aren't in the basal list

## exclude ES, food web only
mat_results <- collapse_wrap_given(N = mat.spp, basal = mat_basal_given_total)


# ADJUST Y AXIS
mat_results$nodes_susc <- tot.susc
mat_results$nontarget_lost_each <- stri_count_fixed(mat_results$name_lost,";") # number species LOST each step
mat_results$cum_nontarget_lost <- cumsum(mat_results$nontarget_lost_each)

mat_results$Y <- (mat_results$nodes_susc - mat_results$cum_nontarget_lost) / mat_results$nodes_susc


# ADJUST X AXIS
target_remove <- length(mat_basal_given)
mat_results$prop_target_removed <- mat_results$num_removed_tot / target_remove


# we want to get the robustness of the system for ONLY the sequence we're interested in, not the full removal
mat_results <- subset(mat_results, basal_removed_each %in% mat_basal_given | prop_remain == 1)
write.csv(mat_results, "CSM_MatOutput_ESPsl_Xtarget.csv")


root <- ggplot(mat_results, aes(x=prop_target_removed, y=Y))
(root + geom_point(shape=19) 
  + geom_line() + theme_bw(base_size = 14) 
  + xlim(0,1) + ylim(0,1)
  + xlab("Proportion of target species removed")
  + ylab("Proportion of susceptible species remaining")
)

# calc robustness (AUC)
robust_auc(x = mat_results$prop_target_removed, y = mat_results$Y) 


## include ES
mat_results_ES <- collapse_wrap_given(N = mat.ES, basal = mat_basal_given_total) # run ext seq like normal

# we want to get the robustness of the system for ONLY the sequence we're interested in, not the full removal
mat_results_ES <- subset(mat_results_ES, basal_removed_each %in% mat_basal_given | prop_remain == 1)

# make a column to track number of services lost at each step
mat_results_ES$service_lost <- ifelse(str_count(string = mat_results_ES$name_lost, pattern = "350|450|550|650|750|850|950"), str_count(mat_results_ES$name_lost,".50") , NA)

# track number services remain to use in propES_remain column
mat_results_ES[["service_lost"]][is.na(mat_results_ES[["service_lost"]])] <- 0
mat_results_ES$cum_nontarget_lost <- cumsum(mat_results_ES$service_lost)
mat_results_ES$ES_remain <- c(7-cumsum(mat_results_ES$service_lost))
mat_results_ES$propES_remain <- c(mat_results_ES$ES_remain/7)

# ADJUST X AXIS
target_remove <- length(mat_basal_given)
mat_results_ES$prop_removed <- mat_results_ES$num_removed_tot / target_remove

# revise output to get ES AUC
mat_results_ES %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "350|450|550|650|750|850|950")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / length(service_nodes))) %>%
  auc_wrapper_given()

root <- ggplot(mat_results_ES, aes(x=prop_removed, y=propES_remain))
(root + geom_point(shape=19) 
  + geom_line() + theme_bw(base_size = 14) 
  + xlim(0,1) + ylim(0,1)
  + xlab("Proportion of target species removed")
  + ylab("Proportion of ecosystem services remaining")
)

write.csv(mat_results_ES, "CSM_MatOutput_ESPsl_ES_Xtarget.csv")

## HL Supporting  --------------------------------------------------------------------------
## set sequence
mat_basal_given <-  indirect.seq 
mat_basal_given_total <- c(mat_basal_given[which(mat_basal_given %in%
                                                   as.numeric(colnames(mat.spp)))],# only the given basal species that ARE in the matrix
                           as.numeric(colnames(mat.spp)[which(!colnames(mat.spp) %in% mat_basal_given)]))# and also the rest of the species in the matrix that aren't in the basal list

## exclude ES, food web only
mat_results <- collapse_wrap_given(N = mat.spp, basal = mat_basal_given_total)


# ADJUST Y AXIS
mat_results$nodes_susc <- tot.susc
mat_results$nontarget_lost_each <- stri_count_fixed(mat_results$name_lost,";") # number species LOST each step
mat_results$cum_nontarget_lost <- cumsum(mat_results$nontarget_lost_each)

mat_results$Y <- (mat_results$nodes_susc - mat_results$cum_nontarget_lost) / mat_results$nodes_susc


# ADJUST X AXIS
target_remove <- length(mat_basal_given)
mat_results$prop_target_removed <- mat_results$num_removed_tot / target_remove



# we want to get the robustness of the system for ONLY the sequence we're interested in, not the full removal
mat_results <- subset(mat_results, basal_removed_each %in% mat_basal_given | prop_remain == 1)
write.csv(mat_results, "CSM_MatOutput_SuppSppHL_Xtarget.csv")


root <- ggplot(mat_results, aes(x=prop_target_removed, y=Y))
(root + geom_point(shape=19) 
  + geom_line() + theme_bw(base_size = 14) 
  + xlim(0,1) + ylim(0,1)
  + xlab("Proportion of target species removed")
  + ylab("Proportion of susceptible species remaining")
)

# calc robustness (AUC)
robust_auc(x = mat_results$prop_target_removed, y = mat_results$Y) 


## include ES
mat_results_ES <- collapse_wrap_given(N = mat.ES, basal = mat_basal_given) # run ext seq like normal

# we want to get the robustness of the system for ONLY the sequence we're interested in, not the full removal
mat_results_ES <- subset(mat_results_ES, basal_removed_each %in% mat_basal_given | prop_remain == 1)

# make a column to track number of services lost at each step
mat_results_ES$service_lost <- ifelse(str_count(string = mat_results_ES$name_lost, pattern = "350|450|550|650|750|850|950"), str_count(mat_results_ES$name_lost,".50") , NA)

# track number services remain to use in propES_remain column
mat_results_ES[["service_lost"]][is.na(mat_results_ES[["service_lost"]])] <- 0
mat_results_ES$cum_nontarget_lost <- cumsum(mat_results_ES$service_lost)
mat_results_ES$ES_remain <- c(7-cumsum(mat_results_ES$service_lost))
mat_results_ES$propES_remain <- c(mat_results_ES$ES_remain/7)

# ADJUST X AXIS
target_remove <- length(mat_basal_given)
mat_results_ES$prop_removed <- mat_results_ES$num_removed_tot / target_remove

# revise output to get ES AUC
mat_results_ES %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "350|450|550|650|750|850|950")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / length(service_nodes))) %>%
  auc_wrapper_given()

root <- ggplot(mat_results_ES, aes(x=prop_removed, y=propES_remain))
(root + geom_point(shape=19) 
  + geom_line() + theme_bw(base_size = 14) 
  + xlim(0,1) + ylim(0,1)
  + xlab("Proportion of target species removed")
  + ylab("Proportion of ecosystem services remaining")
)

write.csv(mat_results_ES, "CSM_MatOutput_SuppSppHL_ES_Xtarget.csv")

## LH Supporting  --------------------------------------------------------------------------
## set sequence
mat_basal_given <-  rev(indirect.seq) 
mat_basal_given_total <- c(mat_basal_given[which(mat_basal_given %in%
                                                   as.numeric(colnames(mat.spp)))],# only the given basal species that ARE in the matrix
                           as.numeric(colnames(mat.spp)[which(!colnames(mat.spp) %in% mat_basal_given)]))# and also the rest of the species in the matrix that aren't in the basal list

## exclude ES, food web only
mat_results <- collapse_wrap_given(N = mat.spp, basal = mat_basal_given_total)


# ADJUST Y AXIS
mat_results$nodes_susc <- tot.susc
mat_results$nontarget_lost_each <- stri_count_fixed(mat_results$name_lost,";") # number species LOST each step
mat_results$cum_nontarget_lost <- cumsum(mat_results$nontarget_lost_each)

mat_results$Y <- (mat_results$nodes_susc - mat_results$cum_nontarget_lost) / mat_results$nodes_susc


# ADJUST X AXIS
target_remove <- length(mat_basal_given)
mat_results$prop_target_removed <- mat_results$num_removed_tot / target_remove





# we want to get the robustness of the system for ONLY the sequence we're interested in, not the full removal
mat_results <- subset(mat_results, basal_removed_each %in% mat_basal_given | prop_remain == 1)
write.csv(mat_results, "CSM_MatOutput_SuppSppLH_Xtarget.csv")


root <- ggplot(mat_results, aes(x=prop_target_removed, y=Y))
(root + geom_point(shape=19) 
  + geom_line() + theme_bw(base_size = 14) 
  + xlim(0,1) + ylim(0,1)
  + xlab("Proportion of target species removed")
  + ylab("Proportion of susceptible species remaining")
)

# calc robustness (AUC)
robust_auc(x = mat_results$prop_target_removed, y = mat_results$Y) 



## include ES
mat_results_ES <- collapse_wrap_given(N = mat.ES, basal = mat_basal_given) # run ext seq like normal

# we want to get the robustness of the system for ONLY the sequence we're interested in, not the full removal
mat_results_ES <- subset(mat_results_ES, basal_removed_each %in% mat_basal_given | prop_remain == 1)

# make a column to track number of services lost at each step
mat_results_ES$service_lost <- ifelse(str_count(string = mat_results_ES$name_lost, pattern = "350|450|550|650|750|850|950"), str_count(mat_results_ES$name_lost,".50") , NA)

# track number services remain to use in propES_remain column
mat_results_ES[["service_lost"]][is.na(mat_results_ES[["service_lost"]])] <- 0
mat_results_ES$cum_nontarget_lost <- cumsum(mat_results_ES$service_lost)
mat_results_ES$ES_remain <- c(7-cumsum(mat_results_ES$service_lost))
mat_results_ES$propES_remain <- c(mat_results_ES$ES_remain/7)

# ADJUST X AXIS
target_remove <- length(mat_basal_given)
mat_results_ES$prop_removed <- mat_results_ES$num_removed_tot / target_remove

# revise output to get ES AUC
mat_results_ES %>%
  # add new column that tracks whenever a service node is lost
  # input the service node names sep by a vertical bar: "|" ("or")
  ## NOTE: cannot have any spaces between the species names and the |
  mutate(service_lost = str_count(string = name_lost, pattern = "350|450|550|650|750|850|950")) %>%
  # add column with the cumulative sum of the services lost
  mutate(service_lost_c = cumsum(service_lost)) %>%
  # overwrite original column of proportion species remaining, for now
  mutate(prop_remain = 1 - (service_lost_c / length(service_nodes))) %>%
  auc_wrapper_given()

root <- ggplot(mat_results_ES, aes(x=prop_removed, y=propES_remain))
(root + geom_point(shape=19) 
  + geom_line() + theme_bw(base_size = 14) 
  + xlim(0,1) + ylim(0,1)
  + xlab("Proportion of target species removed")
  + ylab("Proportion of ecosystem services remaining")
)

write.csv(mat_results_ES, "CSM_MatOutput_SuppSppLH_ES_Xtarget.csv")


## All Random --------------------------------------------------------------------------
## This one runs until ALL species are lost, so no need to adjust prop_removed

mat_basal_given <- colnames(mat.spp)

# create empty data frame to fill with robustness values
rand <- data.frame(matrix(ncol=3,nrow=1000))
x <- c("Randomization","R_web","R_ES")
colnames(rand) <- x
rand$Randomization <- seq(1:1000)

# create empty list and fill it with 1000 randomizations of the direct esps
n <- c(1:1000)
sequences <- vector(mode="list",length=1000)
for(i in n){
  sequences[[i]] <- sample(mat_basal_given)
}


## * FOOD WEB ----

# create empty list to store robustness output for 1000 randomizations
mat_output_list <- vector(mode="list",length=1000)
for(i in n) {
  mat_output_list[[i]] <- collapse_wrap_given(N = mat.spp, basal = sequences[[i]])
}

mat_results_list <- mat_output_list

## need to add columns to calc new REMAIN proportions
for( i in seq_along(mat_results_list)){
  mat_results_list[[i]]$nodes_susc <- tot.susc
}


for( i in seq_along(mat_results_list)){
  mat_results_list[[i]]$nontarget_lost_each <- stri_count_fixed(mat_results_list[[i]]$name_lost,";") #count num species LOST
}


for( i in seq_along(mat_results_list)){
  mat_results_list[[i]]$cum_nontarget_lost <- cumsum(mat_results_list[[i]]$nontarget_lost_each) # calculate the cumulative prop nontarget LOST
}

for( i in seq_along(mat_results_list)){
  mat_results_list[[i]]$Y <- (mat_results_list[[i]]$nodes_susc - mat_results_list[[i]]$cum_nontarget_lost)/mat_results_list[[i]]$nodes_susc
}

##### calculate and store robustness for each of the randomizations
for (i in n){
  rand$R_web[i] <- robust_auc(x = mat_results_list[[i]]$prop_removed, y = mat_results_list[[i]]$Y)
}

min(rand$R_web)
mean(rand$R_web)
max(rand$R_web)

## * ES ----
# create empty list to store robustness output for 1000 randomizations
mat_output_list_ES <- vector(mode="list",length=1000)
for(i in n) {
  mat_output_list_ES[[i]] <- collapse_wrap_given(N = mat.ES, basal = sequences[[i]])
}



# calculate and store AUC for each of the randomizations
for(i in n){
  rand$ES_AUC[i] <- mat_output_list_ES[[i]] %>%
    # add new column that tracks whenever a service node is lost
    # input the service node names sep by a vertical bar: "|" ("or")
    ## NOTE: cannot have any spaces between the species names and the |
    mutate(service_lost = str_count(string = name_lost, pattern = "350|450|550|650|750|850|950")) %>%
    # add column with the cumulative sum of the services lost
    mutate(service_lost_c = cumsum(service_lost)) %>%
    # overwrite original column of proportion species remaining, for now
    mutate(prop_remain = 1 - (service_lost_c / length(service_nodes))) %>%
    auc_wrapper_given()
}


# unlist AUC then calculate robustness and store
rand$R_ES <- unlist(rand$ES_AUC)
rand <- rand[,-4]

min(rand$R_ES)
mean(rand$R_ES)
max(rand$R_ES)

write.csv(rand, "RobustValues_1000RAND_CSM.csv")

##### DIRECT RAND ----
# for this one, you will need to RErun the 100 robustness simulations on the networks

# set sequence
mat_basal_given <- direct.rand.seq

# create empty data frame to fill with robustness values
direct.rand <- data.frame(matrix(ncol=3,nrow=100))
x <- c("Randomization","R_web","R_ES")
colnames(direct.rand) <- x
direct.rand$Randomization <- seq(1:100)

# create empty list and fill it with 100 randomizations of the direct esps
n <- c(1:100)
sequences <- vector(mode="list",length=100)
for(i in n){
  sequences[[i]] <- sample(direct.rand.seq)
}

# create new empty list and fill it with the 100 randomized vectors + extra spp that we'll remove later. Need full web to run the robustness! 
sequences.total <- vector(mode="list",length=100)
for(i in n){
  sequences.total[[i]] <- c(sequences[[i]][which(mat_basal_given %in%
                                                   as.numeric(colnames(mat.spp)))],# only the given basal species that ARE in the matrix
                            as.numeric(colnames(mat.spp)[which(!colnames(mat.spp) %in% mat_basal_given)]))# and also the rest of the species in the matrix that aren't in the basal list
}

## * FOOD WEB ----
# create empty list to store robustness output for 100 randomizations
mat_output_list <- vector(mode="list",length=100)
for(i in n) {
  mat_output_list[[i]] <- collapse_wrap_given(N = mat.spp, basal = sequences.total[[i]])
}

# select only the rows that are in the original sequence (only direct ESPs!)
mat_results_list <- vector(mode="list",length=100)
for(i in n){
  mat_results_list[[i]] <- subset(mat_output_list[[i]], basal_removed_each %in% mat_basal_given | prop_remain == 1)
}


## need to ADJUST prop_removed to be prop_target_removed for EACH of the 100 randomizations
target_remove <- length(direct.rand.seq)

for( i in seq_along(mat_results_list)){
  mat_results_list[[i]]$prop_target_removed <- mat_results_list[[i]]$num_removed_tot / target_remove
}

## need to add columns to calc new REMAIN proportions
for( i in seq_along(mat_results_list)){
  mat_results_list[[i]]$nodes_susc <- tot.susc
}


for( i in seq_along(mat_results_list)){
  mat_results_list[[i]]$nontarget_lost_each <- stri_count_fixed(mat_results_list[[i]]$name_lost,";") #count num species LOST
}


for( i in seq_along(mat_results_list)){
  mat_results_list[[i]]$cum_nontarget_lost <- cumsum(mat_results_list[[i]]$nontarget_lost_each) # calculate the cumulative prop nontarget LOST
}

for( i in seq_along(mat_results_list)){
  mat_results_list[[i]]$Y <- (mat_results_list[[i]]$nodes_susc - mat_results_list[[i]]$cum_nontarget_lost)/mat_results_list[[i]]$nodes_susc
}

##### calculate and store robustness for each of the randomizations
for (i in n){
  direct.rand$R_web[i] <- robust_auc(x = mat_results_list[[i]]$prop_target_removed, y = mat_results_list[[i]]$Y)
}

# CALCULATE THE MEAN AUC TO USE IN ANALYSIS
mean(direct.rand$R_web)


## * ES ----
# create empty list to store robustness output for 100 randomizations
mat_output_list_ES <- vector(mode="list",length=100)
for(i in n) {
  mat_output_list_ES[[i]] <- collapse_wrap_given(N = mat.ES, basal = sequences.total[[i]])
}

# select only the rows that are in the original sequence (only direct ESPs!)
mat_results_list_ES <- vector(mode="list",length=100)
for(i in n){
  mat_results_list_ES[[i]] <- subset(mat_output_list_ES[[i]], basal_removed_each %in% mat_basal_given | prop_remain == 1)
}

for( i in seq_along(mat_results_list_ES)){
  mat_results_list_ES[[i]]$prop_removed <- mat_results_list_ES[[i]]$num_removed_tot / target_remove
}

# calculate and store AUC for each of the randomizations
for(i in n){
  direct.rand$ES_AUC[i] <- mat_results_list_ES[[i]] %>%
    # add new column that tracks whenever a service node is lost
    # input the service node names sep by a vertical bar: "|" ("or")
    ## NOTE: cannot have any spaces between the species names and the |
    mutate(service_lost = str_count(string = name_lost, pattern = "350|450|550|650|750|850|950")) %>%
    # add column with the cumulative sum of the services lost
    mutate(service_lost_c = cumsum(service_lost)) %>%
    # overwrite original column of proportion species remaining, for now
    mutate(prop_remain = 1 - (service_lost_c / length(service_nodes))) %>%
    auc_wrapper_given()
}

# unlist AUC then calculate robustness and store
direct.rand$R_ES <- unlist(direct.rand$ES_AUC)


# calculate mean robustness 
mean(direct.rand$R_ES)

# not really sure how you'd write a list of 100 dataframes of 9 columns each into a csv...so skipped it on this one. 
direct.rand <- direct.rand[,-4] 
write.csv(direct.rand, "DirectRandomVALUES_CSM.csv")



##### SUPP RANDOM ---------------------------------------------------------------------------
mat_basal_given <- indirect.seq

# create empty data frame to fill with robustness values
indirect.rand <- data.frame(matrix(ncol=3,nrow=100))
x <- c("Randomization","R_web","R_ES")
colnames(indirect.rand) <- x
indirect.rand$Randomization <- seq(1:100)

# create empty list and fill it with 100 randomizations of the direct esps
n <- c(1:100)
sequences <- vector(mode="list",length=100)
for(i in n){
  sequences[[i]] <- sample(indirect.seq)
}

# create new empty list and fill it with the 100 randomized vectors + extra spp that we'll remove later. Need full web to run the robustness! 
sequences.total <- vector(mode="list",length=100)
for(i in n){
  sequences.total[[i]] <- c(sequences[[i]][which(mat_basal_given %in%
                                                   as.numeric(colnames(mat.spp)))],# only the given basal species that ARE in the matrix
                            as.numeric(colnames(mat.spp)[which(!colnames(mat.spp) %in% mat_basal_given)]))# and also the rest of the species in the matrix that aren't in the basal list
}

## * FOOD WEB ----
# create empty list to store robustness output for 100 randomizations
mat_output_list <- vector(mode="list",length=100)
for(i in n) {
  mat_output_list[[i]] <- collapse_wrap_given(N = mat.spp, basal = sequences.total[[i]])
}

# select only the rows that are in the original sequence (only direct ESPs!)
mat_results_list <- vector(mode="list",length=100)
for(i in n){
  mat_results_list[[i]] <- subset(mat_output_list[[i]], basal_removed_each %in% mat_basal_given | prop_remain == 1)
}

## need to ADJUST prop_removed to be prop_target_removed for EACH of the 100 randomizations
target_remove <- length(indirect.seq)

for( i in seq_along(mat_results_list)){
  mat_results_list[[i]]$prop_target_removed <- mat_results_list[[i]]$num_removed_tot / target_remove
}

## need to add columns to calc new REMAIN proportions
for( i in seq_along(mat_results_list)){
  mat_results_list[[i]]$nodes_susc <- tot.susc
}


for( i in seq_along(mat_results_list)){
  mat_results_list[[i]]$nontarget_lost_each <- stri_count_fixed(mat_results_list[[i]]$name_lost,";") #count num species LOST
}


for( i in seq_along(mat_results_list)){
  mat_results_list[[i]]$cum_nontarget_lost <- cumsum(mat_results_list[[i]]$nontarget_lost_each) # calculate the cumulative prop nontarget LOST
}

for( i in seq_along(mat_results_list)){
  mat_results_list[[i]]$Y <- (mat_results_list[[i]]$nodes_susc - mat_results_list[[i]]$cum_nontarget_lost)/mat_results_list[[i]]$nodes_susc
}


##### calculate and store robustness for each of the randomizations
for (i in n){
  indirect.rand$R_web[i] <- robust_auc(x = mat_results_list[[i]]$prop_target_removed, y = mat_results_list[[i]]$Y)
}

mean(indirect.rand$R_web)


## * ES ----
# create empty list to store robustness output for 100 randomizations
mat_output_list_ES <- vector(mode="list",length=100)
for(i in n) {
  mat_output_list_ES[[i]] <- collapse_wrap_given(N = mat.ES, basal = sequences.total[[i]])
}

# select only the rows that are in the original sequence (only direct ESPs!)
mat_results_list_ES <- vector(mode="list",length=100)
for(i in n){
  mat_results_list_ES[[i]] <- subset(mat_output_list_ES[[i]], basal_removed_each %in% mat_basal_given | prop_remain == 1)
}

for( i in seq_along(mat_results_list_ES)){
  mat_results_list_ES[[i]]$prop_removed <- mat_results_list_ES[[i]]$num_removed_tot / target_remove
}

# calculate and store AUC for each of the randomizations
for(i in n){
  indirect.rand$ES_AUC[i] <- mat_results_list_ES[[i]] %>%
    # add new column that tracks whenever a service node is lost
    # input the service node names sep by a vertical bar: "|" ("or")
    ## NOTE: cannot have any spaces between the species names and the |
    mutate(service_lost = str_count(string = name_lost, pattern = "350|450|550|650|750|850|950")) %>%
    # add column with the cumulative sum of the services lost
    mutate(service_lost_c = cumsum(service_lost)) %>%
    # overwrite original column of proportion species remaining, for now
    mutate(prop_remain = 1 - (service_lost_c / length(service_nodes))) %>%
    auc_wrapper_given()
}


# unlist AUC then calculate robustness and store
indirect.rand$R_ES <- unlist(indirect.rand$ES_AUC)
indirect.rand <- indirect.rand[,-4]

# calculate mean robustness and store
mean(indirect.rand$R_ES)

write.csv(indirect.rand, "IndirectRandomVALUES_CSM.csv")






















