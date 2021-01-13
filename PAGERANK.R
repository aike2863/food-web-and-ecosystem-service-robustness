####################################################
## Pagerank for important spp ######################
# Written by: Aislyn Keyes #########################
# Last edited: Feb 25, 2020  #######################


#load packages
library(igraph)
library(dplyr)

# load and process data
nodes.ES <- read.csv("Final_EPB_Nodes_ES.csv") 
edges.ES <- read.csv("Final.EPB.edges.filtered.csv")

# want to flip direction of edges so that services are pointing to spp
# this allows us to do random walks from ES, specifying ES as important nodes
edges.flip <- edges.ES[,c(2,1,3)]  

# create graph object using edges.ES
net.ES <- graph.data.frame(edges.flip,
                           directed = T,
                           vertices = nodes.ES)
plot(net.ES) # check to make sure arrows are flipped

##########################################################################
##########################################################################
##########################################################################

#PAGERANK
# we want to run pagerank for each service individually
# this allows us to look at species importance for each service
# we can average the probabilities at the end to look at total importance

# we won't run a page rank for 350 or 450 (wave atten. and shoreline stab.)
# b/c there are no supporting spp. here, both services rely only on basal spp.
# for each of the services, we'll run a pagerank algorithm
# to account for spp directly providing, we will replace prob. with 0 or NA
# if they are directly linked to a service

# 550, carbon sequestration
carbon.seq <- data.frame(SpeciesID=nodes.ES$SpeciesID, prob=page_rank(graph=net.ES, damping = 0.85, directed = T, 
                  personalized = c(rep(0,times=nrow(nodes.ES)-nrow(nodes.ES[nodes.ES$NodeType=="ES",])),
                                   0,0,1,0,0,0,0))$vector)

carbon.seq.SPP <- carbon.seq[!(carbon.seq$SpeciesID>300),] # drop ES from ranking

carbon.seq.direct <- data.frame(SpeciesID = ifelse(edges.ES$ConsumerSpeciesID=="550",
                                                edges.ES$ResourceSpeciesID,NA)) # pull the SpeciesID for direct ESP
carbon.seq.direct <- c(na.omit(carbon.seq.direct)) # remove NAs 
carbon.seq.support <- carbon.seq.SPP[!(carbon.seq.SPP$SpeciesID %in% carbon.seq.direct$SpeciesID),] # keep only supporting spp! 

attach(carbon.seq.support)
carbon.seq.support <- carbon.seq.support[order(-prob),] # sort by probability (highest to lowest)
detach(carbon.seq.support)
write.csv(carbon.seq.support,"EPB_Indirect_CarbonSeq.csv")


# 650, water filtration
water.filt <- data.frame(SpeciesID=nodes.ES$SpeciesID, prob=page_rank(graph=net.ES, damping = 0.85, directed = T, 
                        personalized = c(rep(0,times=nrow(nodes.ES)-nrow(nodes.ES[nodes.ES$NodeType=="ES",])),
                                         0,0,0,1,0,0,0))$vector)

water.filt.SPP <- water.filt[!(water.filt$SpeciesID>300),] # drop ES from ranking

water.filt.direct <- data.frame(SpeciesID = ifelse(edges.ES$ConsumerSpeciesID=="650",
                                                edges.ES$ResourceSpeciesID,NA)) # pull the SpeciesID for direct ESP
water.filt.direct <- c(na.omit(water.filt.direct)) # remove NAs 
water.filt.support <- water.filt.SPP[!(water.filt.SPP$SpeciesID %in% water.filt.direct$SpeciesID),] # keep only supporting spp! 

attach(water.filt.support)
water.filt.support <- water.filt.support[order(-prob),] # sort by probability (highest to lowest)
detach(water.filt.support)
write.csv(water.filt.support,"EPB_Indirect_Waterfilt.csv")


# 750, Fishery
fishery <- data.frame(SpeciesID=nodes.ES$SpeciesID, prob=page_rank(graph=net.ES, damping = 0.85, directed = T, 
                        personalized = c(rep(0,times=nrow(nodes.ES)-nrow(nodes.ES[nodes.ES$NodeType=="ES",])),
                                         0,0,0,0,1,0,0))$vector)

fishery.SPP <- fishery[!(fishery$SpeciesID>300),] # drop ES from ranking

fishery.direct <- data.frame(SpeciesID = ifelse(edges.ES$ConsumerSpeciesID=="750",
                                                edges.ES$ResourceSpeciesID,NA)) # pull the SpeciesID for direct ESP
fishery.direct <- c(na.omit(fishery.direct)) # remove NAs 
fishery.support <- fishery.SPP[!(fishery.SPP$SpeciesID %in% fishery.direct$SpeciesID),] # keep only supporting spp! 

attach(fishery.support)
fishery.support <- fishery.support[order(-prob),] # sort by probability (highest to lowest)
detach(fishery.support)
write.csv(fishery.support,"EPB_Indirect_Fishery.csv")

# 850, bird watching
bird.watch <- data.frame(SpeciesID=nodes.ES$SpeciesID, prob=page_rank(graph=net.ES, damping = 0.85, directed = T, 
                        personalized = c(rep(0,times=nrow(nodes.ES)-nrow(nodes.ES[nodes.ES$NodeType=="ES",])),
                                         0,0,0,0,0,1,0))$vector)

birdwatch.SPP <- bird.watch[!(bird.watch$SpeciesID>300),] # drop ES from ranking

birdwatch.direct <- data.frame(SpeciesID = ifelse(edges.ES$ConsumerSpeciesID=="850",
                                                edges.ES$ResourceSpeciesID,NA)) # pull the SpeciesID for direct ESP
birdwatch.direct <- c(na.omit(birdwatch.direct)) # remove NAs 
birdwatch.support <- birdwatch.SPP[!(birdwatch.SPP$SpeciesID %in% birdwatch.direct$SpeciesID),] # keep only supporting spp! 

attach(birdwatch.support)
birdwatch.support <- birdwatch.support[order(-prob),] # sort by probability (highest to lowest)
detach(birdwatch.support)
write.csv(birdwatch.support,"EPB_Indirect_Birdwatch.csv")


# 950, waterfowl hunting
wat.hunt <- data.frame(SpeciesID=nodes.ES$SpeciesID, prob=page_rank(graph=net.ES, damping = 0.85, directed = T, 
                        personalized = c(rep(0,times=nrow(nodes.ES)-nrow(nodes.ES[nodes.ES$NodeType=="ES",])),
                                         0,0,0,0,0,0,1))$vector)

wathunt.SPP <- wat.hunt[!(wat.hunt$SpeciesID>300),] # drop ES from ranking

wathunt.direct <- data.frame(SpeciesID = ifelse(edges.ES$ConsumerSpeciesID=="950",
                                                edges.ES$ResourceSpeciesID,NA)) # pull the SpeciesID for direct ESP
wathunt.direct <- c(na.omit(wathunt.direct)) # remove NAs 
wathunt.support <- wathunt.SPP[!(wathunt.SPP$SpeciesID %in% wathunt.direct$SpeciesID),] # keep only supporting spp, WE DON'T WANT DIRECT ESP 

attach(wathunt.support)
wathunt.support <- wathunt.support[order(-prob),] # sort by probability (highest to lowest)
detach(wathunt.support)
write.csv(wathunt.support,"EPB_Indirect_WaterfowlHunt.csv")


##########################################################################
##########################################################################
##########################################################################


# we have already exported .csv file for each service individually
# we want to get an mean prob for each species across all services for the main ES robustness

all.support <- rbind(water.filt.support,wathunt.support,carbon.seq.support,fishery.support,birdwatch.support)
all.support$SpeciesID <- as.factor(all.support$SpeciesID)
mean.support <- aggregate(prob ~ SpeciesID, all.support, mean) # calculate mean prob. for each spp.

attach(mean.support)
mean.support <- mean.support[order(-prob),]
detach(mean.support)

write.csv(mean.support, "EPB_IndirectAll.csv")








