library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(Rmisc)
library(TreePar)

tree_size <- 24

############################################################
############################################################
###SIMULATE_TREES###########################################
############################################################
############################################################

############################################################
###CONSTANT_DIV_RATE_TREES##################################
############################################################

###TIMETREE#################################################

entire_tree <- sim.bd.taxa(tree_size, 1, lambda = 3, mu = 0, frac = 1, complete = TRUE, stochsampling = FALSE)
entire_tree <- entire_tree[[1]]
outgroup <- list(edge=matrix(c(2,1),1,2), tip.label="t1001", edge.length=max(node.depth.edgelength(entire_tree)), Nnode=1)
class(outgroup) <- "phylo"
entire_tree <- bind.tree(entire_tree, outgroup, position=0)

###MOLECULAR_TREES##########################################

entire_clock_tree <- entire_tree
entire_uncorrelated_relaxed_tree <- entire_tree

entire_clock_tree$edge.length <- entire_clock_tree$edge.length*0.05
entire_random_relaxed_tree <- entire_clock_tree

for (i in 1:length(entire_uncorrelated_relaxed_tree$edge.length)){
entire_uncorrelated_relaxed_tree$edge.length[[i]] <- entire_uncorrelated_relaxed_tree$edge.length[[i]]*rlnorm(1, log(0.05) - ((0.2937025^2)/2), sdlog=0.2937025)
}

for (i in 1:length(entire_random_relaxed_tree$edge.length)){
jump <- sample(seq(1, 10, 1), 1)
print(jump)
if (jump > 9){
multiplier <- rgamma(1,2,2)
entire_random_relaxed_tree$edge.length[[i]] <- entire_random_relaxed_tree$edge.length[[i]]*multiplier
if (entire_random_relaxed_tree[[1]][,2][[i]] > length(entire_random_relaxed_tree$tip.label)){
descendant_nodes <- getDescendants(entire_random_relaxed_tree, i)
for (b in 1:nrow(entire_random_relaxed_tree[[1]])){
for (c in 1:length(descendant_nodes)){
if (entire_random_relaxed_tree[[1]][,2][[b]] %in% descendant_nodes){
entire_random_relaxed_tree$edge.length[[b]] <- entire_random_relaxed_tree$edge.length[[b]]*multiplier
}
}
}
}
}
}

############################################################
###VARIABLE_DIV_RATE_TREES##################################
############################################################

###TIMETREE#################################################

continue <- 1
while (continue == 1){
entire_div_shift_tree <- sim.bd.taxa(tree_size, 1, lambda = 3, mu = 0, frac = 1, complete = TRUE, stochsampling = FALSE)
entire_div_shift_tree <- entire_div_shift_tree[[1]] 
if ((length(extract.clade(entire_div_shift_tree, length(entire_div_shift_tree$tip.label)+2)$tip.label) >= 18)
&
(length(extract.clade(entire_div_shift_tree, length(entire_div_shift_tree$tip.label)+2)$tip.label) < 22)){
continue <- 0
}
}


retainer_tree <- entire_div_shift_tree

descendant_nodes <- getDescendants(entire_div_shift_tree, length(entire_div_shift_tree$tip.label)+2)
for (j in 1:length(entire_div_shift_tree$edge.length)){
for (i in 1:length(descendant_nodes)){
if (entire_div_shift_tree[[1]][,2][[j]] == descendant_nodes[[i]]){
entire_div_shift_tree$edge.length[[j]] <- entire_div_shift_tree$edge.length[[j]] * 0.2
}
}
}

entire_div_shift_tree$edge.length[[which(entire_div_shift_tree[[1]][,2] == length(entire_div_shift_tree$tip.label)+2)]] <- 
entire_div_shift_tree$edge.length[[which(entire_div_shift_tree[[1]][,2] == length(entire_div_shift_tree$tip.label)+2)]] + 
(max(node.depth.edgelength(entire_div_shift_tree)) - node.depth.edgelength(entire_div_shift_tree)[min(descendant_nodes)])

outgroup <- list(edge=matrix(c(2,1),1,2), tip.label="t1001", edge.length=max(node.depth.edgelength(entire_div_shift_tree)), Nnode=1)
class(outgroup) <- "phylo"
entire_div_shift_tree <- bind.tree(entire_div_shift_tree, outgroup, position=0)

###MOLECULAR_TREES##############################################

entire_div_shift_clock_tree <- entire_div_shift_tree
entire_div_shift_uncorrelated_relaxed_tree <- entire_div_shift_tree

entire_div_shift_clock_tree$edge.length <- entire_div_shift_clock_tree$edge.length*0.05
entire_div_shift_random_relaxed_tree <- entire_div_shift_clock_tree

for (i in 1:length(entire_div_shift_uncorrelated_relaxed_tree$edge.length)){
entire_div_shift_uncorrelated_relaxed_tree$edge.length[[i]] <- entire_div_shift_uncorrelated_relaxed_tree$edge.length[[i]]*rlnorm(1, log(0.05) - ((0.2937025^2)/2), sdlog=0.2937025)
}

for (i in 1:length(entire_div_shift_random_relaxed_tree$edge.length)){
jump <- sample(seq(1, 10, 1), 1)
print(jump)
if (jump > 9){
multiplier <- rgamma(1, 2, 2)
entire_div_shift_random_relaxed_tree$edge.length[[i]] <- entire_div_shift_random_relaxed_tree$edge.length[[i]]*multiplier
if (entire_div_shift_random_relaxed_tree[[1]][,2][[i]] > length(entire_div_shift_random_relaxed_tree$tip.label)){
descendant_nodes <- getDescendants(entire_div_shift_random_relaxed_tree, i)
for (b in 1:nrow(entire_div_shift_random_relaxed_tree[[1]])){
for (c in 1:length(descendant_nodes)){
if (entire_div_shift_random_relaxed_tree[[1]][,2][[b]] %in% descendant_nodes){
entire_div_shift_random_relaxed_tree$edge.length[[b]] <- entire_div_shift_random_relaxed_tree$edge.length[[b]]*multiplier
}
}
}
}
}
}




################################################################
################################################################
###SIMULATE_SEQUENCES###########################################
################################################################
################################################################

################################################################
###STANDARD_TREE_SEQUENCES######################################
################################################################

sequences_clock <- simSeq(entire_clock_tree, l = 50000, type = "DNA", rate = 1)
sequences_uncorrelated_relaxed <- simSeq(entire_uncorrelated_relaxed_tree, l = 50000, type = "DNA", rate = 1)
sequences_random_relaxed <- simSeq(entire_random_relaxed_tree, l = 50000, type = "DNA", rate = 1)

################################################################
###DIV_SHIFT_SEQUENCES##########################################
################################################################

sequences_div_shift_clock <- simSeq(entire_div_shift_clock_tree, l = 50000, type = "DNA", rate = 1)
sequences_div_shift_uncorrelated_relaxed <- simSeq(entire_div_shift_uncorrelated_relaxed_tree, l = 50000, type = "DNA", rate = 1)
sequences_div_shift_random_relaxed <- simSeq(entire_div_shift_random_relaxed_tree, l = 50000, type = "DNA", rate = 1)

################################################################
###STANDARD_SUBSET_SEQUENCES####################################
################################################################

sequences_clock_small <- subset(sequences_clock, select = seq(1, 1000, 1), site.pattern = FALSE)
sequences_uncorrelated_relaxed_small <- subset(sequences_uncorrelated_relaxed, select = seq(1, 1000, 1), site.pattern = FALSE)
sequences_random_relaxed_small <- subset(sequences_random_relaxed, select = seq(1, 1000, 1), site.pattern = FALSE)

################################################################
###DIV_SHIFT_SUBSET_SEQUENCES###################################
################################################################

sequences_div_shift_clock_small <- subset(sequences_div_shift_clock, select = seq(1, 1000, 1), site.pattern = FALSE)
sequences_div_shift_uncorrelated_relaxed_small <- subset(sequences_div_shift_uncorrelated_relaxed, select = seq(1, 1000, 1), site.pattern = FALSE)
sequences_div_shift_random_relaxed_small <- subset(sequences_div_shift_random_relaxed, select = seq(1, 1000, 1), site.pattern = FALSE)

################################################################
################################################################
###GET_STATS####################################################
################################################################
################################################################

################################################################
###FOR_KNOWN_ANALYSIS###########################################
################################################################

min_rates_clock_tree <- 0.05 
max_rates_clock_tree <- 0.05

min_rates_uncorrelated_relaxed_tree <- min(entire_uncorrelated_relaxed_tree$edge.length/entire_tree$edge.length)
max_rates_uncorrelated_relaxed_tree <- max(entire_uncorrelated_relaxed_tree$edge.length/entire_tree$edge.length)
mean_rates_uncorrelated_relaxed_tree <- mean(entire_uncorrelated_relaxed_tree$edge.length/entire_tree$edge.length)
log_mean_rates_uncorrelated_relaxed_tree <- mean(log(entire_uncorrelated_relaxed_tree$edge.length/entire_tree$edge.length))
log_sd_rates_uncorrelated_relaxed_tree <- sd(log(entire_uncorrelated_relaxed_tree$edge.length/entire_tree$edge.length))

min_rates_random_relaxed_tree <- min(entire_random_relaxed_tree$edge.length/entire_tree$edge.length)
max_rates_random_relaxed_tree <- max(entire_random_relaxed_tree$edge.length/entire_tree$edge.length)
mean_rates_random_relaxed_tree <- mean(entire_random_relaxed_tree$edge.length/entire_tree$edge.length)
log_mean_rates_random_relaxed_tree <- mean(log(entire_random_relaxed_tree$edge.length/entire_tree$edge.length))
log_sd_rates_random_relaxed_tree <- sd(log(entire_random_relaxed_tree$edge.length/entire_tree$edge.length))

min_rates_div_shift_clock_tree <- 0.05
max_rates_div_shift_clock_tree <- 0.05

min_rates_div_shift_uncorrelated_relaxed_tree <- min(entire_div_shift_uncorrelated_relaxed_tree$edge.length/entire_div_shift_tree$edge.length)
max_rates_div_shift_uncorrelated_relaxed_tree <- max(entire_div_shift_uncorrelated_relaxed_tree$edge.length/entire_div_shift_tree$edge.length)
mean_rates_div_shift_uncorrelated_relaxed_tree <- mean(entire_div_shift_uncorrelated_relaxed_tree$edge.length/entire_div_shift_tree$edge.length)
log_mean_rates_div_shift_uncorrelated_relaxed_tree <- mean(log(entire_div_shift_uncorrelated_relaxed_tree$edge.length/entire_div_shift_tree$edge.length))
log_sd_rates_div_shift_uncorrelated_relaxed_tree <- sd(log(entire_div_shift_uncorrelated_relaxed_tree$edge.length/entire_div_shift_tree$edge.length))

min_rates_div_shift_random_relaxed_tree <- min(entire_div_shift_random_relaxed_tree$edge.length/entire_div_shift_tree$edge.length)
max_rates_div_shift_random_relaxed_tree <- max(entire_div_shift_random_relaxed_tree$edge.length/entire_div_shift_tree$edge.length)
mean_rates_div_shift_random_relaxed_tree <- mean(entire_div_shift_random_relaxed_tree$edge.length/entire_div_shift_tree$edge.length)
log_mean_rates_div_shift_random_relaxed_tree <- mean(log(entire_div_shift_random_relaxed_tree$edge.length/entire_div_shift_tree$edge.length))
log_sd_rates_div_shift_random_relaxed_tree <- sd(log(entire_div_shift_random_relaxed_tree$edge.length/entire_div_shift_tree$edge.length))

###################################################################
###################################################################
###WRITE_OUTPUTS_TO_FILE###########################################
###################################################################
###################################################################

dir.create("sims_output")

###################################################################
###TREES###########################################################
###################################################################

write.tree(entire_tree, "sims_output/entire_tree.tre")
write.tree(drop.tip(entire_tree, "t1001"), "sims_output/entire_tree_no_og.tre") 
write.tree(entire_clock_tree, "sims_output/entire_clock_tree.tre")
write.tree(entire_uncorrelated_relaxed_tree, "sims_output/entire_uncorrelated_relaxed_tree.tre")
write.tree(entire_random_relaxed_tree, "sims_output/entire_random_relaxed_tree.tre")
write.tree(entire_div_shift_tree, "sims_output/entire_div_shift_tree.tre")
write.tree(drop.tip(entire_div_shift_tree, "t1001"), "sims_output/entire_div_shift_tree_no_og.tre") 
write.tree(entire_div_shift_clock_tree, "sims_output/entire_div_shift_clock_tree.tre")
write.tree(entire_div_shift_uncorrelated_relaxed_tree, "sims_output/entire_div_shift_uncorrelated_relaxed_tree.tre")
write.tree(entire_div_shift_random_relaxed_tree, "sims_output/entire_div_shift_random_relaxed_tree.tre")

###################################################################
###SEQUENCES#######################################################
###################################################################

###CONSTANT_SEQUENCES##############################################

write.phyDat(sequences_clock, file = "sims_output/clock_fixed.nexus", format="nexus")
write.phyDat(sequences_uncorrelated_relaxed, file = "sims_output/uncorrelated_relaxed_fixed.nexus", format="nexus")
write.phyDat(sequences_random_relaxed, file = "sims_output/random_relaxed_fixed.nexus", format="nexus")

write.phyDat(subset(sequences_clock, drop.tip(entire_tree, "t1001")$tip.label), file = "sims_output/clock_fixed_no_og.nexus", format="nexus")
write.phyDat(subset(sequences_uncorrelated_relaxed, drop.tip(entire_tree, "t1001")$tip.label), file = "sims_output/uncorrelated_relaxed_fixed_no_og.nexus", format="nexus")
write.phyDat(subset(sequences_random_relaxed, drop.tip(entire_tree, "t1001")$tip.label), file = "sims_output/random_relaxed_fixed_no_og.nexus", format="nexus")

###DIV_SHIFT_SEQUENCES##############################################

write.phyDat(sequences_div_shift_clock, file = "sims_output/clock_div_shift.nexus", format="nexus")
write.phyDat(sequences_div_shift_uncorrelated_relaxed, file = "sims_output/uncorrelated_relaxed_div_shift.nexus", format="nexus")
write.phyDat(sequences_div_shift_random_relaxed, file = "sims_output/random_relaxed_div_shift.nexus", format="nexus")

write.phyDat(subset(sequences_div_shift_clock, drop.tip(entire_tree, "t1001")$tip.label), file = "sims_output/clock_div_shift_no_og.nexus", format="nexus")
write.phyDat(subset(sequences_div_shift_uncorrelated_relaxed, drop.tip(entire_tree, "t1001")$tip.label), file = "sims_output/uncorrelated_relaxed_div_shift_no_og.nexus", format="nexus")
write.phyDat(subset(sequences_div_shift_random_relaxed, drop.tip(entire_tree, "t1001")$tip.label), file = "sims_output/random_relaxed_div_shift_no_og.nexus", format="nexus")

###CONSTANT_SUBSET_SEQUENCES########################################

write.phyDat(sequences_clock_small, file = "sims_output/clock_fixed_small.nexus", format="nexus")
write.phyDat(sequences_uncorrelated_relaxed_small, file = "sims_output/uncorrelated_relaxed_fixed_small.nexus", format="nexus")
write.phyDat(sequences_random_relaxed_small, file = "sims_output/random_relaxed_fixed_small.nexus", format="nexus")

write.phyDat(subset(sequences_clock_small, drop.tip(entire_tree, "t1001")$tip.label), file = "sims_output/clock_fixed_no_og_small.nexus", format="nexus")
write.phyDat(subset(sequences_uncorrelated_relaxed_small, drop.tip(entire_tree, "t1001")$tip.label), file = "sims_output/uncorrelated_relaxed_fixed_no_og_small.nexus", format="nexus")
write.phyDat(subset(sequences_random_relaxed_small, drop.tip(entire_tree, "t1001")$tip.label), file = "sims_output/random_relaxed_fixed_no_og_small.nexus", format="nexus")

###DIV_SHIFT_SUBSET_SEQUENCES#######################################

write.phyDat(sequences_div_shift_clock_small, file = "sims_output/clock_div_shift_small.nexus", format="nexus")
write.phyDat(sequences_div_shift_uncorrelated_relaxed_small, file = "sims_output/uncorrelated_relaxed_div_shift_small.nexus", format="nexus")
write.phyDat(sequences_div_shift_random_relaxed_small, file = "sims_output/random_relaxed_div_shift_small.nexus", format="nexus")

write.phyDat(subset(sequences_div_shift_clock_small, drop.tip(entire_tree, "t1001")$tip.label), file = "sims_output/clock_div_shift_no_og_small.nexus", format="nexus")
write.phyDat(subset(sequences_div_shift_uncorrelated_relaxed_small, drop.tip(entire_tree, "t1001")$tip.label), file = "sims_output/uncorrelated_relaxed_div_shift_no_og_small.nexus", format="nexus")
write.phyDat(subset(sequences_div_shift_random_relaxed_small, drop.tip(entire_tree, "t1001")$tip.label), file = "sims_output/random_relaxed_div_shift_no_og_small.nexus", format="nexus")

####################################################################
###STATS############################################################
####################################################################

fileConn <- file("sims_output/rate_info_known.txt")
writeLines(c(
paste("min_rates_clock_tree<-",min_rates_clock_tree, sep=""),
paste("max_rates_clock_tree<-",max_rates_clock_tree, sep=""), 
paste("min_rates_uncorrelated_relaxed_tree<-",min_rates_uncorrelated_relaxed_tree, sep=""), 
paste("max_rates_uncorrelated_relaxed_tree<-",max_rates_uncorrelated_relaxed_tree, sep=""),
paste("mean_rates_uncorrelated_relaxed_tree<-",mean_rates_uncorrelated_relaxed_tree, sep=""),
paste("log_mean_rates_uncorrelated_relaxed_tree<-",log_mean_rates_uncorrelated_relaxed_tree, sep=""),
paste("log_sd_rates_uncorrelated_relaxed_tree<-",log_sd_rates_uncorrelated_relaxed_tree, sep=""),
paste("min_rates_random_relaxed_tree<-",min_rates_random_relaxed_tree, sep=""), 
paste("max_rates_random_relaxed_tree<-",max_rates_random_relaxed_tree, sep=""),
paste("mean_rates_random_relaxed_tree<-",mean_rates_random_relaxed_tree, sep=""),
paste("log_mean_rates_random_relaxed_tree<-",log_mean_rates_random_relaxed_tree, sep=""),
paste("log_sd_rates_random_relaxed_tree<-",log_sd_rates_random_relaxed_tree, sep=""),
paste("min_rates_div_shift_clock_tree<-",min_rates_div_shift_clock_tree, sep=""),
paste("max_rates_div_shift_clock_tree<-",max_rates_div_shift_clock_tree, sep=""), 
paste("min_rates_div_shift_uncorrelated_relaxed_tree<-",min_rates_div_shift_uncorrelated_relaxed_tree, sep=""), 
paste("max_rates_div_shift_uncorrelated_relaxed_tree<-",max_rates_div_shift_uncorrelated_relaxed_tree, sep=""),
paste("mean_rates_div_shift_uncorrelated_relaxed_tree<-",mean_rates_div_shift_uncorrelated_relaxed_tree, sep=""),
paste("log_mean_rates_div_shift_uncorrelated_relaxed_tree<-",log_mean_rates_div_shift_uncorrelated_relaxed_tree, sep=""),
paste("log_sd_rates_div_shift_uncorrelated_relaxed_tree<-",log_sd_rates_div_shift_uncorrelated_relaxed_tree, sep=""),
paste("min_rates_div_shift_random_relaxed_tree<-",min_rates_div_shift_random_relaxed_tree, sep=""), 
paste("max_rates_div_shift_random_relaxed_tree<-",max_rates_div_shift_random_relaxed_tree, sep=""),
paste("mean_rates_div_shift_random_relaxed_tree<-",mean_rates_div_shift_random_relaxed_tree, sep=""),
paste("log_mean_rates_div_shift_random_relaxed_tree<-",log_mean_rates_div_shift_random_relaxed_tree, sep=""),
paste("log_sd_rates_div_shift_random_relaxed_tree<-",log_sd_rates_div_shift_random_relaxed_tree, sep="")
),fileConn)
close(fileConn)
