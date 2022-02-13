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

entire_uncorrelated_relaxed_tree <- entire_tree

for (i in 1:length(entire_uncorrelated_relaxed_tree$edge.length)){
entire_uncorrelated_relaxed_tree$edge.length[[i]] <- entire_uncorrelated_relaxed_tree$edge.length[[i]]*rlnorm(1, log(0.05) - ((0.2937025^2)/2), sdlog=0.2937025)
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

entire_div_shift_uncorrelated_relaxed_tree <- entire_div_shift_tree

for (i in 1:length(entire_div_shift_uncorrelated_relaxed_tree$edge.length)){
entire_div_shift_uncorrelated_relaxed_tree$edge.length[[i]] <- entire_div_shift_uncorrelated_relaxed_tree$edge.length[[i]]*rlnorm(1, log(0.05) - ((0.2937025^2)/2), sdlog=0.2937025)
}

################################################################
################################################################
###SIMULATE_SEQUENCES###########################################
################################################################
################################################################

################################################################
###STANDARD_TREE_SEQUENCES######################################
################################################################

sequences_uncorrelated_relaxed <- simSeq(entire_uncorrelated_relaxed_tree, l = 50000, type = "DNA", rate = 1)

################################################################
###DIV_SHIFT_SEQUENCES##########################################
################################################################

sequences_div_shift_uncorrelated_relaxed <- simSeq(entire_div_shift_uncorrelated_relaxed_tree, l = 50000, type = "DNA", rate = 1)

###################################################################
###################################################################
###WRITE_OUTPUTS_TO_FILE###########################################
###################################################################
###################################################################

dir.create("sims_output")

###################################################################
###TREES###########################################################
###################################################################

write.tree(entire_tree, "sims_output/entire_tree.tre") # constant time tree with outgroup
write.tree(drop.tip(entire_tree, "t1001"), "sims_output/entire_tree_no_og.tre") # constant time tree without outgroup 
write.tree(entire_uncorrelated_relaxed_tree, "sims_output/entire_uncorrelated_relaxed_tree.tre") # constant tree with ucln
write.tree(entire_div_shift_tree, "sims_output/entire_div_shift_tree.tre") # div shift time tree with outgroup
write.tree(drop.tip(entire_div_shift_tree, "t1001"), "sims_output/entire_div_shift_tree_no_og.tre") #div shift time tree without outgroup 
write.tree(entire_div_shift_uncorrelated_relaxed_tree, "sims_output/entire_div_shift_uncorrelated_relaxed_tree.tre") # div shift tree with ucln

###################################################################
###SEQUENCES#######################################################
###################################################################

###CONSTANT_SEQUENCES##############################################

write.phyDat(sequences_uncorrelated_relaxed, file = "sims_output/uncorrelated_relaxed_fixed.nexus", format="nexus")

###DIV_SHIFT_SEQUENCES##############################################

write.phyDat(sequences_div_shift_uncorrelated_relaxed, file = "sims_output/uncorrelated_relaxed_div_shift.nexus", format="nexus")
