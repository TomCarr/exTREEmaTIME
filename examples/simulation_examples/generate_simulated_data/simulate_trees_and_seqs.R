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

tree <- sim.bd.taxa(tree_size, 1, lambda = 3, mu = 0, frac = 1, complete = TRUE, stochsampling = FALSE)
tree <- tree[[1]]
outgroup <- list(edge=matrix(c(2,1),1,2), tip.label="t1001", edge.length=max(node.depth.edgelength(tree)), Nnode=1)
class(outgroup) <- "phylo"
tree <- bind.tree(tree, outgroup, position=0)

###MOLECULAR_TREES##########################################

###UNCORRELATED_RELAXED

uncorrelated_relaxed_tree <- tree

for (i in 1:length(uncorrelated_relaxed_tree$edge.length)){
uncorrelated_relaxed_tree$edge.length[[i]] <- uncorrelated_relaxed_tree$edge.length[[i]]*rlnorm(1, log(0.05) - ((0.2937025^2)/2), sdlog=0.2937025)
}

############################################################
###VARIABLE_DIV_RATE_TREES##################################
############################################################

###TIMETREE#################################################

continue <- 1
while (continue == 1){
div_shift_tree <- sim.bd.taxa(tree_size, 1, lambda = 3, mu = 0, frac = 1, complete = TRUE, stochsampling = FALSE)
div_shift_tree <- div_shift_tree[[1]] 
if ((length(extract.clade(div_shift_tree, length(div_shift_tree$tip.label)+2)$tip.label) >= 18)
&
(length(extract.clade(div_shift_tree, length(div_shift_tree$tip.label)+2)$tip.label) < 22)){
continue <- 0
}
}

retainer_tree <- div_shift_tree

descendant_nodes <- getDescendants(div_shift_tree, length(div_shift_tree$tip.label)+2)
for (j in 1:length(div_shift_tree$edge.length)){
for (i in 1:length(descendant_nodes)){
if (div_shift_tree[[1]][,2][[j]] == descendant_nodes[[i]]){
div_shift_tree$edge.length[[j]] <- div_shift_tree$edge.length[[j]] * 0.2
}
}
}

div_shift_tree$edge.length[[which(div_shift_tree[[1]][,2] == length(div_shift_tree$tip.label)+2)]] <- 
div_shift_tree$edge.length[[which(div_shift_tree[[1]][,2] == length(div_shift_tree$tip.label)+2)]] + 
(max(node.depth.edgelength(div_shift_tree)) - node.depth.edgelength(div_shift_tree)[min(descendant_nodes)])

outgroup <- list(edge=matrix(c(2,1),1,2), tip.label="t1001", edge.length=max(node.depth.edgelength(div_shift_tree)), Nnode=1)
class(outgroup) <- "phylo"
div_shift_tree <- bind.tree(div_shift_tree, outgroup, position=0)

###MOLECULAR_TREES##############################################

div_shift_uncorrelated_relaxed_tree <- div_shift_tree

for (i in 1:length(div_shift_uncorrelated_relaxed_tree$edge.length)){
div_shift_uncorrelated_relaxed_tree$edge.length[[i]] <- div_shift_uncorrelated_relaxed_tree$edge.length[[i]]*rlnorm(1, log(0.05) - ((0.2937025^2)/2), sdlog=0.2937025)
}

################################################################
################################################################
###SIMULATE_SEQUENCES###########################################
################################################################
################################################################

################################################################
###STANDARD_TREE_SEQUENCES######################################
################################################################

sequences_uncorrelated_relaxed <- simSeq(uncorrelated_relaxed_tree, l = 50000, type = "DNA", rate = 1)

################################################################
###DIV_SHIFT_SEQUENCES##########################################
################################################################

sequences_div_shift_uncorrelated_relaxed <- simSeq(div_shift_uncorrelated_relaxed_tree, l = 50000, type = "DNA", rate = 1)

###################################################################
###################################################################
###WRITE_OUTPUTS_TO_FILE###########################################
###################################################################
###################################################################

dir.create("sims_output")

###################################################################
###TREES###########################################################
###################################################################

write.tree(tree, "sims_output/constant_tree.tre") # constant time tree with outgroup
write.tree(drop.tip(tree, "t1001"), "sims_output/constant_tree_no_og.tre") # constant time tree without outgroup 
write.tree(uncorrelated_relaxed_tree, "sims_output/constant_uncorrelated_relaxed_tree.tre") # constant tree with ucln
write.tree(div_shift_tree, "sims_output/div_shift_tree.tre") # div shift time tree with outgroup
write.tree(drop.tip(div_shift_tree, "t1001"), "sims_output/div_shift_tree_no_og.tre") #div shift time tree without outgroup 
write.tree(div_shift_uncorrelated_relaxed_tree, "sims_output/div_shift_uncorrelated_relaxed_tree.tre") # div shift tree with ucln

###################################################################
###SEQUENCES#######################################################
###################################################################

###CONSTANT_SEQUENCES##############################################

write.phyDat(sequences_uncorrelated_relaxed, file = "sims_output/constant_uncorrelated_relaxed.nexus", format="nexus")
write.phyDat(subset(sequences_uncorrelated_relaxed, drop.tip(tree, "t1001")$tip.label), file = "sims_output/constant_uncorrelated_relaxed_no_og.nexus", format="nexus")

###DIV_SHIFT_SEQUENCES##############################################

write.phyDat(sequences_div_shift_uncorrelated_relaxed, file = "sims_output/div_shift_uncorrelated_relaxed.nexus", format="nexus")
write.phyDat(subset(sequences_div_shift_uncorrelated_relaxed, drop.tip(tree, "t1001")$tip.label), file = "sims_output/div_shift_uncorrelated_relaxed_no_og.nexus", format="nexus")
