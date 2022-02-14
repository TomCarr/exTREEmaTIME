library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(Rmisc)
library(TreePar)
library(NELSI)

root_time <- 1.897911 # set to correct value
noise_level <- 0.005

##############################
###BRANCH_INDEX_INFORMATION###
##############################

index_tree <- read.annotated.nexus("estimate_molecular_branch_lengths/output/div_shift_uncorrelated_relaxed_index_tree.tre") # change to relevant tree

branch_indices <- vector(mode="numeric", length=0)
for (i in 1:length(index_tree[[6]])){
branch_indices <- append(branch_indices, unlist(index_tree[[6]][[i]]))
}

n_branches <- length(index_tree$tip.label)*2 - 3

###############################
###POSTEIOR_TREE_INFORMATION###
###############################

burnin = 0.25

###Log### #nb tree files with constrained topology doesn't print properly, so have to extract bl estimate from log file

run_one <- read.table("estimate_molecular_branch_lengths/output/div_shift_uncorrelated_relaxed_bl.log", header=F) # change to relevant directory

###Trees###

trees_table <- read.table("estimate_molecular_branch_lengths/output/div_shift_uncorrelated_relaxed_bl.trees", header=F) # change to relevant directory
trees <- vector ("list", nrow(trees_table)-round(burnin*nrow(trees_table), digits = 0))
for (k in seq(round(burnin*nrow(trees_table), digits = 0)+1, nrow(trees_table), 1)){
write(as.character(trees_table[,5][[k]]), "tree.tre")
trees[[k - (round(burnin*nrow(trees_table), digits = 0))]] <- read.tree("tree.tre")
}

######################
###DO_FULL_ANALYSIS###
######################

final_max_diff <- vector(mode="numeric", length=0)
final_min_diff <- vector(mode="numeric", length=0)
final_mean <- vector(mode="numeric", length=0)
final_sd <- vector(mode="numeric", length=0)
trees_edges_storage <- vector("list", length(trees[[1]]$edge.length))

for (z in 1:length(trees)){

tree <- trees[[z]]
trees_edges <- vector("list", n_branches)

###get_branch_lengths_from_posterior###

branch_values <- vector("list", n_branches)
for (i in 5:(n_branches+4)){
branch_values[[i-4]] <- as.numeric(run_one[,i][[round(burnin*nrow(trees_table), digits = 0) + z]])
}

branch_data_frame <- data.frame(index_tree[[1]][,1], index_tree[[1]][,2], seq(1, n_branches, 1), branch_indices[-length(branch_indices)])
branch_data_frame <- branch_data_frame[order(branch_data_frame[,4]),] 
branch_data_frame <- cbind(branch_data_frame, unlist(branch_values))
branch_data_frame <- branch_data_frame[order(branch_data_frame[,3]),]

for (i in 1:nrow(branch_data_frame)){
tree$edge.length[[i]] <- branch_data_frame[,5][[i]]
trees_edges_storage[[i]] <- append(trees_edges_storage[[i]], branch_data_frame[,5][[i]])
}

###process_tree###

tree <- root(tree, "t1001")
tree <- drop.tip(tree, "t1001")

clades_tips <- vector("list", length(tree$tip.label) - 1)
for (i in 1:length(clades_tips)){
clades_tips[[i]] <- extract.clade(tree, length(tree$tip.label) + i)$tip.label
}

retention <- vector(mode="numeric", length=0)
for (i in 1:length(clades_tips)){
if (length(clades_tips[[i]]) == 2){
retention <- append(retention, i)
}
}
clades_tips <- clades_tips[retention]

clades <- vector("list", length(clades_tips))
for (i in 1:length(clades)){
clades[[i]] <- extract.clade(tree, findMRCA(tree, clades_tips[[i]], "node"))
}

###get_rate_mean_from_posterior_replicate###

mean <- mean(node.depth.edgelength(tree)[seq(1, length(tree$tip.label), 1)])/root_time

###get_rate_sd_from_posterior_replicate###

sd_vals <- vector(mode="numeric", length=0)
for (i in 1:length(retention)){
if ((tree$edge.length[[which(tree[[1]][,2] == which(tree$tip.label == clades_tips[[i]][[1]]))]] > noise_level)
&
(tree$edge.length[[which(tree[[1]][,2] == which(tree$tip.label == clades_tips[[i]][[2]]))]] > noise_level)){
sd_vals <- append(sd_vals, (tree$edge.length[[which(tree[[1]][,2] == which(tree$tip.label == clades_tips[[i]][[1]]))]]/mean(clades[[i]]$edge.length)))
sd_vals <- append(sd_vals, (tree$edge.length[[which(tree[[1]][,2] == which(tree$tip.label == clades_tips[[i]][[2]]))]]/mean(clades[[i]]$edge.length)))
}
}

###get_min_max_from_posterior_replicate### #for rate_max_means and rate_min_means in rev shell script. This defines boundaries of clock rate multiplier in random local clock

remove <- vector(mode="numeric", length=0)
min_max <- vector("list", length(clades))
for (i in 1:length(retention)){
if ((tree$edge.length[[which(tree[[1]][,2] == which(tree$tip.label == clades_tips[[i]][[1]]))]] > noise_level)
&
(tree$edge.length[[which(tree[[1]][,2] == which(tree$tip.label == clades_tips[[i]][[2]]))]] > noise_level)){
min_max[[i]] <- append(min_max[[i]], (tree$edge.length[[which(tree[[1]][,2] == which(tree$tip.label == clades_tips[[i]][[1]]))]]))
min_max[[i]] <- append(min_max[[i]], (tree$edge.length[[which(tree[[1]][,2] == which(tree$tip.label == clades_tips[[i]][[2]]))]]))
min_max[[i]] <- max(min_max[[i]])/min(min_max[[i]])
}
}
if (sum(unlist(min_max))>0){
max_diff <- min_max[[which(unlist(min_max) == max(unlist(min_max)))]]
}

###final_storage###

final_mean <- append(final_mean, mean)

if (length(var_vals)>0){
final_sd <- append(final_sd, sd(log(sd_vals)))
}

if (sum(unlist(min_max)>0)){
final_max_diff <- append(final_max_diff, max_diff)
final_min_diff <- append(final_min_diff, 1/max_diff)
}

}

###########################################################
###WRITE_BRANCH_LENGTH_TREE_TO_FILE_FOR_FURTHER_ANALYSIS###
###########################################################

branch_length_tree_edge_lengths <- vector(mode="numeric", length=0)
for (i in 1:length(trees_edges_storage)){
branch_length_tree_edge_lengths <- append(branch_length_tree_edge_lengths, mean(trees_edges_storage[[i]]))
}

branch_length_tree <- trees[[1]]
branch_length_tree$edge.length <- branch_length_tree_edge_lengths

write.tree(branch_length_tree, "estimate_molecular_branch_lengths/output/div_shift_uncorrelated_relaxed_bl_tree.tre") # change to relevant directory