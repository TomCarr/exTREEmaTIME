SetAutoRates <- function(tree, root_min, root_max, unit, auto_rates_type, noise_level){

############################
###VERSION_ONE##############
############################

if(auto_rates_type == 1){

rmax_calc <<- (max(node.depth.edgelength(tree)[seq(1, length(tree$tip.label), 1)])*unit)/(root_min + ((root_max-root_min)/2))
print(rmax)
rmin_calc <<- (min(node.depth.edgelength(tree)[seq(1, length(tree$tip.label), 1)])/unit)/(root_min + ((root_max-root_min)/2))
print(rmin)

}

############################
###VERSION_TWO##############
############################

if(auto_rates_type == 2){

###define_clades_and_tips###

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

###define_diff_values######

var_vals <- vector("list", length(clades))
for (i in 1:length(retention)){
if ((tree$edge.length[[which(tree[[1]][,2] == which(tree$tip.label == clades_tips[[i]][[1]]))]] > noise_level)
&
(tree$edge.length[[which(tree[[1]][,2] == which(tree$tip.label == clades_tips[[i]][[2]]))]] > noise_level)){
var_vals[[i]] <- append(var_vals[[i]], (tree$edge.length[[which(tree[[1]][,2] == which(tree$tip.label == clades_tips[[i]][[1]]))]]))
var_vals[[i]] <- append(var_vals[[i]], (tree$edge.length[[which(tree[[1]][,2] == which(tree$tip.label == clades_tips[[i]][[2]]))]]))
var_vals[[i]] <- max(var_vals[[i]])/min(var_vals[[i]])
}
}

if (sum(unlist(var_vals))>0){
max_diff <- unlist(var_vals)[[which(unlist(var_vals) == max(unlist(var_vals)))]]
}

###set_r_min_and_r_max#####

if (sum(unlist(var_vals))>0){
rmax_calc <<- (mean(node.depth.edgelength(tree)[seq(1, length(tree$tip.label), 1)])/root_min)*(max_diff*unit)
print(rmax_calc)
rmin_calc <<- (mean(node.depth.edgelength(tree)[seq(1, length(tree$tip.label), 1)])/root_max)/(max_diff*unit)
print(rmin_calc)
rate_set <<- 1
} else {
print("SPECIFIED NOISE LEVEL EXCEEDED")
return()
}

}

}
