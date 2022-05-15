#' performs core analysis and writes max age tree and min age tree to file
#' writes branch specific rates in max age tree and min age tree to file
#' @export

exTREEmaTIMEcore <- function(tree, rmax, rmin, noise_level, n_max_constraints, max_constraints_clade, max_constraints_ages, n_min_constraints, min_constraints_clade, min_constraints_ages, tip, sample_time){

###GET_KEY_TREE_VARIABLES###
############################

max_rate <- rmax
min_rate <- rmin

tip_vector <- seq(1, length(tree$tip.label), 1)
node_numbers <- seq(length(tip_vector)+1, (length(tip_vector)+2) + (length(tip_vector)-3), 1) 

clades <- vector("list", (length(tip_vector)-1))
for (i in 1:length(clades)){
clades[[i]] <- extract.clade(tree, node_numbers[[i]])$tip.label
}
clades <- clades[order(sapply(clades,length),decreasing=F)]

print("clades defined")

###GET_CORRECT_TIP_LABELS_FOR_POOL_AND_CONSTRAINTS###
#####################################################

tree_data_frame <- data.frame(tree[[1]][,1], tree[[1]][,2], tree$edge.length)
tree_data_frame <- tree_data_frame[order(tree_data_frame[,2]),]

pool <- vector("list", length(tree$tip.label))
for (i in 1:length(pool)){
pool[[i]] <- list(edge=matrix(c(2,1),1,2), tip.label=tree$tip.label[[i]], edge.length=tree_data_frame[,3][[i]], Nnode=1)
class(pool[[i]]) <- "phylo"
}

if (n_max_constraints > 0){
max_constraints_tip_label <- vector("list", n_max_constraints)
for (i in 1:length(max_constraints_tip_label)){
max_constraints_tip_label[[i]] <- extract.clade(tree, findMRCA(tree, max_constraints_clade[[i]], "node"))$tip.label
}
}

if (n_min_constraints > 0){
min_constraints_tip_label <- vector("list", n_min_constraints)
for (i in 1:length(min_constraints_tip_label)){
min_constraints_tip_label[[i]] <- extract.clade(tree, findMRCA(tree, min_constraints_clade[[i]], "node"))$tip.label
}
}

print("calibrations (if present) and other variables defined")

##############################################################################################
###CORE_CLADE_AGE_CALCULATION_FOR_MAX_AGE#####################################################
##############################################################################################

print("CALCULATING MAX AGE TREE")

print("initial generation of max age tree started")

for (i in 1:length(clades)){

###DETERMINE_WHAT_TO_JOIN###
############################

pool_to_include <- vector(mode="numeric", length=0)
for (a in 1:length(pool)){
counter <- vector(mode="numeric", length=0)
for (b in 1:length(pool[[a]]$tip.label)){
if (pool[[a]]$tip.label[[b]] %in% clades[[i]]){
counter <- append(counter, 1)
}
}
if (length(counter) == length(pool[[a]]$tip.label)){
pool_to_include <- append(pool_to_include, a)
}	
}

###JOIING_TWO_TIPS###
#####################

if ((length(pool[[pool_to_include[[1]]]]$tip.label) == 1) & (length(pool[[pool_to_include[[2]]]]$tip.label) == 1)){

if (pool[[pool_to_include[[1]]]]$edge.length[[1]] > noise_level & pool[[pool_to_include[[2]]]]$edge.length[[1]] > noise_level){

pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length/min_rate
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length/min_rate
if (pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]][[1]] < pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]][[1]]){
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]][[1]] - sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]][[1]]
} else {
pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]][[1]] - sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]][[1]]	
}

} else {

pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length/min_rate
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length/min_rate
if (pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]][[1]] > pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]][[1]]){
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]][[1]] - sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]][[1]]
} else {
pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]][[1]] - sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]][[1]]	
}

}

new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]][[1]]

###CONGRUIFY_TO_NODE_CALIBRATION######

if (n_max_constraints != 0){
for (a in 1:length(max_constraints_clade)){ 
if (setequal(new$tip.label, max_constraints_tip_label[[a]])){
if (true_age > max_constraints_ages[[a]]){
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] - (true_age - max_constraints_ages[[a]])
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] - (true_age - max_constraints_ages[[a]])
}
}
}
}

}

###JOINING_TWO_CLADES###
########################

else if ((length(pool[[pool_to_include[[1]]]]$tip.label) > 1) & (length(pool[[pool_to_include[[2]]]]$tip.label) > 1)){

if ((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]]) > noise_level & (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]]) > noise_level){

if ((((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]])/min_rate) + node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]])
 < 
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]])/min_rate) + node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]])){
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]]/min_rate))
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]]) - (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]])) 
} else {
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]]/min_rate))
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]]) - (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]])) 
}

} else {

if ((((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]])/min_rate) + node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]])
 > 
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]])/min_rate) + node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]])){
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]]/min_rate))
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]]) - (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]])) 
} else {
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]]/min_rate))
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]]) - (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]])) 
}

}

new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]][[1]]

###CONGRUIFY_TO_NODE_CALIBRATION######

if (n_max_constraints != 0){
for (a in 1:length(max_constraints_clade)){ 
if (setequal(new$tip.label, max_constraints_tip_label[[a]])){
if (true_age > max_constraints_ages[[a]]){
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] - (true_age - max_constraints_ages[[a]])
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] - (true_age - max_constraints_ages[[a]])
}
}
}
}

} 

###JOINING_ONE_TIP_WITH_A_CLADE###
##################################

else if ((length(pool[[pool_to_include[[1]]]]$tip.label) == 1) & (length(pool[[pool_to_include[[2]]]]$tip.label) > 1)){

if (pool[[pool_to_include[[1]]]]$edge.length > noise_level & tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]] > noise_level){

if (((pool[[pool_to_include[[1]]]]$edge.length/min_rate) + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]]) 
<
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]])/min_rate) + node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]])){
pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length/min_rate
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]]) - (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]]))
} else {
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]]/min_rate))
pool[[pool_to_include[[1]]]]$edge.length <- (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]]) - sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]]
}

} else {

if (((pool[[pool_to_include[[1]]]]$edge.length/min_rate) + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]]) 
>
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]])/min_rate) + node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]])){
pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length/min_rate
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]]) - (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]]))
} else {
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]]/min_rate))
pool[[pool_to_include[[1]]]]$edge.length <- (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]]) - sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]]
}

}

new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]][[1]]

###CONGRUIFY_TO_NODE_CALIBRATION######

if (n_max_constraints != 0){
for (a in 1:length(max_constraints_clade)){ 
if (setequal(new$tip.label, max_constraints_tip_label[[a]])){
if (true_age > max_constraints_ages[[a]]){
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] - (true_age - max_constraints_ages[[a]])
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] - (true_age - max_constraints_ages[[a]])
}
}
}
}

}

###ALTERNATE############################ 

else if ((length(pool[[pool_to_include[[2]]]]$tip.label) == 1) & (length(pool[[pool_to_include[[1]]]]$tip.label) > 1)){

if (pool[[pool_to_include[[2]]]]$edge.length > noise_level & tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]] > noise_level){

if (((pool[[pool_to_include[[2]]]]$edge.length/min_rate) + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]]) 
<
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]])/min_rate) + node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]])){
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length/min_rate
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]]) - (node.depth.edgelength(pool[[pool_to_include[[1]]]]$tip.label[[1]]) + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]]))
} else {
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]]/min_rate))
pool[[pool_to_include[[2]]]]$edge.length <- (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]]) - sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]]
} 

} else {

if (((pool[[pool_to_include[[2]]]]$edge.length/min_rate) + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]]) 
>
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]])/min_rate) + node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]])){
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length/min_rate
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]]) - (node.depth.edgelength(pool[[pool_to_include[[1]]]]$tip.label[[1]]) + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]]))
} else {
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]]/min_rate))
pool[[pool_to_include[[2]]]]$edge.length <- (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[1]]) - sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[1]]
}

}

new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]][[1]]

###CONGRUIFY_TO_NODE_CALIBRATION######

if (n_max_constraints != 0){
for (a in 1:length(max_constraints_clade)){ 
if (setequal(new$tip.label, max_constraints_tip_label[[a]])){
if (true_age > max_constraints_ages[[a]]){
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] - (true_age - max_constraints_ages[[a]])
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] - (true_age - max_constraints_ages[[a]])
}
}
}
}

}

pool <- pool[-c(pool_to_include[[1]], pool_to_include[[2]])]
pool[[length(pool)+1]] <- new

}

print("initial generation of max age tree complete")

###FILTER_OUT_POOR_BRANCHES##########

print("correcting branches in max age tree")

corresponding_tree <- extract.clade(tree, findMRCA(tree, pool[[1]]$tip.label, "node"))
min_lengths <- vector(mode="numeric", length=0)
for (a in 1:nrow(pool[[1]][[1]])){
if (pool[[1]][[1]][,2][[a]] > length(pool[[1]]$tip.label)){
if (corresponding_tree$edge.length[[which(corresponding_tree[[1]][,2] == findMRCA(corresponding_tree, extract.clade(pool[[1]], pool[[1]][[1]][,2][[a]])$tip.label, "node"))]] > noise_level){
min_lengths <- append(min_lengths, corresponding_tree$edge.length[[which(corresponding_tree[[1]][,2] == findMRCA(corresponding_tree, extract.clade(pool[[1]], pool[[1]][[1]][,2][[a]])$tip.label, "node"))]]/max_rate)
} else {
min_lengths <- append(min_lengths, 0)
}
}
if (pool[[1]][[1]][,2][[a]] <= length(pool[[1]]$tip.label)){
if (corresponding_tree$edge.length[[which(corresponding_tree[[1]][,2] == which(corresponding_tree$tip.label == pool[[1]]$tip.label[[pool[[1]][[1]][,2][[a]]]]))]] > noise_level){
min_lengths <- append(min_lengths, corresponding_tree$edge.length[[which(corresponding_tree[[1]][,2] == which(corresponding_tree$tip.label == pool[[1]]$tip.label[[pool[[1]][[1]][,2][[a]]]]))]]/max_rate)
} else {
min_lengths <- append(min_lengths, 0)
}
}

}

###CORRECT_THE_BRANCH_LENGTHS

continue_correction <- 0
while (continue_correction == 0){
negative_branch <- vector(mode="numeric", length=0)
minimum_permissable <- vector(mode="numeric", length=0)
for (a in 1:length(min_lengths)){
if (pool[[1]]$edge.length[[a]] < min_lengths[[a]]){
negative_branch <- append(negative_branch, a)
minimum_permissable <- append(minimum_permissable, min_lengths[[a]])
break
}
}
print(paste((a/nrow(pool[[1]][[1]]))*100, "% of branches confirmed as correct", sep=""))
if (length(negative_branch) > 0){
negative_branch <- negative_branch[[1]]
minimum_permissable <- minimum_permissable[[1]]
if (pool[[1]][[1]][,2][[negative_branch]] <= length(pool[[1]]$tip.label)){
keep <<-pool[[1]]
print(paste(paste("ISSUE: INCOMPATIBLE ASSUMPTIONS, MAX AGE ", a, sep=""), " TOO YOUNG OR MAX RATE TOO LOW", sep=""))
max_rate_unit_number <<- append(max_rate_unit_number, 1)
print(pool[[1]])
return()
}
descendant <- extract.clade(pool[[1]], pool[[1]][[1]][,2][[negative_branch]])
pool[[1]]$edge.length[[which(pool[[1]][[1]][,1] == findMRCA(pool[[1]], descendant$tip.label, "node"))[[1]]]] <- pool[[1]]$edge.length[[which(pool[[1]][[1]][,1] == findMRCA(pool[[1]], descendant$tip.label, "node"))[[1]]]] - (minimum_permissable - pool[[1]]$edge.length[[negative_branch]])
pool[[1]]$edge.length[[which(pool[[1]][[1]][,1] == findMRCA(pool[[1]], descendant$tip.label, "node"))[[2]]]] <- pool[[1]]$edge.length[[which(pool[[1]][[1]][,1] == findMRCA(pool[[1]], descendant$tip.label, "node"))[[2]]]] - (minimum_permissable - pool[[1]]$edge.length[[negative_branch]])
pool[[1]]$edge.length[[negative_branch]] <- pool[[1]]$edge.length[[negative_branch]] + (minimum_permissable - pool[[1]]$edge.length[[negative_branch]])
} else {
continue_correction <- 1
}
}

max_ages_tree <<- pool[[1]]

print("correction of branches in max age tree completed")

print("MAX AGE TREE CALCULATED")

############################################################################################
###INITIAL_PROCESSING_TWO###################################################################
############################################################################################

###RE-POOL###
#############

pool <- vector("list", length(tree$tip.label))
for (i in 1:length(pool)){
pool[[i]] <- list(edge=matrix(c(2,1),1,2), tip.label=tree$tip.label[[i]], edge.length=tree_data_frame[,3][[i]], Nnode=1)
class(pool[[i]]) <- "phylo"
}

##############################################################################################
###CORE_CLADE_AGE_CALCULATION_FOR_MIN_AGE#####################################################
##############################################################################################

print("CALCULATING MIN AGE TREE")

print("initial generation of min age tree started")

for (i in 1:length(clades)){

###DETERMINE_WHAT_TO_JOIN###
############################

pool_to_include <- vector(mode="numeric", length=0)
for (a in 1:length(pool)){
counter <- vector(mode="numeric", length=0)
for (b in 1:length(pool[[a]]$tip.label)){
if (pool[[a]]$tip.label[[b]] %in% clades[[i]]){
counter <- append(counter, 1)
}
}
if (length(counter) == length(pool[[a]]$tip.label)){
pool_to_include <- append(pool_to_include, a)
}	
}

###JOIING_TWO_TIPS###
#####################

if ((length(pool[[pool_to_include[[1]]]]$tip.label) == 1) & (length(pool[[pool_to_include[[2]]]]$tip.label) == 1)){
pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length/max_rate
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length/max_rate
if (pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]][[2]] > pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]][[2]]){
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]][[2]] - sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]][[2]]
} else {
pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]][[2]] - sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]][[2]]	
}

new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]][[2]]

###CONGRUIFY_TO_NODE_CALIBRATION######

if (n_min_constraints != 0){
for (a in 1:length(min_constraints_clade)){ 
if (setequal(new$tip.label, min_constraints_tip_label[[a]])){
if (true_age < min_constraints_ages[[a]]){
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] + (min_constraints_ages[[a]]-true_age)
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] + (min_constraints_ages[[a]]-true_age)
}
}
}
}

}

###JOINING_TWO_CLADES###
########################

else if ((length(pool[[pool_to_include[[1]]]]$tip.label) > 1) & (length(pool[[pool_to_include[[2]]]]$tip.label) > 1)){ 
if ((((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]])/max_rate) + node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[2]])
 > 
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]])/max_rate) + node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[2]])){
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]]/max_rate))
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[2]]) - (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[2]])) 
} else {
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]]/max_rate))
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[2]]) - (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[2]])) 
}

new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]][[2]]

###CONGRUIFY_TO_NODE_CALIBRATION######

if (n_min_constraints != 0){
for (a in 1:length(min_constraints_clade)){ 
if (setequal(new$tip.label, min_constraints_tip_label[[a]])){
if (true_age < min_constraints_ages[[a]]){
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] + (min_constraints_ages[[a]]-true_age)
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] + (min_constraints_ages[[a]]-true_age)
}
}
}
}

} 

###JOINING_ONE_TIP_WITH_A_CLADE###
##################################

else if ((length(pool[[pool_to_include[[1]]]]$tip.label) == 1) & (length(pool[[pool_to_include[[2]]]]$tip.label) > 1)){
if (((pool[[pool_to_include[[1]]]]$edge.length/max_rate) + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[2]]) 
>
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]])/max_rate) + node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[2]])){
pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length/max_rate
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[2]]) - (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[2]]))
} else {
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]]/max_rate))
pool[[pool_to_include[[1]]]]$edge.length <- (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[2]]) - sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[2]]
} 

new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]][[2]]

###CONGRUIFY_TO_NODE_CALIBRATION######

if (n_min_constraints != 0){
for (a in 1:length(min_constraints_clade)){ 
if (setequal(new$tip.label, min_constraints_tip_label[[a]])){
if (true_age < min_constraints_ages[[a]]){
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] + (min_constraints_ages[[a]]-true_age)
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] + (min_constraints_ages[[a]]-true_age)
}
}
}
}

}

###ALTERNATE 

else if ((length(pool[[pool_to_include[[2]]]]$tip.label) == 1) & (length(pool[[pool_to_include[[1]]]]$tip.label) > 1)){
if (((pool[[pool_to_include[[2]]]]$edge.length/max_rate) + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[2]]) 
>
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]])/max_rate) + node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[2]])){
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length/max_rate
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[2]]) - (node.depth.edgelength(pool[[pool_to_include[[1]]]]$tip.label[[1]]) + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[2]]))
} else {
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]]/max_rate))
pool[[pool_to_include[[2]]]]$edge.length <- (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]][[2]]) - sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]][[2]]
} 
new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]][[2]]

###CONGRUIFY_TO_NODE_CALIBRATION######

if (n_min_constraints != 0){
for (a in 1:length(min_constraints_clade)){ 
if (setequal(new$tip.label, min_constraints_tip_label[[a]])){
if (true_age < min_constraints_ages[[a]]){
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] + (min_constraints_ages[[a]]-true_age)
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] + (min_constraints_ages[[a]]-true_age)
}
}
}
}

}

pool <- pool[-c(pool_to_include[[1]], pool_to_include[[2]])]
pool[[length(pool)+1]] <- new

}

print("initial generation of min age tree complete")

###FILTER_OUT_POOR_BRANCHES##########

print("correcting branches in min age tree")

corresponding_tree <- extract.clade(tree, findMRCA(tree, pool[[1]]$tip.label, "node"))
max_lengths <- vector(mode="numeric", length=0)
for (a in 1:nrow(pool[[1]][[1]])){
if (pool[[1]][[1]][,2][[a]] > length(pool[[1]]$tip.label)){
if (corresponding_tree$edge.length[[which(corresponding_tree[[1]][,2] == findMRCA(corresponding_tree, extract.clade(pool[[1]], pool[[1]][[1]][,2][[a]])$tip.label, "node"))]] > noise_level){
max_lengths <- append(max_lengths, corresponding_tree$edge.length[[which(corresponding_tree[[1]][,2] == findMRCA(corresponding_tree, extract.clade(pool[[1]], pool[[1]][[1]][,2][[a]])$tip.label, "node"))]]/min_rate)
} else {
max_lengths <- append(max_lengths, 1000000000)
}
}
if (pool[[1]][[1]][,2][[a]] <= length(pool[[1]]$tip.label)){
if (corresponding_tree$edge.length[[which(corresponding_tree[[1]][,2] == which(corresponding_tree$tip.label == pool[[1]]$tip.label[[pool[[1]][[1]][,2][[a]]]]))]] > noise_level){
max_lengths <- append(max_lengths, corresponding_tree$edge.length[[which(corresponding_tree[[1]][,2] == which(corresponding_tree$tip.label == pool[[1]]$tip.label[[pool[[1]][[1]][,2][[a]]]]))]]/min_rate)
} else {
max_lengths <- append(max_lengths, 1000000000)
}
}

}

###CORRECT_THE_BRANCH_LENGTHS

continue_correction <- 0
while (continue_correction == 0){
long_branch <- vector(mode="numeric", length=0)
maximum_permissable <- vector(mode="numeric", length=0)
for (a in 1:length(max_lengths)){
if (pool[[1]]$edge.length[[a]] > max_lengths[[a]]){
long_branch <- append(long_branch, a)
maximum_permissable <- append(maximum_permissable, max_lengths[[a]])
break
}
}
print(paste((a/nrow(pool[[1]][[1]]))*100, "% of branches confirmed as correct", sep=""))
if (length(long_branch) > 0){
long_branch <- long_branch[[1]]
maximum_permissable <- maximum_permissable[[1]]
if (pool[[1]][[1]][,2][[long_branch]] <= length(pool[[1]]$tip.label)){
keep <<-pool[[1]]
print(paste(paste("ISSUE: INCOMPATIBLE ASSUMPTIONS, MIN AGE ", a, sep=""), " TOO OLD OR MIN RATE TOO HIGH", sep=""))
min_rate_unit_number <<- append(min_rate_unit_number, 1)
print(pool[[1]])
return()
}
descendant <- extract.clade(pool[[1]], pool[[1]][[1]][,2][[long_branch]])
pool[[1]]$edge.length[[which(pool[[1]][[1]][,1] == findMRCA(pool[[1]], descendant$tip.label, "node"))[[1]]]] <- pool[[1]]$edge.length[[which(pool[[1]][[1]][,1] == findMRCA(pool[[1]], descendant$tip.label, "node"))[[1]]]] + (pool[[1]]$edge.length[[long_branch]]-maximum_permissable)
pool[[1]]$edge.length[[which(pool[[1]][[1]][,1] == findMRCA(pool[[1]], descendant$tip.label, "node"))[[2]]]] <- pool[[1]]$edge.length[[which(pool[[1]][[1]][,1] == findMRCA(pool[[1]], descendant$tip.label, "node"))[[2]]]] + (pool[[1]]$edge.length[[long_branch]]-maximum_permissable)
pool[[1]]$edge.length[[long_branch]] <- pool[[1]]$edge.length[[long_branch]] - (pool[[1]]$edge.length[[long_branch]]-maximum_permissable)
} else {
continue_correction <- 1
}
}

min_ages_tree <<- pool[[1]]

print("correction of branches in max age tree completed")

print("MAX AGE TREE CALCULATED")

#################################
###PLOT_BRANCH_RATES#############
#################################

print("EXTRACTING_RATES")

max_ages_tree_data_frame <- data.frame(max_ages_tree[[1]], max_ages_tree$edge.length)
min_ages_tree_data_frame <- data.frame(min_ages_tree[[1]], min_ages_tree$edge.length)

max_ages_tree_branch_rates <- vector(mode="numeric", length=0)
for (j in 1:nrow(max_ages_tree_data_frame)){
if (max_ages_tree_data_frame[,2][[j]] > length(max_ages_tree$tip.label)){
max_ages_tree_branch_rates <- append(max_ages_tree_branch_rates, tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, extract.clade(max_ages_tree, max_ages_tree_data_frame[,2][[j]])$tip.label, "node"))]]/max_ages_tree_data_frame[,3][[j]])
}
if (max_ages_tree_data_frame[,2][[j]] <= length(max_ages_tree$tip.label)){
max_ages_tree_branch_rates <- append(max_ages_tree_branch_rates, tree_data_frame[,3][[which(tree_data_frame[,2] == which(tree$tip.label == max_ages_tree$tip.label[[max_ages_tree_data_frame[,2][[j]]]]))]]/max_ages_tree_data_frame[,3][[j]])
}
}

min_ages_tree_branch_rates <- vector(mode="numeric", length=0)
for (j in 1:nrow(min_ages_tree_data_frame)){
if (min_ages_tree_data_frame[,2][[j]] > length(min_ages_tree$tip.label)){
min_ages_tree_branch_rates <- append(min_ages_tree_branch_rates, tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, extract.clade(min_ages_tree, min_ages_tree_data_frame[,2][[j]])$tip.label, "node"))]]/min_ages_tree_data_frame[,3][[j]])
}
if (min_ages_tree_data_frame[,2][[j]] <= length(min_ages_tree$tip.label)){
min_ages_tree_branch_rates <- append(min_ages_tree_branch_rates, tree_data_frame[,3][[which(tree_data_frame[,2] == which(tree$tip.label == min_ages_tree$tip.label[[min_ages_tree_data_frame[,2][[j]]]]))]]/min_ages_tree_data_frame[,3][[j]])
}
}

continue_through_units <<-0
rates <<- data.frame(max_ages_tree_branch_rates, min_ages_tree_branch_rates)
min_ages_tree <<- min_ages_tree
max_ages_tree <<- max_ages_tree
print("FINISHED")
}
