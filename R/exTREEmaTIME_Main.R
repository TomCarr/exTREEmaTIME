#' performs core analysis and writes max age tree and min age tree to file
#' writes branch specific rates in max age tree and min age tree to file
#' @export
exTREEmaTIMEcalc <- function(tree, rmax, rmin, n_max_constraints, max_constraint_clade, max_constraint_age, n_min_constraints, min_constraint_clade, min_constraint_age, calibration_implementation_precision, sample_time){

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

##############################################################################################
###CORE_CLADE_AGE_CALCULATION_FOR_MAX_AGE#####################################################
##############################################################################################

print("CALCULATING MAX AGE TREE")

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

print(paste(paste("maximum age estimate", i, ""), " initiated", "")) 

###JOIING_TWO_TIPS###
#####################

if ((length(pool[[pool_to_include[[1]]]]$tip.label) == 1) & (length(pool[[pool_to_include[[2]]]]$tip.label) == 1)){
print("joining two tips")
pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length/min_rate
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length/min_rate
if (pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]] < pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]]){
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]] - sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]]
} else {
pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]] - sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]]	
}

new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]]

###CONGRUIFY_TO_NODE_CALIBRATION######

if (n_max_constraints != 0){
for (a in 1:length(max_constraints_clade)){
if (setequal(new$tip.label, max_constraints_tip_label[[a]])){
print("node calibration being implemented")
while (true_age > max_constraints_age[[a]]){
true_age <- true_age - calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] - calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] - calibration_implementation_precision
}
}
}
}

###FILTER_OUT_POOR_BRANCHES##########

corresponding_tree <- extract.clade(tree, findMRCA(tree, new$tip.label, "node"))
min_lengths <- vector(mode="numeric", length=0)
for (a in 1:nrow(new[[1]])){
for (b in 1:nrow(corresponding_tree[[1]])){
if (new[[1]][,2][[a]] > length(new$tip.label) & corresponding_tree[[1]][,2][[b]] > length(corresponding_tree$tip.label)){
if (setequal(extract.clade(new, new[[1]][,2][[a]])$tip.label, extract.clade(corresponding_tree, corresponding_tree[[1]][,2][[b]])$tip.label)){
min_lengths <- append(min_lengths, corresponding_tree$edge.length[[b]]/max_rate)
}
}
if (new[[1]][,2][[a]] <= length(new$tip.label) & corresponding_tree[[1]][,2][[b]] <= length(corresponding_tree$tip.label)){
if (new$tip.label[[new[[1]][,2][[a]]]] == corresponding_tree$tip.label[[corresponding_tree[[1]][,2][[b]]]]){
min_lengths <- append(min_lengths, corresponding_tree$edge.length[[b]]/max_rate)
}
}
}
}

continue_correction <- 0
while (continue_correction == 0){
negative_branch <- vector(mode="numeric", length=0)
for (a in 1:length(min_lengths)){
if (new$edge.length[[a]] < min_lengths[[a]]){
negative_branch <- append(negative_branch, a)
}
}
if (length(negative_branch) > 0){
negative_branch <- negative_branch[[1]]
if (new[[1]][,2][[negative_branch]] <= length(new$tip.label)){
print(paste(paste("ERROR: INCOMPATIBLE ASSUMPTIONS, MAX AGE ", a, sep=""), " TOO YOUNG OR MAX RATE TOO LOW", sep=""))
quit()
}
descendant <- extract.clade(new, new[[1]][,2][[negative_branch]])
ancestral <- extract.clade(new, new[[1]][,1][[negative_branch]])
min_diff <- tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, descendant$tip.label, "node"))]]/max_rate
keep_reducing <- 0
while (keep_reducing == 0){
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] - calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] - calibration_implementation_precision
new$edge.length[[negative_branch]] <- new$edge.length[[negative_branch]] + calibration_implementation_precision
if (max(node.depth.edgelength(extract.clade(new, findMRCA(new, descendant$tip.label, "node")))) <= max(node.depth.edgelength(extract.clade(new, findMRCA(new, ancestral$tip.label, "node")))) - min_diff){
keep_reducing <- 1
}
}
} else {
continue_correction <- 1
}
}

}

###JOINING_TWO_CLADES###
########################

else if ((length(pool[[pool_to_include[[1]]]]$tip.label) > 1) & (length(pool[[pool_to_include[[2]]]]$tip.label) > 1)){
print("joining two clades")
if ((((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]])/min_rate) + node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]])
 < 
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]])/min_rate) + node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]])){
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]]/min_rate))
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]]) - (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]])) 
} else {
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]]/min_rate))
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]]) - (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]])) 
}

new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]]

###CONGRUIFY_TO_NODE_CALIBRATION######

if (n_max_constraints != 0){
for (a in 1:length(max_constraints_clade)){
if (setequal(new$tip.label, max_constraints_tip_label[[a]])){
print("age shift due to node calibration")
while (true_age > max_constraints_age[[a]]){
true_age <- true_age - calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] - calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] - calibration_implementation_precision
}
}
}
}

###FILTER_OUT_POOR_BRANCHES##########

corresponding_tree <- extract.clade(tree, findMRCA(tree, new$tip.label, "node"))
min_lengths <- vector(mode="numeric", length=0)
for (a in 1:nrow(new[[1]])){
for (b in 1:nrow(corresponding_tree[[1]])){
if (new[[1]][,2][[a]] > length(new$tip.label) & corresponding_tree[[1]][,2][[b]] > length(corresponding_tree$tip.label)){
if (setequal(extract.clade(new, new[[1]][,2][[a]])$tip.label, extract.clade(corresponding_tree, corresponding_tree[[1]][,2][[b]])$tip.label)){
min_lengths <- append(min_lengths, corresponding_tree$edge.length[[b]]/max_rate)
}
}
if (new[[1]][,2][[a]] <= length(new$tip.label) & corresponding_tree[[1]][,2][[b]] <= length(corresponding_tree$tip.label)){
if (new$tip.label[[new[[1]][,2][[a]]]] == corresponding_tree$tip.label[[corresponding_tree[[1]][,2][[b]]]]){
min_lengths <- append(min_lengths, corresponding_tree$edge.length[[b]]/max_rate)
}
}
}
}

continue_correction <- 0
while (continue_correction == 0){
negative_branch <- vector(mode="numeric", length=0)
for (a in 1:length(min_lengths)){
if (new$edge.length[[a]] < min_lengths[[a]]){
negative_branch <- append(negative_branch, a)
}
}
if (length(negative_branch) > 0){
negative_branch <- negative_branch[[1]]
if (new[[1]][,2][[negative_branch]] <= length(new$tip.label)){
print(paste(paste("ERROR: INCOMPATIBLE ASSUMPTIONS, MAX AGE ", a, sep=""), " TOO YOUNG OR MAX RATE TOO LOW", sep=""))
quit()
}
descendant <- extract.clade(new, new[[1]][,2][[negative_branch]])
ancestral <- extract.clade(new, new[[1]][,1][[negative_branch]])
min_diff <- tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, descendant$tip.label, "node"))]]/max_rate
keep_reducing <- 0
while (keep_reducing == 0){
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] - calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] - calibration_implementation_precision
new$edge.length[[negative_branch]] <- new$edge.length[[negative_branch]] + calibration_implementation_precision
if (max(node.depth.edgelength(extract.clade(new, findMRCA(new, descendant$tip.label, "node")))) <= max(node.depth.edgelength(extract.clade(new, findMRCA(new, ancestral$tip.label, "node")))) - min_diff){
keep_reducing <- 1
}
}
} else {
continue_correction <- 1
}
}

} 

###JOINING_ONE_TIP_WITH_A_CLADE###
##################################

else if ((length(pool[[pool_to_include[[1]]]]$tip.label) == 1) & (length(pool[[pool_to_include[[2]]]]$tip.label) > 1)){
print("joining one tip with a clade")
if (((pool[[pool_to_include[[1]]]]$edge.length/min_rate) + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]]) 
<
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]])/min_rate) + node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]])){
pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length/min_rate
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]]) - (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]]))
} else {
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]]/min_rate))
pool[[pool_to_include[[1]]]]$edge.length <- (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]]) - sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]]
}
 
new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]]

###CONGRUIFY_TO_NODE_CALIBRATION######

if (n_max_constraints != 0){
for (a in 1:length(max_constraints_clade)){
if (setequal(new$tip.label, max_constraints_tip_label[[a]])){
while (true_age > max_constraints_age[[a]]){
true_age <- true_age - calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] - calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] - calibration_implementation_precision
}
}
}
}

###FILTER_OUT_POOR_BRANCHES##########

corresponding_tree <- extract.clade(tree, findMRCA(tree, new$tip.label, "node"))
min_lengths <- vector(mode="numeric", length=0)
for (a in 1:nrow(new[[1]])){
for (b in 1:nrow(corresponding_tree[[1]])){
if (new[[1]][,2][[a]] > length(new$tip.label) & corresponding_tree[[1]][,2][[b]] > length(corresponding_tree$tip.label)){
if (setequal(extract.clade(new, new[[1]][,2][[a]])$tip.label, extract.clade(corresponding_tree, corresponding_tree[[1]][,2][[b]])$tip.label)){
min_lengths <- append(min_lengths, corresponding_tree$edge.length[[b]]/max_rate)
}
}
if (new[[1]][,2][[a]] <= length(new$tip.label) & corresponding_tree[[1]][,2][[b]] <= length(corresponding_tree$tip.label)){
if (new$tip.label[[new[[1]][,2][[a]]]] == corresponding_tree$tip.label[[corresponding_tree[[1]][,2][[b]]]]){
min_lengths <- append(min_lengths, corresponding_tree$edge.length[[b]]/max_rate)
}
}
}
}

continue_correction <- 0
while (continue_correction == 0){
negative_branch <- vector(mode="numeric", length=0)
for (a in 1:length(min_lengths)){
if (new$edge.length[[a]] < min_lengths[[a]]){
negative_branch <- append(negative_branch, a)
}
}
if (length(negative_branch) > 0){
negative_branch <- negative_branch[[1]]
if (new[[1]][,2][[negative_branch]] <= length(new$tip.label)){
print(paste(paste("ERROR: INCOMPATIBLE ASSUMPTIONS, MAX AGE ", a, sep=""), " TOO YOUNG OR MAX RATE TOO LOW", sep=""))
quit()
}
descendant <- extract.clade(new, new[[1]][,2][[negative_branch]])
ancestral <- extract.clade(new, new[[1]][,1][[negative_branch]])
min_diff <- tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, descendant$tip.label, "node"))]]/max_rate
keep_reducing <- 0
while (keep_reducing == 0){
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] - calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] - calibration_implementation_precision
new$edge.length[[negative_branch]] <- new$edge.length[[negative_branch]] + calibration_implementation_precision
if (max(node.depth.edgelength(extract.clade(new, findMRCA(new, descendant$tip.label, "node")))) <= max(node.depth.edgelength(extract.clade(new, findMRCA(new, ancestral$tip.label, "node")))) - min_diff){
keep_reducing <- 1
}
}
} else {
continue_correction <- 1
}
}

}

###ALTERNATE############################ 

else if ((length(pool[[pool_to_include[[2]]]]$tip.label) == 1) & (length(pool[[pool_to_include[[1]]]]$tip.label) > 1)){
print("joining one tip with a clade")
if (((pool[[pool_to_include[[2]]]]$edge.length/min_rate) + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]]) 
<
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]])/min_rate) + node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]])){
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length/min_rate
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]]) - (node.depth.edgelength(pool[[pool_to_include[[1]]]]$tip.label[[1]]) + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]]))
} else {
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]]/min_rate))
pool[[pool_to_include[[2]]]]$edge.length <- (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]]) - sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]]
} 

new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]]

###CONGRUIFY_TO_NODE_CALIBRATION######

if (n_max_constraints != 0){
for (a in 1:length(max_constraints_clade)){
if (setequal(new$tip.label, max_constraints_tip_label[[a]])){
print("age shift due to node calibration")
while (true_age > max_constraints_age[[a]]){
true_age <- true_age - calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] - calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] - calibration_implementation_precision
}
}
}
}

###FILTER_OUT_POOR_BRANCHES##########

corresponding_tree <- extract.clade(tree, findMRCA(tree, new$tip.label, "node"))
min_lengths <- vector(mode="numeric", length=0)
for (a in 1:nrow(new[[1]])){
for (b in 1:nrow(corresponding_tree[[1]])){
if (new[[1]][,2][[a]] > length(new$tip.label) & corresponding_tree[[1]][,2][[b]] > length(corresponding_tree$tip.label)){
if (setequal(extract.clade(new, new[[1]][,2][[a]])$tip.label, extract.clade(corresponding_tree, corresponding_tree[[1]][,2][[b]])$tip.label)){
min_lengths <- append(min_lengths, corresponding_tree$edge.length[[b]]/max_rate)
}
}
if (new[[1]][,2][[a]] <= length(new$tip.label) & corresponding_tree[[1]][,2][[b]] <= length(corresponding_tree$tip.label)){
if (new$tip.label[[new[[1]][,2][[a]]]] == corresponding_tree$tip.label[[corresponding_tree[[1]][,2][[b]]]]){
min_lengths <- append(min_lengths, corresponding_tree$edge.length[[b]]/max_rate)
}
}
}
}

continue_correction <- 0
while (continue_correction == 0){
negative_branch <- vector(mode="numeric", length=0)
for (a in 1:length(min_lengths)){
if (new$edge.length[[a]] < min_lengths[[a]]){
negative_branch <- append(negative_branch, a)
}
}
if (length(negative_branch) > 0){
negative_branch <- negative_branch[[1]]
if (new[[1]][,2][[negative_branch]] <= length(new$tip.label)){
print(paste(paste("ERROR: INCOMPATIBLE ASSUMPTIONS, MAX AGE ", a, sep=""), " TOO YOUNG OR MAX RATE TOO LOW", sep=""))
quit()
}
descendant <- extract.clade(new, new[[1]][,2][[negative_branch]])
ancestral <- extract.clade(new, new[[1]][,1][[negative_branch]])
min_diff <- tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, descendant$tip.label, "node"))]]/max_rate
keep_reducing <- 0
while (keep_reducing == 0){
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] - calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] - calibration_implementation_precision
new$edge.length[[negative_branch]] <- new$edge.length[[negative_branch]] + calibration_implementation_precision
if (max(node.depth.edgelength(extract.clade(new, findMRCA(new, descendant$tip.label, "node")))) <= max(node.depth.edgelength(extract.clade(new, findMRCA(new, ancestral$tip.label, "node")))) - min_diff){
keep_reducing <- 1
}
}
} else {
continue_correction <- 1
}
}

}

pool <- pool[-c(pool_to_include[[1]], pool_to_include[[2]])]
pool[[length(pool)+1]] <- new
print(i)

}

max_ages_tree <<- pool[[1]]

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

print(paste(paste("minimum age estimate", i, ""), " initiated", "")) 

###JOIING_TWO_TIPS###
#####################

if ((length(pool[[pool_to_include[[1]]]]$tip.label) == 1) & (length(pool[[pool_to_include[[2]]]]$tip.label) == 1)){
print("joining two tips")
pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length/max_rate
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length/max_rate
if (pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]] > pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]]){
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]] - sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]]
} else {
pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label)]] - sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label)]]	
}

new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]]

###CONGRUIFY_TO_NODE_CALIBRATION#######

if (n_min_constraints != 0){
for (a in 1:length(min_constraints_clade)){
if (setequal(new$tip.label, min_constraints_tip_label[[a]])){
print("node calibration being implemented")
while (true_age < min_constraints_age[[a]]){
true_age <- true_age + calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] + calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] + calibration_implementation_precision
}
}
}
}

###FILTER_OUT_POOR_BRANCHES##########

corresponding_tree <- extract.clade(tree, findMRCA(tree, new$tip.label, "node"))
max_lengths <- vector(mode="numeric", length=0)
for (a in 1:nrow(new[[1]])){
for (b in 1:nrow(corresponding_tree[[1]])){
if (new[[1]][,2][[a]] > length(new$tip.label) & corresponding_tree[[1]][,2][[b]] > length(corresponding_tree$tip.label)){
if (setequal(extract.clade(new, new[[1]][,2][[a]])$tip.label, extract.clade(corresponding_tree, corresponding_tree[[1]][,2][[b]])$tip.label)){
max_lengths <- append(max_lengths, corresponding_tree$edge.length[[b]]/min_rate)
}
}
if (new[[1]][,2][[a]] <= length(new$tip.label) & corresponding_tree[[1]][,2][[b]] <= length(corresponding_tree$tip.label)){
if (new$tip.label[[new[[1]][,2][[a]]]] == corresponding_tree$tip.label[[corresponding_tree[[1]][,2][[b]]]]){
max_lengths <- append(max_lengths, corresponding_tree$edge.length[[b]]/min_rate)
}
}
}
}

continue_correction <- 0
while (continue_correction == 0){
long_branch <- vector(mode="numeric", length=0)
for (a in 1:length(max_lengths)){
if (new$edge.length[[a]] > max_lengths[[a]]){
long_branch <- append(long_branch, a)
}
}
if (length(long_branch) > 0){
long_branch <- long_branch[[1]]
if (new[[1]][,2][[long_branch]] <= length(new$tip.label)){
print(paste(paste("ERROR: INCOMPATIBLE ASSUMPTIONS, MIN AGE ", a, sep=""), " TOO OLD OR MIN RATE TOO HIGH", sep=""))
quit()
}
descendant <- extract.clade(new, new[[1]][,2][[long_branch]])
ancestral <- extract.clade(new, new[[1]][,1][[long_branch]])
max_diff <- tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, descendant$tip.label, "node"))]]/min_rate
keep_reducing <- 0
while (keep_reducing == 0){
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] + calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] + calibration_implementation_precision
new$edge.length[[long_branch]] <- new$edge.length[[long_branch]] - calibration_implementation_precision
if (max(node.depth.edgelength(extract.clade(new, findMRCA(new, descendant$tip.label, "node")))) >= max(node.depth.edgelength(extract.clade(new, findMRCA(new, ancestral$tip.label, "node")))) - max_diff){
keep_reducing <- 1
}
}
} else {
continue_correction <- 1
}
}

}

###JOINING_TWO_CLADES###
########################

else if ((length(pool[[pool_to_include[[1]]]]$tip.label) > 1) & (length(pool[[pool_to_include[[2]]]]$tip.label) > 1)){ 
print("joining two clades")
if ((((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]])/max_rate) + node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]])
 > 
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]])/max_rate) + node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]])){
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]]/max_rate))
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]]) - (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]])) 
} else {
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]]/max_rate))
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]]) - (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]])) 
}

new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]]

###CONGRUIFY_TO_NODE_CALIBRATION#######

if (n_min_constraints != 0){
for (a in 1:length(min_constraints_clade)){
if (setequal(new$tip.label, min_constraints_tip_label[[a]])){
print("node calibration being implemented")
while (true_age < min_constraints_age[[a]]){
true_age <- true_age + calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] + calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] + calibration_implementation_precision
}
}
}
}

###FILTER_OUT_POOR_BRANCHES##########

corresponding_tree <- extract.clade(tree, findMRCA(tree, new$tip.label, "node"))
max_lengths <- vector(mode="numeric", length=0)
for (a in 1:nrow(new[[1]])){
for (b in 1:nrow(corresponding_tree[[1]])){
if (new[[1]][,2][[a]] > length(new$tip.label) & corresponding_tree[[1]][,2][[b]] > length(corresponding_tree$tip.label)){
if (setequal(extract.clade(new, new[[1]][,2][[a]])$tip.label, extract.clade(corresponding_tree, corresponding_tree[[1]][,2][[b]])$tip.label)){
max_lengths <- append(max_lengths, corresponding_tree$edge.length[[b]]/min_rate)
}
}
if (new[[1]][,2][[a]] <= length(new$tip.label) & corresponding_tree[[1]][,2][[b]] <= length(corresponding_tree$tip.label)){
if (new$tip.label[[new[[1]][,2][[a]]]] == corresponding_tree$tip.label[[corresponding_tree[[1]][,2][[b]]]]){
max_lengths <- append(max_lengths, corresponding_tree$edge.length[[b]]/min_rate)
}
}
}
}

continue_correction <- 0
while (continue_correction == 0){
long_branch <- vector(mode="numeric", length=0)
for (a in 1:length(max_lengths)){
if (new$edge.length[[a]] > max_lengths[[a]]){
long_branch <- append(long_branch, a)
}
}
if (length(long_branch) > 0){
long_branch <- long_branch[[1]]
if (new[[1]][,2][[long_branch]] <= length(new$tip.label)){
print(paste(paste("ERROR: INCOMPATIBLE ASSUMPTIONS, MIN AGE ", a, sep=""), " TOO OLD OR MIN RATE TOO HIGH", sep=""))
quit()
}
descendant <- extract.clade(new, new[[1]][,2][[long_branch]])
ancestral <- extract.clade(new, new[[1]][,1][[long_branch]])
max_diff <- tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, descendant$tip.label, "node"))]]/min_rate
keep_reducing <- 0
while (keep_reducing == 0){
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] + calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] + calibration_implementation_precision
new$edge.length[[long_branch]] <- new$edge.length[[long_branch]] - calibration_implementation_precision
if (max(node.depth.edgelength(extract.clade(new, findMRCA(new, descendant$tip.label, "node")))) >= max(node.depth.edgelength(extract.clade(new, findMRCA(new, ancestral$tip.label, "node")))) - max_diff){
keep_reducing <- 1
}
}
} else {
continue_correction <- 1
}
}

} 

###JOINING_ONE_TIP_WITH_A_CLADE###
##################################

else if ((length(pool[[pool_to_include[[1]]]]$tip.label) == 1) & (length(pool[[pool_to_include[[2]]]]$tip.label) > 1)){
print("joining one tip with a clade")
if (((pool[[pool_to_include[[1]]]]$edge.length/max_rate) + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]]) 
>
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]])/max_rate) + node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]])){
pool[[pool_to_include[[1]]]]$edge.length <- pool[[pool_to_include[[1]]]]$edge.length/max_rate
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (pool[[pool_to_include[[1]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]]) - (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]]))
} else {
pool[[pool_to_include[[2]]]] <- addroot(pool[[pool_to_include[[2]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[2]]]]$tip.label, "node"))]]/max_rate))
pool[[pool_to_include[[1]]]]$edge.length <- (node.depth.edgelength(pool[[pool_to_include[[2]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]]) - sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]]
} 

new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]]

###CONGRUIFY_TO_NODE_CALIBRATION#######

if (n_min_constraints != 0){
for (a in 1:length(min_constraints_clade)){
if (setequal(new$tip.label, min_constraints_tip_label[[a]])){
print("node calibration being implemented")
while (true_age < min_constraints_age[[a]]){
true_age <- true_age + calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] + calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] + calibration_implementation_precision
}
}
}
}

###FILTER_OUT_POOR_BRANCHES##########

corresponding_tree <- extract.clade(tree, findMRCA(tree, new$tip.label, "node"))
max_lengths <- vector(mode="numeric", length=0)
for (a in 1:nrow(new[[1]])){
for (b in 1:nrow(corresponding_tree[[1]])){
if (new[[1]][,2][[a]] > length(new$tip.label) & corresponding_tree[[1]][,2][[b]] > length(corresponding_tree$tip.label)){
if (setequal(extract.clade(new, new[[1]][,2][[a]])$tip.label, extract.clade(corresponding_tree, corresponding_tree[[1]][,2][[b]])$tip.label)){
max_lengths <- append(max_lengths, corresponding_tree$edge.length[[b]]/min_rate)
}
}
if (new[[1]][,2][[a]] <= length(new$tip.label) & corresponding_tree[[1]][,2][[b]] <= length(corresponding_tree$tip.label)){
if (new$tip.label[[new[[1]][,2][[a]]]] == corresponding_tree$tip.label[[corresponding_tree[[1]][,2][[b]]]]){
max_lengths <- append(max_lengths, corresponding_tree$edge.length[[b]]/min_rate)
}
}
}
}

continue_correction <- 0
while (continue_correction == 0){
long_branch <- vector(mode="numeric", length=0)
for (a in 1:length(max_lengths)){
if (new$edge.length[[a]] > max_lengths[[a]]){
long_branch <- append(long_branch, a)
}
}
if (length(long_branch) > 0){
long_branch <- long_branch[[1]]
if (new[[1]][,2][[long_branch]] <= length(new$tip.label)){
print(paste(paste("ERROR: INCOMPATIBLE ASSUMPTIONS, MIN AGE ", a, sep=""), " TOO OLD OR MIN RATE TOO HIGH", sep=""))
quit()
}
descendant <- extract.clade(new, new[[1]][,2][[long_branch]])
ancestral <- extract.clade(new, new[[1]][,1][[long_branch]])
max_diff <- tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, descendant$tip.label, "node"))]]/min_rate
keep_reducing <- 0
while (keep_reducing == 0){
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] + calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] + calibration_implementation_precision
new$edge.length[[long_branch]] <- new$edge.length[[long_branch]] - calibration_implementation_precision
if (max(node.depth.edgelength(extract.clade(new, findMRCA(new, descendant$tip.label, "node")))) >= max(node.depth.edgelength(extract.clade(new, findMRCA(new, ancestral$tip.label, "node")))) - max_diff){
keep_reducing <- 1
}
}
} else {
continue_correction <- 1
}
}

}

###ALTERNATE 

else if ((length(pool[[pool_to_include[[2]]]]$tip.label) == 1) & (length(pool[[pool_to_include[[1]]]]$tip.label) > 1)){
print("joining one tip with a clade")
if (((pool[[pool_to_include[[2]]]]$edge.length/max_rate) + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]]) 
>
(((tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]])/max_rate) + node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]])){
pool[[pool_to_include[[2]]]]$edge.length <- pool[[pool_to_include[[2]]]]$edge.length/max_rate
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (pool[[pool_to_include[[2]]]]$edge.length + sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]]) - (node.depth.edgelength(pool[[pool_to_include[[1]]]]$tip.label[[1]]) + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]]))
} else {
pool[[pool_to_include[[1]]]] <- addroot(pool[[pool_to_include[[1]]]], (tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, pool[[pool_to_include[[1]]]]$tip.label, "node"))]]/max_rate))
pool[[pool_to_include[[2]]]]$edge.length <- (node.depth.edgelength(pool[[pool_to_include[[1]]]])[[1]] + sample_time[[which(tip == pool[[pool_to_include[[1]]]]$tip.label[[1]])]]) - sample_time[[which(tip == pool[[pool_to_include[[2]]]]$tip.label[[1]])]]
} 
new <- bind.tree(pool[[pool_to_include[[1]]]], pool[[pool_to_include[[2]]]])
true_age <- node.depth.edgelength(new)[[1]] + sample_time[[which(tip == new$tip.label[[1]])]]

###CONGRUIFY_TO_NODE_CALIBRATION#######

if (n_min_constraints != 0){
for (a in 1:length(min_constraints_clade)){
if (setequal(new$tip.label, min_constraints_tip_label[[a]])){
print("node calibration being implemented")
while (true_age < min_constraints_age[[a]]){
true_age <- true_age + calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[1]]]] + calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] <- new$edge.length[[which(new[[1]][,1] == length(new$tip.label) +1)[[2]]]] + calibration_implementation_precision
}
}
}
}

###FILTER_OUT_POOR_BRANCHES##########

corresponding_tree <- extract.clade(tree, findMRCA(tree, new$tip.label, "node"))
max_lengths <- vector(mode="numeric", length=0)
for (a in 1:nrow(new[[1]])){
for (b in 1:nrow(corresponding_tree[[1]])){
if (new[[1]][,2][[a]] > length(new$tip.label) & corresponding_tree[[1]][,2][[b]] > length(corresponding_tree$tip.label)){
if (setequal(extract.clade(new, new[[1]][,2][[a]])$tip.label, extract.clade(corresponding_tree, corresponding_tree[[1]][,2][[b]])$tip.label)){
max_lengths <- append(max_lengths, corresponding_tree$edge.length[[b]]/min_rate)
}
}
if (new[[1]][,2][[a]] <= length(new$tip.label) & corresponding_tree[[1]][,2][[b]] <= length(corresponding_tree$tip.label)){
if (new$tip.label[[new[[1]][,2][[a]]]] == corresponding_tree$tip.label[[corresponding_tree[[1]][,2][[b]]]]){
max_lengths <- append(max_lengths, corresponding_tree$edge.length[[b]]/min_rate)
}
}
}
}

continue_correction <- 0
while (continue_correction == 0){
long_branch <- vector(mode="numeric", length=0)
for (a in 1:length(max_lengths)){
if (new$edge.length[[a]] > max_lengths[[a]]){
long_branch <- append(long_branch, a)
}
}
if (length(long_branch) > 0){
long_branch <- long_branch[[1]]
if (new[[1]][,2][[long_branch]] <= length(new$tip.label)){
print(paste(paste("ERROR: INCOMPATIBLE ASSUMPTIONS, MIN AGE ", a, sep=""), " TOO OLD OR MIN RATE TOO HIGH", sep=""))
quit()
}
descendant <- extract.clade(new, new[[1]][,2][[long_branch]])
ancestral <- extract.clade(new, new[[1]][,1][[long_branch]])
max_diff <- tree_data_frame[,3][[which(tree_data_frame[,2] == findMRCA(tree, descendant$tip.label, "node"))]]/min_rate
keep_reducing <- 0
while (keep_reducing == 0){
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[1]]]] + calibration_implementation_precision
new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] <- new$edge.length[[which(new[[1]][,1] == findMRCA(new, descendant$tip.label, "node"))[[2]]]] + calibration_implementation_precision
new$edge.length[[long_branch]] <- new$edge.length[[long_branch]] - calibration_implementation_precision
if (max(node.depth.edgelength(extract.clade(new, findMRCA(new, descendant$tip.label, "node")))) >= max(node.depth.edgelength(extract.clade(new, findMRCA(new, ancestral$tip.label, "node")))) - max_diff){
keep_reducing <- 1
}
}
} else {
continue_correction <- 1
}
}


}

pool <- pool[-c(pool_to_include[[1]], pool_to_include[[2]])]
pool[[length(pool)+1]] <- new
print(i)

}

min_ages_tree <<- pool[[1]]

print("MIN AGE TREE CALCULATED")

#################################
###PLOT_BRANCH_RATES#############
#################################

print("EXTRACTING_RATES")

max_ages_tree_data_frame <- data.frame(max_ages_tree[[1]], max_ages_tree$edge.length)
min_ages_tree_data_frame <- data.frame(min_ages_tree[[1]], min_ages_tree$edge.length)

max_ages_tree_branch_rates <- vector(mode="numeric", length=0)
for (j in 1:nrow(tree_data_frame)){
for (i in 1:nrow(max_ages_tree_data_frame)){
if (tree_data_frame[,2][[j]] > length(tree$tip.label) & max_ages_tree_data_frame[,2][[i]] > length(max_ages_tree$tip.label)){
if (setequal(extract.clade(tree, tree_data_frame[,2][[j]])$tip.label, extract.clade(max_ages_tree, max_ages_tree_data_frame[,2][[i]])$tip.label)){
max_ages_tree_branch_rates <- append(max_ages_tree_branch_rates, tree_data_frame[,3][[j]]/max_ages_tree_data_frame[,3][[i]])
}
}
if (tree_data_frame[,2][[j]] < length(tree$tip.label) & max_ages_tree_data_frame[,2][[i]] < length(max_ages_tree$tip.label)){
if (tree$tip.label[[tree_data_frame[,2][[j]]]] == max_ages_tree$tip.label[[max_ages_tree_data_frame[,2][[i]]]]){
max_ages_tree_branch_rates <- append(max_ages_tree_branch_rates, tree_data_frame[,3][[j]]/max_ages_tree_data_frame[,3][[i]])
}
}
}
}

min_ages_tree_branch_rates <- vector(mode="numeric", length=0)
for (j in 1:nrow(tree_data_frame)){
for (i in 1:nrow(min_ages_tree_data_frame)){
if (tree_data_frame[,2][[j]] > length(tree$tip.label) & min_ages_tree_data_frame[,2][[i]] > length(min_ages_tree$tip.label)){
if (setequal(extract.clade(tree, tree_data_frame[,2][[j]])$tip.label, extract.clade(min_ages_tree, min_ages_tree_data_frame[,2][[i]])$tip.label)){
min_ages_tree_branch_rates <- append(min_ages_tree_branch_rates, tree_data_frame[,3][[j]]/min_ages_tree_data_frame[,3][[i]])
}
}
if (tree_data_frame[,2][[j]] < length(tree$tip.label) & min_ages_tree_data_frame[,2][[i]] < length(min_ages_tree$tip.label)){
if (tree$tip.label[[tree_data_frame[,2][[j]]]] == min_ages_tree$tip.label[[min_ages_tree_data_frame[,2][[i]]]]){
min_ages_tree_branch_rates <- append(min_ages_tree_branch_rates, tree_data_frame[,3][[j]]/min_ages_tree_data_frame[,3][[i]])
}
}
}
}

rates <- data.frame(max_ages_tree_branch_rates, min_ages_tree_branch_rates)
write.table(rates, "rates.tsv")
write.tree(min_ages_tree, "min_age_tree.tre")
write.tree(max_ages_tree, "max_age_tree.tre")
}
