exTREEmaTIMEmain <- function(tree, 
auto_rates, 
auto_rates_type,
noise_level,
unit, 
root_max,
root_min,
rmax, 
rmin, 
n_max_constraints, 
max_constraints_clade, 
max_constraints_ages, 
n_min_constraints, 
min_constraints_clade, 
min_constraints_ages, 
calibration_implementation_precision, 
tip, 
sample_time){

#######################################
###UNSET_RATES#########################
#######################################

if (auto_rates == 1){
print("option 1: unset rates")

###BEGIN_LOOP_THROUGH_TREES############

continue_through_trees <<- 1 
#tree_number <- rep(1, 76)
tree_number <- vector(mode="numeric", length=0)
while(continue_through_trees == 1){
tree_number <- append(tree_number, 1)
print("doing next tree")

###BEGIN_LOOP_THROUGH_UNITS############

rate_set <<- 0
continue_through_units <<- 1
unit_number <- vector(mode="numeric", length=0)
while(continue_through_units == 1){
unit_number <- append(unit_number, 1)
print("trying new rate unit")

###RUN_ANALYSIS_WITH_NUMBER_COUNTER###

SetAutoRates(tree[[length(tree_number)]], root_min, root_max, unit[[length(unit_number)]], auto_rates_type, noise_level)

if (rate_set == 1){
exTREEmaTIMEcore(tree[[length(tree_number)]], rmax_calc, rmin_calc, n_max_constraints, max_constraints_clade, max_constraints_ages, n_min_constraints, min_constraints_clade, min_constraints_ages, calibration_implementation_precision, tip, sample_time)

if (continue_through_units == 0){
write.table(rates, paste(length(tree_number), "rates.tsv", sep=""))
write.tree(min_ages_tree, paste(length(tree_number), "min_age_tree.tre", sep=""))
write.tree(max_ages_tree, paste(length(tree_number), "max_age_tree.tre", sep=""))
print(paste("CALCULATED FROM TREE", length(tree_number), sep=""))
}

if ((length(unit_number) == length(unit)) & (continue_through_units == 1)){
continue_through_trees <<- 0
stop("DID NOT REACH COMPLETION WITH SPECIFIED UNITS")
}

} else {
continue_through_units <<- 0
}

}

if (length(tree_number) == length(tree)){
continue_through_trees <<- 0
}

}

}

######################################
######################################
######################################

######################################
###SET_RATES##########################
######################################

if (auto_rates == 0){
print("option_2: set rates")

###BEGIN_LOOP_THROUGH_TREES###########

continue_through_trees <<- 1
tree_number <- vector(mode="numeric", length=0)
while(continue_through_trees == 1){
tree_number <- append(tree_number, 1)
print("doing next tree")

###SET_UP_RESIDUAL_UNIT_LOOP##########

continue_through_units <<- 1

###RUN_ANALYSIS_WITH_NUMBER_COUNTER###

exTREEmaTIMEcore(tree[[length(tree_number)]], rmax, rmin, n_max_constraints, max_constraints_clade, max_constraints_ages, n_min_constraints, min_constraints_clade, min_constraints_ages, calibration_implementation_precision, tip, sample_time)
write.table(rates, paste(length(tree_number), "rates.tsv", sep=""))
write.tree(min_ages_tree, paste(length(tree_number), "min_age_tree.tre", sep=""))
write.tree(max_ages_tree, paste(length(tree_number), "max_age_tree.tre", sep=""))
print(paste("CALCULATED FROM TREE", length(tree_number), sep=""))

if (continue_through_units == 1){
continue_through_trees <<- 0
stop("DID NOT REACH COMPLETION WITH SPECIFIED RATE")
}

if (length(tree_number) == length(tree)){
continue_through_trees <<- 0
}

}

}

}