#' main input function 
#' enables parameterisation of rmin and rmax (including SetAutoRates) and enables full exTREEmaTIME analysis
#' @export

exTREEmaTIMEmain <- function(tree, 
auto_rates, 
auto_rates_type,
max_rate_unit, 
min_rate_unit, 
rmax,
rmin, 
noise_level,
root_max,
root_min,
n_max_constraints, 
max_constraints_clade, 
max_constraints_ages, 
n_min_constraints, 
min_constraints_clade, 
min_constraints_ages, 
tip, 
sample_time,
output_directory){

for (t in 1:length(n_max_constraints)){

#######################################
###UNSET_RATES#########################
#######################################

if (auto_rates == 1){
print("option 1: unset rates")

###BEGIN_LOOP_THROUGH_TREES############

continue_through_trees <<- 1 
tree_number <<- vector(mode="numeric", length=0)
while(continue_through_trees == 1){
tree_number <<- append(tree_number, 1)
print("doing next tree")

###BEGIN_LOOP_THROUGH_UNITS############

rate_set <<- 0
continue_through_units <<- 1
min_rate_unit_number <<- 1
max_rate_unit_number <<- 1
while(continue_through_units == 1){
print("trying new rate unit")

###RUN_ANALYSIS_WITH_NUMBER_COUNTER###

SetAutoRates(tree[[length(tree_number)]], auto_rates_type, max_rate_unit[[length(max_rate_unit_number)]], min_rate_unit[[length(min_rate_unit_number)]], noise_level, root_min[[t]], root_max[[t]])

if (rate_set == 1){
exTREEmaTIMEcore(tree[[length(tree_number)]], rmax_calc, rmin_calc, noise_level, n_max_constraints[[t]], max_constraints_clade[[t]], max_constraints_ages[[t]], n_min_constraints[[t]], min_constraints_clade[[t]], min_constraints_ages[[t]], tip, sample_time)

if (continue_through_units == 0){
write.table(rates, paste(output_directory, "/", t, "_constraint_config_", length(tree_number), "rates.tsv", sep=""))
write.tree(min_ages_tree, paste(output_directory, "/", t, "_constraint_config_", length(tree_number), "min_age_tree.tre", sep=""))
write.tree(max_ages_tree, paste(output_directory, "/", t, "_constraint_config_", length(tree_number), "max_age_tree.tre", sep=""))
print(paste("CALCULATED FROM TREE", length(tree_number), sep=""))
}

if ((length(max_rate_unit_number) > length(max_rate_unit)) || (length(min_rate_unit_number) > length(min_rate_unit))){
if (continue_through_units == 1){
continue_through_units <<- 0
print("DID NOT REACH COMPLETION WITH SPECIFIED UNITS")
}
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
tree_number <<- vector(mode="numeric", length=0)
while(continue_through_trees == 1){
tree_number <<- append(tree_number, 1)
print("doing next tree")

###SET_UP_RESIDUAL_UNIT_LOOP##########

continue_through_units <<- 1
min_rate_unit_number <<- 1
max_rate_unit_number <<- 1

###RUN_ANALYSIS_WITH_NUMBER_COUNTER###

exTREEmaTIMEcore(tree[[length(tree_number)]], rmax, rmin, noise_level, n_max_constraints[[t]], max_constraints_clade[[t]], max_constraints_ages[[t]], n_min_constraints[[t]], min_constraints_clade[[t]], min_constraints_ages[[t]], tip, sample_time)
write.table(rates, paste(output_directory, "/", t, "_constraint_config_", "rates.tsv", sep=""))
write.tree(min_ages_tree, paste(output_directory, "/", t, "_constraint_config_", "min_age_tree.tre", sep=""))
write.tree(max_ages_tree, paste(output_directory, "/", t, "_constraint_config_", "max_age_tree.tre", sep=""))
print(paste("CALCULATED FROM TREE", length(tree_number), sep=""))

if (continue_through_units == 1){
continue_through_units <<- 0
print("DID NOT REACH COMPLETION WITH SPECIFIED RATE")
}

if (length(tree_number) == length(tree)){
continue_through_trees <<- 0
}

}

}

######################################
######################################
######################################

}

}
