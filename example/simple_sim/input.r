output_direc <- c("example_output")

tree <- vector("list", 1)
tree[[1]] <- drop.tip(read.tree("example_tree.tre"), "t1001")

auto_rates <- 0

auto_rates_type <- 0

max_rates_unit <- 0

min_rates_unit <- 0

rmax <- 0.0529547629430849

rmin <- 0.0183985869781416

noise_level <- 0.0000000001

root_max <- 0.9970

root_min <- 0.9969

n_min_constraints <- c(7)
n_max_constraints <- c(7)

max_constraints_clade <- vector("list", length(n_max_constraints))
max_constraints_clade[[1]] <- vector("list", n_max_constraints)
max_constraints_clade[[1]][[1]] <- tree[[1]]$tip.label
max_constraints_clade[[1]][[2]] <- c("t18","t22","t9","t2")
max_constraints_clade[[1]][[3]] <- c("t19","t3","t7","t13","t8","t20","t4","t16","t15")
max_constraints_clade[[1]][[4]] <- c("t20","t4","t16","t15")
max_constraints_clade[[1]][[5]] <- c("t17","t21","t6","t5","t11","t1","t12","t14","t10","t24","t23")
max_constraints_clade[[1]][[6]] <- c("t17","t21")
max_constraints_clade[[1]][[7]] <- c("t11","t1","t12","t14","t10","t24","t23")

max_constraints_ages <- vector("list", length(n_max_constraints))
max_constraints_ages[[1]] <- vector("list", n_max_constraints)
max_constraints_ages[[1]][[1]] <- 1.0967
max_constraints_ages[[1]][[2]] <- 1.0967
max_constraints_ages[[1]][[3]] <- 1.0967
max_constraints_ages[[1]][[4]] <- 1.0967
max_constraints_ages[[1]][[5]] <- 1.0967
max_constraints_ages[[1]][[6]] <- 1.0967
max_constraints_ages[[1]][[7]] <- 1.0967

min_constraints_clade <- vector("list", length(n_min_constraints))
min_constraints_clade[[1]] <- vector("list", n_min_constraints)
min_constraints_clade[[1]][[1]] <- tree[[1]]$tip.label
min_constraints_clade[[1]][[2]] <- c("t18","t22","t9","t2")
min_constraints_clade[[1]][[3]] <- c("t19","t3","t7","t13","t8","t20","t4","t16","t15")
min_constraints_clade[[1]][[4]] <- c("t20","t4","t16","t15")
min_constraints_clade[[1]][[5]] <- c("t17","t21","t6","t5","t11","t1","t12","t14","t10","t24","t23")
min_constraints_clade[[1]][[6]] <- c("t17","t21")
min_constraints_clade[[1]][[7]] <- c("t11","t1","t12","t14","t10","t24","t23")

min_constraints_ages <- vector("list", length(n_min_constraints))
min_constraints_ages[[1]] <- vector("list", n_min_constraints)
min_constraints_ages[[1]][[1]] <- 0.973533933364968
min_constraints_ages[[1]][[2]] <- 0.430614845745634
min_constraints_ages[[1]][[3]] <- 0.071614484241418
min_constraints_ages[[1]][[4]] <- 0.0332609084620804
min_constraints_ages[[1]][[5]] <- 0.125563353695493
min_constraints_ages[[1]][[6]] <- 0.0800184920956898
min_constraints_ages[[1]][[7]] <- 0.0452547472362872

calibration_implementation_precision <- 0.0001

tip <- tree[[1]]$tip.label

sample_time <- vector("list", length(tree[[1]]$tip.label))
for (i in 1:length(sample_time)){
sample_time[[i]] <- rep(0, 2)
}

output_directory <- output_direc[[1]]

exTREEmaTIMEmain(tree, 
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
calibration_implementation_precision, 
tip, 
sample_time,
output_directory) 


   