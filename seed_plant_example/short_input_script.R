library(phangorn)
library(devtools)
library(phytools)
library(phylobase)
library(TreePar)

tree <- read.tree("rooted.tre")

root_min <- 318 
root_max <- 350

n_max_constraints <- 1
max_constraints_clade <- vector("list", n_max_constraints)
max_constraints_age <- vector("list", n_max_constraints)
max_constraints_clade[[1]] <- c("Acalypha", "Ginkgo")
max_constraints_age[[1]] <- 350

n_min_constraints <- 113
source("fossil_info_two.R")

calibration_implementation_precision <- 0.01
sample_time <- rep(0, length(tree$tip.label))










