# exTREEmaTIME
exTREEmaTIME: a simple method for incorporating uncertainty into divergence time estimates

## Installation

requires r packages phangorn, devtools, phytools, phylobase, TreePar

to download use install_github("TomCarr/exTREEmaTIME")

## Functions

### To perform the core divergence time analyses use:

exTREEmaTIMEmain(tree,\
auto_rates,\
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

#### arguments:

tree: input phylogenies with molecular branch lengths. A list of phylo objects. Can be a single tree, or many trees (eg. bootstrap replicates)

auto_rates: 1 or 0. 1 will use the SetAutoRates function to define rmax and rmin, 0 will not.

auto_rates_type: 1 or 2. Specifies the equation to use with SetAutoRates. 1 is based on variance of root to tip distance, 2 is based on variance of in length of sister branches.

max_rate_unit: 

rmax: the maximum possible substitution rate

rmin: the minimum possible substitution rate

n_max_constraints: the number of specified maximum age constraints

max_constraints_clade: a list containing the taxa that belong to each clade with a maximum age constraint. Each list element contains each set of taxa.

max_constraints_age: a vector containing the age of each maximum constraint, must be in order that corresponds to max_constraints_clade

n_min_constraints: the number of specified minimum age constraints

min_constraints_clade: a list containing the taxa that belong to each clade with a minimum age constraint. Each list element contains each set of taxa.

min_constraints_age: a vector containing the age of each minimum constraint, must be in order that corresponds to min_constraints_clade

calibration_implementation_precision: numeric value that describes the scale of each round of age transformations when rescaling branches. Larger values will make the analysis run more quickly

tip: a vector of all tip labels in the input phylogeny

sample_time: vector of the sampling time for each tip. Must be in the same order as the tip labels in the input phylogeny

#### outputs:

min_ages_tree, max_ages_tree, rates


## To estimate reasonable values for rmin and rmax, as set out in Appendix 1, use:

SetAutoRates(tree, root_min, root_max, unit)

#### arguments:

tree: input phylogeny with branch lengths in units of n

root_min: minimum possible age of the root node

root_max: maximum possible age of the root node

unit: Increases the range between estimated values for rmin and rmax. Sequentially increase this value until analysis reaches completion.

#### outputs:

rmin, rmax
 

