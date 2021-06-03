# exTREEmaTIME
exTREEmaTIME: a simple method for incorporating uncertainty into divergence time estimates

## Installation

requires r packages phangorn, devtools, phytools, phylobase, TreePar

to download use install_github("TomCarr/exTREEmaTIME")

## Functions

### To perform the core divergence time analyses use:

exTREEmaTIMEmain(tree,\
auto_rates,\
auto_rates_type,\
max_rate_unit,\
min_rate_unit,\
rmax,\
rmin,\ 
noise_level,\
root_max,\
root_min,\
n_max_constraints,\ 
max_constraints_clade,\ 
max_constraints_ages,\ 
n_min_constraints,\ 
min_constraints_clade,\ 
min_constraints_ages,\ 
calibration_implementation_precision,\ 
tip,\ 
sample_time,\
output_directory)

#### arguments:

tree: input phylogenies with molecular branch lengths. A list of phylo objects. Can be a single tree, or many trees (eg. bootstrap replicates)

auto_rates: 1 or 0. 1 will use SetAutoRates to define rmax and rmin, 0 will not.

auto_rates_type: 1 or 2. Specifies the equation to use with SetAutoRates. 1 is based on variance of root to tip distance (equation 2 in Carruthers and Scotland (2021), 2 is based on variance in length of terminal sister branches (equation 1 in Carruthers and Scotland (2021)).

max_rate_unit: numeric vector with which to scale calculations of rmax in SetAutoRates (Carruthers and Scotland 2021). Analysis will work from smallest to highest value until divergence time estimation can reach completion. A higher max_rate_unit means a higher rmax is calculated by SetAutoRates.   

min_rate_unit: numeric vector with which to scale calculations of rmin in SetAutoRates (Carruthers and Scotland 2021). Analysis will work from smallest to highest value until divergence time estimation can reach completion. A higher min_rate_unit means a lower rmin is calculated by SetAutoRates.   

rmax: specified value for the maximum substitution rate. This value will be used if auto_rates is set to 0. Set this parameter to 0 if auto_rates is set to 1. 

rmin: specified value for the minimum substitution rate. This value will be used if auto_rates is set to 0. Set this parameter to 0 if auto_rates is set to 1.

n_max_constraints: the number of maximum age constraints specified as a numeric vector. Several configurations for fossil calibrations can be inputted into an analysis, with each value in this vector specifying the number of maximum constraints in each configuration. 

max_constraints_clade: lists containing the taxa for which the most recent common ancestor defines each clade with a maximum age constraint

max_constraints_age: list containing the ages of maximum constraints, in same order as max_constraints_clade 

n_min_constraints: the number of minimum age constraints specified as a numeric vector. Several configurations for fossil calibrations can be inputted into an analysis, with each value in this vector specifying the number of minimum constraints in each configuration. 

min_constraints_clade: lists containing the taxa for which the most recent common ancestor defines each clade with a minimum age constraint

min_constraints_age: list containing the ages of minimum constraints, in same order as min_constraints_clade 

calibration_implementation_precision: numeric value defining the percision with which temporal calibrations are implemented. Smaller values are more precise, but will cause the analysis to run somewhat more slowly. 

tip: a vector of tip labels in the input tree. 

sample_time: a list, with each list containing the minimum and maximum sampling time for each tip label

output_directory: directory in which to write output files

#### outputs:

min_ages_tree, max_ages_tree, a summary of table of substitution rates for each branch
