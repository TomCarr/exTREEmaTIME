# exTREEmaTIME
exTREEmaTIME: incorporating uncertainty into divergence time estimates

This package comprises two functions.

# To perform the core divergence time analyses use:

exTREEmaTIMEcalc(tree, rmax, rmin, n_max_constraints, max_constraints_clade, max_constraints_age, n_min_constraints, min_constraints_clade, min_constraints_age, calibration_implementation_precision, sample_time)

ARGUEMENTS:

tree: input phylogeny with branch lengths in units of n

rmax: the maximum possible substitution rate

rmin: the minimum possible substitution rate

n_max_constraints: the number of specified maximum age constraints

max_constraints_clade: a list containing the taxa that belong to each clade with a maximum age constraint. Each list element contains each set of taxa.

max_constraints_age: a list containing the age of each maximum constraint, must be in order that corresponds to max_constraints_clade

n_min_constraints: the number of specified minimum age constraints

min_constraints_clade: a list containing the taxa that belong to each clade with a minimum age constraint. Each list element contains each set of taxa.

min_constraints_age: a list containing the age of each minimum constraint, must be in order that corresponds to min_constraints_clade

calibration_implementation_precision: numeric value that describes the scale of each round of age transformations when rescaling branches. Larger values will make the analysis run more quickly

tip: a list of all tip labels in the input phylogeny

sample_time: list of the sampling time for each tip. Must be in the same order as the tip labels in the input phylogeny

OUTPUTS:

min_ages_tree, max_ages_tree, rates


# To estimate reasonable values for rmin and rmax, as set out in appendix 1, use:

SetAutoRates(tree, root_min, root_max, unit)

tree: input phylogeny with branch lengths in units of n

root_min: minimum possible age of the root node

root_max: maximum possible age of the root node

unit: Increases the range between estimated values for rmin and rmax. Sequentially increase this value until analysis reaches completion.

OUTPUT:

rmin, rmax
 

