#' calculates plausible values for rmax and rmin for subsequent use in UncertainTree
#' increase the unit until UncertainTreeCalc reaches completion
#' @export
SetAutoRates <- function(tree, root_min, root_max, unit){

rmax <<- (max(node.depth.edgelength(tree)[seq(1, length(tree$tip.label), 1)])*unit)/(root_min + ((root_max-root_min)/2))
rmin <<- (min(node.depth.edgelength(tree)[seq(1, length(tree$tip.label), 1)])/unit)/(root_min + ((root_max-root_min)/2))

}