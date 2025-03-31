#' Kidney-Related R Functions for Clinical and Epidemiological Research
#' 
#' Package contains different functions for use in the field of kidney disease and general epidemiology
#' (calculation of estimated GFR by different equations, calculation of KDPI and KDRI for kidney transplant donors, etc).
#' 
#' @name kidney.epi 
#'
#' @noRd
#' @keywords internal 
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("The ", pkgname, " package is made with care by the research consultancy Scientific-Tools.Org.\nContact us at https://Scientific-Tools.Org or via 'maintainer(\"", pkgname, "\")' for data analysis or software development."))
}