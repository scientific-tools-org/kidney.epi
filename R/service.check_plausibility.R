#' Service functions for data check on biological plausibility and biochemistry conversion which could be applied in any function of the package or externally
#'
##################################################################
# FUNCTION: BEGIN
#' Check and modify if necessary the age values. 
#' @noRd
#' @details Service function which check whether age is in biologically plausible boundaries, shows to user warnings if any, and substitute unplausable values.
#' 
#' @param age Numeric. The value to be checked.
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric Vector with controlled values.
#' @name service.check_plausibility.age
service.check_plausibility.age <- function(age, max_age = 100) {
  # logic for this and other similar functions to check plausible biologic boundaries is the following:
  #   check and inform user whether any values out of boundaries were substituted by NA
  #   after the check change the value to boundaries in the possible range (i.e. age > 0 and < 100)

  suspiciosly_low <- service.count_lowerequal_threshold(age, 0)
  if(suspiciosly_low > 0) warning(service.output_message(suspiciosly_low, "negative values for age", "NA"))
  suspiciosly_high <- service.count_greater_threshold(age, max_age)
  if(suspiciosly_high > 0) warning(service.output_message(suspiciosly_high, paste0("age >", max_age, " years"), "NA"))
  age <- service.strict_to_numeric_threshold_lower(age, 0)
  age <- service.strict_to_numeric_threshold_greater(age, max_age)

return (age)
}
# FUNCTION: END
##################################################################



##################################################################
# FUNCTION: BEGIN
#' Check and modify if necessary the creatinine values. 
#' @noRd
#' @details Service function which check whether creatinine is in biologically plausible boundaries, shows to user warnings if any, and substitute unplausable values.
#' 
#' @param creatinine Numeric. The value to be checked.
#' @return numeric Vector with controlled values.
#' @name service.check_plausibility.creatinine
service.check_plausibility.creatinine <- function(creatinine) {

  suspiciosly_low <- service.count_lowerequal_threshold(creatinine, 0)
  if(suspiciosly_low > 0) warning(service.output_message(suspiciosly_low, "creatinine <=0", "NA"))
  creatinine <- service.strict_to_numeric_threshold_lower(creatinine, 0)

return (creatinine)
}
# FUNCTION: END
##################################################################


##################################################################
# FUNCTION: BEGIN
#' Check and modify if necessary albuminuria values. 
#' @noRd
#' @details Service function which check whether albuminuria is in biologically plausible boundaries, shows to user warnings if any, and substitute unplausable values.
#' 
#' @param albuminuria Numeric. The value to be checked.
#' @return numeric Vector with controlled values.
#' @name service.check_plausibility.albuminuria
service.check_plausibility.albuminuria <- function(albuminuria) {

  suspiciosly_low <- service.count_lower_threshold(albuminuria, 0)
  if(suspiciosly_low > 0) warning(service.output_message(suspiciosly_low, "albuminuria <0", "NA"))
  albuminuria <- service.strict_to_numeric_threshold_lower(albuminuria, 0)

return (albuminuria)
}
# FUNCTION: END
##################################################################




##################################################################
# FUNCTION: BEGIN
#' Check and modify if necessary gfr values. 
#' @noRd
#' @details Service function which check whether gfr is in biologically plausible boundaries, shows to user warnings if any, and substitute unplausable values.
#' 
#' @param gfr Numeric. The value to be checked.
#' @return numeric Vector with controlled values.
#' @name service.check_plausibility.gfr
service.check_plausibility.gfr <- function(gfr) {

  suspiciosly_low <- service.count_lowerequal_threshold(gfr, 0)
  if(suspiciosly_low > 0) warning(service.output_message(suspiciosly_low, "gfr <=0", "NA"))
  gfr <- service.strict_to_numeric_threshold_lower(gfr, 0)

return (gfr)
}
# FUNCTION: END
##################################################################



##################################################################
# FUNCTION: BEGIN
#' Convert creatinine values between measurement units. 
#' @details Service function which check measurement units and convert creatinine values to selected by user.
#' 
#' @param creatinine Numeric. The creatinine value from a data set.
#' @param creatinine_units Character. Creatinine mesurement units in a data set.
#' @param creatinine_reference_units Character. Creatinine measurement units as a desired output (mg/dl by default).
#' @return numeric Creatinine values converted into reference measurement units.
#' @export
#' @name service.convert_creatinine
service.convert_creatinine <- function(creatinine, creatinine_units, creatinine_reference_units = "mg/dl") {

  creatinine_units <- unique(creatinine_units)
  if(creatinine_units != creatinine_reference_units){
    if(creatinine_reference_units == "mg/dl"){
      if(creatinine_units == "micromol/l") creatinine <- creatinine / 88.4 
      if(creatinine_units == "mmol/l") creatinine <- 1000 * creatinine / 88.4
    }else if(creatinine_reference_units == "micromol/l"){
      if(creatinine_units == "mg/dl") creatinine <- creatinine * 88.4
      if(creatinine_units == "mmol/l") creatinine <- 1000 * creatinine
    }else if(creatinine_reference_units == "mmol/l"){
      if(creatinine_units == "mg/dl") creatinine <- (creatinine * 88.4) / 1000
      if(creatinine_units == "micromol/l") creatinine <- creatinine / 1000
    }
  }

return (creatinine)
}
# FUNCTION: END
##################################################################



##################################################################
# FUNCTION: BEGIN
#' Convert cystatin C values between measurement units. 
#' @details Service function which check measurement units and convert cystatin C values to mg/l.
#' 
#' @param cystatin Numeric. The cystatin C values from a data set.
#' @param cystatin_units Character. Cystatin C measurement units in a data set.
#' @return numeric Cystatin C values converted in mg/l.
#' @export
#' @name service.convert_cystatin
service.convert_cystatin <- function(cystatin, cystatin_units) {
  cystatin_units <- unique(cystatin_units)
  # by default cystatin C is measured in mg/l and it is the reference in all kidney.epi functions
  # convert if data are in nanomol/l
  if(cystatin_units == "nanomol/l") cystatin <- cystatin / 74.9

  return (cystatin)
}
# FUNCTION: END
##################################################################



##################################################################
# FUNCTION: BEGIN
#' Check string values in data. 
#' @noRd
#' @details Service function which check whether string values in data are limited to possible values for this variable, shows to user warnings if any, and substitute values out of range.
#' 
#' @param input_values string. The values to be checked.
#' @param possible_values string. Possible values for this variable.
#' @return string Vector with controlled values.
#' @name service.check_string_values
service.check_string_values <- function(input_values, possible_values) {

  values_out_of_range <- input_values[!input_values %in% possible_values]
  
  if(length(values_out_of_range) > 0) warning("Provided data contain ", length(values_out_of_range), " values out of the possible range (", paste(possible_values, sep = "", collapse = ", "), "). They are not included in the analysis. The list of these values: ", paste(unique(values_out_of_range), sep = "", collapse = ", "), ".")
}
# FUNCTION: END
##################################################################

