#' Service functions for data check which could be applied in any function of the package or externally
#'
#' @noRd
#' @details Verifies whether the single value is among the values of the vector. 
#' Function is useful to check whether the argument of the function defined by the user is among the possible arguments recognized inside the function.
#' 
#' @param param2check Numeric value or character string. The single value to be verified.
#' @param possible_params Vector. The vector of values which contains all possible values.
#' @return logic returns TRUE if argument param2check is foundin possible values possible_params, and FALSE if it is not.
#' @name service.is.param_possible
#' @examples
#' possible_params = c("KDPI", " KDRI_Rao", "KDRI_median")
#' service.is.param_possible("KDZO", possible_params) # return FALSE
#' service.is.param_possible("KDPI", possible_params) # return TRUE

##################################################################
# FUNCTION: BEGIN
service.is.param_possible <- function(param2check, possible_params) {
  if(param2check %in% possible_params){
    # OK, given value is among possible values
    return(TRUE)
  }else{
    # given value is NOT among possible values
    return(FALSE)
  }
}
# FUNCTION: END
##################################################################





#' Check the argument of a given parameter which set by user and stop function if the value set by user is not among the possible values of the argument 
#' @noRd
#' @details Check the argument of a given parameter which set by user and stop function if the value set by user is not among the possible values of the argument.
#' Service function that will not be exported to user.
#' 
#' @param param2check List, Character string, Number. Parameter used in a function.
#' @param possible_params List. List of possible values of the parameter arguments
#' @param custom_message Character string. Custom message to be output. If not defined, the standart output message is provided.
#' nothing to return
#' @name service.check_param_arguments
# @examples
# service.check_param_arguments(creatinine_method, possible_params)
# 
##################################################################
# FUNCTION: BEGIN
service.check_param_arguments <- function(param2check, possible_params, custom_message = "") {

  # set possible value to lowercase
  possible_params <- tolower(possible_params)
  
  if (!service.is.param_possible(tolower(param2check), possible_params)){
    if( length(custom_message) > 0){
        warning(custom_message)
    }
    # use deparse and substitute to get the name of a function argument
    stop("The defined by user value '", param2check, "' for parameter '", deparse(substitute(param2check)), "' is not among possible values of the parameter. ", "The execution of the function is interrupted.", "\n")
  }  

}
# FUNCTION: END
##################################################################




#' Select only numeric values greater than defined threshold.
#' @noRd
#' @details Select only numeric values greater than defined threshold, and substitute other values with NA. 
#'
#' @param x the vector to be checked.
#' @param threshold numeric the threshold to compare with.
#' @return numeric returns only numeric values greater than threshold.
#' @name service.strict_to_numeric_threshold_lower
#' @examples
#' myvals <- c(1, 8, -5, "oggi", NA)
#' # return to myvals2 only numeric values greater than defined threshold (0 in this case)
#' #    and susbstitute non-numeric or negative values with NA
#' myvals2 <- service.strict_to_numeric_threshold_lower(myvals, 0)
#' myvals2 # 1, 8, NA, NA, NA
#' 
##################################################################
# FUNCTION: BEGIN
service.strict_to_numeric_threshold_lower <- function(x, threshold) {
  y <- ifelse( is.numeric(x) & x >= threshold,
              x,
              NA)
return(y)
}
# FUNCTION: END
##################################################################



#' Select only numeric values lower than defined threshold
#' @noRd
#' @details Select only numeric values lower than defined threshold, and substitute other values with NA. 
#' 
#' @param x the vector to be checked.
#' @param threshold numeric the threshold to compare with.
#' @return numeric returns only numeric values lower than threshold.
#' @name service.strict_to_numeric_threshold_greater
#' @examples
#' myvals <- c(1, 8, -5, "oggi", NA)
#' # ruturn to myvals2 only numeric values lower than threshold  (3 in this case)
#' #   susbstitute non-numeric or negative values with NA
#' myvals2 <- service.strict_to_numeric_threshold_greater(myvals, 3)
#' myvals2 # 1, NA, -5, NA, NA
#' 
##################################################################
# FUNCTION: BEGIN
service.strict_to_numeric_threshold_greater <- function(x, threshold) {
  y <- ifelse( is.numeric(x) & x < threshold,
              x,
              NA)
return(y)
}
# FUNCTION: END
##################################################################




#' Count how many values are less or equal than the defined threshold.
#' @noRd
#' @details Count how many values are less or equal than the defined threshold. 
#' 
#' @param x the vector to be checked.
#' @param threshold numeric the threshold to compare with.
#' @return numeric returns number of numeric values less or equal to the threshold.
#' @name service.count_lowerequal_threshold
#' @examples
#' myvals <- c(1, 8, -5, "oggi", NA)
#' myvals2 <- service.count_lowerequal_threshold(myvals, 0)
#' myvals2 # 1
##################################################################
# FUNCTION: BEGIN
service.count_lowerequal_threshold <- function(x, threshold) {
  mycounter <- sum (is.numeric(x) & x <= threshold & !is.na(x))
return (mycounter)
}
# FUNCTION: END
##################################################################


#' Count how many values are less than the defined threshold.
#' @noRd
#' @details Count how many values are less than the defined threshold. 
#' 
#' @param x the vector to be checked.
#' @param threshold numeric the threshold to compare with.
#' @return numeric returns number of numeric values less the threshold.
#' @name service.count_lower_threshold
#' @examples
#' myvals <- c(1, 8, -5, "oggi", NA)
#' myvals2 <- service.count_lower_threshold(myvals, 0)
#' myvals2 # 1
##################################################################
# FUNCTION: BEGIN
service.count_lower_threshold <- function(x, threshold) {
  mycounter <- sum (is.numeric(x) & x < threshold & !is.na(x))
  return (mycounter)
}
# FUNCTION: END
##################################################################


#' Count how many values are greater than the defined threshold. 
#' @noRd
#' @details Count how many values are greater than the defined threshold. 
#' 
#' @param x the vector to be checked.
#' @param threshold numeric the threshold to compare with.
#' @return numeric returns number of numeric values greater or equal to the threshold.
#' @name service.count_greater_threshold
#' @examples
#' myvals <- c(1, 8, -5, "oggi", NA)
#' myvals2 <- service.count_greater_threshold(myvals, 0)
#' myvals2 # 2
##################################################################
# FUNCTION: BEGIN
service.count_greater_threshold <- function(x, threshold) {
  mycounter <- sum (is.numeric(x) & x > threshold & !is.na(x))
return (mycounter)
}
# FUNCTION: END
##################################################################




#' Check whether a vector is numeric.
#' @noRd
#' @details Check whether a vector is numeric. 
#' 
#' @param x the vector to be checked.
#' @return logic whether vector x is numeric or not.
#' @name service.is_numeric
# @examples
# myvals <- c(1, 8, -5, "oggi", NA)
# service.is_numeric(myvals) # FALSE
##################################################################
# FUNCTION: BEGIN
service.is_numeric <- function(x) {
  return(is.numeric(x))
}
# FUNCTION: END
##################################################################




#' Form output message in singular or plural.
#' @noRd
#' @details Provide different output for constructing messages in singular or plural. 
#' 
#' @param x Numeric. The value to be checked (usualy a counter of some variable).
#' @param singular Character string. The value to be returned in case of singular form (usualy a string, but could be any type).
#' @param plural Character string. The value to be returned in case of plural form (usualy a string, but could be any type).
#' @return Character string. Returns a value for constructing messages in singular or plural form.
#' @name service.singular_or_plural
#' @examples
#' service.singular_or_plural(1, "This value was", "These values were") # "This value was"
#' service.singular_or_plural(99, "This value was", "These values were") # "These values were"
##################################################################
# FUNCTION: BEGIN
service.singular_or_plural <- function(x, singular, plural) {
  if(x == 1){
    return (singular)
  }else{
    return (plural)
  }
}
# FUNCTION: END
##################################################################





#' Produce message for warning or cat
#' @noRd
#' @details Produce message that is used by warning or cat in the ktx.kdpi.optn function. 
#' Service function that will not be exported to user, and used only in the ktx.kdpi.optn function.
#' 
#' @param x Numeric. The value to be checked (usualy a counter of some variable).
#' @param custom_phrase Character string. Custom message to be inserted in the middle of standard message.
#' @param warning_type Character string. The type of message: warning (with substitution to NA) or cat (with leave as is).
#' @return Character string. Returns a phrase.
#' @name service.output_message
# @examples
# service.output_message(suspiciosly_high, "age >100 years", "NA")
# 
##################################################################
# FUNCTION: BEGIN
service.output_message <- function(x, custom_phrase, warning_type) {
  if(warning_type == "NA"){
    last_sentence = paste(" ", service.singular_or_plural(x, "This value was", "These values were"), " substituted to NA.", sep = "")
  }else if(warning_type == "as is"){
    last_sentence = paste(" ", service.singular_or_plural(x, "This value was", "These values were"), " kept as is.", sep = "")
  }
  
  whole_phrase = paste("There", service.singular_or_plural(x, " is ", " are "), x, " patient", service.singular_or_plural(x, "", "s"), " with ", custom_phrase, ". ", last_sentence, "\n", sep = "")
  
return(whole_phrase)
}
# FUNCTION: END
##################################################################






#' Check whether all obligatory paramenters of a given function are present.
#' @noRd
#' @details Check whether all obligatory paramenters of a given function are present. 
#' 
#' @param fx_params List. List of parameters required by function.
#' @param args List. Arguments transferred to the function upon user call.
#' @param predefined_result Logical. Required only in case if other checks were performed in the main script and the result of this check has to be processed to the function.
#'    For example, if in the parent script I've checked the presence of height parameter, and it is absent (while is obligatory), I transfer this info in the "predefined_result = FALSE", so in the function the fx_params_resulting become False and will lead to stop().
#' @return Character string. Returns a messages and stops function if any of the obligatory parameters are absent.
#' @name service.check_obligatory_params
#' @examples
#' # could be run only inside function wich receives some parameters 
#' # fx_params <- c("creatinine", "age", "ethnicity", "sex")
#' # args <- names(as.list(match.call())[-1])
#' # service.check_obligatory_params(fx_params, args)
##################################################################
# FUNCTION: BEGIN
service.check_obligatory_params <- function(fx_params, args, predefined_result = TRUE) {
  # at the beginning assume that all function arguments are present, and in the following code change it to FALSE if any of the obligatory argument(s) is(are) absent
  fx_params_resulting <- TRUE
  err_num <- 0
  
  # if the predefined_result is already suggesting about a missing parameter revealed in the parent script
  if(predefined_result == FALSE){
    fx_params_resulting <- FALSE
    err_num <- 1
  }
  
  # Check whether all necessary params are defined by user
  for(i in 1:length(fx_params)){
    fx_param_local <- fx_params[i] %in% args # whether param is found among the obligatory 
    if( fx_param_local == FALSE ){
      warning("Obligatory argument ", fx_params[i], " is not defined by user in the function arguments", "\n")
      err_num <- err_num + 1
    }
    fx_params_resulting <- fx_params_resulting && fx_param_local
  }
  
  # final message if any of the Obligatory arguments are missing
  if(fx_params_resulting == FALSE){
    stop("Obligatory argument", service.singular_or_plural(err_num, "", "s"), service.singular_or_plural(err_num, " is ", " are "), " not defined by user", "\n", "The execution of the function is interrupted.", "\n")
  }
  
  #return TRUE
}
# FUNCTION: END
##################################################################








#' Check number of parameters and stop function if it exceeds the expected number of parameters 
#' @noRd
#' @details Check number of parameters and stop function if it exceeds the expected number of parameters.
#' Service function that will not be exported to user.
#' 
#' @param param2check List, Character string, Number. Parameter used in a function.
#' @param acceptable_number Numeric. Acceptable number of arguments in the list param2check (by default is "1")
#' @param custom_message Character string. Custom message to be output. If not defined, the standart output message is provided.
#' nothing to return
#' @name service.check_param_number
# @examples
# service.check_param_number(creatinine_method)
# service.check_param_number(param2check = creatinine_method, acceptable_number = 1)
# 
##################################################################
# FUNCTION: BEGIN
service.check_param_number <- function(param2check, acceptable_number = 1, custom_message = "") {

  if( length(param2check) != acceptable_number ){
    if( length(custom_message) > 0){
        warning(custom_message)
    }
	# use deparse and substitute to get the name of a function argument
    stop("The value for '", deparse(substitute(param2check)), "' has to be a single character string. ", "The execution of the function is interrupted.", "\n")
  }
}
# FUNCTION: END
##################################################################












#' Check whether the following variables are numeric and stop function if at least one of them is not numeric 
#' @noRd
#' @details Check whether the following variables are numeric and stop function if at least one of them is not numeric.
#' Service function that will not be exported to user.
#' 
#' @param ... Argument list. Argument list (arbitruary number of valiables) with data to check.
#' nothing to return
#' @name service.check_params_numeric
# @examples
# service.check_params_numeric(creatinine, albumin)
# 
##################################################################
# FUNCTION: BEGIN
service.check_params_numeric <- function(...) {

  # check and inform user whether argument contains non-numeric values
  numeric_resulting = TRUE # assume that all arguments are numeric

  # construct object from the variables defined by user
  # to avoid conversion of numeric values to factor necessary to use rbind.data.frame and cbind.data.frame instead of rbind and cbind
  dta <- cbind.data.frame(...)
  # define number of variables
  n <- length(dta)

  for(i in 1:n){
    # things are rather complicated
    # dta is a list which consists from several variables. to get all values of a single variable I need [[]], but the name of the variable should be taken from the list (so use [] only)
    values2check <- dta[[i]]
    varname <- names(dta[i])
    numeric_local <- service.is_numeric(values2check)
    # numeric_local <- !is.na(as.numeric(values2check)) - doesn't work good
    # numeric_local <- sapply(values2check, is.numeric) - doesn't work good
    if(!numeric_local) warning(varname, " is non-numeric argument.", "\n")
    numeric_resulting <- numeric_local && numeric_resulting # if any of the argument is non-numeric, the numeric_resulting become FALSE
  }

  # resulting message for numeric check
  if(!numeric_resulting){
    stop("At least one of the defined by user arguments is not numeric.", " The execution of the function is interrupted.", "\n")
  }
}
# FUNCTION: END
##################################################################



#' Round decimal for nice output
#' @noRd
#' @details Rounds decimal for nice output. 
#' 
#' Programming: Boris Bikbov \email{boris.bikbov@scientific-tools.org}.
#'
#' @param value numeric. The single value to be formatted.
#' @param rounding_factor numeric. Rounding factor for percentages.
#' @return numeric. Returns formatted decimal.
#' @name service.str.round

##################################################################
# FUNCTION: BEGIN
service.str.round <- function(value, rounding_factor = 1) {
  return(trimws(format(round(value, digits=rounding_factor), nsmall = rounding_factor)))
}


