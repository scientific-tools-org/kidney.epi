#' Calculate eGFR by Schwartz equation (for children only, both "classic" and "quadratic")
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by Schwartz equation.
#'
#' Reference to the equation: Gao A, Cachat F, Faouzi M et al. Comparison of the glomerular filtration rate in children by the new revised Schwartz formula and a new generalized formula. Kidney International 2013;83:524–30.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param age Numeric vector. Age, in years. Age does not accounted in Schwartz equation, but used in the function to check whether Schwartz equation could be applied to a given patient.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female. Required only in case of quadratic Schwartz equation. 
#' @param height_cm Numeric vector. Could be defined either as height_cm if is measured in cm, or as height_ft and height_inch if is measured in feet and inches.
#'    If the parameter height_cm is greater than 0, the function uses cm, otherwise - feet and inches.
#' @param height_ft see height_cm
#' @param height_inch see height_cm
#' @param creatinine_units Character string. Units in which serum creatinne is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param equation_type Character string. Define whether to calculate eGFR either by classic Schwartz or quadratic Schwartz equation. Could be one of the following: "classic", "quadratic". If not explicitly defined by user, the default assumption is "classic".
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.schwartz
#' @examples
#' # for a single patient
#' egfr.schwartz (creatinine = 1.4, age = 10, height_cm = 90, sex = "Male",
#'   creatinine_units = "mg/dl")
#' egfr.schwartz (creatinine = 1.4, age = 10, height_cm = 90, sex = "Male",
#'   creatinine_units = "mg/dl", equation_type = "quadratic")
#' # for a dataset - see vignettes for details
#' # egfr.schwartz (creatinine = dta$scr, age = dta$age, height_cm = dta$ht,
#' #  sex = dta$sex, creatinine_units = "mg/dl")

egfr.schwartz <- function(
  # variables for calculation of eGFR
  creatinine, age, sex, height_cm = 0, height_ft = 0, height_inch = 0, 
  # creatinine measurement units
  creatinine_units = "micromol/l",
  # Schwartz equation type
  equation_type = "classic",
  # custom labels for factor parameters and their unknown values - more explanations are available in the vignette
    # label for definition male sex in data set
    label_sex_male = c ("Male", 1),
    # label for definition female sex in data set
    label_sex_female = c ("Female", 0)
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  # check whether all obligatory argument(s) is(are) defined by user
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  # at the beginning assume that all function arguments are present, and in the following code change it to FALSE if any of the obligatory argument(s) is(are) absent
  fx_params_resulting <- TRUE
  # introduce variable for understanding how user defined height
  local_height_type <- "unknown"
  # Check of height is complex because they could be expressed in different variables by different dimentions
  if ("height_cm" %in% args){
    # do nothing
  local_height_type <- "cm"
  }else{
    # if height_cm not presented check whether height_ft and height_inch are presented
    if( ("height_ft" %in% args) && ("height_inch" %in% args)){
      # do nothing
    local_height_type <- "ft"
    }else{
      # if neither height_cm nor height_ft and height_inch are presented - there is no data for height at all
      warning("Obligatory argument for height is not defined by user in the function arguments", "\n")
      err_num <- err_num + 1
      fx_params_resulting <- FALSE
    }
  }
  
  # List of obligatory function params which have to be defined by user at the function call
  if(equation_type == "classic"){fx_params <- c("creatinine")} 
  if(equation_type == "quadratic"){fx_params <- c("creatinine", "sex")} 
 
  service.check_obligatory_params(fx_params, args, fx_params_resulting)

  # Check the type of arguments inputed by user
  # check that user defined a single creatinine_units
  service.check_param_number(creatinine_units)

  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)

  # check the range of creatinine_units
  service.check_param_arguments(creatinine_units, c("mg/dl", "micromol/l", "mmol/l"))

  # Check whether numerical arguments inputed by user are fine (weight, height etc have to be positive numbers)
  service.check_params_numeric(age, creatinine)
  if(equation_type == "quadratic"){
    if(local_height_type == "cm"){service.check_params_numeric(height_cm)}
    if(local_height_type == "ft"){service.check_params_numeric(height_ft, height_inch)}
  } 
  
  
  # check plausible biologic boundaries (by functions in the service.check_plausibility.R):
  #   check and inform user whether any values out of boundaries were substituted by NA
  #   after the check change the value to boundaries in the possible range (i.e. age > 0 and < 100)
  # age
  # first: general check and tidy: age <0 OR age >100
  age <- service.check_plausibility.age(age)
  # second: since this eGFR equation was developed and validated only for children older 2 years, notify user if younger children were found, and exclude them from calculation
  suspiciosly_low <- service.count_lowerequal_threshold(age, 2)
  if(suspiciosly_low > 0){
    cat(service.output_message(suspiciosly_low, "age <=2 years", "NA"))
    cat("eGFR for children younger than 2 years has to be defined according nomograms\n")
  }
  age <- service.strict_to_numeric_threshold_lower(age, 2)
  # third: since this eGFR equation was developed and validated only for children, notify user if adults were found, and exclude them from calculation
  suspiciosly_high <- service.count_greater_threshold(age, 18)
  if(suspiciosly_high > 0) warning(service.output_message(suspiciosly_high, "age >18 years", "NA"))
  age <- service.strict_to_numeric_threshold_greater(age, 18)

  # creatinine
  service.check_plausibility.creatinine(creatinine)
  
  # I don't check height, since it should be done by user
  

  # CHECK FUNCTION INPUT: END
  ##################################################################


  #
  # repeat creatinine_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  creatinine_units <- rep(creatinine_units, length(creatinine))  
  
  #
  # convert creatinine units if necessary
  #
  creatinine <- service.convert_creatinine(creatinine, creatinine_units)

  #
  # convert height units if necessary
  #
  height <-   ifelse(  height_cm > 0,
            height_cm, # if height is defined in cm, just take this value
            ifelse(height_ft > 0 | height_inch > 0, # if height is defined NOT in cm, check whether height in feets or inches is defined
            height_ft * 30.48 + height_inch * 2.54, # if they contain number than convert to cm
            NA) # if height in feets or inches is not defined, assume NA
        )
  

  # if user required classic equation
  if(equation_type == "classic"){
    eGFR <- 41.3 * (0.01 * height / creatinine )
  }

  # if user required quadratic equation
  if(equation_type == "quadratic"){
  eGFR <- 0.68 * ( height / creatinine ) -0.0008 * (height / creatinine ) * (height / creatinine ) + 0.48 * age
  eGFR <- ifelse( sex %in% label_sex_female,
                        eGFR -25.68,
                        ifelse( sex %in% label_sex_male,
                                eGFR  -21.53,
                                NA
                                )
                )
  
  }


return (round(eGFR, 2))

}


# FUNCTION: END
##################################################################
