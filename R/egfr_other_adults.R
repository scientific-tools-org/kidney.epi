#' Calculate eGFR by the revised Lund-Malmö creatinine-based equation
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by the revised Lund-Malmö creatinine-based equation.
#'
#' Reference to the equation: Björk J, Grubb A, Sterner G, Nyman U. Revised equations for estimating glomerular filtration rate based on the Lund-Malmö Study cohort. Scand J Clin Lab Invest. 2011;71: 232-239.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param age Numeric vector. Age, in years.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female.
#' @param creatinine_units Character string. Units in which serum creatinne is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.lm.cr
#' @examples
#' # for a single patient
#' egfr.lm.cr (creatinine = 1.4, age = 60, sex = "Male", 
#'   creatinine_units = "mg/dl")
#' # for a dataset - see vignettes for details
#' # egfr.lm.cr (creatinine = dta$scr, age = dta$age, sex = dta$sex, 
#' #  creatinine_units = "mg/dl")

egfr.lm.cr <- function(
  # variables for calculation of eGFR
  creatinine, age, sex,
  # creatinine measurement units
  creatinine_units = "micromol/l",
  # custom labels for factor parameters and their unknown values - more eexplanations are available in the vignette
    # label for definition male sex in data set
    label_sex_male = c ("Male", 1),
    # label for definition female sex in data set
    label_sex_female = c ("Female", 0),
  max_age = 100
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("creatinine", "age", "sex") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # check that user defined a single creatinine_units
  service.check_param_number(creatinine_units)

  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)

  # check the range of creatinine_units
  service.check_param_arguments(creatinine_units, c("mg/dl", "micromol/l", "mmol/l"))

  # Check whether numerical arguments inputed by user are fine (weight, height etc have to be positive numbers)
  service.check_params_numeric(age, creatinine)
  
  
  # check plausible biologic boundaries (by functions in the service.check_plausibility.R):
  #   check and inform user whether any values out of boundaries were substituted by NA
  #   after the check change the value to boundaries in the possible range (i.e. age > 0 and < 100)
  # age
  # first: general check and tidy: age <0 OR age >100
  age <- service.check_plausibility.age(age, max_age)
  # second: since this eGFR equation was developed and validated for adults only, notify user if any children were found, and exclude them from calculation
  suspiciosly_low <- service.count_lower_threshold(age, 18)
  if(suspiciosly_low > 0) cat(service.output_message(suspiciosly_low, "age <18 years", "NA"))
  age <- service.strict_to_numeric_threshold_lower(age, 18)
  # creatinine
  service.check_plausibility.creatinine(creatinine)

  # CHECK FUNCTION INPUT: END
  ##################################################################


  #
  # repeat creatinine_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  creatinine_units <- rep(creatinine_units, length(creatinine))
  
  #
  # convert creatinine units
  # since the equation is developed using micromol/L, use them as a reference
  creatinine <- service.convert_creatinine(creatinine, creatinine_units, creatinine_reference_units = "micromol/l")

  # x coefficient
  x <- ifelse(sex %in% label_sex_female, 
    ifelse(creatinine < 150, 2.5 + 0.0121 * (150 - creatinine), 2.5 - 0.926 * log(creatinine / 150)),
  ifelse(sex %in% label_sex_male, 
    ifelse(creatinine < 180, 2.56 + 0.00968 * (180 - creatinine), 2.56 - 0.926 * log(creatinine / 180)), NA
  )
  )

  eGFR <- exp(x - 0.0158 * age + 0.438 * log(age))

return (round(eGFR, 2))

}


# FUNCTION: END
##################################################################

