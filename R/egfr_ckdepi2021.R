#' @title Calculate eGFR by the CKD-EPI 2021 creatinine-based equation
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by the CKD-EPI 2021 creatinine-based equation.
#'
#' Reference to the equation: Inker LA, Eneanya ND, Coresh J, et al. New creatinine- and cystatin C–based equations to estimate GFR without race. N Engl J Med. 2021;385:1737-1749.
#'
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#'
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
#' @name egfr.ckdepi.cr.2021
#' @examples
#' # for a single patient
#' egfr.ckdepi.cr.2021 (creatinine = 1.4, age = 60, sex = "Male", 
#'   creatinine_units = "mg/dl")
#' # for a dataset - see vignettes for details
#' # egfr.ckdepi.cr.2021 (creatinine = dta$scr, age = dta$age, sex = dta$sex, 
#' #  creatinine_units = "mg/dl")

egfr.ckdepi.cr.2021 <- function(
  # variables for calculation of eGFR
  creatinine, age, sex, 
  # creatinine measurement units
  creatinine_units = "micromol/l",
  # custom labels for factor parameters and their unknown values - more explanations are available in the vignette
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
  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)
  # check that user defined a single creatinine_units
  service.check_param_number(creatinine_units)

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
  # convert creatinine units if necessary
  #
  creatinine <- service.convert_creatinine(creatinine, creatinine_units)

  # sex-specific equation coefficients
  k <- rep(NA, length(creatinine))
  k [sex %in% label_sex_female] <- 0.7
  k [sex %in% label_sex_male] <- 0.9
  a <- rep(NA, length(creatinine))
  a [sex %in% label_sex_female] <- -0.241
  a [sex %in% label_sex_male] <- -0.302

  cr_over_k <- creatinine/k
  
  # apply coefficients
  eGFR <- 142 * 
  (pmin(cr_over_k, 1)^a) *
  (pmax(cr_over_k, 1)^(-1.2)) *
  (0.9938^age)
  # apply final sex coefficient
  eGFR <- ifelse( sex %in% label_sex_female, eGFR * 1.012, eGFR)

return (round(eGFR, 2))

}


# FUNCTION: END
##################################################################




#' Calculate eGFR by the CKD-EPI 2021 creatinine-cystatin-based equation
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by the CKD-EPI 2021 creatinine-cystatin-based equation.
#'
#' Reference to the equation: Inker LA, Eneanya ND, Coresh J, et al. New creatinine- and cystatin C–based equations to estimate GFR without race. N Engl J Med. 2021;385:1737-1749.
#'
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#'
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param cystatin Numeric vector. Serum cystatin, could be expressed in "mg/L" or "nanomol/L". Units of measurement should be defined in variable cystatin_units (if not defined explicitly by user, the default value is "mg/L").
#' @param age Numeric vector. Age, in years.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female.
#' @param creatinine_units Character string. Units in which serum creatinne is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param cystatin_units Character string. Units in which serum cystatin is expressed. Could be one of the following: "mg/L" or "nanomol/L"
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.ckdepi.cr_cys.2021
#' @examples
#' # for a single patient
#' egfr.ckdepi.cr_cys.2021 (creatinine = 1.4, cystatin = 0.8, age = 60,
#'    sex = "Male", creatinine_units = "mg/dl")
#' # for a dataset - see vignettes for details
#' # egfr.ckdepi.cr_cys.2021 (creatinine = dta$scr, cystatin = dta$cys,
#' #  age = dta$age, sex = dta$sex, creatinine_units = "mg/dl")

egfr.ckdepi.cr_cys.2021 <- function(
  # variables for calculation of eGFR
  creatinine, cystatin, age, sex, 
  # measurement units
  creatinine_units = "micromol/l", cystatin_units = "mg/L",
  # custom labels for factor parameters and their unknown values - more explanations are available in the vignette
	  # label for definition male sex in data set
	  label_sex_male = c ("Male", 1),
	  # label for definition female sex in data set
	  label_sex_female = c ("Female", 0),
  max_age = 100
) {
  

  
  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("creatinine", "cystatin", "age", "sex") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # check that user defined a single creatinine_units
  service.check_param_number(creatinine_units)
  service.check_param_number(cystatin_units)

  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)
  cystatin_units <- tolower(cystatin_units)
  # check the range of creatinine_units
  service.check_param_arguments(creatinine_units, c("mg/dl", "micromol/l", "mmol/l"))
  service.check_param_arguments(cystatin_units, c("mg/l", "nanomol/l"))

  # Check whether numerical arguments inputed by user are fine (weight, height etc have to be positive numbers)
  service.check_params_numeric(age, creatinine, cystatin)
  
  
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
  service.check_plausibility.creatinine(cystatin)

  # CHECK FUNCTION INPUT: END
  ##################################################################


  #
  # repeat creatinine_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  creatinine_units <- rep(creatinine_units, length(creatinine))
  cystatin_units <- rep(cystatin_units, length(cystatin))
  
  #
  # convert creatinine units if necessary
  #
  creatinine <- service.convert_creatinine(creatinine, creatinine_units)
  cystatin <- service.convert_cystatin(cystatin, cystatin_units)

  # sex-specific equation coefficients
  k <- rep(NA, length(creatinine))
  k [sex %in% label_sex_female] <- 0.7
  k [sex %in% label_sex_male] <- 0.9
  a <- rep(NA, length(creatinine))
  a [sex %in% label_sex_female] <- -0.219
  a [sex %in% label_sex_male] <- -0.144

  cr_over_k <- creatinine / k
  
  # apply coefficients
  eGFR <- 135 * 
    (pmin(cr_over_k, 1)^a) *
    (pmax(cr_over_k, 1)^(-0.544)) *
    (pmin(cystatin / 0.8, 1)^(-0.323)) *
    (pmax(cystatin / 0.8, 1)^(-0.778)) *
    (0.9961^age)
  # apply final sex coefficient
  eGFR <- ifelse( sex %in% label_sex_female, eGFR * 0.963, eGFR)

return (round(eGFR, 2))

}


# FUNCTION: END
##################################################################



#' Alias to the latest eGFR CKD-EPI creatinine-based equation
##################################################################
# FUNCTION: BEGIN
#' @details The function is just an alias to the latest eGFR CKD-EPI creatinine-based equation.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @param ... all arguments for the egfr.ckdepi.cr.2021 function.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.ckdepi.cr
egfr.ckdepi.cr <- function(...) {
  cat("The latest 2021 CKD-EPI creatinine-based equation is used", "\r\n")
  egfr.ckdepi.cr.2021(...)
}
# FUNCTION: END
##################################################################
