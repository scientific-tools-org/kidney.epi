#' Calculate eGFR by the Full age spectrum (FAS) creatinine-based equation
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by the Full age spectrum (FAS) creatinine-based equation.
#'
#' Reference to the equation: Pottel H, Hoste L, Dubourg L et al. An estimating glomerular filtration rate equation for the full age spectrum. Nephrol Dial Transplant 2016; 31:798–806 doi:10.1093/ndt/gfv454.
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
#' @name egfr.fas.cr
#' @examples
#' # for a single patient
#' egfr.fas.cr (creatinine = 1.4, age = 60, sex = "Male", 
#'   creatinine_units = "mg/dl")
#' # for a dataset - see vignettes for details
#' # egfr.fas.cr (creatinine = dta$scr, age = dta$age, sex = dta$sex, 
#' #  creatinine_units = "mg/dl")

egfr.fas.cr <- function(
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
  suspiciosly_low <- service.count_lower_threshold(age, 2)
  if(suspiciosly_low > 0) cat(service.output_message(suspiciosly_low, "age <2 years", "NA"))
  age <- service.strict_to_numeric_threshold_lower(age, 2)
  # creatinine
  service.check_plausibility.creatinine(creatinine)

  # CHECK FUNCTION INPUT: END
  ##################################################################


  #
  # repeat creatinine_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  creatinine_units <- rep(creatinine_units, length(creatinine))
  
  # define q for children in micromol/L, age is the index of value, so q_boys[1] is 23 for 1-year old, q_boys[7] is 39 for 7-years old, etc
  q_boys <- c(23, 26, 27, 30, 34, 36, 39, 41, 43, 45, 47, 50, 52, 54, 64, 69, 72, 75, 78)
  q_girls <- c(23, 26, 27, 30, 34, 36, 39, 41, 43, 45, 47, 50, 52, 54, 57, 59, 61, 61, 62)
  #
  # convert creatinine units if necessary
  #
  creatinine <- service.convert_creatinine(creatinine, creatinine_units, creatinine_reference_units = "micromol/l")
  # no need to convert q_boys and q_girls because creatinine was forced to be in micromol/l

  # q for all ages
  q_cr <- ifelse(sex %in% label_sex_female,
    ifelse(age >= 20,
      62,
      q_girls[round(age, 0)]
    ),
      ifelse(sex %in% label_sex_male,
        ifelse(age >= 20,
          80,
          q_boys[round(age, 0)]
        ),
      NA # If sex is not in any known label
      )
  )
  
  # apply coefficients
  eGFR <- 107.3 / (creatinine / q_cr)
  eGFR <- ifelse(age <= 40, eGFR, eGFR*(0.988^(age-40)))

return (round(eGFR, 2))

}


# FUNCTION: END
##################################################################




#' Calculate eGFR by the Full age spectrum (FAS) cystatin-based equation
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by the Full age spectrum (FAS) cystatin-based equation.
#'
#' Reference to the equation: Pottel H, Delanaye P, Schaeffner E et al. Estimating glomerular filtration rate for the full age spectrum from serum creatinine and cystatin C. Nephrol Dial Transplant 2017; 32:497–507 doi:10.1093/ndt/gfw425.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param cystatin Numeric vector. Serum cystatin, could be expressed in "mg/L" or "nanomol/L". Units of measurement should be defined in variable cystatin_units (if not defined explicitly by user, the default value is "mg/L").
#' @param age Numeric vector. Age, in years.
#' @param cystatin_units Character string. Units in which serum cystatin is expressed. Could be one of the following: "mg/L" or "nanomol/L"
#' @param equation_type Character string. Whether to use "precise" or "simplified" equation.
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.fas.cys
#' @examples
#' # for a single patient
#' egfr.fas.cys (cystatin = 0.8, age = 60)
#' # for a dataset - see vignettes for details
#' # egfr.fas.cys (cystatin = dta$cys, age = dta$age)

egfr.fas.cys <- function(
  # variables for calculation of eGFR
  cystatin, age,
  # measurement units
  cystatin_units = "mg/L",
  equation_type = "precise",
  max_age = 100
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("cystatin", "age") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # check that user defined a single cystatin_units
  service.check_param_number(cystatin_units)

  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  cystatin_units <- tolower(cystatin_units)

  # check the range of cystatin_units
  service.check_param_arguments(cystatin_units, c("mg/l", "nanomol/l") )

  service.check_param_arguments(equation_type, c("precise", "simplified") )

  # Check whether numerical arguments inputed by user are fine (weight, height etc have to be positive numbers)
  service.check_params_numeric(age, cystatin)
  
  
  # check plausible biologic boundaries (by functions in the service.check_plausibility.R):
  #   check and inform user whether any values out of boundaries were substituted by NA
  #   after the check change the value to boundaries in the possible range (i.e. age > 0 and < 100)
  # age
  # first: general check and tidy: age <0 OR age >100
  age <- service.check_plausibility.age(age, max_age)
  # second: since this eGFR equation was developed and validated for adults only, notify user if any children were found, and exclude them from calculation
  suspiciosly_low <- service.count_lower_threshold(age, 2)
  if(suspiciosly_low > 0) cat(service.output_message(suspiciosly_low, "age <2 years", "NA"))
  age <- service.strict_to_numeric_threshold_lower(age, 2)
  service.check_plausibility.creatinine(cystatin)

  # CHECK FUNCTION INPUT: END
  ##################################################################


  #
  # repeat creatinine_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  cystatin_units <- rep(cystatin_units, length(cystatin))
  
  #
  # convert creatinine units if necessary
  #
  cystatin <- service.convert_cystatin(cystatin, cystatin_units)

  # q for all ages 
  # Based on this analysis, and for the sake of simplicity, we fixed QcysC to 0.82 mg/L for all ages <70 years and to 0.95 mg/L for older ages.
  # If we use the linear function QcysC = 0.863 + 0.01704 × (Age – 70) to normalize ScysC in the FAScysC and FAScombi equations, then the performance results (data not shown) are not significantly different than when QcysC = 0.95 is used to normalize ScysC in the FAS equation. 
  if (equation_type == "precise") {
    q_cys <- ifelse(age < 70, 0.82, 0.863 + 0.01704 * (age - 70))
  }
  
  if (equation_type == "simplified") {
    q_cys <- ifelse(age < 70, 0.82, 0.95)
  }
  
  
  # apply coefficients
  eGFR <- 107.3 / (cystatin / q_cys)
  eGFR <- ifelse(age > 40, eGFR * (0.988^(age - 40)), eGFR)

return (round(eGFR, 2))

}


# FUNCTION: END
##################################################################






#' Calculate eGFR by the Full age spectrum (FAS) creatinine-cystatin-based equation
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by the Full age spectrum (FAS) creatinine-cystatin-based equation.
#'
#' Reference to the equation: Pottel H, Delanaye P, Schaeffner E et al. Estimating glomerular filtration rate for the full age spectrum from serum creatinine and cystatin C. Nephrol Dial Transplant 2017; 32:497–507 doi:10.1093/ndt/gfw425.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param cystatin Numeric vector. Serum cystatin, could be expressed in "mg/L" or "nanomol/L". Units of measurement should be defined in variable cystatin_units (if not defined explicitly by user, the default value is "mg/L").
#' @param equation_type Character string. Whether to use "precise" or "simplified" equation.
#' @param age Numeric vector. Age, in years.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female.
#' @param alpha Numeric vector. Alpha coefficient for the combined creatinine-cystatin equation. By default is equal to 0.5.
#' @param creatinine_units Character string. Units in which serum creatinne is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param cystatin_units Character string. Units in which serum cystatin is expressed. Could be one of the following: "mg/L" or "nanomol/L"
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.fas.cr_cys
#' @examples
#' # for a single patient
#' egfr.fas.cr_cys (creatinine = 1.4, cystatin = 0.8, age = 60,
#'    sex = "Male", creatinine_units = "mg/dl")
#' # for a dataset - see vignettes for details
#' # egfr.fas.cr_cys (creatinine = dta$scr, cystatin = dta$cys,
#' #  age = dta$age, sex = dta$sex, creatinine_units = "mg/dl")

egfr.fas.cr_cys <- function(
  # variables for calculation of eGFR
  creatinine, cystatin, age, sex, alpha = 0.5,
  # measurement units
  creatinine_units = "micromol/l", cystatin_units = "mg/L",
  equation_type = "precise",
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

  service.check_param_arguments(equation_type, c("precise", "simplified") )

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
  suspiciosly_low <- service.count_lower_threshold(age, 2)
  if(suspiciosly_low > 0) cat(service.output_message(suspiciosly_low, "age <2 years", "NA"))
  age <- service.strict_to_numeric_threshold_lower(age, 2)
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
  alpha <- rep(alpha, length(cystatin))
  
  #
  # convert creatinine units if necessary
  #
  creatinine <- service.convert_creatinine(creatinine, creatinine_units)
  cystatin <- service.convert_cystatin(cystatin, cystatin_units)

  # define q for children in micromol/L, age is the index of value, so q_boys[1] is 23 for 1-year old, q_boys[7] is 39 for 7-years old, etc
  q_boys <- c(23, 26, 27, 30, 34, 36, 39, 41, 43, 45, 47, 50, 52, 54, 64, 69, 72, 75, 78)
  q_girls <- c(23, 26, 27, 30, 34, 36, 39, 41, 43, 45, 47, 50, 52, 54, 57, 59, 61, 61, 62)

  # convert creatinine units
  creatinine <- service.convert_creatinine(creatinine, creatinine_units, creatinine_reference_units = "micromol/l")
  # no need to convert q_boys and q_girls because creatinine was forced to be in micromol/l

  # q for all ages
  q_cr <- ifelse(sex %in% label_sex_female,
    ifelse(age >= 20,
      62,
      q_girls[round(age, 0)]
    ),
      ifelse(sex %in% label_sex_male,
        ifelse(age >= 20,
          80,
          q_boys[round(age, 0)]
        ),
      NA # If sex is not in any known label
      )
  )
 

  if (equation_type == "precise") {
    q_cys <- ifelse(age < 70, 0.82, 0.863 + 0.01704 * (age - 70))
  }
  
  if (equation_type == "simplified") {
    q_cys <- ifelse(age < 70, 0.82, 0.95)
  }
  
  # apply coefficients
  eGFR <-  1 / ( alpha * (creatinine / q_cr) + (1 - alpha) * (cystatin / q_cys) )
  eGFR <- ifelse(age > 40, eGFR * (0.988^(age - 40)), eGFR)
  eGFR <- 107.3 * eGFR

return (round(eGFR, 2))
}


# FUNCTION: END
##################################################################
