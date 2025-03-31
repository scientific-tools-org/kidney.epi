#' Calculate eGFR by the EKFC creatinine-based equation
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by the EKFC creatinine-based equation.
#'
#' References to the equation:  
#' \itemize{
#'   \item Initial creatinine-based equation was reported in Pottel H, Björk J, Courbebaisse M, et al. Development and validation of a modified full age spectrum creatinine-based equation to estimate glomerular filtration rate. a cross-sectional analysis of pooled data. Ann Int Med. 2021;174:183–192 doi:10.7326/M20-4366.
#'   \item Subsequent definition of Q coefficients for African and Black European subjects was reported in Pottel H, Björk J, Rule AD, et al. Cystatin C–based equation to estimate GFR without the inclusion of race and sex. N Engl J Med. 2023;388:333-343 doi: 10.1056/NEJMoa22037.
#' }
#'
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param age Numeric vector. Age, in years.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female.
#' @param ethnicity Vector. Ethnicity. If no ethnicity will be defined, the calculation will use coefficients for White European subjects. Specify ethnicity if a study includes African and Black European subjects, and define the values of variable in the parameter label_african.
#' @param creatinine_units Character string. Units in which serum creatinne is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @param label_african List. Label(s) for African ethnicity.
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.ekfc.cr
#' @examples
#' # for a single patient
#' egfr.ekfc.cr (creatinine = 1.4, age = 60, sex = "Male", 
#'   creatinine_units = "mg/dl")
#' # for a dataset - see vignettes for details
#' # egfr.ekfc.cr (creatinine = dta$scr, age = dta$age, sex = dta$sex, 
#' #  creatinine_units = "mg/dl")

egfr.ekfc.cr <- function(
  # variables for calculation of eGFR
  creatinine, age, sex, ethnicity = NA,
  # creatinine measurement units
  creatinine_units = "micromol/l",
  # custom labels for factor parameters and their unknown values - more eexplanations are available in the vignette
    # label for definition male sex in data set
    label_sex_male = c ("Male", 1),
    # label for definition female sex in data set
    label_sex_female = c ("Female", 0),
    # label for African ethnicity
    label_african = c ("African"),
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
  
  #
  # convert creatinine units if necessary
  #
  # since the equation is developed using micromol/L, use them as a reference
  creatinine <- service.convert_creatinine(creatinine, creatinine_units, creatinine_reference_units = "micromol/l")

  # q coefficient from original publication
  q <- ifelse(
    age <= 25 & sex %in% label_sex_female, 
    exp(3.080 + 0.177 * age - 0.223 * log(age) - 0.00596 * (age^2) + 0.0000686 * (age^3)),
    ifelse(
      age <= 25 & sex %in% label_sex_male,
      exp(3.200 + 0.259 * age - 0.543 * log(age) - 0.00763 * (age^2) + 0.0000790 * (age^3)),
      ifelse(
        age > 25 & sex %in% label_sex_female, 
        62,
        ifelse(
          age > 25 & sex %in% label_sex_male, 
          80,
          NA
        )
      )
    )
  )
  
  # ethnicity coefficients from 2023 publication
  q <- ifelse(
    ethnicity %in% label_african & sex %in% label_sex_female, 
    0.74 * 88.4,
      ifelse(
        ethnicity %in% label_african & sex %in% label_sex_male,
        1.02 * 88.4,
        q
      )
  )

  cr_over_q <- creatinine / q
  
  cr_over_q_degree <- ifelse(cr_over_q < 1, -0.322, -1.132)
  
  age_component <- ifelse(age <= 40, 1, 0.990^(age - 40))
  
  # apply everything
  eGFR <- 107.3 * (cr_over_q^cr_over_q_degree) * age_component

return (round(eGFR, 2))
}


# FUNCTION: END
##################################################################




#' Calculate eGFR by the EKFC cystatin-based equation
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by the EKFC cystatin-based equation.
#'
#' Reference to the equation: Pottel H, Björk J, Rule AD, et al. Cystatin C–based equation to estimate GFR without the inclusion of race and sex. N Engl J Med. 2023;388:333-343 doi: 10.1056/NEJMoa22037.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param cystatin Numeric vector. Serum cystatin, could be expressed in "mg/L" or "nanomol/L". Units of measurement should be defined in variable cystatin_units (if not defined explicitly by user, the default value is "mg/L").
#' @param age Numeric vector. Age, in years.
#' @param cystatin_units Character string. Units in which serum cystatin is expressed. Could be one of the following: "mg/L" or "nanomol/L"
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.ekfc.cys
#' @examples
#' # for a single patient
#' egfr.ekfc.cys (cystatin = 0.8, age = 60)
#' # for a dataset - see vignettes for details
#' # egfr.ekfc.cys (cystatin = dta$cys, age = dta$age)

egfr.ekfc.cys <- function(
  # variables for calculation of eGFR
  cystatin, age,
  # measurement units
  cystatin_units = "mg/L",
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

  cystatin_units <- tolower(cystatin_units)

  # check the range of cystatin_units
  service.check_param_arguments(cystatin_units, c("mg/l", "nanomol/l"))

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
  # creatinine
  service.check_plausibility.creatinine(cystatin)

  # CHECK FUNCTION INPUT: END
  ##################################################################


  #
  # repeat cystatin_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  cystatin_units <- rep(cystatin_units, length(cystatin))
  
  # convert cystatin units if necessary
  #

  cystatin <- service.convert_cystatin(cystatin, cystatin_units)

  # q coefficient
  q <- ifelse(age <= 50, 0.83, 0.83 + 0.005*(age - 50))
  
  cys_over_q <- cystatin/q
  
  cys_over_q_degree <- ifelse(cys_over_q < 1, -0.322, -1.132)
  
  age_component <- ifelse(age <= 40, 1, 0.990^(age-40))
  
  # apply everything
  eGFR <- 107.3 * (cys_over_q^cys_over_q_degree) * age_component

return (round(eGFR, 2))

}


# FUNCTION: END
##################################################################

