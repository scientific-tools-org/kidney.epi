#' Calculate estimated glomerular filtration rate (eGFR) by MDRD equation
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by MDRD equation.
#'
#' Reference to the equation: Levey AS, Coresh J, Greene T, et al. Using standardized serum creatinine values in the modification of diet in renal disease study equation for estimating glomerular filtration rate. Annals of Internal Medicine 2006;145:247–54.
#'
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#'
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param age Numeric vector. Age, in years.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female.
#' @param ethnicity Vector. Ethnicity, specify in case of African-American patients. The value of variable refers to the parameter label_afroamerican.
#' @param creatinine_units Character string. Units in which serum creatinne is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param creatinine_method Character string. Creatinine standartisation method in a laboratory. Could be either "IDMS" or "non-IDMS". If not explicitly defined by user, the default assumption is "non-IDMS".
#' @param label_afroamerican List. Label(s) for Afroamerican ethnicity.
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.mdrd4
#' @examples
#' # for a single patient
#' egfr.mdrd4 (creatinine = 1.4, age = 60, sex = "Male", ethnicity = "White",
#'   creatinine_units = "mg/dl")
#' # for a dataset - see vignettes for details
#' # egfr.mdrd4 (creatinine = dta$scr, age = dta$age, sex = dta$sex,
#' #  ethnicity = dta$race, creatinine_units = "mg/dl")


egfr.mdrd4 <- function(
  # variables for calculation of eGFR
  creatinine, age, sex, ethnicity, 
  # creatinine measurement units
  creatinine_units = "micromol/l",
  # creatinine standartisation method in a laboratory: "IDMS" in case of isotope dilution mass spectrometry standartisation, "non-IDMS" in other cases
  creatinine_method = "non-IDMS",
  # custom labels for factor parameters and their unknown values
    # introduce variables' labels which any user can adapt to their own labeling in the data file
    # all labels could be a character, numeric, vector, or list - whatever you prefer
    # you have to change the labels according your data file
    #  for example:
    #    if males in your data file defined as "Male", you have to change below the variable label_sex_male = "Male";
    #    if you use 1 for definition of male sex, you have to change below the variable label_sex_male = 1;
    #    if males in your data file defined either as "Male" or as "MALE" or as "Males" or as "M", you have to change below the variable label_sex_male = c ("Male", "MALE", "Males", "M"), but in this case it is worth to normalize your data dictionary;
    #
    # label for Afroamerican ethnicity
    label_afroamerican = c ("Afroamerican"),
    # label for definition male sex in data set
    label_sex_male = c ("Male", 1),
    # label for definition female sex in data set
    label_sex_female = c ("Female", 0),
  max_age = 100
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("creatinine", "age", "ethnicity", "sex") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # check that user defined a single creatinine_units
  service.check_param_number(creatinine_units)
  # check that user defined a single creatinine_method
  service.check_param_number(creatinine_method)

  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)
  creatinine_method <- tolower(creatinine_method)

  # check the range of creatinine_units
  service.check_param_arguments(creatinine_units, c("mg/dl", "micromol/l", "mmol/l"))
  # check the range of creatinine_method
  # it is completely out of logic and has no explanations, but check for creatinine_method produce error in any case, while the check for creatinine_units works fine
  service.check_param_arguments(creatinine_method, c("idms", "non-idms"))

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
  # repeat creatinine_method
  creatinine_method <- rep(creatinine_method, length(creatinine))  
  
  #
  # convert creatinine units if necessary
  #
  creatinine <- service.convert_creatinine(creatinine, creatinine_units)

  # start eGFR calculation
  eGFR <- creatinine^(-1.154)*age^(-0.203)
  

  # apply sex coefficients
  eGFR <- ifelse( sex %in% label_sex_male,
                eGFR, # males
                ifelse( sex %in% label_sex_female,
                        eGFR * 0.742, # females
                        NA # if sex value is not corresponding neither to male nor female labels)
                       )
                )

  # apply ethnicity coefficient
  eGFR <- ifelse(  ethnicity %in% label_afroamerican,
                eGFR * 1.212,
                eGFR
                )

  # apply creatinine standartisation method coefficient
  eGFR <- ifelse(  creatinine_method == "idms",
                eGFR * 175,
                eGFR * 186
                )


return (round(eGFR, 2))

}


# FUNCTION: END
##################################################################


