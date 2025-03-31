#' Calculate KDRI and KDPI for deceased kidney donor
#'
#' @details Calculate Kidney Donor Risk Index (KDRI) and Kidney Donor Profile Index (KDPI) based on the algorithm of US Organ Procurement and Transplantation Network. 
#' The Kidney Donor Profile Index (KDPI) is a numerical measure that combines ten donor factors to summarize into a single number the quality of deceased donor kidneys relative to other recovered kidneys. 
#' \emph{KDRI could be calculated only for a deceased donor}!
#' 
#' More reading:
#' \itemize{
#'   \item \href{https://optn.transplant.hrsa.gov/resources/allocation-calculators/kdpi-calculator/}{OPTN web-based calculator}
#'   \item \href{https://optn.transplant.hrsa.gov/data/allocation-calculators/kdpi-calculator/learn-about-kdpi/}{Guide to calculating and interpreting KDPI}
#'   \item \href{https://optn.transplant.hrsa.gov/media/2150/kdpi_mapping_table.pdf}{Latest data for mapping table, scaling factor, etc}
#' }
#' 
#' Citation: Bikbov B. R open source programming code for calculation of the Kidney Donor Profile Index and Kidney Donor Risk Index. Kidney Diseases, 2018; 4:269–272. DOI: 10.1159/000492427
#'
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#'
#' @param age Numeric vector. Age, in years.
#' @param height_cm Numeric vector. Could be defined either as height_cm if is measured in cm, or as height_ft and height_inch if is measured in feet and inches.
#'    If the parameter height_cm is greater than 0, the function uses cm, otherwise - feet and inches.
#' @param height_ft see height_cm
#' @param height_inch see height_cm
#' @param weight_kg Numeric vector. Could be defined either as weight_kg if measured in kg, or as weight_lb if is measured in pounds.
#'    If the parameter weight_kg is greater than 0, the function uses kg, otherwise - pounds.
#' @param weight_lb see weight_kg
#' @param ethnicity Vector. Ethnicity, specify in case of African-American donors which have special coefficient in the regression equation. The value of variable refers to the parameter label_afroamerican.
#' @param hypertension Vector. History of hypertension, specify in case of hypertensive donors which have special coefficient in the regression equation. The value of variable refers to the parameters label_hypertension_positive and label_hypertension_unknown.
#' @param diabetes Vector. History of diabetes, specify in case of donors with diabetes which have special coefficient in the regression equation. The value of variable refers to the parameters label_diabetes_positive and label_diabetes_unknown.
#' @param causeofdeath Vector. Cause of death, specify whether death was due to cerebrovascular disease, or other reasons.
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param hcv Vector. Hepatitis C virus status. The value of variable refers to the parameters label_hcv_positive and label_hcv_unknown.
#' @param dcdstatus Vector. Donation after circulatory death status. Specify whether organ was from a donor after circulatory death or not. The value of variable refers to the parameter label_dcdstatus.
#' @param creatinine_units Character string. Units in which serum creatinne is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param return_output_type Character string. Specify which calculated parameter to return from the function: "KDRI_Rao" - Raw Kidney Donor Risk Index, "KDRI_median" - scaled to the median Kidney Donor Risk Index, or "KDPI" - Kidney Donor Profile Index.
#' @param mapping_values_year Numeric value or character string.  Specify which year to take for the OPTN mapping table, as well as KDRI scaling factor and chances of hypertension and diabetes in case if they were unknown for donor.
#'    
#'    By default the value is "latest", and the function takes the latest available OPTN mapping table and coefficients from the internal dataframes ktx.kdpi_mapping_table and ktx.kdpi_coefficients_table.
#'    
#'    But if necessary, a user could define the exact year (i.e. mapping_values_year = 2015).
#'    
#'    For a list of available years run the following: ktx.kdpi.optn.show.years().
#' @param label_afroamerican List. Label(s) for Afroamerican ethnicity.
#' @param label_hypertension_positive List. Label(s) for a positive history of hypertension.
#' @param label_hypertension_unknown List. Label(s) for donors with unknown history of hypertension.
#' @param label_diabetes_positive List. Label(s) for a positive history of diabetes.
#' @param label_diabetes_unknown List. Label(s) for donors with unknown history of diabetes.
#' @param label_causeofdeath List. Label(s) for a cause of death due to cerebrovascular/stroke.
#' @param label_hcv_positive List. Label(s) for a positive HCV status.
#' @param label_hcv_unknown List. Label(s) for an unknown, not done, indeterminate, or pending HCV status.
#' @param label_dcdstatus List. Label(s) for a donor after circulatory death status.
#' @return numeric One of the following values based on the return_output_type argument: Raw Kidney Donor Risk Index (KDRI), Scaled to the median Kidney Donor Risk Index (KDRI), or Kidney Donor Profile Index (KDPI).
#' @export
#' @name ktx.kdpi.optn
#' @examples
#' # for a single patient
#' ktx.kdpi.optn (age = 60, height_cm = 168, weight_kg = 93, ethnicity = "White",
#'   hypertension = "yes", diabetes = "no", causeofdeath = "roadinjury",
#'   creatinine = 1.4, hcv = "negative", dcdstatus = "no",
#'   creatinine_units = "mg/dl", return_output_type = "KDRI_Rao")
#' ktx.kdpi.optn (age = 30, height_cm = 176, weight_kg = 82, ethnicity = "White", 
#'   hypertension = "NA", diabetes = "no", causeofdeath = "roadinjury", # note that NA is in the quotes
#'   creatinine = 150, hcv = "negative", dcdstatus = "no", return_output_type = "KDPI")
#' # for a dataset - see vignettes for details

##################################################################
# FUNCTION: BEGIN

ktx.kdpi.optn <- function(
  # variables for calculation of KDPI/KDRI
  age, height_cm = 0, height_ft = 0, height_inch = 0, weight_kg = 0, weight_lb = 0, ethnicity, hypertension, diabetes, causeofdeath, creatinine, hcv, dcdstatus, creatinine_units = "micromol/l",
  # which calculated parameter to return from the function
  return_output_type = "KDPI",
  # which year for coefficients and KDPI mapping values to use - by default the latest available in the tables ktx.kdpi_mapping_table and ktx.kdpi_coefficients_table
  mapping_values_year = "latest",
  # custom labels for factor parameters and their unknown values
    # introduce variables' labels which any user can adapt to their own labeling in the data file
    # all labels could be a character, numeric, vector, or list - whatever you prefer
    # you have to change the labels according your data file
    #  for example:
    #    if donors with a history of hypertension in your data file defined as "hypertensive", you have to change below the variable label_hypertension_positive = "hypertensive";
    #    if you use 1 for donors with a history of hypertension, you have to change below the variable label_hypertension_positive = 1
    # be careful with labeling of missing/unknown data for a history of hypertension, diabetes, HCV
    #    if donors with an unknown history of hypertension in your data file has a missing value, set the variable label_hypertension_unknown = "NA" (with quotes!)
    #    if donors with an unknown history of hypertension in your data file defined as "no data", set the variable label_hypertension_unknown = "no data"
    #    if donors with an unknown history of hypertension in your data file defined as 999, set the variable label_hypertension_unknown = 999
    #
    # label for Afroamerican ethnicity
    label_afroamerican = c ("Afroamerican"),
    # label for a positive history of hypertension
    label_hypertension_positive = c ("yes"),
    # label for an unknown history of hypertension
    label_hypertension_unknown = "NA", # if missing values defined unknown history then use "NA" (with quotes!)
    # label for a positive history of diabetes
    label_diabetes_positive = c ("yes"),
    # label for an unknown history of diabetes
    label_diabetes_unknown = "NA", # if missing values defined unknown history then use "NA" (with quotes!)
    # label for a cause of death due to cerebrovascular/stroke
    label_causeofdeath = c ("cva"),
    # label for a positive hcv status
    label_hcv_positive = c ("positive"),
    # label for an unknown, not done, indeterminate, or pending hcv status
    label_hcv_unknown = "NA", # if missing values defined unknown history then use "NA" (with quotes!)
    # label for a donation after circulatory death status
    label_dcdstatus = c ("yes")
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  # at the beginning assume that all function arguments are present, and in the following code change it to FALSE if any of the obligatory argument(s) is(are) absent
  fx_params_resulting <- TRUE
  err_num <- 0
  # Check whether all necessary params are defined by user
  fx_params <- c("age", "ethnicity", "hypertension", "diabetes", "causeofdeath", "creatinine", "hcv", "dcdstatus") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  for(i in 1:length(fx_params)){
    fx_param_local <- fx_params[i] %in% args # whether param is found among the obligatory 
    if( fx_param_local == FALSE ){
      warning("Obligatory argument ", fx_params[i], " is not defined by user in the function arguments", "\n")
      err_num <- err_num + 1
    }
    fx_params_resulting <- fx_params_resulting && fx_param_local
  }
  
  # Check of height and weight is more complex because they could be expressed in different variables by different dimentions
  if ("height_cm" %in% args){
    # do nothing
  }else{
    # if height_cm not presented check whether height_ft and height_inch are presented
    if( ("height_ft" %in% args) && ("height_inch" %in% args)){
      # do nothing
    }else{
      # if neither height_cm nor height_ft and height_inch are presented - there is no data for height at all
      warning("Obligatory argument for height is not defined by user in the function arguments", "\n")
      err_num <- err_num + 1
      fx_params_resulting <- FALSE
    }
  }
  # do the same check for weight
  if ("weight_kg" %in% args){
    # do nothing
  }else{
    # if weight_kg not presented check whether weight_lb are presented
    if( "weight_lb" %in% args ){
      # do nothing
    }else{
      # if neither weight_kg nor weight_lb are presented - there is no data for weight at all
      warning("Obligatory argument for weight is not defined by user in the function arguments", "\n")
      err_num <- err_num + 1
      fx_params_resulting <- FALSE
    }
  }
  # final message if any of the Obligatory arguments are missing
  if(fx_params_resulting == FALSE){
    stop("Obligatory argument", service.singular_or_plural(err_num, "", "s"), service.singular_or_plural(err_num, " is ", " are "), " not defined by user", "\n", "The execution of the function is interrupted.", "\n")
  }

  
  # Check the type of arguments inputed by user

  # check that user defined a single output param
  if( length(return_output_type) != 1 ){
    stop("The value for 'return_output_type' has to be a single character string. ", "The execution of the function is interrupted.", "\n")
  }
  # check that user defined a single creatinine_units
  if( length(creatinine_units) != 1 ){
    stop("The value for 'creatinine_units' has to be a single character string. ", "The execution of the function is interrupted.", "\n")
  }

  # check the range of output params
  possible_params = c("KDPI", "KDRI_Rao", "KDRI_median")
  if (!service.is.param_possible(return_output_type, possible_params)){
    stop("The defined by user value '", return_output_type, "' for 'return_output_type' is not among possible values of the parameter. ", "The execution of the function is interrupted.", "\n")
  }
  # check the range of creatinine_units
  possible_params = c("mg/dl", "micromol/l", "mmol/l")
  if (!service.is.param_possible(tolower(creatinine_units), possible_params)){
    stop("The defined by user value '", creatinine_units, "'for 'creatinine_units' is not among possible values of the parameter. ", "The execution of the function is interrupted.", "\n")
  }  
  
  # Check whether numericl arguments inputed by user are fine (weight, height etc have to be positive numbers)
  #   check and inform user whether argument contains non-numeric values
  numeric_resulting = TRUE # assume that all arguments are numeric
  numeric_local <- service.is_numeric(age)
  if(!numeric_local) warning("Age is non-numeric argument.", "\n")
  numeric_resulting <- numeric_local && numeric_resulting # if any of the argument is non-numeric, the numeric_resulting become FALSE
  if ("height_cm" %in% args){
    numeric_local <- service.is_numeric(height_cm)
    if(!numeric_local) warning("Height in cm is non-numeric argument.", "\n")
    numeric_resulting <- numeric_local && numeric_resulting
  }
  if ("height_ft" %in% args){
    numeric_local <- service.is_numeric(height_ft)
    if(!numeric_local) warning("Height in feet is non-numeric argument.", "\n")
    numeric_resulting <- numeric_local && numeric_resulting
  }
  if ("height_inch" %in% args){
    numeric_local <- service.is_numeric(height_inch)
    if(!numeric_local) warning("Height in inches is non-numeric argument.", "\n")
    numeric_resulting <- numeric_local && numeric_resulting
  }
  if ("weight_kg" %in% args){
    numeric_local <- service.is_numeric(weight_kg)
    if(!numeric_local) warning("Weight in kg is non-numeric argument.", "\n")
    numeric_resulting <- numeric_local && numeric_resulting
  }
  if ("weight_lb" %in% args){
    numeric_local <- service.is_numeric(weight_lb)
    if(!numeric_local) warning("Weight in lb is non-numeric argument.", "\n")
    numeric_resulting <- numeric_local && numeric_resulting
  }
  numeric_local <- service.is_numeric(creatinine)
  if(!numeric_local) warning("Creatinine is non-numeric argument.", "\n")
  numeric_resulting <- numeric_local && numeric_resulting # if any of the argument is non-numeric, the numeric_resulting become FALSE
  # resulting message for numeric check
  if(!numeric_resulting){
    stop("At least one of the defined by user arguments is not numeric.", "The execution of the function is interrupted.", "\n")
  }
  
  
  # check plausible biologic boundaries:
  #   check and inform user whether any values out of boundaries were substituted by NA
  #   after the check chnge the value to boundaries in the possible range (i.e. age > 0 and < 100)
  # age
  # for age calculate suspiciosly_low twice: first time for warning if <0, seond time for cat if <17
  suspiciosly_low <- service.count_lowerequal_threshold(age, 0)
  if(suspiciosly_low > 0) warning(service.output_message(suspiciosly_low, "negative values for age", "NA"))
  suspiciosly_high <- service.count_greater_threshold(age, 100)
  if(suspiciosly_high > 0) warning(service.output_message(suspiciosly_high, "age >100 years", "NA"))
  age <- service.strict_to_numeric_threshold_lower(age, 0)
  age <- service.strict_to_numeric_threshold_greater(age, 100)
  suspiciosly_low <- service.count_lowerequal_threshold(age, 17)
  if(suspiciosly_low > 0) cat(service.output_message(suspiciosly_low, "age <=17 years", "as is"))
  # height_cm
  if ("height_cm" %in% args){
    suspiciosly_low <- service.count_lowerequal_threshold(height_cm, 100)
    suspiciosly_high <- service.count_greater_threshold(height_cm, 230)
    if(suspiciosly_low > 0) warning(service.output_message(suspiciosly_low, "height <= 100 cm", "NA"))
    if(suspiciosly_high > 0) warning(service.output_message(suspiciosly_high, "height > 230 cm", "NA"))
    height_cm <- service.strict_to_numeric_threshold_lower(height_cm, 100)
    height_cm <- service.strict_to_numeric_threshold_greater(height_cm, 230)
  }
  # height_ft
  if ("height_ft" %in% args){
    suspiciosly_low <- service.count_lowerequal_threshold(height_ft, 3)
    suspiciosly_high <- service.count_greater_threshold(height_ft, 8)
    if(suspiciosly_low > 0) warning(service.output_message(suspiciosly_low, "height <= 3 feet", "NA"))
    if(suspiciosly_high > 0) warning(service.output_message(suspiciosly_high, "height > 8 feet", "NA"))
    height_ft <- service.strict_to_numeric_threshold_lower(height_ft, 3)
    height_ft <- service.strict_to_numeric_threshold_greater(height_ft, 8)
  }
  # height_inch
  if ("height_inch" %in% args){
    suspiciosly_low <- service.count_lowerequal_threshold(height_inch, -0.001)
    if(suspiciosly_low > 0) warning(service.output_message(suspiciosly_low, "negative values of inches", "NA"))
    height_inch <- service.strict_to_numeric_threshold_lower(height_inch, -0.001) # since 0 is possible value for inches
  }
  # weight_kg
  if ("weight_kg" %in% args){
    suspiciosly_low <- service.count_lowerequal_threshold(weight_kg, 30)
    suspiciosly_high <- service.count_greater_threshold(weight_kg, 300)
    if(suspiciosly_low > 0) warning(service.output_message(suspiciosly_low, "weight <= 30 kg", "NA"))
    if(suspiciosly_high > 0) warning(service.output_message(suspiciosly_high, "weight > 300 kg", "NA"))
    weight_kg <- service.strict_to_numeric_threshold_lower(weight_kg, 30)
    weight_kg <- service.strict_to_numeric_threshold_greater(weight_kg, 300)
  }
  # weight_lb
  if ("weight_lb" %in% args){
    suspiciosly_low <- service.count_lowerequal_threshold(weight_lb, 65)
    suspiciosly_high <- service.count_greater_threshold(weight_lb, 650)
    if(suspiciosly_low > 0) warning(service.output_message(suspiciosly_low, "weight <= 65 lb", "NA"))
    if(suspiciosly_high > 0) warning(service.output_message(suspiciosly_high, "weight > 650 lb", "NA"))
    weight_kg <- service.strict_to_numeric_threshold_lower(weight_lb, 65)
    weight_kg <- service.strict_to_numeric_threshold_greater(weight_lb, 650)
  }
  # creatinine
  suspiciosly_low <- service.count_lowerequal_threshold(creatinine, 0)
  if(suspiciosly_low > 0) warning(service.output_message(suspiciosly_low, "creatinine <=0", "NA"))
  creatinine <- service.strict_to_numeric_threshold_lower(creatinine, 0)

  # CHECK FUNCTION INPUT: END
  ##################################################################


  
  
  # The following variables (listed below under CUSTOMIZE_YEAR heading) are defined according OPTN reference year:
  # 1. KDRI scaling factor for appropriate year
  # 2. Chances of hypertension and diabetes in case they were unknown for donor
  #    OPTN currently assumes the chances of HCV in patient with unknown HCV status as negative, and I entered the parameter assumed_chance_of_hcv_for_unknown which is set to 0 (i.e. negative). However, in case of future modification of chances of HCV in a patient with unknown HCV status, it could be set to a probability.
  # 3. KDRI to KDPI Mapping Table.
  
  ##################################################################
  # CUSTOMIZE_YEAR: BEGIN

  # year for which OPTN coefficients are reported
  # set calculator_year to NA
  calculator_year <- NA
  # if mapping_values_year = "latest", lets define in the table ktx.kdpi_mapping_table which year is the latest
  if (mapping_values_year == "latest"){
    # take max yer from internal dataframe ktx.kdpi_mapping_table
    maxyear4ktx.kdpi_mapping_table <- max(ktx.kdpi_mapping_table$mapping_year)
    # take max yer from internal dataframe ktx.kdpi_coefficients_table
    maxyear4ktx.kdpi_coefficients_table <- max(ktx.kdpi_coefficients_table$coefficients_year)
    # now compare if in both service tables the max year is the same
    if (maxyear4ktx.kdpi_mapping_table == maxyear4ktx.kdpi_coefficients_table){
      # in both service tables the max year is the same
      # set calculator_year which will be used for taking data from ktx.kdpi_mapping_table and ktx.kdpi_coefficients_table
      calculator_year <- maxyear4ktx.kdpi_mapping_table
    }else{
      # max year is different between service tables
      cat("WARNING!", "\n")
      cat("The latest year in table ktx.kdpi_mapping_table is not matching ktx.kdpi_coefficients_table!", "\n")
      cat("Please contact the R package developer, since these tables are internal and could not be edited by the user.", "\n")
      stop("The execution of the function is interrupted.")
    }
  }else{
  # in case a user defined its own year
    # check whether the user defined year exists in the service tables ktx.kdpi_mapping_table and ktx.kdpi_coefficients_table
    useryear2ktx.kdpi_mapping_table <- mapping_values_year %in% ktx.kdpi_mapping_table$mapping_year
    if (useryear2ktx.kdpi_mapping_table != TRUE){
      cat("WARNING!", "\n")
      cat("The year ", mapping_values_year, " you defined is not presented in table ktx.kdpi_mapping_table!", "\n")
      cat("To look at the available years for OPTN mapping tables and coefficients run the following: ktx.kdpi.optn.show.years()", "\n")
      stop("The execution of the function is interrupted.")
    }

    useryear2ktx.kdpi_coefficients_table <- mapping_values_year %in% ktx.kdpi_coefficients_table$coefficients_year
    if (useryear2ktx.kdpi_coefficients_table != TRUE){
      cat("WARNING!", "\n")
      cat("The year ", mapping_values_year, " you defined is not presented in table ktx.kdpi_coefficients_table!", "\n")
      cat("To look at the available years for OPTN mapping tables and coefficients run the following: ktx.kdpi.optn.show.years()", "\n")
      stop("The execution of the function is interrupted.")
    }
    
    if (useryear2ktx.kdpi_mapping_table == TRUE & useryear2ktx.kdpi_coefficients_table == TRUE){
    # only if the user defined year is present in both tables, set the calculator_year
      calculator_year <- mapping_values_year
    }

  }

  # check whether the calculator_year is set
  if (!is.na(calculator_year)){
  # OK, the calculator_year is set, lets proceed
    # define KDRI scaling factor for the appropriate year
    kdri_scaling_factor <- ktx.kdpi_coefficients_table$kdri_scaling_factor [ktx.kdpi_coefficients_table$coefficients_year == calculator_year]
    # define chances of hypertension, diabetes, HCV in case if they were unknown for donor
    assumed_chance_of_hypertension_for_unknown <- ktx.kdpi_coefficients_table$assumed_chance_of_hypertension_for_unknown [ktx.kdpi_coefficients_table$coefficients_year == calculator_year]
    assumed_chance_of_diabetes_for_unknown <- ktx.kdpi_coefficients_table$assumed_chance_of_diabetes_for_unknown [ktx.kdpi_coefficients_table$coefficients_year == calculator_year]
    assumed_chance_of_hcv_for_unknown <- ktx.kdpi_coefficients_table$assumed_chance_of_hcv_for_unknown [ktx.kdpi_coefficients_table$coefficients_year == calculator_year]
    # define mapping values for calculation of KDPI
    kdri_lower <- ktx.kdpi_mapping_table$mapping_kdri_lower [ktx.kdpi_mapping_table$mapping_year == calculator_year]
    kdri_upper <- ktx.kdpi_mapping_table$mapping_kdri_upper [ktx.kdpi_mapping_table$mapping_year == calculator_year]
    kdpi_percent_all <- ktx.kdpi_mapping_table$mapping_kdpi_percent [ktx.kdpi_mapping_table$mapping_year == calculator_year]
    cat("For the calculations of KDRI and KDPI the ", calculator_year, " year is used for mapping values, KDRI scaling factor, as well as chances of hypertension and diabetes in case if they were unknown for donor. ", "\n")
  }else{
  # the calculator_year is not set, stop everything
    cat("WARNING!", "\n")
    cat("The variable calculator_year is not set!", "\n")
    stop("The execution of the function is interrupted.")
  }
  
  # CUSTOMIZE_YEAR: END
  ##################################################################

  #
  # repeat creatinine_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)
  creatinine_units <- rep(creatinine_units, length(creatinine))  
  
  #
  # convert creatinine units if necessary
  #
  creatinine <- ifelse (  creatinine_units == "mg/dl",
              creatinine, # do nothing 
              ifelse(  creatinine_units == "micromol/l",
                  creatinine / 88.4,
                  ifelse(  creatinine_units == "mmol/l",
                    1000 * creatinine / 88.4,
                    NA # if any other undefined units is used, assume NA
                  )
              )
            )

  #
  # limit maximum creatinine value
  #
  creatinine <- ifelse (creatinine > 8, 8, creatinine)
  

  #
  # convert height units if necessary
  #
  height <-   ifelse(  height_cm > 0,
            height_cm, # if height is defined in cm, just take this value
            ifelse(height_ft > 0 | height_inch > 0, # if height is defined NOT in cm, check whether height in feets or inches is defined
            height_ft * 30.48 + height_inch * 2.54, # if they contain number than convert to cm
            NA) # if height in feets or inches is not defined, assume NA
        )

  #
  # convert weight units if necessary
  #
  weight <-   ifelse(  weight_kg > 0, 
            weight_kg, # if weight is defined in kg, just take this value
            ifelse  (  weight_lb > 0, # if weight is defined NOT in kg, check whether weight in pounds is defined
                  weight_lb * 0.453592, # if it contains number than convert to kg
                  NA) # if weight in pounds is not defined, assume NA
            )


  # 
  # assign xbeta coefficients to all factors
  # 

  # assign xbeta to age
  age_xbeta <- ifelse(  age < 18,
              0.0128 * ( age - 40 ) - 0.0194 * ( age - 18 ),
              ifelse( age > 50,
                  0.0128 * ( age - 40 ) + 0.0107 * ( age - 50 ),
                  0.0128 * ( age - 40 )
                  )
            )
  

  # assign xbeta to height
  height_xbeta <- -0.0464 * ( height - 170 ) / 10

  # assign xbeta to weight
  weight_xbeta <- ifelse(  weight < 80,
              -0.0199 * ( weight - 80 ) / 5,
              0
              )

  # assign xbeta to ethnicity
  ethnicity_xbeta <- ifelse(  ethnicity %in% label_afroamerican,
                0.1790,
                0
              )
                

  # assign xbeta to history of hypertension
  hypertension_xbeta <- ifelse(    hypertension %in% label_hypertension_positive,
                    0.1260,
                    ifelse(  hypertension == label_hypertension_unknown | (label_hypertension_unknown == "NA" & is.na(hypertension)),
                        0.1260 * assumed_chance_of_hypertension_for_unknown,
                        0
                        )
                )

  # assign xbeta to history of diabetes
  diabetes_xbeta <- ifelse(  diabetes %in% label_diabetes_positive,
                0.1300,
                ifelse(  diabetes == label_diabetes_unknown | (label_diabetes_unknown == "NA" & is.na(diabetes)),
                    0.1300 * assumed_chance_of_diabetes_for_unknown,
                    0
                    )
              )

  # assign xbeta to cause of death
  causeofdeath_xbeta <- ifelse(  causeofdeath %in% label_causeofdeath,
                  0.0881, 
                  0
                )

  # assign xbeta to serum creatinine
  creatinine_xbeta <- ifelse(  creatinine > 1.5,
                0.2200 * ( creatinine - 1 ) - 0.2090 * ( creatinine - 1.5 ),
                0.2200 * ( creatinine - 1 ) 
                )
  

  # assign xbeta to HCV status
  hcv_xbeta <- ifelse(  hcv %in% label_hcv_positive,
              0.2400,
              ifelse(  hcv == label_hcv_unknown | (label_hcv_unknown == "NA" & is.na(hcv)),
                  0.2400 * assumed_chance_of_hcv_for_unknown,
                  0
                  )
            )

  # assign xbeta to donation after circulatory death status
  dcd_xbeta <- ifelse(dcdstatus %in% label_dcdstatus,
            0.1330,
            0
            )

  #
  # calculate total xbeta 
  #
  total_xbeta <- age_xbeta + height_xbeta + weight_xbeta + ethnicity_xbeta + hypertension_xbeta + diabetes_xbeta + causeofdeath_xbeta + creatinine_xbeta + hcv_xbeta + dcd_xbeta

  # KDRI Rao is interpreted as the relative risk of post-transplant graft failure for this donor compared to a reference donor (age=40 years, non-African American, etc.) as defined in Rao, et al.
  kdri_rao <- exp (total_xbeta)

  # KDRI MEDIAN is interpreted as the relative risk of post-transplant graft failure (in an average, adult recipient) for this donor compared to the median kidney donor recovered last year. In the descriptive text in DonorNet explaining the KDPI, for example, "The estimated risk of kidney graft failure from this donor is higher than 74% of all kidney donors recovered in 2010 and 1.28 times that of the median donor from 2010", the value 1.28 is the "scaled to the median" KDRI.
  kdri_median <- kdri_rao / kdri_scaling_factor

  # KDPI mapping to the reference OPTN population  
  kdpi_mapping <- rep(NA, length(kdri_median))
  for (i in 1:length(kdri_median)){
    kdpi_mapping[i] <- ifelse  ( is.na(kdri_median[i]), # check whether the kdri_median[i] is NULL, because if so the which will produce "replacement has length zero" error
                  NA, 
                  which( abs (kdri_lower - kdri_median[i]) == min (abs (kdri_upper - kdri_median[i]) ) ) # if not NA - calculate
                  )
  }
#    kdpi_mapping[i] <- which( abs (kdri_lower - kdri_median[i]) == min (abs (kdri_upper - kdri_median[i]) ) )

  # if not all necessary input parameters are present, kdri_rao will be NA, but kdpi_percent will be mapped to 0. Due to this introduce additional check to output kdpi_percent as NA in case not all paramenters are present
  kdpi_percent <- ifelse (is.na(kdri_rao),
              NA,
              kdpi_percent_all[kdpi_mapping]
              )

  # round kdri_rao after all calculations just for the nice output              
  kdri_rao <- round(kdri_rao, 2) # initially I used 3 digits, but the revewers demanded to output all with 2 digits 

  # round kdri_median to two decimals before return
  kdri_median <- round(kdri_median, 2)

  # if the year for which OPTN coefficients are reported is more than 2 years far from the current year, show warning to the user 
  current_year <- format(Sys.Date(), "%Y")
  current_year <- as.integer(current_year)
  if (current_year - calculator_year > 2){
    cat("WARNING!", "\n")
    cat("The year for which OPTN coefficients are reported is more than 2 years far from the current year", "\n")
    cat("Please use updated OPTN coefficients.", "\n")
  }

  # which data to return as the function output
  return_output_type <- rep(return_output_type, length(creatinine))  
  function_output <- ifelse (  return_output_type == "KDPI",
                kdpi_percent,
                ifelse(  return_output_type == "KDRI_Rao",
                    kdri_rao,
                    ifelse(  return_output_type == "KDRI_median",
                        kdri_median,
                        NA # just in case neither KDPI no KDRI is requested by user
                        )
                    )
                )


return (function_output)

}


# FUNCTION: END
##################################################################



##################################################################
# FUNCTION: BEGIN

#' Shows which years are available in the R package for the OPTN mapping table, KDRI scaling factor, etc. 
#' @details Service function which shows for user for which year(s) the OPTN mapping table, as well as KDRI scaling factor and chances of hypertension and diabetes in case if they were unknown for donor in the ktx.kdpi_mapping_table and ktx.kdpi_coefficients_table. This years could be used for the argument \emph{mapping_values_year} of the ktx.kdpi.optn function.
#' 
#' This function has no arguments.
#' @return numeric List of years which could be used for the argument mapping_values_year of the ktx.kdpi.optn function.
# Because list of years has to be the same in the ktx.kdpi_mapping_table and ktx.kdpi_coefficients_table (and this assumption is checked in the ktx.kdpi.optn function), I will look only into the ktx.kdpi_coefficients_table
#' @export
#' @name ktx.kdpi.optn.show.years
ktx.kdpi.optn.show.years <- function() {
  
  list_of_years <- unique(ktx.kdpi_coefficients_table$coefficients_year)
  cat("The following years are presented in tables with OPTN mapping values, KDRI scaling factor, etc.", "\n")
  print (list_of_years)
  cat("You can use any of this years in the mapping_values_year parameter of the ktx.kdpi.optn function for mapping calculations of KDPI and KDRI.", "\n")
}
  
# FUNCTION: END
##################################################################

