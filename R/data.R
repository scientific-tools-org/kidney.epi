#' A sample dataset with kidney transplant patients.
#'
#' A dataset contains 10 records for kidney transplant patients, including information about deceased donors.
#'
#' @format A data frame with 10 rows and 12 variables:
#' \describe{
#'   \item{ptid}{patient identifier}
#'   \item{rec.age}{age of the recipient, in years}
#'   \item{don.age}{age of the donor, in years}
#'   \item{don.height}{height of the donor, in cm}
#'   \item{don.weight}{weight of the donor, in kg}
#'   \item{don.ethnicity}{ethnicity of the donor}
#'   \item{don.hypertension}{history of hypertension for the donor}
#'   \item{don.diabetes}{history of diabetes for the donor}
#'   \item{don.causeofdeath}{cause of death for the donor}
#'   \item{don.creatinine}{serum creatinine of the donor, in mg/dL}
#'   \item{don.hcv}{hepatitis c virus status of the donor}
#'   \item{don.dcdstatus}{donation after circulatory death status of the donor}
#'   \item{don.sex}{sex of the donor}
#' }
#' @source Generation from different patients' records
"ktx.data"



#' A synthetic dataset contains variables for eGFR calculation.
#'
#' A synthetic dataset contains variables for eGFR calculation for 1000 adults and 1000 children.
#'
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#'
#' @format A data frame with 2000 rows (1000 adults and 1000 children/young adults) and 12 variables:
#' \describe{
#'   \item{cr}{Serum creatinine, micromol/L}
#'   \item{cys}{Serum cystatin C, mg/L}
#'   \item{age}{Age, years}
#'   \item{sex}{Sex}
#'   \item{ethnicity}{Ethnicity}
#'   \item{height}{Height, cm}
#'   \item{category}{Indication on whether the generated data refer to adults or children}
#' }
#' @source Synthetic dataset was generated based on two publications:  
#' \itemize{
#'   \item adults: Lamb EJ, Barratt J, Brettell EA et al. Accuracy of glomerular filtration rate estimation using creatinine and cystatin C for identifying and monitoring moderate chronic kidney disease: the eGFR-C study.  Health Technol Assess 2024;28(35), doi:10.3310/HYHN1078.  
#'   \item children/young adults: Pierce CB, Muñoz A, Ng DK et al. Age- and sex-dependent clinical equations to estimate glomerular filtration rates in children and young adults with chronic kidney disease. Kidney International. 2021;99(4):948–956, doi:10.1016/j.kint.2020.10.047.
#' }
"ckd.data"