# Kidney.epi R package

This is the source code for the "kidney.epi" R package that contains kidney-related functions for clinical and epidemiological research.  

The kidney.epi package is made with care by the research consultancy [Scientific-Tools.Org](https://Scientific-Tools.Org).
Contact us for data analysis or software development at [Scientific-Tools.Org](https://Scientific-Tools.Org/contact) or via 'maintainer("kidney.epi")', connect with the [author on LinkedIn](https://www.linkedin.com/in/boris-bikbov).  
Support this and other Scientific-Tools.Org projects: https://Scientific-Tools.Org/support-us/  
Home page of the package: https://kidney.scientific-tools.org/r  
Connect with us on social platforms: [LinkedIn](https://www.linkedin.com/company/scientific-tools-org/) [X/Twitter](https://twitter.com/SciToolsOrg)


## Citation for journal publications
	citation("kidney.epi")
Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. doi:10.32614/CRAN.package.kidney.epi.

## Agreement on prefix of function names and internal data frames
nephro. - related to nephrology in general  
ktx. - related to kidney transplantation  
hd. - related to hemodialysis  
pd. - related to peritoneal dialysis  
egfr. - related to equations for calculation of estimated glomerular filtration rate  
epi. - related to epidemiology in general  
service. - service functions or internal purposes of the package (data check, etc)  
matrix. - functions for working with matrix  

The package gets posted to the comprehensive R archive (CRAN) at intervals, each such posting preceded a thorough test. In general, each new push to CRAN with new function(s) will update the second term of the version number, e.g. 1.2.0 to 1.3.0. Updates only to the code of existing functions increment the third term of the version number, e.g. 1.2.0 to 1.2.1.
	
## Internal datafarames
Internal datafarames are contained in the R/sysdata.rda of the source code, used in the R package functions *but not accessible to the user*.  

- ktx.kdpi_mapping_table - contains data with mapping KDPI and KDRI reported by OPTN for the years 2014-2024  
- ktx.kdpi_coefficients - contains data with coefficients used by OPTN for the calculation (KDRI scaling factor, chances of hypertension and diabetes in case if they were unknown for donor, etc)  

## External datafarames 
- ckd.data - contains synthetic data for eGFR calculation in 1000 adults and 1000 children.  
- ktx.data - contains data for 10 kidney transplant patients.
	
## Vignettes
List of vignettes are available via  

	browseVignettes(package = "kidney.epi").

## License

This project is licensed under the [CC BY-NC 4.0 License](https://creativecommons.org/licenses/by-nc/4.0/).
