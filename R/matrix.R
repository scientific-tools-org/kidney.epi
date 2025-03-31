#' Access matrix cells by row and column names
##################################################################
# FUNCTION: BEGIN
#' @details Access matrix cells by row and column names.
#' @param matrix_data matrix. Matrix with data.
#' @param row_name character. Row name in the matrix.
#' @param col_name character. Column name in the matrix.
#' @return vector. Matrix value.
#' @export
#' @name matrix.get_named_matrix_value
matrix.get_named_matrix_value <- function(matrix_data, row_name, col_name) {
  if (row_name %in% rownames(matrix_data) && col_name %in% colnames(matrix_data)) {
    return(matrix_data[row_name, col_name])
  } else {
    stop("Specified row or column name not found in matrix")
  }
}

# FUNCTION: END
##################################################################


#' Read Excel file and convert it to matrix with row and column names.
##################################################################
# FUNCTION: BEGIN
#' @import readxl
#' @import purrr
#' @details Read Excel file and convert it to matrix with row and column names.
#' @param file_path character. Path to Excel file.
#' @param sheet_name character. Name of Excel sheet. Optional, if there is only one sheet with data, the function will read it with no need to specifying the sheet name.
#' @return Excel file saved to a specified folder.
#' @export
#' @name matrix.read_excel_to_named_matrix
matrix.read_excel_to_named_matrix <- function(file_path, sheet_name = NULL) {
  # Get list of all sheet names in the file
  sheet_names <- excel_sheets(file_path)
  
  # Check if sheet_name is provided
  if (!is.null(sheet_name)) {
    # Check if the specified sheet exists
    if (sheet_name %in% sheet_names) {
      data <- readxl::read_excel(file_path, sheet = sheet_name, .name_repair = "minimal")
    } else {
      message("The specified sheet '", sheet_name, "' does not exist in the file.")
      return(NULL)
    }
  } else {
    # sheet_name is NULL, automatically decide which sheet to read
    # Read data from each sheet and check if it has content
    non_empty_sheets <- purrr::keep(sheet_names, function(sheet) {
      data <- readxl::read_excel(file_path, sheet = sheet, .name_repair = "minimal")
      any(!is.na(data))  # Return TRUE if the sheet has non-NA data
    })
    
    if (length(non_empty_sheets) == 1) {
      # If only one non-empty sheet, read that sheet
      data <- readxl::read_excel(file_path, sheet = non_empty_sheets[[1]], .name_repair = "minimal")
    } else if (length(non_empty_sheets) > 1) {
      # If multiple non-empty sheets, print a message and do not load the data
      message("Multiple sheets with data found. Please specify a sheet using the 'sheet_name' argument.")
      return(NULL)
    } else {
      # No non-empty sheets found
      message("No data found in any sheet.")
      return(NULL)
    }
  }
  
  # Set the first column as row names and remove it from the data
  row_names <- data[[1]]
  data <- data[,-1]
  data_matrix <- as.matrix(data)
  rownames(data_matrix) <- row_names
  
  return(data_matrix)
}

# FUNCTION: END
##################################################################




#' Save a named matrix as an Excel file.
##################################################################
# FUNCTION: BEGIN
#' @import openxlsx
#' @details Save a named matrix as an Excel file.
#' @param matrix_data matrix. Matrix for saving.
#' @param file_path character. Path to the Excel file.
#' @param sheet_name character. Name of the Excel sheet.
#' @param save_type character. Defines whether the Excel file should be created or overwritten (save_type = "new"), or new sheet should be added to the existing Excel file (save_type = "add").
#' @return Excel file saved to a specified folder.
#' @export
#' @name matrix.save_named_matrix_to_excel
matrix.save_named_matrix_to_excel <- function(matrix_data, file_path, sheet_name = "Sheet1", save_type = "new") {
  # Convert the matrix to a data frame, adding row names as the first column
  df <- as.data.frame(matrix_data)
  df <- cbind(RowNames = rownames(matrix_data), df)
  
  if (save_type == "new") {
    # Create a new workbook and add the data as a new sheet
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, df)
    
    # Save workbook as a new file
    openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
    message("Matrix saved successfully to a new file at ", file_path)
    
  } else if (save_type == "add") {
    # Check if the file exists
    if (file.exists(file_path)) {
      # Load existing workbook
      wb <- openxlsx::loadWorkbook(file_path)
    } else {
      # If file doesn't exist, create a new workbook
      wb <- openxlsx::createWorkbook()
    }
    
    # Add the new sheet with data
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, df)
    
    # Save workbook, updating the file
    openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
    message("Matrix added to existing file at ", file_path)
    
  } else {
    stop("Invalid save_type. Choose either 'new' or 'add'.")
  }
}


# FUNCTION: END
##################################################################




#' Creates a named matrix from two variables.
##################################################################
# FUNCTION: BEGIN
#' @details Creates a named matrix from two variables.
#' @param var1 Character vector. Values representing first variable.
#' @param var2 Character vector. Values representing second variable.
#' @param predefined_levels Character vector. Levels for var1 and var2. If omitted, the variables just coded according to the levels they have. If contains vector, the variables are coded according to predefined_levels values, the latter could be useful if var1 and var2 contain not all levels of interest
#' @return matrix with cross-tabulation of var1 and var2.
#' @export
#' @name matrix.cross_table
matrix.cross_table <- function(var1, var2, predefined_levels = NA) {

  var1_name <- deparse(substitute(var1))
  var2_name <- deparse(substitute(var2))

  if (length(var1) != length(var2)) {
    warning("var1 and var2 do not have the same length.")
  }
  
  # Convert variables to factors with predefined levels if provided
  if (!is.na(predefined_levels[1])) {
    missing_levels_var1 <- setdiff(unique(var1), predefined_levels)
    missing_levels_var2 <- setdiff(unique(var2), predefined_levels)
    
    if (length(missing_levels_var1) > 0 || length(missing_levels_var2) > 0) {
      warning("Some values in var1 or var2 are not in predefined_levels: ", 
              paste(unique(c(missing_levels_var1, missing_levels_var2)), collapse = ", "))
    }
    
    var1 <- factor(var1, levels = predefined_levels)
    var2 <- factor(var2, levels = predefined_levels)
  } else {
    var1 <- factor(var1)
    var2 <- factor(var2)
  }
  
  # Create a cross-table
  cross_tab <- table(var1, var2)
  
  # Convert cross-table to matrix
  cross_matrix <- unclass(as.matrix(cross_tab))

  # Print variable names in console
  cat(var1_name, "represents rows, ", var2_name, "represents columns\n")
    
  return(cross_matrix)
}

# FUNCTION: END
##################################################################
