library(tidyverse)


#' clean read-in of coverage file
#'
#' @param file_name: coverage file 
#'
#' @return
#' @export
#'
#' @examples
read_and_mark_cov <- function(file_name){
  
  df <- readr::read_table(file_name,
                          skip = 1,
                          col_names=c("ref", "Position", "Coverage"))
  df$barcode <- str_extract(file_name, "barcode\\d+")
  
  return(df)
}
