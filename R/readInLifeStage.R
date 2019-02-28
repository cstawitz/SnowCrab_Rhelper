#'Function to read in individuals of a certain lifestage from Results.csv file
#' @param resultsDir - string, the directory results from the model run are found
#' @param LifeStage - string, the name of the life history type
#' @param numAttributes - integer, the number of attributes specific to your life stage  
readInLifeStage <- function(resultsDir, LifeStage, numAttributes, resultsName){
  require(dplyr)
  require(readr)
  setwd(resultsDir)
  #select which columns to read
  colspec <- cols(
    .default = col_integer(),
    typeName = col_character(),
    time = col_datetime(format = ""),
    horizPos1 = col_double(),
    horizPos2 = col_double(),
    vertPos = col_double(),
    number = col_double(),
    size = col_double(),
    weight = col_double(),
    temperature = col_double(),
    age = col_double(),
    ageInInstar = col_double()
  )
  
  #Read in csv
  f <- function(x, pos) subset(x, (typeName==LifeStage))
  dfLifeStage <- read_csv_chunked(file=resultsName, DataFrameCallback$new(f),col_types = colspec)
  
  
  return(dfLifeStage)
}