install.packages("drake")
require(drake)
setwd("C:\\Users\\chris\\Dropbox\\Postdoc\\Code")
source("BioenergeticsFunctions.R")

plan <- drake_plan(
  raw_data = read.csv("SnowCrabParams.csv")
  
)