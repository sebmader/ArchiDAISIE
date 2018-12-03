### pipeline to analyse output of ArchiDAISIE ###
### written by Sebastian Mader on the 20.11.2018, s.mader@student.rug.nl ###
library(DAISIE)
library(readr)

col_names = c("Clade_name", "Status", "Missing_species", "Branching_times")
rep_1_branching <- read_csv("~/Studium/Master/RU Groningen/Research Project 1/DAISIE_Archipelago/ArchiDAISIE/ArchiDAISIE/test_sims/sim_1/rep_1_branching.txt", 
                            col_names = col_names)
View(rep_1_branching)
prepData <- DAISIE_dataprep(rep_1_branching,10,100)
#estimates <- DAISIE_ML(prepData)
  #which starting points for ML? -> the once I gave it?
