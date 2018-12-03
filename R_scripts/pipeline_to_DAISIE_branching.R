### pipeline to analyse branching output of ArchiDAISIE ###
### written by Sebastian Mader on the 20.11.2018, s.mader@student.rug.nl ###

library(DAISIE)
library(readr)

rep_1_branching <- read_csv("ArchiDAISIE/ArchiDAISIE/sims/sim_1/rep_1_branching.txt", 
                            col_names = TRUE)
View(rep_1_branching)
prepData <- DAISIE_dataprep(rep_1_branching,10,100)
DAISIE_plot_island(prepData,10)
#estimates <- DAISIE_ML(prepData)
  #which starting points for ML? -> the once I gave it => YES!
