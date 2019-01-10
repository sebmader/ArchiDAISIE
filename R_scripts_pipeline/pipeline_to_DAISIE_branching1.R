### pipeline to analyse branching output of ArchiDAISIE ###
### written by Sebastian Mader on the 20.11.2018, s.mader@student.rug.nl ###

rm(list = ls())

library(DAISIE)
library(readr)

dir_name_sims <- "sims/"
sim_names <- dir(path = paste(dir_name_sims, sep=""))
parameter_sets <- read.csv(file = paste("parameter_sets.txt",sep = ""),header = T,sep = ",")
length(sim_names)/3
sim_names1 <- sim_names[1:length(sim_names)/3]

for(sim_name in sim_names1) {
  sim_name
  replicate_estimates = list()
  parameters <- as.vector(parameter_sets[which(parameter_sets[,1]==sim_name),])
  parameters
  assertthat::are_equal(sim_name, as.character(parameters[[1]]))
  archiImmi <- parameters[[2]]
  archiCladoL <- parameters[[4]]
  archiAnaL <- parameters[[5]]
  archiExtL <- parameters[[6]]
  archiCladoG <- parameters[[7]]
  archiAnaG <- parameters[[8]]
  archiExtG <- parameters[[9]]
  archiKPerIsl <- parameters[[10]]
  archiNIsl <- parameters[[11]]
  archiAge <- parameters[[12]]
  archiMainlandSp <- parameters[[13]]
  archiReplicates <- parameters[[14]]
  immi <- archiImmi * archiNIsl
  clado <- archiCladoL + archiCladoG
  ana <- archiAnaL + archiAnaG
  ext <- archiExtL + archiExtG
  K <- archiKPerIsl * archiNIsl
  initOptPars <- c(clado,ext,K,immi,ana)
  for(rep in 1:archiReplicates) {  
    branching <- read_csv(paste(dir_name_sims, sim_name, "/rep_", rep, "_branching.txt", sep = ""), 
                          col_names = TRUE)
    if(!is.na.data.frame(branching)) {
      prepData <- DAISIE_dataprep(branching,island_age = archiAge, M = archiMainlandSp)
      for(i in 2:length(prepData)) {
        prepData[[i]]$missing_species <- as.numeric(prepData[[i]]$missing_species)
      }
      replicate_estimates[rep] <- DAISIE_ML(datalist = prepData,
                                            datatype = "single",
                                            initparsopt = initOptPars,
                                            idparsopt = c(1,2,3,4,5),
                                            parsfix = c(),
                                            idparsfix = c())
    }
    else {
      replicate_estimates[rep] <- data.frame()
    }
  }
  assertthat::are_equal(length(replicate_estimates),archiReplicates)
  replicate_estimates
  save(replicate_estimates, parameters, initOptPars, 
       file = paste("estimates/",sim_name, "_estimates.Rdata", sep = ""))
}

