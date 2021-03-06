### pipeline to plot STT output of ArchiDAISIE ###
### per parameter text file ###
### written by Sebastian Mader on the 20.11.2018, s.mader@student.rug.nl ###

rm(list = ls())

.libPaths("/data/s3503704/R/3.5")

library(DAISIE)
library(readr)


replicates <- 1000
y_lims <- c(20,50,150)

main_dir <- ""
parameter_file <- "even_more_pars.txt"
dir_name_sims <- "sims/"
parameter_sets <- read.csv(file = paste(main_dir,parameter_file,sep = ""),header = T, sep = ",")
sim_names <- parameter_sets$name
num_species <- list()
for(sim in 1:length(sim_names)) {
  num_species[[sim]] = list(sim_name = sim_names[[sim]],
                            median_endemic = vector(),
                            median_nonend = vector(),
                            median_total = vector())
}
num_species

for(y_lim in y_lims) {
  for(sim_name in sim_names) {
    island_replicates = list()
    for(rep in 1 : replicates) {
      island_replicates[[rep]] = list()
      stt_all <- as.matrix(read_csv(paste(main_dir, dir_name_sims, sim_name, "/rep_", as.character(rep), "_STT.txt", sep = "")))
      island_replicates[[rep]][[1]] = list(stt_all = stt_all)
    }
    island_replicates
    length(island_replicates)
    png(filename = paste(main_dir,"figures/ylim",y_lim,"/",sim_name,".png",sep=""), type = "cairo", units = "in",
        width = 6, height = 6, res = 300)
    
    ### part from DAISIE_plot_sims
    time<-max(island_replicates[[1]][[1]]$stt_all[,1])
    
    ###  STT ALL species
    s_freq<-length(island_replicates[[1]][[1]]$stt_all[,1])
    complete_arr <-array(dim=c(s_freq,6,replicates))
    
    for(x in 1:replicates)
    {
      sum_endemics<-island_replicates[[x]][[1]]$stt_all[,"nA"]+island_replicates[[x]][[1]]$stt_all[,"nC"]
      total<-island_replicates[[x]][[1]]$stt_all[,"nA"]+island_replicates[[x]][[1]]$stt_all[,"nC"]+island_replicates[[x]][[1]]$stt_all[,"nI"]
      complete_arr[,,x]<-cbind(island_replicates[[x]][[1]]$stt_all[,c('Time',"nI","nA","nC")],sum_endemics,total)
    }
    
    stt_average_all<-apply(complete_arr,c(1,2),median)
    stt_q0.025_all<-apply(complete_arr,c(1,2),quantile,0.025)
    stt_q0.25_all<-apply(complete_arr,c(1,2),quantile,0.25)
    stt_q0.75_all<-apply(complete_arr,c(1,2),quantile,0.75)
    stt_q0.975_all<-apply(complete_arr,c(1,2),quantile,0.975)
    
    colnames(stt_average_all)<-c("Time","nI","nA","nC","Endemic","Total")
    colnames(stt_q0.025_all)<-c("Time","nI","nA","nC","Endemic","Total")
    colnames(stt_q0.25_all)<-c("Time","nI","nA","nC","Endemic","Total")
    colnames(stt_q0.75_all)<-c("Time","nI","nA","nC","Endemic","Total")
    colnames(stt_q0.975_all)<-c("Time","nI","nA","nC","Endemic","Total")
    
    par(mfrow=c(1,1)) 
    suppressWarnings(plot(NULL,NULL,xlim=rev(c(0, time)),ylim=c(1,y_lim),ylab="No of species + 1",bty="l", xaxs="i",xlab="Time before present",main="Species-through-time - All species",log='y',cex.lab=1.2,cex.main=1.2,cex.axis=1.2))
    polygon(c(stt_average_all[,"Time"],rev(stt_average_all[,"Time"])),c(stt_q0.025_all[,"Total"]+1,rev(stt_q0.975_all[,"Total"]+1)),col="light grey",border=NA)
    polygon(c(stt_average_all[,"Time"],rev(stt_average_all[,"Time"])),c(stt_q0.25_all[,"Total"]+1,rev(stt_q0.75_all[,"Total"]+1)),col="dark grey",border=NA)
    lines(stt_average_all[,"Time"],stt_average_all[,"Total"]+1,lwd=2)
    lines(stt_average_all[,"Time"],stt_average_all[,"nI"]+1,lwd=2,col='cyan3')
    lines(stt_average_all[,"Time"],stt_average_all[,"Endemic"]+1,lwd=2,col='dodgerblue1')
    
    legend(time,max(stt_q0.975_all),c("Total","Non-endemic","Endemic"),lty=1,lwd=2,col=c("black","cyan3","dodgerblue1"),cex=1.2,border=NA,bty='n')
    
    dev.off()
  }
}