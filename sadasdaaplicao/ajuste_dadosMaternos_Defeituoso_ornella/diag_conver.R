graphs_diag <- function(chain,par.names,archive){
  npar <- dim(chain)[2]
  #autocorrelation
  pdf(file=paste("Atcor_",archive,".pdf",sep=""))
   for (v in 1:npar) {
    acf(chain[,v],main=" ",xlab=" ")
   }
  dev.off() 
  #time series
  pdf(file=paste("TimeSer_",archive,".pdf",sep=""))
  for (v in 1:npar) {
     plot(ts(chain[,v]),main=" ",xlab=" ",ylab=" ")
  }
  dev.off() 
  #histogram
  pdf(file=paste("Hist_",archive,".pdf",sep=""))
  for (v in 1:npar) {
    hist(chain[,v],prob=TRUE,ylab="Density",xlab=" ",main=" ")
    lines(density(chain[,v]),col="red")      
  } 
  dev.off() 
}


