CE <- function(data,N_max=10,eps=0.01,rho=0.05,M=200,h=5, a=0.8,parallel=FALSE){
  
  if(is.data.frame(data)=="FALSE"| is.null(dim(data)[2])) {print("Error in data : dataframe only")
                                                           
  } else if(dim(data)[2]!=1) {print("Error in data : single column dataframe only") 
                              
  } else { 
    create.folder("CE")
    Nelite<-M*rho
    L<-length(data[,1])
    L0<-1
    
    k<-seq(0,N_max,1)
    
    if(parallel==TRUE & .Platform$OS.type == "windows"){
      
      suppressPackageStartupMessages(pkgs <- require(doSNOW))           
      if (!pkgs) {
        stop("doSNOW and snow packages are not available for parallel computation!")
      } else {
        cl<-makeCluster(parallel::detectCores(), type="SOCK") 
        clusterExport(cl, c("ce.beta","betarand", "fun_alpha","fun_beta","modi.BICnorm.new"), envir=environment())
        clusterExport(cl, c("data","rho","M","h","eps","Nelite","L","L0","a"), envir=environment())  
        registerDoSNOW(cl)
        
        write.table(system.time(sim<-foreach(k=k, .errorhandling=c('pass')) %dopar% ce.beta(k,data,h,L0,L,M,Nelite,eps,a))[[3]],file="Processing Time.txt",sep="\t")
        stopCluster(cl)
        
      }                            
      
    }  else if (parallel==TRUE & .Platform$OS.type == "unix"){
      
      suppressPackageStartupMessages(pkgs <- require(doMC))        
      if (!pkgs) {
        stop("doMC and parallel packages are not available for parallel computation!")
      } else {
        registerDoMC(parallel::detectCores()) 
        write.table(system.time(sim<-foreach(k=k, .errorhandling=c('pass')) %dopar% ce.beta(k,data,h,L0,L,M,Nelite,eps,a))[[3]],file="Processing Time.txt",sep="\t")
      }                            
      
    } else write.table(system.time(sim<-foreach(k=k, .errorhandling=c('pass')) %do% ce.beta(k,data,h,L0,L,M,Nelite,eps,a))[[3]],file="Processing Time.txt",sep="\t")
    
    k<-seq(0,N_max,1)
    
    mBic_summary<-c()
    for(i in 1:length(sim)){
      mBic_summary[i]<-sim[[i]]$mBIC
    }
    
    loci.mBIC<-sim[[k[which(mBic_summary==max(mBic_summary))]+1]]$locis
    
    xaxis<-seq(1,L,1)
    dd.mbic<-plott(data,loci.mBIC,L)[,1]
    
    plt<-qplot(xaxis,data[1:L,1],data=as.data.frame(xaxis),xlab="Data Sequence", ylab="Value",xlim=c(1,L))  + geom_line(aes(xaxis,dd.mbic), colour="red", size=1) + theme(axis.line=element_line(colour="black"), panel.grid.minor=element_blank(), panel.grid.major=element_blank(),panel.border=element_blank())
    print(plt)     
    dev.copy(pdf, 'Mean Profile plot.pdf')
    dev.off()
    
    if(length(loci.mBIC)>=3) {
      write.table(loci.mBIC[2:(length(loci.mBIC)-1)],file=paste("Break-Points.txt"),sep="\t")
      return(list("No.BPs"=length(loci.mBIC)-2,"BP.Loc" =loci.mBIC[2:(length(loci.mBIC)-1)]))
    } else {return(paste("No Break-Points are Estimated")) }
    
  }
}