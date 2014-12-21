CE.Normal <-
function(data, Nmax=10, eps=0.01, rho=0.05, M=200, h=5, a=0.8, b=0.8, distyp = 1, parallel=FALSE){
  
  if(is.data.frame(data) == "FALSE"| is.null(dim(data)[2])) {
    print("Error in data : dataframe only")                                                           
  } else if(dim(data)[2] != 1) {
    print("Error in data : single column dataframe only")                               
  } else { 
    
    if(distyp == 1){
      Melite <- M * rho
      L <- length(data[, 1])
      L0 <- 1      
      k <- seq(0, Nmax, 1)
      
        if(parallel == TRUE & .Platform$OS.type == "windows"){
            
            cl <- makeCluster(parallel::detectCores(), type="SOCK") 
            clusterExport(cl, c("ce.sim4beta", "betarand", "fun.alpha", "fun.beta", "mBIC"), envir=environment())
            clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a"), envir=environment())  
            registerDoParallel(cl)
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta(k, data, h, L0, L, M, Melite, eps, a)
            stopCluster(cl)   

            } else if (parallel == TRUE & .Platform$OS.type == "unix"){         
            
            registerDoParallel(parallel::detectCores()) 
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta(k, data, h, L0, L, M, Melite, eps, a)                           
             
        } else sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.sim4beta(k, data, h, L0, L, M, Melite, eps, a)
           
    k <- seq(0, Nmax, 1)    
    mBic.summary <- c()
      for(i in 1:length(sim)){
        mBic.summary[i] <- sim[[i]]$mBIC
      }
    
    loci.mBIC <- sim[[k[which(mBic.summary == max(mBic.summary))] + 1]]$loci
    
    if(length(loci.mBIC) >= 3) {
      return(list("No.BPs" = length(loci.mBIC) - 2, "BP.Loc" = loci.mBIC[2:(length(loci.mBIC) - 1)]))
    } else {
      return(paste("No Break-Points are Estimated")) 
      }
    
    } else if(distyp == 2){
      
      Melite <- M*rho
      L <- length(data[, 1])
      L0 <- 1      
      k <- seq(0, Nmax, 1)
      
      if(parallel == TRUE & .Platform$OS.type == "windows"){

            cl <- makeCluster(parallel::detectCores(), type="SOCK") 
            clusterExport(cl, c("ce.simnormal", "normrand", "mBIC"), envir = environment())
            clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "b"), envir = environment())  
            registerDoParallel(cl)
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal(k, data, h, L0, L, M, Melite, eps, a, b)
            stopCluster(cl)   
            
        } else if (parallel == TRUE & .Platform$OS.type == "unix"){
            
            registerDoParallel(parallel::detectCores()) 
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal(k, data, h, L0, L, M, Melite, eps, a, b)
          
        } else { 
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.simnormal(k, data, h, L0, L, M, Melite, eps, a, b)
        }      
      
        k <- seq(0, Nmax, 1)    
        mBic.summary <- c()
        for(i in 1:length(sim)){
          mBic.summary[i] <- sim[[i]]$mBIC
        }
    
      loci.mBIC <- sim[[k[which(mBic.summary==max(mBic.summary))] + 1]]$loci
      
        if(length(loci.mBIC) >= 3) {
          return(list("No.BPs" = length(loci.mBIC) - 2, "BP.Loc" = loci.mBIC[2:(length(loci.mBIC) - 1)]))
        } else {
          return(paste("No Break-Points are Estimated")) 
        }          
      } 
  }
}
