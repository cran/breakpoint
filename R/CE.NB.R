CE.NB <-
function(data, Nmax=10, eps=0.01, rho=0.05, M=200, h=5, a=0.8, b=0.8, distyp = 1, parallel=FALSE){
  
  if(is.data.frame(data) == "FALSE" | is.null(dim(data)[2])) {
    print("Error in data: dataframe only")                                                           
  } else if(dim(data)[2] != 1) {
    print("Error in data: single column dataframe only")                               
  } else { 
    
    if(distyp == 1){
      Melite <- M*rho
      L <- length(data[, 1])
      L0 <- 1      
      k <- seq(0, Nmax, 1)
      r <- suppressWarnings(try(fitdistr(data[, 1], "negative binomial")[[1]][[1]], silent = T)) 
      
      if(parallel == TRUE & .Platform$OS.type == "windows"){

          cl <- makeCluster(parallel::detectCores(), type = "SOCK") 
          clusterExport(cl, c("ce.4betaNB", "betarand", "fun.alpha", "fun.beta", "BICnb", "llhoodnb", "logliknb"), envir = environment())
          clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "r"), envir = environment())  
          registerDoParallel(cl)
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.4betaNB(k, data, h, L0, L, M, Melite, eps, a, r)
          stopCluster(cl)   
          
        } else if (parallel == TRUE & .Platform$OS.type == "unix"){

          registerDoParallel(parallel::detectCores()) 
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.4betaNB(k, data, h, L0, L, M, Melite, eps, a, r)
                  
        } else { 
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.4betaNB(k, data, h, L0, L, M, Melite, eps, a, r)
          }
      
      k <- seq(0, Nmax, 1)    
      Bic.summary <- c()
      for(i in 1 : length(sim)){
        Bic.summary[i] <- sim[[i]]$BIC
      }
      
      loci.BIC <- sim[[k[which(Bic.summary == min(Bic.summary))] + 1]]$loci
      
      if(length(loci.BIC) >= 3) {
        return(list("No.BPs" = length(loci.BIC) - 2,"BP.Loc" = loci.BIC[2 : (length(loci.BIC) - 1)]))
      } else {
        return(paste("No Break-Points are Estimated")) 
      }
    
    } else if(distyp == 2){
      
      Melite <- M*rho
      L <- length(data[, 1])
      L0 <- 1      
      k <- seq(0, Nmax, 1)
      r <- suppressWarnings(try(fitdistr(data[, 1], "negative binomial")[[1]][[1]], silent = T)) 
      
      if(parallel == TRUE & .Platform$OS.type == "windows"){

          cl <- makeCluster(parallel::detectCores(), type = "SOCK") 
          clusterExport(cl, c("ce.simNormalNB", "normrand", "BICnb", "llhoodnb", "logliknb"), envir = environment())
          clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "b", "r"), envir=environment())  
          registerDoParallel(cl)
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simNormalNB(k, data, h, L0, L, M, Melite, eps, a, b, r)
          stopCluster(cl)   

          } else if (parallel == TRUE & .Platform$OS.type == "unix"){
            
            registerDoParallel(parallel::detectCores()) 
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simNormalNB(k, data, h, L0, L, M, Melite, eps, a, b, r)
          
        } else {
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.simNormalNB(k, data, h, L0, L, M, Melite, eps, a, b, r)
        }
      
      k <- seq(0, Nmax, 1)    
      Bic.summary <- c()
      for(i in 1 : length(sim)){
        Bic.summary[i] <- sim[[i]]$BIC
      }
      
      loci.BIC <- sim[[k[which(Bic.summary == min(Bic.summary))] + 1]]$loci
      
      if(length(loci.BIC) >= 3) {
        return(list("No.BPs" = length(loci.BIC) - 2,"BP.Loc" = loci.BIC[2 : (length(loci.BIC) - 1)]))
      } else {
        return(paste("No Break-Points are Estimated")) 
      }          
      }
  }
}
