ce.beta <-
function(N,data,h,L0,L,M,Nelite,eps,a){
  
  if (N==0){
    seql<-c(0,L)
    mBic.full<-modi.BICnorm.new(seql,data,0,L,h)
    
    return(list(locis=c(1,L+1),mBIC=mBic.full))
    rm(mBic.full,seql)
    
  } else {
    
  ########################Parameter initialization######################################################  
  new_para<-array(1,dim=c(2,N))
  ######################################################################################################  
  modbic<-c()
  
  k<-0
  repeat
  {
    k<-k+1
    ch<-array(0,dim=c(M,N+2))       
    ch[,1]<-c(1)                    
    ch[,N+2]<-c(L+1)    
    ch[,(2:(N+1))]<-apply(new_para,2,betarand,L0,L,M)   
    ch<-t(apply(ch,1,sort))                       
    mod_bic<-apply(ch,1,modi.BICnorm.new,data,N,L,h)
    ch<-cbind(ch,mod_bic)                         
    ch<-ch[order(ch[,(N+3)],decreasing=TRUE),]  
    nelitesmpl<-ch[1:Nelite,]                     
    modbic[k]<-nelitesmpl[1,(N+3)] 
    
    new_par_n<-array(0,dim=c(2,N))
    new_par_n[1,]<-apply(as.matrix(nelitesmpl[,(2:(N+1))]),2,mean)
    new_par_n[2,]<-apply(as.matrix(nelitesmpl[,(2:(N+1))]),2,var)   
    
    new_par_new<-array(0,dim=c(2,N))
    new_par_new[1,]<-apply(new_par_n,2,fun_alpha,L0,L)
    new_par_new[2,]<-apply(new_par_n,2,fun_beta,L0,L)
    new_para<-a*new_par_new + (1-a)*new_para
    
    mad<-apply(as.matrix(nelitesmpl[,(2:(N+1))]),2,mad)
    
    if(max(mad)<=eps){break}
  }
  return(list(locis=ch[1,(1:(N+2))], mBIC=modbic[k]))
  }
}