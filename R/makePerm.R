#' Check convergence criteria
#'
#' @param phi shapley values
#' @param convTol tolerance
#'
#' @return TRUE if criteria met, FALSE otherwise
#' @export
convCriteria<-function(phi,convTol=0.05){
  if(length(phi)<=100){
    return(FALSE)
  }else{
    t<-length(phi)
    return(sum(abs(phi[[t]]-phi[[t-100]])/(1e-5+abs(phi[[t]])))< convTol)
  }
 return(TRUE)
}

convCriteriaErr<-function(phi,convTol=0.05,gap=2){
  if(length(phi)<=100){
    return(FALSE)
  }else{
    r<-getPhi(phi)
    m<-r$phi
    sd<-r$sd
    t<-r$N
    mN<-m/sum(m)
    idx<-which(abs(mN)>(1e-1/length(m)))
    Z<-qnorm(convTol,lower.tail = FALSE)
    e<-sqrt((Z^2*sd)/(t))
    ide<-which(e==0)
    e[ide] <- 1e-3/length(m)
    cat(format(Sys.time(), "%b %d %X"),'convergence: n=',length(idx),', nonconv=',length(which(abs(m[idx]/e[idx])>gap)),'\n')
    cat(format(Sys.time(), "%b %d %X"),'convergence: mN=[',summary(mN),'] \n')
    cat(format(Sys.time(), "%b %d %X"),'convergence: ratio=[',summary(abs(m[idx]/e[idx])),']\n')
    return(all(abs(m[idx]/e[idx])>gap))
  }
  return(TRUE)
}


getPhi<-function(val){
  t<-length(val)
  p<-do.call(rbind,val)
  rS<-rowSums(p)
  m<-apply(p, 2, mean)
  sd<-apply(p,2,sd)
  return(list(phi=m,sd=sd,N=t))
}

#' Make data point permutations.
#'
#' THe function is extracted to be able to apply group permutations later.
#'
#' @param N length of datasets
#'
#' @return indices of permuted dataset
#' @export
makePerm<-function(N){
 return(sample.int(N))
}

#' @importFrom stats IQR median sd
NULL

tolMeanScore<-function(model,V,T){
  res<-c()
  N<-dim(T)[1]
  for(i in 1:100){
    perm<-sample.int(N,size = N,replace = TRUE)
    res[i]<-V(model,T[perm,])
  }
  return(list(min=min(res),mean=mean(res),std=sd(res),median=median(res),iqr=IQR(res),max=max(res)))
}
