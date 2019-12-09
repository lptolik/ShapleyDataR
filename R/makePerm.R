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

tolMeanScore<-function(model,V,T){
  res<-c()
  N<-dim(T)[1]
  for(i in 1:100){
    perm<-sample.int(N,size = N,replace = TRUE)
    res[i]<-V(model,T[perm,])
  }
  return(list(bin=min(res),mean=mean(res),std=sd(res),median=median(res),iqr=IQR(res),max=max(res)))
}
