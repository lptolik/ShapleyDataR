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


#' Combine results obtained from separate workers
#'
#' @param results list of combined results
#' @param x list of results from one worker
#'
#' @return new list of combined results
#' @export
combResults <-function(results, x){
  i <- x$i
  perm <- x$perm
  v <- x$v
  results$val[[i]] <- x$val
  results$val[[i]][perm] <- v
  results$permL[[i]] <- perm
  results$phi[[i]] <- x$phi
  results$m2[[i]] <- x$m2
  return(results)

}

