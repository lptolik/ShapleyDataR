#' Truncated Monte Carlo Shapley
#'
#' @param D dataset of N points
#' @param A learning algorithm
#' @param V performance score
#'
#' @return Shapley value of training points
#' @export
#'
dataShapley<-function(D,A,V){
  N<-dim(D)[1]
  phi<-rep(0.0,N)
  vTot<-V(D,A)
  t<-0
  while(!convCriteria()){
    t<-t+1
    perm<-makePerm(D)
    v<-list()
    v[[1]]<-V(perm[FALSE,],A)
    for (j in (1:N)){
      if(abs(vTot-v[j])< perfTolerance){
        v[[j+1]]<-v[[j]]
      }else{
        v[[j+1]]<-V(perm[1:j,],A)
      }
    }

  }
}