#' Truncated Monte Carlo Shapley
#'
#' @param D dataset of N points
#' @param A learning algorithm
#' @param V performance score
#' @param T test dataset
#' @param tol trancation tolerance
#' @param convTol convergence tolerance
#'
#' @return Shapley value of training points
#' @export
#'
dataShapley<-function(D,A,V,T,tol=0.05,convTol=tol){
  N<-dim(D)[1]
  phi<-list()
  vTot<-V(D,A,T)
  perfTolerance<-tol*vTot
  t<-1
  phi[[t]]<-rep(0.0,N)
  while(!convCriteria(phi,convTol)){
    t<-t+1
    if(t<=101){
    cat(format(Sys.time(), "%b %d %X"),'t=',t,'\n')
    }else{
      tolV<-sum(abs(phi[[t-1]]-phi[[t-101]])/(1e-5+abs(phi[[t-1]])))
      cat(format(Sys.time(), "%b %d %X"),'t=',t,'tol=',tolV,'\n')
    }
    perm<-makePerm(N)
    vNull<-V(D[FALSE,],A,T)
    v<-c()
    v[1]<-vNull
    for (j in (1:N)){
      if(abs(vTot-v[j])< perfTolerance){
        v[j+1]<-v[j]
      }else{
        v[j+1]<-V(D[perm[1:j],],A,T)
      }
    }
    phi[[t]]<-rep(0.0,N)
    phi[[t]][perm]<-phi[[t-1]][perm]*(t-1)/t+(v[2:N+1]-v[1:N])/t
  }
  return(phi)
}

#' Data Shapley with different truncation rule
#'
#' @param D dataset of N points
#' @param A learning algorithm
#' @param V performance score
#' @param T test dataset
#' @param tol trancation tolerance
#' @param convTol convergence tolerance
#'
#' @return Shapley value of training points
#' @export
#'
dataShapleyI5<-function(D,A,V,T,tol=0.01,convTol=tol*5){
  N<-dim(D)[1]
  phi<-list()
  vTot<-V(D,A,T)
  v<-rep(0.0,N)
  vNull<-V(D[FALSE,],A,T)
  perfTolerance<-tol*vTot
  t<-1
  phi[[t]]<-rep(0.0,N)
  while(!convCriteria(phi,convTol)){
    t<-t+1
    if(t<=101){
      cat(format(Sys.time(), "%b %d %X"),'t=',t,'\n')
    }else if(t%%1000==0){
      tolV<-sum(abs(phi[[t-1]]-phi[[t-101]])/(1e-5+abs(phi[[t-1]])))
      cat(format(Sys.time(), "%b %d %X"),'t=',t,'tol=',tolV,'\n')
      save(phi,t,N,vTot,v,perm,perfTolerance,vNull,file = 'tmpShapley.RData')
    }
    perm<-makePerm(N)
    vNull<-V(D[FALSE,],A,T)
    v<-rep(0.0,N)
    newRes<-vNull
    belowIdx<-0
    for (j in (1:N)){
      oldRes<-newRes
      newRes<-V(D[perm[1:j],],A,T)
      if(abs(vTot-newRes)< perfTolerance){
        belowIdx<-belowIdx+1
      }else{
        belowIdx<-0
      }
      if(belowIdx>5){
        v[j+1:N]<-0
        if(t%%100==0){
        cat(format(Sys.time(), "%b %d %X"),t,"Tolerance break:",j,"\n")
        }
        break()
      }
      v[j]<-newRes-oldRes
    }
    phi[[t]]<-rep(0.0,N)
    phi[[t]][perm]<-phi[[t-1]][perm]*(t-1)/t+v/t
  }
  return(phi)
}

