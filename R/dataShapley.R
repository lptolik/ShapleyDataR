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
  model<-A(D)
  vTot<-V(model,T)
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
    model<-A(D[FALSE,])
    vNull<-V(model,T)
    v<-c()
    v[1]<-vNull
    for (j in (1:N)){
      if(abs(vTot-v[j])< perfTolerance){
        v[j+1]<-v[j]
      }else{
        model<-A(D[perm[1:j],])
        v[j+1]<-V(model,T)
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
  model<-A(D)
  tolMS<-tolMeanScore(model,V,T)
  vTot<-tolMS$mean #V(D,model,T)
  v<-rep(0.0,N)
  vNull<-V(NULL,T)
  perfTolerance<-tol*vTot
  t<-1
  phi[[t]]<-rep(0.0,N)
  while(!convCriteria(phi,convTol)){
    t<-t+1
    if(t<=101){
      cat(format(Sys.time(), "%b %d %X"),'t=',t,'\n')
    }else if(t%%100==0){
      tolV<-sum(abs(phi[[t-1]]-phi[[t-101]])/(1e-5+abs(phi[[t-1]])))
      cat(format(Sys.time(), "%b %d %X"),'t=',t,'tol=',tolV,'\n')
      save(phi,t,N,vTot,v,val,perm,perfTolerance,vNull,tolMS,file = 'tmpShapley.RData')
    }
    perm<-makePerm(N)
    #model<-A(D[FALSE,])
    vNull<-V(NULL,T)
    v<-rep(0.0,N)
    newRes<-vNull
    belowIdx<-0
    for (j in (1:N)){
      oldRes<-newRes
      model<-A(D[perm[1:j],])
      newRes<-V(model,T)
      if(abs(vTot-newRes)< perfTolerance){
        belowIdx<-belowIdx+1
        #cat(format(Sys.time(), "%b %d %X"),t,"Instance:",j,'',belowIdx,"\n")
      }else{
        belowIdx<-0
      }
      if(belowIdx>5){
        v[j:N]<-0
        #if(t%%100==0){
        cat(format(Sys.time(), "%b %d %X"),t,"Tolerance break:",j,"\n")
        #}
        break()
      }
      v[j]<-newRes-oldRes
    }
    val[[t]]<-rep(0.0,N)
    val[[t]][perm]<-v
    phi[[t]]<-rep(0.0,N)
    phi[[t]][perm]<-phi[[t-1]][perm]*(t-1)/t+v/t
  }
  tolV<-sum(abs(phi[[t-1]]-phi[[t-101]])/(1e-5+abs(phi[[t-1]])))
  cat(format(Sys.time(), "%b %d %X"),'t=',t,'tol=',tolV,'\n')
  save(phi,t,N,vTot,v,val,perm,perfTolerance,vNull,tolMS,file = 'tmpShapley.RData')
  return(phi)
}

