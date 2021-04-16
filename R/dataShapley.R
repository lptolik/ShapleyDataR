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
  val<-list()
  sd<-list()
  permL<-list()
  model<-A(D)
  vTot<-V(model,T)
  perfTolerance<-tol*vTot
  t<-1
  phi[[1]]<-rep(0.0,N)
  sd[[1]]<-rep(0.0,N)
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
    sd[[t]]<-rep(0.0,N)
    sd[[t]][perm]<-sd[[t]][perm]+(v-phi[[t-1]][perm])*(v-phi[[t]][perm])
    val[[t]]<-v
    permL[[t]]<-perm
  }
  return(list(phi=phi,val=val,perm=permL))
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
  sd<-list()
  val<-list()
  alph<-c(0.01,0.05,0.1)
  Z<-qnorm(alph,lower.tail = FALSE)
  m2 <- list()
  permL<-list()
  model<-A(D)
  tolMS<-tolMeanScore(model,V,T)
  vTot<-tolMS$mean #V(D,model,T)
  v<-rep(0.0,N)
  vNull<-V(NULL, T)
  perfTolerance<-tol*vTot
  t<-1
  phi[[t]]<-rep(0.0,N)
  val[[t]]<-rep(0.0,N)
  sd[[t]]<-rep(0.0,N)
  m2[[t]]<-rep(0.0,N)
  while(!convCriteria(phi,convTol)){
    t<-t+1
    if(t<=101){
      cat(format(Sys.time(), "%b %d %X"),'t=',t,'\n')
    }else if(t%%100==0){
      sd<-m2[[t-1]]/(t-2)
      e<-sapply(Z,function(.x)sqrt((.x^2*sd)/(t-1)))
      tolV<-sum(abs(phi[[t-1]]-phi[[t-101]])/(1e-5+abs(phi[[t-1]])))
      cat(format(Sys.time(), "%b %d %X"),'t=',t,'tol=',tolV,'\n')
      save(phi,t,N,vTot,v,val,permL,sd,perfTolerance,vNull,tolMS,m2,e,file = 'tmpShapley.RData')
    }
    perm<-makePerm(N)
    v<-rep(0.0,N)
    newRes<-vNull
    belowIdx<-0
    for (j in (1:N)){
      oldRes<-newRes
      model<-A(D[perm[1:j],])
      if (is.null(model)) {
        newRes <- vNull
        cat(format(Sys.time(), "%b %d %X"), "Model is null, j =", j, "\n")
      } else {
        newRes<-V(model,T)
      }
      if(abs(vTot-newRes)< perfTolerance){
        belowIdx<-belowIdx+1
      }else{
        belowIdx<-0
      }
      if(belowIdx>5){
        v[j:N]<-0
        cat(format(Sys.time(), "%b %d %X"),t,"Tolerance break:",j,"\n")
        break()
      }
      v[j]<-newRes-oldRes
    }
    val[[t]]<-rep(0.0,N)
    val[[t]][perm]<-v
    permL[[t]]<-perm
    phi[[t]]<-rep(0.0,N)
    phi[[t]][perm]<-phi[[t-1]][perm]+(v-phi[[t-1]][perm])/t
    m2[[t]]<-rep(0.0,N)
    m2[[t]][perm]<-m2[[t-1]][perm]+(v-phi[[t-1]][perm])*(v-phi[[t]][perm])
  }
  sd<-m2[[t]]/(t-1)
  e<-sapply(Z,function(.x)sqrt((.x^2*sd)/(t)))
  tolV<-sum(abs(phi[[t-1]]-phi[[t-101]])/(1e-5+abs(phi[[t-1]])))
  cat(format(Sys.time(), "%b %d %X"),'t=',t,'tol=',tolV,'\n')
  save(phi,t,N,vTot,v,val,permL,perfTolerance,vNull,tolMS,m2,e,file = 'tmpShapley.RData')
  return(list(phi=phi,
              val=val,
              perm=permL,
              sd=sd,
              err01=e[,1],
              err05=e[,2],
              err10=e[,3]))
}

