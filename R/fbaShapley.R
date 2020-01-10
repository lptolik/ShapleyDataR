#' Shapley values for FBA model.
#'
#' \code{fbaShapley} return object of class \code{shapley}, which contains
#' \describe{
#'   \item{reId}{Reaction IDs as provided by argument}
#'   \item{shapley}{Shapley values, vector with reId as names}
#'   \item{sd}{Standard deviation for the Shapley values}
#'   \item{p}{Probability that the \code{shapley} contains real Shapley values.}
#' }
#'
#'
#'
#' @param mod FBA model in Sybil format.
#' @param reId list of reaction Ids to be used in calculations.
#' @param tol tolerance for Shapley value convergence
#' @param logName name of the file to append log records to
#' @param tmpSave frequency of tmp results saving
#' @param perfTolerance tolerance for fail early check
#' @param cacheDepth how many first knockout have to be cached
#'
#' @return object of class \code{shapley}
#' @export
#'
#' @importFrom sybil optimizeProb
#' @importFrom stats qnorm
#' @import hash
fbaShapley<-function(mod,reId,tol=0.05,logName='fbaShapley.log',tmpSave=500,perfTolerance=1e-5,cacheDepth=3){
  convTol<-0.1
  alph<-c(0.01,0.05,0.1)
  Z<-qnorm(alph,lower.tail = FALSE)
  N<-length(reId)
  tmpFile<-paste0(sub('.log','',logName),format(Sys.time(), "%Y%m%d%H%M"),'.RData')
  phi<-list()
  m2<-list()
  val<-list()
  permL<-list()
  jStop<-c()
  mTot<- optimizeProb(mod,algorithm = "fba",retOptSol = FALSE)
  vTot<-mTot$obj
  v<-rep(0.0,N)
  vNull<-0
  t<-1
  phi[[1]]<-rep(0.0,N)
  val[[1]]<-rep(0.0,N)
  m2[[1]]<-rep(0.0,N)
  objHash<-hash()
  fmt<-paste0('%0',floor(log10(length(reId))+1),'d')
  while(!convCriteria(phi,convTol)){
    t<-t+1
    if(t<=100){
      cat(format(Sys.time(), "%b %d %X"),'t=',t,'\n',append = TRUE,file = logName)
    }else{
      tolV<-sum(abs(phi[[t-1]]-phi[[t-100]])/(1e-5+abs(phi[[t-1]])))
      if(t%%tmpSave==0){
        sd<-m2[[t-1]]/(t-2)
        e<-sapply(Z,function(.x)sqrt((.x^2*sd)/(t-1)))
        save(phi,t,N,vTot,v,val,permL,m2,sd,e,Z,alph,perfTolerance,vNull,tolV,mod,reId,jStop,objHash,file = tmpFile)
        cat(format(Sys.time(), "%b %d %X"),'t=',t,'tol=',tolV,' err (1%,5%,10%)=',apply(e,2,max),'\n',append = TRUE,file = logName)
      }else{
        cat(format(Sys.time(), "%b %d %X"),'t=',t,'tol=',tolV,'\n',append = TRUE,file = logName)
      }
    }
    perm<-makePerm(N)
    v<-rep(0.0,N)
    newRes<-vTot
    belowIdx<-0
    for (j in (1:N)){
      oldRes<-newRes
      if(j<=cacheDepth){
        key<-paste(sprintf(fmt,sort(perm[1:j])),collapse = '')
        if(!has.key(key,objHash)){
      objHash[[key]]<-calcObj(mod,reId[perm[1:j]])
      }
      newRes<-objHash[[key]]
      }else{
        newRes<-calcObj(mod,reId[perm[1:j]])
      }
      if(abs(vNull-newRes)< perfTolerance){
        belowIdx<-belowIdx+1
      }else{
        belowIdx<-0
      }
      if(belowIdx>5){
        v[j:N]<-0
        jStop[t]<-j
        break()
      }
      v[j]<-newRes-oldRes
    }
    val[[t]]<-rep(0.0,N)
    val[[t]]<-v
    permL[[t]]<-perm
    phi[[t]]<-rep(0.0,N)
    phi[[t]][perm]<-phi[[t-1]][perm]+(v-phi[[t-1]][perm])/t
    m2[[t]]<-rep(0.0,N)
    m2[[t]][perm]<-m2[[t]][perm]+(v-phi[[t-1]][perm])*(v-phi[[t]][perm]) # don't forget to divide by T at the end
  }
  sd<-m2[[t]]/(t-1)
  e<-sapply(Z,function(.x)sqrt((.x^2*sd)/(t)))
  save(phi,t,N,vTot,v,val,permL,m2,sd,e,Z,alph,perfTolerance,vNull,tolV,mod,reId,jStop,objHash,file = tmpFile)
  idx<-match(reId,mod@react_id)
  shpl<--1*phi[[t]]
  res<-list(reId=reId,
            name=mod@react_name[idx],
            shapley=shpl,
            sd=m2[[t]]/(t-1),
            err01=e[,1],
            err05=e[,2],
            err10=e[,3])
  class(res)<-'shapley'
  clear(objHash)
  rm(objHash)
  return(res)
}

calcObj<-function(mod,react){
  bound<-rep(0,length(react))
  model<-optimizeProb(mod,algorithm = "fba",retOptSol = FALSE,react=react,lb=bound,ub=bound)
  return(model$obj)
}