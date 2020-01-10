library(glmnet)
library(ShapleyDataR)
library(sybil)
library(sybilcycleFreeFlux)
logFile<-'/Users/lptolik/Documents/Projects/ShapleyDataR/FBAfilesReportPengFull.log'
mod<-readRDS(file='/Users/lptolik/Documents/Projects/Penicillium/pen_model/cleanModel/pen_model.rds')
ex <- findExchReact(mod)
idxUp<-grep('IN',ex@react_id)
modU <- changeBounds(mod, ex[idxUp], lb = rep(0,length(idxUp)), ub = rep(0,length(idxUp)))
#Penicillin medium
modMedP <- changeBounds(modU, ex[c("o2IN", "glcIN",'nh3IN','slfIN','paaIN')], ub = rep(1000,5))
#Growth medium
modMedGr <- changeBounds(modU, ex[c('o2IN','glcIN','piIN','slfIN','pimIN','thmIN','nh3IN')], ub = rep(1000,7))

modPeng<-changeObjFunc(mod,c('pengOUT'),c(1))
cat(format(Sys.time(), "%b %d %X"),'Start\n',append = TRUE,file = logFile)
idxRe<-which((uppbnd(modPeng)-lowbnd(modPeng)>0)&!grepl('^(pengOUT|bmOUT)$',modPeng@react_id))
reId<-modPeng@react_id[idxRe]
res<-fbaShapley(modPeng,reId = reId,logName = logFile)
# N<-modPeng@react_num
# convTol<-0.1
# tol<-0.01
# phi<-list()
# sd<-list()
# val<-list()
# permL<-list()
# mTot<- optimizeProb(modPeng,algorithm = "fba",retOptSol = FALSE)
# vTot<-mTot$obj
# v<-rep(0.0,N)
# vNull<-0
# perfTolerance<-tol*vTot
# t<-1
# phi[[t]]<-rep(0.0,N)
# val[[t]]<-rep(0.0,N)
# sd[[t]]<-rep(0.0,N)
# while(!convCriteria(phi,convTol)){
#   t<-t+1
#   if(t<=101){
#     cat(format(Sys.time(), "%b %d %X"),'t=',t,'\n')
#   }else{
#     tolV<-sum(abs(phi[[t-1]]-phi[[t-101]])/(1e-5+abs(phi[[t-1]])))
#     cat(format(Sys.time(), "%b %d %X"),'t=',t,'tol=',tolV,'\n',append = TRUE,file = logFile)
#   if(t%%10==0){
#     save(phi,t,N,vTot,v,val,permL,sd,perfTolerance,vNull,tolV,file = 'fbaShapley.RData')
#   }
#   }
#   perm<-makePerm(N)
#   v<-rep(0.0,N)
#   newRes<-vNull
#   belowIdx<-0
#   for (j in (1:N)){
#     oldRes<-newRes
#     bound<-rep(0,j)
#     model<-optimizeProb(modPeng,algorithm = "fba",retOptSol = FALSE,react=perm[1:j],lb=bound,ub=bound)
#     newRes<-model$obj
#     #cat(format(Sys.time(), "%b %d %X"),t,"Instance:",j,'; newRes=',newRes,"\n")
#     if(abs(vNull-newRes)< perfTolerance){
#       belowIdx<-belowIdx+1
#       #cat(format(Sys.time(), "%b %d %X"),t,"Instance:",j,'',belowIdx,"\n")
#     }else{
#       belowIdx<-0
#     }
#     if(belowIdx>5){
#       v[j:N]<-0
#       #if(t%%100==0){
#       cat(format(Sys.time(), "%b %d %X"),t,"Tolerance break:",j,"\n")
#       #}
#       break()
#     }
#     v[j]<-newRes-oldRes
#   }
#   val[[t]]<-rep(0.0,N)
#   val[[t]]<-v
#   permL[[t]]<-perm
#   phi[[t]]<-rep(0.0,N)
#   phi[[t]][perm]<-phi[[t-1]][perm]+(v-phi[[t-1]][perm])/t
#   sd[[t]]<-rep(0.0,N)
#   sd[[t]][perm]<-sd[[t]][perm]+(v-phi[[t-1]][perm])*(v-phi[[t]][perm])
# }
#
# save.image(file='fba.shapley.rdata')
