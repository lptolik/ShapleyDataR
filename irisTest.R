library(glmnet)
library(datasets)
library(ShapleyDataR)
library(doParallel)
data(iris)

y <- as.numeric(iris[,5])
X <- iris[y!=1, 1:4]
y <- y[y!=1]-2

data<-iris[y!=1,]

testIdx<-sample.int(dim(data)[1],dim(data)[1]*0.1)

trainData<-data[-testIdx,]
testData<-data[testIdx,]

alg<-function(data){
  if(dim(data)[1]<=0){
    return(NULL)
  }
  if(length(unique(data[,5]))<2|
     min(table(as.numeric(data[,5])))<2){
    return(NULL)
  }
  x<-data[,1:4]
  y<-as.numeric(data[,5])
  model_lambda <- glmnet(as.matrix(x), as.factor(y),alpha=1, family="binomial", type.measure="class")
return(model_lambda)
}

validF<-function(model,test){
  if(is.null(model)){
    return(max(table(test[,5])/dim(test)[1])) #(runif(1))
  }
  model_lambda<-model
  best_s  <- model_lambda$lambda.1se
  x_test<-test[,1:4]
  y_test<-as.numeric(test[,5])
  pred <- as.numeric(predict(model_lambda, newx=as.matrix(x_test), type="class" , s=best_s))
  res=sum(y_test==pred)/NROW(pred)
  return(res)
}

res<-dataShapleyParallel(data,alg,validF,testData,numCores = 8)
