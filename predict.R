#############################################
########### funcio pred MDP #################

predict.mdp<-function(modelFit,newdata){
  data=as.matrix(modelFit$input[["data"]])
  newdata=as.matrix(newdata)
  k=modelFit$k
  g=modelFit$g
  D=modelFit$D
  lambda=modelFit$input[["lambda"]]
  comb=modelFit$comb
  ind=modelFit$ind
  class.assign=modelFit$class.assign
  start=dim(D)[1]
  vec.new=(start+1):(start+dim(newdata)[1])
  C.pool=modelFit$C.pool
  
  if(is.null(modelFit$input[["distance"]]) && is.null(modelFit$input[["kernel"]])){
    modelFit$input[["distance"]]="mahalanobis"
  } else if(!is.null(modelFit$input[["distance"]])){
    distance=modelFit$input[["distance"]]
    if(distance=="mahalanobis"){
      for(i in 1:dim(newdata)[1]){
        d=c()
        D=cbind(D,c(0))
        D=rbind(D,c(0))
        for(j in 1:dim(data)[1]){
          d=c(d,mahalanobis(newdata[i,], data[j,],C.pool))
          D[vec.new[i],j]=d[j]
        }
      }
    } else if(distance=="euclidean"){
      for(i in 1:dim(newdata)[1]){
        d=c()
        D=cbind(D,c(0))
        D=rbind(D,c(0))
        for(j in 1:dim(data)[1]){
          d=c(d,as.numeric(dist(rbind(newdata[i,], data[j,]))))
          D[vec.new[i],j]=d[j]
        }
      }
    } else if(distance=="gower"){
      if(!require(cluster)){
        install.packages("cluster")
        library(cluster)
      }
      data2=rbind(data,newdata)
      D=as.matrix(daisy(data2, metric = "gower"))
    } 
  } else if(!is.null(modelFit$input[["kernel"]])){
    if(!require(kernlab)){
      install.packages("kernlab")} 
    kernel=modelFit$input[["kernel"]]
    if(substr(kernel, 1, 1)=="g"){
      sigma=as.numeric(substr(kernel, 5, nchar(kernel)))
      K.rbfdot<-function(x,y){
        rbf <- rbfdot(sigma = sigma)
        dist = 2*(1-rbf(x,y))
        sqrt(dist)
      }      
      for(i in 1:dim(newdata)[1]){
        d=c()
        D=cbind(D,c(0))
        D=rbind(D,c(0))
        for(j in 1:dim(data)[1]){
          d=c(d,K.rbfdot(newdata[i,], data[j,]))
          D[vec.new[i],j]=d[j]
        }
      }
    } else if(substr(kernel, 1, 1)=="l"){
      sigma=as.numeric(substr(kernel, 5, nchar(kernel)))
      K.laplacedot<-function(x,y){
        lap <- laplacedot(sigma = sigma)
        dist = 2*(1-lap(x,y))
        sqrt(dist)
      }      
      for(i in 1:dim(newdata)[1]){
        d=c()
        D=cbind(D,c(0))
        D=rbind(D,c(0))
        for(j in 1:dim(data)[1]){
          d=c(d,K.laplacedot(newdata[i,], data[j,]))
          D[vec.new[i],j]=d[j]
        }
      }
      
    } else if(substr(kernel, 1, 1)=="p"){
      degree=as.numeric(substr(kernel, 5, nchar(kernel)))
      K.polydot<-function(x,y){
        lap <- polydot(degree = degree, scale = 1, offset = 1)
        dist = lap(x,x)+lap(y,y)-2*lap(x,y)
        sqrt(dist)
      } 
      for(i in 1:dim(newdata)[1]){
        d=c()
        D=cbind(D,c(0))
        D=rbind(D,c(0))
        for(j in 1:dim(data)[1]){
          d=c(d,K.polydot(newdata[i,], data[j,]))
          D[vec.new[i],j]=d[j]
        }
      }
      
    }
  }

  r=numeric(dim(newdata)[1])
  R=rep(0,k)
  ER=rep(0,k)
  mat=matrix(data=0,nrow=k,ncol=k)
  t=1
  for(row in vec.new){
    h=matrix(data=0,nrow=start,ncol=k)
    col=1
    while(col<=k){
      for(i in 1:dim(h)[1]){
        for(j in which(g[,col]!=0)){
          if(D[row,i]<D[row,j]){h[i,col]=h[i,col]+1}
        }
      }
      col<-col+1
    }

    prod=0
    vec.q<-rep(0,k)
    col=1
    
    while(col<=k){
      for(i in which(g[,col]!=0)){
        vec.k<-rep(0,k)
        for(j in 1:dim(h)[2]){
          if(h[i,j]>=lambda-ind(i,j)){vec.k[j]=comb(h[i,j],lambda-ind(i,j))}
        }
        
        if(sum(vec.k==0)==0){
          prod = prod+prod(vec.k)
        }
      }
      vec.q[col]=prod; prod=0
      col<-col+1
    }
    
    j=which.max(vec.q)
    if(sum(vec.q[j]==vec.q[-j],na.rm = T)!=0){
      w = which(vec.q[j]==vec.q)
      j=sample(w,1)
    }
    
    r[t]=class.assign[j]
    t<-t+1
  }
  if(!is.null(modelFit$cm.levels)){
    n=length(modelFit$cm.levels[[1]])
    for(i in 1:n){
      r[which(r==i)]<-modelFit$cm.levels[[1]][i]
    }
    r=factor(r)
  }
  return(r)
}
  

