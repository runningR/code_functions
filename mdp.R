mdp<-function(Class,data,lambda,distance=NULL,prior=NULL,info.pred=NULL,kernel=NULL){
  input=mget(names(formals()),sys.frame(sys.nframe()))
  if(is.numeric(Class)){Class=as.numeric(Class); lev=FALSE
  } else { lev=TRUE
  if(!is.factor(Class)){
    Class=as.matrix(Class)
    Class=as.factor(Class)
  }
  cm.levels=list(levels(Class),levels(Class))
  Class=as.numeric(Class)
  }
  data = as.matrix(data)
  k=nrow(table(Class))
  K.laplacedot=NULL; K.polydot=NULL; K.rbfdot=NULL
  C.pool=NULL
  if(missing(distance) && missing(kernel)){distance="mahalanobis"}
  
  if(!is.null(distance)){
    if(distance=="mahalanobis"){
      D=matrix(data=0,nrow=dim(data)[1],ncol=dim(data)[1])
      Data=data.frame(data); Data$class=Class
      d=c()
      L.cov=vector("list",k)
      for(cl in 1:k){
        L.cov[[cl]]=cov(subset(Data,Data$class==cl,select = -c(class)))
      }
      C.pool=matrix(0,nrow=dim(data)[2],ncol=dim(data)[2])
      for (i in 1:ncol(data)){
        for(j in 1:ncol(data)){
          for(cl in 1:k){
            C.pool[i,j]=C.pool[i,j]+L.cov[[cl]][i,j]*(as.numeric(table(Class)[cl])-1)/(dim(data)[1]-k)
          }
        }
      }
      
      for(i in 1:dim(data)[1]){
        d=c()
        for(j in 1:dim(data)[1]){
          d=c(d,mahalanobis(data[i,], data[j,],C.pool))
          D[i,j]=d[j]
        }
      }
    
      } else if(distance=="euclidean"){
      D=matrix(data=NA,nrow=dim(data)[1],ncol=dim(data)[1])
      d=c()
      for(i in 1:dim(data)[1]){
        d=c()
        for(j in 1:dim(data)[1]){
          d=c(d,as.numeric(dist(rbind(data[i,],data[j,]))))
          D[i,j]=d[j]
        }
      }
    
      } else if(distance=="gower"){
      if(!require(cluster)){
        install.packages("cluster")
        library(cluster)
      }
      D=as.matrix(daisy(data, metric = "gower"))
    }
  }
  
  if(!is.null(kernel)){
    if(!require(kernlab)){
      install.packages("kernlab")
      library(kernlab)
    } 
    D=matrix(data=NA,nrow=dim(data)[1],ncol=dim(data)[1])
    d=c()
    if(substr(kernel, 1, 1)=="g"){
      sigma=as.numeric(substr(kernel, 5, nchar(kernel)))

      K.rbfdot<-function(x,y){
        rbf <- rbfdot(sigma = sigma)
        dist = 2*(1-rbf(x,y))
        sqrt(dist)
      }
      for(i in 1:dim(data)[1]){
        d=c()
        for(j in 1:dim(data)[1]){
          d=c(d,K.rbfdot(data[i,],data[j,]))
          D[i,j]=d[j]
        }
      }
    } else if(substr(kernel, 1, 1)=="p"){
      degree=as.numeric(substr(kernel, 5, nchar(kernel)))
      K.polydot<-function(x,y){
        lap <- polydot(degree = degree, scale = 1, offset = 1)
        dist = lap(x,x)+lap(y,y)-2*lap(x,y)
        sqrt(dist)
      }
      for(i in 1:dim(data)[1]){
        d=c()
        for(j in 1:dim(data)[1]){
          d=c(d,K.polydot(data[i,],data[j,]))
          D[i,j]=d[j]
        }
      }
    } else if(substr(kernel, 1, 1)=="l"){
      sigma=as.numeric(substr(kernel, 5, nchar(kernel)))
      
      K.laplacedot<-function(x,y){
        lap <- laplacedot(sigma = sigma)
        dist = 2*(1-lap(x,y))
        sqrt(dist)
      }
      for(i in 1:dim(data)[1]){
        d=c()
        for(j in 1:dim(data)[1]){
          d=c(d,K.laplacedot(data[i,],data[j,]))
          D[i,j]=d[j]
        }
      }
    }
  }
    
  comb = function(n, x) {
    return(factorial(n) / (factorial(x) * factorial(n-x)))
  }
  
  ind<-function(x,j){
    delta=0;
    if(x%in%g[,j]){delta=1}
    delta
  }
  
  h=matrix(data=0,nrow=dim(data)[1],ncol=k)
  g=matrix(data=0,nrow=dim(data)[1],ncol=k)
  s<-rep(0,k)
  for(c in 1:k){
    for(i in 1:dim(data)[1]){
      if(Class[i]==c){g[i,c]=i; s[c]<-s[c]+1}
    }
  }
  #S = prod(comb(s,lambda))
  
  row=1
  R=rep(0,k)
  ER=rep(0,k)
  mat=matrix(data=0,nrow=k,ncol=k)
  while(row<=(dim(h)[1])){
    h=matrix(data=0,nrow=dim(data)[1],ncol=k)
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
          #if(sum(vec.k==0)>0){vec.k[which(vec.k==0)]=1}
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
    
    R[j]=R[j]+1
    
    if(row%in%g[,j]){ER[j]=ER[j]+1}
    else{
      for(f in (1:k)[-j]){
        if(row%in%g[,f]){mat[j,f]=mat[j,f]+1}
      }
    }
    
    row<-row+1  
  }
  cm=diag(ER)+mat
  if(lev){dimnames(cm)=cm.levels}
  if(missing(prior)){prior="1/k"}
  
  if(class(prior)=="character"){
    
    if(prior=="1/k"){
      prob=matrix(ncol=k,nrow=k)
      pRE=matrix(ncol=k,nrow=k)
      for(i in 1:k){
        for(j in 1:k){
          pRE[i,j]=cm[i,j]/sum(cm[,j])
        }
      }
      for(i in 1:k){
        for(j in 1:k){
          prob[i,j]=pRE[i,j]/sum(pRE[i,])
        }
      }
      assignment = numeric(k)
      class.assign = numeric(k)
      for(i in 1:k){
        assignment[i] = max(prob[i,])
        class.assign[i]=which.max(prob[i,])
      }
    
      } else if(prior=="ni/N"){
      pi=as.numeric(table(Class))/dim(data)[1]
      prob=matrix(ncol=k,nrow=k)
      for(i in 1:k){
        for(j in 1:k){
          prob[i,j]=cm[i,j]/sum(cm[i,])
        }
      }
      assignment = numeric(k)
      class.assign = numeric(k)
      for(i in 1:k){
        assignment[i] = max(prob[i,])
        class.assign[i]=which.max(prob[i,])
      }
    }
  }
  

  if(class(prior)=="numeric"){
    pER=matrix(ncol=k,nrow=k)
    pRE=matrix(ncol=k,nrow=k)
    for(i in 1:k){
      for(j in 1:k){
        pRE[i,j]=cm[i,j]/sum(cm[,j])
        pRE[i,j]=pRE[i,j]*prior[j]
      }
    }
    
    for(i in 1:k){
      for(j in 1:k){
        pER[i,j]=pRE[i,j]/sum(pRE[i,])
      }
    }
    assignment = numeric(k)
    class.assign = numeric(k)
    for(i in 1:k){
      assignment[i] = max(pER[i,])
      class.assign[i]=which.max(pER[i,])
    }
  }
  
  
  
  if(missing(info.pred)||info.pred=="no"){
    return(list(s=s,lambda=lambda,Ri=R,Confusion_matrix=cm,prob_post=prob,assignment=assignment,class.assign=class.assign))
    
  } else if(lev){return(list(input=input,D=D,k=k,g=g,comb=comb,ind=ind,
                      K.rbfdot=K.rbfdot,K.laplacedot=K.laplacedot,K.polydot=K.polydot,
                      C.pool=C.pool,Ri=R,Confusion_matrix=cm,cm.levels=cm.levels,prob_post=prob,assignment=assignment,class.assign=class.assign))
  } else{
    return(list(input=input,D=D,k=k,g=g,comb=comb,ind=ind,
                K.rbfdot=K.rbfdot,K.laplacedot=K.laplacedot,K.polydot=K.polydot,
                C.pool=C.pool,Ri=R,Confusion_matrix=cm,prob_post=prob,assignment=assignment,class.assign=class.assign))
    }
  
} 


