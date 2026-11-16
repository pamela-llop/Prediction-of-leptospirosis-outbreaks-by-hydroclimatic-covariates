

###############################################################################
# ARIMA     
###############################################################################


get.best.arima <- function(x.ts, maxord)
{
  best.aic <- 1e8
  n <- length(x.ts)
  for (p in 0:maxord[1]) for(d in 0:maxord[2]) for(q in 0:maxord[3])
    for (P in 0:maxord[4]) for(D in 0:maxord[5]) for(Q in 0:maxord[6])
    {
      fit <- arima(x.ts, order = c(p,d,q),
                   seas = list(order = c(P,D,Q),
                               frequency(x.ts)), method = "CSS")
      fit.aic <- -2 * fit$loglik + (log(n) + 1) * length(fit$coef)
      if (fit.aic < best.aic)
      {
        best.aic <- fit.aic
        best.fit <- fit
        best.model <- c(p,d,q,P,D,Q)
      }
    }
  list(best.aic, best.fit, best.model)
}



###############################################################################
# SFPLR    
##############################################################################

sfplr=function(Z,x1,x2,x3=NULL,month,hvec,semimetric){
  #month indicates which we wats to predict 
  #Z matrix of functional response data
  #x1, x2 y x3 matrix of covariates (same size of Z)
  #h bandwith
  
  library(fda.usc)
  
  n=dim(Z)[1] 
  tau=dim(Z)[2] 
  j=month
  hh=length(hvec)
  
  G=Z[1:n,] 
  if (is.null(x3)){
    X=cbind(x1[1:n,j],x2[1:n,j]) 
  }else{
    X=cbind(x1[1:n,j],x2[1:n,j],x3[1:n,j]) 
  }
  
  if (semimetric=="PCA"){
    qvec=1:12
    qq=length(qvec)
    CV=matrix(NA,hh,qq)
    
    for (l in 1:qq){
      for(k in 1:hh){
        W=pesos(Z[1:n,], ker = 'gaussian', hvec[k],loocv=TRUE,semimetric,q=qvec[l]) 
        
        Xpelito=(diag(1,n)-W[1:n,1:n])%*%X 
        Gpelito=(diag(1,n)-W[1:n,1:n])%*%G[,j]
        
        inv=solve(t(Xpelito)%*%Xpelito) 
        beta=inv%*%t(Xpelito)%*%Gpelito 
        
        diag(W)=0
        m=W%*%(G[,j]-X%*%beta) 
        pred=X%*%beta+m 
        CV[k,l]=mean((G[,j]-pred)^2)
        
      }
      plot(hvec,CV[,l])
    }
    coor=which(CV==min(CV),arr.ind=TRUE)
    hopt=hvec[coor[1,1]] 
    qopt=qvec[coor[1,2]] 
    
    W=pesos(Z[1:n,], ker = 'gaussian', hopt, loocv=TRUE, semimetric, q=qopt) 
    
    Xpelito=(diag(1,n)-W[1:n,1:n])%*%X 
    Gpelito=(diag(1,n)-W[1:n,1:n])%*%G[,j] 
    
    inv=solve(t(Xpelito)%*%Xpelito) 
    beta=inv%*%t(Xpelito)%*%Gpelito 
    
    m=W[n,]%*%(G[,j]-X%*%beta)
    
    pred=X[n,]%*%beta+m 
    
    return(list(pred, hopt, qopt))
  }
  else {
    CV=NULL
    
    for(k in 1:hh){
      W=pesos(Z[1:n,], ker = 'gaussian', hvec[k],loocv=TRUE,semimetric) 
      
      Xpelito=(diag(1,n)-W[1:n,1:n])%*%X 
      Gpelito=(diag(1,n)-W[1:n,1:n])%*%G[,j] 
      
      inv=solve(t(Xpelito)%*%Xpelito) 
      beta=inv%*%t(Xpelito)%*%Gpelito 
      
      diag(W)=0
      m=W%*%(G[,j]-X%*%beta) 
      pred=X%*%beta+m 
      CV[k]=mean((G[,j]-pred)^2)
      
    }
    hopt=hvec[which.min(CV)]
    plot(hvec,CV)
    
    
    W=pesos(Z[1:n,], ker = 'gaussian', hopt, loocv=TRUE, semimetric)
    
    Xpelito=(diag(1,n)-W[1:n,1:n])%*%X 
    Gpelito=(diag(1,n)-W[1:n,1:n])%*%G[,j] 
    
    inv=solve(t(Xpelito)%*%Xpelito) 
    beta=inv%*%t(Xpelito)%*%Gpelito
    
    m=W[n,]%*%(G[,j]-X%*%beta)
    
    pred=X[n,]%*%beta+m 
    
    return(list((pred, hopt, qopt))
  }
}



###############################################################################
# NW-1K     
##############################################################################


# kernels
epanechnikov = function(x){
  result = (3/4) * (1 - x^2) * (abs(x) <= 1) 
  result}
biweight = function(x){
  result = (15/16) * (1 - x^2)^2 * (abs(x) <= 1) 
  result}
triweight = function(x){
  result = (35/32) * (1 - x^2)^3 * (abs(x) <= 1) 
  result}
triangular = function(x){
  result = (1 - abs(x)) * (abs(x) <= 1) 
  result}
uniform = function(x){
  result = (1/2) * (abs(x) <= 1) 
  result}
gaussian = function(x){
  result = (1/sqrt(2*pi)) * exp(-0.5*x^2)
  result}


pesos = function(X, ker, h,loocv=TRUE,semimetric,q=NULL){
  # X matrix de datos nxp
  # ker elije un tipo de nucleo de la lista
  # h ancho de banda
  # loocv=TRUE indica leave one out
  # q cantidad de componentes prncipales en la semimetrica PCA
  
  if (semimetric=="PCA"){
    DX=semimetric.pca(X,q=q) 
  }
  else if (semimetric=="L2"){
    DX=as.matrix(dist(X)) 
  }
  
  K = gaussian(DX/h)
  
  if(loocv==TRUE){
    diag(K) = 0
  }
  
  Den = apply(K, 2, sum)
  K[Den==0] = 0
  W = K/Den
  
  return(W)
}

norm_vec <- function(x) sqrt(sum(x^2))



one.kernel <- function(serie, k, r, ker = gaussian, h = ((4/3)^{1/5})*sd(serie)*length(serie)^{-1/5}){
  # serie vector of data lenght T
  # k lag to predict (hat Y_{T+k}) 
  # r Markovian parameter 
  # ker choose one of the list
  # h bandwith, real number
  
  # Predice Y_{T+k}
  
  TT = length(serie)
  Y_T = serie[seq(TT-r+1,TT)] 
  
  Y = serie[seq(r+k,TT)] 
  
  # distances
  DX = rep(NA, TT-k-r+1)
  cont = 0 
  for (t in r:(TT-k)){
    cont = cont + 1
    
    Y_t = serie[seq(t-r+1,t)]
    
    DX[cont] = norm_vec(Y_T - Y_t)
  }
  
  K = ker(DX/h)
  
  Den = sum(K) 
  Num = sum(K*Y)
  yhat = Num/Den
  yhat[Den==0] = 0
  
  return(yhat)
}



###############################################################################
# NW-2K    
##############################################################################


two.kernel <- function(serie, ts, k_vec, r, ker1 = gaussian, ker2 = gaussian, h1, h2){
  # serie vector of data lenght T
  # k lag to predict (hat Y_{T+k}) can be a vector
  # r Markovian parameter 
  # ker1 and ker2 choose one of the list
  # h bandwith, real number
  
  
  TT = length(serie)
  Y_T = serie[seq(TT-r+1,TT)]
  kk = length(k_vec)
  
  yhat = rep(NA, kk)
  for (k in 1:kk){
    Y = serie[seq(r+k,TT)]
    
    # distances
    DX = rep(NA, TT-k-r+1)
    Dt = rep(NA, TT-k-r+1)
    cont = 0 
    for (t in r:(TT-k)){
      cont = cont + 1
      
      Y_t = serie[seq(t-r+1,t)]
      
      DX[cont] = norm_vec(Y_T - Y_t)
      
      Dt[cont] = TT - t
      
    }
    
    K1 = ker1(DX/h1)
    K2 = ker2(Dt/h2)
    
    Den = sum(K1*K2) 
    Num = sum(K1*K2*Y)  
    yhat_aux = Num/Den
    yhat_aux[Den==0] = 0
    
    yhat[k] = yhat_aux
  }
  
  return(yhat)
}
