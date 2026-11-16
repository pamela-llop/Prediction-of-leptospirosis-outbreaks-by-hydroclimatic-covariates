

# Librarys and sources

library(forecast)
library(tseries)
library(TSA)
library(readxl)
library(tidyverse)
library(reshape2)
library(Hmisc)

# files that are required

source("~/CODES/models.R") 


# load tha data: copy the dataset file into CODES folder
# data has four columns: 
# column 1: the response (leptospirosis cases)
# columns 2 to 4: covariates (precipitation, hydrometric level, ONI)

source("~/CODES/data")

# lags of precipitation (lp), hydrometric level (lhl) and ONI (lo)
# after choosing the optimal lags, use them here (no lags is 0)
lp = # to complete
lhl =  # to complete
lo = # to complete

data[,2] = Lag(data[,2],lp)
data[,3] = Lag(data[,3],lhl)
data[,4] = Lag(data[,4],lo)

maxlag = max(lp,lhl,lo)
N = dim(data)[1]



# in what follows it is assumed that the data is available from 
# 2009 to 2020

# 2009 to 2019 is used to train the methods
# 2020 is used to compute the error



# created the functional data 

# when the covariates are lagged, the first year 2009 is deleted
if(maxlag!=0){
dat = data[-c(1:12),]
} else {
dat = data
}

tau=12 # months
n=dim(dat)[1]/tau # years in the data set
NN=n*tau
Z=matrix(0,n,tau) 
x1=matrix(0,n,tau) 
x2=matrix(0,n,tau) 
x3=matrix(0,n,tau) 

for (i in 1:n){
  Z[i,]=dat$Cases[(tau*(i-1)+1):(tau*i)]
  x1[i,]=dat$Precip[(tau*(i-1)+1):(tau*i)]
  x2[i,]=dat$Hydro_Level[(tau*(i-1)+1):(tau*i)]
  x3[i,]=dat$ONI[(tau*(i-1)+1):(tau*i)]
}


# Predictions ----------------------------------------------------------------------

# SFPLR 

pred_SFPLR=NULL
hopt=NULL
qopt=NULL
hvec=seq(1,10,0.05)

for (i in 1:tau){
  pred=sfplr(Z[2:(n-1),], x1[2:(n-1),], x3[2:(n-1),], 
                 month=i,hvec=hvec, semimetric = "PCA")
  pred_SFPLR[i]=pred[[1]]
  hopt[i]=pred[[2]]
  qopt[i]=pred[[3]]
}

pred_SFPLR[pred_SFPLR<0]=0

pred_SFPLR = round(pred_SFPLR)



# ARIMA

train = 2019 # use the last year of your data as training sample
serie=ts(data$Cases[1:(N-12)],frequency = 12, start =c(2009,01), 
         end = c(train,12)) 
plot(serie)
decompose(serie, type = c("additive", "multiplicative"), filter = NULL) #descomposicion usando medias moviles
plot(decompose(serie, type = c("additive", "multiplicative"), filter = NULL),col="grey")

acf(serie)
pacf(serie)

best.arima.elec=get.best.arima( serie,maxord = c(2,2,2,2,2,2))
best.fit.elec=best.arima.elec[[2]]
acf(resid(best.fit.elec) )
orders=best.arima.elec [[3]]
ts.plot(cbind( window(serie,start = c(2009,01)),predict(best.fit.elec,12)$pred ), lty = 1:2)

ar_pred=round(as.numeric(predict(best.fit.elec,12)$pred))
ar_pred[ar_pred<0]=0


# ARIMAX 

# to running ARIMAX method, the covariates values must be predicted (we use ARIMA)

# Precipitation prediction

# when the covariates are lagged, the first year 2009 is deleted
prec_serie=ts(dat$Precip[1:(NN-12)],frequency = 12, start = c(2010,1), 
              end = c(train,12)) 
best.arima.elec=get.best.arima(prec_serie,maxord = c(2,2,2,2,2,2))
best.fit.elec=best.arima.elec[[2]]
acf( resid(best.fit.elec) )

ts.plot( cbind( window(prec_serie,start = c(2010,1)),
                predict(best.fit.elec,12)$pred ), lty = 1:2)

prec_predic=as.numeric(predict(best.fit.elec,12)$pred)
prec_predic[prec_predic<0]=0


# Hydrometric Level prediction

# when the covariates are lagged, the first year 2009 is deleted
nivel_serie=ts(dat$Hydro_Level[1:(NN-12)],frequency = 12, start = c(2010,1), 
               end = c(train,12)) 
best.arima.elec=get.best.arima(nivel_serie,maxord = c(2,2,2,2,2,2))
best.fit.elec=best.arima.elec[[2]]
acf( resid(best.fit.elec) )

ts.plot( cbind(window(nivel_serie,start = c(2010,1)),
               predict(best.fit.elec,12)$pred ), lty = 1:2)

nivel_predic=as.numeric(predict(best.fit.elec,12)$pred)
nivel_predic[nivel_predic<0]=0


# ONI prediction 

# when the covariates are lagged, the first year 2009 is deleted
ONI_serie=ts(dat$ONI[1:(NN-12)],frequency = 12, start = c(2010,1), 
              end = c(train,12)) 
best.arima.elec=get.best.arima(ONI_serie,maxord = c(2,2,2,2,2,2))
best.fit.elec=best.arima.elec[[2]]
acf( resid(best.fit.elec))

ts.plot( cbind( window(ONI_serie,start = c(2010,1)),
                predict(best.fit.elec,12)$pred ), lty = 1:2)

ONI_predic=as.numeric(predict(best.fit.elec,12)$pred)


# ARIMAX prediction

# when the covariates are lagged, the first year 2009 is deleted
serie=ts(dat$Cases[1:(NN-12)],frequency = 12, start = c(2010,1), 
         end = c(train,12))

xreg1=data.frame(dat$Precip[1:(NN-12)],dat$Hydro_Level[1:(NN-12)],
                 dat$ONI[1:(NN-12)])
ajustex1=arimax(serie, order = orders[1:3], xreg=xreg1, 
                include.mean = TRUE, seasonal = list(order = orders[4:6], period = 12))
newxreg1=data.frame(prec_predic,nivel_predic,ONI_predic)
predx1=predict(ajustex1, n.ahead = 12, newxreg=newxreg1)

arx_pred=round(as.numeric(predx1$pred))
arx_pred[arx_pred<0]=0


# NW-1K 

# the series to predict is the number of leptospirosis cases (column 1 of data)
y_cases = as.numeric(data$Cases)

TT=N-12
k=1:12
kk=length(k)
r=1:12
rr=length(r)
h=seq(0.1,30,0.5)
hh=length(h)
ECM=matrix(0,rr,hh)
r_opt=NULL
h_opt=NULL

for (l in 1:kk){
  for (i in 1:rr){
    for (j in 1:hh){
      hat=NULL
      cont=0
      for (t in (r[i]+k[l]): (TT-k[l])){
        cont=cont+1
        serie=y_cases[1:t] 
        hat[cont]=one.kernel(serie, k=k[l], r=r[i], ker = gaussian,h=h[j])
      }
      ECM[i,j]=mean((round(hat)-y_cases[(r[i]+2*k[l]):TT])^2)
    }
  }
  r_opt[l]=r[which(ECM==min(ECM),arr.ind=TRUE)[1,1]]
  h_opt[l]=h[which(ECM==min(ECM),arr.ind=TRUE)[1,2]]
  print(l)
}

h_opt
r_opt


hat=NULL
for (j in 1:kk){
  hat[j]=one.kernel(y_cases[1:(N-12)], k[j], r = r_opt[j], ker = gaussian, h = h_opt[j])
}
kpred=round(hat)



# NW-2K 

TT=N-12
k=1:12
kk=length(k)
r=1:12
rr=length(r)
h1=seq(0.1,20,1) 
hh1=length(h1)
h2=seq(1,12,1)
hh2=length(h2)
ECM=array(NA,c(hh1,hh2,rr))
r_opt=NULL
h1_opt=NULL
h2_opt=NULL
ts=1:TT

for (l in 1:kk){
  for (i in 1:rr){
    for (j in 1:hh1){
      for (p in 1:hh2){
        hat=NULL
        cont=0
        for (t in (r[i]+k[l]): (TT-k[l])){
          cont=cont+1
          serie=y_cases[1:t] 
          hat[cont]=two.kernel(serie, ts, k[l], r[i], ker1 = gaussian, 
                               ker2 = gaussian, h1=h1[j], h2=h2[p]) 
        }
        ECM[j,p,i]=mean((round(hat)-y_cases[(r[i]+2*k[l]):TT])^2)
      }
    }
  }
  r_opt[l]=r[which(ECM==min(ECM),arr.ind=TRUE)[1,3]]
  h1_opt[l]=h1[which(ECM==min(ECM),arr.ind=TRUE)[1,1]]
  h2_opt[l]=h2[which(ECM==min(ECM),arr.ind=TRUE)[1,2]]
  print(l)
}


kkpred=NULL
for (j in 1:kk){
  kkpred[j]=two.kernel(y_cases[1:(N-12)], ts, k[j], r_opt[j], ker1 = gaussian, 
                          ker2 = gaussian, h1 = h1_opt[j], h2=h2_opt[j])
}

kkpred=round(kkpred)



# save the results 

real=Z[n,] # testing dataset
real 
save_list = list(pred_SFPLR,ar_pred, arx_pred, kpred, kkpred, real)
save(save_list, file = "results.RData")


# RMSE -----------------------------------------------------------------------

# compute the RMSE using the predictions and the testing dataset

er2_sfplr=mean((round(pred_SFPLR)-real)^2)
sqrt(er2_sfplr)

er2_arima=mean((round(ar_pred)-real)^2)
sqrt(er2_arima) 

er2_arimax=mean((round(arx_pred)-real)^2)
sqrt(er2_arimax)

er2_k=mean((round(kpred)-real)^2)
sqrt(er2_k)

er2_kk=mean((round(kkpred)-real)^2)
sqrt(er2_kk)



